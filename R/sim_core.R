# ===========================================================================
# sim_core.R — Gompertz-Makeham lifespan simulator
# ===========================================================================
# Lambert W closed-form solver replaces numerically unstable Newton-Raphson.
# Bisection+Newton fallback for edge cases and MGG model.
# SR Euler-Maruyama solver compiled via Rcpp for ~5x speedup.
# ===========================================================================

# Rcpp SR solver: lazy-compile on first use per process
.sr_cpp_env <- new.env(parent = emptyenv())
.sr_cpp_env$loaded <- FALSE

.ensure_sr_cpp <- function() {
  if (!.sr_cpp_env$loaded) {
    cpp_file <- file.path(getwd(), "src", "sim_lifespan_sr.cpp")
    if (file.exists(cpp_file)) {
      Rcpp::sourceCpp(cpp_file, verbose = FALSE, rebuild = FALSE)
    }
    .sr_cpp_env$loaded <- TRUE
  }
}

#' Simulate lifespans from individual Gompertz-Makeham hazards (vectorized)
#'
#' Hazard:  mu_i(t) = c_i + exp(theta_i) * exp(b*t)
#' Cumul:   H_i(t) = c_i*t + exp(theta_i)/b * (exp(b*t) - 1)
#' Solve H_i(t) = -log(U_i) via Lambert W closed-form.
#'
#' @param theta Log-Gompertz intercept (vector)
#' @param b Gompertz slope (scalar)
#' @param c_ext Extrinsic hazard rate (scalar or vector)
#' @param t_max Maximum age (censoring boundary)
#' @param u Pre-generated uniform draws for CRN (NULL = generate fresh)
#' @return Vector of lifespans with attributes n_censored, max_error
sim_lifespan_gm <- function(theta, b, c_ext, t_max = 120,
                            u = NULL) {
  n <- length(theta)
  a <- exp(theta)
  c_ext <- rep_len(c_ext, n)
  if (is.null(u)) u <- runif(n)
  y <- -log(u)

  # Check feasibility: max cumulative hazard at t_max
  H_max <- c_ext * t_max + (a / b) * (exp(pmin(b * t_max, 18)) - 1)
  infeasible <- y > H_max

  # --- Lambert W closed-form for c > 0 ---
  # H(t) = c*t + (a/b)*(exp(bt)-1) = y
  # Solution: t = y/c + a/(bc) - (1/b) * W( (a/c) * exp(a/c + b*y/c) )
  #
  # Handle c ≈ 0 with pure Gompertz formula: t = (1/b)*log(1 + b*y/a).
  # For c > 0, use Lambert W with log-space Newton for overflow safety.
  is_pure_gomp <- c_ext < 1e-12
  t <- numeric(n)

  if (any(is_pure_gomp & !infeasible)) {
    idx <- which(is_pure_gomp & !infeasible)
    t[idx] <- (1 / b) * log1p(b * y[idx] / a[idx])
  }

  if (any(!is_pure_gomp & !infeasible)) {
    idx <- which(!is_pure_gomp & !infeasible)
    # Log-stable computation of Lambert W argument:
    # W_arg = (a/c) * exp(a/c + b*y/c)
    # log(W_arg) = log(a/c) + a/c + b*y/c
    log_W_arg <- log(a[idx] / c_ext[idx]) +
      a[idx] / c_ext[idx] +
      b * y[idx] / c_ext[idx]

    # Compute W(exp(z)) using two paths:
    # - Standard: exp() + lamW for z <= 500 (safe from overflow)
    # - Log-space: solve w + log(w) = z via Newton for z > 500
    W_val <- numeric(length(log_W_arg))
    large <- log_W_arg > 500

    if (any(!large)) {
      W_val[!large] <- lamW::lambertW0(exp(log_W_arg[!large]))
    }
    if (any(large)) {
      # Log-space Newton: W(exp(z)) satisfies w + log(w) = z
      # Avoids materializing exp(z), works for z up to ~1e15
      z <- log_W_arg[large]
      z[!is.finite(z)] <- 700  # guard against Inf from extreme a/c
      w <- z - log(z)                                  # asymptotic init
      for (iter in 1:3) {                               # 3 Newton steps
        w <- w - (w + log(w) - z) / (1 + 1 / w)        # quadratic convergence
      }
      W_val[large] <- w
    }

    t[idx] <- y[idx] / c_ext[idx] +
      a[idx] / (b * c_ext[idx]) -
      W_val / b
  }

  # Censoring
  t[infeasible] <- t_max
  t <- pmax(t, 0)
  t <- pmin(t, t_max)

  # Verification: check H(t) ≈ y for feasible cases
  # Only check non-extreme cases (bt < 18) where exp() doesn't overflow.
  max_err <- 0
  feasible <- which(!infeasible)
  if (length(feasible) > 0) {
    bt_raw <- b * t[feasible]
    checkable <- bt_raw < 18
    if (any(checkable)) {
      ci <- feasible[checkable]
      bt <- b * t[ci]
      H_check <- c_ext[ci] * t[ci] + (a[ci] / b) * (exp(bt) - 1)
      max_err <- max(abs(H_check - y[ci]))
      if (max_err > 0.01) {
        warning(sprintf("Lambert W max error: %.6f (n_checked=%d/%d)",
                        max_err, sum(checkable), length(feasible)))
      }
    }
  }

  attr(t, "n_censored") <- sum(infeasible)
  attr(t, "max_error") <- max_err
  t
}


#' Bisection+Newton hybrid solver (fallback)
#'
#' Used when Lambert W is unavailable or for models (like MGG) where the
#' cumulative hazard doesn't have a Lambert W inverse.
#'
#' @param H_fn Function(t) returning cumulative hazard H(t) (vectorized over t)
#' @param mu_fn Function(t) returning instantaneous hazard mu(t) (vectorized)
#' @param y Target values -log(U) (vector)
#' @param t_max Maximum age
#' @param n_bisect Number of bisection iterations for bracketing
#' @param n_newton Number of Newton refinement iterations
#' @return Vector of lifespans
solve_cumhaz_hybrid <- function(H_fn, mu_fn, y, t_max = 120,
                                n_bisect = 30, n_newton = 10) {
  n <- length(y)

  # Check feasibility
  H_max <- H_fn(rep(t_max, n))
  infeasible <- y > H_max

  # Bisection phase: bracket the root

  lo <- rep(0, n)
  hi <- rep(t_max, n)

  for (iter in seq_len(n_bisect)) {
    mid <- (lo + hi) / 2
    H_mid <- H_fn(mid)
    below <- H_mid < y
    lo[below] <- mid[below]
    hi[!below] <- mid[!below]
  }

  t <- (lo + hi) / 2

  # Newton refinement phase
  for (iter in seq_len(n_newton)) {
    H_t <- H_fn(t)
    mu_t <- mu_fn(t)
    mu_t <- pmax(mu_t, 1e-15)
    delta <- (H_t - y) / mu_t
    t <- t - delta
    t <- pmax(pmin(t, t_max), 0)
  }

  # Final residual check
  H_final <- H_fn(t)
  resid <- abs(H_final - y)
  feasible_idx <- which(!infeasible)
  max_resid <- if (length(feasible_idx) > 0) max(resid[feasible_idx]) else 0

  t[infeasible] <- t_max
  t <- pmax(t, 0)
  t <- pmin(t, t_max)

  attr(t, "n_censored") <- sum(infeasible)
  attr(t, "max_residual") <- max_resid

  if (max_resid > 0.01) {
    warning(sprintf("Hybrid solver max residual: %.6f", max_resid))
  }

  t
}


#' Simulate lifespans from MGG model via bisection+Newton
#'
#' Shenhar et al. Table S1 (Makeham-Gamma-Gompertz):
#'   Hazard: mu_i(t) = m_ex_i + a_i * exp(b_i * t) * s(t)
#'   where s(t) = exp(c) / (exp(c) + exp(b*t) - 1)
#'   SM coupling: b_i = q_i * b0
#'   Parameterization flag controls a_i mapping (see below).
#'
#' Cumulative: H(t) = m_ex*t + (a*exp(c)/b) * log((exp(c)+exp(bt)-1)/exp(c))
#'
#' @param q Individual scale parameters (vector)
#' @param m_ex Individual extrinsic hazard rates (scalar or vector)
#' @param a0 Baseline Gompertz intercept
#' @param b0 Baseline Gompertz slope
#' @param c_param Deceleration parameter
#' @param t_max Maximum age
#' @param u Pre-generated uniform draws for CRN (NULL = generate fresh)
#' @param sm_mapping SM coupling for a: "compensatory" (a=a0^q, correct)
#'   or "paper" (a=a0^(1/q), as written in paper but inverted exponent).
#'   Default "compensatory" matches Shenhar's published twin correlations.
#' @return Vector of lifespans
sim_lifespan_mgg <- function(q, m_ex,
                             a0 = NULL, b0 = NULL, c_param = NULL,
                             t_max = 120, u = NULL,
                             sm_mapping = "compensatory") {
  if (is.null(a0)) a0 <- 1e-5
  if (is.null(b0)) b0 <- 0.115
  if (is.null(c_param)) c_param <- 30

  n <- length(q)
  m_ex <- rep_len(m_ex, n)
  q <- pmax(q, 0.01)

  # Strehler-Mildvan coupling: b_i = q_i * b0
  # Paper states log(a) = log(a0)/q => a = a0^(1/q), but this is INVERTED.
  # Correct SM-compensatory form: log(a) = q*log(a0) => a = a0^q
  # With a0^q: higher q => lower a, higher b (true compensation, curves cross at t*)
  # With a0^(1/q): higher q => higher a AND higher b (no compensation)
  b <- q * b0
  if (sm_mapping == "compensatory") {
    a <- a0^q
  } else {
    a <- a0^(1 / q)
  }

  ec <- exp(c_param)

  if (is.null(u)) u <- runif(n)
  y <- -log(u)

  # Use hybrid solver with MGG-specific cumulative hazard
  H_fn <- function(t) {
    bt <- pmin(b * t, 700)
    ebt <- exp(bt)
    D <- ec + ebt - 1
    m_ex * t + (a * ec / b) * log1p((ebt - 1) / ec)
  }

  mu_fn <- function(t) {
    bt <- pmin(b * t, 700)
    ebt <- exp(bt)
    D <- ec + ebt - 1
    m_ex + a * ebt * ec / D
  }

  solve_cumhaz_hybrid(H_fn, mu_fn, y, t_max = t_max)
}


#' Simulate lifespans from SR damage model via Euler-Maruyama (Rcpp)
#'
#' dx = (eta*t - beta*x/(x+kappa)) dt + sqrt(2*epsilon) dW
#' Death when x > Xc OR extrinsic Poisson event
#'
#' @param Xc Individual death thresholds (vector)
#' @param m_ex Individual extrinsic hazard rates (scalar or vector)
#' @param eta Damage production rate
#' @param beta Max removal rate
#' @param epsilon Noise amplitude
#' @param kappa Half-saturation constant
#' @param dt Euler-Maruyama timestep
#' @param t_max Maximum age
#' @return Vector of lifespans
sim_lifespan_sr <- function(Xc, m_ex,
                            eta = NULL, beta = NULL,
                            epsilon = NULL, kappa = NULL,
                            dt = NULL, t_max = 120) {
  .ensure_sr_cpp()  # lazy-compile Rcpp on first call per process

  if (is.null(eta)) eta <- 0.66
  if (is.null(beta)) beta <- 62.96
  if (is.null(epsilon)) epsilon <- 51.83
  if (is.null(kappa)) kappa <- 0.50
  if (is.null(dt)) dt <- 1 / 52

  n <- length(Xc)
  m_ex <- rep_len(m_ex, n)
  Xc <- pmax(Xc, 0.01)

  # Derive a seed from R's current RNG state for the fast C++ path
  fast_seed <- sample.int(.Machine$integer.max, 1L)
  sim_lifespan_sr_cpp(Xc, m_ex, eta, beta, epsilon, kappa, dt, t_max,
                      rng_seed = fast_seed)
}
