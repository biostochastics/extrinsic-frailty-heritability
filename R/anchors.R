# ===========================================================================
# anchors.R — σ_γ bridge, m_ex splitting, negative ρ sensitivity
# ===========================================================================

#' σ_γ bridge experiment (Vulnerability #2)
#'
#' Simulates competing-risk infection deaths, computes tetrachoric correlation,
#' and finds σ_γ that reproduces Obel et al.'s h²_L ≈ 0.40.
#'
#' @param sigma_theta_true True sigma_theta
#' @param params Parameter list
#' @return List with bridge_sigma_gamma and curve data
run_sigma_gamma_bridge <- function(sigma_theta_true, params) {
  sigma_gamma_grid <- seq(0.1, 1.5, by = 0.025)

  results <- lapply(sigma_gamma_grid, function(sg) {
    compute_infection_death_h2(
      sigma_gamma = sg,
      rho = params$RHO_DEF,
      sigma_theta = sigma_theta_true,
      m_ex = params$M_EX_HIST,
      n_pairs = params$N_PAIRS,
      cutoff = params$CUTOFF_AGE,
      seed = params$MASTER_SEED + round(sg * 10000)
    )
  })

  h2_vals <- sapply(results, `[[`, "h2_L")
  curve <- data.frame(
    sigma_gamma = sigma_gamma_grid,
    h2_infection = h2_vals,
    stringsAsFactors = FALSE
  )

  # Interpolate to find sigma_gamma matching h²_L = 0.40
  bridge_sg <- tryCatch({
    approx(h2_vals, sigma_gamma_grid, xout = 0.40)$y
  }, error = function(e) NA_real_)

  list(
    bridge_sigma_gamma = bridge_sg,
    curve = curve,
    target_h2_L = 0.40
  )
}

#' Compute infection-death liability h² for a given σ_γ
#'
#' Uses competing risks: T_intrinsic (Gompertz) vs T_extrinsic (Exponential).
#' Computes tetrachoric correlation of binary (died_extrinsic) for MZ and DZ.
#'
#' @param sigma_gamma SD of log-extrinsic frailty
#' @param rho Pleiotropy correlation
#' @param sigma_theta True sigma_theta
#' @param m_ex Historical extrinsic hazard
#' @param n_pairs Number of twin pairs
#' @param cutoff Minimum age
#' @param seed RNG seed
#' @return List with h2_L, r_tet_mz, r_tet_dz
compute_infection_death_h2 <- function(sigma_gamma, rho, sigma_theta,
                                        m_ex, n_pairs, cutoff, seed,
                                        b = PARAMS$B_GOMP) {
  set.seed(seed)
  shift <- sigma_gamma^2 / 2

  compute_for_zygosity <- function(r_g) {
    params_tw <- gen_twin_params_gm(n_pairs, r_g, sigma_theta, sigma_gamma,
                                     rho, gamma_heritable = TRUE)

    c1 <- m_ex * exp(params_tw$gamma1 - shift)
    c2 <- m_ex * exp(params_tw$gamma2 - shift)

    # Intrinsic-only lifespans (c_ext = 0)
    T_int1 <- sim_lifespan_gm(params_tw$theta1, b, c_ext = 0)
    T_int2 <- sim_lifespan_gm(params_tw$theta2, b, c_ext = 0)

    # Extrinsic event times
    T_ext1 <- rexp(n_pairs, rate = pmax(c1, 1e-20))
    T_ext2 <- rexp(n_pairs, rate = pmax(c2, 1e-20))

    # Actual death = min
    T_death1 <- pmin(T_int1, T_ext1)
    T_death2 <- pmin(T_int2, T_ext2)

    # Both survived cutoff
    keep <- T_death1 > cutoff & T_death2 > cutoff

    # Died from extrinsic cause (among those surviving cutoff who eventually die)
    died_ext1 <- T_ext1 < T_int1
    died_ext2 <- T_ext2 < T_int2

    list(
      died_ext1 = died_ext1[keep],
      died_ext2 = died_ext2[keep],
      n_keep = sum(keep)
    )
  }

  mz <- compute_for_zygosity(1.0)
  dz <- compute_for_zygosity(0.5)

  # Tetrachoric correlation from 2x2 table
  r_tet_mz <- tetrachoric_2x2(mz$died_ext1, mz$died_ext2)
  r_tet_dz <- tetrachoric_2x2(dz$died_ext1, dz$died_ext2)

  list(
    h2_L = 2 * (r_tet_mz - r_tet_dz),
    r_tet_mz = r_tet_mz,
    r_tet_dz = r_tet_dz,
    n_mz = mz$n_keep,
    n_dz = dz$n_keep
  )
}

#' Tetrachoric correlation from two binary vectors
#'
#' Uses maximum likelihood via optim on the bivariate normal threshold model.
#'
#' @param x Binary vector (logical or 0/1)
#' @param y Binary vector (logical or 0/1)
#' @return Tetrachoric correlation estimate
tetrachoric_2x2 <- function(x, y) {
  x <- as.integer(x)
  y <- as.integer(y)
  n <- length(x)

  # 2x2 cell counts
  a <- sum(x == 1 & y == 1)
  b <- sum(x == 1 & y == 0)
  cc <- sum(x == 0 & y == 1)
  d <- sum(x == 0 & y == 0)

  # Edge cases
  if (a == 0 || b == 0 || cc == 0 || d == 0) {
    # Add 0.5 continuity correction
    a <- a + 0.5; b <- b + 0.5; cc <- cc + 0.5; d <- d + 0.5
    n <- n + 2
  }

  # Thresholds
  p1 <- (a + b) / n
  p2 <- (a + cc) / n
  t1 <- qnorm(1 - p1)
  t2 <- qnorm(1 - p2)

  # MLE: find rho that maximizes bivariate normal probability
  neg_loglik <- function(rho) {
    rho <- max(min(rho, 0.999), -0.999)
    Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    # P(X>t1, Y>t2) using pmvnorm
    p11 <- tryCatch(
      mvtnorm::pmvnorm(lower = c(t1, t2), upper = c(Inf, Inf), sigma = Sigma)[1],
      error = function(e) NA
    )
    if (is.na(p11) || p11 <= 0) return(1e10)
    p10 <- pnorm(-t1) - p11  # P(X>t1, Y<t2)
    p01 <- pnorm(-t2) - p11  # P(X<t1, Y>t2)
    p00 <- 1 - p11 - p10 - p01

    # Clamp
    p10 <- max(p10, 1e-15)
    p01 <- max(p01, 1e-15)
    p00 <- max(p00, 1e-15)
    p11 <- max(p11, 1e-15)

    -(a * log(p11) + b * log(p10) + cc * log(p01) + d * log(p00))
  }

  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for tetrachoric correlation.")
  }
  opt <- optim(par = 0, fn = neg_loglik, method = "Brent",
               lower = -0.999, upper = 0.999)
  opt$par
}

#' m_ex heritable fraction sensitivity (Vulnerability #6)
#'
#' Splits m_ex = m_inf + m_other. Only infection fraction gets γ.
#'
#' @param sigma_theta_true True sigma_theta
#' @param params Parameter list
#' @return Data frame with frac_heritable, bias, sigma_infl_pct
run_mex_split_sensitivity <- function(sigma_theta_true, params) {
  frac_grid <- seq(0, 1, by = 0.05)
  oracle <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0,
                         rng_seed = params$MASTER_SEED)

  results <- lapply(frac_grid, function(f_inf) {
    seed_base <- params$MASTER_SEED + 50000L + round(f_inf * 1000)
    res <- sim_twin_h2_split_mex(
      sigma_theta_true,
      m_inf = f_inf * params$M_EX_HIST,
      m_other = (1 - f_inf) * params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF,
      rho = params$RHO_DEF,
      rng_seed = seed_base
    )

    # Calibrate with total m_ex (Shenhar's model sees total m_ex)
    sigma_fit <- calibrate_sigma_theta(res$r_mz, params$M_EX_HIST)
    corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
                         rng_seed = seed_base + 5000L)

    data.frame(
      frac_heritable = f_inf,
      bias = corr$h2 - oracle$h2,
      bias_pp = 100 * (corr$h2 - oracle$h2),
      sigma_infl_pct = 100 * (sigma_fit / sigma_theta_true - 1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Simulate twin h² with split m_ex (heritable + constant components)
#'
#' c_i = m_other + m_inf * exp(gamma_i - sigma^2/2)
#'
#' @param sigma_theta True sigma_theta
#' @param m_inf Heritable infection-like component
#' @param m_other Constant non-heritable component
#' @param sigma_gamma SD of log-extrinsic frailty
#' @param rho Pleiotropy correlation
#' @param rng_seed RNG seed
#' @return List with r_mz, r_dz, h2
sim_twin_h2_split_mex <- function(sigma_theta, m_inf, m_other,
                                   sigma_gamma = PARAMS$SIGMA_GAMMA_DEF,
                                   rho = PARAMS$RHO_DEF,
                                   rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)

  n <- PARAMS$N_PAIRS
  b <- PARAMS$B_GOMP
  t_max <- PARAMS$T_MAX
  cutoff <- PARAMS$CUTOFF_AGE
  shift <- sigma_gamma^2 / 2

  # MZ
  params_mz <- gen_twin_params_gm(n, r_g = 1.0, sigma_theta, sigma_gamma,
                                    rho, gamma_heritable = TRUE)
  c_mz1 <- m_other + m_inf * exp(params_mz$gamma1 - shift)
  c_mz2 <- m_other + m_inf * exp(params_mz$gamma2 - shift)
  L_mz1 <- sim_lifespan_gm(params_mz$theta1, b, c_mz1, t_max)
  L_mz2 <- sim_lifespan_gm(params_mz$theta2, b, c_mz2, t_max)

  # DZ
  params_dz <- gen_twin_params_gm(n, r_g = 0.5, sigma_theta, sigma_gamma,
                                    rho, gamma_heritable = TRUE)
  c_dz1 <- m_other + m_inf * exp(params_dz$gamma1 - shift)
  c_dz2 <- m_other + m_inf * exp(params_dz$gamma2 - shift)
  L_dz1 <- sim_lifespan_gm(params_dz$theta1, b, c_dz1, t_max)
  L_dz2 <- sim_lifespan_gm(params_dz$theta2, b, c_dz2, t_max)

  keep_mz <- (L_mz1 > cutoff) & (L_mz2 > cutoff)
  keep_dz <- (L_dz1 > cutoff) & (L_dz2 > cutoff)

  r_mz <- cor(L_mz1[keep_mz], L_mz2[keep_mz])
  r_dz <- cor(L_dz1[keep_dz], L_dz2[keep_dz])

  list(r_mz = r_mz, r_dz = r_dz, h2 = 2 * (r_mz - r_dz))
}

#' Negative ρ sensitivity (Task 10)
#'
#' Shows that negative ρ → downward bias (falsifiable prediction).
#'
#' @param sigma_theta_true True sigma_theta
#' @param params Parameter list
#' @return Data frame with rho, bias, bias_pp
run_negative_rho_sensitivity <- function(sigma_theta_true, params) {
  rho_grid <- seq(-0.6, 0.8, by = 0.05)
  oracle <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0,
                         rng_seed = params$MASTER_SEED)

  results <- lapply(rho_grid, function(rho_val) {
    seed_base <- params$MASTER_SEED + 60000L + round((rho_val + 1) * 1000)

    obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = rho_val,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = seed_base
    )

    sigma_fit <- calibrate_sigma_theta(obs$r_mz, params$M_EX_HIST)

    corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
      rng_seed = seed_base + 5000L
    )

    data.frame(
      rho = rho_val,
      bias = corr$h2 - oracle$h2,
      bias_pp = 100 * (corr$h2 - oracle$h2),
      sigma_fit = sigma_fit,
      sigma_infl_pct = 100 * (sigma_fit / sigma_theta_true - 1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}
