# ===========================================================================
# calibrate.R — Bisection calibration with CRN support
# ===========================================================================

#' Calibrate sigma_theta to match a target r_MZ (Shenhar's model)
#'
#' Uses seed-per-iteration approach (stochastic bisection).
#' For deterministic bisection, use calibrate_sigma_theta_crn().
#'
#' @param target_r_mz Target MZ correlation
#' @param m_ex Extrinsic hazard rate for calibration sims
#' @param lo Lower bound for search
#' @param hi Upper bound for search
#' @param n_iter Number of bisection iterations
#' @param cutoff Minimum age for inclusion
#' @return Calibrated sigma_theta
calibrate_sigma_theta <- function(target_r_mz, m_ex,
                                  lo = 0.01, hi = 5.0,
                                  n_iter = PARAMS$N_CALIB_ITER,
                                  cutoff = PARAMS$CUTOFF_AGE) {
  for (i in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    result <- sim_twin_h2(mid, m_ex,
      sigma_gamma = 0, rho = 0,
      individual_ext = FALSE,
      rng_seed = PARAMS$MASTER_SEED + i * 100,
      cutoff = cutoff
    )
    if (result$r_mz < target_r_mz) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}

#' Calibrate sigma_theta to match a target h² at m_ex=0 (oracle calibration)
#'
#' @param target_h2 Target heritability
#' @param lo Lower bound
#' @param hi Upper bound
#' @param n_iter Number of bisection iterations
#' @return Calibrated sigma_theta
calibrate_oracle <- function(target_h2 = PARAMS$TARGET_H2,
                             lo = 0.5, hi = 4.0,
                             n_iter = PARAMS$N_CALIB_ITER) {
  for (i in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    test <- sim_twin_h2(mid, 0, sigma_gamma = 0,
                        rng_seed = PARAMS$MASTER_SEED + i * 50)
    if (test$h2 < target_h2) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}

#' Calibrate dispersion for SR/MGG models
#'
#' @param target_r_mz Target MZ correlation
#' @param model "mgg" or "sr"
#' @param m_ex Extrinsic hazard rate
#' @param lo Lower bound
#' @param hi Upper bound
#' @param n_iter Number of bisection iterations
#' @param n_pairs Number of twin pairs per evaluation
#' @return Calibrated dispersion (sigma_q or sigma_Xc)
calibrate_dispersion <- function(target_r_mz, model = "mgg",
                                 m_ex = PARAMS$M_EX_HIST,
                                 lo = 0.01, hi = NULL,
                                 n_iter = PARAMS$N_CALIB_ITER,
                                 n_pairs = PARAMS$N_PAIRS) {
  if (is.null(hi)) {
    hi <- if (model == "mgg") 0.80 else 12.0
  }

  sim_fn <- if (model == "mgg") sim_twin_h2_mgg else sim_twin_h2_sr

  for (i in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    if (model == "mgg") {
      res <- sim_fn(sigma_q = mid, m_ex = m_ex, n_pairs = n_pairs,
                    rng_seed = PARAMS$MASTER_SEED + 3000 + i * 100)
    } else {
      res <- sim_fn(sigma_Xc = mid, m_ex = m_ex, n_pairs = n_pairs,
                    rng_seed = PARAMS$MASTER_SEED + 3000 + i * 100)
    }
    if (res$r_mz < target_r_mz) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  (lo + hi) / 2
}

#' Calibrate with uncertainty quantification
#'
#' Runs calibration n_reps times with different master seeds,
#' reporting mean ± SD of all key outputs.
#'
#' @param target_r_mz Target MZ correlation (or target_h2 for oracle)
#' @param m_ex Extrinsic hazard rate
#' @param calibration_type "r_mz" or "oracle_h2"
#' @param n_reps Number of replications
#' @param base_seed Base seed (each rep uses base_seed + rep*1000)
#' @return List with means, SDs, and CI for sigma_fit and h2
calibrate_with_uncertainty <- function(target_r_mz = NULL,
                                       target_h2 = NULL,
                                       m_ex = PARAMS$M_EX_HIST,
                                       calibration_type = "r_mz",
                                       n_reps = 20,
                                       base_seed = PARAMS$MASTER_SEED) {
  results <- lapply(seq_len(n_reps), function(rep) {
    rep_seed <- base_seed + rep * 1000L

    if (calibration_type == "oracle_h2") {
      # Calibrate oracle: find sigma_theta giving target h2 at m_ex=0
      sigma_fit <- calibrate_oracle_crn(
        target_h2 = target_h2,
        crn_seed = rep_seed
      )
      check <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
                           rng_seed = rep_seed + 500L)
      list(sigma_fit = sigma_fit, h2 = check$h2,
           r_mz = check$r_mz, r_dz = check$r_dz)
    } else {
      # Calibrate sigma_theta to match r_MZ
      sigma_fit <- calibrate_sigma_theta_crn(
        target_r_mz = target_r_mz,
        m_ex = m_ex,
        crn_seed = rep_seed
      )
      check <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
                           rng_seed = rep_seed + 500L)
      list(sigma_fit = sigma_fit, h2 = check$h2,
           r_mz = check$r_mz, r_dz = check$r_dz)
    }
  })

  sigma_fits <- sapply(results, `[[`, "sigma_fit")
  h2s <- sapply(results, `[[`, "h2")

  list(
    sigma_fit_mean = mean(sigma_fits),
    sigma_fit_sd = sd(sigma_fits),
    h2_mean = mean(h2s),
    h2_sd = sd(h2s),
    h2_ci95 = quantile(h2s, c(0.025, 0.975)),
    n_reps = n_reps,
    all_sigma_fits = sigma_fits,
    all_h2s = h2s
  )
}

#' CRN-based sigma_theta calibration (deterministic bisection)
#'
#' Pre-generates all random draws, making the objective function
#' deterministic and monotone within a calibration run.
#'
#' @param target_r_mz Target MZ correlation
#' @param m_ex Extrinsic hazard rate
#' @param crn_seed Seed for pre-generating draws
#' @param lo Lower bound
#' @param hi Upper bound
#' @param n_iter Number of bisection iterations
#' @param n Number of twin pairs
#' @return Calibrated sigma_theta
calibrate_sigma_theta_crn <- function(target_r_mz, m_ex,
                                      crn_seed = PARAMS$MASTER_SEED,
                                      lo = 0.01, hi = 5.0,
                                      n_iter = PARAMS$N_CALIB_ITER,
                                      n = PARAMS$N_PAIRS) {
  set.seed(crn_seed)
  # Pre-generate all needed draws
  # For GM model: need normals for twin params + uniforms for survival inversion
  # MZ: 2 normals (z1, z2) + 2 uniforms (u_mz1, u_mz2)
  # DZ: 2 normals (z1, z2) + 2 uniforms (u_dz1, u_dz2)
  z_mz1 <- rnorm(n)
  z_mz2 <- rnorm(n)
  u_mz1 <- runif(n)
  u_mz2 <- runif(n)
  z_dz1 <- rnorm(n)
  z_dz2 <- rnorm(n)
  u_dz1 <- runif(n)
  u_dz2 <- runif(n)

  b <- PARAMS$B_GOMP
  t_max <- PARAMS$T_MAX
  cutoff <- PARAMS$CUTOFF_AGE
  mu_theta <- PARAMS$MU_THETA

  for (i in seq_len(n_iter)) {
    mid <- (lo + hi) / 2

    # MZ twins (r_g = 1): theta1 = theta2
    theta_mz <- mu_theta + mid * z_mz1
    L_mz1 <- sim_lifespan_gm(theta_mz, b, m_ex, t_max, u = u_mz1)
    L_mz2 <- sim_lifespan_gm(theta_mz, b, m_ex, t_max, u = u_mz2)

    # DZ twins (r_g = 0.5)
    g_dz1 <- z_dz1
    g_dz2 <- 0.5 * z_dz1 + sqrt(0.75) * z_dz2
    theta_dz1 <- mu_theta + mid * g_dz1
    theta_dz2 <- mu_theta + mid * g_dz2
    L_dz1 <- sim_lifespan_gm(theta_dz1, b, m_ex, t_max, u = u_dz1)
    L_dz2 <- sim_lifespan_gm(theta_dz2, b, m_ex, t_max, u = u_dz2)

    keep_mz <- (L_mz1 > cutoff) & (L_mz2 > cutoff)
    keep_dz <- (L_dz1 > cutoff) & (L_dz2 > cutoff)

    r_mz <- cor(L_mz1[keep_mz], L_mz2[keep_mz])

    if (r_mz < target_r_mz) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}

#' CRN-based oracle calibration
#'
#' @param target_h2 Target heritability at m_ex=0
#' @param crn_seed Seed for pre-generating draws
#' @param lo Lower bound
#' @param hi Upper bound
#' @param n_iter Number of bisection iterations
#' @param n Number of twin pairs
#' @return Calibrated sigma_theta
calibrate_oracle_crn <- function(target_h2 = PARAMS$TARGET_H2,
                                 crn_seed = PARAMS$MASTER_SEED,
                                 lo = 0.5, hi = 4.0,
                                 n_iter = PARAMS$N_CALIB_ITER,
                                 n = PARAMS$N_PAIRS) {
  set.seed(crn_seed)
  z_mz1 <- rnorm(n)
  z_mz2 <- rnorm(n)
  u_mz1 <- runif(n)
  u_mz2 <- runif(n)
  z_dz1 <- rnorm(n)
  z_dz2 <- rnorm(n)
  u_dz1 <- runif(n)
  u_dz2 <- runif(n)

  b <- PARAMS$B_GOMP
  t_max <- PARAMS$T_MAX
  cutoff <- PARAMS$CUTOFF_AGE
  mu_theta <- PARAMS$MU_THETA

  for (i in seq_len(n_iter)) {
    mid <- (lo + hi) / 2

    theta_mz <- mu_theta + mid * z_mz1
    L_mz1 <- sim_lifespan_gm(theta_mz, b, 0, t_max, u = u_mz1)
    L_mz2 <- sim_lifespan_gm(theta_mz, b, 0, t_max, u = u_mz2)

    g_dz1 <- z_dz1
    g_dz2 <- 0.5 * z_dz1 + sqrt(0.75) * z_dz2
    theta_dz1 <- mu_theta + mid * g_dz1
    theta_dz2 <- mu_theta + mid * g_dz2
    L_dz1 <- sim_lifespan_gm(theta_dz1, b, 0, t_max, u = u_dz1)
    L_dz2 <- sim_lifespan_gm(theta_dz2, b, 0, t_max, u = u_dz2)

    keep_mz <- (L_mz1 > cutoff) & (L_mz2 > cutoff)
    keep_dz <- (L_dz1 > cutoff) & (L_dz2 > cutoff)

    r_mz <- cor(L_mz1[keep_mz], L_mz2[keep_mz])
    r_dz <- cor(L_dz1[keep_dz], L_dz2[keep_dz])
    h2 <- 2 * (r_mz - r_dz)

    if (h2 < target_h2) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}

#' Calibrate sigma_theta under correctly specified two-component model
#'
#' Oracle fix arm: holds (sigma_gamma, rho) at their true DGP values
#' and calibrates only sigma_theta to match observed r_MZ.
#' If the bias is from omitting gamma, this should recover the true sigma_theta.
#'
#' @param target_r_mz Target MZ correlation (from misspecified Arm 2 data)
#' @param m_ex Extrinsic hazard rate
#' @param sigma_gamma True sigma_gamma from DGP
#' @param rho True rho from DGP
#' @param lo Lower bound for search
#' @param hi Upper bound for search
#' @param n_iter Number of bisection iterations
#' @param cutoff Minimum age for inclusion
#' @return Calibrated sigma_theta
calibrate_sigma_theta_fix <- function(target_r_mz, m_ex,
                                      sigma_gamma, rho,
                                      lo = 0.01, hi = 5.0,
                                      n_iter = PARAMS$N_CALIB_ITER,
                                      cutoff = PARAMS$CUTOFF_AGE) {
  for (i in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    result <- sim_twin_h2(mid, m_ex,
      sigma_gamma = sigma_gamma, rho = rho,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = PARAMS$MASTER_SEED + 90000L + i * 100L,
      cutoff = cutoff
    )
    if (result$r_mz < target_r_mz) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  (lo + hi) / 2
}

#' Joint (r_MZ, r_DZ) calibration using optim
#'
#' Minimizes weighted SSE over both MZ and DZ correlations.
#'
#' @param target_r_mz Target MZ correlation
#' @param target_r_dz Target DZ correlation
#' @param m_ex Extrinsic hazard rate
#' @param rng_seed RNG seed
#' @return Calibrated sigma_theta
calibrate_joint <- function(target_r_mz, target_r_dz, m_ex,
                            rng_seed = PARAMS$MASTER_SEED) {
  # Fisher-z transform for better optimization landscape
  z_target_mz <- atanh(target_r_mz)
  z_target_dz <- atanh(target_r_dz)

  clamp_cor <- function(r) max(min(r, 0.999999), -0.999999)

  optim(par = 1.5, fn = function(sigma) {
    res <- sim_twin_h2(sigma, m_ex, sigma_gamma = 0,
                       rng_seed = rng_seed)
    if (!is.finite(res$r_mz) || !is.finite(res$r_dz)) return(1e10)
    z_mz <- atanh(clamp_cor(res$r_mz))
    z_dz <- atanh(clamp_cor(res$r_dz))
    (z_mz - z_target_mz)^2 + (z_dz - z_target_dz)^2
  }, method = "Brent", lower = 0.1, upper = 5.0)$par
}

#' Multi-target calibration: fit sigma_theta to multiple summary statistics
#'
#' Minimizes weighted SSE over r_MZ, mean age, SD age, q25, q75.
#' With 1 free parameter and 5 targets the system is overdetermined,
#' so the optimizer finds a best-compromise sigma_theta.
#'
#' @param targets Named list: r_mz, mean_age, sd_age, q25, q75
#' @param m_ex Extrinsic hazard rate
#' @param rng_seed RNG seed for reproducibility
#' @param n Number of twin pairs per evaluation
#' @param weights Named vector of relative weights
#' @return List with sigma_fit, loss, fitted_targets, residuals, gof_rmse
calibrate_multi_target <- function(targets, m_ex,
                                    rng_seed = PARAMS$MASTER_SEED,
                                    n = PARAMS$N_PAIRS,
                                    weights = c(r_mz = 1, mean_age = 1,
                                                sd_age = 1, q25 = 1, q75 = 1)) {
  z_target_mz <- atanh(max(min(targets$r_mz, 0.999), -0.999))

  obs_vec <- c(
    r_mz = z_target_mz,
    mean_age = targets$mean_age,
    sd_age = targets$sd_age,
    q25 = targets$q25,
    q75 = targets$q75
  )
  # Standardize by target magnitude so all components are O(1)
  scales <- pmax(abs(obs_vec), 0.1)

  eval_fn <- function(sigma) {
    # Simulate MZ twin pairs with concordant-survivor conditioning
    # (matching the conditioning in run_multi_target_arm)
    set.seed(rng_seed + 99L)
    pars <- gen_twin_params_gm(n, r_g = 1.0, sigma, 0, 0)
    L1 <- sim_lifespan_gm(pars$theta1, PARAMS$B_GOMP, m_ex, PARAMS$T_MAX)
    L2 <- sim_lifespan_gm(pars$theta2, PARAMS$B_GOMP, m_ex, PARAMS$T_MAX)
    keep_pair <- (L1 > PARAMS$CUTOFF_AGE) & (L2 > PARAMS$CUTOFF_AGE)
    if (sum(keep_pair) < 100) return(list(loss = 1e10, sim_vec = rep(NA, 5)))

    r_mz_sim <- cor(L1[keep_pair], L2[keep_pair])
    if (!is.finite(r_mz_sim)) return(list(loss = 1e10, sim_vec = rep(NA, 5)))

    # Marginals from concordant pairs (same conditioning as observed targets)
    Lk <- c(L1[keep_pair], L2[keep_pair])

    sim_vec <- c(
      r_mz = atanh(max(min(r_mz_sim, 0.999), -0.999)),
      mean_age = mean(Lk),
      sd_age = sd(Lk),
      q25 = unname(quantile(Lk, 0.25)),
      q75 = unname(quantile(Lk, 0.75))
    )

    resid <- (sim_vec - obs_vec) / scales
    loss <- sum(weights[names(resid)] * resid^2)
    list(loss = loss, sim_vec = sim_vec)
  }

  opt <- optim(par = 1.5, fn = function(s) eval_fn(s)$loss,
               method = "Brent", lower = 0.1, upper = 5.0)

  # Final evaluation for diagnostics
  final <- eval_fn(opt$par)
  fitted <- list(
    r_mz = unname(tanh(final$sim_vec["r_mz"])),
    mean_age = unname(final$sim_vec["mean_age"]),
    sd_age = unname(final$sim_vec["sd_age"]),
    q25 = unname(final$sim_vec["q25"]),
    q75 = unname(final$sim_vec["q75"])
  )
  residuals <- mapply(function(f, t) f - t,
                       fitted, targets[names(fitted)])

  list(
    sigma_fit = opt$par,
    loss = opt$value,
    fitted_targets = fitted,
    residuals = residuals,
    gof_rmse = sqrt(mean(unlist(residuals)^2))
  )
}
