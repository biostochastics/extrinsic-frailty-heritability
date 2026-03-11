# ===========================================================================
# diagnostics.R — Variance decomposition, joint diagnostic, α/β estimation,
#                  MGG hazard curves, SR dt-sensitivity
# ===========================================================================

#' Variance decomposition diagnostic
#'
#' Tests: sigma_fit^2 ≈ sigma_true^2 + alpha*sigma_gamma^2
#'        + 2*beta*rho*sigma_theta*sigma_gamma
#'
#' @param sweep_df Data frame with sigma_fit, sigma_gamma, rho columns
#' @param sigma_theta_true True sigma_theta
#' @return List with predicted_var, actual_var, residuals
variance_decomposition <- function(sweep_df, sigma_theta_true) {
  sweep_df$predicted_var <- sigma_theta_true^2 + sweep_df$sigma_gamma^2 +
    2 * sweep_df$rho * sigma_theta_true * sweep_df$sigma_gamma
  sweep_df$actual_var <- sweep_df$sigma_fit^2
  sweep_df$residual <- sweep_df$actual_var - sweep_df$predicted_var
  sweep_df$residual_pct <- 100 * sweep_df$residual / sweep_df$predicted_var

  list(
    data = sweep_df,
    mar = mean(abs(sweep_df$residual)),
    mar_pct = mean(abs(sweep_df$residual_pct)),
    max_resid = max(abs(sweep_df$residual)),
    max_resid_pct = max(abs(sweep_df$residual_pct))
  )
}

#' Joint (r_MZ, r_DZ) calibration diagnostic (Vulnerability #1)
#'
#' After fitting sigma_theta to match r_MZ, report the predicted r_DZ
#' under the misspecified model and compare to true DGP r_DZ.
#'
#' @param sigma_theta_true True sigma_theta
#' @param params Parameter list
#' @return List with true vs predicted correlations and discrepancies
run_joint_rmz_rdz_diagnostic <- function(sigma_theta_true, params) {
  # 1. Generate data under true DGP (heritable gamma)
  true_dgp <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
    sigma_gamma = params$SIGMA_GAMMA_DEF,
    rho = params$RHO_DEF,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = params$MASTER_SEED + 42001L
  )

  # 2. Shenhar's calibration: fit sigma_theta to match r_MZ
  sigma_fit <- calibrate_sigma_theta(true_dgp$r_mz, params$M_EX_HIST)

  # 3. Predict r_DZ under misspecified model with fitted sigma
  misspec_pred <- sim_twin_h2(sigma_fit, params$M_EX_HIST,
    sigma_gamma = 0,
    rng_seed = params$MASTER_SEED + 42002L
  )

  # 4. Also do joint calibration via optim
  sigma_joint <- calibrate_joint(
    true_dgp$r_mz, true_dgp$r_dz, params$M_EX_HIST,
    rng_seed = params$MASTER_SEED + 42003L
  )
  joint_pred <- sim_twin_h2(sigma_joint, 0, sigma_gamma = 0,
    rng_seed = params$MASTER_SEED + 42004L
  )

  list(
    true_r_mz = true_dgp$r_mz,
    true_r_dz = true_dgp$r_dz,
    misspec_pred_r_mz = misspec_pred$r_mz,
    misspec_pred_r_dz = misspec_pred$r_dz,
    r_dz_discrepancy = misspec_pred$r_dz - true_dgp$r_dz,
    true_h2 = true_dgp$h2,
    misspec_h2 = misspec_pred$h2,
    sigma_fit_rmz_only = sigma_fit,
    sigma_fit_joint = sigma_joint,
    joint_h2 = joint_pred$h2,
    joint_bias = joint_pred$h2 - sim_twin_h2(sigma_theta_true, 0,
      sigma_gamma = 0, rng_seed = params$MASTER_SEED)$h2
  )
}

#' Empirical α, β estimation from sweep results (Task 11)
#'
#' Fits: sigma_fit^2 = sigma_true^2 + alpha*sigma_gamma^2
#'       + 2*beta*rho*sigma_true*sigma_gamma
#'
#' @param sweep_df Data frame from sweep_results
#' @param sigma_theta_true True sigma_theta
#' @return List with alpha, beta, R², residuals
estimate_alpha_beta <- function(sweep_df, sigma_theta_true) {
  y <- sweep_df$sigma_fit^2 - sigma_theta_true^2
  x1 <- sweep_df$sigma_gamma^2
  x2 <- 2 * sweep_df$rho * sigma_theta_true * sweep_df$sigma_gamma

  fit <- lm(y ~ 0 + x1 + x2)
  list(
    alpha = unname(coef(fit)[1]),
    beta = unname(coef(fit)[2]),
    r_squared = summary(fit)$r.squared,
    residuals_pct = 100 * residuals(fit) / (sigma_theta_true^2 + fitted(fit)),
    max_resid_pct = max(abs(100 * residuals(fit) / (sigma_theta_true^2 + fitted(fit))))
  )
}

#' MGG hazard curve comparison figure data (Task 8)
#'
#' Generates hazard curves for q ∈ {0.7, 1.0, 1.3} under both SM mappings.
#'
#' @param params Parameter list
#' @return Data frame with t, hazard, q, mapping columns
run_mgg_hazard_curves <- function(params) {
  q_vals <- c(0.7, 1.0, 1.3)
  t_grid <- seq(0, 110, by = 0.5)
  ec <- exp(params$MGG_C)

  curves_list <- lapply(c("compensatory", "paper"), function(mapping) {
    do.call(rbind, lapply(q_vals, function(q) {
      b <- q * params$MGG_B0
      if (mapping == "compensatory") {
        a <- params$MGG_A0^q
      } else {
        a <- params$MGG_A0^(1 / q)
      }
      ebt <- exp(pmin(b * t_grid, 700))
      D <- ec + ebt - 1
      mu_int <- a * ebt * ec / D
      data.frame(
        t = t_grid, hazard = mu_int,
        q = factor(q), mapping = mapping,
        stringsAsFactors = FALSE
      )
    }))
  })

  result <- do.call(rbind, curves_list)
  # SM convergence age
  attr(result, "t_star") <- -log(params$MGG_A0) / params$MGG_B0
  result
}

#' SR dt-sensitivity check (Task 12)
#'
#' Shows that halving dt doesn't materially change Misspecified bias.
#'
#' @param params Parameter list
#' @return Data frame with dt, arm2_bias columns
run_sr_dt_check <- function(params) {
  dt_values <- c(1 / 26, 1 / 52, 1 / 104)  # biweekly, weekly, half-weekly
  n_small <- min(params$N_PAIRS, 10000L)

  results <- lapply(dt_values, function(dt_val) {
    p <- params
    p$SR_DT <- dt_val
    p$N_PAIRS <- n_small
    p$N_CALIB_ITER <- 25L
    res <- run_model_analysis("sr", p)
    data.frame(
      dt = dt_val,
      dt_label = sprintf("1/%d", round(1 / dt_val)),
      oracle_h2 = res$oracle_h2,
      arm2_h2 = res$arm2_h2,
      arm2_bias_pp = round(100 * res$arm2_bias, 1),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}


#' MGG parameterization comparison (a=a0^q vs a=a0^{1/q})
#'
#' Side-by-side comparison showing that the paper's stated formula
#' (a=a0^{1/q}) breaks Strehler-Mildvan compensation and produces
#' unrealistic survival patterns and twin correlations.
#'
#' @param params Parameter list
#' @return List with hazard_df, survival_df, twin_stats, demo_stats
run_mgg_param_comparison <- function(params) {
  mappings <- c("compensatory", "paper")
  q_vals <- c(0.7, 1.0, 1.3)
  t_grid <- seq(0, 110, by = 0.5)
  ec <- exp(params$MGG_C)

  # --- Panel A data: Hazard curves ---
  hazard_list <- lapply(mappings, function(mapping) {
    do.call(rbind, lapply(q_vals, function(q) {
      b <- q * params$MGG_B0
      a <- if (mapping == "compensatory") params$MGG_A0^q else params$MGG_A0^(1 / q)
      ebt <- exp(pmin(b * t_grid, 700))
      D <- ec + ebt - 1
      mu_int <- a * ebt * ec / D
      data.frame(t = t_grid, hazard = mu_int,
                 q = factor(q), mapping = mapping,
                 stringsAsFactors = FALSE)
    }))
  })
  hazard_df <- do.call(rbind, hazard_list)
  t_star <- -log(params$MGG_A0) / params$MGG_B0

  # --- Panel B data: Cohort survival curves ---
  set.seed(params$MASTER_SEED + 55000L)
  n_cohort <- 20000L
  q_cohort <- rnorm(n_cohort, mean = params$MGG_MU_Q, sd = params$MGG_CV_Q)
  q_cohort <- pmax(q_cohort, 0.01)
  u_cohort <- runif(n_cohort)

  surv_list <- lapply(mappings, function(mapping) {
    lifespans <- sim_lifespan_mgg(q_cohort, m_ex = 0,
      a0 = params$MGG_A0, b0 = params$MGG_B0, c_param = params$MGG_C,
      u = u_cohort, sm_mapping = mapping)

    # Empirical survival function on a grid
    ages <- seq(0, 120, by = 1)
    surv <- sapply(ages, function(a) mean(lifespans > a))
    data.frame(age = ages, survival = surv, mapping = mapping,
               mean_age = mean(lifespans), median_age = median(lifespans),
               stringsAsFactors = FALSE)
  })
  survival_df <- do.call(rbind, surv_list)

  # --- Panel C data: Twin h² under each mapping ---
  seed_base <- params$MASTER_SEED + 56000L
  twin_list <- lapply(mappings, function(mapping) {
    res <- sim_twin_h2_mgg(
      sigma_q = params$MGG_CV_Q, mu_q = params$MGG_MU_Q,
      m_ex = 0, n_pairs = params$N_PAIRS,
      rng_seed = seed_base, sm_mapping = mapping)
    data.frame(
      mapping = mapping,
      r_mz = res$r_mz, r_dz = res$r_dz,
      h2 = res$h2, n_mz = res$n_mz, n_dz = res$n_dz,
      stringsAsFactors = FALSE)
  })
  twin_stats <- do.call(rbind, twin_list)

  # --- Demographic stats ---
  demo_stats <- data.frame(
    mapping = mappings,
    mean_age = c(
      survival_df$mean_age[survival_df$mapping == "compensatory"][1],
      survival_df$mean_age[survival_df$mapping == "paper"][1]),
    median_age = c(
      survival_df$median_age[survival_df$mapping == "compensatory"][1],
      survival_df$median_age[survival_df$mapping == "paper"][1]),
    stringsAsFactors = FALSE)

  list(
    hazard_df = hazard_df,
    survival_df = survival_df,
    twin_stats = twin_stats,
    demo_stats = demo_stats,
    t_star = t_star
  )
}


# ===========================================================================
# Bivariate survival model check
# ===========================================================================

#' Simulation-based model check via bivariate survival dependence
#'
#' Compares age-specific dependence structure under three DGPs:
#' (1) True DGP with heritable extrinsic frailty,
#' (2) Misspecified model (theta-only, calibrated to match r_MZ),
#' (3) Recovery: two-component refit (correctly specified, positive control).
#'
#' Uses held-out seeds (200000+ block, well separated from calibration
#' seeds 1–100000) for evaluation to prevent information leakage.
#' Returns 1D dependence-ratio curves R(t) and 2D bivariate surface
#' residuals for MZ twins.
#'
#' Expected R(t) signature: the misspecified model (inflated sigma_theta,
#' no gamma) should show excess late-life dependence and deficit in
#' early-life dependence relative to the true DGP, because it forces
#' all concordance through intrinsic aging rather than shared extrinsic
#' frailty. The recovery model should track the true DGP closely.
#'
#' @param sigma_theta_true True sigma_theta (oracle calibrated)
#' @param sigma_theta_fit Calibrated sigma_theta from misspecified model
#' @param sigma_theta_fix Calibrated sigma_theta from correctly specified model
#' @param params Parameter list
#' @return List with dep_curves (data frame), biv_diff (data frame),
#'   IAE scalars, peak deviation age, sample sizes
run_bivariate_survival_check <- function(sigma_theta_true, sigma_theta_fit,
                                          sigma_theta_fix, params) {
  n <- params$N_PAIRS
  b <- params$B_GOMP
  t_max <- params$T_MAX
  cutoff <- params$CUTOFF_AGE
  m_ex <- params$M_EX_HIST
  sg <- params$SIGMA_GAMMA_DEF
  rho <- params$RHO_DEF

  # Guard: if recovery model calibration failed, return NA structure
  if (!is.finite(sigma_theta_fix)) {
    warning("sigma_theta_fix is NA; bivariate check skipped")
    na_df <- data.frame(age = numeric(), R = numeric(),
                        zygosity = character(), model = character(),
                        stringsAsFactors = FALSE)
    na_biv <- data.frame(t1 = numeric(), t2 = numeric(),
                         delta_R = numeric(), model = character(),
                         stringsAsFactors = FALSE)
    return(list(
      dep_curves = na_df, biv_diff = na_biv, age_2d = numeric(),
      iae_misspec_mz = NA_real_, iae_misspec_dz = NA_real_,
      iae_fix_mz = NA_real_, iae_fix_dz = NA_real_,
      peak_dev_age_mz = NA_real_, peak_dev_age_dz = NA_real_,
      n_mz_true = 0L, n_dz_true = 0L, n_pairs = n
    ))
  }

  # Held-out seed block (200000+, separated from calibration 1–100000)
  seed_eval <- params$MASTER_SEED + 200000L

  # 1D age grid (starts above cutoff; conditional on both twins surviving)
  age_grid <- seq(cutoff + 5, 100, by = 1)

  # Minimum joint count for stable R(t) estimates
  min_joint_count <- 30L

  # Helper: dependence ratio R(t) = log(P(T1>t, T2>t) / P(T>t)^2)
  compute_dep_ratio <- function(L1, L2, t_grid, keep) {
    l1 <- L1[keep]
    l2 <- L2[keep]
    n_eff <- length(l1)
    if (n_eff == 0) return(rep(NA_real_, length(t_grid)))
    sapply(t_grid, function(t) {
      n_joint <- sum(l1 > t & l2 > t)
      marg <- mean(c(l1, l2) > t)
      # Mask: require min count AND marginal support > 3%
      if (n_joint < min_joint_count || marg < 0.03) return(NA_real_)
      joint <- n_joint / n_eff
      log(joint / marg^2)
    })
  }

  # Helper: generate MZ + DZ lifespans under specified DGP
  gen_lifespans <- function(sigma_theta, sigma_gamma, rho_val,
                            gamma_heritable, individual_ext,
                            m_ex_val, seed) {
    set.seed(seed)

    # MZ
    pars_mz <- gen_twin_params_gm(n, r_g = 1.0, sigma_theta,
                                   sigma_gamma, rho_val, gamma_heritable)
    if (individual_ext && sigma_gamma > 0) {
      shift <- sigma_gamma^2 / 2
      c_mz1 <- m_ex_val * exp(pars_mz$gamma1 - shift)
      c_mz2 <- m_ex_val * exp(pars_mz$gamma2 - shift)
    } else {
      c_mz1 <- c_mz2 <- m_ex_val
    }
    L_mz1 <- sim_lifespan_gm(pars_mz$theta1, b, c_mz1, t_max)
    L_mz2 <- sim_lifespan_gm(pars_mz$theta2, b, c_mz2, t_max)

    # DZ
    pars_dz <- gen_twin_params_gm(n, r_g = 0.5, sigma_theta,
                                   sigma_gamma, rho_val, gamma_heritable)
    if (individual_ext && sigma_gamma > 0) {
      shift <- sigma_gamma^2 / 2
      c_dz1 <- m_ex_val * exp(pars_dz$gamma1 - shift)
      c_dz2 <- m_ex_val * exp(pars_dz$gamma2 - shift)
    } else {
      c_dz1 <- c_dz2 <- m_ex_val
    }
    L_dz1 <- sim_lifespan_gm(pars_dz$theta1, b, c_dz1, t_max)
    L_dz2 <- sim_lifespan_gm(pars_dz$theta2, b, c_dz2, t_max)

    keep_mz <- (L_mz1 > cutoff) & (L_mz2 > cutoff)
    keep_dz <- (L_dz1 > cutoff) & (L_dz2 > cutoff)

    list(L_mz1 = L_mz1, L_mz2 = L_mz2, keep_mz = keep_mz,
         L_dz1 = L_dz1, L_dz2 = L_dz2, keep_dz = keep_dz)
  }

  # --- Generate lifespans under three DGPs (held-out seeds) ---
  true <- gen_lifespans(sigma_theta_true, sg, rho,
                        gamma_heritable = TRUE, individual_ext = TRUE,
                        m_ex, seed_eval)

  misspec <- gen_lifespans(sigma_theta_fit, 0, 0,
                           gamma_heritable = FALSE, individual_ext = FALSE,
                           m_ex, seed_eval + 1000L)

  fix <- gen_lifespans(sigma_theta_fix, sg, rho,
                       gamma_heritable = TRUE, individual_ext = TRUE,
                       m_ex, seed_eval + 2000L)

  # --- 1D dependence ratio curves R(t) ---
  R_true_mz    <- compute_dep_ratio(true$L_mz1,    true$L_mz2,    age_grid, true$keep_mz)
  R_true_dz    <- compute_dep_ratio(true$L_dz1,    true$L_dz2,    age_grid, true$keep_dz)
  R_misspec_mz <- compute_dep_ratio(misspec$L_mz1, misspec$L_mz2, age_grid, misspec$keep_mz)
  R_misspec_dz <- compute_dep_ratio(misspec$L_dz1, misspec$L_dz2, age_grid, misspec$keep_dz)
  R_fix_mz     <- compute_dep_ratio(fix$L_mz1,     fix$L_mz2,     age_grid, fix$keep_mz)
  R_fix_dz     <- compute_dep_ratio(fix$L_dz1,     fix$L_dz2,     age_grid, fix$keep_dz)

  dep_curves <- rbind(
    data.frame(age = age_grid, R = R_true_mz,    zygosity = "MZ", model = "True DGP",
               stringsAsFactors = FALSE),
    data.frame(age = age_grid, R = R_true_dz,    zygosity = "DZ", model = "True DGP",
               stringsAsFactors = FALSE),
    data.frame(age = age_grid, R = R_misspec_mz, zygosity = "MZ", model = "Misspecified fit",
               stringsAsFactors = FALSE),
    data.frame(age = age_grid, R = R_misspec_dz, zygosity = "DZ", model = "Misspecified fit",
               stringsAsFactors = FALSE),
    data.frame(age = age_grid, R = R_fix_mz,     zygosity = "MZ", model = "Recovery: two-component refit",
               stringsAsFactors = FALSE),
    data.frame(age = age_grid, R = R_fix_dz,     zygosity = "DZ", model = "Recovery: two-component refit",
               stringsAsFactors = FALSE)
  )

  # --- IAE (integrated absolute error) scalars ---
  .iae <- function(R_ref, R_test) {
    ok <- !is.na(R_ref) & !is.na(R_test)
    if (sum(ok) == 0) return(NA_real_)
    mean(abs(R_test[ok] - R_ref[ok]))
  }
  iae_misspec_mz <- .iae(R_true_mz, R_misspec_mz)
  iae_misspec_dz <- .iae(R_true_dz, R_misspec_dz)
  iae_fix_mz     <- .iae(R_true_mz, R_fix_mz)
  iae_fix_dz     <- .iae(R_true_dz, R_fix_dz)

  # Peak deviation age (where |R_misspec - R_true| is maximized)
  .peak_age <- function(R_ref, R_test, t_grid) {
    ok <- !is.na(R_ref) & !is.na(R_test)
    if (sum(ok) == 0) return(NA_real_)
    diffs <- abs(R_test[ok] - R_ref[ok])
    t_grid[ok][which.max(diffs)]
  }
  peak_dev_age_mz <- .peak_age(R_true_mz, R_misspec_mz, age_grid)
  peak_dev_age_dz <- .peak_age(R_true_dz, R_misspec_dz, age_grid)

  # --- 2D bivariate surface residuals (MZ only, for appendix) ---
  # Coarser grid than 1D (computational cost vs resolution trade-off)
  age_2d <- seq(cutoff + 5, 95, by = 2)

  compute_biv_R <- function(L1, L2, keep, t_grid) {
    l1 <- L1[keep]
    l2 <- L2[keep]
    n_eff <- length(l1)
    ng <- length(t_grid)
    R_mat <- matrix(NA_real_, nrow = ng, ncol = ng)
    if (n_eff == 0) return(R_mat)

    # Precompute marginals once per grid point
    marg1 <- sapply(t_grid, function(t) mean(l1 > t))
    marg2 <- sapply(t_grid, function(t) mean(l2 > t))

    for (i in seq_len(ng)) {
      if (marg1[i] < 0.03) next
      s1 <- l1 > t_grid[i]
      for (j in seq_len(ng)) {
        if (marg2[j] < 0.03) next
        n_joint <- sum(s1 & l2 > t_grid[j])
        if (n_joint < min_joint_count) next
        joint <- n_joint / n_eff
        R_mat[i, j] <- log(joint / (marg1[i] * marg2[j]))
      }
    }
    R_mat
  }

  biv_true_mz    <- compute_biv_R(true$L_mz1,    true$L_mz2,    true$keep_mz,    age_2d)
  biv_misspec_mz <- compute_biv_R(misspec$L_mz1, misspec$L_mz2, misspec$keep_mz, age_2d)
  biv_fix_mz     <- compute_biv_R(fix$L_mz1,     fix$L_mz2,     fix$keep_mz,     age_2d)

  delta_R_misspec <- biv_misspec_mz - biv_true_mz
  delta_R_fix     <- biv_fix_mz - biv_true_mz

  biv_to_df <- function(mat, t_grid, label) {
    idx <- which(!is.na(mat), arr.ind = TRUE)
    if (nrow(idx) == 0) return(data.frame(t1 = numeric(), t2 = numeric(),
                                           delta_R = numeric(), model = character(),
                                           stringsAsFactors = FALSE))
    data.frame(
      t1 = t_grid[idx[, 1]],
      t2 = t_grid[idx[, 2]],
      delta_R = mat[idx],
      model = label,
      stringsAsFactors = FALSE
    )
  }

  biv_diff <- rbind(
    biv_to_df(delta_R_misspec, age_2d, "Misspecified"),
    biv_to_df(delta_R_fix,     age_2d, "Recovery: two-component refit")
  )

  list(
    dep_curves = dep_curves,
    biv_diff = biv_diff,
    age_2d = age_2d,
    iae_misspec_mz = iae_misspec_mz,
    iae_misspec_dz = iae_misspec_dz,
    iae_fix_mz = iae_fix_mz,
    iae_fix_dz = iae_fix_dz,
    peak_dev_age_mz = peak_dev_age_mz,
    peak_dev_age_dz = peak_dev_age_dz,
    n_mz_true = sum(true$keep_mz),
    n_dz_true = sum(true$keep_dz),
    n_pairs = n
  )
}
