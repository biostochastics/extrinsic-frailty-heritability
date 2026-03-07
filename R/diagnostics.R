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
#' Shows that halving dt doesn't materially change Arm 2 bias.
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
