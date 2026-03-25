# ===========================================================================
# variance_components.R — Latent variance component inflation reporting
# ===========================================================================
# Primary estimand: Δσ_θ = σ_θ_fit - σ_θ_true and percent inflation.
# Reports fitted latent variance components rather than Falconer h².
# ===========================================================================

#' Compute σ_θ inflation metrics
#'
#' @param sigma_fit Calibrated (possibly inflated) sigma_theta
#' @param sigma_true True sigma_theta from oracle calibration
#' @return Named list with delta, pct, variance ratio
compute_sigma_inflation <- function(sigma_fit, sigma_true) {
  delta <- sigma_fit - sigma_true
  pct <- 100 * delta / sigma_true
  list(
    sigma_true = sigma_true,
    sigma_fit = sigma_fit,
    delta_sigma = delta,
    sigma_infl_pct = pct,
    var_ratio = sigma_fit^2 / sigma_true^2,
    var_infl_pct = 100 * (sigma_fit^2 / sigma_true^2 - 1)
  )
}

#' Build variance components table across all arms and models
#'
#' @param sigma_theta_true True σ_θ (oracle calibrated)
#' @param sigma_theta_arm1 Calibrated σ_θ for Arm1 (null)
#' @param sigma_theta_fit Calibrated σ_θ for Arm2 (misspecified)
#' @param sigma_theta_fix Recovery σ_θ (two-component refit)
#' @param oracle Oracle sim result (list with h2)
#' @param arm1_corr Arm1 corrected sim result
#' @param arm2_corr Arm2 corrected sim result
#' @param oracle_fix Recovery result (list with fix_h2)
#' @param mgg_validation MGG model analysis result (optional)
#' @param sr_validation SR model analysis result (optional)
#' @return Data frame with one row per condition, variance component columns
build_variance_components_table <- function(sigma_theta_true,
                                             sigma_theta_arm1,
                                             sigma_theta_fit,
                                             sigma_theta_fix,
                                             oracle, arm1_corr, arm2_corr,
                                             oracle_fix,
                                             mgg_validation = NULL,
                                             sr_validation = NULL) {
  # Helper to build one row
  make_row <- function(model, condition, s_true, s_fit, h2_falc) {
    data.frame(
      model = model, condition = condition,
      sigma_true = s_true, sigma_fit = s_fit,
      delta_sigma = s_fit - s_true,
      sigma_infl_pct = 100 * (s_fit / s_true - 1),
      var_infl_pct = 100 * (s_fit^2 / s_true^2 - 1),
      h2_falconer = h2_falc,
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    make_row("GM", "Oracle",
             sigma_theta_true, sigma_theta_true, oracle$h2),
    make_row("GM", "Baseline (correctly specified)",
             sigma_theta_true, sigma_theta_arm1, arm1_corr$h2),
    make_row("GM", "Misspecified (omitted extrinsic)",
             sigma_theta_true, sigma_theta_fit, arm2_corr$h2)
  )

  if (is.finite(sigma_theta_fix)) {
    rows[[length(rows) + 1]] <- make_row(
      "GM", "Recovery (two-component refit)",
      sigma_theta_true, sigma_theta_fix, oracle_fix$fix_h2)
  }

  # MGG
  if (!is.null(mgg_validation)) {
    td <- mgg_validation$true_disp
    rows[[length(rows) + 1]] <- make_row(
      "MGG", "Oracle", td, td, mgg_validation$oracle_h2)
    rows[[length(rows) + 1]] <- make_row(
      "MGG", "Baseline (correctly specified)",
      td, mgg_validation$arm1_disp, mgg_validation$arm1_h2)
    rows[[length(rows) + 1]] <- make_row(
      "MGG", "Misspecified (omitted extrinsic)",
      td, mgg_validation$arm2_disp, mgg_validation$arm2_h2)
  }

  # SR
  if (!is.null(sr_validation)) {
    td <- sr_validation$true_disp
    rows[[length(rows) + 1]] <- make_row(
      "SR", "Oracle", td, td, sr_validation$oracle_h2)
    rows[[length(rows) + 1]] <- make_row(
      "SR", "Baseline (correctly specified)",
      td, sr_validation$arm1_disp, sr_validation$arm1_h2)
    rows[[length(rows) + 1]] <- make_row(
      "SR", "Misspecified (omitted extrinsic)",
      td, sr_validation$arm2_disp, sr_validation$arm2_h2)
  }

  do.call(rbind, rows)
}

#' Compute variance absorption decomposition across sweep
#'
#' For each sweep cell: σ_fit² = σ_true² + absorbed_variance
#' Decomposes absorbed variance into σ_γ² contribution and cross-term.
#'
#' @param sweep_results Data frame from sweep_results target
#' @param sigma_theta_true True σ_θ
#' @return Data frame with decomposition columns added
compute_sweep_variance_decomposition <- function(sweep_results, sigma_theta_true) {
  df <- sweep_results
  df$sigma_true_sq <- sigma_theta_true^2
  df$sigma_fit_sq <- df$sigma_fit^2
  df$absorbed_var <- df$sigma_fit_sq - df$sigma_true_sq
  df$gamma_var_term <- df$sigma_gamma^2
  df$cross_term <- 2 * df$rho * sigma_theta_true * df$sigma_gamma
  df$predicted_absorbed <- df$gamma_var_term + df$cross_term
  df$residual_var <- df$absorbed_var - df$predicted_absorbed
  df$residual_pct <- 100 * df$residual_var / df$sigma_true_sq
  df
}

#' Compute replicated variance component inflation with uncertainty
#'
#' Runs multiple seeds of observe→calibrate→extrapolate and reports
#' σ_θ inflation statistics with Monte Carlo uncertainty.
#'
#' @param sigma_theta_true True σ_θ
#' @param params Parameter list
#' @param n_seeds Number of MC replication seeds
#' @return List with per_seed data frame and summary statistics
run_variance_component_uncertainty <- function(sigma_theta_true, params,
                                                n_seeds = 20L) {
  seed_base <- params$MASTER_SEED + 300000L
  results <- lapply(seq_len(n_seeds), function(s) {
    seed <- seed_base + (s - 1L) * 1000L

    # Oracle
    oracle_s <- sim_twin_h2(sigma_theta_true, 0,
      sigma_gamma = 0, rng_seed = seed)

    # Arm1: correctly specified
    obs_null <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = 0, rng_seed = seed + 100L)
    sigma_arm1 <- calibrate_sigma_theta(obs_null$r_mz, params$M_EX_HIST)
    arm1_s <- sim_twin_h2(sigma_arm1, 0, sigma_gamma = 0,
      rng_seed = seed + 200L)

    # Arm2: misspecified
    true_obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = seed + 300L)
    sigma_fit <- calibrate_sigma_theta(true_obs$r_mz, params$M_EX_HIST)
    arm2_s <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
      rng_seed = seed + 400L)

    data.frame(
      seed = s,
      sigma_arm1 = sigma_arm1,
      sigma_fit = sigma_fit,
      delta_sigma_arm1 = sigma_arm1 - sigma_theta_true,
      delta_sigma_arm2 = sigma_fit - sigma_theta_true,
      sigma_infl_arm1_pct = 100 * (sigma_arm1 / sigma_theta_true - 1),
      sigma_infl_arm2_pct = 100 * (sigma_fit / sigma_theta_true - 1),
      var_infl_arm2_pct = 100 * (sigma_fit^2 / sigma_theta_true^2 - 1),
      h2_oracle = oracle_s$h2,
      h2_arm1 = arm1_s$h2,
      h2_arm2 = arm2_s$h2,
      bias_pp_arm1 = 100 * (arm1_s$h2 - oracle_s$h2),
      bias_pp_arm2 = 100 * (arm2_s$h2 - oracle_s$h2),
      stringsAsFactors = FALSE
    )
  })

  per_seed <- do.call(rbind, results)

  # Summary statistics
  t_crit <- qt(0.975, df = n_seeds - 1)
  make_summary <- function(x) {
    m <- mean(x); se <- sd(x) / sqrt(length(x))
    data.frame(mean = m, se = se,
               lo95 = m - t_crit * se, hi95 = m + t_crit * se)
  }

  metrics <- c("sigma_infl_arm2_pct", "var_infl_arm2_pct",
               "delta_sigma_arm2", "bias_pp_arm2",
               "sigma_infl_arm1_pct", "bias_pp_arm1")
  summary_rows <- lapply(metrics, function(m) {
    s <- make_summary(per_seed[[m]])
    s$metric <- m
    s
  })
  summary <- do.call(rbind, summary_rows)

  list(per_seed = per_seed, summary = summary, n_seeds = n_seeds)
}
