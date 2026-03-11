# ===========================================================================
# export.R — Results export pipeline (CSV/JSON for Typst consumption)
# ===========================================================================
# All manuscript numbers flow through this module.
# Typst reads tables/*.csv and tables/scalars.json — no hardcoded values.
# ===========================================================================

#' Write all result tables to CSV files in tables/
#'
#' @param summary_table Data frame with main arms results
#' @param model_table Data frame with SR/MGG model validation
#' @param controls_table Data frame with additional controls
#' @param sweep_table Data frame with sensitivity sweep
#' @param anchored_table Data frame with anchored sweep
#' @param scalars Named list of scalar results
#' @param outdir Directory for output files
write_all_tables <- function(summary_table = NULL,
                             model_table = NULL,
                             controls_table = NULL,
                             sweep_table = NULL,
                             anchored_table = NULL,
                             model_controls = NULL,
                             scalars = NULL,
                             outdir = "tables") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  if (!is.null(summary_table)) {
    write.csv(summary_table, file.path(outdir, "summary.csv"),
              row.names = FALSE)
  }
  if (!is.null(model_table)) {
    write.csv(model_table, file.path(outdir, "models.csv"),
              row.names = FALSE)
  }
  if (!is.null(controls_table)) {
    write.csv(controls_table, file.path(outdir, "controls.csv"),
              row.names = FALSE)
  }
  if (!is.null(sweep_table)) {
    write.csv(sweep_table, file.path(outdir, "sweep.csv"),
              row.names = FALSE)
  }
  if (!is.null(anchored_table)) {
    write.csv(anchored_table, file.path(outdir, "anchored.csv"),
              row.names = FALSE)
  }
  if (!is.null(model_controls)) {
    write.csv(model_controls, file.path(outdir, "model_controls.csv"),
              row.names = FALSE)
  }
  if (!is.null(scalars)) {
    write_scalars(scalars, outdir)
  }

  invisible(outdir)
}

#' Write scalar results to JSON for Typst consumption
#'
#' @param scalars Named list of scalar values
#' @param outdir Output directory
write_scalars <- function(scalars, outdir = "tables") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  json_path <- file.path(outdir, "scalars.json")
  # Write JSON manually to avoid jsonlite dependency
  lines <- "{\n"
  nms <- names(scalars)
  for (i in seq_along(scalars)) {
    val <- scalars[[i]]
    if (is.character(val)) {
      val_str <- sprintf('"%s"', gsub('"', '\\\\"', val))
    } else if (is.logical(val)) {
      val_str <- if (val) "true" else "false"
    } else if (is.numeric(val) && !is.finite(val)) {
      val_str <- "null"
    } else {
      val_str <- format(val, digits = 6, scientific = FALSE)
    }
    comma <- if (i < length(scalars)) "," else ""
    lines <- paste0(lines, sprintf('  "%s": %s%s\n', nms[i], val_str, comma))
  }
  lines <- paste0(lines, "}\n")
  writeLines(lines, json_path)
  invisible(json_path)
}

#' Build the main summary table from simulation results
#'
#' @param oracle Oracle simulation result (list)
#' @param arm1 Baseline result
#' @param arm2 Misspecified result
#' @param ctrl_a Control: nonfamilial result
#' @param ctrl_b Control: vanishing result
#' @param arm3 Comparator (Hamilton) result
#' @param sigma_theta_true True sigma_theta
#' @param sigma_theta_arm1 Calibrated sigma for Baseline
#' @param sigma_theta_arm2 Calibrated sigma for Misspecified
#' @param sigma_theta_ctrl_a Calibrated sigma for Control: nonfamilial
#' @param sigma_theta_ctrl_b Calibrated sigma for Control: vanishing
#' @param oracle_fix Recovery (two-component refit) result
#' @param pleiotropy_isolation Check: pleiotropy only result
#' @param irrelevant_trait Control: irrelevant trait result
#' @param robustness_cutoff0 Check: no survivor cutoff result
#' @return Data frame
build_summary_table <- function(oracle, arm1, arm2, ctrl_a, ctrl_b, arm3,
                                sigma_theta_true,
                                sigma_theta_arm1, sigma_theta_arm2,
                                sigma_theta_ctrl_a, sigma_theta_ctrl_b,
                                oracle_fix = NULL,
                                pleiotropy_isolation = NULL,
                                irrelevant_trait = NULL,
                                robustness_cutoff0 = NULL) {
  conditions <- c(
    "Oracle",
    "Baseline: correctly specified",
    "Misspecified: omitted familial extrinsic",
    "Control: nonfamilial extrinsic",
    "Control: vanishing extrinsic"
  )
  sigma_vals <- c(
    sigma_theta_true, sigma_theta_arm1, sigma_theta_arm2,
    sigma_theta_ctrl_a, sigma_theta_ctrl_b
  )
  h2_vals <- c(oracle$h2, arm1$h2, arm2$h2, ctrl_a$h2, ctrl_b$h2)
  bias_vals <- c(0,
                 arm1$h2 - oracle$h2,
                 arm2$h2 - oracle$h2,
                 ctrl_a$h2 - oracle$h2,
                 ctrl_b$h2 - oracle$h2)

  # Check: pleiotropy only
  if (!is.null(pleiotropy_isolation)) {
    conditions <- c(conditions, "Check: pleiotropy only")
    sigma_vals <- c(sigma_vals, pleiotropy_isolation$sigma_fit)
    h2_vals <- c(h2_vals, pleiotropy_isolation$h2_corrected)
    bias_vals <- c(bias_vals, pleiotropy_isolation$bias)
  }

  # Control: irrelevant trait
  if (!is.null(irrelevant_trait)) {
    conditions <- c(conditions, "Control: irrelevant trait")
    sigma_vals <- c(sigma_vals, irrelevant_trait$sigma_fit)
    h2_vals <- c(h2_vals, irrelevant_trait$h2_corrected)
    bias_vals <- c(bias_vals, irrelevant_trait$bias)
  }

  # Recovery — bias relative to the GLOBAL oracle (h2_vals[1]), not local
  if (!is.null(oracle_fix)) {
    conditions <- c(conditions, "Recovery: two-component refit")
    sigma_vals <- c(sigma_vals, oracle_fix$sigma_theta_fix)
    h2_vals <- c(h2_vals, oracle_fix$fix_h2)
    bias_vals <- c(bias_vals, oracle_fix$fix_h2 - h2_vals[1])
  }

  # Check: no survivor cutoff
  if (!is.null(robustness_cutoff0)) {
    conditions <- c(conditions, "Check: no survivor cutoff")
    sigma_vals <- c(sigma_vals, robustness_cutoff0$sigma_fit)
    h2_vals <- c(h2_vals, robustness_cutoff0$corr_h2)
    bias_vals <- c(bias_vals, robustness_cutoff0$bias_pp / 100)
  }

  # Comparator always last
  conditions <- c(conditions, "Comparator: Hamilton")
  sigma_vals <- c(sigma_vals, NA)
  h2_vals <- c(h2_vals, arm3$h2)
  bias_vals <- c(bias_vals, arm3$h2 - oracle$h2)

  df <- data.frame(
    Condition = conditions,
    sigma_theta = sigma_vals,
    h2 = h2_vals,
    bias = bias_vals,
    stringsAsFactors = FALSE
  )

  # Canonical presentation order (matches Results prose)
  canonical <- c(
    "Oracle",
    "Baseline: correctly specified",
    "Misspecified: omitted familial extrinsic",
    "Recovery: two-component refit",
    "Check: no survivor cutoff",
    "Check: pleiotropy only",
    "Control: nonfamilial extrinsic",
    "Control: irrelevant trait",
    "Control: vanishing extrinsic",
    "Comparator: Hamilton"
  )
  ord <- match(df$Condition, canonical)
  df <- df[order(ord), ]
  rownames(df) <- NULL
  df
}

#' Build the model validation table
#'
#' @param oracle_gm GM oracle result
#' @param arm1_gm GM Baseline corrected result
#' @param arm2_gm GM Misspecified corrected result
#' @param mgg_result MGG model analysis result (list with oracle_h2, arm1_h2, etc.)
#' @param sr_result SR model analysis result
#' @return Data frame
build_model_table <- function(oracle_gm, arm1_gm, arm2_gm,
                              mgg_result, sr_result) {
  gm_arm1_bias <- arm1_gm$h2 - oracle_gm$h2
  gm_arm2_bias <- arm2_gm$h2 - oracle_gm$h2
  data.frame(
    Model = c("Gompertz-Makeham", "MGG (Shenhar)", "SR (Shenhar)"),
    Oracle_h2 = c(oracle_gm$h2, mgg_result$oracle_h2, sr_result$oracle_h2),
    Baseline_h2 = c(arm1_gm$h2, mgg_result$arm1_h2, sr_result$arm1_h2),
    Baseline_bias_pp = c(round(100 * gm_arm1_bias, 1),
                     round(100 * mgg_result$arm1_bias, 1),
                     round(100 * sr_result$arm1_bias, 1)),
    Misspec_h2 = c(arm2_gm$h2, mgg_result$arm2_h2, sr_result$arm2_h2),
    Misspec_bias_pp = c(round(100 * gm_arm2_bias, 1),
                     round(100 * mgg_result$arm2_bias, 1),
                     round(100 * sr_result$arm2_bias, 1)),
    stringsAsFactors = FALSE
  )
}

#' Collect ALL scalars needed by the manuscript into a named list
#'
#' Every number in the Typst manuscript should come from this list.
#' The Typst document loads tables/scalars.json and references these keys.
#'
#' @return Named list
collect_scalars <- function(oracle, arm1, arm2, ctrl_a, ctrl_b, arm3,
                            sigma_theta_true, sigma_theta_arm1,
                            sigma_theta_fit, sigma_theta_ctrl_a,
                            sigma_theta_ctrl_b,
                            var_decomp, sweep_results, anchored_results,
                            robustness_cutoff0, dose_response,
                            pleiotropy_isolation, irrelevant_trait,
                            mgg_validation, sr_validation,
                            empirical_alpha_beta,
                            joint_diagnostic,
                            sigma_gamma_bridge, mex_split,
                            oracle_fix,
                            multi_target_arm = NULL,
                            mgg_param_comparison = NULL,
                            alt_ext_forms = NULL,
                            extended_sigma_gamma = NULL,
                            mc_uncertainty = NULL,
                            mc_uncertainty_multimodel = NULL,
                            mex_split_hires = NULL,
                            bivariate_check = NULL,
                            mz_concordance = NULL,
                            mz_concordance_decoupled = NULL,
                            bridge_uncertainty = NULL,
                            main_arms_replicated = NULL) {

  arm2_bias <- arm2$h2 - oracle$h2
  arm1_bias <- arm1$h2 - oracle$h2

  # Sweep extremes (rho=0.8, sigma_gamma=0.55)
  sweep_max <- sweep_results[
    abs(sweep_results$rho - 0.8) < 0.01 &
    abs(sweep_results$sigma_gamma - 0.55) < 0.01, ]
  sweep_rho0_sg02 <- sweep_results[
    abs(sweep_results$rho) < 0.01 &
    abs(sweep_results$sigma_gamma - 0.2) < 0.01, ]

  # Anchored stats
  anch_biases <- anchored_results$bias
  anch_point <- anchored_results[
    abs(anchored_results$rho - 0.35) < 0.01 &
    abs(anchored_results$sigma_gamma - 0.65) < 0.01, ]

  # Dose-response at historical and 3x
  dr_hist <- dose_response[dose_response$m_ex == 0.004, ]
  dr_3x <- dose_response[dose_response$m_ex == 0.012, ]

  c(list(
    # --- Oracle ---
    oracle_h2 = round(oracle$h2, 3),
    oracle_r_mz = round(oracle$r_mz, 3),
    oracle_r_dz = round(oracle$r_dz, 3),

    # --- Baseline ---
    arm1_h2 = round(arm1$h2, 3),
    arm1_bias_pp = round(100 * arm1_bias, 1),
    sigma_theta_arm1 = round(sigma_theta_arm1, 3),

    # --- Misspecified ---
    arm2_h2 = round(arm2$h2, 3),
    arm2_r_mz = round(arm2$r_mz, 3),
    arm2_bias = round(arm2_bias, 3),
    arm2_bias_pp = round(100 * arm2_bias, 1),
    arm2_vs_arm1_pp = round(100 * (arm2$h2 - arm1$h2), 1),
    sigma_theta_true = round(sigma_theta_true, 4),
    sigma_theta_fit = round(sigma_theta_fit, 3),
    sigma_infl_pct = round(100 * (sigma_theta_fit / sigma_theta_true - 1), 1),

    # --- Control: nonfamilial ---
    ctrl_a_h2 = round(ctrl_a$h2, 3),
    ctrl_a_bias_pp = round(100 * (ctrl_a$h2 - oracle$h2), 1),
    sigma_theta_ctrl_a = round(sigma_theta_ctrl_a, 3),

    # --- Control: vanishing ---
    ctrl_b_h2 = round(ctrl_b$h2, 3),
    ctrl_b_bias_pp = round(100 * (ctrl_b$h2 - oracle$h2), 1),
    sigma_theta_ctrl_b = round(sigma_theta_ctrl_b, 3),

    # --- Comparator (Hamilton) ---
    arm3_h2 = round(arm3$h2, 3),
    arm3_bias_pp = round(100 * (arm3$h2 - oracle$h2), 1),

    # --- Variance decomposition ---
    var_decomp_mar_pct = round(var_decomp$mar_pct, 1),
    var_decomp_alpha = round(empirical_alpha_beta$alpha, 2),
    var_decomp_beta = round(empirical_alpha_beta$beta, 2),
    var_decomp_r2 = round(empirical_alpha_beta$r_squared, 3),

    # --- Sweep extremes ---
    sweep_max_bias_pp = round(100 * sweep_max$bias[1], 1),
    sweep_max_h2 = round(sweep_max$h2_shenhar[1], 2),
    sweep_max_sigma_infl = round(sweep_max$sigma_infl_pct[1], 0),

    # --- Anchored sweep ---
    anchored_min_pp = round(100 * min(anch_biases), 1),
    anchored_max_pp = round(100 * max(anch_biases), 1),
    anchored_mean_pp = round(100 * mean(anch_biases), 0),
    anchored_point_pp = round(100 * anch_point$bias[1], 1),

    # --- Check: no survivor cutoff (20-seed replicated) ---
    cutoff0_bias_pp = round(robustness_cutoff0$bias_pp, 1),
    cutoff0_se_pp = round(robustness_cutoff0$se_pp, 2),
    cutoff0_lo95_pp = round(robustness_cutoff0$lo95_pp, 1),
    cutoff0_hi95_pp = round(robustness_cutoff0$hi95_pp, 1),
    cutoff0_sigma_fit = round(robustness_cutoff0$sigma_fit, 3),
    cutoff0_h2 = round(robustness_cutoff0$corr_h2, 3),

    # --- Dose-response ---
    dose_hist_bias_pp = round(dr_hist$bias_pp[1], 1),
    dose_3x_bias_pp = round(dr_3x$bias_pp[1], 1),

    # --- Check: pleiotropy only ---
    pleiotropy_iso_h2 = round(pleiotropy_isolation$h2_corrected, 3),
    pleiotropy_iso_sigma_fit = round(pleiotropy_isolation$sigma_fit, 3),
    pleiotropy_iso_bias_pp = round(pleiotropy_isolation$bias_pp, 1),
    pleiotropy_iso_sigma_infl_pct = round(
      100 * (pleiotropy_isolation$sigma_fit / sigma_theta_true - 1), 1),

    # --- Control: irrelevant trait ---
    irrel_h2 = round(irrelevant_trait$h2_corrected, 3),
    irrel_sigma_fit = round(irrelevant_trait$sigma_fit, 3),
    irrel_bias_pp = round(irrelevant_trait$bias_pp, 1),

    # --- Dose-response 3x (with uncertainty) ---
    dose_3x_se_pp = round(dr_3x$se_pp[1], 2),
    dose_3x_lo95_pp = round(dr_3x$lo95_pp[1], 1),
    dose_3x_hi95_pp = round(dr_3x$hi95_pp[1], 1),

    # --- Model validation (MGG) ---
    mgg_oracle_h2 = round(mgg_validation$oracle_h2, 3),
    mgg_arm1_h2 = round(mgg_validation$arm1_h2, 3),
    mgg_arm1_bias_pp = round(100 * mgg_validation$arm1_bias, 1),
    mgg_arm2_h2 = round(mgg_validation$arm2_h2, 3),
    mgg_arm2_bias_pp = round(100 * mgg_validation$arm2_bias, 1),

    # --- Joint diagnostic ---
    joint_rdz_discrep = round(abs(joint_diagnostic$r_dz_discrepancy), 4),
    joint_sigma_fit = round(joint_diagnostic$sigma_fit_joint, 3),
    joint_bias_pp = round(100 * joint_diagnostic$joint_bias, 1),

    # --- Model validation (SR) ---
    sr_oracle_h2 = round(sr_validation$oracle_h2, 3),
    sr_arm1_h2 = round(sr_validation$arm1_h2, 3),
    sr_arm1_bias_pp = round(100 * sr_validation$arm1_bias, 1),
    sr_arm2_h2 = round(sr_validation$arm2_h2, 3),
    sr_arm2_bias_pp = round(100 * sr_validation$arm2_bias, 1),

    # --- σ_γ bridge ---
    bridge_sigma_gamma = round(sigma_gamma_bridge$bridge_sigma_gamma, 2),

    # --- m_ex split ---
    mex_split_50pct_pp = round(mex_split$bias_pp[mex_split$frac_heritable == 0.5], 1),
    mex_split_100pct_pp = round(mex_split$bias_pp[mex_split$frac_heritable == 1.0], 1),

    # --- Recovery (two-component refit) arm (bias relative to global oracle, not local) ---
    fix_sigma_theta = round(oracle_fix$sigma_theta_fix, 4),
    fix_sigma_recovery_pct = round(oracle_fix$sigma_recovery_pct, 1),
    fix_h2 = round(oracle_fix$fix_h2, 3),
    fix_bias_pp = round(100 * (oracle_fix$fix_h2 - oracle$h2), 1)
  ),

  # --- Multi-target calibration ---
  if (!is.null(multi_target_arm)) list(
    mt_h2 = round(multi_target_arm$mt_h2, 3),
    mt_bias_pp = multi_target_arm$mt_bias_pp,
    mt_sigma_fit = round(multi_target_arm$mt_sigma_fit, 4),
    mt_gof_rmse = round(multi_target_arm$mt_gof_rmse, 2),
    mt_sigma_infl_pct = round(multi_target_arm$mt_sigma_infl_pct, 1),
    mt_rmz_only_h2 = round(multi_target_arm$rmz_h2, 3),
    mt_rmz_only_bias_pp = multi_target_arm$rmz_bias_pp,
    mt_rmz_only_sigma_fit = round(multi_target_arm$rmz_sigma_fit, 4),
    mt_rmz_only_sigma_infl_pct = round(multi_target_arm$rmz_sigma_infl_pct, 1),
    mt_resid_q25 = round(multi_target_arm$mt_residuals["q25"], 1),
    mt_resid_mean = round(multi_target_arm$mt_residuals["mean_age"], 1),
    mt_resid_sd = round(multi_target_arm$mt_residuals["sd_age"], 1),
    mt_resid_q75 = round(multi_target_arm$mt_residuals["q75"], 1)
  ),

  # --- MGG parameterization comparison ---
  if (!is.null(mgg_param_comparison)) {
    ts <- mgg_param_comparison$twin_stats
    ds <- mgg_param_comparison$demo_stats
    comp <- ts[ts$mapping == "compensatory", ]
    paper <- ts[ts$mapping == "paper", ]
    list(
      mgg_comp_h2 = round(comp$h2, 3),
      mgg_comp_r_mz = round(comp$r_mz, 3),
      mgg_comp_r_dz = round(comp$r_dz, 3),
      mgg_paper_h2 = round(paper$h2, 3),
      mgg_paper_r_mz = round(paper$r_mz, 3),
      mgg_paper_r_dz = round(paper$r_dz, 3),
      mgg_comp_mean_age = round(ds$mean_age[ds$mapping == "compensatory"]),
      mgg_paper_mean_age = round(ds$mean_age[ds$mapping == "paper"]),
      mgg_param_h2_diff_pp = round(100 * (paper$h2 - comp$h2), 1)
    )
  },

  # --- B1: Alternative extrinsic frailty forms (multi-rep) ---
  if (!is.null(alt_ext_forms)) {
    b1_agg <- aggregate(bias_pp ~ ext_form, data = alt_ext_forms, FUN = mean)
    b1_se  <- aggregate(bias_pp ~ ext_form, data = alt_ext_forms,
      FUN = function(x) sd(x) / sqrt(length(x)))
    ln_mean <- b1_agg$bias_pp[b1_agg$ext_form == "lognormal"]
    ad_mean <- b1_agg$bias_pp[b1_agg$ext_form == "additive"]
    ga_mean <- b1_agg$bias_pp[b1_agg$ext_form == "gamma"]
    ln_se <- b1_se$bias_pp[b1_se$ext_form == "lognormal"]
    ad_se <- b1_se$bias_pp[b1_se$ext_form == "additive"]
    ga_se <- b1_se$bias_pp[b1_se$ext_form == "gamma"]
    list(
      b1_ln_bias_pp = round(ln_mean, 1),
      b1_ln_se_pp = round(ln_se, 2),
      b1_add_bias_pp = round(ad_mean, 1),
      b1_add_se_pp = round(ad_se, 2),
      b1_gamma_bias_pp = round(ga_mean, 1),
      b1_gamma_se_pp = round(ga_se, 2),
      b1_min_pp = round(min(b1_agg$bias_pp), 1),
      b1_max_pp = round(max(b1_agg$bias_pp), 1),
      b1_spread_pp = round(diff(range(b1_agg$bias_pp)), 1),
      b1_mean_pp = round(mean(b1_agg$bias_pp), 1),
      b1_n_reps = max(alt_ext_forms$rep)
    )
  },

  # --- B2: Extended σ_γ range (multi-rep) ---
  if (!is.null(extended_sigma_gamma)) {
    b2_agg <- aggregate(bias_pp ~ sigma_gamma, data = extended_sigma_gamma, FUN = mean)
    b2_se  <- aggregate(bias_pp ~ sigma_gamma, data = extended_sigma_gamma,
      FUN = function(x) sd(x) / sqrt(length(x)))
    b2_sinfl <- aggregate(sigma_infl_pct ~ sigma_gamma, data = extended_sigma_gamma,
      FUN = mean)
    sg070_m <- b2_agg$bias_pp[abs(b2_agg$sigma_gamma - 0.70) < 0.01]
    sg100_m <- b2_agg$bias_pp[abs(b2_agg$sigma_gamma - 1.00) < 0.01]
    bridge_m <- b2_agg$bias_pp[abs(b2_agg$sigma_gamma - 1.47) < 0.01]
    bridge_se <- b2_se$bias_pp[abs(b2_se$sigma_gamma - 1.47) < 0.01]
    bridge_sinfl <- b2_sinfl$sigma_infl_pct[abs(b2_sinfl$sigma_gamma - 1.47) < 0.01]
    b2_nreps <- max(extended_sigma_gamma$rep)
    b2_t_crit <- qt(0.975, df = b2_nreps - 1)
    list(
      b2_sg070_pp = round(sg070_m, 1),
      b2_sg100_pp = round(sg100_m, 1),
      b2_bridge_pp = round(bridge_m, 1),
      b2_bridge_se_pp = round(bridge_se, 2),
      b2_bridge_ci_half = round(b2_t_crit * bridge_se, 1),
      b2_bridge_sigma_infl = round(bridge_sinfl, 1),
      b2_n_reps = b2_nreps
    )
  },

  # --- B3: MC uncertainty bounds ---
  if (!is.null(mc_uncertainty)) {
    s <- mc_uncertainty$summary
    seeds <- mc_uncertainty$per_seed
    m <- mean(seeds$arm2_bias_pp)
    se <- sd(seeds$arm2_bias_pp) / sqrt(nrow(seeds))
    list(
      b3_n_seeds = nrow(seeds),
      b3_arm2_mean_pp = round(m, 2),
      b3_arm2_se = round(se, 2),
      b3_arm2_ci_lo = round(m - qt(0.975, df = nrow(seeds) - 1) * se, 1),
      b3_arm2_ci_hi = round(m + qt(0.975, df = nrow(seeds) - 1) * se, 1),
      b3_oracle_mean = round(s$mean[s$statistic == "oracle_h2"], 3),
      b3_oracle_sd = round(s$sd[s$statistic == "oracle_h2"], 3),
      b3_arm1_mean_pp = round(s$mean[s$statistic == "arm1_bias_pp"], 2),
      b3_arm1_sd_pp = round(s$sd[s$statistic == "arm1_bias_pp"], 2)
    )
  },

  # --- B3b: Multi-model MC uncertainty ---
  if (!is.null(mc_uncertainty_multimodel)) {
    mm <- mc_uncertainty_multimodel$summary
    mm_gm <- mm[mm$model == "GM", ]
    mm_mgg <- mm[mm$model == "MGG", ]
    mm_sr <- mm[mm$model == "SR", ]
    list(
      b3b_n_seeds = mm_gm$n_seeds,
      b3b_gm_arm2_mean_pp = round(mm_gm$arm2_bias_mean, 1),
      b3b_gm_arm2_se_pp = round(mm_gm$arm2_bias_se, 2),
      b3b_gm_arm1_mean_pp = round(mm_gm$arm1_bias_mean, 1),
      b3b_gm_arm1_se_pp = round(mm_gm$arm1_bias_se, 2),
      b3b_gm_oracle_h2 = round(mm_gm$oracle_h2_mean, 3),
      b3b_mgg_arm2_mean_pp = round(mm_mgg$arm2_bias_mean, 1),
      b3b_mgg_arm2_se_pp = round(mm_mgg$arm2_bias_se, 2),
      b3b_mgg_arm1_mean_pp = round(mm_mgg$arm1_bias_mean, 1),
      b3b_mgg_arm1_se_pp = round(mm_mgg$arm1_bias_se, 2),
      b3b_mgg_oracle_h2 = round(mm_mgg$oracle_h2_mean, 3),
      b3b_sr_arm2_mean_pp = round(mm_sr$arm2_bias_mean, 1),
      b3b_sr_arm2_se_pp = round(mm_sr$arm2_bias_se, 2),
      b3b_sr_arm1_mean_pp = round(mm_sr$arm1_bias_mean, 1),
      b3b_sr_arm1_se_pp = round(mm_sr$arm1_bias_se, 2),
      b3b_sr_oracle_h2 = round(mm_sr$oracle_h2_mean, 3),
      b3b_min_lo95_pp = {
        t_crit <- qt(0.975, df = mm_gm$n_seeds - 1)
        lo <- c(mm_gm$arm2_bias_mean - t_crit * mm_gm$arm2_bias_se,
                mm_mgg$arm2_bias_mean - t_crit * mm_mgg$arm2_bias_se,
                mm_sr$arm2_bias_mean - t_crit * mm_sr$arm2_bias_se)
        round(min(lo), 1)
      }
    )
  },

  # --- B4: m_ex split hires ---
  if (!is.null(mex_split_hires)) {
    agg <- aggregate(bias_pp ~ frac_heritable, data = mex_split_hires,
                     FUN = mean)
    list(
      b4_n_reps = max(mex_split_hires$rep),
      b4_frac0_pp = round(agg$bias_pp[agg$frac_heritable == 0], 1),
      b4_frac50_pp = round(agg$bias_pp[agg$frac_heritable == 0.5], 1),
      b4_frac100_pp = round(agg$bias_pp[agg$frac_heritable == 1.0], 1)
    )
  },

  # --- B5: Bivariate survival model check ---
  if (!is.null(bivariate_check)) list(
    biv_iae_misspec_mz = round(bivariate_check$iae_misspec_mz, 4),
    biv_iae_misspec_dz = round(bivariate_check$iae_misspec_dz, 4),
    biv_iae_fix_mz = round(bivariate_check$iae_fix_mz, 4),
    biv_iae_fix_dz = round(bivariate_check$iae_fix_dz, 4),
    biv_iae_ratio_mz = round(
      bivariate_check$iae_misspec_mz / max(bivariate_check$iae_fix_mz, 1e-10), 1),
    biv_peak_dev_age_mz = bivariate_check$peak_dev_age_mz,
    biv_peak_dev_age_dz = bivariate_check$peak_dev_age_dz,
    biv_n_mz = bivariate_check$n_mz_true,
    biv_n_dz = bivariate_check$n_dz_true
  ),

  # --- B6: MZ extrinsic concordance sensitivity ---
  if (!is.null(mz_concordance)) {
    d <- mz_concordance[order(mz_concordance$r_gamma_mz), ]
    # Zero crossing via first sign change in ordered data
    diffs <- diff(sign(d$bias_pp))
    sign_changes <- which(diffs != 0)
    zero_cross <- if (length(sign_changes) > 0) {
      i <- sign_changes[1]
      w <- (0 - d$bias_pp[i]) / (d$bias_pp[i + 1] - d$bias_pp[i])
      d$r_gamma_mz[i] + w * (d$r_gamma_mz[i + 1] - d$r_gamma_mz[i])
    } else NA_real_
    # Interpolated bias at selected concordance levels
    get_bias <- function(r) {
      tryCatch(
        approx(d$r_gamma_mz, d$bias_pp, xout = r)$y,
        error = function(e) NA_real_
      )
    }
    list(
      conc_bias_pp_00 = round(get_bias(0.0), 1),
      conc_bias_pp_05 = round(get_bias(0.5), 1),
      conc_bias_pp_07 = round(get_bias(0.7), 1),
      conc_bias_pp_085 = round(get_bias(0.85), 1),
      conc_bias_pp_10 = round(get_bias(1.0), 1),
      conc_zero_crossing = if (is.na(zero_cross)) "none" else round(zero_cross, 2)
    )
  },

  # --- B6b: Decoupled concordance sweep (rho=0) ---
  if (!is.null(mz_concordance_decoupled)) {
    dd <- mz_concordance_decoupled[order(mz_concordance_decoupled$r_gamma_mz), ]
    get_bias_d <- function(r) {
      tryCatch(
        approx(dd$r_gamma_mz, dd$bias_pp, xout = r)$y,
        error = function(e) NA_real_
      )
    }
    list(
      conc_rho0_bias_pp_00 = round(get_bias_d(0.0), 1),
      conc_rho0_bias_pp_05 = round(get_bias_d(0.5), 1),
      conc_rho0_bias_pp_07 = round(get_bias_d(0.7), 1),
      conc_rho0_bias_pp_085 = round(get_bias_d(0.85), 1),
      conc_rho0_bias_pp_10 = round(get_bias_d(1.0), 1)
    )
  },

  # --- B7: Bridge uncertainty ---
  if (!is.null(bridge_uncertainty)) {
    s <- bridge_uncertainty$summary
    central_se <- s$bridge_se[s$target_label == "central"]
    central_n <- s$n_valid[s$target_label == "central"]
    central_t <- if (!is.na(central_n) && central_n > 1) qt(0.975, df = central_n - 1) else NA_real_
    list(
      bridge_unc_central_mean = round(s$bridge_mean[s$target_label == "central"], 2),
      bridge_unc_central_se = round(central_se, 2),
      bridge_unc_central_ci_half = if (!is.na(central_t)) round(central_t * central_se, 2) else NA_real_,
      bridge_unc_lower_mean = round(s$bridge_mean[s$target_label == "lower"], 2),
      bridge_unc_upper_mean = if (is.na(s$bridge_mean[s$target_label == "upper"])) "NA" else round(s$bridge_mean[s$target_label == "upper"], 2),
      bridge_all_above_065 = all(s$all_above_065, na.rm = TRUE)
    )
  },

  # --- B8: Replicated main arms (50-seed canonical numbers) ---
  if (!is.null(main_arms_replicated)) {
    s <- main_arms_replicated$summary
    get_row <- function(metric) {
      row <- s[s$metric == metric, ]
      if (nrow(row) == 0) return(list(mean = NA, se = NA, lo95 = NA, hi95 = NA))
      list(mean = row$mean[1], se = row$se[1], lo95 = row$lo95[1], hi95 = row$hi95[1])
    }
    scalars <- list()
    for (arm in c("arm1", "arm2", "ctrl_a", "ctrl_b", "fix")) {
      r <- get_row(paste0(arm, "_bias_pp"))
      scalars[[paste0("rep_", arm, "_bias_pp")]] <- round(r$mean, 1)
      scalars[[paste0("rep_", arm, "_se_pp")]] <- round(r$se, 2)
      scalars[[paste0("rep_", arm, "_lo95_pp")]] <- round(r$lo95, 1)
      scalars[[paste0("rep_", arm, "_hi95_pp")]] <- round(r$hi95, 1)
    }
    # Replicated oracle and arm2 h²
    ora <- get_row("oracle_h2")
    a2 <- get_row("arm2_h2")
    scalars$rep_oracle_h2 <- sprintf("%.3f", ora$mean)
    scalars$rep_arm2_h2 <- sprintf("%.3f", a2$mean)
    # Sigma inflation
    si <- get_row("arm2_sigma_infl")
    scalars$rep_arm2_sigma_infl <- round(si$mean, 1)
    # Pleiotropy isolation (replicated)
    for (arm in c("pleio_iso", "irrel", "arm3")) {
      r <- get_row(paste0(arm, "_bias_pp"))
      scalars[[paste0("rep_", arm, "_bias_pp")]] <- round(r$mean, 1)
      scalars[[paste0("rep_", arm, "_se_pp")]] <- round(r$se, 2)
      scalars[[paste0("rep_", arm, "_lo95_pp")]] <- round(r$lo95, 1)
      scalars[[paste0("rep_", arm, "_hi95_pp")]] <- round(r$hi95, 1)
    }
    # Replicated h2 values for all arms
    for (arm in c("pleio_iso", "irrel", "arm3")) {
      rr <- get_row(paste0(arm, "_h2"))
      scalars[[paste0("rep_", arm, "_h2")]] <- sprintf("%.3f", rr$mean)
    }
    scalars
  }
  )
}
