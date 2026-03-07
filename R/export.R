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
#' @param arm1 Arm 1 result
#' @param arm2 Arm 2 result
#' @param ctrl_a Control A result
#' @param ctrl_b Control B result
#' @param arm3 Arm 3 (Hamilton) result
#' @param sigma_theta_true True sigma_theta
#' @param sigma_theta_arm1 Calibrated sigma for arm 1
#' @param sigma_theta_arm2 Calibrated sigma for arm 2
#' @param sigma_theta_ctrl_a Calibrated sigma for control A
#' @param sigma_theta_ctrl_b Calibrated sigma for control B
#' @param oracle_fix Oracle fix arm result
#' @param pleiotropy_isolation Pleiotropy isolation result
#' @param irrelevant_trait Irrelevant trait result
#' @param robustness_cutoff0 Cutoff=0 result
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
    "Oracle (intrinsic only)",
    "Arm 1: Correctly specified",
    "Arm 2: Misspecification",
    "Control A: Non-heritable extrinsic",
    "Control B: Vanishing m_ex"
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

  # Pleiotropy isolation
  if (!is.null(pleiotropy_isolation)) {
    conditions <- c(conditions, "Pleiotropy isolation (rho=0)")
    sigma_vals <- c(sigma_vals, pleiotropy_isolation$sigma_fit)
    h2_vals <- c(h2_vals, pleiotropy_isolation$h2_corrected)
    bias_vals <- c(bias_vals, pleiotropy_isolation$bias)
  }

  # Irrelevant trait
  if (!is.null(irrelevant_trait)) {
    conditions <- c(conditions, "Irrelevant trait (zeta)")
    sigma_vals <- c(sigma_vals, irrelevant_trait$sigma_fit)
    h2_vals <- c(h2_vals, irrelevant_trait$h2_corrected)
    bias_vals <- c(bias_vals, irrelevant_trait$bias)
  }

  # Oracle fix — bias relative to the GLOBAL oracle (h2_vals[1]), not local
  if (!is.null(oracle_fix)) {
    conditions <- c(conditions, "Oracle fix: correct model")
    sigma_vals <- c(sigma_vals, oracle_fix$sigma_theta_fix)
    h2_vals <- c(h2_vals, oracle_fix$fix_h2)
    bias_vals <- c(bias_vals, oracle_fix$fix_h2 - h2_vals[1])
  }

  # Cutoff = 0
  if (!is.null(robustness_cutoff0)) {
    conditions <- c(conditions, "Cutoff = 0")
    sigma_vals <- c(sigma_vals, robustness_cutoff0$sigma_fit)
    h2_vals <- c(h2_vals, robustness_cutoff0$corr_h2)
    bias_vals <- c(bias_vals, robustness_cutoff0$bias_pp / 100)
  }

  # Arm 3 always last
  conditions <- c(conditions, "Arm 3: Hamilton conditioning")
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
    "Oracle (intrinsic only)",
    "Arm 1: Correctly specified",
    "Arm 2: Misspecification",
    "Control A: Non-heritable extrinsic",
    "Irrelevant trait (zeta)",
    "Control B: Vanishing m_ex",
    "Pleiotropy isolation (rho=0)",
    "Oracle fix: correct model",
    "Cutoff = 0",
    "Arm 3: Hamilton conditioning"
  )
  ord <- match(df$Condition, canonical)
  df <- df[order(ord), ]
  rownames(df) <- NULL
  df
}

#' Build the model validation table
#'
#' @param oracle_gm GM oracle result
#' @param arm1_gm GM Arm 1 corrected result
#' @param arm2_gm GM Arm 2 corrected result
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
    Arm1_h2 = c(arm1_gm$h2, mgg_result$arm1_h2, sr_result$arm1_h2),
    Arm1_bias_pp = c(round(100 * gm_arm1_bias, 1),
                     round(100 * mgg_result$arm1_bias, 1),
                     round(100 * sr_result$arm1_bias, 1)),
    Arm2_h2 = c(arm2_gm$h2, mgg_result$arm2_h2, sr_result$arm2_h2),
    Arm2_bias_pp = c(round(100 * gm_arm2_bias, 1),
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
                            mgg_param_comparison = NULL) {

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

    # --- Arm 1 ---
    arm1_h2 = round(arm1$h2, 3),
    arm1_bias_pp = round(100 * arm1_bias, 1),
    sigma_theta_arm1 = round(sigma_theta_arm1, 3),

    # --- Arm 2 ---
    arm2_h2 = round(arm2$h2, 3),
    arm2_r_mz = round(arm2$r_mz, 3),
    arm2_bias = round(arm2_bias, 3),
    arm2_bias_pp = round(100 * arm2_bias, 1),
    arm2_vs_arm1_pp = round(100 * (arm2$h2 - arm1$h2), 1),
    sigma_theta_true = round(sigma_theta_true, 4),
    sigma_theta_fit = round(sigma_theta_fit, 3),
    sigma_infl_pct = round(100 * (sigma_theta_fit / sigma_theta_true - 1), 1),

    # --- Control A ---
    ctrl_a_h2 = round(ctrl_a$h2, 3),
    ctrl_a_bias_pp = round(100 * (ctrl_a$h2 - oracle$h2), 1),
    sigma_theta_ctrl_a = round(sigma_theta_ctrl_a, 3),

    # --- Control B ---
    ctrl_b_h2 = round(ctrl_b$h2, 3),
    ctrl_b_bias_pp = round(100 * (ctrl_b$h2 - oracle$h2), 1),
    sigma_theta_ctrl_b = round(sigma_theta_ctrl_b, 3),

    # --- Arm 3 (Hamilton) ---
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

    # --- Cutoff=0 robustness ---
    cutoff0_bias_pp = round(robustness_cutoff0$bias_pp, 1),
    cutoff0_sigma_fit = round(robustness_cutoff0$sigma_fit, 3),
    cutoff0_h2 = round(robustness_cutoff0$corr_h2, 3),

    # --- Dose-response ---
    dose_hist_bias_pp = round(100 * dr_hist$bias[1], 1),
    dose_3x_bias_pp = round(100 * dr_3x$bias[1], 1),

    # --- Pleiotropy isolation ---
    pleiotropy_iso_h2 = round(pleiotropy_isolation$h2_corrected, 3),
    pleiotropy_iso_sigma_fit = round(pleiotropy_isolation$sigma_fit, 3),
    pleiotropy_iso_bias_pp = round(pleiotropy_isolation$bias_pp, 1),
    pleiotropy_iso_sigma_infl_pct = round(
      100 * (pleiotropy_isolation$sigma_fit / sigma_theta_true - 1), 1),

    # --- Irrelevant trait ---
    irrel_h2 = round(irrelevant_trait$h2_corrected, 3),
    irrel_sigma_fit = round(irrelevant_trait$sigma_fit, 3),
    irrel_bias_pp = round(irrelevant_trait$bias_pp, 1),

    # --- Dose-response large m_ex ---
    dose_3x_h2 = round(dr_3x$h2_corrected[1], 3),
    dose_3x_sigma_fit = round(dr_3x$sigma_fit[1], 3),

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

    # --- Oracle fix arm (bias relative to global oracle, not local) ---
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
  }
  )
}
