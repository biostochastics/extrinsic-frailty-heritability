# ===========================================================================
# _targets.R — Shenhar-Hamilton-Kornilov Pipeline
# ===========================================================================
# Reproducible simulation pipeline using {targets}.
# All R/ modules are sourced automatically.
# ===========================================================================

library(targets)
library(tarchetypes)

# Source all R/ modules
tar_source("R")

# Parallel execution via crew
tar_option_set(
  packages = c("MASS", "ggplot2", "patchwork", "lamW",
               "ggpubr", "ggthemes", "viridis", "mvtnorm", "scales", "Rcpp"),
  controller = crew::crew_controller_local(workers = 10L),
  deployment = "worker"
)

# Static grids (evaluated at definition time for tar_map)
sweep_grid_static <- make_sweep_grid()
anchored_grid_static <- make_anchored_grid()

# Mapped targets
sweep_mapped <- tar_map(
  values = sweep_grid_static,
  names = tidyselect::all_of(c("rho_val", "sigma_gamma_val")),
  tar_target(sweep_cell, run_sweep_cell(
    sigma_theta_true, rho_val, sigma_gamma_val, oracle, PARAMS
  ))
)

anchored_mapped <- tar_map(
  values = anchored_grid_static,
  names = tidyselect::all_of(c("rho_val", "sigma_gamma_val")),
  tar_target(anchored_cell, run_sweep_cell(
    sigma_theta_true, rho_val, sigma_gamma_val, oracle, PARAMS
  ))
)

list(
  # =================================================================
  # Stage 0: Calibration
  # =================================================================
  tar_target(sigma_theta_true, calibrate_oracle(PARAMS$TARGET_H2)),

  tar_target(oracle, sim_twin_h2(sigma_theta_true, 0,
    sigma_gamma = 0, rng_seed = PARAMS$MASTER_SEED)),

  # =================================================================
  # Stage 1: Main arms
  # =================================================================
  tar_target(obs_null, sim_twin_h2(sigma_theta_true, PARAMS$M_EX_HIST,
    sigma_gamma = 0, rng_seed = PARAMS$MASTER_SEED + 100L)),

  tar_target(sigma_theta_arm1,
    calibrate_sigma_theta(obs_null$r_mz, PARAMS$M_EX_HIST)),

  tar_target(arm1_corr, sim_twin_h2(sigma_theta_arm1, 0,
    sigma_gamma = 0, rng_seed = PARAMS$MASTER_SEED + 200L)),

  tar_target(true_obs, sim_twin_h2(sigma_theta_true, PARAMS$M_EX_HIST,
    sigma_gamma = PARAMS$SIGMA_GAMMA_DEF, rho = PARAMS$RHO_DEF,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = PARAMS$MASTER_SEED + 300L)),

  tar_target(sigma_theta_fit,
    calibrate_sigma_theta(true_obs$r_mz, PARAMS$M_EX_HIST)),

  tar_target(arm2_corr, sim_twin_h2(sigma_theta_fit, 0,
    sigma_gamma = 0, rng_seed = PARAMS$MASTER_SEED + 400L)),

  # --- Control A: Non-heritable extrinsic ---
  tar_target(ctrl_a_obs, sim_twin_h2(sigma_theta_true, PARAMS$M_EX_HIST,
    sigma_gamma = PARAMS$SIGMA_GAMMA_DEF, rho = 0,
    gamma_heritable = FALSE, individual_ext = TRUE,
    rng_seed = PARAMS$MASTER_SEED + 500L)),

  tar_target(sigma_theta_ctrl_a,
    calibrate_sigma_theta(ctrl_a_obs$r_mz, PARAMS$M_EX_HIST)),

  tar_target(ctrl_a_corr, sim_twin_h2(sigma_theta_ctrl_a, 0,
    sigma_gamma = 0, rng_seed = PARAMS$MASTER_SEED + 600L)),

  # --- Control B: Vanishing m_ex ---
  tar_target(ctrl_b_obs, sim_twin_h2(sigma_theta_true, PARAMS$M_EX_TINY,
    sigma_gamma = PARAMS$SIGMA_GAMMA_DEF, rho = PARAMS$RHO_DEF,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = PARAMS$MASTER_SEED + 700L)),

  tar_target(sigma_theta_ctrl_b,
    calibrate_sigma_theta(ctrl_b_obs$r_mz, PARAMS$M_EX_TINY)),

  tar_target(ctrl_b_corr, sim_twin_h2(sigma_theta_ctrl_b, 0,
    sigma_gamma = 0, rng_seed = PARAMS$MASTER_SEED + 800L)),

  # --- Arm 3: Hamilton ---
  tar_target(arm3_hamilton, run_hamilton_arm(PARAMS)),

  # --- Oracle fix arm: correctly specified two-component model ---
  tar_target(oracle_fix, run_oracle_fix_arm(sigma_theta_true, PARAMS)),

  # --- Multi-target calibration robustness ---
  tar_target(multi_target_arm,
    run_multi_target_arm(sigma_theta_true, oracle, PARAMS)),

  # =================================================================
  # Stage 2: Sensitivity sweep (PARALLEL via static branching)
  # =================================================================
  sweep_mapped,

  tar_combine(sweep_results, sweep_mapped,
    command = do.call(rbind, lapply(list(!!!.x), as.data.frame))),

  # =================================================================
  # Stage 3: Anchored sweep (PARALLEL via static branching)
  # =================================================================
  anchored_mapped,

  tar_combine(anchored_results, anchored_mapped,
    command = do.call(rbind, lapply(list(!!!.x), as.data.frame))),

  # =================================================================
  # Stage 4: Variance decomposition
  # =================================================================
  tar_target(var_decomp,
    variance_decomposition(sweep_results, sigma_theta_true)),

  # =================================================================
  # Stage 5: Model validation (MGG + SR)
  # =================================================================
  tar_target(mgg_validation, run_model_analysis("mgg", PARAMS)),
  tar_target(sr_validation, run_model_analysis("sr", PARAMS)),

  tar_target(model_table, build_model_table(
    oracle, arm1_corr, arm2_corr, mgg_validation, sr_validation)),

  # =================================================================
  # Stage 6: Additional controls
  # =================================================================
  tar_target(dose_response,
    run_dose_response(sigma_theta_true, oracle, PARAMS)),

  tar_target(pleiotropy_isolation,
    run_pleiotropy_isolation(sigma_theta_true, oracle, PARAMS)),

  tar_target(irrelevant_trait,
    run_irrelevant_trait(sigma_theta_true, oracle, PARAMS)),

  # =================================================================
  # Stage 6b: MGG/SR control extensions
  # =================================================================
  tar_target(negative_rho_mgg,
    run_model_control_sweep("mgg", "negative_rho", PARAMS)),

  tar_target(negative_rho_sr,
    run_model_control_sweep("sr", "negative_rho", PARAMS)),

  tar_target(pleiotropy_isolation_mgg,
    run_model_control_sweep("mgg", "pleiotropy_isolation", PARAMS)),

  tar_target(pleiotropy_isolation_sr,
    run_model_control_sweep("sr", "pleiotropy_isolation", PARAMS)),

  tar_target(dose_response_mgg,
    run_model_control_sweep("mgg", "dose_response", PARAMS)),

  tar_target(dose_response_sr,
    run_model_control_sweep("sr", "dose_response", PARAMS)),

  # =================================================================
  # Stage 6c: Robustness (cutoff = 0)
  # =================================================================
  tar_target(robustness_cutoff0,
    run_cutoff0_robustness(sigma_theta_true, PARAMS)),

  # =================================================================
  # Stage 7: New diagnostics (Batch 4)
  # =================================================================
  tar_target(joint_diagnostic,
    run_joint_rmz_rdz_diagnostic(sigma_theta_true, PARAMS)),

  tar_target(sigma_gamma_bridge,
    run_sigma_gamma_bridge(sigma_theta_true, PARAMS)),

  tar_target(mex_split,
    run_mex_split_sensitivity(sigma_theta_true, PARAMS)),

  tar_target(mgg_hazard_curves,
    run_mgg_hazard_curves(PARAMS)),

  tar_target(mgg_param_comparison,
    run_mgg_param_comparison(PARAMS)),

  tar_target(negative_rho,
    run_negative_rho_sensitivity(sigma_theta_true, PARAMS)),

  tar_target(empirical_alpha_beta,
    estimate_alpha_beta(sweep_results, sigma_theta_true)),

  tar_target(sr_dt_sensitivity,
    run_sr_dt_check(PARAMS)),

  # =================================================================
  # Stage 8: Calibration UQ
  # =================================================================
  tar_target(calib_uq_oracle, calibrate_with_uncertainty(
    target_h2 = PARAMS$TARGET_H2, calibration_type = "oracle_h2")),

  tar_target(calib_uq_arm2, calibrate_with_uncertainty(
    target_r_mz = true_obs$r_mz, m_ex = PARAMS$M_EX_HIST,
    calibration_type = "r_mz")),

  # =================================================================
  # Stage 8b: Revision analyses
  # =================================================================

  # B1: Alternative functional forms for extrinsic frailty DGP (10 reps)
  tar_target(alt_ext_forms,
    run_alt_ext_forms(sigma_theta_true, oracle, PARAMS, n_reps = 10L)),

  # B2: Extended sigma_gamma range (above anchored ceiling, 10 reps)
  tar_target(extended_sigma_gamma,
    run_extended_sigma_gamma(sigma_theta_true, oracle, PARAMS, n_reps = 10L)),

  # B3: Monte Carlo uncertainty bounds (20 seeds)
  tar_target(mc_uncertainty,
    run_mc_uncertainty(sigma_theta_true, PARAMS, n_seeds = 20L)),

  # B3b: Multi-model MC uncertainty (split for crew parallelism)
  tar_target(mc_unc_gm,
    run_mc_uncertainty_model("gm", sigma_theta_true, PARAMS, n_seeds = 20L)),
  tar_target(mc_unc_mgg,
    run_mc_uncertainty_model("mgg", sigma_theta_true, PARAMS, n_seeds = 20L)),
  tar_target(mc_unc_sr,
    run_mc_uncertainty_model("sr", sigma_theta_true, PARAMS, n_seeds = 20L)),
  tar_target(mc_uncertainty_multimodel,
    combine_mc_uncertainty_multimodel(mc_unc_gm, mc_unc_mgg, mc_unc_sr)),

  # B4: High-replication m_ex split (5 reps per fraction level)
  tar_target(mex_split_hires,
    run_mex_split_hires(sigma_theta_true, oracle, PARAMS, n_reps = 5L)),

  # =================================================================
  # Stage 9: Tables
  # =================================================================
  tar_target(summary_table, build_summary_table(
    oracle, arm1_corr, arm2_corr, ctrl_a_corr, ctrl_b_corr, arm3_hamilton,
    sigma_theta_true, sigma_theta_arm1, sigma_theta_fit,
    sigma_theta_ctrl_a, sigma_theta_ctrl_b,
    oracle_fix = oracle_fix,
    pleiotropy_isolation = pleiotropy_isolation,
    irrelevant_trait = irrelevant_trait,
    robustness_cutoff0 = robustness_cutoff0)),

  tar_target(scalars, collect_scalars(
    oracle, arm1_corr, arm2_corr, ctrl_a_corr, ctrl_b_corr, arm3_hamilton,
    sigma_theta_true, sigma_theta_arm1,
    sigma_theta_fit, sigma_theta_ctrl_a, sigma_theta_ctrl_b,
    var_decomp, sweep_results, anchored_results,
    robustness_cutoff0, dose_response,
    pleiotropy_isolation, irrelevant_trait,
    mgg_validation, sr_validation,
    empirical_alpha_beta,
    joint_diagnostic,
    sigma_gamma_bridge, mex_split,
    oracle_fix,
    multi_target_arm = multi_target_arm,
    mgg_param_comparison = mgg_param_comparison,
    alt_ext_forms = alt_ext_forms,
    extended_sigma_gamma = extended_sigma_gamma,
    mc_uncertainty = mc_uncertainty,
    mc_uncertainty_multimodel = mc_uncertainty_multimodel,
    mex_split_hires = mex_split_hires)),

  tar_target(controls_table, data.frame(
    Control = c(
      "Arm 2 reference",
      "Pleiotropy isolation (rho=0)",
      "Irrelevant trait (zeta)",
      "Vanishing m_ex (ctrl B)"
    ),
    bias = c(
      arm2_corr$h2 - oracle$h2,
      pleiotropy_isolation$bias,
      irrelevant_trait$bias,
      ctrl_b_corr$h2 - oracle$h2
    ),
    bias_pp = c(
      100 * (arm2_corr$h2 - oracle$h2),
      pleiotropy_isolation$bias_pp,
      irrelevant_trait$bias_pp,
      100 * (ctrl_b_corr$h2 - oracle$h2)
    ),
    stringsAsFactors = FALSE
  )),

  tar_target(model_controls_combined, rbind(
    negative_rho_mgg, negative_rho_sr,
    pleiotropy_isolation_mgg, pleiotropy_isolation_sr,
    dose_response_mgg, dose_response_sr
  )),

  tar_target(write_tables_done, write_all_tables(
    summary_table = summary_table,
    model_table = model_table,
    controls_table = controls_table,
    sweep_table = sweep_results,
    anchored_table = anchored_results,
    model_controls = model_controls_combined,
    scalars = scalars)),

  # =================================================================
  # Stage 10: Figures
  # =================================================================
  tar_target(all_figures, generate_all_figures(
    summary_table = summary_table,
    sweep_results = sweep_results,
    anchored_results = anchored_results,
    model_table = model_table,
    controls_table = controls_table,
    dose_response = dose_response,
    negative_rho = negative_rho,
    mex_split = mex_split,
    sigma_gamma_bridge = sigma_gamma_bridge,
    mgg_hazard_curves = mgg_hazard_curves,
    var_decomp = var_decomp,
    empirical_alpha_beta = empirical_alpha_beta,
    joint_diagnostic = joint_diagnostic,
    sr_dt_sensitivity = sr_dt_sensitivity,
    model_controls = model_controls_combined,
    mgg_param_comparison = mgg_param_comparison,
    alt_ext_forms = alt_ext_forms,
    extended_sigma_gamma = extended_sigma_gamma,
    mc_uncertainty = mc_uncertainty,
    mc_uncertainty_multimodel = mc_uncertainty_multimodel,
    mex_split_hires = mex_split_hires
  ))
)
