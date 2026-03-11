# ===========================================================================
# utils.R — Shared utilities: grid builders, arm runners, Hamilton sim
# ===========================================================================

#' Run a function over multiple independent seeds and summarize
#'
#' Wraps any single-seed simulation cycle in an MC outer loop.
#' Returns mean, SE, and 95% CI for a numeric result vector.
#'
#' @param fn Function that takes `rng_seed` and returns a named numeric vector
#' @param n_reps Number of independent replications
#' @param seed_base Base seed (each rep uses seed_base + rep * 7L)
#' @return List with mean, se, lo95, hi95 (named vectors), and raw matrix
replicate_arm <- function(fn, n_reps, seed_base) {
  raw <- lapply(seq_len(n_reps), function(rep) {
    seed <- seed_base + rep * 7L
    fn(rng_seed = seed)
  })
  mat <- do.call(rbind, raw)
  means <- colMeans(mat)
  ses <- apply(mat, 2, sd) / sqrt(n_reps)
  list(
    mean = means,
    se = ses,
    lo95 = means - qt(0.975, df = n_reps - 1) * ses,
    hi95 = means + qt(0.975, df = n_reps - 1) * ses,
    raw = mat
  )
}

#' Build the main sensitivity sweep grid
#' @return Data frame with rho_val, sigma_gamma_val columns (for tar_map)
make_sweep_grid <- function() {
  rho_grid <- seq(0, 0.8, by = 0.1)
  sgamma_grid <- seq(0.15, 0.55, by = 0.05)
  g <- expand.grid(rho_val = rho_grid, sigma_gamma_val = sgamma_grid)
  g
}

#' Build the anchored parameter sweep grid
#' @return Data frame with rho_val, sigma_gamma_val columns
make_anchored_grid <- function() {
  rho_anchored <- seq(0.20, 0.50, by = 0.05)
  sgamma_anchored <- seq(0.30, 0.65, by = 0.05)
  g <- expand.grid(rho_val = rho_anchored, sigma_gamma_val = sgamma_anchored)
  g
}

#' Run a single sweep cell: observe under true DGP, calibrate, extrapolate
#'
#' @param sigma_theta_true True sigma_theta
#' @param rho_val Pleiotropy correlation for this cell
#' @param sigma_gamma_val Extrinsic heterogeneity SD for this cell
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return Data frame row with results
run_sweep_cell <- function(sigma_theta_true, rho_val, sigma_gamma_val,
                           oracle, params) {
  seed_base <- params$MASTER_SEED +
    round(rho_val * 1000) + round(sigma_gamma_val * 10000)

  obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
    sigma_gamma = sigma_gamma_val, rho = rho_val,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = seed_base
  )

  sigma_fit <- calibrate_sigma_theta_crn(obs$r_mz, params$M_EX_HIST,
    crn_seed = seed_base + 12345L)

  corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
    rng_seed = seed_base + 50000L
  )

  data.frame(
    rho = rho_val,
    sigma_gamma = sigma_gamma_val,
    sigma_fit = sigma_fit,
    h2_shenhar = corr$h2,
    bias = corr$h2 - oracle$h2,
    sigma_infl_pct = 100 * (sigma_fit / sigma_theta_true - 1),
    obs_r_mz = obs$r_mz,
    obs_r_dz = obs$r_dz,
    stringsAsFactors = FALSE
  )
}

#' Run cutoff=0 robustness check (20-seed MC replication)
#'
#' @param sigma_theta_true True sigma_theta
#' @param params Parameter list
#' @param n_seeds Number of independent seeds
#' @return List with replicated mean/SE/CI (backward-compatible field names)
run_cutoff0_robustness <- function(sigma_theta_true, params, n_seeds = 20L) {
  seed_base <- params$MASTER_SEED + 9000L

  results <- lapply(seq_len(n_seeds), function(s) {
    sb <- seed_base + s * 10000L

    oracle_nc <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0,
      cutoff = 0, rng_seed = sb)

    true_obs_nc <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
      gamma_heritable = TRUE, individual_ext = TRUE,
      cutoff = 0, rng_seed = sb + 100L)

    sigma_fit_nc <- calibrate_sigma_theta_crn(true_obs_nc$r_mz, params$M_EX_HIST,
      crn_seed = sb + 150L)

    corr_nc <- sim_twin_h2(sigma_fit_nc, 0, sigma_gamma = 0,
      cutoff = 0, rng_seed = sb + 200L)

    data.frame(
      seed = s,
      oracle_h2 = oracle_nc$h2,
      corr_h2 = corr_nc$h2,
      bias_pp = 100 * (corr_nc$h2 - oracle_nc$h2),
      sigma_fit = sigma_fit_nc,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, results)
  m <- mean(df$bias_pp); se <- sd(df$bias_pp) / sqrt(nrow(df))

  # Use t-distribution for proper CI with n=20
  t_crit <- qt(0.975, df = nrow(df) - 1)

  # SE of corrected h2 (for error bars on the h2 bar itself)
  h2_se <- sd(df$corr_h2) / sqrt(nrow(df))
  h2_mean <- mean(df$corr_h2)

  list(
    per_seed = df,
    n_seeds = n_seeds,
    oracle_h2 = mean(df$oracle_h2),
    corr_h2 = h2_mean,
    corr_h2_se = h2_se,
    corr_h2_lo95 = h2_mean - t_crit * h2_se,
    corr_h2_hi95 = h2_mean + t_crit * h2_se,
    bias = mean(df$bias_pp) / 100,
    bias_pp = m,
    se_pp = se,
    lo95_pp = m - t_crit * se,
    hi95_pp = m + t_crit * se,
    sigma_fit = mean(df$sigma_fit)
  )
}

#' Run Hamilton's concordant-survivor conditioning arm
#'
#' @param params Parameter list
#' @return List with r_mz, r_dz, h2 (full and survivor-conditioned)
run_hamilton_arm <- function(params) {
  set.seed(params$MASTER_SEED + 8000L)
  n <- params$N_PAIRS

  hamilton_sim <- function(n, r_g, g_weight) {
    g1 <- rnorm(n)
    g2 <- r_g * g1 + sqrt(1 - r_g^2) * rnorm(n)
    e1 <- rnorm(n)
    e2 <- rnorm(n)
    f1 <- g_weight * g1 + params$HAM_E_WEIGHT * e1
    f2 <- g_weight * g2 + params$HAM_E_WEIGHT * e2
    p1 <- 1 / (1 + exp(-params$HAM_EXT_SLOPE * (f1 - params$HAM_EXT_SHIFT)))
    p2 <- 1 / (1 + exp(-params$HAM_EXT_SLOPE * (f2 - params$HAM_EXT_SHIFT)))
    d1 <- runif(n) < p1
    d2 <- runif(n) < p2
    L1 <- 80 - 5 * f1 + rnorm(n, 0, 5.8)
    L2 <- 80 - 5 * f2 + rnorm(n, 0, 5.8)
    if (sum(d1) > 0) L1[d1] <- 30 + rnorm(sum(d1), 0, 4.8)
    if (sum(d2) > 0) L2[d2] <- 30 + rnorm(sum(d2), 0, 4.8)
    list(L1 = L1, L2 = L2, d1 = d1, d2 = d2)
  }

  ham_calibrate <- function(target_h2, lo = 0.1, hi = 2.0) {
    for (i in 1:40) {
      mid <- (lo + hi) / 2
      mz <- hamilton_sim(n, 1.0, mid)
      dz <- hamilton_sim(n, 0.5, mid)
      k_mz <- !mz$d1 & !mz$d2
      k_dz <- !dz$d1 & !dz$d2
      h2c <- 2 * (cor(mz$L1[k_mz], mz$L2[k_mz]) -
                   cor(dz$L1[k_dz], dz$L2[k_dz]))
      if (h2c < target_h2) lo <- mid else hi <- mid
    }
    (lo + hi) / 2
  }

  ham_g <- ham_calibrate(params$HAM_TARGET_H2)
  mz_h <- hamilton_sim(n, 1.0, ham_g)
  dz_h <- hamilton_sim(n, 0.5, ham_g)

  r_mz_full <- cor(mz_h$L1, mz_h$L2)
  r_dz_full <- cor(dz_h$L1, dz_h$L2)

  km <- !mz_h$d1 & !mz_h$d2
  kd <- !dz_h$d1 & !dz_h$d2
  r_mz_surv <- cor(mz_h$L1[km], mz_h$L2[km])
  r_dz_surv <- cor(dz_h$L1[kd], dz_h$L2[kd])

  list(
    r_mz = r_mz_surv, r_dz = r_dz_surv,
    h2 = 2 * (r_mz_surv - r_dz_surv),
    r_mz_full = r_mz_full, r_dz_full = r_dz_full,
    h2_full = 2 * (r_mz_full - r_dz_full),
    ext_death_mz = mean(c(mz_h$d1, mz_h$d2)),
    ext_death_dz = mean(c(dz_h$d1, dz_h$d2))
  )
}

#' Oracle fix arm: calibrate sigma_theta using correctly specified model
#'
#' Holds (sigma_gamma, rho) at their true DGP values during calibration,
#' then extrapolates to m_ex=0 to show the bias vanishes.
#'
#' @param sigma_theta_true True sigma_theta (oracle)
#' @param params Parameter list
#' @return List with fix_sigma_theta, fix_h2, fix_bias, fix_bias_pp
run_oracle_fix_arm <- function(sigma_theta_true, params) {
  # Step 1: Generate observed data from the misspecified DGP (same as Arm 2)
  obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
    sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = params$MASTER_SEED + 95000L)

  # Step 2: Calibrate sigma_theta using a correctly specified model
  # (holds sigma_gamma and rho at their true values)
  sigma_theta_fix <- calibrate_sigma_theta_fix(
    target_r_mz = obs$r_mz,
    m_ex = params$M_EX_HIST,
    sigma_gamma = params$SIGMA_GAMMA_DEF,
    rho = params$RHO_DEF)

  # Step 3: Extrapolate to m_ex=0 using true sigma_gamma=0
  # (in intrinsic-only world, extrinsic frailty vanishes too)
  fix_corr <- sim_twin_h2(sigma_theta_fix, 0,
    sigma_gamma = 0, rho = 0,
    rng_seed = params$MASTER_SEED + 96000L)

  oracle <- sim_twin_h2(sigma_theta_true, 0,
    sigma_gamma = 0, rho = 0,
    rng_seed = params$MASTER_SEED + 97000L)

  list(
    obs_r_mz = obs$r_mz,
    obs_r_dz = obs$r_dz,
    sigma_theta_fix = sigma_theta_fix,
    sigma_theta_true = sigma_theta_true,
    sigma_recovery_pct = 100 * (sigma_theta_fix - sigma_theta_true) / sigma_theta_true,
    fix_h2 = fix_corr$h2,
    fix_r_mz = fix_corr$r_mz,
    fix_r_dz = fix_corr$r_dz,
    oracle_h2 = oracle$h2,
    fix_bias = fix_corr$h2 - oracle$h2,
    fix_bias_pp = round(100 * (fix_corr$h2 - oracle$h2), 1)
  )
}

#' Run full model analysis (oracle + arm1 + arm2) for MGG or SR
#'
#' @param model_name "mgg" or "sr"
#' @param params Parameter list
#' @return List with oracle_h2, arm1_h2, arm1_bias, arm2_h2, arm2_bias, etc.
run_model_analysis <- function(model_name, params) {
  if (model_name == "mgg") {
    sim_fn <- sim_twin_h2_mgg
    true_disp <- params$MGG_CV_Q     # sigma_q = 0.27
    disp_name <- "sigma_q"
  } else {
    sim_fn <- sim_twin_h2_sr
    true_disp <- params$SR_CV_XC * params$SR_MU_XC  # sigma_Xc = 3.57
    disp_name <- "sigma_Xc"
  }

  n_pairs_eff <- params$N_PAIRS

  seed_base <- params$MASTER_SEED +
    if (model_name == "mgg") 20000L else 30000L

  # Oracle (m_ex = 0, true dispersion)
  if (model_name == "mgg") {
    oracle <- sim_fn(sigma_q = true_disp, m_ex = 0,
                     n_pairs = n_pairs_eff, rng_seed = seed_base)
  } else {
    oracle <- sim_fn(sigma_Xc = true_disp, m_ex = 0,
                     n_pairs = n_pairs_eff, rng_seed = seed_base)
  }

  # Arm 1: correct specification
  if (model_name == "mgg") {
    obs_arm1 <- sim_fn(sigma_q = true_disp, m_ex = params$M_EX_HIST,
                       n_pairs = n_pairs_eff, rng_seed = seed_base + 1000L)
  } else {
    obs_arm1 <- sim_fn(sigma_Xc = true_disp, m_ex = params$M_EX_HIST,
                       n_pairs = n_pairs_eff, rng_seed = seed_base + 1000L)
  }
  calib_arm1 <- calibrate_dispersion(obs_arm1$r_mz, model = model_name,
                                      m_ex = params$M_EX_HIST,
                                      n_pairs = n_pairs_eff)
  if (model_name == "mgg") {
    extrap_arm1 <- sim_fn(sigma_q = calib_arm1, m_ex = 0,
                          n_pairs = n_pairs_eff, rng_seed = seed_base + 2000L)
  } else {
    extrap_arm1 <- sim_fn(sigma_Xc = calib_arm1, m_ex = 0,
                          n_pairs = n_pairs_eff, rng_seed = seed_base + 2000L)
  }

  # Arm 2: misspecified (heritable gamma)
  if (model_name == "mgg") {
    obs_arm2 <- sim_fn(sigma_q = true_disp, m_ex = params$M_EX_HIST,
                       sigma_gamma = params$SIGMA_GAMMA_DEF,
                       rho = params$RHO_DEF, individual_ext = TRUE,
                       n_pairs = n_pairs_eff, rng_seed = seed_base + 3000L)
  } else {
    obs_arm2 <- sim_fn(sigma_Xc = true_disp, m_ex = params$M_EX_HIST,
                       sigma_gamma = params$SIGMA_GAMMA_DEF,
                       rho = params$RHO_DEF, individual_ext = TRUE,
                       n_pairs = n_pairs_eff, rng_seed = seed_base + 3000L)
  }
  calib_arm2 <- calibrate_dispersion(obs_arm2$r_mz, model = model_name,
                                      m_ex = params$M_EX_HIST,
                                      n_pairs = n_pairs_eff)
  if (model_name == "mgg") {
    extrap_arm2 <- sim_fn(sigma_q = calib_arm2, m_ex = 0,
                          n_pairs = n_pairs_eff, rng_seed = seed_base + 4000L)
  } else {
    extrap_arm2 <- sim_fn(sigma_Xc = calib_arm2, m_ex = 0,
                          n_pairs = n_pairs_eff, rng_seed = seed_base + 4000L)
  }

  list(
    model = model_name,
    oracle_h2 = oracle$h2,
    oracle_r_mz = oracle$r_mz,
    oracle_r_dz = oracle$r_dz,
    arm1_h2 = extrap_arm1$h2,
    arm1_bias = extrap_arm1$h2 - oracle$h2,
    arm1_disp = calib_arm1,
    arm2_h2 = extrap_arm2$h2,
    arm2_bias = extrap_arm2$h2 - oracle$h2,
    arm2_disp = calib_arm2,
    true_disp = true_disp
  )
}

#' Dose-response sweep: bias as a function of m_ex
#'
#' Sweeps m_ex from 0 to 0.012 with heritable gamma.
#' Shows bias scales with extrinsic mortality magnitude.
#'
#' @param sigma_theta_true True sigma_theta
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return Data frame with m_ex, sigma_fit, h2_corrected, bias, sigma_infl_pct
run_dose_response <- function(sigma_theta_true, oracle, params,
                               n_reps = params$N_REPS_SWEEP) {
  m_ex_grid <- c(0, 0.00005, 0.0001, 0.0003, 0.0005, 0.0008,
                  0.001, 0.0015, 0.002, 0.003, 0.004, 0.005,
                  0.006, 0.008, 0.012)

  summary_rows <- list()
  per_seed_rows <- list()

  for (i in seq_along(m_ex_grid)) {
    mx <- m_ex_grid[i]

    if (mx == 0) {
      summary_rows[[i]] <- data.frame(
        m_ex = 0, bias_pp = 0, se_pp = 0,
        lo95_pp = 0, hi95_pp = 0,
        stringsAsFactors = FALSE
      )
      per_seed_rows[[i]] <- data.frame(
        m_ex = 0, rep = seq_len(n_reps), bias_pp = 0,
        stringsAsFactors = FALSE
      )
      next
    }

    seed_base <- params$MASTER_SEED + 70000L + i * 100L

    rep_biases <- vapply(seq_len(n_reps), function(r) {
      sb <- seed_base + r * 100000L
      obs_i <- sim_twin_h2(sigma_theta_true, mx,
        sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
        gamma_heritable = TRUE, individual_ext = TRUE,
        rng_seed = sb)
      sf <- calibrate_sigma_theta(obs_i$r_mz, mx)
      corr_i <- sim_twin_h2(sf, 0, sigma_gamma = 0,
        rng_seed = sb + 5000L)
      100 * (corr_i$h2 - oracle$h2)
    }, numeric(1))

    m <- mean(rep_biases); se <- sd(rep_biases) / sqrt(n_reps)
    t_crit <- qt(0.975, df = n_reps - 1)
    summary_rows[[i]] <- data.frame(
      m_ex = mx,
      bias_pp = m,
      se_pp = se,
      lo95_pp = m - t_crit * se,
      hi95_pp = m + t_crit * se,
      stringsAsFactors = FALSE
    )
    per_seed_rows[[i]] <- data.frame(
      m_ex = mx, rep = seq_len(n_reps), bias_pp = rep_biases,
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, summary_rows),
    per_seed = do.call(rbind, per_seed_rows)
  )
}

#' Pleiotropy isolation control (rho=0, heritable gamma)
#'
#' Tests whether bias requires genetic correlation between intrinsic
#' and extrinsic susceptibility.
#'
#' @param sigma_theta_true True sigma_theta
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return List with obs, corrected, bias, and comparison to default rho
run_pleiotropy_isolation <- function(sigma_theta_true, oracle, params) {
  seed_base <- params$MASTER_SEED + 80000L

  obs_rho0 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
    sigma_gamma = params$SIGMA_GAMMA_DEF, rho = 0,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = seed_base
  )

  sf_rho0 <- calibrate_sigma_theta(obs_rho0$r_mz, params$M_EX_HIST)
  corr_rho0 <- sim_twin_h2(sf_rho0, 0, sigma_gamma = 0,
    rng_seed = seed_base + 5000L
  )

  bias_rho0 <- corr_rho0$h2 - oracle$h2

  list(
    obs_r_mz = obs_rho0$r_mz,
    obs_r_dz = obs_rho0$r_dz,
    obs_h2 = obs_rho0$h2,
    sigma_fit = sf_rho0,
    sigma_infl_pct = 100 * (sf_rho0 / sigma_theta_true - 1),
    h2_corrected = corr_rho0$h2,
    bias = bias_rho0,
    bias_pp = 100 * bias_rho0,
    pleiotropy_required = abs(bias_rho0) < 0.02
  )
}

#' Irrelevant heritable trait control
#'
#' Creates heritable zeta correlated with theta but NOT in hazard.
#' Proves bias requires a mortality-relevant pathway.
#'
#' @param sigma_theta_true True sigma_theta
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return List with obs, corrected, bias (should be ~0)
run_irrelevant_trait <- function(sigma_theta_true, oracle, params) {
  seed_base <- params$MASTER_SEED + 90000L

  obs_zeta <- sim_twin_h2_with_zeta(
    sigma_theta_true, params$M_EX_HIST,
    sigma_zeta = params$SIGMA_GAMMA_DEF,
    rho_tz = params$RHO_DEF,
    n = params$N_PAIRS,
    rng_seed = seed_base
  )

  sf_zeta <- calibrate_sigma_theta(obs_zeta$r_mz, params$M_EX_HIST)
  corr_zeta <- sim_twin_h2(sf_zeta, 0, sigma_gamma = 0,
    rng_seed = seed_base + 5000L
  )

  bias_zeta <- corr_zeta$h2 - oracle$h2

  list(
    obs_r_mz = obs_zeta$r_mz,
    obs_r_dz = obs_zeta$r_dz,
    obs_h2 = obs_zeta$h2,
    sigma_fit = sf_zeta,
    sigma_infl_pct = 100 * (sf_zeta / sigma_theta_true - 1),
    h2_corrected = corr_zeta$h2,
    bias = bias_zeta,
    bias_pp = 100 * bias_zeta,
    trait_irrelevant_confirmed = abs(bias_zeta) < 0.02
  )
}

#' Simulate twin h² with irrelevant heritable zeta (not in hazard)
#'
#' Zeta is correlated with theta genetically but does NOT enter
#' the mortality hazard function.
#'
#' @param sigma_theta SD of log-Gompertz intercept
#' @param m_ex Extrinsic hazard rate
#' @param sigma_zeta SD of irrelevant trait
#' @param rho_tz Genetic correlation between theta and zeta
#' @param n Number of twin pairs
#' @param rng_seed RNG seed
#' @return List with r_mz, r_dz, h2
# ===================================================================
# MGG/SR extensions: dose-response, pleiotropy isolation, negative rho
# ===================================================================

#' Generic control sweep for MGG or SR models
#'
#' Runs observe → calibrate → extrapolate cycle for a sweep parameter.
#'
#' @param model_name "mgg" or "sr"
#' @param sweep_type "negative_rho", "pleiotropy_isolation", or "dose_response"
#' @param params Parameter list
#' @return Data frame with sweep results
run_model_control_sweep <- function(model_name, sweep_type, params,
                                    n_reps = params$N_REPS_SWEEP) {
  if (model_name == "mgg") {
    sim_fn <- sim_twin_h2_mgg
    true_disp <- params$MGG_CV_Q
    disp_name <- "sigma_q"
  } else {
    sim_fn <- sim_twin_h2_sr
    true_disp <- params$SR_CV_XC * params$SR_MU_XC
    disp_name <- "sigma_Xc"
  }

  n_pairs_eff <- params$N_PAIRS

  seed_offset <- if (model_name == "mgg") 40000L else 45000L
  seed_base <- params$MASTER_SEED + seed_offset

  # Oracle (m_ex = 0, true dispersion)
  oracle_args <- list(
    m_ex = 0, sigma_gamma = 0, rho = 0,
    n_pairs = n_pairs_eff, rng_seed = seed_base
  )
  oracle_args[[disp_name]] <- true_disp
  oracle <- do.call(sim_fn, oracle_args)

  # Build sweep grid
  if (sweep_type == "negative_rho") {
    grid <- data.frame(
      rho_val = seq(-0.6, 0.8, by = 0.05),
      m_ex_val = params$M_EX_HIST,
      sg_val = params$SIGMA_GAMMA_DEF,
      stringsAsFactors = FALSE
    )
  } else if (sweep_type == "pleiotropy_isolation") {
    grid <- data.frame(
      rho_val = 0,
      m_ex_val = params$M_EX_HIST,
      sg_val = params$SIGMA_GAMMA_DEF,
      stringsAsFactors = FALSE
    )
  } else if (sweep_type == "dose_response") {
    m_grid <- c(0, 0.00005, 0.0001, 0.0003, 0.0005, 0.0008,
                0.001, 0.0015, 0.002, 0.003, 0.004, 0.005,
                0.006, 0.008, 0.012)
    grid <- data.frame(
      rho_val = params$RHO_DEF,
      m_ex_val = m_grid,
      sg_val = params$SIGMA_GAMMA_DEF,
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown sweep_type: ", sweep_type)
  }

  # Pre-compile SR Rcpp before forking workers (avoids parallel build race)
  if (model_name == "sr") .ensure_sr_cpp()

  # Use parallel::mclapply for SR (Euler-Maruyama is expensive per grid point)
  # Cap at INNER_MC_CORES to avoid oversubscription under {crew} workers
  mc_cores <- if (model_name == "sr" && nrow(grid) > 1) {
    inner_cap <- params$INNER_MC_CORES %||% 1L
    min(nrow(grid), inner_cap)
  } else {
    1L
  }

  run_one_point <- function(i) {
    rho_i <- grid$rho_val[i]
    m_ex_i <- grid$m_ex_val[i]
    sg_i <- grid$sg_val[i]
    sweep_val <- if (sweep_type == "dose_response") m_ex_i else rho_i

    # If m_ex = 0, bias is 0 by definition
    if (m_ex_i == 0) {
      return(data.frame(
        model = model_name, sweep_type = sweep_type,
        sweep_val = sweep_val,
        oracle_h2 = oracle$h2, extrap_h2 = oracle$h2,
        bias = 0, bias_pp = 0, se_pp = 0,
        lo95_pp = 0, hi95_pp = 0,
        stringsAsFactors = FALSE
      ))
    }

    # Run n_reps independent seeds per grid point
    rep_biases <- vapply(seq_len(n_reps), function(r) {
      seed_i <- seed_base + 1000L * i + r * 100000L

      # 1. Simulate observed
      obs_args <- list(
        m_ex = m_ex_i, sigma_gamma = sg_i, rho = rho_i,
        individual_ext = TRUE,
        n_pairs = n_pairs_eff, rng_seed = seed_i
      )
      obs_args[[disp_name]] <- true_disp
      obs <- do.call(sim_fn, obs_args)

      # 2. Calibrate dispersion to match observed r_MZ
      calib_disp <- calibrate_dispersion(obs$r_mz, model = model_name,
                                          m_ex = m_ex_i, n_pairs = n_pairs_eff)

      # 3. Extrapolate to m_ex=0
      extrap_args <- list(
        m_ex = 0, sigma_gamma = 0,
        n_pairs = n_pairs_eff, rng_seed = seed_i + 500L
      )
      extrap_args[[disp_name]] <- calib_disp
      extrap <- do.call(sim_fn, extrap_args)

      100 * (extrap$h2 - oracle$h2)
    }, numeric(1))

    m <- mean(rep_biases)
    se <- if (n_reps > 1) sd(rep_biases) / sqrt(n_reps) else NA_real_
    t_crit <- if (n_reps > 1) qt(0.975, df = n_reps - 1) else NA_real_

    data.frame(
      model = model_name, sweep_type = sweep_type,
      sweep_val = sweep_val,
      oracle_h2 = oracle$h2, extrap_h2 = oracle$h2 + m / 100,
      bias = m / 100, bias_pp = m,
      se_pp = se,
      lo95_pp = m - t_crit * se,
      hi95_pp = m + t_crit * se,
      stringsAsFactors = FALSE
    )
  }

  if (mc_cores > 1L) {
    message(sprintf("  [%s/%s] Running %d grid points on %d cores",
                    model_name, sweep_type, nrow(grid), mc_cores))
    results <- parallel::mclapply(seq_len(nrow(grid)), run_one_point,
                                   mc.cores = mc_cores)
  } else {
    results <- lapply(seq_len(nrow(grid)), run_one_point)
  }

  do.call(rbind, results)
}

#' Multi-target calibration arm: bias persists under richer fitting criterion
#'
#' Generates data from the true DGP (heritable gamma), computes summary
#' statistics (r_MZ + marginal moments/quantiles), fits the misspecified
#' model to all five targets, and extrapolates to m_ex=0.
#'
#' @param sigma_theta_true True sigma_theta (oracle)
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return List with observed targets, multi-target fit, r_MZ-only comparison
run_multi_target_arm <- function(sigma_theta_true, oracle, params) {
  seed_base <- params$MASTER_SEED + 98000L

  # 1. Generate "observed" data from the true DGP (heritable gamma)
  set.seed(seed_base)
  n <- params$N_PAIRS
  pars_mz <- gen_twin_params_gm(n, r_g = 1.0, sigma_theta_true,
    params$SIGMA_GAMMA_DEF, params$RHO_DEF,
    gamma_heritable = TRUE)
  shift <- params$SIGMA_GAMMA_DEF^2 / 2
  c_mz1 <- params$M_EX_HIST * exp(pars_mz$gamma1 - shift)
  c_mz2 <- params$M_EX_HIST * exp(pars_mz$gamma2 - shift)
  L_mz1 <- sim_lifespan_gm(pars_mz$theta1, params$B_GOMP, c_mz1, params$T_MAX)
  L_mz2 <- sim_lifespan_gm(pars_mz$theta2, params$B_GOMP, c_mz2, params$T_MAX)
  keep <- (L_mz1 > params$CUTOFF_AGE) & (L_mz2 > params$CUTOFF_AGE)
  r_mz_obs <- cor(L_mz1[keep], L_mz2[keep])
  L_all_obs <- c(L_mz1[keep], L_mz2[keep])

  obs_targets <- list(
    r_mz = r_mz_obs,
    mean_age = mean(L_all_obs),
    sd_age = sd(L_all_obs),
    q25 = unname(quantile(L_all_obs, 0.25)),
    q75 = unname(quantile(L_all_obs, 0.75))
  )

  # 2. Fit misspecified model to all 5 targets
  mt_fit <- calibrate_multi_target(obs_targets, params$M_EX_HIST,
    rng_seed = seed_base + 1000L, n = n)

  # 3. Fit to r_MZ only (for comparison)
  sigma_rmz_only <- calibrate_sigma_theta_crn(r_mz_obs, params$M_EX_HIST,
    crn_seed = seed_base + 2000L)

  # 4. Extrapolate both to m_ex=0 (CRN: same seed for fair comparison)
  extrap_seed <- seed_base + 3000L
  extrap_mt <- sim_twin_h2(mt_fit$sigma_fit, 0, sigma_gamma = 0,
    rng_seed = extrap_seed)
  extrap_rmz <- sim_twin_h2(sigma_rmz_only, 0, sigma_gamma = 0,
    rng_seed = extrap_seed)

  list(
    obs_targets = obs_targets,
    mt_sigma_fit = mt_fit$sigma_fit,
    mt_h2 = extrap_mt$h2,
    mt_bias = extrap_mt$h2 - oracle$h2,
    mt_bias_pp = round(100 * (extrap_mt$h2 - oracle$h2), 1),
    mt_fitted = mt_fit$fitted_targets,
    mt_residuals = mt_fit$residuals,
    mt_gof_rmse = mt_fit$gof_rmse,
    mt_loss = mt_fit$loss,
    rmz_sigma_fit = sigma_rmz_only,
    rmz_h2 = extrap_rmz$h2,
    rmz_bias = extrap_rmz$h2 - oracle$h2,
    rmz_bias_pp = round(100 * (extrap_rmz$h2 - oracle$h2), 1),
    sigma_theta_true = sigma_theta_true,
    mt_sigma_infl_pct = 100 * (mt_fit$sigma_fit / sigma_theta_true - 1),
    rmz_sigma_infl_pct = 100 * (sigma_rmz_only / sigma_theta_true - 1)
  )
}

sim_twin_h2_with_zeta <- function(sigma_theta, m_ex, sigma_zeta, rho_tz,
                                   n = PARAMS$N_PAIRS,
                                   rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)

  b <- PARAMS$B_GOMP
  t_max <- PARAMS$T_MAX
  cutoff <- PARAMS$CUTOFF_AGE
  mu_theta <- PARAMS$MU_THETA

  Sigma_within <- matrix(c(1, rho_tz, rho_tz, 1), 2, 2)

  # --- MZ twins (r_g = 1) ---
  G1_mz <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_within)
  theta1_mz <- mu_theta + sigma_theta * G1_mz[, 1]
  theta2_mz <- theta1_mz  # MZ: identical genotype

  L_mz1 <- sim_lifespan_gm(theta1_mz, b, m_ex, t_max)
  L_mz2 <- sim_lifespan_gm(theta2_mz, b, m_ex, t_max)

  # --- DZ twins (r_g = 0.5) ---
  G1_dz <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_within)
  Z_dz  <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_within)
  g_theta1_dz <- G1_dz[, 1]
  g_theta2_dz <- 0.5 * g_theta1_dz + sqrt(0.75) * Z_dz[, 1]

  theta1_dz <- mu_theta + sigma_theta * g_theta1_dz
  theta2_dz <- mu_theta + sigma_theta * g_theta2_dz

  L_dz1 <- sim_lifespan_gm(theta1_dz, b, m_ex, t_max)
  L_dz2 <- sim_lifespan_gm(theta2_dz, b, m_ex, t_max)

  keep_mz <- (L_mz1 > cutoff) & (L_mz2 > cutoff)
  keep_dz <- (L_dz1 > cutoff) & (L_dz2 > cutoff)

  r_mz <- cor(L_mz1[keep_mz], L_mz2[keep_mz])
  r_dz <- cor(L_dz1[keep_dz], L_dz2[keep_dz])

  list(r_mz = r_mz, r_dz = r_dz, h2 = 2 * (r_mz - r_dz))
}


# ===========================================================================
# B1: Alternative functional forms for extrinsic frailty DGP
# ===========================================================================

#' Run Arm 2 under alternative extrinsic frailty functional forms
#'
#' Tests whether the bias mechanism is robust to different DGP specifications
#' for how extrinsic frailty enters the hazard:
#'   - "lognormal" (default): c_i = m_ex * exp(gamma_i) (multiplicative log-normal)
#'   - "additive": c_i = m_ex + sigma_gamma * gamma_i (additive Gaussian, truncated at 0)
#'   - "gamma": c_i = m_ex * g_i, g_i ~ Gamma(shape=1/sigma^2, rate=1/sigma^2)
#'
#' @param sigma_theta_true True sigma_theta (oracle)
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return Data frame with one row per functional form, showing bias
run_alt_ext_forms <- function(sigma_theta_true, oracle, params, n_reps = 10L) {
  forms <- c("lognormal", "additive", "gamma")
  seed_base <- params$MASTER_SEED + 110000L

  results <- lapply(seq_along(forms), function(i) {
    form <- forms[i]

    rep_results <- lapply(seq_len(n_reps), function(r) {
      seed_ir <- seed_base + i * 1000L + r * 100L

      # Step 1: Observe under true DGP with this functional form
      obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
        sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
        gamma_heritable = TRUE, individual_ext = TRUE,
        rng_seed = seed_ir, ext_form = form)

      # Step 2: Calibrate misspecified model using CRN for stability
      sigma_fit <- calibrate_sigma_theta_crn(obs$r_mz, params$M_EX_HIST,
        crn_seed = seed_ir + 500L)

      # Step 3: Extrapolate to m_ex=0
      corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
        rng_seed = seed_ir + 900L)

      data.frame(
        ext_form = form,
        rep = r,
        obs_r_mz = obs$r_mz,
        sigma_fit = sigma_fit,
        sigma_infl_pct = 100 * (sigma_fit / sigma_theta_true - 1),
        h2 = corr$h2,
        bias = corr$h2 - oracle$h2,
        bias_pp = 100 * (corr$h2 - oracle$h2),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rep_results)
  })

  do.call(rbind, results)
}


# ===========================================================================
# B2: Extended sigma_gamma range
# ===========================================================================

#' Run bias sweep at high sigma_gamma values (above anchored ceiling)
#'
#' Explores bias behavior at sigma_gamma between 0.65 and 1.47,
#' documenting how the bias scales above the anchored regime.
#'
#' @param sigma_theta_true True sigma_theta
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @return Data frame with one row per sigma_gamma value
run_extended_sigma_gamma <- function(sigma_theta_true, oracle, params,
                                     n_reps = 10L) {
  sg_grid <- c(0.70, 0.80, 0.90, 1.00, 1.20, 1.47)
  rho_mid <- 0.35  # midpoint of anchored range
  seed_base <- params$MASTER_SEED + 120000L

  results <- lapply(seq_along(sg_grid), function(i) {
    sg <- sg_grid[i]

    rep_results <- lapply(seq_len(n_reps), function(r) {
      seed_ir <- seed_base + i * 1000L + r * 100L

      obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
        sigma_gamma = sg, rho = rho_mid,
        gamma_heritable = TRUE, individual_ext = TRUE,
        rng_seed = seed_ir)

      # Use CRN calibration for stability at extreme parameters
      sigma_fit <- calibrate_sigma_theta_crn(obs$r_mz, params$M_EX_HIST,
        crn_seed = seed_ir + 500L)

      corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
        rng_seed = seed_ir + 900L)

      data.frame(
        sigma_gamma = sg,
        rho = rho_mid,
        rep = r,
        obs_r_mz = obs$r_mz,
        sigma_fit = sigma_fit,
        sigma_infl_pct = 100 * (sigma_fit / sigma_theta_true - 1),
        h2 = corr$h2,
        bias = corr$h2 - oracle$h2,
        bias_pp = 100 * (corr$h2 - oracle$h2),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rep_results)
  })

  do.call(rbind, results)
}


# ===========================================================================
# B3: Monte Carlo uncertainty bounds
# ===========================================================================

#' Run primary arms across multiple seeds for MC uncertainty estimation
#'
#' Runs Oracle, Arm 1, and Arm 2 at n_seeds independent seeds.
#' Reports MC SE on h², r_MZ, bias.
#'
#' @param sigma_theta_true True sigma_theta
#' @param params Parameter list
#' @param n_seeds Number of independent replications (default: 20)
#' @return List with per-seed results and summary statistics
run_mc_uncertainty <- function(sigma_theta_true, params, n_seeds = 20L) {
  seed_base <- params$MASTER_SEED + 130000L

  results <- lapply(seq_len(n_seeds), function(s) {
    seed_s <- seed_base + s * 10000L

    # Oracle
    ora <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0,
      rng_seed = seed_s)

    # Arm 1: correctly specified
    obs1 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = 0, rng_seed = seed_s + 100L)
    sf1 <- calibrate_sigma_theta_crn(obs1$r_mz, params$M_EX_HIST,
      crn_seed = seed_s + 150L)
    corr1 <- sim_twin_h2(sf1, 0, sigma_gamma = 0,
      rng_seed = seed_s + 200L)

    # Arm 2: misspecified
    obs2 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = seed_s + 300L)
    sf2 <- calibrate_sigma_theta_crn(obs2$r_mz, params$M_EX_HIST,
      crn_seed = seed_s + 350L)
    corr2 <- sim_twin_h2(sf2, 0, sigma_gamma = 0,
      rng_seed = seed_s + 400L)

    data.frame(
      seed = s,
      oracle_h2 = ora$h2,
      oracle_r_mz = ora$r_mz,
      arm1_h2 = corr1$h2,
      arm1_bias_pp = 100 * (corr1$h2 - ora$h2),
      arm2_h2 = corr2$h2,
      arm2_bias_pp = 100 * (corr2$h2 - ora$h2),
      arm2_sigma_fit = sf2,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, results)

  list(
    per_seed = df,
    summary = data.frame(
      statistic = c("oracle_h2", "arm1_bias_pp", "arm2_bias_pp",
                     "arm2_h2", "arm2_sigma_fit"),
      mean = c(mean(df$oracle_h2), mean(df$arm1_bias_pp),
               mean(df$arm2_bias_pp), mean(df$arm2_h2),
               mean(df$arm2_sigma_fit)),
      sd = c(sd(df$oracle_h2), sd(df$arm1_bias_pp),
             sd(df$arm2_bias_pp), sd(df$arm2_h2),
             sd(df$arm2_sigma_fit)),
      se = c(sd(df$oracle_h2), sd(df$arm1_bias_pp),
             sd(df$arm2_bias_pp), sd(df$arm2_h2),
             sd(df$arm2_sigma_fit)) / sqrt(nrow(df)),
      stringsAsFactors = FALSE
    )
  )
}


# ===========================================================================
# B3b: Multi-model MC uncertainty (GM + MGG + SR)
# ===========================================================================

#' Run primary arms across multiple seeds for all three mortality models
#'
#' Extends B3 to MGG and SR. For each seed, runs Oracle/Arm1/Arm2 under
#' each model and records bias. Returns long-format data frame with model column.
#'
#' @param sigma_theta_true True sigma_theta (for GM)
#' @param params Parameter list
#' @param n_seeds Number of independent replications (default: 20)
#' @return List with per_seed data frame and summary data frame
#' Run MC uncertainty for a single mortality model
#'
#' @param model_name One of "gm", "mgg", "sr"
#' @param sigma_theta_true True sigma_theta (used for GM)
#' @param params Parameter list
#' @param n_seeds Number of independent replications
#' @return Data frame with one row per seed
run_mc_uncertainty_model <- function(model_name, sigma_theta_true, params,
                                     n_seeds = 20L) {
  seed_base <- params$MASTER_SEED + 150000L
  n_pairs <- params$N_PAIRS
  soff <- switch(model_name, gm = 0L, mgg = 1000000L, sr = 2000000L)

  mc_cores <- min(n_seeds, params$INNER_MC_CORES %||% 1L)
  if (mc_cores > 1L) {
    message(sprintf("  [mc_unc_%s] Running %d seeds on %d cores",
                    model_name, n_seeds, mc_cores))
  }

  # Pre-compile SR Rcpp before forking workers (avoids parallel build race)
  if (model_name == "sr") .ensure_sr_cpp()

  results <- parallel::mclapply(seq_len(n_seeds), mc.cores = mc_cores, function(s) {
    seed_s <- seed_base + soff + s * 10000L

    if (model_name == "gm") {
      ora <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0,
        rng_seed = seed_s)
      obs1 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
        sigma_gamma = 0, rng_seed = seed_s + 100L)
      sf1 <- calibrate_sigma_theta_crn(obs1$r_mz, params$M_EX_HIST,
        crn_seed = seed_s + 150L)
      corr1 <- sim_twin_h2(sf1, 0, sigma_gamma = 0,
        rng_seed = seed_s + 200L)
      obs2 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
        sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
        gamma_heritable = TRUE, individual_ext = TRUE,
        rng_seed = seed_s + 300L)
      sf2 <- calibrate_sigma_theta_crn(obs2$r_mz, params$M_EX_HIST,
        crn_seed = seed_s + 350L)
      corr2 <- sim_twin_h2(sf2, 0, sigma_gamma = 0,
        rng_seed = seed_s + 400L)
    } else {
      sim_fn <- if (model_name == "mgg") sim_twin_h2_mgg else sim_twin_h2_sr
      dname <- if (model_name == "mgg") "sigma_q" else "sigma_Xc"
      true_d <- if (model_name == "mgg") params$MGG_CV_Q
        else params$SR_CV_XC * params$SR_MU_XC

      mk_args <- function(disp, ...) {
        a <- list(...)
        a[[dname]] <- disp
        a
      }

      ora <- do.call(sim_fn, mk_args(true_d,
        m_ex = 0, n_pairs = n_pairs, rng_seed = seed_s))
      obs1 <- do.call(sim_fn, mk_args(true_d,
        m_ex = params$M_EX_HIST, n_pairs = n_pairs,
        rng_seed = seed_s + 100L))
      sf1 <- calibrate_dispersion(obs1$r_mz, model = model_name,
        m_ex = params$M_EX_HIST, n_pairs = n_pairs)
      corr1 <- do.call(sim_fn, mk_args(sf1,
        m_ex = 0, n_pairs = n_pairs, rng_seed = seed_s + 200L))
      obs2 <- do.call(sim_fn, mk_args(true_d,
        m_ex = params$M_EX_HIST,
        sigma_gamma = params$SIGMA_GAMMA_DEF,
        rho = params$RHO_DEF, individual_ext = TRUE,
        n_pairs = n_pairs, rng_seed = seed_s + 300L))
      sf2 <- calibrate_dispersion(obs2$r_mz, model = model_name,
        m_ex = params$M_EX_HIST, n_pairs = n_pairs)
      corr2 <- do.call(sim_fn, mk_args(sf2,
        m_ex = 0, n_pairs = n_pairs, rng_seed = seed_s + 400L))
    }

    data.frame(
      model = toupper(model_name),
      seed = s,
      oracle_h2 = ora$h2,
      arm1_h2 = corr1$h2,
      arm1_bias_pp = 100 * (corr1$h2 - ora$h2),
      arm2_h2 = corr2$h2,
      arm2_bias_pp = 100 * (corr2$h2 - ora$h2),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Combine per-model MC uncertainty results into the multimodel structure
#'
#' @param ... Data frames from run_mc_uncertainty_model() calls
#' @return List with per_seed data frame and summary data frame
combine_mc_uncertainty_multimodel <- function(...) {
  df <- do.call(rbind, list(...))

  summ <- do.call(rbind, lapply(split(df, df$model), function(d) {
    data.frame(
      model = d$model[1],
      n_seeds = nrow(d),
      oracle_h2_mean = mean(d$oracle_h2),
      oracle_h2_se = sd(d$oracle_h2) / sqrt(nrow(d)),
      arm1_bias_mean = mean(d$arm1_bias_pp),
      arm1_bias_se = sd(d$arm1_bias_pp) / sqrt(nrow(d)),
      arm2_bias_mean = mean(d$arm2_bias_pp),
      arm2_bias_se = sd(d$arm2_bias_pp) / sqrt(nrow(d)),
      stringsAsFactors = FALSE
    )
  }))

  list(per_seed = df, summary = summ)
}

# ===========================================================================
# B4: High-replication m_ex split
# ===========================================================================

#' Run m_ex split at higher replication to resolve MC noise
#'
#' Re-runs the heritable fraction sweep at n_reps replications
#' to clarify whether the 50% dip is real.
#'
#' @param sigma_theta_true True sigma_theta
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @param n_reps Number of replications per fraction level (default: 5)
#' @return Data frame with per-replication results (frac_heritable, rep, h2, bias, bias_pp)
run_mex_split_hires <- function(sigma_theta_true, oracle, params,
                                 n_reps = 5L) {
  fractions <- seq(0, 1, by = 0.05)
  seed_base <- params$MASTER_SEED + 140000L

  results <- lapply(seq_along(fractions), function(i) {
    frac <- fractions[i]

    rep_results <- lapply(seq_len(n_reps), function(r) {
      seed_ir <- seed_base + i * 1000L + r * 100L
      m_inf <- frac * params$M_EX_HIST
      m_other <- (1 - frac) * params$M_EX_HIST

      if (frac == 0) {
        # No heritable fraction: all extrinsic is population-constant
        obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
          sigma_gamma = 0, rng_seed = seed_ir)
      } else {
        # Split: c_i = m_other + m_inf * exp(gamma_i - sigma^2/2)
        obs <- sim_twin_h2_split_mex(sigma_theta_true, m_inf, m_other,
          sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
          rng_seed = seed_ir)
      }

      sf <- calibrate_sigma_theta_crn(obs$r_mz, params$M_EX_HIST,
        crn_seed = seed_ir + 500L)
      corr <- sim_twin_h2(sf, 0, sigma_gamma = 0,
        rng_seed = seed_ir + 50L)

      data.frame(
        frac_heritable = frac,
        rep = r,
        h2 = corr$h2,
        bias = corr$h2 - oracle$h2,
        bias_pp = 100 * (corr$h2 - oracle$h2),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rep_results)
  })

  do.call(rbind, results)
}
