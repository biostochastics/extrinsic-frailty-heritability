# ===========================================================================
# replication.R — 50-seed MC replication of all GM arms
# ===========================================================================

#' Run all GM arms across n_reps independent seeds
#'
#' For each seed, runs the full observe -> calibrate -> extrapolate cycle
#' for each arm, recording h2 and bias. Returns per-seed data frame and
#' summary with mean/SE/CI.
#'
#' @param sigma_theta_true True sigma_theta (from oracle calibration)
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @param n_reps Number of seeds (default: params$N_REPS_MAIN)
#' @return List with per_seed data frame and summary data frame
run_main_arms_replicated <- function(sigma_theta_true, oracle, params,
                                      n_reps = params$N_REPS_MAIN) {
  seed_base <- params$MASTER_SEED + 200000L

  results <- lapply(seq_len(n_reps), function(s) {
    sb <- seed_base + s * 100000L

    # Oracle (re-run per seed for unbiased comparison)
    ora <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0, rng_seed = sb)

    # Arm 1: correctly specified (no gamma)
    obs1 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = 0, rng_seed = sb + 100L)
    sf1 <- calibrate_sigma_theta_crn(obs1$r_mz, params$M_EX_HIST,
      crn_seed = sb + 150L)
    corr1 <- sim_twin_h2(sf1, 0, sigma_gamma = 0, rng_seed = sb + 200L)

    # Arm 2: misspecified (heritable gamma)
    obs2 <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = sb + 300L)
    sf2 <- calibrate_sigma_theta_crn(obs2$r_mz, params$M_EX_HIST,
      crn_seed = sb + 350L)
    corr2 <- sim_twin_h2(sf2, 0, sigma_gamma = 0, rng_seed = sb + 400L)

    # Control A: non-heritable extrinsic
    obs_a <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = 0,
      gamma_heritable = FALSE, individual_ext = TRUE,
      rng_seed = sb + 500L)
    sf_a <- calibrate_sigma_theta_crn(obs_a$r_mz, params$M_EX_HIST,
      crn_seed = sb + 550L)
    corr_a <- sim_twin_h2(sf_a, 0, sigma_gamma = 0, rng_seed = sb + 600L)

    # Control B: vanishing m_ex
    obs_b <- sim_twin_h2(sigma_theta_true, params$M_EX_TINY,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = sb + 700L)
    sf_b <- calibrate_sigma_theta_crn(obs_b$r_mz, params$M_EX_TINY,
      crn_seed = sb + 750L)
    corr_b <- sim_twin_h2(sf_b, 0, sigma_gamma = 0, rng_seed = sb + 800L)

    # Oracle fix: correctly specified two-component calibration
    obs_fix <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = params$RHO_DEF,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = sb + 900L)
    sf_fix <- calibrate_sigma_theta_fix(
      target_r_mz = obs_fix$r_mz,
      m_ex = params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF,
      rho = params$RHO_DEF)
    corr_fix <- sim_twin_h2(sf_fix, 0, sigma_gamma = 0, rng_seed = sb + 1000L)

    # Pleiotropy isolation (rho=0, heritable gamma)
    obs_pi <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
      sigma_gamma = params$SIGMA_GAMMA_DEF, rho = 0,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = sb + 1100L)
    sf_pi <- calibrate_sigma_theta_crn(obs_pi$r_mz, params$M_EX_HIST,
      crn_seed = sb + 1150L)
    corr_pi <- sim_twin_h2(sf_pi, 0, sigma_gamma = 0, rng_seed = sb + 1200L)

    # Irrelevant trait (zeta correlated with theta but not in hazard)
    obs_irr <- sim_twin_h2_with_zeta(
      sigma_theta_true, params$M_EX_HIST,
      sigma_zeta = params$SIGMA_GAMMA_DEF,
      rho_tz = params$RHO_DEF,
      n = params$N_PAIRS,
      rng_seed = sb + 1300L)
    sf_irr <- calibrate_sigma_theta_crn(obs_irr$r_mz, params$M_EX_HIST,
      crn_seed = sb + 1350L)
    corr_irr <- sim_twin_h2(sf_irr, 0, sigma_gamma = 0, rng_seed = sb + 1400L)

    # Arm 3: Hamilton conditioning
    ham_h2 <- run_hamilton_arm_single(params, rng_seed = sb + 1500L)

    data.frame(
      seed = s,
      oracle_h2 = ora$h2,
      arm1_h2 = corr1$h2,
      arm1_bias_pp = 100 * (corr1$h2 - ora$h2),
      arm2_h2 = corr2$h2,
      arm2_bias_pp = 100 * (corr2$h2 - ora$h2),
      arm2_sigma_fit = sf2,
      arm2_sigma_infl = 100 * (sf2 / sigma_theta_true - 1),
      ctrl_a_h2 = corr_a$h2,
      ctrl_a_bias_pp = 100 * (corr_a$h2 - ora$h2),
      ctrl_b_h2 = corr_b$h2,
      ctrl_b_bias_pp = 100 * (corr_b$h2 - ora$h2),
      fix_h2 = corr_fix$h2,
      fix_bias_pp = 100 * (corr_fix$h2 - ora$h2),
      pleio_iso_h2 = corr_pi$h2,
      pleio_iso_bias_pp = 100 * (corr_pi$h2 - ora$h2),
      irrel_h2 = corr_irr$h2,
      irrel_bias_pp = 100 * (corr_irr$h2 - ora$h2),
      arm3_h2 = ham_h2$h2,
      arm3_bias_pp = 100 * (ham_h2$h2 - ora$h2),
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, results)

  # Summary: mean, SE, 95% CI for each metric (t-distribution, NA-safe)
  summarize_col <- function(col) {
    vals <- df[[col]]
    vals <- vals[is.finite(vals)]
    n_used <- length(vals)
    m <- mean(vals); se <- if (n_used > 1) sd(vals) / sqrt(n_used) else NA_real_
    t_crit <- if (n_used > 1) qt(0.975, df = n_used - 1) else NA_real_
    data.frame(metric = col, mean = m, se = se,
               lo95 = m - t_crit * se, hi95 = m + t_crit * se,
               n_used = n_used,
               stringsAsFactors = FALSE)
  }

  cols <- c("oracle_h2", "arm1_h2", "arm1_bias_pp",
            "arm2_h2", "arm2_bias_pp", "arm2_sigma_fit", "arm2_sigma_infl",
            "ctrl_a_h2", "ctrl_a_bias_pp",
            "ctrl_b_h2", "ctrl_b_bias_pp",
            "fix_h2", "fix_bias_pp",
            "pleio_iso_h2", "pleio_iso_bias_pp",
            "irrel_h2", "irrel_bias_pp",
            "arm3_h2", "arm3_bias_pp")

  summary <- do.call(rbind, lapply(cols, summarize_col))

  list(per_seed = df, summary = summary)
}

#' Run a single Hamilton arm simulation (for replication loop)
#'
#' Self-contained Hamilton conditioning arm with explicit seed.
#'
#' @param params Parameter list
#' @param rng_seed RNG seed
#' @return List with h2, r_mz, r_dz
run_hamilton_arm_single <- function(params, rng_seed) {
  set.seed(rng_seed)
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
      if (sum(k_mz) < 3 || sum(k_dz) < 3) next
      h2c <- 2 * (cor(mz$L1[k_mz], mz$L2[k_mz]) -
                   cor(dz$L1[k_dz], dz$L2[k_dz]))
      if (!is.finite(h2c)) next
      if (h2c < target_h2) lo <- mid else hi <- mid
    }
    (lo + hi) / 2
  }

  ham_g <- ham_calibrate(params$HAM_TARGET_H2)
  mz_h <- hamilton_sim(n, 1.0, ham_g)
  dz_h <- hamilton_sim(n, 0.5, ham_g)

  km <- !mz_h$d1 & !mz_h$d2
  kd <- !dz_h$d1 & !dz_h$d2
  r_mz <- if (sum(km) >= 3) cor(mz_h$L1[km], mz_h$L2[km]) else NA_real_
  r_dz <- if (sum(kd) >= 3) cor(dz_h$L1[kd], dz_h$L2[kd]) else NA_real_

  list(r_mz = r_mz, r_dz = r_dz, h2 = 2 * (r_mz - r_dz))
}
