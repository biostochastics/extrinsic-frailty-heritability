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
run_negative_rho_sensitivity <- function(sigma_theta_true, params,
                                          n_reps = params$N_REPS_SWEEP) {
  rho_grid <- seq(-0.6, 0.8, by = 0.05)
  oracle <- sim_twin_h2(sigma_theta_true, 0, sigma_gamma = 0,
                         rng_seed = params$MASTER_SEED)

  results <- lapply(seq_along(rho_grid), function(i) {
    rho_val <- rho_grid[i]
    seed_base <- params$MASTER_SEED + 60000L + round((rho_val + 1) * 1000)

    rep_biases <- vapply(seq_len(n_reps), function(r) {
      sb <- seed_base + r * 100000L
      obs <- sim_twin_h2(sigma_theta_true, params$M_EX_HIST,
        sigma_gamma = params$SIGMA_GAMMA_DEF, rho = rho_val,
        gamma_heritable = TRUE, individual_ext = TRUE,
        rng_seed = sb)
      sigma_fit <- calibrate_sigma_theta(obs$r_mz, params$M_EX_HIST)
      corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
        rng_seed = sb + 5000L)
      100 * (corr$h2 - oracle$h2)
    }, numeric(1))

    m <- mean(rep_biases); se <- sd(rep_biases) / sqrt(n_reps)
    t_crit <- qt(0.975, df = n_reps - 1)
    data.frame(
      rho = rho_val,
      bias_pp = m,
      se_pp = se,
      lo95_pp = m - t_crit * se,
      hi95_pp = m + t_crit * se,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

# ===========================================================================
# MZ extrinsic concordance sensitivity sweep
# ===========================================================================

#' Sweep within-MZ gamma concordance from 0 to 1
#'
#' Decomposes extrinsic frailty into a genetic component (shared at r_g)
#' and individual-specific exposure noise. gamma_concordance is the fraction
#' of gamma variance that is genetically shared, so
#' Corr(gamma_1, gamma_2) = gamma_concordance for MZ twins (0.5× for DZ).
#'
#' NOTE: Reducing gamma_concordance also reduces the effective within-person
#' pleiotropy Corr(theta, gamma_total) = sqrt(gamma_concordance) * rho,
#' because only the genetic component of gamma is pleiotropic with theta.
#' This is the physically correct model: if only a fraction of extrinsic
#' susceptibility is genetic, both twin concordance AND pleiotropy decrease.
#'
#' @param sigma_theta_true True sigma_theta from oracle calibration
#' @param oracle Oracle simulation result
#' @param params Parameter list
#' @param rho Pleiotropy correlation (default: params$RHO_DEF; set 0 for decoupled)
#' @param n_reps Number of replications per concordance level
#' @return Data frame with r_gamma_mz, bias_pp, se_pp columns
run_mz_concordance_sweep <- function(sigma_theta_true, oracle, params,
                                     rho = params$RHO_DEF,
                                     n_reps = params$N_REPS_SWEEP) {
  r_gamma_grid <- sort(unique(c(seq(0, 1, by = 0.1), 0.85)))

  results <- lapply(r_gamma_grid, function(r_gam) {
    biases <- vapply(seq_len(n_reps), function(rep) {
      seed <- params$MASTER_SEED + 90000L + round(r_gam * 1000) + rep * 7L

      obs <- sim_twin_h2_concordance(
        sigma_theta_true,
        m_ex = params$M_EX_HIST,
        sigma_gamma = params$SIGMA_GAMMA_DEF,
        rho = rho,
        gamma_concordance = r_gam,
        rng_seed = seed
      )

      sigma_fit <- calibrate_sigma_theta(obs$r_mz, params$M_EX_HIST)

      corr <- sim_twin_h2(sigma_fit, 0, sigma_gamma = 0,
        rng_seed = seed + 50000L
      )
      corr$h2 - oracle$h2
    }, numeric(1))

    data.frame(
      r_gamma_mz = r_gam,
      bias_mean = mean(biases),
      bias_se = sd(biases) / sqrt(n_reps),
      bias_pp = 100 * mean(biases),
      se_pp = 100 * sd(biases) / sqrt(n_reps),
      n_reps = n_reps,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

# ===========================================================================
# Bridge uncertainty quantification
# ===========================================================================

#' Run sigma_gamma bridge with MC uncertainty + Obel CI targets
#'
#' Wraps the bridge computation over multiple seeds and inverts at
#' Obel et al.'s central estimate (0.40) plus CI endpoints (0.12, 0.50).
#'
#' @param sigma_theta_true True sigma_theta from oracle calibration
#' @param params Parameter list
#' @param n_seeds Number of independent seeds
#' @return List with bridge_table and per-target summaries
run_bridge_uncertainty <- function(sigma_theta_true, params, n_seeds = 10L) {
  sigma_gamma_grid <- seq(0.1, 1.5, by = 0.025)
  obel_targets <- c(lower = 0.12, central = 0.40, upper = 0.50)

  seed_results <- lapply(seq_len(n_seeds), function(s) {
    base_seed <- params$MASTER_SEED + 70000L + s * 100000L

    h2_vals <- vapply(sigma_gamma_grid, function(sg) {
      compute_infection_death_h2(
        sigma_gamma = sg,
        rho = params$RHO_DEF,
        sigma_theta = sigma_theta_true,
        m_ex = params$M_EX_HIST,
        n_pairs = params$N_PAIRS,
        cutoff = params$CUTOFF_AGE,
        seed = base_seed + round(sg * 10000)
      )$h2_L
    }, numeric(1))

    # Invert via first sign-change bracket (robust to MC noise)
    bridges <- vapply(obel_targets, function(target) {
      diffs <- h2_vals - target
      sign_changes <- which(diff(sign(diffs)) != 0)
      if (length(sign_changes) == 0) return(NA_real_)
      # Use first crossing (lowest sigma_gamma)
      i <- sign_changes[1]
      # Linear interpolation within bracket
      w <- (target - h2_vals[i]) / (h2_vals[i + 1] - h2_vals[i])
      sigma_gamma_grid[i] + w * (sigma_gamma_grid[i + 1] - sigma_gamma_grid[i])
    }, numeric(1))

    data.frame(
      seed = s,
      target_h2_L = obel_targets,
      target_label = names(obel_targets),
      bridge_sigma_gamma = bridges,
      stringsAsFactors = FALSE
    )
  })

  all_bridges <- do.call(rbind, seed_results)

  summary_df <- do.call(rbind, lapply(names(obel_targets), function(lab) {
    vals <- all_bridges$bridge_sigma_gamma[all_bridges$target_label == lab]
    vals <- vals[!is.na(vals)]
    nv <- length(vals)
    se <- if (nv > 1) sd(vals) / sqrt(nv) else NA_real_
    data.frame(
      target_label = lab,
      target_h2_L = obel_targets[lab],
      bridge_mean = if (nv > 0) mean(vals) else NA_real_,
      bridge_se = se,
      bridge_lo95 = if (nv > 1) mean(vals) - qt(0.975, df = nv - 1) * se else NA_real_,
      bridge_hi95 = if (nv > 1) mean(vals) + qt(0.975, df = nv - 1) * se else NA_real_,
      n_valid = nv,
      all_above_065 = if (nv > 0) all(vals > 0.65) else NA,
      stringsAsFactors = FALSE
    )
  }))

  list(
    all_bridges = all_bridges,
    summary = summary_df,
    obel_targets = obel_targets,
    n_seeds = n_seeds
  )
}
