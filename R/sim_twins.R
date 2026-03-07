# ===========================================================================
# sim_twins.R — Twin parameter generation and Falconer h² computation
# ===========================================================================

#' Generate correlated twin parameters (theta, gamma) for GM model
#'
#' @param n Number of twin pairs
#' @param r_g Genetic correlation (1.0 for MZ, 0.5 for DZ)
#' @param sigma_theta SD of log-Gompertz intercept
#' @param sigma_gamma SD of log-extrinsic frailty (0 = no extrinsic heterogeneity)
#' @param rho Genetic correlation between theta and gamma (pleiotropy)
#' @param gamma_heritable Whether gamma is shared between twins
#' @param mu_theta Mean log-Gompertz intercept
#' @return List with theta1, theta2, gamma1, gamma2
gen_twin_params_gm <- function(n, r_g, sigma_theta, sigma_gamma, rho,
                               gamma_heritable = TRUE,
                               mu_theta = PARAMS$MU_THETA) {
  if (sigma_gamma == 0 || !gamma_heritable) {
    z1 <- rnorm(n)
    z2 <- rnorm(n)
    g_theta1 <- z1
    g_theta2 <- r_g * z1 + sqrt(1 - r_g^2) * z2

    theta1 <- mu_theta + sigma_theta * g_theta1
    theta2 <- mu_theta + sigma_theta * g_theta2

    if (sigma_gamma > 0 && !gamma_heritable) {
      gamma1 <- sigma_gamma * rnorm(n)
      gamma2 <- sigma_gamma * rnorm(n)
    } else {
      gamma1 <- rep(0, n)
      gamma2 <- rep(0, n)
    }
  } else {
    Sigma_within <- matrix(c(1, rho, rho, 1), 2, 2)
    G1 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_within)
    g_theta1 <- G1[, 1]
    g_gamma1 <- G1[, 2]

    Z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma_within)
    g_theta2 <- r_g * g_theta1 + sqrt(1 - r_g^2) * Z[, 1]
    g_gamma2 <- r_g * g_gamma1 + sqrt(1 - r_g^2) * Z[, 2]

    theta1 <- mu_theta + sigma_theta * g_theta1
    theta2 <- mu_theta + sigma_theta * g_theta2
    gamma1 <- sigma_gamma * g_gamma1
    gamma2 <- sigma_gamma * g_gamma2
  }

  list(
    theta1 = theta1, theta2 = theta2,
    gamma1 = gamma1, gamma2 = gamma2
  )
}

#' Generate correlated twin parameters for SR/MGG models
#'
#' @param n_pairs Number of twin pairs
#' @param sigma_intrinsic SD of the intrinsic genetic parameter
#' @param mu_intrinsic Mean of intrinsic parameter
#' @param sigma_gamma SD of log-extrinsic frailty
#' @param rho Pleiotropy correlation
#' @param r_g Genetic correlation (1 for MZ, 0.5 for DZ)
#' @param intrinsic_sign Direction: -1 if higher parameter = longer life
#' @return List with intrinsic1, intrinsic2, gamma1, gamma2
gen_twin_params_model <- function(n_pairs, sigma_intrinsic, mu_intrinsic,
                                  sigma_gamma = 0, rho = 0, r_g = 1,
                                  intrinsic_sign = 1) {
  rho_eff <- rho * intrinsic_sign
  Sigma <- matrix(c(1, rho_eff, rho_eff, 1), 2, 2)

  g_shared <- MASS::mvrnorm(n_pairs, mu = c(0, 0), Sigma = Sigma)

  if (r_g >= 1) {
    z_int1 <- g_shared[, 1]
    z_int2 <- g_shared[, 1]
    z_gam1 <- g_shared[, 2]
    z_gam2 <- g_shared[, 2]
  } else {
    g_indep <- MASS::mvrnorm(n_pairs, mu = c(0, 0), Sigma = Sigma)
    z_int1 <- g_shared[, 1]
    z_int2 <- r_g * g_shared[, 1] + sqrt(1 - r_g^2) * g_indep[, 1]
    z_gam1 <- g_shared[, 2]
    z_gam2 <- r_g * g_shared[, 2] + sqrt(1 - r_g^2) * g_indep[, 2]
  }

  list(
    intrinsic1 = mu_intrinsic + sigma_intrinsic * z_int1,
    intrinsic2 = mu_intrinsic + sigma_intrinsic * z_int2,
    gamma1 = sigma_gamma * z_gam1,
    gamma2 = sigma_gamma * z_gam2
  )
}

#' Simulate twin pair lifespans and compute Falconer h² (GM model)
#'
#' CRN pattern: if rng_seed is provided, set it once at call site.
#' No internal set.seed() calls.
#'
#' @param sigma_theta SD of log-Gompertz intercept
#' @param m_ex Extrinsic hazard rate
#' @param sigma_gamma SD of log-extrinsic frailty
#' @param rho Pleiotropy correlation
#' @param gamma_heritable Whether gamma is shared between twins
#' @param individual_ext Whether to use individual extrinsic rates
#' @param cutoff Minimum age for inclusion
#' @param n Number of twin pairs
#' @param rng_seed RNG seed (NULL = don't set)
#' @return List with r_mz, r_dz, h2, n_mz, n_dz, diagnostics
sim_twin_h2 <- function(sigma_theta, m_ex, sigma_gamma = 0, rho = 0,
                        gamma_heritable = TRUE, individual_ext = TRUE,
                        cutoff = PARAMS$CUTOFF_AGE,
                        n = PARAMS$N_PAIRS,
                        rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)

  b <- PARAMS$B_GOMP
  t_max <- PARAMS$T_MAX

  # --- MZ twins ---
  params_mz <- gen_twin_params_gm(n,
    r_g = 1.0, sigma_theta, sigma_gamma, rho, gamma_heritable
  )
  if (individual_ext && sigma_gamma > 0) {
    shift <- sigma_gamma^2 / 2
    c_mz1 <- m_ex * exp(params_mz$gamma1 - shift)
    c_mz2 <- m_ex * exp(params_mz$gamma2 - shift)
  } else {
    c_mz1 <- c_mz2 <- m_ex
  }
  L_mz1 <- sim_lifespan_gm(params_mz$theta1, b, c_mz1, t_max)
  L_mz2 <- sim_lifespan_gm(params_mz$theta2, b, c_mz2, t_max)

  n_censored_mz <- attr(L_mz1, "n_censored") + attr(L_mz2, "n_censored")
  max_err_mz <- max(attr(L_mz1, "max_error"), attr(L_mz2, "max_error"))

  # --- DZ twins ---
  params_dz <- gen_twin_params_gm(n,
    r_g = 0.5, sigma_theta, sigma_gamma, rho, gamma_heritable
  )
  if (individual_ext && sigma_gamma > 0) {
    shift <- sigma_gamma^2 / 2
    c_dz1 <- m_ex * exp(params_dz$gamma1 - shift)
    c_dz2 <- m_ex * exp(params_dz$gamma2 - shift)
  } else {
    c_dz1 <- c_dz2 <- m_ex
  }
  L_dz1 <- sim_lifespan_gm(params_dz$theta1, b, c_dz1, t_max)
  L_dz2 <- sim_lifespan_gm(params_dz$theta2, b, c_dz2, t_max)

  n_censored_dz <- attr(L_dz1, "n_censored") + attr(L_dz2, "n_censored")
  max_err_dz <- max(attr(L_dz1, "max_error"), attr(L_dz2, "max_error"))

  # --- Apply cutoff age ---
  keep_mz <- (L_mz1 > cutoff) & (L_mz2 > cutoff)
  keep_dz <- (L_dz1 > cutoff) & (L_dz2 > cutoff)

  if (sum(keep_mz) < 100 || sum(keep_dz) < 100) {
    warning("Very few pairs survived cutoff. Results may be unstable.")
  }

  r_mz <- cor(L_mz1[keep_mz], L_mz2[keep_mz])
  r_dz <- cor(L_dz1[keep_dz], L_dz2[keep_dz])
  h2 <- 2 * (r_mz - r_dz)

  list(
    r_mz = r_mz, r_dz = r_dz, h2 = h2,
    n_mz = sum(keep_mz), n_dz = sum(keep_dz),
    frac_mz = mean(keep_mz), frac_dz = mean(keep_dz),
    mean_age_mz = mean(c(L_mz1[keep_mz], L_mz2[keep_mz])),
    n_censored_mz = n_censored_mz,
    n_censored_dz = n_censored_dz,
    max_err_mz = max_err_mz,
    max_err_dz = max_err_dz
  )
}

#' Simulate twin h² for MGG model
#'
#' @param sigma_q SD of scale parameter q
#' @param mu_q Mean of q
#' @param m_ex Extrinsic hazard
#' @param sigma_gamma SD of log-extrinsic frailty
#' @param rho Pleiotropy correlation
#' @param individual_ext Whether to use individual extrinsic rates
#' @param n_pairs Number of twin pairs
#' @param cutoff_age Minimum age
#' @param rng_seed RNG seed
#' @return List with r_mz, r_dz, h2, diagnostics
sim_twin_h2_mgg <- function(sigma_q, mu_q = PARAMS$MGG_MU_Q,
                            m_ex = PARAMS$M_EX_HIST,
                            sigma_gamma = 0, rho = 0,
                            individual_ext = FALSE,
                            n_pairs = PARAMS$N_PAIRS,
                            cutoff_age = PARAMS$CUTOFF_AGE,
                            rng_seed = NULL,
                            sm_mapping = "compensatory") {
  if (!is.null(rng_seed)) set.seed(rng_seed)

  # With a=a0^q, higher q = longer life, so intrinsic_sign = -1
  params_mz <- gen_twin_params_model(n_pairs, sigma_q, mu_q,
    sigma_gamma, rho, r_g = 1, intrinsic_sign = -1
  )
  if (individual_ext && sigma_gamma > 0) {
    shift <- sigma_gamma^2 / 2
    m_ex_mz1 <- m_ex * exp(params_mz$gamma1 - shift)
    m_ex_mz2 <- m_ex * exp(params_mz$gamma2 - shift)
  } else {
    m_ex_mz1 <- m_ex_mz2 <- m_ex
  }
  t_mz1 <- sim_lifespan_mgg(params_mz$intrinsic1, m_ex_mz1,
    a0 = PARAMS$MGG_A0, b0 = PARAMS$MGG_B0, c_param = PARAMS$MGG_C,
    sm_mapping = sm_mapping
  )
  t_mz2 <- sim_lifespan_mgg(params_mz$intrinsic2, m_ex_mz2,
    a0 = PARAMS$MGG_A0, b0 = PARAMS$MGG_B0, c_param = PARAMS$MGG_C,
    sm_mapping = sm_mapping
  )

  params_dz <- gen_twin_params_model(n_pairs, sigma_q, mu_q,
    sigma_gamma, rho, r_g = 0.5, intrinsic_sign = -1
  )
  if (individual_ext && sigma_gamma > 0) {
    shift <- sigma_gamma^2 / 2
    m_ex_dz1 <- m_ex * exp(params_dz$gamma1 - shift)
    m_ex_dz2 <- m_ex * exp(params_dz$gamma2 - shift)
  } else {
    m_ex_dz1 <- m_ex_dz2 <- m_ex
  }
  t_dz1 <- sim_lifespan_mgg(params_dz$intrinsic1, m_ex_dz1,
    a0 = PARAMS$MGG_A0, b0 = PARAMS$MGG_B0, c_param = PARAMS$MGG_C,
    sm_mapping = sm_mapping
  )
  t_dz2 <- sim_lifespan_mgg(params_dz$intrinsic2, m_ex_dz2,
    a0 = PARAMS$MGG_A0, b0 = PARAMS$MGG_B0, c_param = PARAMS$MGG_C,
    sm_mapping = sm_mapping
  )

  keep_mz <- t_mz1 > cutoff_age & t_mz2 > cutoff_age
  keep_dz <- t_dz1 > cutoff_age & t_dz2 > cutoff_age

  r_mz <- cor(t_mz1[keep_mz], t_mz2[keep_mz])
  r_dz <- cor(t_dz1[keep_dz], t_dz2[keep_dz])

  list(
    r_mz = r_mz, r_dz = r_dz,
    h2 = 2 * (r_mz - r_dz),
    n_mz = sum(keep_mz), n_dz = sum(keep_dz)
  )
}

#' Simulate twin h² for SR model
#'
#' @param sigma_Xc SD of death threshold
#' @param mu_Xc Mean death threshold
#' @param m_ex Extrinsic hazard
#' @param sigma_gamma SD of log-extrinsic frailty
#' @param rho Pleiotropy correlation
#' @param individual_ext Whether to use individual extrinsic rates
#' @param n_pairs Number of twin pairs
#' @param cutoff_age Minimum age
#' @param rng_seed RNG seed
#' @return List with r_mz, r_dz, h2, diagnostics
sim_twin_h2_sr <- function(sigma_Xc, mu_Xc = PARAMS$SR_MU_XC,
                           m_ex = PARAMS$M_EX_HIST,
                           sigma_gamma = 0, rho = 0,
                           individual_ext = FALSE,
                           n_pairs = PARAMS$N_PAIRS,
                           cutoff_age = PARAMS$CUTOFF_AGE,
                           rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)

  params_mz <- gen_twin_params_model(n_pairs, sigma_Xc, mu_Xc,
    sigma_gamma, rho, r_g = 1, intrinsic_sign = -1
  )
  if (individual_ext && sigma_gamma > 0) {
    shift <- sigma_gamma^2 / 2
    m_ex_mz1 <- m_ex * exp(params_mz$gamma1 - shift)
    m_ex_mz2 <- m_ex * exp(params_mz$gamma2 - shift)
  } else {
    m_ex_mz1 <- m_ex_mz2 <- m_ex
  }
  t_mz1 <- sim_lifespan_sr(params_mz$intrinsic1, m_ex_mz1,
    eta = PARAMS$SR_ETA, beta = PARAMS$SR_BETA,
    epsilon = PARAMS$SR_EPSILON, kappa = PARAMS$SR_KAPPA,
    dt = PARAMS$SR_DT
  )
  t_mz2 <- sim_lifespan_sr(params_mz$intrinsic2, m_ex_mz2,
    eta = PARAMS$SR_ETA, beta = PARAMS$SR_BETA,
    epsilon = PARAMS$SR_EPSILON, kappa = PARAMS$SR_KAPPA,
    dt = PARAMS$SR_DT
  )

  params_dz <- gen_twin_params_model(n_pairs, sigma_Xc, mu_Xc,
    sigma_gamma, rho, r_g = 0.5, intrinsic_sign = -1
  )
  if (individual_ext && sigma_gamma > 0) {
    shift <- sigma_gamma^2 / 2
    m_ex_dz1 <- m_ex * exp(params_dz$gamma1 - shift)
    m_ex_dz2 <- m_ex * exp(params_dz$gamma2 - shift)
  } else {
    m_ex_dz1 <- m_ex_dz2 <- m_ex
  }
  t_dz1 <- sim_lifespan_sr(params_dz$intrinsic1, m_ex_dz1,
    eta = PARAMS$SR_ETA, beta = PARAMS$SR_BETA,
    epsilon = PARAMS$SR_EPSILON, kappa = PARAMS$SR_KAPPA,
    dt = PARAMS$SR_DT
  )
  t_dz2 <- sim_lifespan_sr(params_dz$intrinsic2, m_ex_dz2,
    eta = PARAMS$SR_ETA, beta = PARAMS$SR_BETA,
    epsilon = PARAMS$SR_EPSILON, kappa = PARAMS$SR_KAPPA,
    dt = PARAMS$SR_DT
  )

  keep_mz <- t_mz1 > cutoff_age & t_mz2 > cutoff_age
  keep_dz <- t_dz1 > cutoff_age & t_dz2 > cutoff_age

  r_mz <- cor(t_mz1[keep_mz], t_mz2[keep_mz])
  r_dz <- cor(t_dz1[keep_dz], t_dz2[keep_dz])

  list(
    r_mz = r_mz, r_dz = r_dz,
    h2 = 2 * (r_mz - r_dz),
    n_mz = sum(keep_mz), n_dz = sum(keep_dz)
  )
}
