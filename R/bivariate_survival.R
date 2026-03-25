# ===========================================================================
# bivariate_survival.R — Bivariate survival surface analysis
# ===========================================================================
# Computes S(t1,t2) = P(T1>t1, T2>t2) for twin pairs and derives
# formal misfit metrics between true DGP and fitted models.
# ===========================================================================

#' Compute empirical bivariate survival surface S(t1,t2)
#'
#' @param L1 Lifespan vector for twin 1
#' @param L2 Lifespan vector for twin 2
#' @param keep Logical mask for inclusion (both survive cutoff)
#' @param t_grid Age grid for evaluation
#' @return List with S (matrix), S1/S2 (marginals), grid, n
compute_joint_survival_surface <- function(L1, L2, keep, t_grid) {
  l1 <- L1[keep]
  l2 <- L2[keep]
  n <- length(l1)

  # Vectorized: indicator matrices (n_pairs × n_grid)
  ind1 <- outer(l1, t_grid, FUN = ">")
  ind2 <- outer(l2, t_grid, FUN = ">")

  # Joint survival: S(t1,t2) = (ind1' %*% ind2) / n
  S <- crossprod(ind1, ind2) / n

  list(S = unclass(S),
       S1 = colMeans(ind1),
       S2 = colMeans(ind2),
       grid = t_grid,
       n = n)
}

#' Compute Integrated Squared Error between two survival surfaces
#'
#' ISE = ∫∫ (S_true - S_fit)² dt1 dt2 via Riemann sum on uniform grid
#'
#' @param S_true Reference surface (from compute_joint_survival_surface)
#' @param S_fit Comparison surface
#' @return Scalar ISE value
compute_surface_ise <- function(S_true, S_fit) {
  delta <- S_fit$S - S_true$S
  dt <- mean(diff(S_true$grid))
  sum(delta^2) * dt^2
}

#' Compute supremum distance between two survival surfaces
#'
#' @param S_true Reference surface
#' @param S_fit Comparison surface
#' @return Scalar sup|S_true - S_fit|
compute_surface_sup <- function(S_true, S_fit) {
  max(abs(S_fit$S - S_true$S))
}

#' Compute Cramér-von Mises type statistic for survival surfaces
#'
#' CvM = n × Σ (S_true - S_fit)² × dF_true
#' where dF_true are bivariate CDF increments used as weights.
#'
#' @param S_true Reference surface
#' @param S_fit Comparison surface
#' @return Scalar CvM statistic
compute_surface_cvm <- function(S_true, S_fit) {
  delta <- S_fit$S - S_true$S
  ng <- nrow(S_true$S)

  # Bivariate CDF increments as weights
  weights <- matrix(0, ng, ng)
  for (i in 2:ng) {
    for (j in 2:ng) {
      # Rectangle probability mass from bivariate CDF increments;
      # pmax clips sampling noise that can make small increments negative
      weights[i, j] <- pmax(
        S_true$S[i - 1, j - 1] - S_true$S[i, j - 1] -
        S_true$S[i - 1, j] + S_true$S[i, j], 0
      )
    }
  }
  S_true$n * sum(delta^2 * weights)
}

#' Compute Integrated Absolute Error (for comparability with existing B5)
#'
#' @param S_true Reference surface
#' @param S_fit Comparison surface
#' @return Scalar IAE value
compute_surface_iae <- function(S_true, S_fit) {
  delta <- S_fit$S - S_true$S
  dt <- mean(diff(S_true$grid))
  sum(abs(delta)) * dt^2
}

#' Compute age-specific concordance curve C(t) = P(T2>t | T1>t)
#'
#' Uses pair bootstrap for SE estimation: resamples twin pairs (not
#' individuals) to correctly account for within-pair correlation.
#' Falls back to DEFF-adjusted binomial SE when bootstrap is disabled.
#'
#' @param L1 Lifespan vector for twin 1
#' @param L2 Lifespan vector for twin 2
#' @param keep Logical mask
#' @param t_grid Age grid
#' @param min_n Minimum number of conditioning pairs
#' @param n_boot Number of pair bootstrap replicates (0 = use DEFF approx)
#' @return Data frame with age, concordance, se, lo95, hi95, n_cond
compute_concordance_curve <- function(L1, L2, keep, t_grid, min_n = 50L,
                                      n_boot = 200L) {
  l1 <- L1[keep]
  l2 <- L2[keep]
  n_pairs <- length(l1)

  results <- lapply(t_grid, function(t) {
    cond <- l1 > t
    n_cond <- sum(cond)
    if (n_cond < min_n) {
      return(data.frame(age = t, concordance = NA_real_,
                        se = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
                        n_cond = n_cond, stringsAsFactors = FALSE))
    }
    p <- mean(l2[cond] > t)

    if (n_boot > 0 && n_cond >= min_n) {
      # Pair bootstrap: resample twin PAIRS, not individuals
      boot_p <- vapply(seq_len(n_boot), function(b) {
        idx <- sample.int(n_pairs, n_pairs, replace = TRUE)
        b_cond <- l1[idx] > t
        if (sum(b_cond) < 10L) return(NA_real_)
        mean(l2[idx][b_cond] > t)
      }, numeric(1))
      boot_p <- boot_p[!is.na(boot_p)]
      se <- if (length(boot_p) > 10) sd(boot_p) else sqrt(1.3 * p * (1 - p) / n_cond)
      # Percentile bootstrap CI
      lo95 <- if (length(boot_p) > 10) quantile(boot_p, 0.025) else p - 1.96 * se
      hi95 <- if (length(boot_p) > 10) quantile(boot_p, 0.975) else p + 1.96 * se
    } else {
      # DEFF-adjusted binomial SE as fallback
      deff <- 1.3
      se <- sqrt(deff * p * (1 - p) / n_cond)
      lo95 <- p - 1.96 * se
      hi95 <- p + 1.96 * se
    }

    data.frame(age = t, concordance = p, se = se,
               lo95 = unname(lo95), hi95 = unname(hi95),
               n_cond = n_cond, stringsAsFactors = FALSE)
  })
  do.call(rbind, results)
}

#' Run full bivariate survival surface analysis
#'
#' Generates paired lifespans under oracle, misspecified, and recovery DGPs,
#' computes joint survival surfaces, and returns all misfit metrics.
#'
#' @param sigma_theta_true True σ_θ
#' @param sigma_theta_fit Misspecified σ_θ
#' @param sigma_theta_fix Recovery σ_θ
#' @param params Parameter list
#' @return List with surfaces, metrics, concordance curves
run_bivariate_surface_analysis <- function(sigma_theta_true,
                                            sigma_theta_fit,
                                            sigma_theta_fix,
                                            params) {
  seed_base <- params$MASTER_SEED + 400000L
  sg <- params$SIGMA_GAMMA_DEF
  rho <- params$RHO_DEF
  m_ex <- params$M_EX_HIST

  # Generate paired lifespans under three DGPs
  true_data <- sim_twin_lifespans(sigma_theta_true, m_ex,
    sigma_gamma = sg, rho = rho,
    gamma_heritable = TRUE, individual_ext = TRUE,
    rng_seed = seed_base)

  misspec_data <- sim_twin_lifespans(sigma_theta_fit, m_ex,
    sigma_gamma = 0, rho = 0,
    gamma_heritable = FALSE, individual_ext = FALSE,
    rng_seed = seed_base + 1000L)

  has_fix <- is.finite(sigma_theta_fix)
  fix_data <- if (has_fix) {
    sim_twin_lifespans(sigma_theta_fix, m_ex,
      sigma_gamma = sg, rho = rho,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = seed_base + 2000L)
  }

  # Bivariate surface grid (2-year steps, 20-95)
  t_grid <- seq(20, 95, by = 2)
  conc_grid <- seq(20, 90, by = 1)

  results <- list(t_grid = t_grid,
                  n_mz = true_data$n_mz,
                  n_dz = true_data$n_dz)

  for (zyg in c("mz", "dz")) {
    L1_key <- paste0("L_", zyg, "1")
    L2_key <- paste0("L_", zyg, "2")
    keep_key <- paste0("keep_", zyg)

    S_true <- compute_joint_survival_surface(
      true_data[[L1_key]], true_data[[L2_key]],
      true_data[[keep_key]], t_grid)

    S_misspec <- compute_joint_survival_surface(
      misspec_data[[L1_key]], misspec_data[[L2_key]],
      misspec_data[[keep_key]], t_grid)

    # Store surfaces for plotting
    results[[paste0("surface_true_", zyg)]] <- S_true
    results[[paste0("surface_misspec_", zyg)]] <- S_misspec

    # Misfit metrics: misspecified vs true
    results[[paste0("ise_misspec_", zyg)]] <-
      compute_surface_ise(S_true, S_misspec)
    results[[paste0("sup_misspec_", zyg)]] <-
      compute_surface_sup(S_true, S_misspec)
    results[[paste0("cvm_misspec_", zyg)]] <-
      compute_surface_cvm(S_true, S_misspec)
    results[[paste0("iae_misspec_", zyg)]] <-
      compute_surface_iae(S_true, S_misspec)

    if (has_fix) {
      S_fix <- compute_joint_survival_surface(
        fix_data[[L1_key]], fix_data[[L2_key]],
        fix_data[[keep_key]], t_grid)

      results[[paste0("surface_fix_", zyg)]] <- S_fix
      results[[paste0("ise_fix_", zyg)]] <-
        compute_surface_ise(S_true, S_fix)
      results[[paste0("sup_fix_", zyg)]] <-
        compute_surface_sup(S_true, S_fix)
      results[[paste0("cvm_fix_", zyg)]] <-
        compute_surface_cvm(S_true, S_fix)
      results[[paste0("iae_fix_", zyg)]] <-
        compute_surface_iae(S_true, S_fix)
    }

    # Concordance curves
    results[[paste0("conc_true_", zyg)]] <- compute_concordance_curve(
      true_data[[L1_key]], true_data[[L2_key]],
      true_data[[keep_key]], conc_grid)
    results[[paste0("conc_misspec_", zyg)]] <- compute_concordance_curve(
      misspec_data[[L1_key]], misspec_data[[L2_key]],
      misspec_data[[keep_key]], conc_grid)
    if (has_fix) {
      results[[paste0("conc_fix_", zyg)]] <- compute_concordance_curve(
        fix_data[[L1_key]], fix_data[[L2_key]],
        fix_data[[keep_key]], conc_grid)
    }
  }

  results
}
