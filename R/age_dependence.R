# ===========================================================================
# age_dependence.R — Age-specific twin dependence measures
# ===========================================================================
# Conditional correlation, concordance, Kendall's tau, and cross-ratio
# as functions of age threshold. Shows structural misfit in dependence.
# ===========================================================================

#' Compute conditional Pearson/Spearman correlation given both twins survive to age t
#'
#' @param L1 Twin 1 lifespans
#' @param L2 Twin 2 lifespans
#' @param keep Logical mask (both survived cutoff)
#' @param t_grid Age grid
#' @param min_n Minimum pairs for stable estimate
#' @return Data frame with age, cor_pearson, cor_spearman, n_pairs
compute_conditional_correlation <- function(L1, L2, keep, t_grid,
                                            min_n = 100L) {
  l1 <- L1[keep]
  l2 <- L2[keep]

  results <- lapply(t_grid, function(t) {
    idx <- l1 > t & l2 > t
    n_pairs <- sum(idx)
    if (n_pairs < min_n) {
      return(data.frame(age = t, cor_pearson = NA_real_,
                        cor_spearman = NA_real_, n_pairs = n_pairs,
                        stringsAsFactors = FALSE))
    }
    data.frame(
      age = t,
      cor_pearson = cor(l1[idx], l2[idx], method = "pearson"),
      cor_spearman = cor(l1[idx], l2[idx], method = "spearman"),
      n_pairs = n_pairs,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#' Compute Kendall's tau by age threshold (subsampled for speed)
#'
#' @param L1 Twin 1 lifespans
#' @param L2 Twin 2 lifespans
#' @param keep Logical mask
#' @param t_grid Age grid
#' @param min_n Minimum pairs
#' @param max_n Maximum pairs for Kendall computation (O(n²))
#' @return Data frame with age, tau, n_pairs
compute_kendall_by_age <- function(L1, L2, keep, t_grid,
                                   min_n = 100L, max_n = 5000L) {
  l1 <- L1[keep]
  l2 <- L2[keep]

  results <- lapply(t_grid, function(t) {
    idx <- which(l1 > t & l2 > t)
    n_pairs <- length(idx)
    if (n_pairs < min_n) {
      return(data.frame(age = t, tau = NA_real_, n_pairs = n_pairs,
                        stringsAsFactors = FALSE))
    }
    if (n_pairs > max_n) {
      # Deterministic subsample: use age-derived seed for reproducibility
      set.seed(as.integer(t * 1000))
      idx <- sample(idx, max_n)
      n_pairs <- max_n
    }
    data.frame(
      age = t,
      tau = cor(l1[idx], l2[idx], method = "kendall"),
      n_pairs = n_pairs,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#' Compute cross-ratio on the diagonal θ(t,t) from paired lifespans
#'
#' The cross-ratio θ(t₁,t₂) = S(t₁,t₂)·f₁₂(t₁,t₂) / [f₁(t₁)·f₂(t₂)]
#' measures local dependence. On the diagonal, estimated via numerical
#' differentiation of the empirical joint survival function.
#'
#' @param L1 Twin 1 lifespans
#' @param L2 Twin 2 lifespans
#' @param keep Logical mask
#' @param t_grid Age grid (should be evenly spaced)
#' @return Data frame with age, cross_ratio, S_joint, dS_dt1, dS_dt2, d2S
compute_cross_ratio_diagonal <- function(L1, L2, keep, t_grid) {
  l1 <- L1[keep]
  l2 <- L2[keep]
  n <- length(l1)
  ng <- length(t_grid)
  dt <- t_grid[2] - t_grid[1]

  # Helper: empirical joint survival at arbitrary (t1, t2)
  S_hat <- function(t1, t2) mean(l1 > t1 & l2 > t2)

  # Empirical joint survival on diagonal and marginals
  S_joint <- sapply(t_grid, function(t) S_hat(t, t))
  S1 <- sapply(t_grid, function(t) mean(l1 > t))
  S2 <- sapply(t_grid, function(t) mean(l2 > t))

  # Oakes (1989) cross-ratio uses PARTIAL DERIVATIVES of the JOINT survival:
  #   θ(t₁,t₂) = S(t₁,t₂) · [∂²S/∂t₁∂t₂] / [(∂S/∂t₁)(∂S/∂t₂)]
  # NOTE: ∂S/∂t₁ at (t,t) is NOT the marginal density f₁(t).
  # It is the conditional sub-density ∫_{t}^∞ f(t₁,v)dv evaluated at t₁=t.

  dS_dt1 <- rep(NA_real_, ng)  # ∂S(t₁,t₂)/∂t₁ at (t,t), holding t₂=t fixed
  dS_dt2 <- rep(NA_real_, ng)  # ∂S(t₁,t₂)/∂t₂ at (t,t), holding t₁=t fixed
  d2S    <- rep(NA_real_, ng)  # ∂²S/∂t₁∂t₂ = bivariate density f(t,t)

  for (i in 2:(ng - 1)) {
    ti_m <- t_grid[i - 1]; ti <- t_grid[i]; ti_p <- t_grid[i + 1]

    # Partial derivatives via central differences on the joint surface
    dS_dt1[i] <- (S_hat(ti_p, ti) - S_hat(ti_m, ti)) / (2 * dt)
    dS_dt2[i] <- (S_hat(ti, ti_p) - S_hat(ti, ti_m)) / (2 * dt)

    # Mixed partial via 2D central difference
    d2S[i] <- (S_hat(ti_p, ti_p) - S_hat(ti_p, ti_m) -
               S_hat(ti_m, ti_p) + S_hat(ti_m, ti_m)) / (4 * dt^2)
  }

  # Oakes cross-ratio: θ(t,t) = S(t,t) · d²S(t,t) / [dS/dt₁(t,t) · dS/dt₂(t,t)]
  # Both dS/dt₁ and dS/dt₂ are negative; their product is positive.
  denom <- dS_dt1 * dS_dt2
  theta <- ifelse(abs(denom) > 1e-16 & !is.na(d2S),
                  S_joint * d2S / denom, NA_real_)

  data.frame(age = t_grid, cross_ratio = theta,
             S_joint = S_joint, dS_dt1 = dS_dt1, dS_dt2 = dS_dt2, d2S = d2S,
             stringsAsFactors = FALSE)
}

#' Run full age-specific dependence analysis
#'
#' @param sigma_theta_true True σ_θ
#' @param sigma_theta_fit Misspecified σ_θ
#' @param sigma_theta_fix Recovery σ_θ
#' @param params Parameter list
#' @return List with age-specific curves for all arms and zygosities
run_age_dependence_analysis <- function(sigma_theta_true,
                                         sigma_theta_fit,
                                         sigma_theta_fix,
                                         params) {
  seed_base <- params$MASTER_SEED + 600000L
  sg <- params$SIGMA_GAMMA_DEF
  rho <- params$RHO_DEF
  m_ex <- params$M_EX_HIST

  # Generate paired lifespans under three DGPs
  dgps <- list(
    true_dgp = sim_twin_lifespans(sigma_theta_true, m_ex,
      sigma_gamma = sg, rho = rho,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = seed_base),
    misspecified = sim_twin_lifespans(sigma_theta_fit, m_ex,
      sigma_gamma = 0, rho = 0,
      gamma_heritable = FALSE, individual_ext = FALSE,
      rng_seed = seed_base + 1000L)
  )
  if (is.finite(sigma_theta_fix)) {
    dgps$recovery <- sim_twin_lifespans(sigma_theta_fix, m_ex,
      sigma_gamma = sg, rho = rho,
      gamma_heritable = TRUE, individual_ext = TRUE,
      rng_seed = seed_base + 2000L)
  }

  cor_grid <- seq(20, 90, by = 2)
  cr_grid <- seq(25, 85, by = 2)

  cor_list <- list()
  tau_list <- list()
  cr_list <- list()

  for (dgp_name in names(dgps)) {
    d <- dgps[[dgp_name]]
    for (zyg in c("mz", "dz")) {
      L1 <- d[[paste0("L_", zyg, "1")]]
      L2 <- d[[paste0("L_", zyg, "2")]]
      keep <- d[[paste0("keep_", zyg)]]

      cc <- compute_conditional_correlation(L1, L2, keep, cor_grid)
      cc$model <- dgp_name
      cc$zygosity <- toupper(zyg)
      cor_list[[paste0(dgp_name, "_", zyg)]] <- cc

      kt <- compute_kendall_by_age(L1, L2, keep, cor_grid)
      kt$model <- dgp_name
      kt$zygosity <- toupper(zyg)
      tau_list[[paste0(dgp_name, "_", zyg)]] <- kt

      # Cross-ratio: MZ only (most informative)
      if (zyg == "mz") {
        cr <- compute_cross_ratio_diagonal(L1, L2, keep, cr_grid)
        cr$model <- dgp_name
        cr_list[[dgp_name]] <- cr
      }
    }
  }

  cor_combined <- do.call(rbind, cor_list)
  rownames(cor_combined) <- NULL
  tau_combined <- do.call(rbind, tau_list)
  rownames(tau_combined) <- NULL
  cr_combined <- do.call(rbind, cr_list)
  rownames(cr_combined) <- NULL

  # Summary: max absolute correlation difference between misspec and true
  summary_metrics <- list()
  for (zyg in c("mz", "dz")) {
    true_cc <- cor_list[[paste0("true_dgp_", zyg)]]
    mis_cc <- cor_list[[paste0("misspecified_", zyg)]]
    ok <- !is.na(true_cc$cor_pearson) & !is.na(mis_cc$cor_pearson)
    if (any(ok)) {
      diff_cor <- mis_cc$cor_pearson[ok] - true_cc$cor_pearson[ok]
      summary_metrics[[paste0("max_cor_diff_", zyg)]] <- max(abs(diff_cor))
      summary_metrics[[paste0("peak_diff_age_", zyg)]] <-
        cor_grid[ok][which.max(abs(diff_cor))]
      summary_metrics[[paste0("mean_abs_cor_diff_", zyg)]] <-
        mean(abs(diff_cor))
    }
  }

  list(
    correlations = cor_combined,
    kendall_tau = tau_combined,
    cross_ratio = cr_combined,
    summary = summary_metrics
  )
}
