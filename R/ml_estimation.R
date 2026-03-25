# ===========================================================================
# ml_estimation.R — Full maximum likelihood estimation for misspecified GM model
# ===========================================================================
# Demonstrates that ML converges to the same pseudo-true parameter as
# moment-based calibration under omitted extrinsic frailty misspecification.
# Uses Gauss-Hermite quadrature to integrate over shared frailty theta.
# ===========================================================================

# ---------------------------------------------------------------------------
# Per-pair log-likelihood (for debugging / unit verification)
# ---------------------------------------------------------------------------

#' Log-likelihood for one MZ twin pair under shared-frailty GM model
#'
#' Integrates over shared theta via Gauss-Hermite quadrature.
#' Both twins assumed uncensored (exact death times observed).
#'
#' @param t1 Lifespan of twin 1
#' @param t2 Lifespan of twin 2
#' @param mu_theta Mean of theta distribution
#' @param sigma_theta SD of theta distribution
#' @param m_ex Population extrinsic hazard (constant)
#' @param b Gompertz slope
#' @param gh_nodes Pre-computed GH quadrature nodes (standard normal scale)
#' @param gh_weights Pre-computed GH quadrature weights
#' @return Scalar log-likelihood
loglik_mz_pair <- function(t1, t2, mu_theta, sigma_theta, m_ex, b,
                           gh_nodes, gh_weights) {
  theta <- mu_theta + sqrt(2) * sigma_theta * gh_nodes
  a <- exp(theta)

  bt1 <- b * t1
  bt2 <- b * t2
  ebt1 <- exp(bt1)
  ebt2 <- exp(bt2)

  h1 <- m_ex + a * ebt1
  h2 <- m_ex + a * ebt2

  H1 <- m_ex * t1 + (a / b) * (ebt1 - 1)
  H2 <- m_ex * t2 + (a / b) * (ebt2 - 1)

  log_integrand <- log(h1) + log(h2) - H1 - H2

  log_w <- log(gh_weights)
  log_terms <- log_w + log_integrand
  max_log <- max(log_terms)
  log(sum(exp(log_terms - max_log))) + max_log - 0.5 * log(pi)
}

# ---------------------------------------------------------------------------
# Vectorized total log-likelihood (production)
# ---------------------------------------------------------------------------

#' Total MZ log-likelihood across all pairs
#'
#' Vectorized over pairs with chunking to bound memory.
#' Uses log-sum-exp for numerical stability.
#'
#' @param sigma_theta SD of theta (the parameter being estimated)
#' @param t1_vec Vector of twin-1 lifespans
#' @param t2_vec Vector of twin-2 lifespans
#' @param mu_theta Mean of theta distribution (fixed)
#' @param m_ex Extrinsic hazard (fixed)
#' @param b Gompertz slope (fixed)
#' @param n_quad Number of GH quadrature points
#' @param cutoff Age cutoff for left-truncation (0 = none)
#' @param t_max Maximum age (right-censoring boundary)
#' @return Scalar total log-likelihood
loglik_mz_total <- function(sigma_theta, t1_vec, t2_vec,
                            mu_theta = PARAMS$MU_THETA,
                            m_ex = PARAMS$M_EX_HIST,
                            b = PARAMS$B_GOMP,
                            n_quad = 64L,
                            cutoff = PARAMS$CUTOFF_AGE,
                            t_max = PARAMS$T_MAX,
                            constrain_mean_A = FALSE,
                            mu_theta_oracle = PARAMS$MU_THETA,
                            sigma_theta_oracle = NULL) {
  if (sigma_theta <= 0) return(-1e15)

  # --- Mean-frailty constraint (Gemini review fix) ---
  # When constrain_mean_A = TRUE, adjust mu_theta so that E[A] = E[exp(theta)]
  # stays constant as sigma_theta varies:
  #   E[A] = exp(mu_theta + sigma_theta^2 / 2) = constant
  #   => mu_theta(sigma_theta) = C - sigma_theta^2 / 2
  #   where C = mu_theta_oracle + sigma_theta_oracle^2 / 2
  if (constrain_mean_A && !is.null(sigma_theta_oracle)) {
    C_mean <- mu_theta_oracle + sigma_theta_oracle^2 / 2
    mu_theta <- C_mean - sigma_theta^2 / 2
  }

  gh <- statmod::gauss.quad(n_quad, "hermite")

  theta <- mu_theta + sqrt(2) * sigma_theta * gh$nodes
  a <- exp(theta)
  log_w <- log(gh$weights)
  nq <- length(gh$nodes)

  # --- Left-truncation correction ---
  log_trunc_prob <- 0
  if (cutoff > 0) {
    ebc <- exp(b * cutoff)
    H_c <- m_ex * cutoff + (a / b) * (ebc - 1)
    log_S2_c <- -2 * H_c
    log_trunc_terms <- log_w + log_S2_c
    max_lt <- max(log_trunc_terms)
    log_trunc_prob <- log(sum(exp(log_trunc_terms - max_lt))) + max_lt - 0.5 * log(pi)
  }

  # --- Identify censoring status (exact comparison, GPT review fix) ---
  d1 <- as.numeric(t1_vec < t_max)  # 1 = event, 0 = censored at t_max
  d2 <- as.numeric(t2_vec < t_max)

  n_pairs <- length(t1_vec)
  ll <- 0

  chunk_size <- 5000L
  for (start in seq(1L, n_pairs, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n_pairs)
    idx <- start:end
    nc <- length(idx)

    t1_mat <- matrix(t1_vec[idx], nrow = nc, ncol = nq)
    t2_mat <- matrix(t2_vec[idx], nrow = nc, ncol = nq)
    a_mat  <- matrix(a, nrow = nc, ncol = nq, byrow = TRUE)
    d1_vec <- d1[idx]
    d2_vec <- d2[idx]

    bt1 <- b * t1_mat
    bt2 <- b * t2_mat
    ebt1 <- exp(bt1)
    ebt2 <- exp(bt2)

    # Hazards
    h1 <- m_ex + a_mat * ebt1
    h2 <- m_ex + a_mat * ebt2

    # Cumulative hazards
    H1 <- m_ex * t1_mat + (a_mat / b) * (ebt1 - 1)
    H2 <- m_ex * t2_mat + (a_mat / b) * (ebt2 - 1)

    # Log-integrand with censoring:
    #   uncensored: log h(t) - H(t)
    #   censored:   -H(t)   (survival contribution only)
    # Per-twin: d * log(h) - H
    d1_mat <- matrix(d1_vec, nrow = nc, ncol = nq)
    d2_mat <- matrix(d2_vec, nrow = nc, ncol = nq)
    log_integ <- d1_mat * log(h1) + d2_mat * log(h2) - H1 - H2

    log_w_mat <- matrix(log_w, nrow = nc, ncol = nq, byrow = TRUE)
    log_terms <- log_w_mat + log_integ

    max_log <- apply(log_terms, 1, max)
    ll_chunk <- log(rowSums(exp(log_terms - max_log))) + max_log - 0.5 * log(pi)

    ll_chunk[!is.finite(ll_chunk)] <- -1e10
    ll <- ll + sum(ll_chunk)
  }

  # Left-truncation correction
  ll <- ll - n_pairs * log_trunc_prob

  ll
}

# ---------------------------------------------------------------------------
# ML optimizer wrapper
# ---------------------------------------------------------------------------

#' Fit misspecified one-component GM model via maximum likelihood (MZ only)
#'
#' Estimates sigma_theta by maximizing the bivariate MZ log-likelihood.
#' When constrain_mean_A = TRUE, mu_theta adjusts with sigma_theta to keep
#' E[A] = E[exp(theta)] constant — matching the implicit normalization of
#' the moment estimator.
#'
#' @param t1_vec Vector of twin-1 lifespans (after cutoff filtering)
#' @param t2_vec Vector of twin-2 lifespans (after cutoff filtering)
#' @param mu_theta Fixed mean of theta (used when constrain_mean_A = FALSE)
#' @param m_ex Fixed extrinsic hazard
#' @param b Fixed Gompertz slope
#' @param n_quad Number of GH quadrature points
#' @param constrain_mean_A If TRUE, adjust mu_theta to keep E[A] constant
#' @param sigma_theta_oracle Oracle sigma_theta (needed for mean-A constraint)
#' @param lo Lower bound for sigma_theta
#' @param hi Upper bound for sigma_theta
#' @return List with sigma_theta_ml, loglik, se, convergence
fit_ml_misspecified <- function(t1_vec, t2_vec,
                                mu_theta = PARAMS$MU_THETA,
                                m_ex = PARAMS$M_EX_HIST,
                                b = PARAMS$B_GOMP,
                                n_quad = 64L,
                                cutoff = PARAMS$CUTOFF_AGE,
                                t_max = PARAMS$T_MAX,
                                constrain_mean_A = TRUE,
                                sigma_theta_oracle = NULL,
                                lo = 0.1, hi = 4.0) {
  neg_ll <- function(sigma_theta) {
    -loglik_mz_total(sigma_theta, t1_vec, t2_vec,
                     mu_theta = mu_theta, m_ex = m_ex, b = b,
                     n_quad = n_quad, cutoff = cutoff, t_max = t_max,
                     constrain_mean_A = constrain_mean_A,
                     mu_theta_oracle = mu_theta,
                     sigma_theta_oracle = sigma_theta_oracle)
  }

  opt <- optim(
    par = 1.3,
    fn = neg_ll,
    method = "Brent",
    lower = lo, upper = hi,
    hessian = TRUE
  )

  se <- tryCatch(sqrt(1 / opt$hessian[1, 1]), error = function(e) NA_real_)

  list(
    sigma_theta_ml = opt$par,
    loglik = -opt$value,
    se = se,
    convergence = opt$convergence
  )
}

# ---------------------------------------------------------------------------
# ML-vs-moment comparison (single seed)
# ---------------------------------------------------------------------------

#' Run one ML-vs-moment comparison for a given DGP
#'
#' @param sigma_theta_true True (oracle) sigma_theta
#' @param sigma_gamma Extrinsic heterogeneity SD
#' @param rho Pleiotropy correlation
#' @param m_ex Extrinsic hazard
#' @param n_pairs Number of MZ twin pairs
#' @param seed RNG seed
#' @param n_quad GH quadrature points
#' @return Data frame row with ML and moment estimates
run_ml_vs_moment_single <- function(sigma_theta_true, sigma_gamma, rho,
                                     m_ex = PARAMS$M_EX_HIST,
                                     n_pairs = PARAMS$N_PAIRS,
                                     seed = 1L,
                                     n_quad = 64L) {
  data <- sim_twin_lifespans(
    sigma_theta = sigma_theta_true,
    m_ex = m_ex,
    sigma_gamma = sigma_gamma,
    rho = rho,
    gamma_heritable = (sigma_gamma > 0),
    individual_ext = (sigma_gamma > 0),
    n = n_pairs,
    rng_seed = seed
  )

  t1 <- data$L_mz1[data$keep_mz]
  t2 <- data$L_mz2[data$keep_mz]
  r_mz_obs <- cor(t1, t2)

  ml_fit <- fit_ml_misspecified(t1, t2, n_quad = n_quad,
    constrain_mean_A = TRUE, sigma_theta_oracle = sigma_theta_true)
  sigma_theta_moment <- calibrate_sigma_theta(r_mz_obs, m_ex)

  data.frame(
    sigma_gamma = sigma_gamma,
    rho = rho,
    n_pairs = n_pairs,
    seed = seed,
    r_mz_obs = r_mz_obs,
    sigma_theta_true = sigma_theta_true,
    sigma_theta_ml = ml_fit$sigma_theta_ml,
    sigma_theta_moment = sigma_theta_moment,
    ml_se = ml_fit$se,
    ml_loglik = ml_fit$loglik,
    ml_convergence = ml_fit$convergence,
    infl_ml_pct = 100 * (ml_fit$sigma_theta_ml / sigma_theta_true - 1),
    infl_moment_pct = 100 * (sigma_theta_moment / sigma_theta_true - 1),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# ML-vs-moment comparison (multi-seed, parallelized)
# ---------------------------------------------------------------------------

#' Run ML-vs-moment comparison across multiple seeds
#'
#' @param sigma_theta_true True sigma_theta
#' @param sigma_gamma Extrinsic heterogeneity SD
#' @param rho Pleiotropy correlation
#' @param n_pairs Number of pairs
#' @param n_seeds Number of MC replicates
#' @param seed_base Base seed
#' @return Data frame with one row per seed
run_ml_vs_moment <- function(sigma_theta_true, sigma_gamma, rho,
                              n_pairs = PARAMS$N_PAIRS,
                              n_seeds = 20L,
                              seed_base = PARAMS$MASTER_SEED + 700000L) {
  seeds <- seed_base + seq_len(n_seeds)
  mc_cores <- getOption("mc.cores", PARAMS$INNER_MC_CORES %||% 4L)

  results <- parallel::mclapply(seeds, function(s) {
    tryCatch(
      run_ml_vs_moment_single(
        sigma_theta_true = sigma_theta_true,
        sigma_gamma = sigma_gamma,
        rho = rho,
        n_pairs = n_pairs,
        seed = s
      ),
      error = function(e) {
        warning(sprintf("ML fit failed for seed %d: %s", s, e$message))
        NULL
      }
    )
  }, mc.cores = mc_cores)

  results <- results[!vapply(results, is.null, logical(1))]
  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# Full experimental design
# ---------------------------------------------------------------------------

#' Run the full ML-vs-moment experimental design
#'
#' 4 DGP conditions x 3 sample sizes. Replicates vary by sample size.
#'
#' @param sigma_theta_true True sigma_theta from oracle calibration
#' @return List with results data frame and summary
run_ml_experiment <- function(sigma_theta_true) {
  conditions <- data.frame(
    sigma_gamma = c(0,    0.40, 0.40, 0.65),
    rho         = c(0,    0.40, 0.00, 0.40),
    label       = c("sanity", "default", "no_pleiotropy", "extreme"),
    stringsAsFactors = FALSE
  )

  sample_sizes <- data.frame(
    n_pairs = c(5000L,  20000L, 50000L),
    n_seeds = c(50L,    20L,    10L),
    stringsAsFactors = FALSE
  )

  all_results <- list()

  for (i in seq_len(nrow(conditions))) {
    cond <- conditions[i, ]
    for (j in seq_len(nrow(sample_sizes))) {
      ss <- sample_sizes[j, ]
      cat(sprintf("[ML experiment] %s, n=%dk, %d seeds...\n",
                  cond$label, ss$n_pairs / 1000, ss$n_seeds))

      res <- run_ml_vs_moment(
        sigma_theta_true = sigma_theta_true,
        sigma_gamma = cond$sigma_gamma,
        rho = cond$rho,
        n_pairs = ss$n_pairs,
        n_seeds = ss$n_seeds,
        seed_base = PARAMS$MASTER_SEED + 700000L +
          (i - 1L) * 10000L + (j - 1L) * 1000L
      )
      res$label <- cond$label
      all_results[[length(all_results) + 1L]] <- res
    }
  }

  combined <- do.call(rbind, all_results)

  summary_df <- do.call(rbind, lapply(
    split(combined, interaction(combined$label, combined$n_pairs)),
    function(d) {
      data.frame(
        label         = d$label[1],
        sigma_gamma   = d$sigma_gamma[1],
        rho           = d$rho[1],
        n_pairs       = d$n_pairs[1],
        n_seeds       = nrow(d),
        ml_mean_infl  = mean(d$infl_ml_pct),
        ml_se_infl    = sd(d$infl_ml_pct) / sqrt(nrow(d)),
        mom_mean_infl = mean(d$infl_moment_pct),
        mom_se_infl   = sd(d$infl_moment_pct) / sqrt(nrow(d)),
        ml_minus_mom  = mean(d$sigma_theta_ml - d$sigma_theta_moment),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(summary_df) <- NULL

  list(results = combined, summary = summary_df)
}

# ---------------------------------------------------------------------------
# Export for manuscript
# ---------------------------------------------------------------------------

#' Export ML comparison results for manuscript
#'
#' Writes tables/ml_vs_moment.csv and tables/ml_scalars.json.
#'
#' @param ml_experiment Output of run_ml_experiment()
#' @return Invisible list of key scalars
export_ml_results <- function(ml_experiment) {
  write.csv(ml_experiment$summary,
            "tables/ml_vs_moment.csv", row.names = FALSE)

  summ <- ml_experiment$summary
  default_50k <- summ[summ$label == "default" & summ$n_pairs == 50000, ]

  scalars <- list(
    ml_infl_pct     = round(default_50k$ml_mean_infl, 1),
    ml_infl_se_pct  = round(default_50k$ml_se_infl, 1),
    mom_infl_pct    = round(default_50k$mom_mean_infl, 1),
    mom_infl_se_pct = round(default_50k$mom_se_infl, 1),
    ml_minus_mom    = round(default_50k$ml_minus_mom, 3)
  )
  jsonlite::write_json(scalars, "tables/ml_scalars.json",
                        auto_unbox = TRUE, pretty = TRUE)

  invisible(scalars)
}

# ---------------------------------------------------------------------------
# Figure generation wrapper (for targets pipeline)
# ---------------------------------------------------------------------------

#' Generate ML comparison figures
#'
#' @param ml_experiment Output of run_ml_experiment()
#' @return Invisible list of plot objects
generate_ml_figures <- function(ml_experiment) {
  plots_dir <- "figures"

  p1 <- plot_ml_vs_moment_scatter(ml_experiment)
  ggplot2::ggsave(file.path(plots_dir, "fig_ml_scatter.png"),
                  p1, width = 8, height = 7, dpi = 300)

  p2 <- plot_ml_convergence(ml_experiment, "default")
  ggplot2::ggsave(file.path(plots_dir, "fig_ml_convergence.png"),
                  p2, width = 7, height = 8, dpi = 300)

  invisible(list(scatter = p1, convergence = p2))
}
