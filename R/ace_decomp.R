# ===========================================================================
# ace_decomp.R — ACE/ADE twin decomposition via OpenMx
# ===========================================================================
# Classical structural equation model decomposing twin lifespan variance
# into A (additive genetic), C (shared environment), E (unique environment).
# Used as sensitivity analysis alongside Falconer h².
# ===========================================================================

#' Fit univariate ACE twin model using OpenMx
#'
#' @param L_mz1 MZ twin 1 lifespans (after cutoff filtering)
#' @param L_mz2 MZ twin 2 lifespans
#' @param L_dz1 DZ twin 1 lifespans
#' @param L_dz2 DZ twin 2 lifespans
#' @param transform "none", "log", or "rank" for normality handling
#' @return List with A, C, E proportions, fit statistics
fit_ace_model <- function(L_mz1, L_mz2, L_dz1, L_dz2,
                          transform = c("none", "log", "rank")) {
  transform <- match.arg(transform)

  na_result <- list(A = NA, C = NA, E = NA, total_var = NA,
                    a2 = NA, c2 = NA, e2 = NA, h2_ace = NA,
                    minus2ll = NA, status = NA,
                    converged = FALSE, transform = transform)

  if (!requireNamespace("OpenMx", quietly = TRUE)) {
    warning("OpenMx not available; returning NA results")
    return(na_result)
  }

  # Apply transformation
  if (transform == "rank") {
    # Joint ranking across all four twin vectors (van der Waerden scores)
    all_vals <- c(L_mz1, L_mz2, L_dz1, L_dz2)
    n_all <- length(all_vals)
    all_ranks <- rank(all_vals, ties.method = "average")
    all_scores <- qnorm((all_ranks - 0.5) / n_all)
    n1 <- length(L_mz1); n2 <- length(L_mz2)
    n3 <- length(L_dz1); n4 <- length(L_dz2)
    L_mz1 <- all_scores[1:n1]
    L_mz2 <- all_scores[(n1 + 1):(n1 + n2)]
    L_dz1 <- all_scores[(n1 + n2 + 1):(n1 + n2 + n3)]
    L_dz2 <- all_scores[(n1 + n2 + n3 + 1):n_all]
  } else if (transform == "log") {
    L_mz1 <- log(pmax(L_mz1, 1))
    L_mz2 <- log(pmax(L_mz2, 1))
    L_dz1 <- log(pmax(L_dz1, 1))
    L_dz2 <- log(pmax(L_dz2, 1))
  }

  all_tf <- c(L_mz1, L_mz2, L_dz1, L_dz2)
  total_var <- var(all_tf)
  m_start <- mean(all_tf)

  # Starting values: reasonable split
  sv_a <- sqrt(total_var * 0.4)
  sv_c <- sqrt(total_var * 0.1)
  sv_e <- sqrt(total_var * 0.5)

  mz_data <- data.frame(twin1 = L_mz1, twin2 = L_mz2)
  dz_data <- data.frame(twin1 = L_dz1, twin2 = L_dz2)

  # Define variance component path coefficients
  # Note: matrix names must differ from parameter labels to avoid OpenMx clash
  pathA <- OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = sv_a,
                             labels = "par_a", name = "pathA", lbound = 1e-6)
  pathC <- OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = sv_c,
                             labels = "par_c", name = "pathC", lbound = 0)
  pathE <- OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = sv_e,
                             labels = "par_e", name = "pathE", lbound = 1e-6)
  meanG <- OpenMx::mxMatrix("Full", 1, 2, free = TRUE, values = m_start,
                             labels = "par_mean", name = "meanG")

  # Variance components (squared path coefficients)
  covA <- OpenMx::mxAlgebra(pathA * pathA, name = "A")
  covC <- OpenMx::mxAlgebra(pathC * pathC, name = "C")
  covE <- OpenMx::mxAlgebra(pathE * pathE, name = "E")
  covV <- OpenMx::mxAlgebra(A + C + E, name = "V")

  # Expected covariance matrices
  # MZ: Cov(twin1, twin2) = A + C (share all A and all C)
  # DZ: Cov(twin1, twin2) = 0.5*A + C (share half A, all C)
  covMZ <- OpenMx::mxAlgebra(rbind(cbind(V, A + C),
                                    cbind(A + C, V)), name = "expCovMZ")
  covDZ <- OpenMx::mxAlgebra(rbind(cbind(V, 0.5 * A + C),
                                    cbind(0.5 * A + C, V)), name = "expCovDZ")

  datMZ <- OpenMx::mxData(mz_data, type = "raw")
  datDZ <- OpenMx::mxData(dz_data, type = "raw")

  objMZ <- OpenMx::mxExpectationNormal(covariance = "expCovMZ",
                                        means = "meanG",
                                        dimnames = c("twin1", "twin2"))
  objDZ <- OpenMx::mxExpectationNormal(covariance = "expCovDZ",
                                        means = "meanG",
                                        dimnames = c("twin1", "twin2"))

  fitML <- OpenMx::mxFitFunctionML()

  # Submodels
  modelMZ <- OpenMx::mxModel("MZ",
    pathA, pathC, pathE, meanG,
    covA, covC, covE, covV, covMZ,
    datMZ, objMZ, fitML)

  modelDZ <- OpenMx::mxModel("DZ",
    pathA, pathC, pathE, meanG,
    covA, covC, covE, covV, covDZ,
    datDZ, objDZ, fitML)

  # Multi-group model
  fitACE <- OpenMx::mxFitFunctionMultigroup(c("MZ", "DZ"))
  aceModel <- OpenMx::mxModel("ACE", modelMZ, modelDZ, fitACE)

  aceFit <- tryCatch(
    OpenMx::mxRun(aceModel, silent = TRUE),
    error = function(e) {
      warning("OpenMx ACE failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(aceFit)) return(na_result)

  # Extract results
  a_val <- aceFit$MZ$pathA$values[1, 1]^2
  c_val <- aceFit$MZ$pathC$values[1, 1]^2
  e_val <- aceFit$MZ$pathE$values[1, 1]^2
  total <- a_val + c_val + e_val

  list(
    A = a_val, C = c_val, E = e_val,
    total_var = total,
    a2 = a_val / total,
    c2 = c_val / total,
    e2 = e_val / total,
    h2_ace = a_val / total,
    minus2ll = aceFit$output$fit,
    status = aceFit$output$status$code,
    converged = (aceFit$output$status$code == 0),
    transform = transform
  )
}

#' Fit AE sub-model (drop C, compare with ACE via LRT)
#'
#' @param L_mz1 MZ twin 1 lifespans
#' @param L_mz2 MZ twin 2 lifespans
#' @param L_dz1 DZ twin 1 lifespans
#' @param L_dz2 DZ twin 2 lifespans
#' @param transform Transformation to apply
#' @return List with A, E proportions, fit statistics
fit_ae_model <- function(L_mz1, L_mz2, L_dz1, L_dz2,
                         transform = c("none", "log", "rank")) {
  transform <- match.arg(transform)

  na_result <- list(A = NA, E = NA, total_var = NA,
                    a2 = NA, e2 = NA, h2_ae = NA,
                    minus2ll = NA, status = NA,
                    converged = FALSE, transform = transform)

  if (!requireNamespace("OpenMx", quietly = TRUE)) return(na_result)

  if (transform == "rank") {
    all_vals <- c(L_mz1, L_mz2, L_dz1, L_dz2)
    n_all <- length(all_vals)
    all_ranks <- rank(all_vals, ties.method = "average")
    all_scores <- qnorm((all_ranks - 0.5) / n_all)
    n1 <- length(L_mz1); n2 <- length(L_mz2)
    n3 <- length(L_dz1); n4 <- length(L_dz2)
    L_mz1 <- all_scores[1:n1]
    L_mz2 <- all_scores[(n1 + 1):(n1 + n2)]
    L_dz1 <- all_scores[(n1 + n2 + 1):(n1 + n2 + n3)]
    L_dz2 <- all_scores[(n1 + n2 + n3 + 1):n_all]
  } else if (transform == "log") {
    L_mz1 <- log(pmax(L_mz1, 1))
    L_mz2 <- log(pmax(L_mz2, 1))
    L_dz1 <- log(pmax(L_dz1, 1))
    L_dz2 <- log(pmax(L_dz2, 1))
  }

  all_tf <- c(L_mz1, L_mz2, L_dz1, L_dz2)
  total_var <- var(all_tf)
  m_start <- mean(all_tf)

  sv_a <- sqrt(total_var * 0.5)
  sv_e <- sqrt(total_var * 0.5)

  mz_data <- data.frame(twin1 = L_mz1, twin2 = L_mz2)
  dz_data <- data.frame(twin1 = L_dz1, twin2 = L_dz2)

  pathA <- OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = sv_a,
                             labels = "par_a", name = "pathA", lbound = 1e-6)
  pathE <- OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = sv_e,
                             labels = "par_e", name = "pathE", lbound = 1e-6)
  meanG <- OpenMx::mxMatrix("Full", 1, 2, free = TRUE, values = m_start,
                             labels = "par_mean", name = "meanG")

  covA <- OpenMx::mxAlgebra(pathA * pathA, name = "A")
  covE <- OpenMx::mxAlgebra(pathE * pathE, name = "E")
  covV <- OpenMx::mxAlgebra(A + E, name = "V")

  covMZ <- OpenMx::mxAlgebra(rbind(cbind(V, A),
                                    cbind(A, V)), name = "expCovMZ")
  covDZ <- OpenMx::mxAlgebra(rbind(cbind(V, 0.5 * A),
                                    cbind(0.5 * A, V)), name = "expCovDZ")

  datMZ <- OpenMx::mxData(mz_data, type = "raw")
  datDZ <- OpenMx::mxData(dz_data, type = "raw")

  objMZ <- OpenMx::mxExpectationNormal(covariance = "expCovMZ",
                                        means = "meanG",
                                        dimnames = c("twin1", "twin2"))
  objDZ <- OpenMx::mxExpectationNormal(covariance = "expCovDZ",
                                        means = "meanG",
                                        dimnames = c("twin1", "twin2"))

  fitML <- OpenMx::mxFitFunctionML()

  modelMZ <- OpenMx::mxModel("MZ", pathA, pathE, meanG,
    covA, covE, covV, covMZ, datMZ, objMZ, fitML)
  modelDZ <- OpenMx::mxModel("DZ", pathA, pathE, meanG,
    covA, covE, covV, covDZ, datDZ, objDZ, fitML)

  fitAE <- OpenMx::mxFitFunctionMultigroup(c("MZ", "DZ"))
  aeModel <- OpenMx::mxModel("AE", modelMZ, modelDZ, fitAE)

  aeFit <- tryCatch(
    OpenMx::mxRun(aeModel, silent = TRUE),
    error = function(e) { warning("OpenMx AE failed: ", e$message); NULL }
  )

  if (is.null(aeFit)) return(na_result)

  a_val <- aeFit$MZ$pathA$values[1, 1]^2
  e_val <- aeFit$MZ$pathE$values[1, 1]^2
  total <- a_val + e_val

  list(
    A = a_val, E = e_val,
    total_var = total,
    a2 = a_val / total,
    e2 = e_val / total,
    h2_ae = a_val / total,
    minus2ll = aeFit$output$fit,
    status = aeFit$output$status$code,
    converged = (aeFit$output$status$code == 0),
    transform = transform
  )
}

#' Run ACE analysis across oracle, misspecified, and recovery arms
#'
#' Fits ACE and AE models under multiple DGPs and transformations.
#' Compares how the additive genetic component A inflates under
#' misspecification.
#'
#' @param sigma_theta_true True σ_θ
#' @param sigma_theta_fit Misspecified σ_θ
#' @param sigma_theta_fix Recovery σ_θ
#' @param params Parameter list
#' @return List with per-arm ACE results and summary data frame
run_ace_analysis <- function(sigma_theta_true, sigma_theta_fit,
                              sigma_theta_fix, params) {
  seed_base <- params$MASTER_SEED + 500000L
  sg <- params$SIGMA_GAMMA_DEF
  rho <- params$RHO_DEF
  m_ex <- params$M_EX_HIST

  # Define DGPs for each arm
  arms <- list(
    oracle = list(
      sigma = sigma_theta_true, m_ex = 0,
      sg = 0, rho = 0,
      heritable = FALSE, individual = FALSE,
      seed = seed_base
    ),
    arm1_null = list(
      sigma = sigma_theta_true, m_ex = m_ex,
      sg = 0, rho = 0,
      heritable = FALSE, individual = FALSE,
      seed = seed_base + 1000L
    ),
    arm2_true_dgp = list(
      sigma = sigma_theta_true, m_ex = m_ex,
      sg = sg, rho = rho,
      heritable = TRUE, individual = TRUE,
      seed = seed_base + 2000L
    ),
    arm2_misspec_fit = list(
      sigma = sigma_theta_fit, m_ex = m_ex,
      sg = 0, rho = 0,
      heritable = FALSE, individual = FALSE,
      seed = seed_base + 3000L
    )
  )

  if (is.finite(sigma_theta_fix)) {
    arms$recovery <- list(
      sigma = sigma_theta_fix, m_ex = m_ex,
      sg = sg, rho = rho,
      heritable = TRUE, individual = TRUE,
      seed = seed_base + 4000L
    )
  }

  transforms <- c("none", "log", "rank")
  results <- list()

  for (arm_name in names(arms)) {
    arm <- arms[[arm_name]]
    data <- sim_twin_lifespans(
      arm$sigma, arm$m_ex,
      sigma_gamma = arm$sg, rho = arm$rho,
      gamma_heritable = arm$heritable,
      individual_ext = arm$individual,
      rng_seed = arm$seed
    )

    mz1 <- data$L_mz1[data$keep_mz]
    mz2 <- data$L_mz2[data$keep_mz]
    dz1 <- data$L_dz1[data$keep_dz]
    dz2 <- data$L_dz2[data$keep_dz]

    for (tf in transforms) {
      ace <- fit_ace_model(mz1, mz2, dz1, dz2, transform = tf)
      ace$arm <- arm_name
      ace$h2_falconer <- data$h2
      ace$r_mz <- data$r_mz
      ace$r_dz <- data$r_dz
      results[[paste0(arm_name, "_ace_", tf)]] <- ace

      ae <- fit_ae_model(mz1, mz2, dz1, dz2, transform = tf)
      ae$arm <- arm_name
      ae$h2_falconer <- data$h2
      results[[paste0(arm_name, "_ae_", tf)]] <- ae
    }
  }

  # Build summary data frame (ACE models only)
  ace_keys <- grep("_ace_", names(results), value = TRUE)
  ace_df <- do.call(rbind, lapply(ace_keys, function(nm) {
    r <- results[[nm]]
    data.frame(
      arm = r$arm, transform = r$transform, model_type = "ACE",
      a2 = r$a2, c2 = r$c2, e2 = r$e2,
      h2_ace = r$h2_ace, h2_falconer = r$h2_falconer,
      r_mz = r$r_mz, r_dz = r$r_dz,
      converged = r$converged, minus2ll = r$minus2ll,
      stringsAsFactors = FALSE
    )
  }))

  ae_keys <- grep("_ae_", names(results), value = TRUE)
  ae_df <- do.call(rbind, lapply(ae_keys, function(nm) {
    r <- results[[nm]]
    data.frame(
      arm = r$arm, transform = r$transform, model_type = "AE",
      a2 = r$a2, c2 = NA_real_, e2 = r$e2,
      h2_ace = r$h2_ae, h2_falconer = r$h2_falconer,
      r_mz = NA_real_, r_dz = NA_real_,
      converged = r$converged, minus2ll = r$minus2ll,
      stringsAsFactors = FALSE
    )
  }))

  summary <- rbind(ace_df, ae_df)
  rownames(summary) <- NULL

  list(results = results, summary = summary)
}
