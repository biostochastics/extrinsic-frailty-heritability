# ===========================================================================
# export_estimands.R — Export pipeline for alternative estimand analyses
# ===========================================================================
# Standalone export (does not modify existing tables/scalars.json).
# Produces: variance_components.csv, sweep_variance_decomp.csv,
#           bivariate_metrics.csv, ace_results.csv, age_dependence_cor.csv,
#           estimand_scalars.json
# ===========================================================================

# Null-coalesce operator
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Write all alternative estimand tables and scalars
#'
#' @return Invisible output directory path
write_estimand_tables <- function(variance_components,
                                   sweep_variance_decomp,
                                   variance_component_uq,
                                   bivariate_surface,
                                   ace_analysis,
                                   age_dependence) {
  outdir <- "tables"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # --- Variance components table ---
  write.csv(variance_components,
            file.path(outdir, "variance_components.csv"),
            row.names = FALSE)

  # --- Sweep variance decomposition ---
  write.csv(sweep_variance_decomp,
            file.path(outdir, "sweep_variance_decomp.csv"),
            row.names = FALSE)

  # --- Bivariate surface metrics ---
  biv_metrics <- data.frame(
    metric = character(), zygosity = character(),
    model = character(), value = numeric(),
    stringsAsFactors = FALSE
  )
  for (zyg in c("mz", "dz")) {
    for (stat in c("ise", "sup", "cvm", "iae")) {
      for (mdl in c("misspec", "fix")) {
        key <- paste0(stat, "_", mdl, "_", zyg)
        val <- bivariate_surface[[key]]
        if (!is.null(val)) {
          biv_metrics <- rbind(biv_metrics, data.frame(
            metric = stat, zygosity = toupper(zyg),
            model = mdl, value = val,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  write.csv(biv_metrics, file.path(outdir, "bivariate_metrics.csv"),
            row.names = FALSE)

  # --- ACE results ---
  write.csv(ace_analysis$summary,
            file.path(outdir, "ace_results.csv"),
            row.names = FALSE)

  # --- Age-specific dependence ---
  if (!is.null(age_dependence$correlations)) {
    write.csv(age_dependence$correlations,
              file.path(outdir, "age_dependence_cor.csv"),
              row.names = FALSE)
  }

  # --- Estimand scalars JSON ---
  scalars <- list()

  # Variance component inflation (GM misspecified)
  vc_gm_mis <- variance_components[
    variance_components$model == "GM" &
    grepl("Misspecified", variance_components$condition), ]
  if (nrow(vc_gm_mis) > 0) {
    scalars$sigma_infl_pct_gm <- round(vc_gm_mis$sigma_infl_pct[1], 1)
    scalars$var_infl_pct_gm <- round(vc_gm_mis$var_infl_pct[1], 1)
    scalars$delta_sigma_gm <- round(vc_gm_mis$delta_sigma[1], 4)
    scalars$sigma_true_gm <- round(vc_gm_mis$sigma_true[1], 4)
    scalars$sigma_fit_gm <- round(vc_gm_mis$sigma_fit[1], 4)
  }

  # Recovery
  vc_gm_fix <- variance_components[
    variance_components$model == "GM" &
    grepl("Recovery", variance_components$condition), ]
  if (nrow(vc_gm_fix) > 0) {
    scalars$sigma_infl_pct_fix <- round(vc_gm_fix$sigma_infl_pct[1], 1)
  }

  # MGG/SR misspecified
  for (mdl in c("MGG", "SR")) {
    row <- variance_components[
      variance_components$model == mdl &
      grepl("Misspecified", variance_components$condition), ]
    if (nrow(row) > 0) {
      prefix <- tolower(mdl)
      scalars[[paste0("sigma_infl_pct_", prefix)]] <-
        round(row$sigma_infl_pct[1], 1)
    }
  }

  # UQ
  uq_s <- variance_component_uq$summary
  si_row <- uq_s[uq_s$metric == "sigma_infl_arm2_pct", ]
  if (nrow(si_row) > 0) {
    scalars$sigma_infl_uq_mean <- round(si_row$mean[1], 1)
    scalars$sigma_infl_uq_se <- round(si_row$se[1], 2)
    scalars$sigma_infl_uq_lo95 <- round(si_row$lo95[1], 1)
    scalars$sigma_infl_uq_hi95 <- round(si_row$hi95[1], 1)
  }

  # Bivariate surface metrics
  scalars$biv_ise_misspec_mz <-
    round(bivariate_surface$ise_misspec_mz %||% NA, 6)
  scalars$biv_sup_misspec_mz <-
    round(bivariate_surface$sup_misspec_mz %||% NA, 4)
  scalars$biv_iae_misspec_mz <-
    round(bivariate_surface$iae_misspec_mz %||% NA, 4)
  if (!is.null(bivariate_surface$ise_fix_mz)) {
    scalars$biv_ise_fix_mz <- round(bivariate_surface$ise_fix_mz, 6)
    scalars$biv_ise_ratio_mz <- round(
      bivariate_surface$ise_misspec_mz /
      max(bivariate_surface$ise_fix_mz, 1e-10), 1)
  }

  # ACE (rank-normalized, converged)
  ace_rank <- ace_analysis$summary[
    ace_analysis$summary$transform == "rank" &
    ace_analysis$summary$model_type == "ACE" &
    ace_analysis$summary$converged, ]
  for (arm in unique(ace_rank$arm)) {
    row <- ace_rank[ace_rank$arm == arm, ]
    if (nrow(row) > 0) {
      scalars[[paste0("ace_", arm, "_a2")]] <- round(row$a2[1], 3)
      scalars[[paste0("ace_", arm, "_c2")]] <- round(row$c2[1], 3)
      scalars[[paste0("ace_", arm, "_e2")]] <- round(row$e2[1], 3)
      scalars[[paste0("ace_", arm, "_h2")]] <- round(row$h2_ace[1], 3)
    }
  }

  # Age dependence summary
  for (nm in names(age_dependence$summary)) {
    val <- age_dependence$summary[[nm]]
    if (is.numeric(val)) {
      scalars[[paste0("agedep_", nm)]] <- round(val, 3)
    }
  }

  write_estimand_scalars(scalars, outdir)

  invisible(outdir)
}

#' Write estimand scalars to JSON (separate from main scalars.json)
#'
#' @param scalars Named list
#' @param outdir Output directory
write_estimand_scalars <- function(scalars, outdir = "tables") {
  json_path <- file.path(outdir, "estimand_scalars.json")
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
