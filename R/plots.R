# ===========================================================================
# plots.R — All figure generation
# ===========================================================================
# Uses ggpubr/ggthemes with theme_light() base + refined palette.
# Viridis for heatmaps.
# ===========================================================================

# Refined publication theme based on theme_light()
theme_pub <- function(base_size = 14, base_family = "") {
  ggplot2::theme_light(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.1),
        margin = ggplot2::margin(b = 6), color = "gray10"
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.82), color = "gray45",
        margin = ggplot2::margin(b = 12)
      ),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(
        size = ggplot2::rel(0.92), color = "gray25"
      ),
      axis.text = ggplot2::element_text(
        size = ggplot2::rel(0.85), color = "gray35"
      ),
      legend.title = ggplot2::element_text(
        size = ggplot2::rel(0.88), face = "bold"
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.82)),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      plot.margin = ggplot2::margin(14, 16, 10, 12),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Refined palette
PAL <- c(
  Oracle    = "#1B6B93",
  Null      = "#2E8B57",
  Biased    = "#C73E3A",
  ControlA  = "#7B68AE",
  ControlB  = "#5B8FA8",
  Pleiotropy = "#8B6BAE",
  Irrelevant = "#6B8EAE",
  Fix       = "#3D9970",
  Cutoff0   = "#A0A0A0",
  Hamilton  = "#D68910"
)

# Categorical line palette (colorblind-safe)
LINE_PAL <- c("#1B6B93", "#2E8B57", "#C73E3A", "#D68910", "#7B68AE")

# ===================================================================
# Figure 1: Main arms bar chart
# ===================================================================

#' Bar chart of h² across all experimental arms
#'
#' @param summary_table Data frame from build_summary_table()
#' @return ggplot object
plot_main_arms <- function(summary_table) {
  df <- summary_table
  df <- df[!is.na(df$h2), ]
  df$Condition <- factor(df$Condition, levels = rev(df$Condition))

  # Semantic color mapping by condition name
  condition_colors <- setNames(
    unname(PAL[c("Oracle", "Null", "Biased", "ControlA", "Irrelevant",
                 "ControlB", "Pleiotropy", "Fix", "Cutoff0", "Hamilton")]),
    c("Oracle (intrinsic only)", "Arm 1: Correctly specified",
      "Arm 2: Misspecification", "Control A: Non-heritable extrinsic",
      "Irrelevant trait (zeta)", "Control B: Vanishing m_ex",
      "Pleiotropy isolation (rho=0)", "Oracle fix: correct model",
      "Cutoff = 0", "Arm 3: Hamilton conditioning")
  )

  oracle_h2 <- df$h2[grepl("Oracle.*intrinsic", df$Condition)]

  ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = h2, fill = Condition)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = oracle_h2,
                        linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%.3f\n(%+.1f pp)", h2, bias * 100)),
      hjust = -0.05, size = 3.0, color = "gray30") +
    ggplot2::coord_flip(ylim = c(0, max(df$h2, na.rm = TRUE) * 1.25)) +
    ggplot2::scale_fill_manual(values = condition_colors, guide = "none") +
    ggplot2::labs(
      title = "Falconer h\u00b2 across experimental arms",
      subtitle = "Dashed line = oracle (true intrinsic-only) h\u00b2",
      x = NULL, y = expression(hat(h)^2)
    ) +
    theme_pub()
}

# ===================================================================
# Figure 2: Sensitivity sweep heatmap
# ===================================================================

#' Heatmap of bias across (rho, sigma_gamma) sweep
#'
#' @param sweep_df Data frame from sweep_results target
#' @return ggplot object
plot_sweep_heatmap <- function(sweep_df) {
  df <- sweep_df
  df$bias_pp <- 100 * df$bias

  ggplot2::ggplot(df, ggplot2::aes(x = sigma_gamma, y = rho, fill = bias_pp)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(
      name = "Bias (pp)", option = "viridis", direction = -1
    ) +
    ggplot2::coord_fixed(
      ratio = diff(range(df$sigma_gamma)) / diff(range(df$rho))
    ) +
    ggplot2::labs(
      title = "Heritability bias (pp) across parameter space",
      subtitle = expression(paste("Bias = ", hat(h)^2 - h[oracle]^2,
                                  ", pooled over ", sigma[gamma], " and ", rho)),
      x = expression(sigma[gamma]), y = expression(rho)
    ) +
    theme_pub() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# ===================================================================
# Figure 3: Anchored sweep heatmap
# ===================================================================

#' Heatmap of bias in the anchored parameter region
#'
#' @param anchored_df Data frame from anchored_results target
#' @return ggplot object
plot_anchored_heatmap <- function(anchored_df) {
  df <- anchored_df
  df$bias_pp <- 100 * df$bias

  ggplot2::ggplot(df, ggplot2::aes(x = sigma_gamma, y = rho, fill = bias_pp)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(
      name = "Bias (pp)", option = "viridis", direction = -1
    ) +
    ggplot2::coord_fixed(
      ratio = diff(range(df$sigma_gamma)) / diff(range(df$rho))
    ) +
    ggplot2::labs(
      title = "Anchored region: heritability bias (pp)",
      subtitle = expression(paste(
        sigma[gamma] %in% group("[", list(0.30, 0.65), "]"),
        ", ", rho %in% group("[", list(0.20, 0.50), "]")
      )),
      x = expression(sigma[gamma]), y = expression(rho)
    ) +
    theme_pub() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# ===================================================================
# Figure 4: Model validation comparison
# ===================================================================

#' Grouped bar chart comparing bias across GM, MGG, SR models
#'
#' @param model_table Data frame from build_model_table()
#' @return ggplot object
plot_model_comparison <- function(model_table) {
  df_long <- data.frame(
    Model = rep(model_table$Model, 2),
    Arm = rep(c("Arm 1 (null)", "Arm 2 (misspec)"), each = nrow(model_table)),
    bias_pp = c(model_table$Arm1_bias_pp, model_table$Arm2_bias_pp),
    stringsAsFactors = FALSE
  )
  df_long$Model <- factor(df_long$Model, levels = model_table$Model)

  ggplot2::ggplot(df_long, ggplot2::aes(x = Model, y = bias_pp, fill = Arm)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.1f", bias_pp)),
      position = ggplot2::position_dodge(width = 0.7),
      vjust = -0.5, size = 3.5
    ) +
    ggplot2::scale_fill_manual(
      values = c("Arm 1 (null)" = PAL["Null"], "Arm 2 (misspec)" = PAL["Biased"])
    ) +
    ggplot2::labs(
      title = "Bias across mortality models",
      subtitle = "All three models show positive bias under misspecification",
      x = NULL, y = "Bias (percentage points)",
      fill = NULL
    ) +
    theme_pub() +
    ggplot2::theme(legend.position = "top")
}

# ===================================================================
# Figure 5: Controls bar chart
# ===================================================================

#' Bar chart of control experiment biases
#'
#' @param controls_table Data frame from controls_table target
#' @return ggplot object
plot_controls <- function(controls_table) {
  df <- controls_table
  df$Control <- factor(df$Control, levels = rev(df$Control))

  ggplot2::ggplot(df, ggplot2::aes(x = Control, y = bias_pp, fill = bias_pp > 0)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f pp", bias_pp)),
                       hjust = ifelse(df$bias_pp > 0, -0.1, 1.1),
                       size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("TRUE" = PAL["Biased"],
                                          "FALSE" = PAL["Oracle"]),
                               guide = "none") +
    ggplot2::labs(
      title = "Control experiments: bias (pp)",
      subtitle = "Only heritable, mortality-relevant, pleiotropic traits produce substantial bias",
      x = NULL, y = "Bias (percentage points)"
    ) +
    theme_pub()
}

# ===================================================================
# Diagnostic: Dose-response
# ===================================================================

#' Bias vs m_ex dose-response curve
#'
#' @param dose_response_df Data frame from dose_response target
#' @return ggplot object
plot_dose_response <- function(dose_response_df) {
  ggplot2::ggplot(dose_response_df,
                  ggplot2::aes(x = m_ex * 1000, y = bias * 100)) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      title = "Dose-response: bias scales with extrinsic mortality",
      subtitle = expression(paste(sigma[gamma], " = 0.40, ", rho, " = 0.4")),
      x = expression(m[ex] %*% 10^3),
      y = "Bias (pp)"
    ) +
    theme_pub()
}

# ===================================================================
# Diagnostic: Negative rho sensitivity
# ===================================================================

#' Bias vs rho showing sign reversal
#'
#' @param neg_rho_df Data frame from negative_rho target
#' @return ggplot object
plot_negative_rho <- function(neg_rho_df) {
  ggplot2::ggplot(neg_rho_df, ggplot2::aes(x = rho, y = bias_pp)) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::labs(
      title = expression(paste("Bias vs ", rho, ": sign reversal at ", rho, " < 0")),
      subtitle = "Falsifiable prediction: negative pleiotropy causes downward bias",
      x = expression(rho ~ (genetic ~ correlation)),
      y = "Bias (pp)"
    ) +
    theme_pub()
}

# ===================================================================
# Diagnostic: m_ex split sensitivity
# ===================================================================

#' Bias vs fraction of m_ex that is heritable
#'
#' @param mex_split_df Data frame from mex_split target
#' @return ggplot object
plot_mex_split <- function(mex_split_df) {
  ggplot2::ggplot(mex_split_df,
                  ggplot2::aes(x = frac_heritable, y = bias_pp)) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      title = "Bias vs heritable fraction of extrinsic mortality",
      subtitle = expression(paste(m[ex], " = ", m[inf], " + ", m[other],
                                  "; only ", m[inf], " gets ", gamma)),
      x = expression(f[heritable] == m[inf] / m[ex]),
      y = "Bias (pp)"
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_pub()
}

# ===================================================================
# Diagnostic: sigma_gamma bridge
# ===================================================================

#' Bridge curve: h²_infection vs sigma_gamma
#'
#' @param bridge_result List from sigma_gamma_bridge target
#' @return ggplot object
plot_sigma_gamma_bridge <- function(bridge_result) {
  df <- bridge_result$curve
  bridge_sg <- bridge_result$bridge_sigma_gamma

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sigma_gamma, y = h2_infection)) +
    ggplot2::geom_line(color = PAL["Oracle"], linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0.40, linetype = "dashed",
                        color = PAL["Biased"], linewidth = 0.6) +
    ggplot2::annotate("text", x = 0.2, y = 0.42,
                      label = "Target: h\u00b2_L = 0.40",
                      color = PAL["Biased"], size = 3.5, hjust = 0)

  if (!is.na(bridge_sg)) {
    p <- p +
      ggplot2::geom_vline(xintercept = bridge_sg, linetype = "dotted",
                          color = PAL["Null"], linewidth = 0.6) +
      ggplot2::annotate("text", x = bridge_sg + 0.05, y = 0.1,
                        label = sprintf("\u03c3_\u03b3 = %.2f", bridge_sg),
                        color = PAL["Null"], size = 3.5, hjust = 0)
  }

  p +
    ggplot2::labs(
      title = expression(paste(sigma[gamma], " bridge: infection-death heritability")),
      subtitle = "Tetrachoric MLE via bivariate normal threshold model",
      x = expression(sigma[gamma]),
      y = expression(hat(h)[L]^2 ~ "(infection death)")
    ) +
    theme_pub()
}

# ===================================================================
# Diagnostic: MGG hazard curves
# ===================================================================

#' MGG intrinsic hazard curves under two SM mappings
#'
#' @param hazard_df Data frame from mgg_hazard_curves target
#' @return ggplot object
plot_mgg_hazard_curves <- function(hazard_df) {
  t_star <- attr(hazard_df, "t_star")
  df <- hazard_df[hazard_df$t > 30 & hazard_df$t < 110, ]
  df$mapping_label <- ifelse(df$mapping == "compensatory",
                             "Correct: a = a0^q",
                             "Paper: a = a0^(1/q)")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = hazard,
                                         color = q, linetype = mapping_label)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = LINE_PAL[1:3]) +
    ggplot2::facet_wrap(~ mapping_label, scales = "free_y")

  if (!is.null(t_star) && !is.na(t_star)) {
    p <- p + ggplot2::geom_vline(xintercept = t_star, linetype = "dotted",
                                 color = "gray50")
  }

  p +
    ggplot2::labs(
      title = "MGG intrinsic hazard: compensatory vs paper parameterization",
      subtitle = expression(paste("SM convergence at t* = -ln(", a[0], ") / ", b[0])),
      x = "Age (years)", y = expression(mu[int](t) ~ "(log scale)"),
      color = "q", linetype = "Mapping"
    ) +
    theme_pub() +
    ggplot2::theme(legend.position = "bottom")
}

# ===================================================================
# MGG parameterization comparison (a=a0^q vs a=a0^{1/q})
# ===================================================================

#' Side-by-side MGG parameterization comparison
#'
#' Three-panel figure: (A) hazard curves overlaid, (B) survival curves,
#' (C) twin h² under each mapping. Directly demonstrates that the paper's
#' stated formula breaks SM compensation.
#'
#' @param comparison List from run_mgg_param_comparison()
#' @return ggplot (patchwork composite)
plot_mgg_param_comparison <- function(comparison) {
  hazard_df <- comparison$hazard_df
  survival_df <- comparison$survival_df
  twin_stats <- comparison$twin_stats
  demo_stats <- comparison$demo_stats
  t_star <- comparison$t_star

  mapping_labels <- c(
    "compensatory" = "Correct: a = a0^q",
    "paper" = "Paper: a = a0^(1/q)")

  mapping_colors <- c(
    "compensatory" = "#2166AC",
    "paper" = "#B2182B")

  # --- Panel A: Hazard curves overlaid ---
  hdf <- hazard_df[hazard_df$t > 30 & hazard_df$t < 110, ]
  hdf$mapping_label <- mapping_labels[hdf$mapping]

  p_hazard <- ggplot2::ggplot(hdf, ggplot2::aes(
      x = t, y = hazard, color = mapping,
      linetype = q, group = interaction(mapping, q))) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(
      values = mapping_colors, labels = mapping_labels) +
    ggplot2::scale_linetype_manual(
      values = c("0.7" = "dashed", "1" = "solid", "1.3" = "dotdash")) +
    ggplot2::geom_vline(xintercept = t_star, linetype = "dotted",
                        color = "gray50", linewidth = 0.5) +
    ggplot2::annotate("text", x = t_star + 1, y = 0.003,
      label = paste0("t* = ", round(t_star)), hjust = 0, size = 3) +
    ggplot2::labs(
      title = "A. Intrinsic hazard curves",
      x = "Age (years)", y = expression(mu[int](t)),
      color = "SM mapping", linetype = "q") +
    theme_pub() +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical",
                   legend.margin = ggplot2::margin(t = -5))

  # --- Panel B: Survival curves ---
  sdf <- survival_df
  sdf$mapping_label <- mapping_labels[sdf$mapping]

  mean_comp <- demo_stats$mean_age[demo_stats$mapping == "compensatory"]
  mean_paper <- demo_stats$mean_age[demo_stats$mapping == "paper"]

  p_surv <- ggplot2::ggplot(sdf, ggplot2::aes(
      x = age, y = survival, color = mapping)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(
      values = mapping_colors, labels = mapping_labels) +
    ggplot2::geom_vline(xintercept = mean_comp, linetype = "dashed",
                        color = mapping_colors["compensatory"], linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = mean_paper, linetype = "dashed",
                        color = mapping_colors["paper"], linewidth = 0.4) +
    ggplot2::annotate("text", x = mean_comp + 1, y = 0.85,
      label = sprintf("mean = %.0f", mean_comp),
      color = mapping_colors["compensatory"], hjust = 0, size = 3) +
    ggplot2::annotate("text", x = mean_paper + 1, y = 0.75,
      label = sprintf("mean = %.0f", mean_paper),
      color = mapping_colors["paper"], hjust = 0, size = 3) +
    ggplot2::labs(
      title = "B. Cohort survival (intrinsic only)",
      x = "Age (years)", y = "S(t)",
      color = "SM mapping") +
    theme_pub() +
    ggplot2::theme(legend.position = "bottom")

  # --- Panel C: Twin h² comparison ---
  tdf <- twin_stats
  tdf$mapping_label <- factor(mapping_labels[tdf$mapping],
    levels = mapping_labels)

  # Reshape for grouped bars
  tdf_long <- rbind(
    data.frame(mapping_label = tdf$mapping_label,
               metric = "r_MZ", value = tdf$r_mz,
               stringsAsFactors = FALSE),
    data.frame(mapping_label = tdf$mapping_label,
               metric = "r_DZ", value = tdf$r_dz,
               stringsAsFactors = FALSE),
    data.frame(mapping_label = tdf$mapping_label,
               metric = "Falconer h\u00B2", value = tdf$h2,
               stringsAsFactors = FALSE)
  )
  tdf_long$metric <- factor(tdf_long$metric,
    levels = c("r_MZ", "r_DZ", "Falconer h\u00B2"))

  p_twin <- ggplot2::ggplot(tdf_long, ggplot2::aes(
      x = metric, y = value, fill = mapping_label)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
      position = ggplot2::position_dodge(width = 0.7),
      vjust = -0.5, size = 3) +
    ggplot2::scale_fill_manual(values = setNames(
      mapping_colors, mapping_labels)) +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed",
                        color = "gray50", linewidth = 0.4) +
    ggplot2::annotate("text", x = 3.3, y = 1.03, label = "h\u00B2 = 1",
                      size = 2.8, color = "gray40") +
    ggplot2::labs(
      title = "C. Twin correlations & h\u00B2",
      x = NULL, y = "Value", fill = "SM mapping") +
    ggplot2::ylim(NA, max(tdf_long$value) * 1.15) +
    theme_pub() +
    ggplot2::theme(legend.position = "bottom")

  # Combine with patchwork
  patchwork::wrap_plots(p_hazard, p_surv, p_twin, ncol = 3) +
    patchwork::plot_annotation(
      title = "MGG Strehler-Mildvan parameterization: compensatory (a = a0^q) vs paper (a = a0^{1/q})")
}

# ===================================================================
# Diagnostic: Variance decomposition
# ===================================================================

#' Predicted vs actual variance scatter
#'
#' @param var_decomp List from variance_decomposition()
#' @param alpha_beta List from estimate_alpha_beta()
#' @return ggplot object
plot_variance_decomp <- function(var_decomp, alpha_beta) {
  df <- var_decomp$data

  lbl <- sprintf("alpha = %.2f, beta = %.2f\nR^2 = %.3f\nMAR = %.1f%%",
                 alpha_beta$alpha, alpha_beta$beta,
                 alpha_beta$r_squared, var_decomp$mar_pct)

  ggplot2::ggplot(df, ggplot2::aes(x = predicted_var, y = actual_var)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "gray60") +
    ggplot2::geom_point(ggplot2::aes(color = factor(rho)), size = 2.5) +
    ggplot2::scale_color_manual(
      values = viridis::viridis(length(unique(df$rho)))
    ) +
    ggplot2::annotate("text", x = min(df$predicted_var),
                      y = max(df$actual_var) * 0.95,
                      label = lbl, hjust = 0, size = 3.5,
                      family = "mono") +
    ggplot2::labs(
      title = expression(paste("Variance decomposition: ",
        sigma[fit]^2, " vs ", sigma[true]^2 + alpha * sigma[gamma]^2 +
          2 * beta * rho * sigma[theta] * sigma[gamma])),
      x = expression("Predicted " ~ sigma[fit]^2),
      y = expression("Actual " ~ sigma[fit]^2),
      color = expression(rho)
    ) +
    theme_pub()
}

# ===================================================================
# Diagnostic: Joint r_MZ / r_DZ
# ===================================================================

#' Joint calibration diagnostic scatter
#'
#' @param joint_diag List from run_joint_rmz_rdz_diagnostic()
#' @return ggplot object
plot_joint_diagnostic <- function(joint_diag) {
  df <- data.frame(
    Source = c("True DGP", "Misspec (r_MZ-only)", "Joint calibration"),
    r_MZ = c(joint_diag$true_r_mz, joint_diag$misspec_pred_r_mz, NA),
    r_DZ = c(joint_diag$true_r_dz, joint_diag$misspec_pred_r_dz, NA),
    h2 = c(joint_diag$true_h2, joint_diag$misspec_h2, joint_diag$joint_h2),
    stringsAsFactors = FALSE
  )

  lbl <- sprintf(
    "r_DZ discrepancy = %.4f\nJoint h2 = %.3f (bias = %.3f)",
    joint_diag$r_dz_discrepancy, joint_diag$joint_h2, joint_diag$joint_bias
  )

  ggplot2::ggplot(df[1:2, ], ggplot2::aes(x = r_MZ, y = r_DZ, color = Source)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_segment(
      x = df$r_MZ[1], xend = df$r_MZ[2],
      y = df$r_DZ[1], yend = df$r_DZ[2],
      arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches")),
      color = "gray50", linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(values = c(PAL["Oracle"], PAL["Biased"])) +
    ggplot2::annotate("text", x = min(df$r_MZ[1:2]),
                      y = max(df$r_DZ[1:2]) + 0.01,
                      label = lbl, hjust = 0, size = 3.2,
                      family = "mono") +
    ggplot2::labs(
      title = "Joint (r_MZ, r_DZ) calibration diagnostic",
      subtitle = "Arrow: true DGP -> misspecified model prediction",
      x = expression(r[MZ]), y = expression(r[DZ]),
      color = NULL
    ) +
    theme_pub()
}

# ===================================================================
# Diagnostic: SR dt sensitivity
# ===================================================================

#' SR model dt sensitivity bar chart
#'
#' @param sr_dt_df Data frame from sr_dt_sensitivity target
#' @return ggplot object
plot_sr_dt_sensitivity <- function(sr_dt_df) {
  ggplot2::ggplot(sr_dt_df, ggplot2::aes(x = dt_label, y = arm2_bias_pp)) +
    ggplot2::geom_col(fill = PAL["Biased"], width = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f pp", arm2_bias_pp)),
                       vjust = -0.5, size = 3.8) +
    ggplot2::labs(
      title = "SR model: Arm 2 bias stability across time steps",
      subtitle = "Euler-Maruyama convergence check",
      x = expression(Delta * t ~ "(fraction of year)"),
      y = "Arm 2 bias (pp)"
    ) +
    theme_pub()
}

# ===================================================================
# Multi-model comparison: Negative rho across GM, MGG, SR
# ===================================================================

#' Multi-model negative rho comparison
#'
#' @param neg_rho_gm GM negative rho data frame
#' @param model_controls Combined model controls data frame
#' @return ggplot object
plot_negative_rho_multimodel <- function(neg_rho_gm, model_controls) {
  gm_df <- data.frame(
    model = "GM", rho = neg_rho_gm$rho,
    bias_pp = neg_rho_gm$bias_pp, stringsAsFactors = FALSE
  )
  mc_rho <- model_controls[model_controls$sweep_type == "negative_rho", ]
  mc_df <- data.frame(
    model = toupper(mc_rho$model),
    rho = mc_rho$sweep_val,
    bias_pp = mc_rho$bias_pp, stringsAsFactors = FALSE
  )
  df <- rbind(gm_df, mc_df)
  df$model <- factor(df$model, levels = c("GM", "MGG", "SR"))

  ggplot2::ggplot(df, ggplot2::aes(x = rho, y = bias_pp,
                                    color = model, shape = model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::scale_color_manual(values = c(GM = PAL["Oracle"],
                                           MGG = PAL["Biased"],
                                           SR = PAL["Hamilton"])) +
    ggplot2::labs(
      title = expression(paste("Bias vs ", rho, " across all three models")),
      subtitle = "All models show sign reversal — falsifiable prediction holds universally",
      x = expression(rho ~ (genetic ~ correlation)),
      y = "Bias (pp)", color = "Model", shape = "Model"
    ) +
    theme_pub() +
    ggplot2::theme(legend.position = c(0.15, 0.85))
}

# ===================================================================
# Multi-model comparison: Dose-response across GM, MGG, SR
# ===================================================================

#' Multi-model dose-response comparison
#'
#' @param dose_response_gm GM dose-response data frame
#' @param model_controls Combined model controls data frame
#' @return ggplot object
plot_dose_response_multimodel <- function(dose_response_gm, model_controls) {
  gm_df <- data.frame(
    model = "GM", m_ex = dose_response_gm$m_ex,
    bias_pp = dose_response_gm$bias * 100, stringsAsFactors = FALSE
  )
  mc_dr <- model_controls[model_controls$sweep_type == "dose_response", ]
  mc_df <- data.frame(
    model = toupper(mc_dr$model),
    m_ex = mc_dr$sweep_val,
    bias_pp = mc_dr$bias_pp, stringsAsFactors = FALSE
  )
  df <- rbind(gm_df, mc_df)
  df$model <- factor(df$model, levels = c("GM", "MGG", "SR"))

  ggplot2::ggplot(df, ggplot2::aes(x = m_ex * 1000, y = bias_pp,
                                    color = model, shape = model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = c(GM = PAL["Oracle"],
                                           MGG = PAL["Biased"],
                                           SR = PAL["Hamilton"])) +
    ggplot2::labs(
      title = "Dose-response: bias scales with extrinsic mortality (all models)",
      subtitle = expression(paste(sigma[gamma], " = 0.40, ", rho, " = 0.4")),
      x = expression(m[ex] %*% 10^3),
      y = "Bias (pp)", color = "Model", shape = "Model"
    ) +
    theme_pub() +
    ggplot2::theme(legend.position = c(0.15, 0.85))
}

# ===================================================================
# Multi-model comparison: Pleiotropy isolation
# ===================================================================

#' Multi-model pleiotropy isolation bar chart
#'
#' @param controls_table GM controls table
#' @param model_controls Combined model controls data frame
#' @return ggplot object
plot_pleiotropy_multimodel <- function(controls_table, model_controls) {
  mc_pl <- model_controls[model_controls$sweep_type == "pleiotropy_isolation", ]

  # GM arm2 reference and rho=0 from controls_table
  gm_arm2 <- controls_table$bias_pp[controls_table$Control == "Arm 2 reference"]
  gm_rho0 <- controls_table$bias_pp[controls_table$Control == "Pleiotropy isolation (rho=0)"]

  df <- data.frame(
    Model = c("GM", "GM", "MGG", "MGG", "SR", "SR"),
    Condition = rep(c("Arm 2 (rho=0.4)", "rho=0"), 3),
    bias_pp = c(
      gm_arm2, gm_rho0,
      mc_pl$bias_pp[mc_pl$model == "mgg"], NA,
      mc_pl$bias_pp[mc_pl$model == "sr"], NA
    ),
    stringsAsFactors = FALSE
  )

  # Fill in MGG/SR arm2 from model_table (not available here, use combined)
  # Actually use the model_controls for a simpler approach
  # We only have rho=0 from model_controls for MGG/SR
  # Let's just show the rho=0 isolation result across models
  df_simple <- data.frame(
    Model = c("GM", "MGG", "SR"),
    bias_pp_rho0 = c(gm_rho0,
                     mc_pl$bias_pp[mc_pl$model == "mgg"],
                     mc_pl$bias_pp[mc_pl$model == "sr"]),
    stringsAsFactors = FALSE
  )
  df_simple$Model <- factor(df_simple$Model, levels = c("GM", "MGG", "SR"))

  ggplot2::ggplot(df_simple, ggplot2::aes(x = Model, y = bias_pp_rho0,
                                           fill = Model)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f pp", bias_pp_rho0)),
                       vjust = -0.5, size = 3.8) +
    ggplot2::scale_fill_manual(values = c(GM = PAL["Oracle"],
                                          MGG = PAL["Biased"],
                                          SR = PAL["Hamilton"]),
                               guide = "none") +
    ggplot2::labs(
      title = expression(paste("Pleiotropy isolation (", rho, "=0) across models")),
      subtitle = "Removing genetic correlation collapses bias in all models",
      x = NULL, y = "Bias (pp)"
    ) +
    theme_pub()
}

# ===================================================================
# Composite: All diagnostics panel
# ===================================================================

#' Composite panel of all diagnostic plots (3x3 or similar)
#'
#' @param plots Named list of ggplot objects
#' @return patchwork composite
plot_all_diagnostics <- function(plots) {
  # Arrange available plots in a grid
  available <- Filter(Negate(is.null), plots)
  patchwork::wrap_plots(available, ncol = 3) +
    patchwork::plot_annotation(
      title = "Diagnostic panel",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16)
      )
    )
}

# ===================================================================
# Save helper
# ===================================================================

#' Save a ggplot to figures/ directory
#'
#' @param p ggplot object
#' @param filename Filename (without path)
#' @param width Width in inches
#' @param height Height in inches
#' @param outdir Output directory
#' @return Invisible file path
save_figure <- function(p, filename, width = 10, height = 7,
                        outdir = "figures") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  path <- file.path(outdir, filename)
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 300,
                  bg = "white")
  invisible(path)
}

# ===================================================================
# Master figure generation target
# ===================================================================

#' Generate all figures from pipeline results
#'
#' @param summary_table Main arms table
#' @param sweep_results Sweep data frame
#' @param anchored_results Anchored sweep data frame
#' @param model_table Model comparison table
#' @param controls_table Controls table
#' @param dose_response Dose-response data frame
#' @param negative_rho Negative rho data frame
#' @param mex_split m_ex split data frame
#' @param sigma_gamma_bridge Bridge result list
#' @param mgg_hazard_curves Hazard curves data frame
#' @param var_decomp Variance decomposition list
#' @param empirical_alpha_beta Alpha/beta estimation list
#' @param joint_diagnostic Joint diagnostic list
#' @param sr_dt_sensitivity SR dt check data frame
#' @param outdir Output directory
#' @return Named list of file paths
generate_all_figures <- function(summary_table, sweep_results, anchored_results,
                                 model_table, controls_table,
                                 dose_response, negative_rho, mex_split,
                                 sigma_gamma_bridge, mgg_hazard_curves,
                                 var_decomp, empirical_alpha_beta,
                                 joint_diagnostic, sr_dt_sensitivity,
                                 model_controls = NULL,
                                 mgg_param_comparison = NULL,
                                 outdir = "figures") {
  paths <- list()

  # Main manuscript figures
  paths$main_arms <- save_figure(
    plot_main_arms(summary_table),
    "fig1_main_arms.png", width = 10, height = 6, outdir = outdir)

  paths$sweep <- save_figure(
    plot_sweep_heatmap(sweep_results),
    "fig2_sweep_heatmap.png", width = 9, height = 6, outdir = outdir)

  paths$anchored <- save_figure(
    plot_anchored_heatmap(anchored_results),
    "fig3_anchored_heatmap.png", width = 9, height = 6, outdir = outdir)

  paths$models <- save_figure(
    plot_model_comparison(model_table),
    "fig4_model_comparison.png", width = 9, height = 6, outdir = outdir)

  paths$controls <- save_figure(
    plot_controls(controls_table),
    "fig5_controls.png", width = 10, height = 5, outdir = outdir)

  # Diagnostic figures
  paths$dose_response <- save_figure(
    plot_dose_response(dose_response),
    "diag_dose_response.png", width = 8, height = 5.5, outdir = outdir)

  paths$negative_rho <- save_figure(
    plot_negative_rho(negative_rho),
    "diag_negative_rho.png", width = 8, height = 5.5, outdir = outdir)

  paths$mex_split <- save_figure(
    plot_mex_split(mex_split),
    "diag_mex_split.png", width = 8, height = 5.5, outdir = outdir)

  paths$bridge <- save_figure(
    plot_sigma_gamma_bridge(sigma_gamma_bridge),
    "diag_sigma_gamma_bridge.png", width = 8, height = 5.5, outdir = outdir)

  paths$mgg_hazard <- save_figure(
    plot_mgg_hazard_curves(mgg_hazard_curves),
    "diag_mgg_hazard_curves.png", width = 12, height = 6, outdir = outdir)

  if (!is.null(mgg_param_comparison)) {
    paths$mgg_param <- save_figure(
      plot_mgg_param_comparison(mgg_param_comparison),
      "fig_mgg_parameterization.png", width = 16, height = 6, outdir = outdir)
  }

  paths$var_decomp <- save_figure(
    plot_variance_decomp(var_decomp, empirical_alpha_beta),
    "diag_variance_decomp.png", width = 8, height = 7, outdir = outdir)

  paths$joint <- save_figure(
    plot_joint_diagnostic(joint_diagnostic),
    "diag_joint_rmz_rdz.png", width = 8, height = 6, outdir = outdir)

  paths$sr_dt <- save_figure(
    plot_sr_dt_sensitivity(sr_dt_sensitivity),
    "diag_sr_dt_sensitivity.png", width = 7, height = 5, outdir = outdir)

  # Multi-model comparison figures (if model_controls available)
  if (!is.null(model_controls) && nrow(model_controls) > 0) {
    paths$neg_rho_multimodel <- save_figure(
      plot_negative_rho_multimodel(negative_rho, model_controls),
      "fig6_negative_rho_multimodel.png", width = 9, height = 6, outdir = outdir)

    paths$dose_response_multimodel <- save_figure(
      plot_dose_response_multimodel(dose_response, model_controls),
      "fig7_dose_response_multimodel.png", width = 9, height = 6, outdir = outdir)

    paths$pleiotropy_multimodel <- save_figure(
      plot_pleiotropy_multimodel(controls_table, model_controls),
      "fig8_pleiotropy_multimodel.png", width = 7, height = 5.5, outdir = outdir)
  }

  paths
}
