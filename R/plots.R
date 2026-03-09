# ===========================================================================
# plots.R — All figure generation
# ===========================================================================
# Uses ggpubr/ggthemes with theme_light() base + refined palette.
# Viridis for heatmaps.
# ===========================================================================

# Refined publication theme based on theme_light()
# Helvetica Neue: clean sans-serif that stays legible when figures are
# composited at half-width in Typst grids (~3.5" effective on A4).
theme_pub <- function(base_size = 16, base_family = "Helvetica Neue") {
  ggplot2::theme_light(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.05),
        margin = ggplot2::margin(b = 6), color = "gray10"
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.78), color = "gray35",
        margin = ggplot2::margin(b = 12)
      ),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(
        size = ggplot2::rel(0.88), color = "gray20"
      ),
      axis.text = ggplot2::element_text(
        size = ggplot2::rel(0.82), color = "gray30"
      ),
      legend.title = ggplot2::element_text(
        size = ggplot2::rel(0.85), face = "bold"
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.80)),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.88)),
      plot.margin = ggplot2::margin(14, 16, 10, 12),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "gray92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "gray80", linewidth = 0.3)
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

# Model-type palette: consistent across all multi-model figures
MODEL_PAL       <- c(GM = "#1B6B93", MGG = "#C73E3A", SR = "#D68910")
MODEL_SHAPES    <- c(GM = 16, MGG = 17, SR = 15)
MODEL_LINETYPES <- c(GM = "solid", MGG = "dashed", SR = "dotted")

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
      hjust = -0.05, size = 3.8, color = "gray30") +
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
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
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
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
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
      vjust = -0.5, size = 4.2
    ) +
    ggplot2::scale_fill_manual(
      values = c("Arm 1 (null)" = unname(PAL["Null"]),
                 "Arm 2 (misspec)" = unname(PAL["Biased"]))
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
                       size = 4.2) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_fill_manual(values = c("TRUE" = unname(PAL["Biased"]),
                                          "FALSE" = unname(PAL["Oracle"])),
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
  # Anchored regime for rho
  rho_lo <- 0.20
  rho_hi <- 0.50

  ggplot2::ggplot(neg_rho_df, ggplot2::aes(x = rho, y = bias_pp)) +
    # Shaded region: negative pleiotropy (sign-reversal zone)
    ggplot2::annotate("rect",
                      xmin = min(neg_rho_df$rho) - 0.05, xmax = 0,
                      ymin = -Inf, ymax = Inf,
                      fill = PAL["Pleiotropy"], alpha = 0.08) +
    ggplot2::annotate("text", x = mean(c(min(neg_rho_df$rho), 0)), y = max(neg_rho_df$bias_pp) * 0.9,
                      label = expression("Negative " * rho * ": deflation"),
                      color = PAL["Pleiotropy"], size = 4.0,
                      fontface = "italic", lineheight = 0.9) +
    # Shaded region: anchored sensitivity range
    ggplot2::annotate("rect",
                      xmin = rho_lo, xmax = rho_hi,
                      ymin = -Inf, ymax = Inf,
                      fill = PAL["Oracle"], alpha = 0.08) +
    ggplot2::annotate("text", x = (rho_lo + rho_hi) / 2, y = max(neg_rho_df$bias_pp) * 0.9,
                      label = "Anchored\nregime",
                      color = PAL["Oracle"], size = 4.0,
                      fontface = "italic") +
    # Data
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50",
                        linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray50",
                        linewidth = 0.6) +
    ggplot2::labs(
      title = expression("Bias vs genetic correlation " * rho * ": sign reversal under negative " * rho),
      subtitle = expression(
        "At " * sigma[gamma] == 0.40 *
        "; negative " * rho * " reverses inflation to deflation"
      ),
      x = expression("Genetic correlation " * rho *
                      " between intrinsic and extrinsic susceptibility"),
      y = "Bias (percentage points)"
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

  # Sensitivity range used in the paper
  sg_lo <- 0.30
  sg_hi <- 0.65

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sigma_gamma, y = h2_infection)) +
    # Shaded band: our conservative sensitivity range
    ggplot2::annotate("rect", xmin = sg_lo, xmax = sg_hi,
                      ymin = -Inf, ymax = Inf,
                      fill = PAL["Oracle"], alpha = 0.10) +
    ggplot2::annotate("text", x = (sg_lo + sg_hi) / 2, y = 0.48,
                      label = "Our sensitivity\nrange",
                      color = PAL["Oracle"], size = 4.0,
                      fontface = "italic", lineheight = 0.9) +
    # Raw simulated curve (no smoothing)
    ggplot2::geom_line(color = PAL["Oracle"], linewidth = 0.9) +
    # Obel target line
    ggplot2::geom_hline(yintercept = 0.40, linetype = "dashed",
                        color = PAL["Biased"], linewidth = 0.7) +
    ggplot2::annotate("text", x = 0.15, y = 0.425,
                      label = expression(Obel~et~al.~target:~h[L]^2 == 0.40),
                      color = PAL["Biased"], size = 4.2, hjust = 0)

  if (!is.na(bridge_sg)) {
    p <- p +
      ggplot2::geom_vline(xintercept = bridge_sg, linetype = "dotted",
                          color = PAL["Null"], linewidth = 0.7) +
      ggplot2::annotate("text", x = bridge_sg - 0.05, y = 0.12,
                        label = sprintf("\u03c3\u03b3 \u2248 %.2f", bridge_sg),
                        color = PAL["Null"], size = 4.5, hjust = 1,
                        fontface = "bold")
  }

  p +
    ggplot2::scale_x_continuous(
      limits = c(0.1, max(df$sigma_gamma) + 0.1),
      breaks = seq(0.2, 1.6, by = 0.2)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.05, 0.55),
      breaks = seq(0, 0.5, by = 0.1)
    ) +
    ggplot2::labs(
      title = "Bridging simulation input \u03c3\u03b3 to empirical infection-death heritability",
      subtitle = expression(
        "Simulated " * hat(h)[L]^2 *
        " from tetrachoric twin correlations under the full competing-risks DGP"
      ),
      x = expression("Log-hazard frailty dispersion " * sigma[gamma] *
                      " (simulation input)"),
      y = expression("Liability-scale heritability " * hat(h)[L]^2 *
                      " (empirical target)")
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
      label = paste0("t* = ", round(t_star)), hjust = 0, size = 3.8) +
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
      color = mapping_colors["compensatory"], hjust = 0, size = 3.8) +
    ggplot2::annotate("text", x = mean_paper + 1, y = 0.75,
      label = sprintf("mean = %.0f", mean_paper),
      color = mapping_colors["paper"], hjust = 0, size = 3.8) +
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
      vjust = -0.5, size = 3.8) +
    ggplot2::scale_fill_manual(values = setNames(
      mapping_colors, mapping_labels)) +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed",
                        color = "gray50", linewidth = 0.4) +
    ggplot2::annotate("text", x = 3.3, y = 1.03, label = "h\u00B2 = 1",
                      size = 3.5, color = "gray40") +
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

  lbl <- sprintf("alpha = %.2f, beta = %.2f\nR^2 = %.3f\nMean abs. residual = %.1f%%",
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
                      label = lbl, hjust = 0, size = 4.0,
                      family = "mono") +
    ggplot2::labs(
      title = "Pseudo-true variance prediction tracks fitted variance",
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
    ggplot2::scale_color_manual(values = c("True DGP" = unname(PAL["Oracle"]),
                                           "Misspec (r_MZ-only)" = unname(PAL["Biased"]))) +
    ggplot2::annotate("text", x = min(df$r_MZ[1:2]),
                      y = max(df$r_DZ[1:2]) + 0.01,
                      label = lbl, hjust = 0, size = 4.0,
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
                       vjust = -0.5, size = 4.5) +
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
                                    color = model, shape = model,
                                    linetype = model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::scale_color_manual(values = MODEL_PAL) +
    ggplot2::scale_shape_manual(values = MODEL_SHAPES) +
    ggplot2::scale_linetype_manual(values = MODEL_LINETYPES) +
    ggplot2::labs(
      title = expression("Bias vs genetic correlation " * rho * " across all three models"),
      subtitle = "Sign reversal appears in all three models",
      x = expression("Genetic correlation " * rho),
      y = "Bias (pp)", color = "Model", shape = "Model", linetype = "Model"
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
                                    color = model, shape = model,
                                    linetype = model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = MODEL_PAL) +
    ggplot2::scale_shape_manual(values = MODEL_SHAPES) +
    ggplot2::scale_linetype_manual(values = MODEL_LINETYPES) +
    ggplot2::labs(
      title = "Dose-response: bias scales with extrinsic mortality (all models)",
      subtitle = expression(paste(sigma[gamma], " = 0.40, ", rho, " = 0.4")),
      x = expression(m[ex] %*% 10^3),
      y = "Bias (pp)", color = "Model", shape = "Model", linetype = "Model"
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
                       vjust = -0.5, size = 4.5) +
    ggplot2::scale_fill_manual(values = MODEL_PAL, guide = "none") +
    ggplot2::labs(
      title = "Pleiotropy isolation (\u03c1=0) across models",
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
# ===================================================================
# B1: Alternative extrinsic functional forms (bar chart)
# ===================================================================

#' Bar chart of bias across alternative extrinsic frailty forms
#'
#' @param alt_ext_df Data frame from run_alt_ext_forms()
#' @return ggplot object
plot_alt_ext_forms <- function(alt_ext_df) {
  # Clean labels
  form_labels <- c(lognormal = "Log-normal", additive = "Additive",
                   gamma = "Gamma")
  alt_ext_df$form_label <- factor(form_labels[alt_ext_df$ext_form],
    levels = c("Log-normal", "Additive", "Gamma"))

  # Aggregate: mean and SE per form
  agg <- do.call(rbind, lapply(split(alt_ext_df, alt_ext_df$form_label), function(d) {
    data.frame(
      form_label = d$form_label[1],
      mean_pp = mean(d$bias_pp),
      se_pp = sd(d$bias_pp) / sqrt(nrow(d)),
      n_reps = nrow(d),
      stringsAsFactors = FALSE
    )
  }))
  agg$form_label <- factor(agg$form_label,
    levels = c("Log-normal", "Additive", "Gamma"))

  form_colors <- c("Log-normal" = unname(PAL["Biased"]),
                    "Additive" = unname(PAL["ControlA"]),
                    "Gamma" = unname(PAL["Hamilton"]))

  ggplot2::ggplot(agg, ggplot2::aes(x = form_label, y = mean_pp, color = form_label)) +
    # Individual replications as jittered points
    ggplot2::geom_jitter(data = alt_ext_df,
      ggplot2::aes(x = form_label, y = bias_pp, color = form_label),
      width = 0.15, size = 1.5, alpha = 0.35) +
    # Mean + 95% CI
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = mean_pp - 1.96 * se_pp,
                   ymax = mean_pp + 1.96 * se_pp),
      size = 0.8, linewidth = 0.9, fatten = 3) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%+.1f \u00b1 %.1f pp", mean_pp, 1.96 * se_pp)),
      vjust = -1.2, size = 4.5, color = "gray25", show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = form_colors, guide = "none") +
    ggplot2::coord_cartesian(ylim = c(
      min(0, min(alt_ext_df$bias_pp) - 1),
      max(alt_ext_df$bias_pp) + 4)) +
    ggplot2::labs(
      title = "Upward bias is robust to extrinsic frailty construction",
      subtitle = expression(
        paste("Point-range = mean \u00b1 95% CI (10 reps); ",
              sigma[gamma], " = 0.40, ", rho, " = 0.4; GM model")
      ),
      x = "Extrinsic frailty construction",
      y = "Net bias (percentage points)"
    ) +
    theme_pub()
}

# ===================================================================
# B2: Extended sigma_gamma range (line plot)
# ===================================================================

#' Bias vs sigma_gamma beyond the anchored ceiling
#'
#' @param ext_sg_df Data frame from run_extended_sigma_gamma()
#' @return ggplot object
plot_extended_sigma_gamma <- function(ext_sg_df) {
  # Anchored ceiling
  sg_ceil <- 0.65

  # Aggregate: mean and SE per sigma_gamma level
  agg <- do.call(rbind, lapply(split(ext_sg_df, ext_sg_df$sigma_gamma), function(d) {
    data.frame(
      sigma_gamma = d$sigma_gamma[1],
      mean_pp = mean(d$bias_pp),
      se_pp = sd(d$bias_pp) / sqrt(nrow(d)),
      n_reps = nrow(d),
      stringsAsFactors = FALSE
    )
  }))
  agg <- agg[order(agg$sigma_gamma), ]

  ggplot2::ggplot(agg, ggplot2::aes(x = sigma_gamma, y = mean_pp)) +
    # Shade anchored regime
    ggplot2::annotate("rect", xmin = 0.30, xmax = sg_ceil,
                      ymin = -Inf, ymax = Inf,
                      fill = PAL["Oracle"], alpha = 0.08) +
    ggplot2::annotate("text", x = 0.475, y = max(agg$mean_pp) * 0.3,
                      label = "Anchored\nrange",
                      color = PAL["Oracle"], size = 4.0,
                      fontface = "italic", lineheight = 0.9) +
    # 95% CI ribbon
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_pp - 1.96 * se_pp,
                   ymax = mean_pp + 1.96 * se_pp),
      fill = PAL["Biased"], alpha = 0.15) +
    # Individual reps as faint points
    ggplot2::geom_point(data = ext_sg_df,
      ggplot2::aes(x = sigma_gamma, y = bias_pp),
      color = PAL["Biased"], alpha = 0.2, size = 1) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", mean_pp)),
                       vjust = -1, size = 4.0, color = "gray30") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = sg_ceil, linetype = "dotted",
                        color = "gray50", linewidth = 0.5) +
    # Bridge-implied upper bound annotation
    ggplot2::annotate("text", x = max(agg$sigma_gamma), y = max(agg$mean_pp) * 0.55,
                      label = expression("bridge-implied " * sigma[gamma] %~~% 1.47),
                      hjust = 1, size = 3.8, color = "gray40", fontface = "italic") +
    ggplot2::labs(
      title = "Bias continues to increase beyond the anchored ceiling",
      subtitle = expression(
        paste(rho, " = 0.35; ribbon = 95% CI (10 reps); ",
              "dotted = anchored ceiling ", sigma[gamma], " = 0.65")
      ),
      x = expression(sigma[gamma]),
      y = "Bias (percentage points)"
    ) +
    theme_pub()
}

# ===================================================================
# B3: Monte Carlo uncertainty (forest plot)
# ===================================================================

#' Forest plot of per-seed bias estimates with CI band
#'
#' @param mc_result List from run_mc_uncertainty() with $per_seed and $summary
#' @return ggplot object
plot_mc_uncertainty <- function(mc_result) {
  df <- mc_result$per_seed
  summ <- mc_result$summary

  arm2_mean <- summ$mean[summ$statistic == "arm2_bias_pp"]
  arm2_se   <- summ$se[summ$statistic == "arm2_bias_pp"]
  arm1_mean <- summ$mean[summ$statistic == "arm1_bias_pp"]
  arm1_se   <- summ$se[summ$statistic == "arm1_bias_pp"]

  # Reshape to long format
  df_long <- rbind(
    data.frame(seed = df$seed, arm = "Arm 2 (misspecified)",
               bias_pp = df$arm2_bias_pp, stringsAsFactors = FALSE),
    data.frame(seed = df$seed, arm = "Arm 1 (correctly specified)",
               bias_pp = df$arm1_bias_pp, stringsAsFactors = FALSE)
  )
  df_long$arm <- factor(df_long$arm,
    levels = c("Arm 1 (correctly specified)", "Arm 2 (misspecified)"))

  # Summary for mean + CI bands
  ci_df <- data.frame(
    arm = factor(c("Arm 1 (correctly specified)", "Arm 2 (misspecified)"),
                 levels = c("Arm 1 (correctly specified)", "Arm 2 (misspecified)")),
    mean = c(arm1_mean, arm2_mean),
    lo = c(arm1_mean - 1.96 * arm1_se, arm2_mean - 1.96 * arm2_se),
    hi = c(arm1_mean + 1.96 * arm1_se, arm2_mean + 1.96 * arm2_se),
    stringsAsFactors = FALSE
  )

  arm_colors <- c("Arm 1 (correctly specified)" = unname(PAL["Null"]),
                   "Arm 2 (misspecified)" = unname(PAL["Biased"]))

  ggplot2::ggplot(df_long, ggplot2::aes(x = bias_pp, y = factor(seed), color = arm)) +
    ggplot2::facet_wrap(~ arm, scales = "free_x", ncol = 2) +
    # 95% CI band (behind points)
    ggplot2::geom_rect(data = ci_df, inherit.aes = FALSE,
                       ggplot2::aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
                       fill = "steelblue", alpha = 0.10) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    # Mean line (labeled)
    ggplot2::geom_vline(data = ci_df, ggplot2::aes(xintercept = mean),
                        linetype = "solid", linewidth = 0.9, color = "gray20") +
    ggplot2::geom_text(data = ci_df, ggplot2::aes(x = mean, y = Inf, label = "mean"),
                       inherit.aes = FALSE, vjust = 1.5, hjust = -0.1,
                       size = 3.8, color = "gray20", fontface = "italic") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = arm_colors, guide = "none") +
    ggplot2::labs(
      title = "Monte Carlo replication: bias is precisely estimated across 20 seeds",
      subtitle = sprintf(
        "Arm 2 mean = %+.1f pp, 95%% CI [%+.1f, %+.1f]; Arm 1 mean = %+.1f pp, 95%% CI [%+.1f, %+.1f]",
        arm2_mean, arm2_mean - 1.96 * arm2_se, arm2_mean + 1.96 * arm2_se,
        arm1_mean, arm1_mean - 1.96 * arm1_se, arm1_mean + 1.96 * arm1_se
      ),
      x = "Bias (percentage points)",
      y = "Seed"
    ) +
    theme_pub() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
}

# ===================================================================
# B4: High-replication m_ex split (line + ribbon)
# ===================================================================

#' Bias vs heritable fraction with SE ribbons from multiple replications
#'
#' @param hires_df Data frame from run_mex_split_hires()
#' @return ggplot object
plot_mex_split_hires <- function(hires_df) {
  # Aggregate: mean and SE per fraction level
  agg <- do.call(rbind, lapply(split(hires_df, hires_df$frac_heritable), function(d) {
    data.frame(
      frac_heritable = d$frac_heritable[1],
      mean_bias_pp = mean(d$bias_pp),
      se_bias_pp = sd(d$bias_pp) / sqrt(nrow(d)),
      n_reps = nrow(d),
      stringsAsFactors = FALSE
    )
  }))
  agg <- agg[order(agg$frac_heritable), ]

  ggplot2::ggplot(agg, ggplot2::aes(x = frac_heritable, y = mean_bias_pp)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_bias_pp - 1.96 * se_bias_pp,
                   ymax = mean_bias_pp + 1.96 * se_bias_pp),
      fill = PAL["Biased"], alpha = 0.15) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    # Individual replications as faint points
    ggplot2::geom_point(data = hires_df,
                        ggplot2::aes(x = frac_heritable, y = bias_pp),
                        color = PAL["Biased"], alpha = 0.2, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      title = "Bias vs heritable fraction of extrinsic mortality (5-rep average)",
      subtitle = expression(
        paste("Ribbon = 95% CI; faint dots = individual replications; ",
              sigma[gamma], " = 0.40, ", rho, " = 0.4")
      ),
      x = expression(f[heritable] == m[inf] / m[ex]),
      y = "Bias (percentage points)"
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_pub()
}

# ===================================================================
# B3b: Multi-model MC uncertainty (forest plot by model)
# ===================================================================

#' Forest plot of Arm 2 bias across GM, MGG, SR with per-seed dots and CIs
#'
#' @param mm_result List from run_mc_uncertainty_multimodel()
#' @return ggplot object
plot_mc_uncertainty_multimodel <- function(mm_result) {
  df <- mm_result$per_seed
  summ <- mm_result$summary

  # Model factor ordering
  df$model <- factor(df$model, levels = c("GM", "MGG", "SR"))
  summ$model <- factor(summ$model, levels = c("GM", "MGG", "SR"))

  model_colors <- c(GM = unname(MODEL_PAL["GM"]),
                     MGG = unname(MODEL_PAL["MGG"]),
                     SR = unname(MODEL_PAL["SR"]))

  # Build summary data for point-range
  summ$ci_lo_arm2 <- summ$arm2_bias_mean - 1.96 * summ$arm2_bias_se
  summ$ci_hi_arm2 <- summ$arm2_bias_mean + 1.96 * summ$arm2_bias_se
  summ$ci_lo_arm1 <- summ$arm1_bias_mean - 1.96 * summ$arm1_bias_se
  summ$ci_hi_arm1 <- summ$arm1_bias_mean + 1.96 * summ$arm1_bias_se

  # Long format for both arms
  df_long <- rbind(
    data.frame(model = df$model, seed = df$seed,
               arm = "Arm 2 (misspecified)", bias_pp = df$arm2_bias_pp,
               stringsAsFactors = FALSE),
    data.frame(model = df$model, seed = df$seed,
               arm = "Arm 1 (correctly specified)", bias_pp = df$arm1_bias_pp,
               stringsAsFactors = FALSE)
  )
  df_long$arm <- factor(df_long$arm,
    levels = c("Arm 1 (correctly specified)", "Arm 2 (misspecified)"))

  summ_long <- rbind(
    data.frame(model = summ$model, arm = "Arm 2 (misspecified)",
               mean_pp = summ$arm2_bias_mean,
               ci_lo = summ$ci_lo_arm2, ci_hi = summ$ci_hi_arm2,
               stringsAsFactors = FALSE),
    data.frame(model = summ$model, arm = "Arm 1 (correctly specified)",
               mean_pp = summ$arm1_bias_mean,
               ci_lo = summ$ci_lo_arm1, ci_hi = summ$ci_hi_arm1,
               stringsAsFactors = FALSE)
  )
  summ_long$arm <- factor(summ_long$arm,
    levels = c("Arm 1 (correctly specified)", "Arm 2 (misspecified)"))

  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = bias_pp, y = model, color = model)) +
    ggplot2::facet_wrap(~ arm, scales = "free_x", ncol = 2) +
    # Per-seed dots
    ggplot2::geom_jitter(height = 0.15, size = 2.2, alpha = 0.4) +
    # Mean + CI
    ggplot2::geom_pointrange(data = summ_long,
      ggplot2::aes(x = mean_pp, y = model,
                   xmin = ci_lo, xmax = ci_hi),
      size = 0.7, linewidth = 0.9, fatten = 3, color = "gray15") +
    # Mean value labels
    ggplot2::geom_text(data = summ_long,
      ggplot2::aes(x = mean_pp, y = model,
                   label = sprintf("%+.1f", mean_pp)),
      vjust = -1.3, size = 3.8, color = "gray25", fontface = "bold") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = model_colors, guide = "none") +
    ggplot2::labs(
      title = "Misspecification-induced upward bias replicates across all three mortality models",
      subtitle = "20 independent seeds per model; black point-range = mean \u00b1 95% CI",
      x = "Bias (percentage points)",
      y = "Mortality model"
    ) +
    theme_pub()
}

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
                                 alt_ext_forms = NULL,
                                 extended_sigma_gamma = NULL,
                                 mc_uncertainty = NULL,
                                 mc_uncertainty_multimodel = NULL,
                                 mex_split_hires = NULL,
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

  # B1-B4 revision figures
  if (!is.null(alt_ext_forms)) {
    paths$alt_ext_forms <- save_figure(
      plot_alt_ext_forms(alt_ext_forms),
      "fig_b1_alt_ext_forms.png", width = 8, height = 5.5, outdir = outdir)
  }

  if (!is.null(extended_sigma_gamma)) {
    paths$extended_sigma_gamma <- save_figure(
      plot_extended_sigma_gamma(extended_sigma_gamma),
      "fig_b2_extended_sigma_gamma.png", width = 8, height = 5.5, outdir = outdir)
  }

  if (!is.null(mc_uncertainty)) {
    paths$mc_uncertainty <- save_figure(
      plot_mc_uncertainty(mc_uncertainty),
      "fig_b3_mc_uncertainty.png", width = 12, height = 7, outdir = outdir)
  }

  if (!is.null(mc_uncertainty_multimodel)) {
    paths$mc_uncertainty_multimodel <- save_figure(
      plot_mc_uncertainty_multimodel(mc_uncertainty_multimodel),
      "fig_b3b_mc_multimodel.png", width = 10, height = 5.5, outdir = outdir)
  }

  if (!is.null(mex_split_hires)) {
    paths$mex_split_hires <- save_figure(
      plot_mex_split_hires(mex_split_hires),
      "fig_b4_mex_split_hires.png", width = 8, height = 5.5, outdir = outdir)
  }

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
