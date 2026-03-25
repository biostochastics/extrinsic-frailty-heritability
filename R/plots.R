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

# Standard annotation text sizes (absolute, not relative to theme base_size).
# All figures use these to ensure consistent rendered font size.
LABEL_SIZE <- 3.8    # geom_text: data-value labels on points/bars
ANNOT_SIZE <- 3.5    # annotate("text"): region labels, reference annotations
ANNOT_EMPH <- 4.0    # annotate("text"): emphasised callouts (bold)

# Patchwork super-title theme (reused across all composite figures)
patchwork_title_theme <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold", size = 16, family = "Helvetica Neue",
      margin = ggplot2::margin(b = 6)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 11, family = "Helvetica Neue", color = "gray35"
    )
  )
}

# ===================================================================
# Figure 1: Two-panel main arms (pointrange + lollipop)
# ===================================================================

#' Two-panel Figure 1(B): pointrange for h², lollipop for bias
#'
#' Left panel: h² pointrange with 95% CI (cloud of per-seed values optional).
#' Right panel: bias relative to oracle as lollipop stems from zero.
#'
#' @param main_arms_replicated List from run_main_arms_replicated()
#' @param robustness_cutoff0 List from run_cutoff0_robustness() (replicated)
#' @return patchwork composite
plot_main_arms <- function(main_arms_replicated, robustness_cutoff0 = NULL) {
  s <- main_arms_replicated$summary
  get_row <- function(metric) {
    r <- s[s$metric == metric, ]
    if (nrow(r) == 0) {
      warning(sprintf("Metric '%s' not found in replicated summary", metric))
      return(list(mean = NA_real_, se = NA_real_, lo95 = NA_real_, hi95 = NA_real_))
    }
    list(mean = r$mean[1], se = r$se[1], lo95 = r$lo95[1], hi95 = r$hi95[1])
  }

  oracle <- get_row("oracle_h2")
  arms <- list(
    list(name = "Oracle", h2 = oracle, bias = list(mean = 0, lo95 = NA, hi95 = NA)),
    list(name = "Baseline: correctly specified",
         h2 = get_row("arm1_h2"), bias = get_row("arm1_bias_pp")),
    list(name = "Misspecified: omitted familial extrinsic",
         h2 = get_row("arm2_h2"), bias = get_row("arm2_bias_pp")),
    list(name = "Control: nonfamilial extrinsic",
         h2 = get_row("ctrl_a_h2"), bias = get_row("ctrl_a_bias_pp")),
    list(name = "Control: irrelevant trait",
         h2 = get_row("irrel_h2"), bias = get_row("irrel_bias_pp")),
    list(name = "Control: vanishing extrinsic",
         h2 = get_row("ctrl_b_h2"), bias = get_row("ctrl_b_bias_pp")),
    list(name = "Check: pleiotropy only",
         h2 = get_row("pleio_iso_h2"), bias = get_row("pleio_iso_bias_pp")),
    list(name = "Recovery: two-component refit",
         h2 = get_row("fix_h2"), bias = get_row("fix_bias_pp"))
  )

  if (!is.null(robustness_cutoff0)) {
    arms <- c(arms, list(list(
      name = "Check: no survivor cutoff",
      h2 = list(mean = robustness_cutoff0$corr_h2,
                lo95 = robustness_cutoff0$corr_h2_lo95,
                hi95 = robustness_cutoff0$corr_h2_hi95),
      bias = list(mean = robustness_cutoff0$bias_pp, lo95 = NA, hi95 = NA)
    )))
  }

  arms <- c(arms, list(list(
    name = "Comparator: Hamilton",
    h2 = get_row("arm3_h2"), bias = get_row("arm3_bias_pp")
  )))

  df <- data.frame(
    Condition = sapply(arms, `[[`, "name"),
    h2 = sapply(arms, function(a) a$h2$mean),
    h2_lo = sapply(arms, function(a) if (!is.null(a$h2$lo95)) a$h2$lo95 else NA),
    h2_hi = sapply(arms, function(a) if (!is.null(a$h2$hi95)) a$h2$hi95 else NA),
    bias_pp = sapply(arms, function(a) a$bias$mean),
    bias_lo = sapply(arms, function(a) if (!is.null(a$bias$lo95)) a$bias$lo95 else NA),
    bias_hi = sapply(arms, function(a) if (!is.null(a$bias$hi95)) a$bias$hi95 else NA),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$h2), ]

  # Display order: Oracle at top after coord_flip (= last factor level)
  display_order <- c("Comparator: Hamilton", "Control: vanishing extrinsic",
                     "Control: irrelevant trait", "Control: nonfamilial extrinsic",
                     "Check: pleiotropy only", "Check: no survivor cutoff",
                     "Recovery: two-component refit",
                     "Misspecified: omitted familial extrinsic", "Baseline: correctly specified", "Oracle")
  display_order <- display_order[display_order %in% df$Condition]
  df$Condition <- factor(df$Condition, levels = display_order)

  cond_names <- c("Oracle", "Baseline: correctly specified", "Misspecified: omitted familial extrinsic",
                   "Control: nonfamilial extrinsic", "Control: irrelevant trait",
                   "Control: vanishing extrinsic", "Check: pleiotropy only",
                   "Recovery: two-component refit", "Check: no survivor cutoff", "Comparator: Hamilton")
  pal_keys <- c("Oracle", "Null", "Biased", "ControlA", "Irrelevant",
                "ControlB", "Pleiotropy", "Fix", "Cutoff0", "Hamilton")
  condition_colors <- setNames(unname(PAL[pal_keys]), cond_names)

  oracle_h2 <- df$h2[df$Condition == "Oracle"]

  # --- Build per-seed long data for jittered dots ---
  ps <- main_arms_replicated$per_seed
  h2_cols <- c(oracle_h2 = "Oracle", arm1_h2 = "Baseline: correctly specified",
               arm2_h2 = "Misspecified: omitted familial extrinsic",
               ctrl_a_h2 = "Control: nonfamilial extrinsic",
               irrel_h2 = "Control: irrelevant trait",
               ctrl_b_h2 = "Control: vanishing extrinsic",
               pleio_iso_h2 = "Check: pleiotropy only",
               fix_h2 = "Recovery: two-component refit", arm3_h2 = "Comparator: Hamilton")
  seed_rows <- lapply(names(h2_cols), function(col) {
    if (col %in% names(ps))
      data.frame(Condition = h2_cols[[col]], h2_seed = ps[[col]],
                 stringsAsFactors = FALSE)
  })
  # Add Cutoff = 0 per-seed data if available
  if (!is.null(robustness_cutoff0) && !is.null(robustness_cutoff0$per_seed)) {
    seed_rows <- c(seed_rows, list(
      data.frame(Condition = "Check: no survivor cutoff",
                 h2_seed = robustness_cutoff0$per_seed$corr_h2,
                 stringsAsFactors = FALSE)
    ))
  }
  seed_long <- do.call(rbind, seed_rows)
  seed_long$Condition <- factor(seed_long$Condition, levels = display_order)

  # --- Left panel: h² pointrange with per-seed jitter ---
  # Find position of Misspecified condition (highlight row)
  misspec_cond <- "Misspecified: omitted familial extrinsic"
  misspec_pos <- which(display_order == misspec_cond)

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = h2)) +
    # Gray background highlight for Misspecified row
    ggplot2::annotate("rect",
      xmin = misspec_pos - 0.5, xmax = misspec_pos + 0.5,
      ymin = 0.37, ymax = 0.66,
      fill = "gray88", color = NA) +
    ggplot2::geom_hline(yintercept = oracle_h2,
                        linetype = "dashed", color = "gray55", linewidth = 0.4) +
    ggplot2::geom_boxplot(data = seed_long,
      ggplot2::aes(x = Condition, y = h2_seed, fill = Condition),
      width = 0.4, alpha = 0.35, color = "gray50",
      outlier.shape = NA, inherit.aes = FALSE) +
    ggplot2::geom_jitter(data = seed_long,
      ggplot2::aes(x = Condition, y = h2_seed, color = Condition),
      width = 0.15, size = 1.2, alpha = 0.3, inherit.aes = FALSE) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = h2_lo, ymax = h2_hi),
      linewidth = 0.7, color = "gray20") +
    ggplot2::geom_point(size = 1.8, color = "gray10") +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(limits = display_order) +
    ggplot2::scale_color_manual(values = condition_colors, guide = "none") +
    ggplot2::scale_fill_manual(values = condition_colors, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0.38, 0.65),
      breaks = seq(0.40, 0.65, by = 0.05)) +
    ggplot2::labs(title = bquote(hat(h)^2 ~ ": mean \u00b1 95% CI across seeds"),
                  subtitle = NULL,
                  x = NULL, y = expression(hat(h)^2)) +
    theme_pub()

  # --- Right panel: bias lollipop ---
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = bias_pp,
                                          color = Condition)) +
    # Gray background highlight for Misspecified row
    ggplot2::annotate("rect",
      xmin = misspec_pos - 0.5, xmax = misspec_pos + 0.5,
      ymin = -10, ymax = 10,
      fill = "gray88", color = NA) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray55", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(xend = Condition, y = 0, yend = bias_pp),
                          linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(limits = display_order) +
    ggplot2::scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 5)) +
    ggplot2::scale_color_manual(values = condition_colors, guide = "none") +
    ggplot2::labs(title = "Bias vs oracle (pp)",
                  x = NULL, y = "Bias (pp)") +
    theme_pub() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())

  p1 + p2 + patchwork::plot_layout(widths = c(1.3, 1))
}

#' Legacy bar chart of h² (single-seed, no CIs) — fallback only
#'
#' @param summary_table Data frame from build_summary_table()
#' @return ggplot object
plot_main_arms_legacy <- function(summary_table) {
  df <- summary_table
  df <- df[!is.na(df$h2), ]
  df$Condition <- factor(df$Condition, levels = rev(df$Condition))
  df$bias_pp <- 100 * df$bias

  condition_colors <- setNames(
    unname(PAL[c("Oracle", "Null", "Biased", "ControlA", "Irrelevant",
                 "ControlB", "Pleiotropy", "Fix", "Cutoff0", "Hamilton")]),
    c("Oracle", "Baseline: correctly specified",
      "Misspecified: omitted familial extrinsic", "Control: nonfamilial extrinsic",
      "Control: irrelevant trait", "Control: vanishing extrinsic",
      "Check: pleiotropy only", "Recovery: two-component refit",
      "Check: no survivor cutoff", "Comparator: Hamilton")
  )

  oracle_h2 <- df$h2[df$Condition == "Oracle"]

  ggplot2::ggplot(df, ggplot2::aes(x = Condition, y = h2, fill = Condition)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = oracle_h2,
                        linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%.3f\n(%+.1f pp)", h2, bias_pp)),
      hjust = -0.05, size = LABEL_SIZE, color = "gray30") +
    ggplot2::coord_flip(ylim = c(0, max(df$h2, na.rm = TRUE) * 1.25)) +
    ggplot2::scale_fill_manual(values = condition_colors, guide = "none") +
    ggplot2::labs(
      title = expression(hat(h)^2 ~ "across experimental arms"),
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
      title = "Bias across parameter space",
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
      title = "Anchored region",
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
    Arm = rep(c("Baseline: correctly specified", "Misspecified: omitted familial extrinsic"), each = nrow(model_table)),
    bias_pp = c(model_table$Baseline_bias_pp, model_table$Misspec_bias_pp),
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
      vjust = -0.5, size = LABEL_SIZE
    ) +
    ggplot2::scale_fill_manual(
      values = c("Baseline: correctly specified" = unname(PAL["Null"]),
                 "Misspecified: omitted familial extrinsic" = unname(PAL["Biased"]))
    ) +
    ggplot2::labs(
      title = "Bias across mortality models",
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
                       hjust = -0.1,
                       size = LABEL_SIZE) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.15))) +
    ggplot2::scale_fill_manual(values = c("TRUE" = unname(PAL["Biased"]),
                                          "FALSE" = unname(PAL["Oracle"])),
                               guide = "none") +
    ggplot2::labs(
      title = "Control experiments",
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
plot_dose_response <- function(dose_response_obj) {
  # Accept either a list (new format: $summary + $per_seed) or a data.frame (legacy)
  if (is.data.frame(dose_response_obj)) {
    dose_response_df <- dose_response_obj
    per_seed <- NULL
  } else {
    dose_response_df <- dose_response_obj$summary
    per_seed <- dose_response_obj$per_seed
  }

  p <- ggplot2::ggplot(dose_response_df,
                  ggplot2::aes(x = m_ex * 1000, y = bias_pp)) +
    {if ("se_pp" %in% names(dose_response_df))
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lo95_pp, ymax = hi95_pp),
        fill = PAL["Biased"], alpha = 0.15)
    }

  # Individual seed dots (faint, behind summary)
  if (!is.null(per_seed)) {
    p <- p + ggplot2::geom_point(data = per_seed,
      ggplot2::aes(x = m_ex * 1000, y = bias_pp),
      color = PAL["Biased"], alpha = 0.2, size = 1,
      inherit.aes = FALSE)
  }

  p +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      title = "Dose-response",
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
                      color = PAL["Pleiotropy"], size = ANNOT_SIZE,
                      fontface = "italic", lineheight = 0.9) +
    # Shaded region: anchored sensitivity range
    ggplot2::annotate("rect",
                      xmin = rho_lo, xmax = rho_hi,
                      ymin = -Inf, ymax = Inf,
                      fill = PAL["Oracle"], alpha = 0.08) +
    ggplot2::annotate("text", x = (rho_lo + rho_hi) / 2, y = max(neg_rho_df$bias_pp) * 0.9,
                      label = "Anchored\nregime",
                      color = PAL["Oracle"], size = ANNOT_SIZE,
                      fontface = "italic") +
    # Uncertainty ribbon (if available)
    {if ("se_pp" %in% names(neg_rho_df))
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lo95_pp, ymax = hi95_pp),
        fill = PAL["Biased"], alpha = 0.15)
    } +
    # Data
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50",
                        linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray50",
                        linewidth = 0.6) +
    ggplot2::labs(
      title = expression("Bias vs " * rho),
      subtitle = expression(sigma[gamma] == 0.40),
      x = expression("Genetic correlation " * rho),
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
      title = "Heritable fraction sensitivity",
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
                      label = "Sensitivity\nrange",
                      color = PAL["Oracle"], size = ANNOT_SIZE,
                      fontface = "italic", lineheight = 0.9) +
    # Raw simulated curve (no smoothing)
    ggplot2::geom_line(color = PAL["Oracle"], linewidth = 0.9) +
    # Obel target line
    ggplot2::geom_hline(yintercept = 0.40, linetype = "dashed",
                        color = PAL["Biased"], linewidth = 0.7) +
    ggplot2::annotate("text", x = 0.15, y = 0.425,
                      label = expression(Obel~et~al.~target:~h[L]^2 == 0.40),
                      color = PAL["Biased"], size = ANNOT_SIZE, hjust = 0)

  if (!is.na(bridge_sg)) {
    p <- p +
      ggplot2::geom_vline(xintercept = bridge_sg, linetype = "dotted",
                          color = PAL["Null"], linewidth = 0.7) +
      ggplot2::annotate("text", x = bridge_sg - 0.05, y = 0.12,
                        label = sprintf("\u03c3\u03b3 \u2248 %.2f", bridge_sg),
                        color = PAL["Null"], size = ANNOT_EMPH, hjust = 1,
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
      title = expression(sigma[gamma] ~ "bridge curve"),
      x = expression(sigma[gamma]),
      y = expression(hat(h)[L]^2)
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
      title = "MGG intrinsic hazard",
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
      label = paste0("t* = ", round(t_star)), hjust = 0, size = ANNOT_SIZE) +
    ggplot2::labs(
      title = "A. Intrinsic hazard",
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
      color = mapping_colors["compensatory"], hjust = 0, size = ANNOT_SIZE) +
    ggplot2::annotate("text", x = mean_paper + 1, y = 0.75,
      label = sprintf("mean = %.0f", mean_paper),
      color = mapping_colors["paper"], hjust = 0, size = ANNOT_SIZE) +
    ggplot2::labs(
      title = "B. Cohort survival",
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
      vjust = -0.5, size = LABEL_SIZE) +
    ggplot2::scale_fill_manual(values = setNames(
      mapping_colors, mapping_labels)) +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed",
                        color = "gray50", linewidth = 0.4) +
    ggplot2::annotate("text", x = 3.3, y = 1.03, label = "h\u00B2 = 1",
                      size = ANNOT_SIZE, color = "gray40") +
    ggplot2::labs(
      title = expression("C. Twin correlations & " * hat(h)^2),
      x = NULL, y = "Value", fill = "SM mapping") +
    ggplot2::ylim(NA, max(tdf_long$value) * 1.15) +
    theme_pub() +
    ggplot2::theme(legend.position = "bottom")

  # Combine with patchwork
  patchwork::wrap_plots(p_hazard, p_surv, p_twin, ncol = 3) +
    patchwork::plot_annotation(
      title = "MGG SM parameterization comparison",
      theme = patchwork_title_theme()
    )
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
                      label = lbl, hjust = 0, size = ANNOT_SIZE,
                      family = "mono") +
    ggplot2::labs(
      title = "Variance decomposition",
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
                      label = lbl, hjust = 0, size = ANNOT_SIZE,
                      family = "mono") +
    ggplot2::labs(
      title = "Joint calibration",
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
                       vjust = -0.5, size = LABEL_SIZE) +
    ggplot2::labs(
      title = "SR time-step sensitivity",
      x = expression(Delta * t ~ "(fraction of year)"),
      y = "Misspecified bias (pp)"
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
    bias_pp = neg_rho_gm$bias_pp,
    se_pp = if ("se_pp" %in% names(neg_rho_gm)) neg_rho_gm$se_pp else NA_real_,
    lo95_pp = if ("lo95_pp" %in% names(neg_rho_gm)) neg_rho_gm$lo95_pp else NA_real_,
    hi95_pp = if ("hi95_pp" %in% names(neg_rho_gm)) neg_rho_gm$hi95_pp else NA_real_,
    stringsAsFactors = FALSE
  )
  mc_rho <- model_controls[model_controls$sweep_type == "negative_rho", ]
  mc_df <- data.frame(
    model = toupper(mc_rho$model),
    rho = mc_rho$sweep_val,
    bias_pp = mc_rho$bias_pp,
    se_pp = if ("se_pp" %in% names(mc_rho)) mc_rho$se_pp else NA_real_,
    lo95_pp = if ("lo95_pp" %in% names(mc_rho)) mc_rho$lo95_pp else NA_real_,
    hi95_pp = if ("hi95_pp" %in% names(mc_rho)) mc_rho$hi95_pp else NA_real_,
    stringsAsFactors = FALSE
  )
  df <- rbind(gm_df, mc_df)
  df$model <- factor(df$model, levels = c("GM", "MGG", "SR"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = rho, y = bias_pp,
                                    color = model, shape = model,
                                    linetype = model, fill = model))
  # Add CI ribbons if available
  if (any(!is.na(df$lo95_pp))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo95_pp, ymax = hi95_pp),
      alpha = 0.15, linewidth = 0, color = NA)
  }
  p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::scale_color_manual(values = MODEL_PAL) +
    ggplot2::scale_fill_manual(values = MODEL_PAL, guide = "none") +
    ggplot2::scale_shape_manual(values = MODEL_SHAPES) +
    ggplot2::scale_linetype_manual(values = MODEL_LINETYPES) +
    ggplot2::labs(
      title = expression("Bias vs " * rho ~ "(all models)"),
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
  # Accept either list (new) or data.frame (legacy)
  gm_summ <- if (is.data.frame(dose_response_gm)) dose_response_gm else dose_response_gm$summary
  gm_df <- data.frame(
    model = "GM", m_ex = gm_summ$m_ex,
    bias_pp = gm_summ$bias_pp,
    se_pp = if ("se_pp" %in% names(gm_summ)) gm_summ$se_pp else NA_real_,
    lo95_pp = if ("lo95_pp" %in% names(gm_summ)) gm_summ$lo95_pp else NA_real_,
    hi95_pp = if ("hi95_pp" %in% names(gm_summ)) gm_summ$hi95_pp else NA_real_,
    stringsAsFactors = FALSE
  )
  mc_dr <- model_controls[model_controls$sweep_type == "dose_response", ]
  mc_df <- data.frame(
    model = toupper(mc_dr$model),
    m_ex = mc_dr$sweep_val,
    bias_pp = mc_dr$bias_pp,
    se_pp = if ("se_pp" %in% names(mc_dr)) mc_dr$se_pp else NA_real_,
    lo95_pp = if ("lo95_pp" %in% names(mc_dr)) mc_dr$lo95_pp else NA_real_,
    hi95_pp = if ("hi95_pp" %in% names(mc_dr)) mc_dr$hi95_pp else NA_real_,
    stringsAsFactors = FALSE
  )
  df <- rbind(gm_df, mc_df)
  df$model <- factor(df$model, levels = c("GM", "MGG", "SR"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = m_ex * 1000, y = bias_pp,
                                    color = model, shape = model,
                                    linetype = model, fill = model))
  if (any(!is.na(df$lo95_pp))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo95_pp, ymax = hi95_pp),
      alpha = 0.15, linewidth = 0, color = NA)
  }
  p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = MODEL_PAL) +
    ggplot2::scale_fill_manual(values = MODEL_PAL, guide = "none") +
    ggplot2::scale_shape_manual(values = MODEL_SHAPES) +
    ggplot2::scale_linetype_manual(values = MODEL_LINETYPES) +
    ggplot2::labs(
      title = "Dose-response (all models)",
      subtitle = expression(paste(sigma[gamma], " = 0.40, ", rho, " = 0.4")),
      x = expression(m[ex] %*% 10^3),
      y = "Bias (pp)", color = "Model", shape = "Model", linetype = "Model"
    ) +
    theme_pub() +
    ggplot2::theme(legend.position = c(0.15, 0.85))
}

# ===================================================================
# Multi-model comparison: Check pleiotropy only
# ===================================================================

#' Multi-model pleiotropy isolation bar chart
#'
#' @param controls_table GM controls table
#' @param model_controls Combined model controls data frame
#' @return ggplot object
plot_pleiotropy_multimodel <- function(controls_table, model_controls) {
  mc_pl <- model_controls[model_controls$sweep_type == "pleiotropy_isolation", ]

  # GM arm2 reference and rho=0 from controls_table
  gm_arm2 <- controls_table$bias_pp[controls_table$Control == "Misspecified: omitted familial extrinsic"]
  gm_rho0 <- controls_table$bias_pp[controls_table$Control == "Check: pleiotropy only"]

  df <- data.frame(
    Model = c("GM", "GM", "MGG", "MGG", "SR", "SR"),
    Condition = rep(c("Misspecified: omitted familial extrinsic", "Check: pleiotropy only"), 3),
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
                       vjust = -0.5, size = LABEL_SIZE) +
    ggplot2::scale_fill_manual(values = MODEL_PAL, guide = "none") +
    ggplot2::labs(
      title = expression("Check: pleiotropy only (" * rho * " = 0)"),
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
      theme = patchwork_title_theme()
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
  agg$t_crit <- qt(0.975, df = agg$n_reps - 1)
  agg$lo95 <- agg$mean_pp - agg$t_crit * agg$se_pp
  agg$hi95 <- agg$mean_pp + agg$t_crit * agg$se_pp

  form_colors <- c("Log-normal" = unname(PAL["Biased"]),
                    "Additive" = unname(PAL["ControlA"]),
                    "Gamma" = unname(PAL["Hamilton"]))

  ggplot2::ggplot(agg, ggplot2::aes(x = form_label, y = mean_pp, color = form_label)) +
    # Individual replications as jittered points
    ggplot2::geom_jitter(data = alt_ext_df,
      ggplot2::aes(x = form_label, y = bias_pp, color = form_label),
      width = 0.15, size = 1.5, alpha = 0.35) +
    # Mean + 95% CI (t-distribution)
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = lo95, ymax = hi95),
      size = 0.8, linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%+.1f \u00b1 %.1f pp", mean_pp, t_crit * se_pp)),
      vjust = -1.2, size = LABEL_SIZE, color = "gray25", show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = form_colors, guide = "none") +
    ggplot2::coord_cartesian(ylim = c(
      min(0, min(alt_ext_df$bias_pp) - 1),
      max(alt_ext_df$bias_pp) + 4)) +
    ggplot2::labs(
      title = "Functional form robustness",
      subtitle = expression(paste(sigma[gamma], " = 0.40, ", rho, " = 0.4")),
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
  agg$t_crit <- qt(0.975, df = agg$n_reps - 1)
  agg$lo95 <- agg$mean_pp - agg$t_crit * agg$se_pp
  agg$hi95 <- agg$mean_pp + agg$t_crit * agg$se_pp

  ggplot2::ggplot(agg, ggplot2::aes(x = sigma_gamma, y = mean_pp)) +
    # Shade anchored regime
    ggplot2::annotate("rect", xmin = 0.30, xmax = sg_ceil,
                      ymin = -Inf, ymax = Inf,
                      fill = PAL["Oracle"], alpha = 0.08) +
    ggplot2::annotate("text", x = 0.475, y = max(agg$mean_pp) * 0.3,
                      label = "Anchored\nrange",
                      color = PAL["Oracle"], size = ANNOT_SIZE,
                      fontface = "italic", lineheight = 0.9) +
    # 95% CI ribbon (t-distribution)
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo95, ymax = hi95),
      fill = PAL["Biased"], alpha = 0.15) +
    # Individual reps as faint points
    ggplot2::geom_point(data = ext_sg_df,
      ggplot2::aes(x = sigma_gamma, y = bias_pp),
      color = PAL["Biased"], alpha = 0.2, size = 1) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", mean_pp)),
                       vjust = -1, size = LABEL_SIZE, color = "gray30") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::geom_vline(xintercept = sg_ceil, linetype = "dotted",
                        color = "gray50", linewidth = 0.5) +
    # Bridge-implied upper bound annotation
    ggplot2::annotate("text", x = max(agg$sigma_gamma), y = max(agg$mean_pp) * 0.55,
                      label = expression("bridge-implied " * sigma[gamma] %~~% 1.47),
                      hjust = 1, size = ANNOT_SIZE, color = "gray40", fontface = "italic") +
    ggplot2::labs(
      title = expression("Extended " * sigma[gamma] * " stress test"),
      subtitle = expression(paste(rho, " = 0.35")),
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

  n_seeds <- nrow(df)
  t_crit <- qt(0.975, df = n_seeds - 1)
  arm2_mean <- summ$mean[summ$statistic == "arm2_bias_pp"]
  arm2_se   <- summ$se[summ$statistic == "arm2_bias_pp"]
  arm1_mean <- summ$mean[summ$statistic == "arm1_bias_pp"]
  arm1_se   <- summ$se[summ$statistic == "arm1_bias_pp"]

  # Reshape to long format
  df_long <- rbind(
    data.frame(seed = df$seed, arm = "Misspecified: omitted familial extrinsic",
               bias_pp = df$arm2_bias_pp, stringsAsFactors = FALSE),
    data.frame(seed = df$seed, arm = "Baseline: correctly specified",
               bias_pp = df$arm1_bias_pp, stringsAsFactors = FALSE)
  )
  df_long$arm <- factor(df_long$arm,
    levels = c("Baseline: correctly specified", "Misspecified: omitted familial extrinsic"))

  # Summary for mean + CI bands
  ci_df <- data.frame(
    arm = factor(c("Baseline: correctly specified", "Misspecified: omitted familial extrinsic"),
                 levels = c("Baseline: correctly specified", "Misspecified: omitted familial extrinsic")),
    mean = c(arm1_mean, arm2_mean),
    lo = c(arm1_mean - t_crit * arm1_se, arm2_mean - t_crit * arm2_se),
    hi = c(arm1_mean + t_crit * arm1_se, arm2_mean + t_crit * arm2_se),
    stringsAsFactors = FALSE
  )

  arm_colors <- c("Baseline: correctly specified" = unname(PAL["Null"]),
                   "Misspecified: omitted familial extrinsic" = unname(PAL["Biased"]))

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
                       size = LABEL_SIZE, color = "gray20", fontface = "italic") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = arm_colors, guide = "none") +
    ggplot2::labs(
      title = "MC uncertainty",
      x = "Bias (percentage points)",
      y = "Seed"
    ) +
    theme_pub() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.72)))
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

  agg$t_crit <- qt(0.975, df = agg$n_reps - 1)
  agg$lo95 <- agg$mean_bias_pp - agg$t_crit * agg$se_bias_pp
  agg$hi95 <- agg$mean_bias_pp + agg$t_crit * agg$se_bias_pp

  ggplot2::ggplot(agg, ggplot2::aes(x = frac_heritable, y = mean_bias_pp)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo95, ymax = hi95),
      fill = PAL["Biased"], alpha = 0.15) +
    ggplot2::geom_line(color = PAL["Biased"], linewidth = 1) +
    ggplot2::geom_point(color = PAL["Biased"], size = 2.5) +
    # Individual replications as faint points
    ggplot2::geom_point(data = hires_df,
                        ggplot2::aes(x = frac_heritable, y = bias_pp),
                        color = PAL["Biased"], alpha = 0.2, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::labs(
      title = "Heritable fraction sensitivity",
      subtitle = expression(paste(sigma[gamma], " = 0.40, ", rho, " = 0.4")),
      x = expression(f[heritable] == m[inf] / m[ex]),
      y = "Bias (percentage points)"
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_pub()
}

# ===================================================================
# B3b: Multi-model MC uncertainty (forest plot by model)
# ===================================================================

#' Forest plot of Misspecified bias across GM, MGG, SR with per-seed dots and CIs
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

  # Build summary data for point-range (t-distribution CIs)
  t_crit <- qt(0.975, df = summ$n_seeds - 1)
  summ$ci_lo_arm2 <- summ$arm2_bias_mean - t_crit * summ$arm2_bias_se
  summ$ci_hi_arm2 <- summ$arm2_bias_mean + t_crit * summ$arm2_bias_se
  summ$ci_lo_arm1 <- summ$arm1_bias_mean - t_crit * summ$arm1_bias_se
  summ$ci_hi_arm1 <- summ$arm1_bias_mean + t_crit * summ$arm1_bias_se

  # Long format for both arms
  df_long <- rbind(
    data.frame(model = df$model, seed = df$seed,
               arm = "Misspecified: omitted familial extrinsic", bias_pp = df$arm2_bias_pp,
               stringsAsFactors = FALSE),
    data.frame(model = df$model, seed = df$seed,
               arm = "Baseline: correctly specified", bias_pp = df$arm1_bias_pp,
               stringsAsFactors = FALSE)
  )
  df_long$arm <- factor(df_long$arm,
    levels = c("Baseline: correctly specified", "Misspecified: omitted familial extrinsic"))

  summ_long <- rbind(
    data.frame(model = summ$model, arm = "Misspecified: omitted familial extrinsic",
               mean_pp = summ$arm2_bias_mean,
               ci_lo = summ$ci_lo_arm2, ci_hi = summ$ci_hi_arm2,
               stringsAsFactors = FALSE),
    data.frame(model = summ$model, arm = "Baseline: correctly specified",
               mean_pp = summ$arm1_bias_mean,
               ci_lo = summ$ci_lo_arm1, ci_hi = summ$ci_hi_arm1,
               stringsAsFactors = FALSE)
  )
  summ_long$arm <- factor(summ_long$arm,
    levels = c("Baseline: correctly specified", "Misspecified: omitted familial extrinsic"))

  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = bias_pp, y = model, color = model)) +
    ggplot2::facet_wrap(~ arm, scales = "free_x", ncol = 2) +
    # Per-seed dots
    ggplot2::geom_jitter(height = 0.15, size = 2.2, alpha = 0.4) +
    # Mean + CI
    ggplot2::geom_pointrange(data = summ_long,
      ggplot2::aes(x = mean_pp, y = model,
                   xmin = ci_lo, xmax = ci_hi),
      size = 0.8, linewidth = 0.7, color = "gray15") +
    # Mean value labels
    ggplot2::geom_text(data = summ_long,
      ggplot2::aes(x = mean_pp, y = model,
                   label = sprintf("%+.1f", mean_pp)),
      vjust = -1.3, size = LABEL_SIZE, color = "gray25", fontface = "bold") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    ggplot2::scale_color_manual(values = model_colors, guide = "none") +
    ggplot2::labs(
      title = "Multi-model MC uncertainty",
      x = "Bias (percentage points)",
      y = "Mortality model"
    ) +
    theme_pub()
}

# ===================================================================
# Bivariate survival model check (main-text + appendix figures)
# ===================================================================

#' Age-specific dependence ratio R(t) for MZ and DZ twins
#'
#' Main-text figure: three lines (True DGP, Misspecified fit, Recovery: two-component refit)
#' showing R(t) = log(P(T1>t, T2>t) / P(T>t)^2) across age t.
#' The misspecified model gets the *timing* of concordance wrong even
#' after matching overall r_MZ; the oracle fix collapses onto the true DGP.
#'
#' @param biv_check List from run_bivariate_survival_check()
#' @return ggplot object (patchwork composite: MZ / DZ panels)
plot_bivariate_dep_curves <- function(biv_check) {
  dep <- biv_check$dep_curves
  dep <- dep[!is.na(dep$R), ]

  dep$model <- factor(dep$model,
    levels = c("True DGP", "Misspecified fit", "Recovery: two-component refit"))

  model_colors <- c("True DGP"         = "#2166ac",
                     "Misspecified fit"  = "#b2182b",
                     "Recovery: two-component refit" = "#ff7f00")
  model_lty    <- c("True DGP"         = "solid",
                     "Misspecified fit"  = "dashed",
                     "Recovery: two-component refit" = "dotdash")
  model_lw     <- c("True DGP"         = 1.4,
                     "Misspecified fit"  = 1.4,
                     "Recovery: two-component refit" = 1.2)

  make_panel <- function(zyg, show_legend = FALSE) {
    d <- dep[dep$zygosity == zyg, ]

    # --- Zoom inset (ages 75–85) ---
    d_zoom <- d[d$age >= 75 & d$age <= 85, ]
    p_zoom <- ggplot2::ggplot(d_zoom,
                ggplot2::aes(x = age, y = R, color = model, linetype = model,
                             linewidth = model)) +
      ggplot2::geom_line() +
      ggplot2::scale_linewidth_manual(values = model_lw, guide = "none") +
      ggplot2::scale_color_manual(values = model_colors, guide = "none") +
      ggplot2::scale_linetype_manual(values = model_lty, guide = "none") +
      ggplot2::labs(title = "Ages 75\u201385", x = NULL, y = NULL) +
      theme_pub() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text = ggplot2::element_text(size = 7),
        axis.ticks = ggplot2::element_line(linewidth = 0.3),
        plot.background = ggplot2::element_rect(fill = "white", color = "grey50",
                                                 linewidth = 0.4),
        plot.margin = ggplot2::margin(3, 4, 3, 4),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "none"
      )
    inset_grob <- ggplot2::ggplotGrob(p_zoom)

    # --- Main panel ---
    p <- ggplot2::ggplot(d,
           ggplot2::aes(x = age, y = R, color = model, linetype = model,
                        linewidth = model)) +
      ggplot2::geom_line() +
      ggplot2::scale_linewidth_manual(values = model_lw, guide = "none") +
      ggplot2::scale_color_manual(values = model_colors, name = NULL) +
      ggplot2::scale_linetype_manual(values = model_lty, name = NULL) +
      ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                          color = "gray80", linewidth = 0.3) +
      # Zoom region indicator
      ggplot2::annotate("rect", xmin = 75, xmax = 85,
                        ymin = -Inf, ymax = Inf,
                        alpha = 0.06, fill = "grey50") +
      ggplot2::labs(
        title = paste0(zyg, " twins"),
        x = "Age (years)",
        y = expression(R(t))
      ) +
      # Place inset in left-center (data coords: x 18–58, y upper ~65%)
      ggplot2::annotation_custom(
        grob = inset_grob,
        xmin = 18, xmax = 58, ymin = 0.35 * max(d$R, na.rm = TRUE),
        ymax = 1.02 * max(d$R, na.rm = TRUE)
      ) +
      theme_pub()
    if (!show_legend) {
      p <- p + ggplot2::theme(legend.position = "none")
    } else {
      p <- p + ggplot2::theme(
        legend.position = "bottom",
        legend.key.width = ggplot2::unit(2, "cm"))
    }
    p
  }

  p_mz <- make_panel("MZ", show_legend = FALSE)
  p_dz <- make_panel("DZ", show_legend = TRUE)

  p_mz / p_dz +
    patchwork::plot_annotation(
      title = "Bivariate dependence R(t)",
      theme = patchwork_title_theme()
    )
}


#' 2D bivariate surface residuals (appendix figure)
#'
#' Shows Delta R(t1, t2) = R_model(t1,t2) - R_true(t1,t2) for MZ twins,
#' comparing misspecified fit vs oracle fix side by side. Diverging color
#' scale centered at zero.
#'
#' @param biv_check List from run_bivariate_survival_check()
#' @return ggplot object (patchwork: Misspecified | Recovery: two-component refit)
plot_bivariate_surface <- function(biv_check) {
  biv_df <- biv_check$biv_diff
  if (nrow(biv_df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  max_abs <- max(abs(biv_df$delta_R), na.rm = TRUE)

  make_panel <- function(model_label, panel_title) {
    d <- biv_df[biv_df$model == model_label, ]
    ggplot2::ggplot(d,
      ggplot2::aes(x = t1, y = t2, fill = delta_R)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#2166ac", mid = "white", high = "#b2182b",
        midpoint = 0, limits = c(-max_abs, max_abs),
        name = expression(Delta * R)
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        title = panel_title,
        x = expression(t[1] ~ "(years)"),
        y = expression(t[2] ~ "(years)")
      ) +
      theme_pub()
  }

  p_misspec <- make_panel("Misspecified", "Misspecified \u2212 True DGP")
  p_fix     <- make_panel("Recovery: two-component refit",   "Recovery: two-component refit \u2212 True DGP")

  (p_misspec | p_fix) +
    patchwork::plot_annotation(
      title = expression("Bivariate dependence residuals " * Delta * R),
      theme = patchwork_title_theme()
    )
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
#' Plot MZ extrinsic concordance sensitivity
#'
#' @param mz_concordance Data frame from run_mz_concordance_sweep
#' @return ggplot object
plot_mz_concordance <- function(mz_concordance, mz_concordance_decoupled = NULL) {
  # Prepare coupled trace
  d <- mz_concordance[order(mz_concordance$r_gamma_mz), ]
  if ("n_reps" %in% names(d)) {
    d$t_crit <- qt(0.975, df = d$n_reps - 1)
  } else {
    d$t_crit <- qt(0.975, df = 19)
  }
  d$lo95 <- d$bias_pp - d$t_crit * d$se_pp
  d$hi95 <- d$bias_pp + d$t_crit * d$se_pp
  d$Trace <- "Coupled (\u03c1 = 0.4)"

  # Combine traces for legend
  plot_df <- d
  if (!is.null(mz_concordance_decoupled) &&
      "se_pp" %in% names(mz_concordance_decoupled)) {
    dd <- mz_concordance_decoupled[order(mz_concordance_decoupled$r_gamma_mz), ]
    if ("n_reps" %in% names(dd)) {
      dd$t_crit <- qt(0.975, df = dd$n_reps - 1)
    } else {
      dd$t_crit <- qt(0.975, df = 19)
    }
    dd$lo95 <- dd$bias_pp - dd$t_crit * dd$se_pp
    dd$hi95 <- dd$bias_pp + dd$t_crit * dd$se_pp
    dd$Trace <- "Decoupled (\u03c1 = 0)"
    plot_df <- rbind(plot_df[, names(plot_df) %in% names(dd)],
                     dd[, names(dd) %in% names(plot_df)])
  }
  plot_df$Trace <- factor(plot_df$Trace,
    levels = c("Coupled (\u03c1 = 0.4)", "Decoupled (\u03c1 = 0)"))

  trace_colors <- c("Coupled (\u03c1 = 0.4)" = "#2166ac",
                     "Decoupled (\u03c1 = 0)" = "#b2182b")
  trace_shapes <- c("Coupled (\u03c1 = 0.4)" = 16,
                     "Decoupled (\u03c1 = 0)" = 17)
  trace_linetypes <- c("Coupled (\u03c1 = 0.4)" = "solid",
                        "Decoupled (\u03c1 = 0)" = "dashed")

  ggplot2::ggplot(plot_df, ggplot2::aes(x = r_gamma_mz, y = bias_pp,
                                         color = Trace, fill = Trace,
                                         shape = Trace, linetype = Trace)) +
    # Plausible range band
    ggplot2::annotate("rect", xmin = 0.7, xmax = 0.9, ymin = -Inf, ymax = Inf,
                      fill = "#2166ac", alpha = 0.08) +
    ggplot2::annotate("text", x = 0.8, y = min(plot_df$bias_pp) + 0.5,
                      label = "Plausible\nrange", size = ANNOT_SIZE, color = "grey40") +
    # Negative-bias zone
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                      fill = "grey92", alpha = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    # CI ribbons
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95),
                         alpha = 0.15, color = NA) +
    # Lines and points
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2.5) +
    # Scales with legend
    ggplot2::scale_color_manual(values = trace_colors) +
    ggplot2::scale_fill_manual(values = trace_colors) +
    ggplot2::scale_shape_manual(values = trace_shapes) +
    ggplot2::scale_linetype_manual(values = trace_linetypes) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(
      title = "MZ concordance sensitivity",
      x = expression("Within-MZ extrinsic concordance (" * r[gamma*","*MZ] * ")"),
      y = "Net bias (pp)",
      color = NULL, fill = NULL, shape = NULL, linetype = NULL
    ) +
    theme_pub(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -5)
    )
}

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
                                 bivariate_check = NULL,
                                 mz_concordance = NULL,
                                 mz_concordance_decoupled = NULL,
                                 bridge_uncertainty = NULL,
                                 main_arms_replicated = NULL,
                                 robustness_cutoff0 = NULL,
                                 outdir = "figures") {
  paths <- list()

  # Main manuscript figures (replicated means with CIs)
  if (!is.null(main_arms_replicated)) {
    paths$main_arms <- save_figure(
      plot_main_arms(main_arms_replicated, robustness_cutoff0),
      "fig1_main_arms.png", width = 10, height = 6, outdir = outdir)
  } else {
    # Fallback to legacy single-seed plot
    paths$main_arms <- save_figure(
      plot_main_arms_legacy(summary_table),
      "fig1_main_arms.png", width = 10, height = 6, outdir = outdir)
  }

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
    "diag_mgg_hazard_curves.png", width = 10, height = 6, outdir = outdir)

  if (!is.null(mgg_param_comparison)) {
    paths$mgg_param <- save_figure(
      plot_mgg_param_comparison(mgg_param_comparison),
      "fig_mgg_parameterization.png", width = 14, height = 5.5, outdir = outdir)
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
      "fig_b3_mc_uncertainty.png", width = 10, height = 6.5, outdir = outdir)
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

  # Bivariate survival model check
  if (!is.null(bivariate_check)) {
    paths$bivariate_dep <- save_figure(
      plot_bivariate_dep_curves(bivariate_check),
      "fig_bivariate_dep_curves.png", width = 9, height = 10, outdir = outdir)

    paths$bivariate_surface <- save_figure(
      plot_bivariate_surface(bivariate_check),
      "fig_bivariate_surface.png", width = 11, height = 6, outdir = outdir)
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

  # MZ concordance sweep (with optional decoupled trace)
  if (!is.null(mz_concordance)) {
    paths$mz_concordance <- save_figure(
      plot_mz_concordance(mz_concordance, mz_concordance_decoupled),
      "fig_mz_concordance.png", width = 7, height = 5, outdir = outdir)
  }

  paths
}


# ===========================================================================
# Alternative Estimand Figures (Stage 11)
# ===========================================================================

# Shorthand for arm colors used in estimand figures
COL_ORACLE <- PAL[["Oracle"]]
COL_NULL   <- PAL[["Null"]]
COL_BIASED <- PAL[["Biased"]]
COL_FIX    <- PAL[["Fix"]]
MODEL_COLORS <- MODEL_PAL

#' σ_θ inflation bar chart (primary estimand figure)
plot_sigma_inflation <- function(vc) {
  gm <- vc[vc$model == "GM", ]
  gm$condition <- factor(gm$condition, levels = rev(gm$condition))

  fill_vals <- c(
    "Oracle" = COL_ORACLE,
    "Baseline (correctly specified)" = COL_NULL,
    "Misspecified (omitted extrinsic)" = COL_BIASED,
    "Recovery (two-component refit)" = COL_FIX
  )

  ggplot2::ggplot(gm, ggplot2::aes(x = condition, y = sigma_infl_pct,
                                    fill = condition)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%+.1f%%", sigma_infl_pct)),
      hjust = ifelse(gm$sigma_infl_pct >= 0, -0.1, 1.1),
      size = LABEL_SIZE) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::labs(
      x = NULL,
      y = expression(sigma[theta] ~ "inflation (%)"),
      title = expression("Primary estimand:" ~ Delta * sigma[theta])
    ) +
    theme_pub()
}

#' Multi-model σ_θ inflation comparison
plot_sigma_inflation_multimodel <- function(vc) {
  misspec <- vc[grepl("Misspecified", vc$condition), ]
  misspec$model <- factor(misspec$model, levels = c("GM", "MGG", "SR"))

  ggplot2::ggplot(misspec, ggplot2::aes(x = model, y = sigma_infl_pct,
                                         fill = model)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%+.1f%%", sigma_infl_pct)),
      vjust = -0.5, size = LABEL_SIZE) +
    ggplot2::scale_fill_manual(values = MODEL_COLORS, guide = "none") +
    ggplot2::labs(
      x = "Survival model",
      y = expression(sigma ~ "inflation (%)"),
      title = "Intrinsic frailty inflation across models"
    ) +
    theme_pub()
}

#' Bivariate survival surface contour plots (3-panel: true, misspec, recovery)
plot_bivariate_surface_contours <- function(surface_result, zygosity = "mz") {
  t_grid <- surface_result$t_grid

  make_panel <- function(S_obj, title) {
    df <- expand.grid(t1 = t_grid, t2 = t_grid)
    df$S <- as.vector(S_obj$S)
    ggplot2::ggplot(df, ggplot2::aes(x = t1, y = t2, fill = S)) +
      ggplot2::geom_raster(interpolate = TRUE) +
      ggplot2::geom_contour(ggplot2::aes(z = S), color = "white", alpha = 0.5,
                             breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)) +
      viridis::scale_fill_viridis(option = "D", limits = c(0, 1)) +
      ggplot2::labs(x = expression(t[1]), y = expression(t[2]),
                    title = title, fill = "S") +
      ggplot2::coord_fixed() +
      theme_pub(base_size = 12)
  }

  p1 <- make_panel(surface_result[[paste0("surface_true_", zygosity)]],
                    "True DGP")
  p2 <- make_panel(surface_result[[paste0("surface_misspec_", zygosity)]],
                    "Misspecified")

  fix_key <- paste0("surface_fix_", zygosity)
  if (!is.null(surface_result[[fix_key]])) {
    p3 <- make_panel(surface_result[[fix_key]], "Recovery")
    p1 + p2 + p3 + patchwork::plot_layout(ncol = 3, guides = "collect")
  } else {
    p1 + p2 + patchwork::plot_layout(ncol = 2, guides = "collect")
  }
}

#' Bivariate surface difference heatmap (misspec − true)
plot_bivariate_surface_diff <- function(surface_result, zygosity = "mz") {
  t_grid <- surface_result$t_grid
  S_true <- surface_result[[paste0("surface_true_", zygosity)]]
  S_misspec <- surface_result[[paste0("surface_misspec_", zygosity)]]

  delta <- S_misspec$S - S_true$S
  df <- expand.grid(t1 = t_grid, t2 = t_grid)
  df$delta_S <- as.vector(delta)
  lim <- max(abs(df$delta_S), na.rm = TRUE)

  ggplot2::ggplot(df, ggplot2::aes(x = t1, y = t2, fill = delta_S)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = delta_S),
                           color = "grey30", alpha = 0.5) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                   high = "#B2182B",
                                   midpoint = 0, limits = c(-lim, lim)) +
    ggplot2::labs(x = expression(t[1]), y = expression(t[2]),
                  title = paste0("Surface misfit (", toupper(zygosity), ")"),
                  fill = expression(Delta * S)) +
    ggplot2::coord_fixed() +
    theme_pub()
}

#' Concordance curves C(t) across models
plot_concordance_curves <- function(surface_result, zygosity = "mz") {
  dfs <- list()
  for (model in c("true", "misspec", "fix")) {
    key <- paste0("conc_", model, "_", zygosity)
    if (!is.null(surface_result[[key]])) {
      df <- surface_result[[key]]
      df$model <- switch(model,
        true = "True DGP", misspec = "Misspecified", fix = "Recovery")
      dfs[[model]] <- df
    }
  }
  combined <- do.call(rbind, dfs)
  combined <- combined[!is.na(combined$concordance), ]

  model_cols <- c("True DGP" = COL_ORACLE,
                  "Misspecified" = COL_BIASED,
                  "Recovery" = COL_FIX)

  ggplot2::ggplot(combined, ggplot2::aes(x = age, y = concordance,
                                          color = model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = if ("lo95" %in% names(combined)) lo95 else concordance - 1.96 * se,
      ymax = if ("hi95" %in% names(combined)) hi95 else concordance + 1.96 * se,
      fill = model),
      alpha = 0.15, color = NA) +
    ggplot2::scale_color_manual(values = model_cols) +
    ggplot2::scale_fill_manual(values = model_cols) +
    ggplot2::labs(
      x = "Age (years)",
      y = expression("C(t) = P(T"[2] * " > t | T"[1] * " > t)"),
      title = paste0("Concordance (", toupper(zygosity), ")"),
      color = NULL, fill = NULL
    ) +
    theme_pub()
}

#' ACE variance component stacked bars
plot_ace_components <- function(ace_analysis, transform_filter = "rank") {
  df <- ace_analysis$summary
  df <- df[df$transform == transform_filter &
           df$model_type == "ACE" &
           df$converged, ]
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  # Reshape to long
  long <- rbind(
    data.frame(arm = df$arm, component = "A (additive genetic)",
               proportion = df$a2, stringsAsFactors = FALSE),
    data.frame(arm = df$arm, component = "C (shared environment)",
               proportion = df$c2, stringsAsFactors = FALSE),
    data.frame(arm = df$arm, component = "E (unique environment)",
               proportion = df$e2, stringsAsFactors = FALSE)
  )
  long$component <- factor(long$component,
    levels = c("E (unique environment)",
               "C (shared environment)",
               "A (additive genetic)"))
  long$arm <- factor(long$arm, levels = unique(df$arm))

  ggplot2::ggplot(long, ggplot2::aes(x = arm, y = proportion,
                                      fill = component)) +
    ggplot2::geom_col(width = 0.7, position = "stack") +
    ggplot2::scale_fill_manual(values = c(
      "A (additive genetic)" = COL_ORACLE,
      "C (shared environment)" = "#D68910",
      "E (unique environment)" = "#999999"
    )) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = NULL, y = "Variance proportion",
                  title = "ACE decomposition by arm",
                  fill = NULL) +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Age-specific conditional correlation curves
plot_age_conditional_correlation <- function(age_dep, zygosity = "MZ") {
  df <- age_dep$correlations
  df <- df[df$zygosity == zygosity & !is.na(df$cor_pearson), ]
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  model_cols <- c(true_dgp = COL_ORACLE,
                  misspecified = COL_BIASED,
                  recovery = COL_FIX)
  model_labs <- c(true_dgp = "True DGP",
                  misspecified = "Misspecified",
                  recovery = "Recovery")

  ggplot2::ggplot(df, ggplot2::aes(x = age, y = cor_pearson,
                                    color = model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = model_cols, labels = model_labs) +
    ggplot2::labs(
      x = "Age threshold (years)",
      y = expression("cor(" * T[1] * "," ~ T[2] * " | both > t)"),
      title = paste0("Conditional correlation (", zygosity, ")"),
      color = NULL
    ) +
    theme_pub()
}

#' Cross-ratio comparison on diagonal (MZ)
plot_cross_ratio <- function(age_dep) {
  df <- age_dep$cross_ratio
  df <- df[!is.na(df$cross_ratio) & is.finite(df$cross_ratio) &
           df$cross_ratio > 0, ]
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  model_cols <- c(true_dgp = COL_ORACLE,
                  misspecified = COL_BIASED,
                  recovery = COL_FIX)
  model_labs <- c(true_dgp = "True DGP",
                  misspecified = "Misspecified",
                  recovery = "Recovery")

  ggplot2::ggplot(df, ggplot2::aes(x = age, y = cross_ratio,
                                    color = model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = model_cols, labels = model_labs) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Age (years)",
      y = expression(Theta(t, t)),
      title = expression("Cross-ratio" ~ Theta(t, t) ~ "(MZ)"),
      color = NULL
    ) +
    theme_pub()
}

#' Generate all alternative estimand figures
generate_estimand_figures <- function(variance_components,
                                      sweep_variance_decomp,
                                      variance_component_uq,
                                      bivariate_surface,
                                      ace_analysis,
                                      age_dependence,
                                      sigma_theta_true) {
  figdir <- "figures"
  if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
  W <- 7; H <- 5

  paths <- list()

  # E1: σ_θ inflation bar chart
  paths$e1 <- save_figure(
    plot_sigma_inflation(variance_components),
    "fig_e1_sigma_inflation.png", width = W, height = H, outdir = figdir)

  # E2: Multi-model σ_θ inflation
  paths$e2 <- save_figure(
    plot_sigma_inflation_multimodel(variance_components),
    "fig_e2_sigma_inflation_multimodel.png",
    width = W, height = 4, outdir = figdir)

  # E3: Raw bivariate surfaces — DROPPED (differences invisible on 0-1 scale;
  #      use E4 difference heatmap instead)

  # E4: Surface difference heatmap (the informative version)
  paths$e4 <- save_figure(
    plot_bivariate_surface_diff(bivariate_surface, "mz"),
    "fig_e4_surface_diff_mz.png",
    width = 6, height = 5, outdir = figdir)

  # E5: Concordance curves (MZ)
  paths$e5 <- save_figure(
    plot_concordance_curves(bivariate_surface, "mz"),
    "fig_e5_concordance_mz.png", width = W, height = H, outdir = figdir)

  # E6: ACE decomposition
  paths$e6 <- save_figure(
    plot_ace_components(ace_analysis),
    "fig_e6_ace_components.png", width = W, height = H, outdir = figdir)

  # E7: Age-specific conditional correlation (MZ)
  paths$e7 <- save_figure(
    plot_age_conditional_correlation(age_dependence, "MZ"),
    "fig_e7_age_cor_mz.png", width = W, height = H, outdir = figdir)

  # E8: Age-specific conditional correlation (DZ)
  paths$e8 <- save_figure(
    plot_age_conditional_correlation(age_dependence, "DZ"),
    "fig_e8_age_cor_dz.png", width = W, height = H, outdir = figdir)

  # E9: Cross-ratio — DROPPED (too noisy from finite-difference estimation
  #      of second derivatives on empirical step surface; conditional
  #      correlation E7/E8 shows the same pattern much more clearly)

  # E10: Composite (primary estimand + concordance + age correlation)
  p_left <- plot_sigma_inflation(variance_components)
  p_mid <- plot_concordance_curves(bivariate_surface, "mz")
  p_right <- plot_age_conditional_correlation(age_dependence, "MZ")
  p_composite <- p_left + p_mid + p_right +
    patchwork::plot_layout(ncol = 3) +
    patchwork::plot_annotation(tag_levels = "A")
  paths$e10 <- save_figure(
    p_composite,
    "fig_e10_estimand_composite.png",
    width = 14, height = 5, outdir = figdir)

  paths
}


# ===========================================================================
# σ-inflation-first sensitivity plots (CTR primary estimand versions)
# ===========================================================================

#' Crosswalk: Falconer h² bias vs σ_θ inflation across sweep grid
plot_crosswalk <- function(sweep_results) {
  df <- sweep_results[!is.na(sweep_results$sigma_infl_pct), ]
  if (nrow(df) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())

  ggplot2::ggplot(df, ggplot2::aes(x = sigma_infl_pct,
                                    y = 100 * bias)) +
    ggplot2::geom_point(ggplot2::aes(color = rho), size = 2, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "grey40",
                          linewidth = 0.8, linetype = "dashed") +
    viridis::scale_color_viridis(option = "C", name = expression(rho)) +
    ggplot2::labs(
      x = expression(sigma[theta] ~ "inflation (%)"),
      y = expression("Falconer" ~ h^2 ~ "bias (pp)"),
      title = expression("Crosswalk:" ~ sigma[theta] ~ "inflation" %->%
                          "Falconer" ~ h^2 ~ "bias")
    ) +
    theme_pub()
}

#' Negative ρ sweep — σ_θ inflation version (with CIs)
plot_negative_rho_sigma <- function(neg_rho_df) {
  if (!"sigma_infl_pct" %in% names(neg_rho_df))
    return(ggplot2::ggplot() + ggplot2::theme_void())

  ggplot2::ggplot(neg_rho_df, ggplot2::aes(x = rho, y = sigma_infl_pct)) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = sigma_infl_lo95, ymax = sigma_infl_hi95),
      fill = COL_BIASED, alpha = 0.15) +
    ggplot2::geom_line(color = COL_BIASED, linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.3) +
    ggplot2::labs(
      x = expression("Genetic correlation" ~ rho),
      y = expression(sigma[theta] ~ "inflation (%)"),
      title = expression("Sign reversal:" ~ sigma[theta] ~ "inflation vs" ~ rho)
    ) +
    theme_pub()
}

#' Negative ρ sweep — multi-model dispersion inflation version
plot_negative_rho_sigma_multimodel <- function(neg_rho_gm, model_controls) {
  gm <- neg_rho_gm[, c("rho", "sigma_infl_pct", "sigma_infl_lo95",
                         "sigma_infl_hi95")]
  gm$model <- "GM"

  mc <- model_controls[model_controls$sweep_type == "negative_rho", ]
  if (nrow(mc) > 0 && "disp_infl_pct" %in% names(mc)) {
    for (mdl in unique(mc$model)) {
      sub <- mc[mc$model == mdl, ]
      # Use t-critical with n_reps (default 20 seeds)
      n_reps_mc <- nrow(sub)
      tc <- if (n_reps_mc > 1) qt(0.975, df = n_reps_mc - 1) else NA_real_
      gm <- rbind(gm, data.frame(
        rho = sub$sweep_val,
        sigma_infl_pct = sub$disp_infl_pct,
        sigma_infl_lo95 = sub$disp_infl_pct - tc * sub$disp_infl_se,
        sigma_infl_hi95 = sub$disp_infl_pct + tc * sub$disp_infl_se,
        model = toupper(sub$model),
        stringsAsFactors = FALSE
      ))
    }
  }
  gm$model <- factor(gm$model, levels = c("GM", "MGG", "SR"))

  ggplot2::ggplot(gm, ggplot2::aes(x = rho, y = sigma_infl_pct,
                                    color = model, fill = model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = sigma_infl_lo95,
                                       ymax = sigma_infl_hi95),
                          alpha = 0.1, color = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_manual(values = MODEL_PAL) +
    ggplot2::scale_fill_manual(values = MODEL_PAL) +
    ggplot2::labs(
      x = expression("Genetic correlation" ~ rho),
      y = "Dispersion inflation (%)",
      title = "Sign reversal across mortality models",
      color = NULL, fill = NULL
    ) +
    theme_pub()
}

#' Dose-response — σ_θ inflation version
plot_dose_response_sigma <- function(dose_response) {
  dr <- if (is.data.frame(dose_response)) dose_response else dose_response$summary
  if (!"sigma_infl_pct" %in% names(dr))
    return(ggplot2::ggplot() + ggplot2::theme_void())

  ggplot2::ggplot(dr, ggplot2::aes(x = m_ex, y = sigma_infl_pct)) +
    ggplot2::geom_line(color = COL_BIASED, linewidth = 1) +
    ggplot2::geom_point(color = COL_BIASED, size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0.004, linetype = "dotted",
                         alpha = 0.3) +
    ggplot2::annotate("text", x = 0.0045, y = max(dr$sigma_infl_pct) * 0.3,
                       label = "historical", size = ANNOT_SIZE,
                       hjust = 0, color = "grey50") +
    ggplot2::labs(
      x = expression("Extrinsic mortality" ~ m[ex]),
      y = expression(sigma[theta] ~ "inflation (%)"),
      title = expression("Dose-response:" ~ sigma[theta] ~ "inflation vs" ~ m[ex])
    ) +
    theme_pub()
}

#' Heritable fraction sweep — σ_θ inflation version
plot_mex_split_sigma <- function(mex_split_hires) {
  if (!"sigma_infl_pct" %in% names(mex_split_hires))
    return(ggplot2::ggplot() + ggplot2::theme_void())

  agg <- aggregate(sigma_infl_pct ~ frac_heritable,
                   data = mex_split_hires, FUN = mean)
  se_agg <- aggregate(sigma_infl_pct ~ frac_heritable,
                      data = mex_split_hires,
                      FUN = function(x) sd(x) / sqrt(length(x)))
  n_agg <- aggregate(sigma_infl_pct ~ frac_heritable,
                     data = mex_split_hires, FUN = length)
  agg$se <- se_agg$sigma_infl_pct
  agg$t_crit <- qt(0.975, df = pmax(n_agg$sigma_infl_pct - 1, 1))

  ggplot2::ggplot(agg, ggplot2::aes(x = frac_heritable,
                                     y = sigma_infl_pct)) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = sigma_infl_pct - t_crit * se,
      ymax = sigma_infl_pct + t_crit * se),
      fill = COL_BIASED, alpha = 0.15) +
    ggplot2::geom_line(color = COL_BIASED, linewidth = 1) +
    ggplot2::geom_point(color = COL_BIASED, size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(
      x = "Heritable fraction of extrinsic mortality",
      y = expression(sigma[theta] ~ "inflation (%)"),
      title = expression("Heritable fraction:" ~ sigma[theta] ~ "inflation")
    ) +
    theme_pub()
}

#' MZ concordance sweep — σ_θ inflation version
plot_mz_concordance_sigma <- function(mz_concordance,
                                       mz_concordance_decoupled = NULL) {
  if (!"sigma_infl_pct" %in% names(mz_concordance))
    return(ggplot2::ggplot() + ggplot2::theme_void())

  p <- ggplot2::ggplot(mz_concordance, ggplot2::aes(
    x = r_gamma_mz, y = sigma_infl_pct)) +
    ggplot2::geom_line(color = COL_ORACLE, linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

  if (!is.null(mz_concordance_decoupled) &&
      "sigma_infl_pct" %in% names(mz_concordance_decoupled)) {
    p <- p + ggplot2::geom_line(
      data = mz_concordance_decoupled,
      ggplot2::aes(x = r_gamma_mz, y = sigma_infl_pct),
      color = COL_BIASED, linewidth = 1, linetype = "dashed")
  }

  p + ggplot2::labs(
    x = expression("Within-MZ extrinsic concordance" ~ r[gamma * ",MZ"]),
    y = expression(sigma[theta] ~ "inflation (%)"),
    title = expression(sigma[theta] ~ "inflation vs MZ concordance")
  ) + theme_pub()
}

#' σ_θ inflation heatmap over (ρ, σ_γ) grid
plot_sigma_inflation_heatmap <- function(sweep_results) {
  df <- sweep_results[!is.na(sweep_results$sigma_infl_pct), ]

  lim <- max(abs(df$sigma_infl_pct), na.rm = TRUE)

  # Axis order matches existing plot_sweep_heatmap (x=sigma_gamma, y=rho)
  ggplot2::ggplot(df, ggplot2::aes(x = sigma_gamma, y = rho,
                                    fill = sigma_infl_pct)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-lim, lim),
      name = expression(sigma[theta] ~ "inflation (%)")) +
    ggplot2::labs(
      x = expression("Extrinsic heterogeneity" ~ sigma[gamma]),
      y = expression("Genetic correlation" ~ rho),
      title = expression(sigma[theta] ~ "inflation across sensitivity grid")
    ) +
    theme_pub()
}

# ---------------------------------------------------------------------------
# ML vs moment estimation comparison plots
# ---------------------------------------------------------------------------

#' ML vs Moment sigma_theta scatter plot
#'
#' Each point is one MC replicate. Points above the identity line mean
#' ML inflates less than moment.
plot_ml_vs_moment_scatter <- function(ml_results) {
  df <- ml_results$results
  df$label_f <- factor(df$label,
    levels = c("sanity", "default", "no_pleiotropy", "extreme"),
    labels = c("Correctly specified", "Misspecified (default)",
               "No pleiotropy", "Extreme heterogeneity"))

  # Use project palette: green for correct, red for misspec, purple/orange for variants
  cond_colors <- c(
    "Correctly specified"    = PAL[["Null"]],
    "Misspecified (default)" = PAL[["Biased"]],
    "No pleiotropy"          = PAL[["Pleiotropy"]],
    "Extreme heterogeneity"  = PAL[["Hamilton"]]
  )
  cond_shapes <- c(
    "Correctly specified" = 16, "Misspecified (default)" = 17,
    "No pleiotropy" = 15, "Extreme heterogeneity" = 18
  )

  true_sigma <- unique(df$sigma_theta_true)[1]

  ggplot2::ggplot(df, ggplot2::aes(x = sigma_theta_moment, y = sigma_theta_ml,
                                     color = label_f, shape = label_f)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "grey60", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = true_sigma, linetype = "dotted",
                         color = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = true_sigma, linetype = "dotted",
                         color = "grey40", linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.6, size = 2.2) +
    ggplot2::scale_color_manual(values = cond_colors) +
    ggplot2::scale_shape_manual(values = cond_shapes) +
    ggplot2::annotate("text", x = true_sigma, y = Inf,
                       label = expression("true" ~ sigma[theta]),
                       hjust = 1.1, vjust = 1.5, size = ANNOT_SIZE,
                       color = "grey40", fontface = "italic") +
    ggplot2::labs(
      x = expression(hat(sigma)[theta]^"moment"),
      y = expression(hat(sigma)[theta]^"ML"),
      color = NULL, shape = NULL,
      title = "ML vs moment calibration"
    ) +
    theme_pub() +
    ggplot2::theme(legend.position = "right")
}

#' Convergence strip plot for ML and moment estimates
#'
#' Shows individual estimates as jittered points with mean + 95% CI bars,
#' matching the style of plot_mc_uncertainty_multimodel.
plot_ml_convergence <- function(ml_results, condition = "default") {
  df <- ml_results$results
  df <- df[df$label == condition, ]

  true_sigma <- unique(df$sigma_theta_true)[1]

  # Reshape to long
  df_long <- rbind(
    data.frame(n_pairs = df$n_pairs, estimator = "ML",
               sigma_hat = df$sigma_theta_ml, stringsAsFactors = FALSE),
    data.frame(n_pairs = df$n_pairs, estimator = "Moment",
               sigma_hat = df$sigma_theta_moment, stringsAsFactors = FALSE)
  )
  df_long$n_label <- paste0("n = ", formatC(df_long$n_pairs,
                                              format = "d", big.mark = ","))
  df_long$n_label <- factor(df_long$n_label,
    levels = paste0("n = ", formatC(sort(unique(df_long$n_pairs)),
                                      format = "d", big.mark = ",")))
  df_long$estimator <- factor(df_long$estimator, levels = c("ML", "Moment"))

  # Compute summary per group
  summ <- do.call(rbind, lapply(
    split(df_long, interaction(df_long$n_label, df_long$estimator)),
    function(d) {
      n <- nrow(d)
      m <- mean(d$sigma_hat)
      se <- sd(d$sigma_hat) / sqrt(n)
      t_crit <- if (n > 1) qt(0.975, n - 1) else 1.96
      data.frame(n_label = d$n_label[1], estimator = d$estimator[1],
                 mean = m, lo = m - t_crit * se, hi = m + t_crit * se,
                 stringsAsFactors = FALSE)
    }
  ))

  est_colors <- c("ML" = PAL[["Oracle"]], "Moment" = PAL[["Biased"]])

  ggplot2::ggplot(df_long, ggplot2::aes(x = sigma_hat, y = estimator,
                                          color = estimator)) +
    ggplot2::facet_wrap(~ n_label, ncol = 1) +
    # True sigma reference
    ggplot2::geom_vline(xintercept = true_sigma, linetype = "dotted",
                         color = "grey40", linewidth = 0.5) +
    # Individual seeds (faint)
    ggplot2::geom_jitter(alpha = 0.4, size = 1.8, height = 0.15, width = 0) +
    # Mean + CI
    ggplot2::geom_pointrange(data = summ,
      ggplot2::aes(x = mean, xmin = lo, xmax = hi, y = estimator),
      color = "grey10", size = 0.5, linewidth = 0.8, fatten = 3) +
    # Mean label
    ggplot2::geom_text(data = summ,
      ggplot2::aes(x = mean, y = estimator,
                    label = sprintf("%.3f", mean)),
      vjust = -1.2, size = LABEL_SIZE, color = "grey20") +
    ggplot2::scale_color_manual(values = est_colors, guide = "none") +
    ggplot2::labs(
      x = expression(hat(sigma)[theta]),
      y = NULL,
      title = "Convergence to pseudo-true",
      subtitle = expression(paste("Default DGP: ", sigma[gamma], " = 0.40, ",
                                    rho, " = 0.4"))
    ) +
    theme_pub() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.95))
    )
}
