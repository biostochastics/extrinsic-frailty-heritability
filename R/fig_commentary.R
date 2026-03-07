# Generate the Science commentary figure
# Three panels: (A) Core bias, (B) Dose-response, (C) Sign reversal
#
# Usage: source("R/fig_commentary.R"); generate_commentary_figure()

generate_commentary_figure <- function(outdir = "figures") {
  library(ggplot2)
  library(patchwork)

  targets::tar_load(c(
    "oracle", "arm1_corr", "arm2_corr",
    "multi_target_arm", "robustness_cutoff0",
    "dose_response", "negative_rho"
  ))

  oracle_h2 <- oracle$h2

  # Theme
  theme_sci <- function() {
    theme_minimal(base_size = 9.5) +
      theme(
        plot.title = element_text(face = "bold", size = 9.5, hjust = 0,
                                  margin = margin(b = 6)),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray92"),
        axis.line = element_line(color = "gray40", linewidth = 0.3),
        axis.ticks = element_line(color = "gray40", linewidth = 0.3),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 7.5),
        plot.margin = margin(6, 12, 4, 4)
      )
  }

  pal_blue <- "#1B6B93"
  pal_green <- "#2E8B57"
  pal_red <- "#C73E3A"
  pal_orange <- "#D68910"

  # ── Panel A: Core bias comparison ──
  bar_df <- data.frame(
    Condition = c("Oracle", "Correct\nspec.", "Misspec.\n(r_MZ)", "Misspec.\n(5-target)"),
    h2 = c(oracle_h2, arm1_corr$h2, arm2_corr$h2, multi_target_arm$mt_h2),
    fill = c("oracle", "correct", "misspec", "multi"),
    stringsAsFactors = FALSE
  )
  bar_df$bias_pp <- round(100 * (bar_df$h2 - oracle_h2), 1)
  bar_df$label <- ifelse(bar_df$bias_pp == 0,
    sprintf("%.2f", bar_df$h2),
    sprintf("%.2f\n(%+.1f pp)", bar_df$h2, bar_df$bias_pp))
  bar_df$Condition <- factor(bar_df$Condition, levels = bar_df$Condition)

  fill_vals <- c(oracle = pal_blue, correct = pal_green,
                 misspec = pal_red, multi = pal_orange)

  pA <- ggplot(bar_df, aes(x = Condition, y = h2, fill = fill)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = oracle_h2, linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    geom_text(aes(label = label), vjust = -0.2, size = 2.6,
              color = "gray25", lineheight = 0.85) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_y_continuous(limits = c(0, 0.73), expand = c(0, 0)) +
    labs(title = "A  Calibration bias",
         x = NULL, y = expression(hat(h)^2)) +
    theme_sci() +
    theme(axis.text.x = element_text(size = 7, lineheight = 0.85))

  # ── Panel B: Dose-response ──
  dr <- dose_response
  dr$bias_pp <- 100 * dr$bias

  pB <- ggplot(dr, aes(x = m_ex * 1000, y = bias_pp)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray60", linewidth = 0.3) +
    geom_line(color = pal_red, linewidth = 0.8) +
    geom_point(color = pal_red, size = 1.8) +
    geom_vline(xintercept = 4, linetype = "dotted",
               color = "gray50", linewidth = 0.4) +
    annotate("text", x = 4.3, y = -3.5,
             label = "historical", hjust = 0, size = 2.4, color = "gray45") +
    labs(title = "B  Dose-response",
         x = expression(m[ex] ~ "(per 1000 / yr)"),
         y = "Bias (pp)") +
    theme_sci()

  # ── Panel C: Sign reversal ──
  nr <- negative_rho

  pC <- ggplot(nr, aes(x = rho, y = bias_pp)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "gray50", linewidth = 0.4) +
    geom_line(color = pal_red, linewidth = 0.8) +
    geom_point(color = pal_red, size = 1.8) +
    annotate("text", x = -0.35, y = -5,
             label = "\u03c1 < 0\nbias \u2193",
             size = 2.5, color = pal_blue, fontface = "bold",
             lineheight = 0.85) +
    annotate("text", x = 0.45, y = 14,
             label = "\u03c1 > 0\nbias \u2191",
             size = 2.5, color = pal_red, fontface = "bold",
             lineheight = 0.85) +
    labs(title = "C  Sign reversal",
         x = expression(rho ~ "(genetic correlation)"),
         y = "Bias (pp)") +
    theme_sci()

  # ── Compose ──
  fig <- pA + pB + pC +
    plot_layout(ncol = 3, widths = c(1.05, 1, 1))

  ggsave(file.path(outdir, "fig_commentary.pdf"),
         fig, width = 18, height = 8.5, units = "cm", device = cairo_pdf)
  ggsave(file.path(outdir, "fig_commentary.png"),
         fig, width = 18, height = 8.5, units = "cm", dpi = 300)

  message("Saved: ", file.path(outdir, "fig_commentary.pdf"))
  message("Saved: ", file.path(outdir, "fig_commentary.png"))
  invisible(fig)
}
