# ===========================================================================
# commentary_figure.R — Single multi-panel commentary figure
# ===========================================================================
# Layout (full page, 3 rows × 2 cols):
#   Row 1: (A) Main results forest plot    | (B) Variance decomposition
#   Row 2: (C) Three-model comparison      | (D) Sign reversal
#   Row 3: (E) Anchored heatmap            | (F) Heritable fraction
# ===========================================================================

library(ggplot2)
library(patchwork)

source("R/params.R")
source("R/plots.R")

# --- Load all needed data from targets store ---
main_rep     <- targets::tar_read(main_arms_replicated)
cutoff0      <- targets::tar_read(robustness_cutoff0)
mc_multi     <- targets::tar_read(mc_uncertainty_multimodel)
anchored     <- targets::tar_read(anchored_results)
mex_hires    <- targets::tar_read(mex_split_hires)
var_decomp   <- readRDS("_targets/objects/var_decomp")
alpha_beta   <- targets::tar_read(empirical_alpha_beta)
neg_rho      <- targets::tar_read(negative_rho)
neg_rho_mgg  <- targets::tar_read(negative_rho_mgg)
neg_rho_sr   <- targets::tar_read(negative_rho_sr)

# ===================================================================
# Panel A: Main results (streamlined forest plot — bias panel only)
# Colored seeds + boxplots per condition (matching fig1 style),
# with subtle highlight on Misspecified row.
# ===================================================================
make_panel_a <- function() {
  s <- main_rep$summary
  ps <- main_rep$per_seed

  conditions <- c("Correctly specified", "Misspecified", "Nonfamilial extrinsic",
                   "Irrelevant trait", "Vanishing extrinsic",
                   "Pleiotropy only (\u03C1=0)", "Two-component recovery", "No survivor cutoff")
  metrics <- c("arm1_bias_pp", "arm2_bias_pp", "ctrl_a_bias_pp", "irrel_bias_pp",
               "ctrl_b_bias_pp", "pleio_iso_bias_pp", "fix_bias_pp", NA)
  pal_keys <- c("Null", "Biased", "ControlA", "Irrelevant",
                "ControlB", "Pleiotropy", "Fix", "Cutoff0")
  colors <- unname(PAL[pal_keys])

  df <- do.call(rbind, lapply(seq_along(conditions), function(i) {
    if (conditions[i] == "No survivor cutoff") {
      if (is.null(cutoff0)) return(NULL)
      return(data.frame(
        Condition = "No survivor cutoff", mean = cutoff0$bias_pp,
        lo95 = cutoff0$lo95_pp, hi95 = cutoff0$hi95_pp,
        stringsAsFactors = FALSE
      ))
    }
    row <- s[s$metric == metrics[i], ]
    if (nrow(row) == 0) return(NULL)
    data.frame(
      Condition = conditions[i],
      mean = row$mean, lo95 = row$lo95, hi95 = row$hi95,
      stringsAsFactors = FALSE
    )
  }))

  # Build color map keyed by condition name
  cond_color_map <- setNames(colors, conditions)
  condition_colors <- cond_color_map[df$Condition]

  display_order <- rev(df$Condition)
  df$Condition <- factor(df$Condition, levels = display_order)

  # Per-seed jitter data
  seed_rows <- list()
  for (i in seq_along(conditions)) {
    if (conditions[i] == "No survivor cutoff") {
      if (!is.null(cutoff0) && !is.null(cutoff0$per_seed) &&
          "bias_pp" %in% names(cutoff0$per_seed)) {
        vals <- cutoff0$per_seed$bias_pp
        if (length(vals) > 0) {
          seed_rows <- c(seed_rows, list(data.frame(
            Condition = "No survivor cutoff", val = vals, stringsAsFactors = FALSE
          )))
        }
      }
      next
    }
    vals <- ps[[metrics[i]]]
    if (!is.null(vals) && length(vals) > 0) {
      seed_rows <- c(seed_rows, list(data.frame(
        Condition = conditions[i], val = vals, stringsAsFactors = FALSE
      )))
    }
  }
  seed_long <- do.call(rbind, seed_rows)
  seed_long$Condition <- factor(seed_long$Condition, levels = display_order)

  # Find misspecified position for highlight rect
  misspec_pos <- which(display_order == "Misspecified")

  ggplot(df, aes(x = Condition, y = mean)) +
    # Subtle gray highlight on Misspecified row
    annotate("rect",
      xmin = misspec_pos - 0.5, xmax = misspec_pos + 0.5,
      ymin = -Inf, ymax = Inf,
      fill = "gray90", color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_boxplot(data = seed_long,
      aes(x = Condition, y = val, fill = Condition),
      width = 0.4, alpha = 0.35, color = "gray50",
      outlier.shape = NA, inherit.aes = FALSE, show.legend = FALSE) +
    geom_jitter(data = seed_long,
      aes(x = Condition, y = val, color = Condition),
      width = 0.15, alpha = 0.3, size = 0.8,
      inherit.aes = FALSE, show.legend = FALSE) +
    geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.7, color = "gray20") +
    geom_point(size = 1.8, color = "gray10") +
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(limits = display_order) +
    coord_flip() +
    labs(title = "Bias localizes to the intended confounding pathway",
         subtitle = "50-seed mean \u00b1 95% CI (No cutoff: 20 seeds)",
         x = NULL, y = "Bias (percentage points)") +
    theme_pub(base_size = 13) +
    theme(plot.margin = margin(6, 12, 6, 6))
}

# ===================================================================
# Panel B: Variance decomposition (pseudo-parameter tracking)
# ===================================================================
make_panel_b <- function() {
  df <- var_decomp$data

  lbl <- sprintf("\u03B1 = %.2f,  \u03B2 = %.2f\nR\u00B2 = %.3f\nMAR = %.1f%%",
                 alpha_beta$alpha, alpha_beta$beta,
                 alpha_beta$r_squared, var_decomp$mar_pct)

  ggplot(df, aes(x = predicted_var, y = actual_var)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(aes(color = factor(rho)), size = 2) +
    scale_color_manual(values = viridis::viridis(length(unique(df$rho)))) +
    annotate("text", x = min(df$predicted_var),
             y = max(df$actual_var) * 0.95,
             label = lbl, hjust = 0, size = ANNOT_SIZE, family = "mono") +
    labs(title = "Fitted variance tracks omitted-variable absorption",
         subtitle = expression("Predicted vs actual " ~ sigma[fit]^2 ~ " across " ~ rho %*% sigma[gamma] ~ " grid"),
         x = expression("Predicted " ~ sigma[fit]^2),
         y = expression("Actual " ~ sigma[fit]^2),
         color = expression(rho)) +
    theme_pub(base_size = 13) +
    theme(legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          plot.margin = margin(6, 6, 6, 6))
}

# ===================================================================
# Panel C: Three-model comparison (forest with per-seed jitter)
# ===================================================================
make_panel_c <- function() {
  mm <- mc_multi$summary
  ps <- mc_multi$per_seed
  t_crit <- qt(0.975, df = mm$n_seeds[1] - 1)

  # Summary data
  df <- data.frame(
    Model = factor(rep(c("GM", "MGG", "SR"), 2), levels = c("SR", "MGG", "GM")),
    Arm = rep(c("Baseline", "Misspecified"), each = 3),
    mean = c(mm$arm1_bias_mean, mm$arm2_bias_mean),
    se = c(mm$arm1_bias_se, mm$arm2_bias_se),
    stringsAsFactors = FALSE
  )
  df$lo <- df$mean - t_crit * df$se
  df$hi <- df$mean + t_crit * df$se

  # Per-seed data for jitter + boxplots
  seed_long <- rbind(
    data.frame(Model = factor(ps$model, levels = c("SR", "MGG", "GM")),
               Arm = "Baseline", val = ps$arm1_bias_pp, stringsAsFactors = FALSE),
    data.frame(Model = factor(ps$model, levels = c("SR", "MGG", "GM")),
               Arm = "Misspecified", val = ps$arm2_bias_pp, stringsAsFactors = FALSE)
  )

  cols <- c("Baseline" = unname(PAL["Null"]), "Misspecified" = unname(PAL["Biased"]))

  ggplot(df, aes(x = Model, y = mean, color = Arm)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_boxplot(data = seed_long,
      aes(x = Model, y = val, fill = Arm, group = interaction(Model, Arm)),
      width = 0.5, alpha = 0.25, color = "gray50",
      outlier.shape = NA, inherit.aes = FALSE,
      position = position_dodge(width = 0.6), show.legend = FALSE) +
    geom_jitter(data = seed_long,
      aes(x = Model, y = val, color = Arm, group = interaction(Model, Arm)),
      size = 0.8, alpha = 0.35,
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6),
      inherit.aes = FALSE, show.legend = FALSE) +
    geom_pointrange(aes(ymin = lo, ymax = hi),
                    position = position_dodge(width = 0.6),
                    size = 0.5, linewidth = 0.6) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(title = "Upward bias replicates across mortality families",
         subtitle = "20 seeds per model; mean \u00b1 95% CI",
         x = NULL, y = "Bias (pp)", color = NULL) +
    theme_pub(base_size = 13) +
    theme(legend.position = c(0.80, 0.08),
          legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
          legend.key.size = unit(0.4, "cm"),
          plot.margin = margin(6, 12, 6, 6))
}

# ===================================================================
# Panel D: Sign reversal (multimodel negative rho)
# Uses MODEL_PAL, MODEL_SHAPES, MODEL_LINETYPES from plots.R
# ===================================================================
make_panel_d <- function() {
  # GM data
  gm_df <- data.frame(
    rho = neg_rho$rho, bias_pp = neg_rho$bias_pp,
    model = "GM", lo = neg_rho$lo95_pp, hi = neg_rho$hi95_pp,
    stringsAsFactors = FALSE
  )

  # MGG data
  mgg_df <- data.frame(
    rho = neg_rho_mgg$sweep_val, bias_pp = neg_rho_mgg$bias_pp,
    model = "MGG", lo = neg_rho_mgg$lo95_pp, hi = neg_rho_mgg$hi95_pp,
    stringsAsFactors = FALSE
  )

  # SR data
  sr_df <- data.frame(
    rho = neg_rho_sr$sweep_val, bias_pp = neg_rho_sr$bias_pp,
    model = "SR", lo = neg_rho_sr$lo95_pp, hi = neg_rho_sr$hi95_pp,
    stringsAsFactors = FALSE
  )

  all_df <- rbind(gm_df, mgg_df, sr_df)
  all_df$model <- factor(all_df$model, levels = c("GM", "MGG", "SR"))

  ggplot(all_df, aes(x = rho, y = bias_pp,
                     color = model, shape = model,
                     linetype = model, fill = model)) +
    # Subtle gray highlight on lower-left sector (ρ < 0, bias < 0)
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = "gray90", color = NA) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                alpha = 0.15, linewidth = 0, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray60", linewidth = 0.3) +
    scale_color_manual(values = MODEL_PAL) +
    scale_fill_manual(values = MODEL_PAL, guide = "none") +
    scale_shape_manual(values = MODEL_SHAPES) +
    scale_linetype_manual(values = MODEL_LINETYPES) +
    labs(title = "Bias reverses under negative genetic correlation",
         subtitle = "Bias vs \u03C1 (20 seeds/point; ribbon = 95% CI)",
         x = expression("Genetic correlation " * rho),
         y = "Bias (pp)", color = "Model", shape = "Model", linetype = "Model") +
    theme_pub(base_size = 13) +
    theme(legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
          legend.key.size = unit(0.4, "cm"),
          plot.margin = margin(6, 12, 6, 6))
}

# ===================================================================
# Panel E: Anchored heatmap (viridis, matching fig3 axis orientation)
# σ_γ on x-axis, ρ on y-axis
# ===================================================================
make_panel_e <- function() {
  df <- as.data.frame(anchored)
  df$bias_pp <- df$bias * 100

  # Dynamic text color: white on dark tiles, black on light tiles
  mid <- (min(df$bias_pp) + max(df$bias_pp)) / 2
  df$text_col <- ifelse(df$bias_pp < mid, "white", "gray10")

  ggplot(df, aes(x = factor(sigma_gamma), y = factor(rho), fill = bias_pp)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%+.1f", bias_pp), color = text_col),
              size = 3, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    viridis::scale_fill_viridis(option = "viridis", name = "Bias\n(pp)",
                                limits = c(min(df$bias_pp), max(df$bias_pp))) +
    labs(title = "Anchored sensitivity regime",
         subtitle = expression(sigma[gamma] ~ "\u2208 [0.30, 0.65], " ~ rho ~ "\u2208 [0.20, 0.50]"),
         x = expression(sigma[gamma]), y = expression(rho)) +
    theme_pub(base_size = 13) +
    theme(legend.position = "right",
          legend.key.width = unit(0.3, "cm"),
          legend.key.height = unit(1, "cm"),
          plot.margin = margin(6, 6, 6, 6),
          panel.grid = element_blank())
}

# ===================================================================
# Panel F: Heritable fraction sensitivity
# ===================================================================
make_panel_f <- function() {
  df <- mex_hires
  agg <- aggregate(bias_pp ~ frac_heritable, data = df, FUN = mean)
  se_df <- aggregate(bias_pp ~ frac_heritable, data = df, FUN = function(x) {
    if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
  })
  names(se_df)[2] <- "se"
  agg <- merge(agg, se_df, by = "frac_heritable")
  n_per <- max(table(df$frac_heritable))
  t_crit <- qt(0.975, df = n_per - 1)
  agg$lo <- agg$bias_pp - t_crit * agg$se
  agg$hi <- agg$bias_pp + t_crit * agg$se

  ggplot(agg, aes(x = frac_heritable, y = bias_pp)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = PAL["Biased"], alpha = 0.2) +
    geom_jitter(data = df, aes(x = frac_heritable, y = bias_pp),
                width = 0.008, alpha = 0.2, size = 0.5,
                color = PAL["Biased"]) +
    geom_line(linewidth = 0.7, color = PAL["Biased"]) +
    geom_point(size = 1.5, color = PAL["Biased"]) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(title = "Bias increases with heritable extrinsic fraction",
         subtitle = paste0(n_per, " seeds/level; ribbon = 95% CI"),
         x = expression(f[heritable] == m[inf] / m[ex]),
         y = "Bias (pp)") +
    theme_pub(base_size = 13) +
    theme(plot.margin = margin(6, 12, 6, 6))
}

# ===================================================================
# Compose and save
# ===================================================================
cat("Building panels...\n")
pA <- make_panel_a(); cat("  A done\n")
pB <- make_panel_b(); cat("  B done\n")
pC <- make_panel_c(); cat("  C done\n")
pD <- make_panel_d(); cat("  D done\n")
pE <- make_panel_e(); cat("  E done\n")
pF <- make_panel_f(); cat("  F done\n")

# 3×2 grid layout (use / for vertical stacking, | for side-by-side)
composite <- (
  (pA | pB) /
  (pC | pD) /
  (pE | pF)
) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

outfile <- "figures/commentary_figure.png"
ggsave(outfile, composite, width = 14, height = 17, dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", outfile))

# ===================================================================
# Alternate version: mechanism diagram + 4 panels (no E/F)
# Row 1: mechanism.jpg (full width)
# Row 2: A | B
# Row 3: C | D
# ===================================================================
cat("Building alternate (mechanism + 4 panels)...\n")

# Wrap mechanism.png as a ggplot panel
# Read image and compute aspect ratio
mech_img <- png::readPNG("figures/unnamed.png")
img_aspect <- ncol(mech_img) / nrow(mech_img)
# Place image left-aligned at native aspect ratio, 20% larger than panels below
# Image occupies ~85% of width, left-aligned with space on the right
img_w <- 0.95
img_h <- img_w / img_aspect
p_mech <- ggplot() +
  annotation_raster(mech_img, xmin = 0, xmax = img_w, ymin = 0, ymax = img_h) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.01, img_h * 0.98), expand = FALSE, clip = "off") +
  labs(title = "Omitted extrinsic frailty is absorbed into the calibrated intrinsic component") +
  theme_pub(base_size = 13) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.line = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(margin = margin(t = 0, b = 2)),
        plot.margin = margin(0, 6, 0, 6))

composite2 <- (
  p_mech /
  (pA | pB) /
  (pC | pD)
) +
  plot_layout(heights = c(1.3, 1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 4, 4, 4)
    )
  )

outfile2 <- "figures/commentary_figure_v2.png"
ggsave(outfile2, composite2, width = 14, height = 15.5, dpi = 300, bg = "white")
cat(sprintf("Saved: %s\n", outfile2))
