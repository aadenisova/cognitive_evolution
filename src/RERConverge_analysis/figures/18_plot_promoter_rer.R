#!/usr/bin/env Rscript
# 02_plot_promoter_rer_correlation.R
# Faceted scatter plot: RER vs innovation for top 10 pos + 10 neg promoter CNEEs.
# Each panel: one CNEE, points = species, regression line, Rho/P/dist in title.
# Negative Rho panels first (by |Rho| desc), then positive Rho (by Rho desc).

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
})

DATA_DIR  <- "/mnt/results/plotting/data"
OUTPUT_DIR <- "/mnt/results/figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Color palette
COL_POS <- "#E94057"
COL_NEG <- "#0279EE"

# ── Load data ───────────────────────────────────────────────────────────────
message("Loading data...")
df <- read.csv(file.path(DATA_DIR, "promoter_rer_data_20.csv"), stringsAsFactors = FALSE)

# Build facet labels: geneSymbol (CNEE)\nRho=X.XXX, P=X.XXe-XX | dist=Xbp
df$facet_label <- sprintf("%s (%s)\nRho=%.3f, P=%.2e | dist=%dbp",
                          df$geneSymbol, df$CNEE,
                          df$Rho, df$P, df$abs_distance_to_TSS)

# Order facets: negative first (by |Rho| desc), then positive (by Rho desc)
facet_order <- unique(df[, c("CNEE", "facet_order", "direction", "Rho", "facet_label")])
facet_order <- facet_order[order(facet_order$facet_order), ]
df$facet_label <- factor(df$facet_label, levels = facet_order$facet_label)

# Direction color mapping
df$dir_color <- ifelse(df$direction == "positive", COL_POS, COL_NEG)

# ── Plot ─────────────────────────────────────────────────────────────────────
message("Building faceted scatter plot...")

p <- ggplot(df, aes(x = RER, y = TECHINNO2025_ResEff)) +
  geom_point(aes(color = direction), alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", fill = "grey80",
              linewidth = 0.5, alpha = 0.3) +
  facet_wrap(~ facet_label, ncol = 4, scales = "free") +
  scale_color_manual(values = c("positive" = COL_POS, "negative" = COL_NEG),
                     guide = "none") +
  labs(x = "Relative Evolutionary Rate (RER)",
       y = "Technological Innovation (TECHINNO2025_ResEff)") +
  theme_prism(base_size = 10, base_family = "sans") +
  theme(
    strip.text = element_text(size = 7, face = "bold", lineheight = 0.9),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.35, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

# ── Save ────────────────────────────────────────────────────────────────────
png_path <- file.path(OUTPUT_DIR, "figure_promoter_rer_correlation.png")
svg_path <- file.path(OUTPUT_DIR, "figure_promoter_rer_correlation.svg")

ggsave(png_path, p, width = 16, height = 12, dpi = 300, bg = "white")
ggsave(svg_path, p, width = 16, height = 12, bg = "white")
message(sprintf("Saved: %s", png_path))
message(sprintf("Saved: %s", svg_path))

message("\n=== Step 02 complete ===")
