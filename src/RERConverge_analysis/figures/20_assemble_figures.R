#!/usr/bin/env Rscript
# 10_assemble_figures.R
# Cowplot assembly of all individual figures into a combined publication figure.
# Sources 09_figures.R to obtain plot objects, then arranges with cowplot.
#
# Layout:
#   +-------------------+-----------+
#   |                   |   B1      |
#   |        A          |  (intron) |
#   |   (tree + innov)  +-----------+
#   |                   |   B2      |
#   |                   | (promoter)|
#   +-------------------+-----------+
#   |        D          |     S     |
#   |   (DAR enrich)    |  (element)|
#   +-------------------+-----------+

# Workaround: ggtree 3.14.0 calls is.waive() removed from ggplot2 4.0.3
is.waive <- function(x) inherits(x, "waiver")

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplot2)
  library(patchwork)
})

OUTPUT_DIR <- "/mnt/results/figures"

# Source individual figures script — creates p_a, p_b1, p_b2, p_d, p_s
# and re-saves individual PNG/SVG files
message("Sourcing 09_figures.R for plot objects...")
source("/mnt/results/plotting/code/09_figures.R")
message("Plot objects loaded.")

# Convert patchwork p_a to a form cowplot can handle
p_a_wrapped <- wrap_elements(full = p_a)

# Right column: B1 (top) + B2 (bottom)
right_col <- plot_grid(
  p_b1, p_b2,
  ncol = 1,
  rel_heights = c(1.6, 1),
  labels = c("B1", "B2"),
  label_size = 14
)

# Top section: A (tree, left, wide) | right column
top_section <- plot_grid(
  p_a_wrapped, right_col,
  ncol = 2,
  rel_widths = c(1.6, 1),
  labels = c("A", ""),
  label_size = 14
)

# Bottom section: D (left) | S (right)
bottom_section <- plot_grid(
  p_d, p_s,
  ncol = 2,
  rel_widths = c(1.1, 1),
  labels = c("D", "S"),
  label_size = 14
)

# Full assembly
combined <- plot_grid(
  top_section, bottom_section,
  ncol = 1,
  rel_heights = c(3, 1.3)
)

# Save combined figure
png_path <- file.path(OUTPUT_DIR, "figure_combined.png")
svg_path <- file.path(OUTPUT_DIR, "figure_combined.svg")
ggsave(png_path, combined, width = 20, height = 22, dpi = 300, bg = "white")
ggsave(svg_path, combined, width = 20, height = 22, bg = "white")
message(sprintf("Saved: %s + %s", png_path, svg_path))

message("\n=== Step 10 complete ===")
