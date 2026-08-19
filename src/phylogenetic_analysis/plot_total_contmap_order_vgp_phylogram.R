#!/usr/bin/env Rscript
# contMap phylogram for TECHINNO2025_ResEff (straight tree).
#
# Analogue of plot_total_contmap_order_vgp.R, but:
#  - type = "phylogram" (rightwards straight tree);
#  - tip names are shown;
#  - VGP logo sits next to the species name (for the 50/80-species overlap set);
#  - ResearchEffort is annotated in parentheses after the name;
#  - tips with ResearchEffort < 15 are dropped;
#  - Order is a coloured vertical strip to the right of the tip labels
#    (one continuous bar per contiguous Order block), with Order names
#    beside the strip - same muted palette as the arc plot.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(png)
})

library(tidyverse)

PHENOTYPE <- "TECHINNO2025_ResEff"
RE_MIN    <- 15

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ROOT_MAIN <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/cognitive_evolution/"
TREE_FILE_MAIN   <- "tree_clootl_overlap_25.nwk"
TRAITS_FILE_MAIN <- "inno_db_25_clootl_overlap.csv"

# Project root: output figure + VGP logo.
# VGP tip set lives with the main cognitive_evolution tables
# (tree_vgp_25_vgp.nwk / inno_db_25_overlap_vgp.csv).
ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/src/Phylogenetic_tests/"
TREE_FILE_B   <- "tree_vgp_25_vgp.nwk"
TRAITS_FILE_B <- "inno_db_25_overlap_vgp.csv"
VGP_LOGO      <- file.path(ROOT, "data", "VGP_tag.png")

BASE <- file.path(ROOT, "results", PHENOTYPE)
OUT  <- file.path(BASE, "figures", "contMap_order_vgp_phylogram.pdf")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Main tree + trait + Order + ResearchEffort
# ---------------------------------------------------------------------------

tree <- read.tree(file.path(ROOT_MAIN, "data", TREE_FILE_MAIN))
xdf <- read.csv(file.path(ROOT_MAIN, "data", TRAITS_FILE_MAIN), stringsAsFactors = FALSE) %>%
  select(Species, all_of(PHENOTYPE), Order, ResearchEffort) %>%
  filter(!is.na(ResearchEffort), ResearchEffort >= RE_MIN)

tree <- keep.tip(tree, intersect(tree$tip.label, xdf$Species))
xdf  <- xdf %>% filter(Species %in% tree$tip.label)

x   <- setNames(xdf[[PHENOTYPE]], xdf$Species)[tree$tip.label]
re  <- setNames(xdf$ResearchEffort, xdf$Species)[tree$tip.label]
ord <- setNames(xdf$Order, xdf$Species)[tree$tip.label]
ord[is.na(ord) | ord == ""] <- "Unknown"

normalize_order <- function(s) {
  s <- tolower(trimws(s))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}
ord <- vapply(ord, normalize_order, character(1))

message(
  "Tips after ResearchEffort >= ", RE_MIN, ": ", Ntip(tree),
  " (dropped from full traits table before keep.tip filter)"
)

# ---------------------------------------------------------------------------
# 2. Overlap set -> which tips get the VGP marker
# ---------------------------------------------------------------------------

tree_b <- read.tree(file.path(ROOT_MAIN, "data", TREE_FILE_B))
xdf_b  <- read.csv(file.path(ROOT_MAIN, "data", TRAITS_FILE_B), stringsAsFactors = FALSE) %>%
  select(Species)
tree_b <- keep.tip(tree_b, intersect(tree_b$tip.label, xdf_b$Species))

is_vgp <- tree$tip.label %in% tree_b$tip.label
message("VGP-marked tips: ", sum(is_vgp))

vgp_logo <- if (file.exists(VGP_LOGO)) png::readPNG(VGP_LOGO) else NULL
if (is.null(vgp_logo)) {
  message("VGP logo not found at ", VGP_LOGO, " - VGP markers will be skipped.")
}

# ---------------------------------------------------------------------------
# 3. Ancestral-state reconstruction
# ---------------------------------------------------------------------------

cm <- contMap(tree, x, plot = FALSE, lut = c("blue", "white", "red"))

# Blank tip labels so we can draw custom name + VGP + (ResearchEffort).
cm$tree$tip.label <- rep("", Ntip(cm$tree))

# ---------------------------------------------------------------------------
# 4. Order colours (same muted HCL palette as the arc plot)
# ---------------------------------------------------------------------------

order_levels <- sort(unique(ord))
n_ord_levels <- length(order_levels)
muted_hues <- seq(15, 375, length.out = n_ord_levels + 1)[1:n_ord_levels]
order_cols <- setNames(
  grDevices::hcl(h = muted_hues, c = 55, l = 58),
  order_levels
)

# ---------------------------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------------------------

ntip    <- Ntip(tree)
tree_h  <- max(phytools::nodeHeights(tree))
cex_lab <- 0.42
cex_re  <- 0.36
cex_ord <- 0.50

sp_labs <- gsub("_", " ", tree$tip.label)
re_labs <- paste0("(", re, ")")

# Tall page so ~1200 tip labels stay readable (same density as rate-shift plots).
FIG_W <- 16
FIG_H <- max(12, 0.13 * ntip + 4)

# Reserve horizontal room: tip labels + optional VGP icon + (RE) + Order strip.
measure_right_edge <- function(fig_w, fig_h, x_right_guess) {
  pdf(NULL, width = fig_w, height = fig_h)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(3.2, 0.4, 0.6, 0.4), xpd = NA)
  plot(
    cm,
    type    = "phylogram",
    fsize   = 0,
    outline = FALSE,
    lwd     = 1.1,
    legend  = FALSE,
    xlim    = c(0, x_right_guess),
    ylim    = c(-max(6, ntip * 0.02), ntip + 0.5)
  )
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_x <- pp$xx[seq_len(ntip)]
  tip_y <- pp$yy[seq_len(ntip)]

  lab_w <- max(strwidth(sp_labs, cex = cex_lab, font = 3))
  re_w  <- max(strwidth(re_labs, cex = cex_re, font = 1))
  icon_w <- if (!is.null(vgp_logo)) {
    # ~height of one tip row, scaled by logo aspect
    tip_row <- diff(range(tip_y)) / max(ntip - 1, 1)
    tip_row * 0.85 * (dim(vgp_logo)[2] / dim(vgp_logo)[1])
  } else {
    0
  }
  gap1 <- strwidth("  ", cex = cex_lab)
  gap2 <- strwidth(" ", cex = cex_lab)
  gap3 <- strwidth("   ", cex = cex_lab)

  # After the farthest tip: name, optional VGP, (RE), then Order strip + label.
  x_label   <- max(tip_x) + strwidth(" ", cex = cex_lab) * 2
  x_after_name <- x_label + lab_w
  x_after_vgp  <- x_after_name + gap1 + icon_w
  x_after_re   <- x_after_vgp + gap2 + re_w
  x_order      <- x_after_re + gap3
  order_bar_w  <- strwidth("W", cex = cex_ord) * 0.7
  order_lab_w  <- max(strwidth(order_levels, cex = cex_ord, font = 2))
  x_right      <- x_order + order_bar_w + gap1 + order_lab_w + strwidth("  ", cex = cex_ord)

  list(
    tip_x = tip_x, tip_y = tip_y,
    x_label = x_label,
    lab_w = lab_w, icon_w = icon_w, re_w = re_w,
    gap1 = gap1, gap2 = gap2, gap3 = gap3,
    x_order = x_order, order_bar_w = order_bar_w,
    x_right = x_right
  )
}

# Two-pass measure so usr↔inches conversion matches the final xlim.
lay <- measure_right_edge(FIG_W, FIG_H, tree_h * 2.5)
lay <- measure_right_edge(FIG_W, FIG_H, lay$x_right * 1.02)
x_lim <- c(-0.02 * lay$x_right, lay$x_right * 1.01)

pdf(OUT, width = FIG_W, height = FIG_H)
par(mar = c(3.2, 0.4, 0.6, 0.4), xpd = NA)

ylim_use <- c(-max(6, ntip * 0.02), ntip + 0.5)
plot(
  cm,
  type    = "phylogram",
  fsize   = 0,
  outline = FALSE,
  lwd     = 1.1,
  legend  = FALSE,
  xlim    = x_lim,
  ylim    = ylim_use
)

pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_x <- pp$xx[seq_len(ntip)]
tip_y <- pp$yy[seq_len(ntip)]

# Recompute placement in the final device coordinates.
lab_w  <- max(strwidth(sp_labs, cex = cex_lab, font = 3))
re_w   <- max(strwidth(re_labs, cex = cex_re, font = 1))
gap1   <- strwidth("  ", cex = cex_lab)
gap2   <- strwidth(" ", cex = cex_lab)
gap3   <- strwidth("   ", cex = cex_lab)
icon_w <- if (!is.null(vgp_logo)) {
  tip_row <- diff(range(tip_y)) / max(ntip - 1, 1)
  tip_row * 0.85 * (dim(vgp_logo)[2] / dim(vgp_logo)[1])
} else {
  0
}
icon_h <- if (!is.null(vgp_logo)) {
  tip_row <- diff(range(tip_y)) / max(ntip - 1, 1)
  tip_row * 0.85
} else {
  0
}

x_label <- max(tip_x) + strwidth(" ", cex = cex_lab) * 2
x_order <- x_label + lab_w + gap1 + icon_w + gap2 + re_w + gap3
order_bar_w <- strwidth("W", cex = cex_ord) * 0.7

## -- tip labels: Species  [VGP]  (ResearchEffort) -------------------------
for (i in seq_len(ntip)) {
  text(
    x_label, tip_y[i],
    labels = sp_labs[i],
    adj = c(0, 0.5),
    cex = cex_lab,
    font = 3,
    col = "grey15"
  )

  x_cursor <- x_label + strwidth(sp_labs[i], cex = cex_lab, font = 3)

  if (!is.null(vgp_logo) && is_vgp[i]) {
    cx <- x_cursor + gap1 + icon_w / 2
    cy <- tip_y[i]
    ndc_x <- grconvertX(cx, from = "user", to = "ndc")
    ndc_y <- grconvertY(cy, from = "user", to = "ndc")
    ndc_w <- abs(
      grconvertX(cx + icon_w / 2, from = "user", to = "ndc") -
        grconvertX(cx - icon_w / 2, from = "user", to = "ndc")
    )
    ndc_h <- abs(
      grconvertY(cy + icon_h / 2, from = "user", to = "ndc") -
        grconvertY(cy - icon_h / 2, from = "user", to = "ndc")
    )
    grid::pushViewport(grid::viewport(
      x = grid::unit(ndc_x, "npc"), y = grid::unit(ndc_y, "npc"),
      width = grid::unit(ndc_w, "npc"), height = grid::unit(ndc_h, "npc")
    ))
    grid::grid.raster(vgp_logo, interpolate = TRUE)
    grid::popViewport()
  }

  # Align all (ResearchEffort) annotations to the same column so they read
  # as a quiet secondary field rather than sticking to variable name lengths.
  x_re <- x_label + lab_w + gap1 + icon_w + gap2
  text(
    x_re, tip_y[i],
    labels = re_labs[i],
    adj = c(0, 0.5),
    cex = cex_re,
    font = 1,
    col = "grey45"
  )
}

## -- Order strip: one coloured bar per contiguous Order block -------------
ord_runs  <- rle(as.character(ord))
run_end   <- cumsum(ord_runs$lengths)
run_start <- run_end - ord_runs$lengths + 1
half_gap  <- 0.35

for (i in seq_along(ord_runs$lengths)) {
  ii <- run_start[i]:run_end[i]
  y1 <- min(tip_y[ii]) - half_gap
  y2 <- max(tip_y[ii]) + half_gap
  ymid <- (y1 + y2) / 2
  ord_name <- ord_runs$values[i]
  ord_col  <- order_cols[[ord_name]]

  rect(
    x_order, y1,
    x_order + order_bar_w, y2,
    col = ord_col, border = NA
  )
  text(
    x_order + order_bar_w + gap1, ymid,
    labels = ord_name,
    adj = c(0, 0.5),
    cex = cex_ord,
    font = 2,
    col = ord_col
  )
}

## -- Innovation-rate colour bar (bottom-left, inside the tree area) --------
lims <- cm$lims
grad_cols <- unname(cm$cols)
n_grad <- length(grad_cols)

bar_x0 <- 0
bar_x1 <- tree_h * 0.45
bar_h  <- abs(ylim_use[1]) * 0.22
bar_y0 <- ylim_use[1] * 0.55
bar_y1 <- bar_y0 + bar_h

x_edges <- seq(bar_x0, bar_x1, length.out = n_grad + 1)
for (i in seq_len(n_grad)) {
  rect(x_edges[i], bar_y0, x_edges[i + 1], bar_y1, col = grad_cols[i], border = NA)
}

lab_vals <- pretty(lims, n = 4)
lab_vals <- lab_vals[lab_vals >= lims[1] & lab_vals <= lims[2]]
tick_h <- bar_h * 0.45
for (v in lab_vals) {
  xx <- bar_x0 + (v - lims[1]) / diff(lims) * (bar_x1 - bar_x0)
  segments(xx, bar_y0 - tick_h, xx, bar_y0, lwd = 0.8)
  text(xx, bar_y0 - tick_h * 1.15, signif(v, 2),
       adj = c(0.5, 1), cex = 0.55, col = "grey25")
}
text(
  (bar_x0 + bar_x1) / 2,
  bar_y0 - tick_h * 2.4,
  "Innovation rate",
  adj = c(0.5, 1), font = 2, cex = 0.7, col = "grey20"
)

dev.off()
message("Saved: ", OUT)
message("Figure size: ", FIG_W, " x ", round(FIG_H, 1), " in; tips: ", ntip)

