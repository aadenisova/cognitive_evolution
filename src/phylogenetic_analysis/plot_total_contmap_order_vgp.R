#!/usr/bin/env Rscript
# contMap-style arc plot for TOTALINNOVATIONS2025_ResEff.
#
# Differences from plot_total_contmap_big.R:
#  - tip names are NOT plotted (tip.label is blanked out before plotting);
#  - no plot title;
#  - xlim/ylim are set explicitly (from the tree's own depth) so the ring,
#    Order labels and callouts all stay within the page instead of spilling
#    past its edges;
#  - Order is shown as one coloured ring around the tree, one continuous
#    wedge per Order with a small gap between neighbouring Orders. All
#    Order labels use one fixed font size and are coloured to match their
#    own wedge: wide wedges get the name curved along the arc right above
#    them, narrow wedges instead get a short pointer line with the name
#    written off to the side - there is no separate Order legend panel;
#  - the built-in contMap colour legend is removed and redrawn as a thin,
#    borderless "Innovation rate" bar tucked inside the arc, near the
#    baseline (no separate legend panel);
#  - species that are also present in the smaller 50-species overlap set
#    (binomial_b) are marked with a small VGP logo right at the tip of the
#    branch, tilted the way a radiating tip label would be (half the tip's
#    own angle, so it never stands fully upright).

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(png)
})

library(tidyverse)

PHENOTYPE <- "TECHINNO2025_ResEff"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# Main (100-species) tree + trait table, same source as plot_total_contmap_big.R
ROOT_MAIN <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/cognitive_evolution/"
TREE_FILE_MAIN   <- "tree_clootl_overlap_25.nwk"
TRAITS_FILE_MAIN <- "inno_db_25_clootl_overlap.csv"

# Project root: used for the output figure and for the 50-species overlap
# set that defines which species get the VGP marker.
ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/src/Phylogenetic_tests/"
TREE_FILE_B   <- "tree_vgp_0_vgp.nwk"
TRAITS_FILE_B <- "inno_db_0_overlap_vgp.csv"
VGP_LOGO      <- file.path(ROOT, "data", "VGP_tag.png")

BASE <- file.path(ROOT, "results", PHENOTYPE)
OUT  <- file.path(BASE, "figures", "contMap_order_vgp.pdf")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Main tree + trait ("Innovation rate") + Order per tip
# ---------------------------------------------------------------------------

tree <- read.tree(file.path(ROOT_MAIN, "data", TREE_FILE_MAIN))
xdf <- read.csv(file.path(ROOT_MAIN, "data", TRAITS_FILE_MAIN), stringsAsFactors = FALSE) %>%
  select(Species, all_of(PHENOTYPE), Order)

tree <- keep.tip(tree, intersect(tree$tip.label, xdf$Species))
xdf  <- xdf %>% filter(Species %in% tree$tip.label)

x   <- setNames(xdf[[PHENOTYPE]], xdf$Species)[tree$tip.label]
ord <- setNames(xdf$Order, xdf$Species)[tree$tip.label]
ord[is.na(ord) | ord == ""] <- "Unknown"

# Normalise casing: "Order" -> first letter upper case, rest lower case.
normalize_order <- function(s) {
  s <- tolower(trimws(s))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}
ord <- vapply(ord, normalize_order, character(1))

# ---------------------------------------------------------------------------
# 2. 50-species overlap set -> which tips get the VGP marker
# ---------------------------------------------------------------------------

tree_b <- read.tree(file.path(ROOT, "data", TREE_FILE_B))
xdf_b  <- read.csv(file.path(ROOT, "data", TRAITS_FILE_B), stringsAsFactors = FALSE) %>%
  select(Species)
tree_b <- keep.tip(tree_b, intersect(tree_b$tip.label, xdf_b$Species))

# binomial_b: the species names of the 50-species overlap tree, used only to
# decide which tips in the main tree get the VGP marker below.
binomial_b <- gsub("_", " ", tree_b$tip.label)
is_vgp <- tree$tip.label %in% tree_b$tip.label

vgp_logo <- if (file.exists(VGP_LOGO)) png::readPNG(VGP_LOGO) else NULL
if (is.null(vgp_logo)) {
  message("VGP logo not found at ", VGP_LOGO, " - VGP markers will be skipped.")
}

# ---------------------------------------------------------------------------
# 3. Ancestral-state reconstruction for the colour map
# ---------------------------------------------------------------------------

cm <- contMap(tree, x, plot = FALSE, lut = c("blue", "white", "red"))

# Blank out the tip labels on the tree actually used for plotting. This is
# more reliable than relying on fsize/ftype alone (which some phytools plot
# methods can ignore for certain tree types) - an empty tip.label guarantees
# no species name ever gets drawn.
cm$tree$tip.label <- rep("", Ntip(cm$tree))

# ---------------------------------------------------------------------------
# 4. Order colours - muted/desaturated qualitative palette (low chroma,
#    medium lightness) so the ring doesn't compete visually with the
#    Innovation-rate colour map.
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
#
# Single panel: tree, Order ring + labels, VGP markers, and the
# "Innovation rate" colour bar tucked inside the arc near the baseline -
# there is no separate legend panel any more.
# ---------------------------------------------------------------------------

# Tree depth (max root-to-tip distance) computed *before* plotting, so we can
# tell plot() to reserve enough room (via xlim/ylim) for the ring and the
# Order labels/callouts that sit outside the tips - this is what keeps the
# whole figure from spilling past the page edges.
tree_depth   <- max(phytools::nodeHeights(tree))
r_max_factor <- 4.1
r_lim <- tree_depth * r_max_factor

# Overall figure scale: this makes the whole page bigger (same 15:9 aspect
# ratio, so nothing about the layout above changes) and scales up line
# widths / label sizes by the same factor, so the figure gets bigger
# without changing any relative proportions - just a uniform zoom. Raise it
# further for an even bigger output.
FIG_SCALE <- 1

pdf(OUT, width = 15 * FIG_SCALE, height = 9 * FIG_SCALE)

## -- tree - no tip labels, no built-in legend, no title -------------------
par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA)
plot(
  cm,
  type    = "arc",
  fsize   = 0,
  outline = FALSE,
  lwd     = 1.3 * FIG_SCALE,
  legend  = FALSE,
  xlim    = c(-r_lim, r_lim),
  ylim    = c(-tree_depth * 0.02, r_lim)
)

pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
ntip <- Ntip(tree)
tip_x <- pp$xx[1:ntip]
tip_y <- pp$yy[1:ntip]
tip_r <- sqrt(tip_x^2 + tip_y^2)
tip_dir_x <- tip_x / tip_r
tip_dir_y <- tip_y / tip_r
tip_ang <- atan2(tip_y, tip_x)
max_r <- max(tip_r)

## -- Order ring: one continuous coloured band per Order, with a small gap
##    left between neighbouring Orders --------------------------------------
draw_order_wedge <- function(a0, a1, r0, r1, col) {
  if (a1 < a0) { tmp <- a0; a0 <- a1; a1 <- tmp }
  n_seg <- max(2, ceiling((a1 - a0) / (pi / 180)))
  a_seq <- seq(a0, a1, length.out = n_seg + 1)
  xs <- c(r0 * cos(a_seq), rev(r1 * cos(a_seq)))
  ys <- c(r0 * sin(a_seq), rev(r1 * sin(a_seq)))
  polygon(xs, ys, col = col, border = col)
}

## -- All Order labels use one fixed font size, are coloured to match their
##    own wedge, and are drawn the same way: a short pointer line from the
##    ring out to a label written to the side. There is no curved,
##    arc-following label any more - a fixed pointer length is far easier to
##    keep from overlapping a neighbour than text curved along the ring,
##    which is what caused labels for adjacent narrow Orders to run into
##    each other. Each callout picks whichever of three radii was used
##    *least recently in angular terms* (the level whose last label is
##    angularly farthest from this one) rather than a simple 1-2-3 round
##    robin - a fixed cycle can still put two close-together Orders on the
##    same radius whenever an unrelated Order happens to sit between them
##    and shifts the count; picking by angular distance instead means
##    several narrow, closely spaced Orders (e.g. Podicipediformes,
##    Otidiformes, Cuculiformes) reliably end up on three different radii
##    regardless of what else is interspersed between them. -----------------
cex_label <- 0.65 * FIG_SCALE

callout_radii <- NULL       # set below, once max_r is known
callout_last_angle <- NULL  # set below: one "last used at angle" per radius

## The automatic angular-distance placement above still left Cuculiformes
## and Podicipediformes too close together, so these two get an explicit
## manual radius instead of the automatic pick: Cuculiformes pulled in
## closer to the ring, Podicipediformes pushed further out.
callout_radius_override <- NULL  # set below, once max_r is known

draw_order_callout <- function(a_mid, r_from, label, cex, col) {
  if (!is.null(callout_radius_override) && label %in% names(callout_radius_override)) {
    r_to <- callout_radius_override[[label]]
  } else {
    gaps <- abs(a_mid - callout_last_angle)
    lvl  <- which.max(gaps)
    callout_last_angle[lvl] <<- a_mid
    r_to <- callout_radii[lvl]
  }

  x0 <- r_from * cos(a_mid); y0 <- r_from * sin(a_mid)
  x1 <- r_to   * cos(a_mid); y1 <- r_to   * sin(a_mid)
  segments(x0, y0, x1, y1, col = col, lwd = FIG_SCALE)

  side <- if (x1 >= 0) 1 else -1
  text(x1 + side * tree_depth * 0.025, y1, label,
       adj = c(if (side > 0) 0 else 1, 0.5), cex = cex, col = col)
}

ring_r1 <- max_r * 1.04
ring_r2 <- max_r * 1.10
callout_radii <- max_r * c(1.20, 1.28, 1.36)
callout_last_angle <- rep(-Inf, length(callout_radii))
callout_radius_override <- c(
  "Cuculiformes"     = max_r * 1.12,  # pulled in closer to the ring
  "Otidiformes" = max_r * 1.22,   # pushed further out
  "Podicipediformes" = max_r * 1.35   # pushed further out
)

ang_sorted <- sort(tip_ang)
half_step <- stats::median(diff(ang_sorted)) / 2
gap_frac  <- 1                        # fraction of half_step left as a gap
wedge_ext <- half_step * (1 - gap_frac)

ord_runs  <- rle(as.character(ord))
run_end   <- cumsum(ord_runs$lengths)
run_start <- run_end - ord_runs$lengths + 1

for (i in seq_along(ord_runs$lengths)) {
  ii <- run_start[i]:run_end[i]
  a_lo <- min(tip_ang[ii]) - wedge_ext
  a_hi <- max(tip_ang[ii]) + wedge_ext
  a_mid <- (a_lo + a_hi) / 2
  ord_name <- ord_runs$values[i]
  ord_col  <- order_cols[[ord_name]]

  draw_order_wedge(a_lo, a_hi, ring_r1, ring_r2, ord_col)
  draw_order_callout(a_mid, ring_r2, ord_name, cex_label, ord_col)
}

## -- VGP marker: small logo right at the tip of the branch, for species
##    also present in the 50-species overlap set (binomial_b). Rotated
##    exactly the way a radiating tip label (species name) would be: angled
##    along the tip's own radial direction, and flipped by 180 degrees
##    whenever that direction points into the left half of the plot - the
##    same rule ape/phytools use for fan/arc tip labels so text (or, here,
##    the logo) never ends up upside-down. rasterImage() can't rotate, so
##    the icon is drawn with grid::grid.raster() instead, positioned via
##    grconvertX/Y (base R) so it lines up with the same tip coordinates.
##    grid.raster() itself has no rotation argument, so the rotation is done
##    by pushing a rotated grid::viewport() (which does have "angle") sized
##    and centred on the icon, then drawing the raster to fill it. ---------
if (!is.null(vgp_logo)) {
  logo_asp <- dim(vgp_logo)[2] / dim(vgp_logo)[1]
  icon_h <- max_r * 0.0110
  icon_w <- icon_h * logo_asp
  icon_r <- max_r * 1.02

  for (t in which(is_vgp)) {
    cx <- tip_dir_x[t] * icon_r
    cy <- tip_dir_y[t] * icon_r

    # Same orientation rule as a radiating tip label: use the tip's own
    # angle, but flip 180 degrees whenever that would point left (into the
    # (90, 270) range) so the logo is never rotated upside-down.
    ang_deg <- (tip_ang[t] * 180 / pi) %% 360
    rot_deg <- if (ang_deg > 90 && ang_deg < 270) ang_deg - 180 else ang_deg

    ndc_x  <- grconvertX(cx, from = "user", to = "ndc")
    ndc_y  <- grconvertY(cy, from = "user", to = "ndc")
    ndc_w  <- abs(grconvertX(cx + icon_w / 2, from = "user", to = "ndc") -
                  grconvertX(cx - icon_w / 2, from = "user", to = "ndc"))
    ndc_h  <- abs(grconvertY(cy + icon_h / 2, from = "user", to = "ndc") -
                  grconvertY(cy - icon_h / 2, from = "user", to = "ndc"))

    grid::pushViewport(grid::viewport(
      x = grid::unit(ndc_x, "npc"), y = grid::unit(ndc_y, "npc"),
      width = grid::unit(ndc_w, "npc"), height = grid::unit(ndc_h, "npc"),
      angle = rot_deg
    ))
    grid::grid.raster(vgp_logo, interpolate = TRUE)
    grid::popViewport()
  }
}

## -- "Innovation rate" colour bar - tucked inside the arc, near the
##    baseline, well clear of the tree branches above it -------------------
lims <- cm$lims
grad_cols <- unname(cm$cols)
n_grad <- length(grad_cols)

bar_w  <- max_r * 0.55
bar_h  <- max_r * 0.022
bar_x0 <- -bar_w / 2
bar_x1 <- bar_w / 2
bar_y0 <- max_r * 0.11
bar_y1 <- bar_y0 + bar_h

x_edges <- seq(bar_x0, bar_x1, length.out = n_grad + 1)
for (i in seq_len(n_grad)) {
  rect(x_edges[i], bar_y0, x_edges[i + 1], bar_y1, col = grad_cols[i], border = NA)
}

lab_vals <- pretty(lims, n = 4)
lab_vals <- lab_vals[lab_vals >= lims[1] & lab_vals <= lims[2]]
for (v in lab_vals) {
  xx <- bar_x0 + (v - lims[1]) / diff(lims) * (bar_x1 - bar_x0)
  segments(xx, bar_y0 - max_r * 0.012, xx, bar_y0, lwd = FIG_SCALE)
  text(xx, bar_y0 - max_r * 0.028, signif(v, 2), adj = c(0.5, 1), cex = 0.65 * FIG_SCALE)
}

text(0, bar_y0 - max_r * 0.06, "Innovation rate",
     adj = c(0.5, 1), font = 2, cex = 0.85 * FIG_SCALE)

dev.off()
message("Saved: ", OUT)

