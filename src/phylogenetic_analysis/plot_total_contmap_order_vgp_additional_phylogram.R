#!/usr/bin/env Rscript
# contMap phylograms for TECHINNO2025_ResEff — two figures from one script:
#   1) VGP ∪ all Additional_genomes_new.tsv species
#   2) VGP ∪ VGP (new) only
#
# VGP logo: species in vgp_birds.csv.
# VGP* logo: Genome == "VGP (new)" in Additional_genomes_new.tsv.
# ResearchEffort >= 15 still applies.

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
ADDITIONAL_FILE  <- "Additional_genomes_new.tsv"
VGP_BIRDS_FILE   <- "vgp_birds.csv"

ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/src/Phylogenetic_tests/"
VGP_LOGO <- file.path(ROOT, "data", "VGP_tag.png")

BASE <- file.path(ROOT, "results", PHENOTYPE)
OUT_ADDITIONAL <- file.path(BASE, "figures", "contMap_order_vgp_additional_phylogram_.pdf")
OUT_VGP_ONLY   <- file.path(BASE, "figures", "contMap_order_vgp_only_phylogram.pdf")
dir.create(dirname(OUT_ADDITIONAL), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

normalize_order <- function(s) {
  s <- tolower(trimws(s))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

normalize_species <- function(x) {
  x <- gsub(" ", "_", trimws(as.character(x)))
  sub("_+$", "", x)
}

binomial_species <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(parts, function(p) {
    if (!length(p)) return("")
    paste(p[seq_len(min(2L, length(p)))], collapse = "_")
  }, character(1))
}

read_additional_genomes <- function(path) {
  raw <- readLines(path, warn = FALSE)
  if (!length(raw)) {
    return(data.frame(
      Species = character(), Genome = character(), stringsAsFactors = FALSE
    ))
  }
  rows <- strsplit(raw[-1], "\t", fixed = TRUE)
  sci <- trimws(vapply(rows, `[`, character(1), 1))
  gen <- trimws(vapply(rows, function(r) {
    if (length(r) >= 2) r[[2]] else ""
  }, character(1)))
  sci <- sub("[[:space:]]*\\(.*\\)[[:space:]]*$", "", sci)
  data.frame(
    Species = normalize_species(sci),
    Genome  = gen,
    stringsAsFactors = FALSE
  )
}

is_vgp_new_label <- function(genome) {
  g <- tolower(gsub("[[:space:]]+", " ", trimws(as.character(genome))))
  g == "vgp (new)"
}

read_vgp_birds <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Species" %in% names(d)) {
    stop("vgp_birds.csv must contain a Species column")
  }
  sp <- unique(normalize_species(d$Species))
  sp <- sp[nzchar(sp)]
  unique(c(sp, binomial_species(sp)))
}

# ---------------------------------------------------------------------------
# Shared data
# ---------------------------------------------------------------------------

vgp_species <- read_vgp_birds(file.path(ROOT_MAIN, "data", VGP_BIRDS_FILE))
add_df <- read_additional_genomes(file.path(ROOT_MAIN, "data", ADDITIONAL_FILE))
add_species <- unique(add_df$Species[nzchar(add_df$Species)])
vgp_new_species <- unique(add_df$Species[is_vgp_new_label(add_df$Genome) & nzchar(add_df$Species)])

tree_full <- read.tree(file.path(ROOT_MAIN, "data", TREE_FILE_MAIN))
xdf_full <- read.csv(file.path(ROOT_MAIN, "data", TRAITS_FILE_MAIN), stringsAsFactors = FALSE) %>%
  select(Species, all_of(PHENOTYPE), Order, ResearchEffort) %>%
  filter(!is.na(ResearchEffort), ResearchEffort >= RE_MIN)

vgp_logo <- if (file.exists(VGP_LOGO)) png::readPNG(VGP_LOGO) else NULL
if (is.null(vgp_logo)) {
  message("VGP logo not found at ", VGP_LOGO, " - VGP markers will be skipped.")
}

prepare_subset <- function(keep_species, label) {
  missing_on_tree <- setdiff(keep_species, tree_full$tip.label)
  if (length(missing_on_tree)) {
    message(
      "[", label, "] Warning: ", length(missing_on_tree),
      " keep-set species absent from main tree (dropped)"
    )
  }

  xdf <- xdf_full %>% filter(Species %in% keep_species)
  tree <- keep.tip(tree_full, intersect(tree_full$tip.label, xdf$Species))
  xdf  <- xdf %>% filter(Species %in% tree$tip.label)

  x   <- setNames(xdf[[PHENOTYPE]], xdf$Species)[tree$tip.label]
  re  <- setNames(xdf$ResearchEffort, xdf$Species)[tree$tip.label]
  ord <- setNames(xdf$Order, xdf$Species)[tree$tip.label]
  ord[is.na(ord) | ord == ""] <- "Unknown"
  ord <- vapply(ord, normalize_order, character(1))

  is_vgp <- tree$tip.label %in% vgp_species
  is_vgp_new <- tree$tip.label %in% vgp_new_species
  is_vgp_new[is_vgp] <- FALSE

  message(
    "[", label, "] Tips plotted: ", Ntip(tree),
    " (VGP: ", sum(is_vgp),
    "; VGP (new)*: ", sum(is_vgp_new),
    "; other: ", sum(!is_vgp & !is_vgp_new), ")"
  )
  if (sum(is_vgp_new)) {
    message("[", label, "] VGP (new): ", paste(tree$tip.label[is_vgp_new], collapse = ", "))
  }

  list(
    tree = tree, x = x, re = re, ord = ord,
    is_vgp = is_vgp, is_vgp_new = is_vgp_new
  )
}

plot_contmap_phylogram <- function(inp, out_pdf) {
  tree <- inp$tree
  x <- inp$x
  re <- inp$re
  ord <- inp$ord
  is_vgp <- inp$is_vgp
  is_vgp_new <- inp$is_vgp_new

  cm <- contMap(tree, x, plot = FALSE, lut = c("blue", "white", "red"))
  cm$tree$tip.label <- rep("", Ntip(cm$tree))

  order_levels <- sort(unique(ord))
  n_ord_levels <- length(order_levels)
  muted_hues <- seq(15, 375, length.out = n_ord_levels + 1)[1:n_ord_levels]
  order_cols <- setNames(
    grDevices::hcl(h = muted_hues, c = 55, l = 58),
    order_levels
  )

  ntip    <- Ntip(tree)
  tree_h  <- max(phytools::nodeHeights(tree))
  cex_lab <- 0.70
  cex_re  <- 0.58
  cex_ord <- 0.72

  sp_labs <- gsub("_", " ", tree$tip.label)
  re_labs <- paste0("(", re, ")")

  FIG_W <- 14
  FIG_H <- max(10, 0.18 * ntip + 3.5)

  measure_right_edge <- function(fig_w, fig_h, x_right_guess) {
    pdf(NULL, width = fig_w, height = fig_h)
    on.exit(dev.off(), add = TRUE)
    par(mar = c(3.2, 0.4, 0.6, 0.4), xpd = NA)
    plot(
      cm,
      type    = "phylogram",
      fsize   = 0,
      outline = FALSE,
      lwd     = 1.6,
      legend  = FALSE,
      xlim    = c(0, x_right_guess),
      ylim    = c(-max(4, ntip * 0.04), ntip + 0.5)
    )
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    tip_x <- pp$xx[seq_len(ntip)]
    tip_y <- pp$yy[seq_len(ntip)]

    lab_w <- max(strwidth(sp_labs, cex = cex_lab, font = 3))
    re_w  <- max(strwidth(re_labs, cex = cex_re, font = 1))
    icon_w <- if (!is.null(vgp_logo)) {
      tip_row <- diff(range(tip_y)) / max(ntip - 1, 1)
      tip_row * 0.85 * (dim(vgp_logo)[2] / dim(vgp_logo)[1])
    } else {
      0
    }
    star_w <- strwidth("*", cex = cex_lab, font = 2)
    gap1 <- strwidth("  ", cex = cex_lab)
    gap2 <- strwidth(" ", cex = cex_lab)
    gap3 <- strwidth("   ", cex = cex_lab)

    x_label  <- max(tip_x) + strwidth(" ", cex = cex_lab) * 2
    x_order  <- x_label + lab_w + gap1 + icon_w + star_w + gap2 + re_w + gap3
    order_bar_w <- strwidth("W", cex = cex_ord) * 0.7
    order_lab_w <- max(strwidth(order_levels, cex = cex_ord, font = 2))
    x_right  <- x_order + order_bar_w + gap1 + order_lab_w + strwidth("  ", cex = cex_ord)

    list(x_right = x_right)
  }

  lay <- measure_right_edge(FIG_W, FIG_H, tree_h * 2.5)
  lay <- measure_right_edge(FIG_W, FIG_H, lay$x_right * 1.02)
  x_lim <- c(-0.02 * lay$x_right, lay$x_right * 1.01)

  pdf(out_pdf, width = FIG_W, height = FIG_H)
  par(mar = c(3.2, 0.4, 0.6, 0.4), xpd = NA)

  ylim_use <- c(-max(4, ntip * 0.04), ntip + 0.5)
  plot(
    cm,
    type    = "phylogram",
    fsize   = 0,
    outline = FALSE,
    lwd     = 1.6,
    legend  = FALSE,
    xlim    = x_lim,
    ylim    = ylim_use
  )

  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_x <- pp$xx[seq_len(ntip)]
  tip_y <- pp$yy[seq_len(ntip)]

  lab_w  <- max(strwidth(sp_labs, cex = cex_lab, font = 3))
  re_w   <- max(strwidth(re_labs, cex = cex_re, font = 1))
  gap1   <- strwidth("  ", cex = cex_lab)
  gap2   <- strwidth(" ", cex = cex_lab)
  gap3   <- strwidth("   ", cex = cex_lab)
  star_w <- strwidth("*", cex = cex_lab, font = 2)
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
  x_order <- x_label + lab_w + gap1 + icon_w + star_w + gap2 + re_w + gap3
  order_bar_w <- strwidth("W", cex = cex_ord) * 0.7

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
    show_badge <- !is.null(vgp_logo) && (is_vgp[i] || is_vgp_new[i])

    if (show_badge) {
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

      if (is_vgp_new[i]) {
        text(
          x_cursor + gap1 + icon_w, tip_y[i],
          labels = "*",
          adj = c(0, 0.35),
          cex = cex_lab * 1.15,
          font = 2,
          col = "grey20"
        )
      }
    }

    x_re <- x_label + lab_w + gap1 + icon_w + star_w + gap2
    text(
      x_re, tip_y[i],
      labels = re_labs[i],
      adj = c(0, 0.5),
      cex = cex_re,
      font = 1,
      col = "grey45"
    )
  }

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
  message("Saved: ", out_pdf, " (", FIG_W, " x ", round(FIG_H, 1), " in; tips: ", ntip, ")")
}

# ---------------------------------------------------------------------------
# Figure 1: VGP ∪ all Additional_genomes_new
# Figure 2: VGP ∪ VGP (new) only
# ---------------------------------------------------------------------------

keep_additional <- union(vgp_species, add_species)
keep_vgp_only   <- union(vgp_species, vgp_new_species)

message(
  "Keep-sets — additional: ", length(keep_additional),
  "; VGP-only: ", length(keep_vgp_only),
  " (vgp_birds: ", length(vgp_species),
  "; Additional rows: ", length(add_species),
  "; VGP (new): ", length(vgp_new_species), ")"
)

inp_additional <- prepare_subset(keep_additional, "additional")
inp_vgp_only   <- prepare_subset(keep_vgp_only, "vgp_only")

plot_contmap_phylogram(inp_additional, OUT_ADDITIONAL)
plot_contmap_phylogram(inp_vgp_only, OUT_VGP_ONLY)
