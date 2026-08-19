#!/usr/bin/env Rscript
# RRphylo + continuous ASR branch colouring (straight phylogram)
# =====================================================================
# Branch colours = ancestral phenotype (fastAnc); lollipops use tip colours.
# Tip-threshold + search.shift variant:
#   phylo_rrphylo_tip_threshold.R
#
# Outputs: results/rate_shifts/

suppressPackageStartupMessages({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) == 1) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
    local_lib <- file.path(script_dir, ".Rlib")
  } else {
    local_lib <- file.path(getwd(), ".Rlib")
  }
  if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

  library(ape)
  library(phytools)
  library(RRphylo)
  library(tidyverse)
})

PROJECT_ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/src/Phylogenetic_tests"
DATA_DIR <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/cognitive_evolution/data"
OUT_DIR <- file.path(PROJECT_ROOT, "results", "rate_shifts")

TREE_FILE <- "tree_clootl_overlap_50.nwk"
TRAITS_FILE <- "inno_db_50_clootl_overlap.csv"

#TREE_FILE <- "tree_astralpro.l6.tre_overlap_0_vgp.nwk"
#TRAITS_FILE <- "inno_db_0_astralpro.l6.tre_overlap_vgp.csv"

SPECIES_COL <- "Species"
PHENOTYPE <- "TOTALINNOVATIONS2025_ResEff"
INNO_THRESHOLD <- 0.5

# TEMP: tip labels show ResearchEffort on the right. Set to FALSE to hide.
SHOW_RESEARCH_EFFORT <- TRUE

CLUS <- 2 / parallel::detectCores()

# Continuous ASR palette for whole-tree branch colours (low → mid → high)
ASR_LUT <- c("#56B4E9", "#B0B0B0", "#E69F00")
COL_BG <- "#666666"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- data ----
load_and_match <- function() {
  tree <- read.tree(file.path(DATA_DIR, TREE_FILE))
  tree$node.label <- NULL
  dat <- read.csv(
    file.path(DATA_DIR, TRAITS_FILE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  dat <- dat %>% filter(ResearchEffort >= 15)
  
  if ("Unnamed: 0" %in% names(dat)) {
    dat <- dat[, setdiff(names(dat), "Unnamed: 0"), drop = FALSE]
  }
  rownames(dat) <- dat[[SPECIES_COL]]
  keep <- intersect(tree$tip.label, dat[[SPECIES_COL]])
  tree <- keep.tip(tree, keep)
  dat <- dat[tree$tip.label, , drop = FALSE]
  y <- dat[[PHENOTYPE]]
  if (!is.numeric(y)) y <- as.numeric(y)
  names(y) <- tree$tip.label
  ok <- !is.na(y)
  if (!all(ok)) {
    tree <- keep.tip(tree, names(y)[ok])
    dat <- dat[tree$tip.label, , drop = FALSE]
    y <- y[tree$tip.label]
  }
  if (!is.binary(tree)) tree <- multi2di(tree)
  tree$node.label <- NULL
  list(tree = tree, dat = dat, y = y[tree$tip.label])
}


extract_branch_rates <- function(RR) {
  tree <- RR$tree
  rates <- RR$rates[, 1]
  mean_abs <- mean(abs(rates))
  edge <- tree$edge
  child <- edge[, 2]
  ntip <- Ntip(tree)
  is_tip <- child <= ntip
  lab <- ifelse(is_tip, tree$tip.label[child], as.character(child))
  rate_vals <- as.numeric(rates[lab])
  data.frame(
    parent = edge[, 1],
    node = child,
    branch = lab,
    species = ifelse(is_tip, tree$tip.label[child], NA_character_),
    branch_type = ifelse(is_tip, "terminal", "internal"),
    branch_length = tree$edge.length,
    evolutionary_rate = rate_vals,
    relative_rate = rate_vals / mean_abs,
    abs_rate = abs(rate_vals),
    stringsAsFactors = FALSE
  )
}

# Continuous ASR colouring for the whole tree (same idea as contMap)
branch_asr_colors <- function(tree, y, lut = ASR_LUT, lwd = 2.2) {
  y <- y[tree$tip.label]
  anc <- phytools::fastAnc(tree, y)
  ntip <- Ntip(tree)
  state_at <- function(node) {
    if (node <= ntip) as.numeric(y[tree$tip.label[node]])
    else as.numeric(anc[as.character(node)])
  }
  edge_val <- vapply(seq_len(nrow(tree$edge)), function(i) {
    mean(c(state_at(tree$edge[i, 1]), state_at(tree$edge[i, 2])))
  }, numeric(1))
  rng <- range(c(as.numeric(y), as.numeric(anc)), na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) rng <- rng + c(-1e-6, 1e-6)
  to_col <- function(v) {
    norm <- (v - rng[1]) / diff(rng)
    norm[!is.finite(norm)] <- 0.5
    norm <- pmin(pmax(norm, 0), 1)
    rgb(colorRamp(lut)(norm), maxColorValue = 255)
  }
  tip_cols <- to_col(as.numeric(y))
  names(tip_cols) <- names(y)
  list(
    edge_col = to_col(edge_val),
    edge_lwd = rep(lwd, length(edge_val)),
    tip_cols = tip_cols,
    edge_val = edge_val,
    anc = anc,
    range = rng,
    lut = lut
  )
}

draw_asr_legend <- function(x_left, y_bottom, width, height, rng, lut, title = "Reconstructed phenotype") {
  n <- 80
  xs <- seq(x_left, x_left + width, length.out = n + 1)
  cols <- rgb(colorRamp(lut)(seq(0, 1, length.out = n)), maxColorValue = 255)
  for (i in seq_len(n)) {
    rect(xs[i], y_bottom, xs[i + 1], y_bottom + height, col = cols[i], border = NA)
  }
  rect(x_left, y_bottom, x_left + width, y_bottom + height, border = "grey30", lwd = 0.8)
  ticks <- pretty(rng, n = 4)
  ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
  for (tv in ticks) {
    tx <- x_left + (tv - rng[1]) / diff(rng) * width
    segments(tx, y_bottom, tx, y_bottom - 0.25 * height, col = "grey30")
    text(tx, y_bottom - 0.55 * height, labels = signif(tv, 3),
         cex = 0.65, col = "grey20", adj = c(0.5, 1))
  }
  text(
    x_left + width / 2, y_bottom + height + 0.55 * height,
    labels = title, cex = 0.75, font = 2, col = "grey20"
  )
}

plot_straight <- function(
  tree, y, out_pdf, out_png,
  research_effort = NULL,
  orders = NULL
) {
  y <- y[tree$tip.label]
  if (!is.null(research_effort)) {
    research_effort <- research_effort[tree$tip.label]
  }
  if (!is.null(orders)) {
    orders <- setNames(as.character(orders[tree$tip.label]), tree$tip.label)
  }

  # Whole-tree branch colours from continuous ASR (replaces tip-threshold colours)
  asr <- branch_asr_colors(tree, y, lut = ASR_LUT, lwd = 2.2)
  ec <- list(col = asr$edge_col, lwd = asr$edge_lwd)
  tip_cols <- asr$tip_cols
  axis_label <- "Technical Innovation"

  format_order_label <- function(x) {
    paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
  }

  # Contiguous Order blocks along tip order (as plotted)
  order_blocks <- function(ord_vec) {
    if (is.null(ord_vec) || !length(ord_vec)) {
      return(data.frame(
        order = character(0), start = integer(0), end = integer(0),
        stringsAsFactors = FALSE
      ))
    }
    o <- as.character(ord_vec)
    r <- rle(o)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    data.frame(
      order = r$values,
      start = starts,
      end = ends,
      stringsAsFactors = FALSE
    )
  }

  draw_one <- function(fig_w = 14, fig_h = 10) {
    ntip <- Ntip(tree)
    tree_h <- max(node.depth.edgelength(tree))
    cex_lab <- 0.52
    labs <- gsub("_", " ", tree$tip.label)

    # --- TEMP: ResearchEffort next to species name --------------------
    if (isTRUE(SHOW_RESEARCH_EFFORT) && !is.null(research_effort)) {
      re <- research_effort[tree$tip.label]
      labs <- paste0(labs, "  (", re, ")")
    }
    # --- end TEMP ----------------------------------------------------

    y_range <- range(c(0, y, INNO_THRESHOLD), na.rm = TRUE)
    if (diff(y_range) == 0) y_range <- y_range + c(-0.1, 0.1)
    y_pad <- diff(y_range) * 0.08
    y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)

    lab_w_in <- max(strwidth(labs, units = "inches", cex = cex_lab, font = 3))
    order_labs <- if (!is.null(orders)) {
      unique(format_order_label(as.character(orders)))
    } else {
      character(0)
    }
    order_lab_w_in <- if (length(order_labs)) {
      max(strwidth(order_labs, units = "inches", cex = 0.55, font = 2))
    } else {
      0
    }

    ylim_use <- c(-6.2, ntip + 1.5)
    plot_mar <- c(6.5, 0.8, 1.5, 0.8)
    # Centered panel inside the page (balances empty side space)
    fig_region <- c(0.05, 0.95, 0.05, 0.97)

    layout_from_plot <- function() {
      pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
      tip_x <- pp$xx[seq_len(ntip)]
      tip_y <- pp$yy[seq_len(ntip)]
      lab_w <- diff(grconvertX(c(0, lab_w_in), from = "inches", to = "user"))
      lolli_width <- diff(grconvertX(c(0, 2.3), from = "inches", to = "user"))
      tip_gap <- diff(grconvertX(c(0, 0.12), from = "inches", to = "user"))
      after_lab_gap <- diff(grconvertX(c(0, 0.18), from = "inches", to = "user"))
      order_gap <- diff(grconvertX(c(0, 0.22), from = "inches", to = "user"))
      order_lab_w <- diff(grconvertX(c(0, order_lab_w_in), from = "inches", to = "user"))
      right_pad <- diff(grconvertX(c(0, 0.15), from = "inches", to = "user"))
      x_label <- max(tip_x) + tip_gap
      lolli_origin <- x_label + lab_w + after_lab_gap
      x_order <- lolli_origin + lolli_width + order_gap
      x_right <- x_order + order_lab_w + right_pad
      list(
        tip_x = tip_x, tip_y = tip_y,
        lab_w = lab_w, lolli_width = lolli_width,
        tip_gap = tip_gap, after_lab_gap = after_lab_gap,
        x_label = x_label, lolli_origin = lolli_origin,
        x_order = x_order, order_lab_w = order_lab_w,
        x_right = x_right
      )
    }

    # Measure on null device with SAME size as final figure
    x_right <- tree_h * 2.2
    cur <- dev.cur()
    pdf(NULL, width = fig_w, height = fig_h)
    par(fig = fig_region, mar = plot_mar, xpd = FALSE, new = FALSE)
    plot.phylo(
      tree,
      type = "phylogram",
      direction = "rightwards",
      show.tip.label = FALSE,
      edge.color = "transparent",
      edge.width = 1,
      x.lim = c(0, x_right),
      y.lim = ylim_use,
      no.margin = FALSE
    )
    lay <- layout_from_plot()
    x_right <- lay$x_right
    x_pad <- 0.04 * x_right
    x_lim_use <- c(-x_pad, x_right + x_pad)
    # Second measure with final x.lim (usr ↔ inches depends on span)
    plot.phylo(
      tree,
      type = "phylogram",
      direction = "rightwards",
      show.tip.label = FALSE,
      edge.color = "transparent",
      edge.width = 1,
      x.lim = x_lim_use,
      y.lim = ylim_use,
      no.margin = FALSE
    )
    lay <- layout_from_plot()
    x_right <- lay$x_right
    x_pad <- 0.04 * x_right
    x_lim_use <- c(-x_pad, x_right + x_pad)
    dev.off()
    if (cur > 1) dev.set(cur)

    # Final draw, centered on the page
    par(fig = fig_region, mar = plot_mar, xpd = FALSE, new = FALSE)
    plot.phylo(
      tree,
      type = "phylogram",
      direction = "rightwards",
      show.tip.label = FALSE,
      edge.color = ec$col,
      edge.width = ec$lwd,
      x.lim = x_lim_use,
      y.lim = ylim_use,
      no.margin = FALSE
    )
    lay <- layout_from_plot()
    tip_x <- lay$tip_x
    tip_y <- lay$tip_y
    tip_gap <- lay$tip_gap
    x_label <- lay$x_label
    lolli_origin <- lay$lolli_origin
    lolli_width <- lay$lolli_width
    x_order <- lay$x_order

    scale_x <- function(v) (v - y_range[1]) / diff(y_range) * lolli_width
    x_zero <- lolli_origin + scale_x(0)

    # Zero baseline for lollipops
    segments(x_zero, min(tip_y) - 0.2, x_zero, max(tip_y) + 0.2, col = "grey55", lwd = 1)

    for (i in seq_len(ntip)) {
      sp <- tree$tip.label[i]
      tc <- tip_cols[[sp]]

      points(tip_x[i], tip_y[i], pch = 16, cex = 0.45, col = tc)

      text(
        x_label, tip_y[i],
        labels = labs[i],
        adj = c(0, 0.5),
        cex = cex_lab,
        font = 3,
        col = "grey15"
      )

      x_lab_end <- x_label + strwidth(labs[i], cex = cex_lab, font = 3)
      x_val <- lolli_origin + scale_x(y[i])

      segments(
        x_lab_end + tip_gap * 0.3, tip_y[i],
        lolli_origin, tip_y[i],
        col = "grey75", lwd = 0.7, lty = 3
      )
      segments(
        x_zero, tip_y[i],
        x_val, tip_y[i],
        col = tc, lwd = 1.7, lend = "butt"
      )
      points(
        x_val, tip_y[i],
        pch = 21,
        bg = tc,
        col = "white",
        cex = 1.0,
        lwd = 0.7
      )
    }

    # ---- Order strips to the right of innovation lollipops ----
    if (!is.null(orders)) {
      blocks <- order_blocks(orders[tree$tip.label])
      tick <- diff(grconvertX(c(0, 0.06), from = "inches", to = "user"))
      label_gap <- diff(grconvertX(c(0, 0.08), from = "inches", to = "user"))
      half <- 0.35
      # Guide lines start near tip leaves and run to the order bar
      x_guide_left <- min(tip_x) + 0.02 * (max(tip_x) - min(tip_x) + 1e-6)
      op <- par(xpd = NA)
      for (b in seq_len(nrow(blocks))) {
        i1 <- blocks$start[b]
        i2 <- blocks$end[b]
        y1 <- min(tip_y[i1], tip_y[i2]) - half
        y2 <- max(tip_y[i1], tip_y[i2]) + half
        ymid <- (y1 + y2) / 2
        # light dotted guides from tree side to order edge
        segments(
          x_guide_left, y1, x_order, y1,
          col = "grey80", lwd = 0.7, lty = 3
        )
        segments(
          x_guide_left, y2, x_order, y2,
          col = "grey80", lwd = 0.7, lty = 3
        )
        # vertical bar spanning the order's tips
        segments(x_order, y1, x_order, y2, col = "grey35", lwd = 2.4, lend = "butt")
        # end ticks toward the plot
        segments(x_order, y1, x_order - tick, y1, col = "grey35", lwd = 2)
        segments(x_order, y2, x_order - tick, y2, col = "grey35", lwd = 2)
        text(
          x_order + label_gap, ymid,
          labels = format_order_label(blocks$order[b]),
          adj = c(0, 0.5),
          cex = 0.55,
          font = 2,
          col = "grey25"
        )
      }
      par(op)
    }

    # Axis closer to the tip band
    axis_y <- -2.0
    tick_h <- 0.45
    segments(
      lolli_origin, axis_y,
      lolli_origin + lolli_width, axis_y,
      col = "grey30", lwd = 1.3
    )
    ticks <- pretty(y_range, n = 5)
    ticks <- ticks[ticks >= y_range[1] & ticks <= y_range[2]]
    for (tv in ticks) {
      tx <- lolli_origin + scale_x(tv)
      segments(tx, axis_y, tx, axis_y - tick_h, col = "grey30")
      text(
        tx, axis_y - 1.15 * tick_h,
        labels = tv, cex = 0.72, col = "grey20", adj = c(0.5, 1)
      )
    }
    text(
      lolli_origin + lolli_width / 2,
      axis_y - 3.6 * tick_h,
      labels = axis_label,
      cex = 1.0,
      font = 3,
      col = "grey15"
    )

    # Continuous ASR colour legend near the tree base
    usr <- par("usr")
    leg_w <- 0.22 * (usr[2] - usr[1])
    leg_h <- 0.9
    draw_asr_legend(
      x_left = usr[1] + 0.03 * (usr[2] - usr[1]),
      y_bottom = -4.6,
      width = leg_w,
      height = leg_h,
      rng = asr$range,
      lut = asr$lut,
      title = "ASR phenotype"
    )
  }

  fig_w <- 15.5
  fig_h <- max(10, 0.14 * Ntip(tree) + 5.5)

  pdf(out_pdf, width = fig_w, height = fig_h)
  draw_one(fig_w = fig_w, fig_h = fig_h)
  dev.off()

  png(
    out_png,
    width = round(fig_w * 220),
    height = round(fig_h * 220),
    res = 220
  )
  draw_one(fig_w = fig_w, fig_h = fig_h)
  dev.off()
}

# ---- main ----
main <- function() {
  message("Loading ", TREE_FILE, " + ", TRAITS_FILE)
  inp <- load_and_match()
  tree <- inp$tree
  y <- inp$y
  message("Matched species: ", Ntip(tree))

  message("Running RRphylo ...")
  RR <- RRphylo(tree = tree, y = y, clus = CLUS)
  saveRDS(RR, file.path(OUT_DIR, "rrphylo_model.rds"))

  branch_df <- extract_branch_rates(RR)
  write.csv(branch_df, file.path(OUT_DIR, "branch_rates.csv"), row.names = FALSE)

  message("Drawing straight phylogram + ASR branch colours ...")
  research_effort <- NULL
  if ("ResearchEffort" %in% names(inp$dat)) {
    research_effort <- setNames(
      as.numeric(inp$dat$ResearchEffort), inp$dat[[SPECIES_COL]]
    )[RR$tree$tip.label]
  }
  orders <- NULL
  if ("Order" %in% names(inp$dat)) {
    orders <- setNames(
      as.character(inp$dat$Order), inp$dat[[SPECIES_COL]]
    )[RR$tree$tip.label]
  }

  asr_fit <- branch_asr_colors(RR$tree, y[RR$tree$tip.label], lut = ASR_LUT)
  write.csv(
    data.frame(
      node = as.integer(names(asr_fit$anc)),
      ancestral_state = as.numeric(asr_fit$anc),
      stringsAsFactors = FALSE
    ),
    file.path(OUT_DIR, "asr_ancestral_states.csv"),
    row.names = FALSE
  )

  plot_straight(
    tree = RR$tree,
    y = y[RR$tree$tip.label],
    out_pdf = file.path(OUT_DIR, "phylogeny_rate_shifts_innovation_cont.pdf"),
    out_png = file.path(OUT_DIR, "phylogeny_rate_shifts_innovation_cont.png"),
    research_effort = research_effort,
    orders = orders
  )

  message("Done → ", OUT_DIR)
}

if (sys.nframe() == 0) main()


