#!/usr/bin/env Rscript
# RRphylo + sparse search.shift + tip-threshold colouring
# =====================================================================
# Previous visualization variant (kept separate from ASR continuous colouring):
#   - tip / terminal branch colours by phenotype thresholds + low ResearchEffort
#   - sparse search.shift (inno / non-inno)
#
# Tip visual classes:
#   high_inno : phenotype >= INNO_THRESHOLD
#   low_inno  : phenotype <= LOW_INNO_THRESHOLD  (branch stays neutral)
#   mid_inno  : between thresholds
#   low_re    : ResearchEffort <= LOW_RE_THRESHOLD (overrides)
#
# Outputs: results/rate_shifts_tip_threshold/

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
OUT_DIR <- file.path(PROJECT_ROOT, "results", "rate_shifts_tip_threshold")

# TREE_FILE <- "tree_clootl_overlap_500.nwk"
# TRAITS_FILE <- "inno_db_500_clootl_overlap.csv"

TREE_FILE <- "tree_astralpro.l6.tre_overlap_0_vgp.nwk"
TRAITS_FILE <- "inno_db_0_astralpro.l6.tre_overlap_vgp.csv"

#TREE_FILE <- "tree_сlootl_overlap_0_vgp.nwk"
#TRAITS_FILE <- "inno_db_0_сlootl_overlap_vgp.csv"

SPECIES_COL <- "Species"
PHENOTYPE <- "TOTALINNOVATIONS2025_ResEff"
INNO_THRESHOLD <- 1.5
# Low-innovation visual class (blue): phenotype at or below this value
LOW_INNO_THRESHOLD <- -0.5

# TEMP: tip labels show ResearchEffort on the right. Set to FALSE to hide.
SHOW_RESEARCH_EFFORT <- TRUE
# Tips with ResearchEffort at or below this are drawn pale gray
LOW_RE_THRESHOLD <- 15

SHIFT_NREP <- 1000
CLUS <- 2 / parallel::detectCores()
P_ACCEL <- 0.975
P_DECEL <- 0.025

# Visual classes for tips (not rate-shift stats)
COL_HIGH_INNO <- "#E69F00"   # yellow/orange: Technical Innovation >= INNO_THRESHOLD
COL_LOW_INNO  <- "#56B4E9"   # soft blue: Technical Innovation <= LOW_INNO_THRESHOLD
COL_MID_INNO  <- "#666666"   # neutral: between thresholds
COL_LOW_RE    <- "#D0D0D0"   # pale gray: ResearchEffort <= LOW_RE_THRESHOLD
COL_BG        <- "#B0B0B0"   # internal / background branches

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
  # dat <- dat %>% filter(ResearchEffort >= 15)
  
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

make_inno_state <- function(y, threshold = INNO_THRESHOLD) {
  state <- ifelse(y >= threshold, "inno", "non-inno")
  names(state) <- names(y)
  state
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

# Interpret binary sparse result for the "inno" group vs "non-inno"
interpret_sparse <- function(SS, state, RR) {
  sr <- SS$state.results
  rn <- rownames(sr)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(sr)))

  # Binary: one row named like "inno_non-inno" or "non-inno_inno"
  # rate.difference = mean|rate|(second) - mean|rate|(first) from name order
  # in RRphylo: names(status.diff) <- paste(combn(names(sta), 2)[2:1], collapse="_")
  # and status.diff = apply(combn(sta,2), 2, diff) = sta[2]-sta[1]
  # so for name "B_A", difference is sta[B]-sta[A].
  p <- as.numeric(sr$p.status.diff[1])
  dlt <- as.numeric(sr$rate.difference[1])
  cmp <- rn[1]

  # Positive dlt with name "inno_non-inno" => |rate|_inno > |rate|_non-inno
  # Positive dlt with name "non-inno_inno" => |rate|_non-inno > |rate|_inno
  inno_higher <- if (cmp == "inno_non-inno") {
    dlt > 0
  } else if (cmp == "non-inno_inno") {
    dlt < 0
  } else {
    # fallback: compare tip mean absolute rates directly
    tip_rates <- abs(RR$rates[names(state), 1])
    mean(tip_rates[state == "inno"]) > mean(tip_rates[state == "non-inno"])
  }

  # p >= 0.975: first group in the randomization ranking is extreme high difference
  # Same convention as clade: p>0.975 => larger rates; p<0.025 => smaller
  if (p >= P_ACCEL) {
    direction <- if (inno_higher) "accelerated" else "decelerated"
  } else if (p <= P_DECEL) {
    # Significant in the opposite direction of the signed difference
    direction <- if (inno_higher) "decelerated" else "accelerated"
  } else {
    direction <- "ns"
  }

  # For binary, the p-value tests the observed signed difference.
  # Cleaner: use sign of (mean_inno - mean_noninno) with p:
  tip_rates <- abs(RR$rates[RR$tree$tip.label, 1])
  st <- state[RR$tree$tip.label]
  mean_inno <- mean(tip_rates[st == "inno"])
  mean_non <- mean(tip_rates[st == "non-inno"])
  inno_minus_non <- mean_inno - mean_non

  if (p >= P_ACCEL) {
    direction <- if (inno_minus_non > 0) "accelerated" else "decelerated"
  } else if (p <= P_DECEL) {
    direction <- if (inno_minus_non > 0) "decelerated" else "accelerated"
  } else {
    direction <- "ns"
  }

  list(
    comparison = cmp,
    rate_difference = dlt,
    p_value = p,
    direction = direction,
    mean_abs_rate_inno = mean_inno,
    mean_abs_rate_non_inno = mean_non,
    inno_minus_non = inno_minus_non
  )
}

build_sparse_tip_table <- function(RR, state, y, sparse_summary) {
  tips_inno <- names(state)[state == "inno"]
  rates <- RR$rates[, 1]
  mean_abs <- mean(abs(rates))
  if (!length(tips_inno)) {
    return(data.frame(
      node = integer(0), branch = character(0), species = character(0),
      shift_type = character(0), evolutionary_rate = numeric(0),
      relative_rate = numeric(0), phenotype = numeric(0),
      rate_difference = numeric(0), p_value = numeric(0),
      direction = character(0), sparse_state = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(tips_inno, function(sp) {
    data.frame(
      node = as.integer(match(sp, RR$tree$tip.label)),
      branch = sp,
      species = sp,
      shift_type = "terminal",
      evolutionary_rate = as.numeric(rates[[sp]]),
      relative_rate = as.numeric(rates[[sp]]) / mean_abs,
      phenotype = as.numeric(y[[sp]]),
      rate_difference = sparse_summary$rate_difference,
      p_value = sparse_summary$p_value,
      direction = sparse_summary$direction,
      sparse_state = "inno",
      stringsAsFactors = FALSE
    )
  }))
}

# Tip visual class: low ResearchEffort overrides innovation colouring
classify_tip_style <- function(y, research_effort = NULL) {
  cls <- rep("mid_inno", length(y))
  names(cls) <- names(y)
  cls[y >= INNO_THRESHOLD] <- "high_inno"
  cls[y <= LOW_INNO_THRESHOLD] <- "low_inno"
  if (!is.null(research_effort)) {
    re <- research_effort[names(y)]
    low_re <- !is.na(re) & as.numeric(re) <= LOW_RE_THRESHOLD
    cls[low_re] <- "low_re"
  }
  cls
}

style_colors <- function(cls) {
  cols <- c(
    high_inno = COL_HIGH_INNO,
    low_inno = COL_LOW_INNO,
    mid_inno = COL_MID_INNO,
    low_re = COL_LOW_RE
  )
  unname(cols[cls])
}

color_tip_terminals <- function(tree, tip_class) {
  ne <- Nedge(tree)
  edge_col <- rep(COL_BG, ne)
  edge_lwd <- rep(1.1, ne)
  tip_class <- tip_class[tree$tip.label]
  tip_cols <- style_colors(tip_class)
  names(tip_cols) <- tree$tip.label
  for (sp in tree$tip.label) {
    tip_i <- match(sp, tree$tip.label)
    e <- which(tree$edge[, 2] == tip_i)
    if (!length(e)) next
    edge_col[e] <- tip_cols[[sp]]
    edge_lwd[e] <- switch(
      tip_class[[sp]],
      high_inno = 2.8,
      low_inno = 1.1,
      low_re = 1.2,
      1.1
    )
    # Low-innovation: keep branch neutral (only label/lollipop are blue)
    if (identical(tip_class[[sp]], "low_inno")) {
      edge_col[e] <- COL_BG
    }
  }
  list(col = edge_col, lwd = edge_lwd, tip_cols = tip_cols)
}

plot_straight <- function(
  tree, y, state, direction, out_pdf, out_png,
  research_effort = NULL,
  orders = NULL
) {
  y <- y[tree$tip.label]
  state <- state[tree$tip.label]
  if (!is.null(research_effort)) {
    research_effort <- research_effort[tree$tip.label]
  }
  if (!is.null(orders)) {
    orders <- setNames(as.character(orders[tree$tip.label]), tree$tip.label)
  }

  tip_class <- classify_tip_style(y, research_effort)
  ec <- color_tip_terminals(tree, tip_class)
  tip_cols <- ec$tip_cols
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

    y_range <- range(c(0, y, INNO_THRESHOLD, LOW_INNO_THRESHOLD), na.rm = TRUE)
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
      is_faint <- tip_class[[sp]] == "low_re"

      points(tip_x[i], tip_y[i], pch = 16, cex = 0.55, col = tc)

      text(
        x_label, tip_y[i],
        labels = labs[i],
        adj = c(0, 0.5),
        cex = cex_lab,
        font = 3,
        col = tc
      )

      x_lab_end <- x_label + strwidth(labs[i], cex = cex_lab, font = 3)
      x_val <- lolli_origin + scale_x(y[i])

      segments(
        x_lab_end + tip_gap * 0.3, tip_y[i],
        lolli_origin, tip_y[i],
        col = if (is_faint) COL_LOW_RE else "grey68",
        lwd = 0.8, lty = 3
      )
      segments(
        x_zero, tip_y[i],
        x_val, tip_y[i],
        col = tc, lwd = if (is_faint) 1.2 else 1.8, lend = "butt"
      )
      points(
        x_val, tip_y[i],
        pch = 21,
        bg = tc,
        col = "white",
        cex = if (is_faint) 0.85 else 1.05,
        lwd = 0.7
      )
    }

    # Innovation threshold marker on the lollipop axis
    x_thr <- lolli_origin + scale_x(INNO_THRESHOLD)
    segments(
      x_thr, min(tip_y) - 0.2, x_thr, max(tip_y) + 0.2,
      col = "grey40", lwd = 1, lty = 2
    )

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

    # Tip-class legend near the tree base
    usr <- par("usr")
    legend(
      x = usr[1] + 0.02 * (usr[2] - usr[1]),
      y = min(tip_y) - 0.8,
      legend = c(
        sprintf("Technical Innovation >= %.1f", INNO_THRESHOLD),
        sprintf("Technical Innovation <= %.1f", LOW_INNO_THRESHOLD),
        "Between thresholds",
        sprintf("ResearchEffort <= %s (faded)", LOW_RE_THRESHOLD),
        "Background branch",
        sprintf("Inno threshold = %.1f", INNO_THRESHOLD)
      ),
      col = c(COL_HIGH_INNO, COL_LOW_INNO, COL_MID_INNO, COL_LOW_RE, COL_BG, "grey40"),
      pch = c(21, 21, 21, 21, NA, NA),
      pt.bg = c(COL_HIGH_INNO, COL_LOW_INNO, COL_MID_INNO, COL_LOW_RE, NA, NA),
      lwd = c(NA, NA, NA, NA, 1.1, 1),
      lty = c(NA, NA, NA, NA, 1, 2),
      pt.cex = c(1.1, 1.1, 1.0, 0.9, NA, NA),
      bty = "n",
      cex = 0.68
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

  state <- make_inno_state(y, INNO_THRESHOLD)
  message(
    "State: inno=", sum(state == "inno"),
    ", non-inno=", sum(state == "non-inno"),
    " (threshold=", INNO_THRESHOLD, ")"
  )
  write.csv(
    data.frame(
      species = names(state),
      phenotype = as.numeric(y[names(state)]),
      sparse_state = as.character(state),
      stringsAsFactors = FALSE
    ),
    file.path(OUT_DIR, "sparse_state_assignment.csv"),
    row.names = FALSE
  )

  message("Running RRphylo ...")
  RR <- RRphylo(tree = tree, y = y, clus = CLUS)
  saveRDS(RR, file.path(OUT_DIR, "rrphylo_model.rds"))

  branch_df <- extract_branch_rates(RR)
  branch_df$sparse_state <- ifelse(
    !is.na(branch_df$species),
    state[branch_df$species],
    NA_character_
  )
  write.csv(branch_df, file.path(OUT_DIR, "branch_rates.csv"), row.names = FALSE)

  message('Running search.shift(RR, status.type = "sparse", state = inno/non-inno) ...')
  SS <- search.shift(
    RR,
    status.type = "sparse",
    state = state,
    nrep = SHIFT_NREP
  )
  saveRDS(SS, file.path(OUT_DIR, "search_shift_sparse.rds"))

  write.csv(
    data.frame(
      comparison = rownames(SS$state.results),
      SS$state.results,
      stringsAsFactors = FALSE
    ),
    file.path(OUT_DIR, "sparse_search_shift_state_results.csv"),
    row.names = FALSE
  )

  sparse_summary <- interpret_sparse(SS, state, RR)
  message(
    "Sparse result: comparison=", sparse_summary$comparison,
    " rate.difference=", round(sparse_summary$rate_difference, 4),
    " p=", sparse_summary$p_value,
    " direction(inno)=", sparse_summary$direction
  )
  write.csv(
    as.data.frame(sparse_summary, stringsAsFactors = FALSE),
    file.path(OUT_DIR, "sparse_shift_summary.csv"),
    row.names = FALSE
  )

  sparse_tab <- build_sparse_tip_table(RR, state, y, sparse_summary)
  write.csv(sparse_tab, file.path(OUT_DIR, "sparse_rate_shifts.csv"), row.names = FALSE)
  message("Inno tips in sparse_rate_shifts.csv: ", nrow(sparse_tab))
  if (nrow(sparse_tab)) {
    print(sparse_tab[, c("species", "phenotype", "evolutionary_rate", "direction", "p_value")])
  }

  message("Drawing straight phylogram + tip-threshold colours ...")
  # ResearchEffort is optional tip annotation only
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

  plot_straight(
    tree = RR$tree,
    y = y[RR$tree$tip.label],
    state = state,
    direction = sparse_summary$direction,
    out_pdf = file.path(OUT_DIR, "phylogeny_rate_shifts_innovation.pdf"),
    out_png = file.path(OUT_DIR, "phylogeny_rate_shifts_innovation.png"),
    research_effort = research_effort,
    orders = orders
  )

  pdf(file.path(OUT_DIR, "sparse_plotShift_states.pdf"), width = 7, height = 5)
  tryCatch({
    ps <- plotShift(RR = RR, SS = SS, state = state)
    ps$plotStates()
  }, error = function(e) {
    message("plotShift: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, conditionMessage(e))
  })
  dev.off()

  message("Done → ", OUT_DIR)
}

if (sys.nframe() == 0) main()


