#!/usr/bin/env Rscript
# Phylogenetic analysis of bird innovativeness (based on phenotype + tree only).

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phylolm)
  library(ggplot2)
  library(ggtree)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(patchwork)
})

# ---- parameters ----
# Run this script with the working directory set to the folder that
# contains data/ (i.e. Phylogenetic_tests). Override PROJECT_ROOT below if
# you need to point at a different location.
PROJECT_ROOT <- normalizePath(".", wins = FALSE)
PROJECT_ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/cognitive_evolution/"

# PROJECT_ROOT <- "/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/src/Phylogenetic_tests"

DATA_DIR <- file.path(PROJECT_ROOT, "data")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results_500")

# Pair each traits CSV with the matching tree (same ResEff threshold).
# tree_clootl_overlap.nwk is the small 76-tip tree (= _500), not the 100 set.
TREE_FILE   <- "tree_clootl_overlap_500.nwk"
TRAITS_FILE <- "inno_db_500_clootl_overlap.csv"

# TREE_FILE   <- "tree_clootl_overlap_100.nwk"
# TRAITS_FILE <- "inno_db_100_clootl_overlap.csv"

# TREE_FILE   <- "tree_clootl_overlap_50.nwk"
# TRAITS_FILE <- "inno_db_50_clootl_overlap.csv"

SPECIES_COL <- "Species"

PHENOTYPES <- c(
  #"TOTALINNOVATIONS2025_ResEff",
  #"FOODINNO2025_ResEff",
  "TECHINNO2025_ResEff"
)
COVARIATES <- c(
  "Brain_size", "Relative_brain_size", "Mass", "DietBreadth",
  "HabitatBreadth", "GenerationLength", "Range.Size",
  "UrbanFULL", "Migration", "Trophic_level"
)
SIGNAL_TRAITS <- c(
  "TECHINNO2025_ResEff","FOODINNO2025_ResEff", "TOTALINNOVATIONS2025_ResEff",
  "Brain_size", "Relative_brain_size", "Mass",
  "DietBreadth", "HabitatBreadth", "GenerationLength", "Range.Size"
)
MULTIVARIATE_COVARIATES <- c(
  "Relative_brain_size", "Mass", "DietBreadth",
  "HabitatBreadth", "GenerationLength"
)
PCA_TRAITS <- c(
  "Relative_brain_size", "Mass", "DietBreadth",
  "HabitatBreadth", "GenerationLength", "Range.Size"
)
CORRELATION_TRAITS <- c(
  PHENOTYPES,
  "Brain_size", "Relative_brain_size", "Mass",
  "DietBreadth", "HabitatBreadth", "GenerationLength", "Range.Size"
)
FOREGROUND_SD <- 1
SIGNAL_ALPHA <- 0.05
# PGLS predictors are z-scaled inside pgls_fit() for numerical stability
# (Range.Size otherwise causes a singular design matrix in phylolm).

TRAIT_LABELS <- c(
  Brain_size = "Brain size",
  Relative_brain_size = "Relative brain size",
  Mass = "Body mass",
  GenerationLength = "Generation length",
  Range.Size = "Range size",
  HabitatBreadth = "Habitat breadth",
  DietBreadth = "Diet breadth",
  UrbanFULL = "Urban tolerance",
  Migration = "Migration",
  Trophic_level = "Trophic level",
  TOTALINNOVATIONS2025_ResEff = "Total innovation",
  FOODINNO2025_ResEff = "Food innovation",
  TECHINNO2025_ResEff = "Technical innovation"
)

labelize <- function(x) {
  out <- unname(TRAIT_LABELS[x])
  ifelse(is.na(out), gsub("_", " ", x), out)
}

# ---- helpers ----
dir_ok <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Columns that must be numeric for phylosig / PGLS / PCA. Larger CSVs can
# force e.g. Range.Size to character when a few cells are " NA" (leading
# space); phylosig then dies with "invC %*% x : requires numeric/...".
NUMERIC_ANALYSIS_COLS <- unique(c(
  PHENOTYPES,
  setdiff(COVARIATES, c("Migration", "Trophic_level")),
  SIGNAL_TRAITS,
  MULTIVARIATE_COVARIATES,
  PCA_TRAITS,
  CORRELATION_TRAITS
))

coerce_numeric_cols <- function(dat, cols = NUMERIC_ANALYSIS_COLS) {
  for (col in cols) {
    if (!col %in% names(dat) || is.numeric(dat[[col]])) next
    raw <- trimws(as.character(dat[[col]]))
    raw[raw %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
    dat[[col]] <- suppressWarnings(as.numeric(raw))
  }
  dat
}

load_inputs <- function() {
  tree <- read.tree(file.path(DATA_DIR, TREE_FILE))
  dat <- read.csv(
    file.path(DATA_DIR, TRAITS_FILE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if ("Unnamed: 0" %in% names(dat)) {
    dat <- dat[, setdiff(names(dat), "Unnamed: 0"), drop = FALSE]
  }
  dat <- coerce_numeric_cols(dat)
  list(tree = tree, dat = dat)
}

prepare_data <- function(tree, dat) {
  rownames(dat) <- dat[[SPECIES_COL]]

  tree_species <- sort(tree$tip.label)
  table_species <- sort(unique(dat[[SPECIES_COL]]))

  in_both <- intersect(tree_species, table_species)
  missing_in_table <- setdiff(tree_species, table_species)
  missing_in_tree <- setdiff(table_species, tree_species)

  report <- data.frame(
    metric = c(
      "species_in_tree", "species_in_table",
      "species_matched", "missing_in_table", "missing_in_tree"
    ),
    value = c(
      length(tree_species), length(table_species),
      length(in_both), length(missing_in_table), length(missing_in_tree)
    ),
    stringsAsFactors = FALSE
  )

  # Restrict tree + table to the matched species and align their order.
  # Without this, `dat` keeps every row from the CSV (e.g. 809) while `tree`
  # only has the matched tips (e.g. 76); downstream code that does
  # tree$tip.label[!complete.cases(dat[, traits])] then indexes the tree's
  # tip labels using a logical vector whose positions refer to CSV row
  # order, not tip order -- silently dropping the wrong tips and eventually
  # producing degenerate (non-numeric) covariance matrices in phylosig().
  tree <- drop.tip(tree, setdiff(tree$tip.label, in_both))
  dat <- dat[tree$tip.label, , drop = FALSE]

  list(
    tree = tree,
    dat = dat,
    report = report,
    missing_in_table = missing_in_table,
    missing_in_tree = missing_in_tree
  )
}

trait_vector <- function(dat, tree, trait) {
  x <- dat[[trait]]
  if (!is.numeric(x)) {
    x <- suppressWarnings(as.numeric(trimws(as.character(x))))
  }
  names(x) <- dat[[SPECIES_COL]]
  x <- x[tree$tip.label]
  if (all(is.na(x))) {
    stop("Trait '", trait, "' is all NA after numeric coercion / tip matching")
  }
  x
}

prune_complete <- function(tree, dat, traits) {
  # Index dat by tree tip labels (name-based) before checking completeness,
  # so `ok` always lines up with tree$tip.label positionally -- regardless
  # of dat's row order or whether it contains extra species not in the tree.
  dat_tip <- dat[tree$tip.label, traits, drop = FALSE]
  ok <- complete.cases(dat_tip)
  drop_tips <- tree$tip.label[!ok]
  tree_pr <- if (length(drop_tips) > 0) drop.tip(tree, drop_tips) else tree
  dat_pr <- dat[tree_pr$tip.label, , drop = FALSE]
  list(tree = tree_pr, dat = dat_pr)
}

compute_phylogenetic_signal <- function(tree, dat, traits) {
  rows <- lapply(traits, function(tr) {
    cc <- prune_complete(tree, dat, tr)
    x <- trait_vector(cc$dat, cc$tree, tr)
    k <- phytools::phylosig(cc$tree, x, method = "K", test = TRUE)
    l <- phytools::phylosig(cc$tree, x, method = "lambda", test = TRUE)
    data.frame(
      Trait = tr,
      K = unname(k$K),
      K_p = k$P,
      lambda = unname(l$lambda),
      lambda_p = l$P,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

compute_ancestral_states <- function(tree, x) {
  anc <- fastAnc(tree, x)
  data.frame(
    node = as.integer(names(anc)),
    state = as.numeric(anc),
    stringsAsFactors = FALSE
  )
}

compute_branch_changes <- function(tree, x, anc) {
  ntip <- Ntip(tree)
  edge <- tree$edge
  tip_vals <- setNames(as.numeric(x), names(x))

  node_value <- function(node) {
    if (node <= ntip) {
      return(tip_vals[tree$tip.label[node]])
    }
    anc[[as.character(node)]]
  }

  out <- data.frame(
    branch_id = paste0("edge_", seq_len(nrow(edge))),
    parent = edge[, 1],
    child = edge[, 2],
    parent_value = vapply(edge[, 1], node_value, numeric(1)),
    child_value = vapply(edge[, 2], node_value, numeric(1)),
    branch_length = tree$edge.length,
    stringsAsFactors = FALSE
  )
  out$delta <- out$child_value - out$parent_value
  out
}

identify_foreground <- function(branch_changes, sd_mult = FOREGROUND_SD) {
  threshold <- mean(branch_changes$delta) + sd_mult * sd(branch_changes$delta)
  fg <- branch_changes[branch_changes$delta > threshold, , drop = FALSE]
  list(foreground = fg, threshold = threshold)
}

# Project the branch-level foreground/background split onto tips: a tip is
# part of a "foreground lineage" if it descends from (or is the child of)
# any foreground branch. This mirrors how foreground/background branches
# are used in molecular selection tests, but here the "selection" criterion
# is a lineage-specific jump in reconstructed innovation score.
foreground_tip_labels <- function(tree, foreground) {
  if (is.null(foreground) || nrow(foreground) == 0) return(character(0))
  ntip <- Ntip(tree)
  desc_nodes <- unique(unlist(lapply(foreground$child, function(nd) {
    if (nd <= ntip) nd else c(nd, phytools::getDescendants(tree, nd))
  })))
  tip_nodes <- desc_nodes[desc_nodes <= ntip]
  tree$tip.label[tip_nodes]
}

pgls_fit <- function(tree, dat, formula, scale_terms = NULL) {
  d <- dat
  if (!is.null(scale_terms)) {
    for (tr in scale_terms) {
      if (tr %in% names(d) && is.numeric(d[[tr]])) {
        d[[tr]] <- as.numeric(scale(d[[tr]]))
      }
    }
  }
  phylolm(
    formula,
    data = d,
    phy = tree,
    model = "lambda",
    boot = 0
  )
}

extract_slope <- function(fit, term) {
  coefs <- summary(fit)$coefficients
  if (!term %in% rownames(coefs)) {
    return(data.frame(beta = NA_real_, SE = NA_real_, p = NA_real_, R2 = NA_real_))
  }
  data.frame(
    beta = coefs[term, "Estimate"],
    SE = coefs[term, "StdErr"],
    p = coefs[term, "p.value"],
    R2 = summary(fit)$adj.r.squared,
    stringsAsFactors = FALSE
  )
}

run_univariate_pgls <- function(tree, dat, phenotype, covariates) {
  rows <- lapply(covariates, function(tr) {
    d <- dat
    if (tr %in% c("Migration", "Trophic_level")) {
      d[[tr]] <- as.numeric(factor(d[[tr]]))
    }
    fit <- tryCatch(
      pgls_fit(tree, d, as.formula(paste(phenotype, "~", tr)), scale_terms = tr),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(Trait = tr, beta = NA_real_, SE = NA_real_, p = NA_real_, R2 = NA_real_))
    }
    out <- extract_slope(fit, tr)
    out$Trait <- tr
    out
  })
  bind_rows(rows) |>
    transmute(Trait, beta, SE, p, R2) |>
    arrange(p)
}

run_multivariate_pgls <- function(tree, dat, phenotype, predictors) {
  rhs <- paste(predictors, collapse = " + ")
  fit <- tryCatch(
    pgls_fit(
      tree, dat,
      as.formula(paste(phenotype, "~", rhs)),
      scale_terms = predictors
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    out <- data.frame(
      term = predictors,
      beta = NA_real_,
      SE = NA_real_,
      p = NA_real_,
      stringsAsFactors = FALSE
    )
    attr(out, "R2") <- NA_real_
    return(out)
  }
  coefs <- summary(fit)$coefficients
  keep <- rownames(coefs) != "(Intercept)"
  out <- data.frame(
    term = rownames(coefs)[keep],
    beta = coefs[keep, "Estimate"],
    SE = coefs[keep, "StdErr"],
    p = coefs[keep, "p.value"],
    stringsAsFactors = FALSE
  )
  attr(out, "R2") <- summary(fit)$adj.r.squared
  out
}

run_pca <- function(dat, traits, phenotype) {
  cc <- complete.cases(dat[, traits, drop = FALSE])
  x <- scale(dat[cc, traits, drop = FALSE])
  pca <- prcomp(x, center = FALSE, scale. = FALSE)
  coords <- data.frame(
    species = dat[[SPECIES_COL]][cc],
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    innovation = dat[[phenotype]][cc],
    stringsAsFactors = FALSE
  )
  list(coords = coords, pca = pca)
}

save_genome_inputs <- function(dat, phenotype, out_dir) {
  write_tsv(
    data.frame(
      species = dat[[SPECIES_COL]],
      innovation = dat[[phenotype]],
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "continuous_phenotype.tsv")
  )
}

# ---- figures (original) ----
plot_contmap <- function(tree, x, out_file) {
  pdf(out_file, width = 12, height = 16)
  contMap(
    tree, x,
    lut = c("blue", "white", "red"),
    outline = FALSE,
    lwd = 3,
    fsize = 0.35
  )
  dev.off()
}

plot_ancestral_tree <- function(tree, anc_df, out_file) {
  anc_df$label <- round(anc_df$state, 2)
  # ggtree's internal fortify.phylo() method declares named ... parameters
  # (as.Date, yscale_mapping, hang, etc.) that it never uses for a plain
  # rectangular layout; ggplot2 >= 3.5's check_dots_used() flags this as a
  # warning even though nothing is actually wrong with the plot. It's an
  # upstream ggtree/ggplot2 compatibility quirk, not a bug in this script,
  # so it's suppressed here rather than left to confuse readers of the log.
  p <- suppressWarnings(
    ggtree(tree, branch.length = "none", layout = "rectangular") %<+% anc_df +
      geom_nodelab(aes(label = label), size = 1.2, color = "darkblue") +
      theme_tree2() +
      ggtitle("Ancestral state reconstruction")
  )
  ggsave(out_file, p, width = 12, height = 16)
}

plot_foreground_tree <- function(tree, branch_changes, foreground, out_file) {
  ne <- Nedge(tree)
  edge_col <- rep("gray60", ne)
  edge_w <- rep(0.6, ne)
  fg_ids <- foreground$branch_id
  idx <- match(fg_ids, branch_changes$branch_id)
  idx <- idx[!is.na(idx)]
  edge_col[idx] <- "red"
  edge_w[idx] <- 2.5

  pdf(out_file, width = 12, height = 16)
  plot(
    tree,
    type = "phylogram",
    edge.color = edge_col,
    edge.width = edge_w,
    cex = 0.35,
    no.margin = TRUE
  )
  title(main = "Foreground branches (independent innovation candidates)")
  dev.off()
}

plot_delta_histogram <- function(branch_changes, threshold, out_file) {
  p <- ggplot(branch_changes, aes(delta)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "Distribution of branch-wise trait changes",
      x = expression(Delta ~ "(child - parent)"),
      y = "Count"
    ) +
    theme_bw()
  ggsave(out_file, p, width = 7, height = 5)
}

plot_parent_child_scatter <- function(branch_changes, foreground, out_file) {
  branch_changes$foreground <- branch_changes$branch_id %in% foreground$branch_id
  p <- ggplot(branch_changes, aes(parent_value, child_value, color = foreground)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
    labs(
      title = "Parent vs child values along branches",
      x = "Parent value",
      y = "Child value",
      color = "Foreground"
    ) +
    theme_bw()
  ggsave(out_file, p, width = 6, height = 5)
}

plot_forest <- function(pgls_uni, out_file) {
  pgls_uni <- pgls_uni |>
    mutate(
      lower = beta - 1.96 * SE,
      upper = beta + 1.96 * SE
    )
  p <- ggplot(pgls_uni, aes(x = beta, y = reorder(Trait, beta))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.2) +
    geom_point(size = 2) +
    labs(
      title = "Univariate PGLS coefficients (95% CI)",
      x = expression(beta),
      y = NULL
    ) +
    theme_bw()
  ggsave(out_file, p, width = 7, height = 5)
}

plot_correlation_matrix <- function(dat, traits, out_file) {
  cor_mat <- cor(dat[, traits, drop = FALSE], use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "r")
  p <- ggplot(cor_df, aes(Var1, Var2, fill = r)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    coord_fixed() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    ) +
    labs(title = "Trait correlation matrix", x = NULL, y = NULL, fill = "r")
  ggsave(out_file, p, width = 9, height = 8)
}

plot_pca <- function(pca_coords, out_file) {
  p <- ggplot(pca_coords, aes(PC1, PC2, color = innovation)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = median(pca_coords$innovation, na.rm = TRUE)
    ) +
    labs(
      title = "PCA of ecological traits",
      x = "PC1",
      y = "PC2",
      color = "Innovation"
    ) +
    theme_bw()
  ggsave(out_file, p, width = 7, height = 6)
}

# ---- figures (new) ----

# 1) Foreground vs background lineage boxplot, per innovation phenotype.
# Tips descending from a foreground branch (fast lineage-specific increase
# in the phenotype, as flagged by identify_foreground()) are compared
# against the rest of the tree for a chosen trait (often Range.Size),
# using a Wilcoxon rank-sum test.
#
# log_transform = TRUE plots log10(trait); useful for geographic range.
plot_foreground_background_box <- function(
  tree, dat, phenotype, foreground, out_file = NULL,
  log_transform = TRUE,
  point_alpha = 0.18,
  ylim = c(3, 8)
) {
  fg_tips <- foreground_tip_labels(tree, foreground)

  d <- dat[!is.na(dat[[phenotype]]), , drop = FALSE]
  if (log_transform) {
    d <- d[is.finite(d[[phenotype]]) & d[[phenotype]] > 0, , drop = FALSE]
    d$y_plot <- log10(d[[phenotype]])
    ylab <- bquote(log[10](.(labelize(phenotype))))
  } else {
    d$y_plot <- d[[phenotype]]
    ylab <- labelize(phenotype)
  }

  d$Lineage <- ifelse(
    d[[SPECIES_COL]] %in% fg_tips,
    "Foreground\n(rapid-innovation lineage)",
    "Background"
  )
  d$Lineage <- factor(
    d$Lineage,
    levels = c("Background", "Foreground\n(rapid-innovation lineage)")
  )

  if (length(unique(d$Lineage)) < 2 || any(table(d$Lineage) < 2)) {
    message("Skipping foreground/background boxplot for ", phenotype, ": not enough tips in one group")
    return(invisible(NULL))
  }

  p <- ggplot(d, aes(Lineage, y_plot, fill = Lineage)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = point_alpha, size = 1.3) +
    scale_fill_manual(values = c(
      "Background" = "grey75",
      "Foreground\n(rapid-innovation lineage)" = "#35B779"
    )) +
    # Match phylogenetic_signal_bar text scale; p-value ~3× smaller,
    # bottom-center so it stays inside the panel.
    stat_compare_means(
      method = "wilcox.test",
      size = 1.2,
      label.x.npc = "center",
      label.y.npc = 0.05
    ) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  if (!is.null(out_file)) {
    dir_ok(dirname(out_file))
    # Same canvas as plot_signal_bars so panels match in the final figure.
    ggsave(out_file, p, width = 8, height = 8)
    message("Saved: ", out_file)
  }
  p
}

# Re-draw only the foreground/background boxplot from saved pipeline outputs
# (cleaned_data.tsv + foreground_branches.tsv + the project tree).
#
# Examples:
#   source("phylo_innovation_analysis_new.R")
#   replot_foreground_background_box("TECHINNO2025_ResEff")
#   replot_foreground_background_box(
#     phenotype = "TECHINNO2025_ResEff",
#     trait = "Range.Size",
#     log_transform = TRUE,
#     ylim = c(3, 8)
#   )
replot_foreground_background_box <- function(
  phenotype = PHENOTYPES[[1]],
  results_dir = RESULTS_DIR,
  trait = "Range.Size",
  log_transform = TRUE,
  point_alpha = 0.18,
  ylim = c(3, 8),
  data_file = NULL,
  foreground_file = NULL,
  tree_file = NULL,
  out_file = NULL
) {
  if (is.null(data_file)) {
    data_file <- file.path(results_dir, phenotype, "cleaned_data.tsv")
  }
  if (is.null(foreground_file)) {
    foreground_file <- file.path(results_dir, phenotype, "foreground_branches.tsv")
  }
  if (is.null(tree_file)) {
    tree_file <- file.path(DATA_DIR, TREE_FILE)
  }
  if (is.null(out_file)) {
    out_file <- file.path(
      results_dir, phenotype, "figures", "foreground_background_boxplot.pdf"
    )
  }
  if (!file.exists(data_file)) stop("cleaned_data.tsv not found: ", data_file)
  if (!file.exists(foreground_file)) {
    stop("foreground_branches.tsv not found: ", foreground_file)
  }
  if (!file.exists(tree_file)) stop("tree not found: ", tree_file)

  dat <- read.delim(data_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!SPECIES_COL %in% names(dat)) {
    stop("Species column '", SPECIES_COL, "' missing in ", data_file)
  }
  if (!trait %in% names(dat)) {
    stop("Trait '", trait, "' missing in ", data_file)
  }
  rownames(dat) <- dat[[SPECIES_COL]]

  foreground <- read.delim(
    foreground_file, stringsAsFactors = FALSE, check.names = FALSE
  )
  tree <- read.tree(tree_file)

  plot_foreground_background_box(
    tree = tree,
    dat = dat,
    phenotype = trait,
    foreground = foreground,
    out_file = out_file,
    log_transform = log_transform,
    point_alpha = point_alpha,
    ylim = ylim
  )
}

# 2) Phylogenetic-signal bar plot (Blomberg's K and Pagel's lambda) across
# all SIGNAL_TRAITS, with BH-adjusted significance stars. This is the
# graphical counterpart of "does this trait carry phylogenetic signal".
#
# When max(K) exceeds k_break[2] (e.g. body mass ~27 while other traits
# sit near 0-2), the K panel uses a single compressed y-axis that skips
# [k_break[1], k_break[2]] while keeping equal unit scale (and equal tick
# spacing) on both sides of the break. The K and lambda panels get equal
# overall height.
plot_signal_bars <- function(
  signal,
  out_file = NULL,
  alpha = SIGNAL_ALPHA,
  k_break = c(3, 24)
) {
  df <- signal
  if (!"Trait" %in% names(df) && "trait" %in% names(df)) {
    df$Trait <- df$trait
  }
  df$K_p_adj <- p.adjust(df$K_p, method = "BH")
  df$lambda_p_adj <- p.adjust(df$lambda_p, method = "BH")
  df$K_sig <- df$K_p_adj <= alpha
  df$lambda_sig <- df$lambda_p_adj <= alpha

  trait_order <- labelize(df$Trait[order(df$K, decreasing = TRUE)])
  df$trait_label <- factor(labelize(df$Trait), levels = trait_order)
  x_lvls <- levels(df$trait_label)
  df$x <- as.numeric(df$trait_label)
  bar_w <- 0.72
  # Half of the previous 19.5 so fonts stay readable when panels are
  # composed into the final multi-panel figure.
  base_sz <- 10
  star_sz <- 4.5
  n_tr <- length(x_lvls)
  # Shared continuous x for both panels (same limits → aligned bars).
  x_limits <- c(0.4, n_tr + 0.6)

  shared_x_scale <- function(show_labels = FALSE) {
    scale_x_continuous(
      breaks = seq_len(n_tr),
      labels = if (show_labels) x_lvls else NULL,
      limits = x_limits,
      expand = expansion(mult = 0, add = 0)
    )
  }

  panel_theme <- function(show_x) {
    theme_bw(base_size = base_sz) +
      theme(
        axis.text.x = if (show_x) {
          element_text(angle = 45, hjust = 1)
        } else {
          element_blank()
        },
        axis.ticks.x = if (show_x) element_line() else element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)
      )
  }

  make_plain_panel <- function(y_col, sig_col, ylab, show_x, yintercept = 1) {
    y_max <- max(df[[y_col]], na.rm = TRUE)
    y_top <- max(1.05, y_max * 1.12)
    ggplot(df, aes(x = x, y = .data[[y_col]])) +
      geom_hline(yintercept = yintercept, linetype = "dashed", color = "gray45", linewidth = 0.7) +
      geom_col(fill = "#35B779", width = bar_w) +
      geom_text(
        data = df[df[[sig_col]], , drop = FALSE],
        aes(y = .data[[y_col]] + 0.03 * y_top, label = "*"),
        vjust = 0.5,
        size = star_sz,
        fontface = "bold"
      ) +
      shared_x_scale(show_labels = show_x) +
      scale_y_continuous(
        limits = c(0, y_top),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(x = NULL, y = ylab) +
      panel_theme(show_x)
  }

  # Single-panel broken K axis: compress [y_lo, y_hi] to a small visual gap
  # so 1 data unit has the same pixel size above and below the break.
  make_broken_k_panel <- function(break_at = k_break, gap = 0.45) {
    y_lo <- break_at[1]
    y_hi <- break_at[2]
    # Fixed top so the Mass star fits without excess empty space.
    y_top <- 28

    to_disp <- function(y) {
      ifelse(
        y <= y_lo, y,
        ifelse(y >= y_hi, y_lo + gap + (y - y_hi), NA_real_)
      )
    }

    # Lower bar segments (clipped at y_lo); upper segments for outliers.
    bars_lo <- data.frame(
      xmin = df$x - bar_w / 2,
      xmax = df$x + bar_w / 2,
      ymin = 0,
      ymax = to_disp(pmin(df$K, y_lo)),
      stringsAsFactors = FALSE
    )
    hi <- df[df$K > y_hi, , drop = FALSE]
    bars_hi <- data.frame(
      xmin = hi$x - bar_w / 2,
      xmax = hi$x + bar_w / 2,
      ymin = to_disp(y_hi),
      ymax = to_disp(hi$K),
      stringsAsFactors = FALSE
    )

    # Upper segment starts at y_hi (24), not one tick above it.
    tick_data <- c(seq(0, y_lo, by = 1), seq(y_hi, y_top, by = 1))
    tick_disp <- to_disp(tick_data)
    y_top_disp <- to_disp(y_top)

    # Stars at bar tops (+ small offset in display units).
    sig <- df[df$K_sig, , drop = FALSE]
    stars <- data.frame(
      x = sig$x,
      y = to_disp(sig$K) + 0.35,
      stringsAsFactors = FALSE
    )
    stars <- stars[is.finite(stars$y), , drop = FALSE]

    ggplot() +
      geom_hline(
        yintercept = to_disp(1),
        linetype = "dashed", color = "gray45", linewidth = 0.7
      ) +
      geom_rect(
        data = bars_lo,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "#35B779", color = NA
      ) +
      geom_rect(
        data = bars_hi,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "#35B779", color = NA
      ) +
      # White band over the compressed gap so bars don't fill it.
      annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = y_lo, ymax = y_lo + gap,
        fill = "white", color = NA
      ) +
      geom_text(
        data = stars,
        aes(x = x, y = y, label = "*"),
        vjust = 0.5, size = star_sz, fontface = "bold"
      ) +
      shared_x_scale(show_labels = FALSE) +
      scale_y_continuous(
        breaks = tick_disp,
        labels = tick_data,
        limits = c(0, y_top_disp),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(x = NULL, y = "Blomberg's K") +
      panel_theme(FALSE) +
      theme(plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5))
  }

  use_break <- !is.null(k_break) &&
    length(k_break) == 2 &&
    is.finite(max(df$K, na.rm = TRUE)) &&
    max(df$K, na.rm = TRUE) > k_break[2]

  p_lambda <- make_plain_panel(
    "lambda", "lambda_sig", expression("Pagel's " * lambda),
    show_x = TRUE
  ) +
    theme(plot.margin = margin(t = 2, r = 5.5, b = 5.5, l = 5.5))

  if (use_break) {
    p_k <- make_broken_k_panel()
    fig <- p_k / p_lambda + plot_layout(heights = c(1, 1))
  } else {
    p_k <- make_plain_panel("K", "K_sig", "Blomberg's K", show_x = FALSE)
    fig <- p_k / p_lambda + plot_layout(heights = c(1, 1))
  }

  if (!is.null(out_file)) {
    dir_ok(dirname(out_file))
    ggsave(out_file, fig, width = 8, height = 8)
    message("Saved: ", out_file)
  }
  fig
}

# Re-draw only the phylogenetic-signal bar figure from a saved TSV
# (no need to re-run phylosig / the full phenotype pipeline).

replot_phylogenetic_signal_bars <- function(
  phenotype = PHENOTYPES[[1]],
  results_dir = RESULTS_DIR,
  signal_file = NULL,
  out_file = NULL,
  alpha = SIGNAL_ALPHA,
  k_break = c(3, 24)
) {
  if (is.null(signal_file)) {
    signal_file <- file.path(results_dir, phenotype, "phylogenetic_signal.tsv")
  }
  if (is.null(out_file)) {
    out_file <- file.path(
      results_dir, phenotype, "figures", "phylogenetic_signal_bar.pdf"
    )
  }
  if (!file.exists(signal_file)) {
    stop("phylogenetic_signal.tsv not found: ", signal_file)
  }
  signal <- read.delim(signal_file, stringsAsFactors = FALSE, check.names = FALSE)
  plot_signal_bars(
    signal,
    out_file = out_file,
    alpha = alpha,
    k_break = k_break
  )
}

# ---- per-phenotype pipeline ----
run_phenotype_analysis <- function(tree, dat, phenotype, prep_info) {
  out_dir <- file.path(RESULTS_DIR, phenotype)
  fig_dir <- file.path(out_dir, "figures")
  dir_ok(out_dir)
  dir_ok(fig_dir)

  message("=== ", phenotype, " ===")

  # Step 1
  write_tsv(dat, file.path(out_dir, "cleaned_data.tsv"))
  write_tsv(prep_info$report, file.path(out_dir, "data_match_report.tsv"))
  write_tsv(
    data.frame(species = prep_info$missing_in_table, stringsAsFactors = FALSE),
    file.path(out_dir, "missing_in_table.tsv")
  )
  write_tsv(
    data.frame(species = prep_info$missing_in_tree, stringsAsFactors = FALSE),
    file.path(out_dir, "missing_in_tree.tsv")
  )
  cat(
    phenotype, ": tree =", prep_info$report$value[1],
    ", table =", prep_info$report$value[2],
    ", matched =", prep_info$report$value[3], "\n"
  )

  x <- trait_vector(dat, tree, phenotype)

  # Step 2: phylogenetic signal + bar-plot figure
  signal <- compute_phylogenetic_signal(tree, dat, SIGNAL_TRAITS)
  write_tsv(signal, file.path(out_dir, "phylogenetic_signal.tsv"))
  signal_bars <- plot_signal_bars(signal, file.path(fig_dir, "phylogenetic_signal_bar.pdf"))

  # Step 3
  anc_df <- compute_ancestral_states(tree, x)
  write_tsv(anc_df, file.path(out_dir, "ancestral_states.tsv"))
  anc_named <- setNames(anc_df$state, anc_df$node)

  # Step 4
  branch_changes <- compute_branch_changes(tree, x, anc_named)
  write_tsv(branch_changes, file.path(out_dir, "branch_changes.tsv"))

  # Step 5: foreground branches + foreground/background lineage boxplot
  fg_info <- identify_foreground(branch_changes, FOREGROUND_SD)
  write_tsv(fg_info$foreground, file.path(out_dir, "foreground_branches.tsv"))
  top20 <- branch_changes |>
    arrange(desc(delta)) |>
    slice_head(n = min(20, nrow(branch_changes)))
  write_tsv(top20, file.path(out_dir, "top20_branch_changes.tsv"))

  # Step 6
  pgls_uni <- run_univariate_pgls(tree, dat, phenotype, COVARIATES)
  write_tsv(pgls_uni, file.path(out_dir, "pgls_univariate.tsv"))

  # Step 7
  dat_tip <- dat[tree$tip.label, c(phenotype, MULTIVARIATE_COVARIATES), drop = FALSE]
  cc <- complete.cases(dat_tip)
  drop_tips_mv <- tree$tip.label[!cc]
  tree_mv <- if (length(drop_tips_mv) > 0) drop.tip(tree, drop_tips_mv) else tree
  dat_mv <- dat[tree_mv$tip.label, , drop = FALSE]
  pgls_mv <- run_multivariate_pgls(tree_mv, dat_mv, phenotype, MULTIVARIATE_COVARIATES)
  pgls_mv_out <- transform(
    pgls_mv,
    model_R2 = attr(pgls_mv, "R2")
  )
  write_tsv(pgls_mv_out, file.path(out_dir, "pgls_multivariate.tsv"))

  # Step 8
  pca_res <- run_pca(dat, PCA_TRAITS, phenotype)
  write_tsv(pca_res$coords, file.path(out_dir, "pca_coordinates.tsv"))

  # Step 9
  save_genome_inputs(dat, phenotype, out_dir)

  # Figures
  plot_contmap(tree, x, file.path(fig_dir, "contMap.pdf"))
  plot_ancestral_tree(tree, anc_df, file.path(fig_dir, "ancestral_states.pdf"))
  plot_foreground_tree(
    tree, branch_changes, fg_info$foreground,
    file.path(fig_dir, "foreground_branches.pdf")
  )
  plot_delta_histogram(
    branch_changes, fg_info$threshold,
    file.path(fig_dir, "delta_histogram.pdf")
  )
  plot_parent_child_scatter(
    branch_changes, fg_info$foreground,
    file.path(fig_dir, "parent_child_scatter.pdf")
  )
  plot_forest(pgls_uni, file.path(fig_dir, "forest_plot.pdf"))
  plot_correlation_matrix(
    dat, CORRELATION_TRAITS,
    file.path(fig_dir, "correlation_matrix.pdf")
  )
  plot_pca(pca_res$coords, file.path(fig_dir, "pca.pdf"))
  
  phenotype <- "Range.Size"

  foreground_background_box <- plot_foreground_background_box(
    tree, dat, phenotype, fg_info$foreground,
    file.path(fig_dir, "foreground_background_boxplot.pdf"),
    log_transform = TRUE,
    point_alpha = 0.18,
    ylim = c(3, 8)
  )

  invisible(out_dir)
  
  return(
    list(
      signal_bars = signal_bars, 
      foreground_background_box = foreground_background_box
    )
  )
  
}

# ---- main ----
main <- function() {
  inputs <- load_inputs()
  prep <- prepare_data(inputs$tree, inputs$dat)
  
  for (ph in PHENOTYPES) {
    plots <- run_phenotype_analysis(prep$tree, prep$dat, ph, prep)
  }

  message("Done. Results in: ", RESULTS_DIR)
  return(plots)
}

if (sys.nframe() == 0) {
  plots <- main()
  print(plots$signal_bars)
  print(plots$foreground_background_box)
}

signal_bars <- replot_phylogenetic_signal_bars("TECHINNO2025_ResEff")
foreground_background_box <- replot_foreground_background_box("TECHINNO2025_ResEff")

signal_bars

