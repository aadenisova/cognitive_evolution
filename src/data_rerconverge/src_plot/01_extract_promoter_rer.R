#!/usr/bin/env Rscript
# 01_extract_promoter_rer.R
# Cluster-only data extraction (not figure plotting).
# Extract per-species RER values for top 40 promoter CNEEs (20 pos + 20 neg).
# Ensures >=5 per direction within 1kb of TSS.
# Strategy: subset treesObj to 40 CNEEs, manually log-transform, coreGetResiduals.
# Aggregate branch-level RER to species-level (mean per species per CNEE).
# Exports long-format CSV for plotting + selection key.
#
# Local copy already in the project: promoter_rer_data_40.csv
# Article figure: Rscript src_plot/figure_article_panel.R

suppressPackageStartupMessages({
  library(RERconverge)
  library(ape)
})

BASE       <- "/mnt/results/rerconverge_cont_v2"
PHENO_PATH <- "/mnt/user-uploads/df_merged_for_rer_no_na.csv"
TREES_PATH <- "/workspace/treesObj_cont.rds"
OUTPUT_DIR <- "/mnt/results/plotting/data"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──────────────────────────────────────────────────────────────
message("Loading robust CNEEs...")
robust <- read.csv(file.path(BASE, "robust_cnees_annotated_v3.csv"), stringsAsFactors = FALSE)
sym_map <- read.csv(file.path(BASE, "gene_symbol_map.csv"), stringsAsFactors = FALSE)
id2sym <- setNames(sym_map$symbol, as.character(sym_map$geneId))

message("Loading phenotype...")
pheno <- read.csv(PHENO_PATH, stringsAsFactors = FALSE, row.names = 1)

message("Loading trees object (3.3 GB)...")
treesObj <- readRDS(TREES_PATH)

# ── Select 40 promoter CNEEs ────────────────────────────────────────────────
prom <- robust[robust$category == "Promoter", ]
prom$abs_dist <- abs(prom$distanceToTSS)

select_top20 <- function(df, dir) {
  d <- df[df$direction == dir, ]
  d <- d[order(d$permpval), ]
  top20 <- head(d, 20)
  n_within_1kb <- sum(top20$abs_dist <= 1000)
  if (n_within_1kb < 5) {
    guaranteed <- head(d, 15)
    remaining <- d[!d$CNEE %in% guaranteed$CNEE & d$abs_dist <= 1000, ]
    remaining <- remaining[order(remaining$permpval), ]
    fill <- head(remaining, 5)
    top20 <- rbind(guaranteed, fill)
    top20 <- top20[order(top20$permpval), ]
    n_within_1kb <- sum(top20$abs_dist <= 1000)
  }
  message(sprintf("  %s: selected %d CNEEs, %d within 1kb of TSS",
                  dir, nrow(top20), n_within_1kb))
  top20
}

message("Selecting top 20 positive + 20 negative promoter CNEEs...")
pos20 <- select_top20(prom, "positive")
neg20 <- select_top20(prom, "negative")
all40 <- rbind(pos20, neg20)

all40$geneSymbol <- sapply(as.character(all40$geneId), function(x) {
  if (is.na(x) || x == "") return("")
  v <- id2sym[as.character(as.integer(x))]
  if (is.na(v)) "" else v
})

# ── Subset treesObj to 40 CNEEs ─────────────────────────────────────────────
all_cnee_names <- rownames(treesObj$report)
cnee_names <- all40$CNEE
row_idx <- match(cnee_names, all_cnee_names)
if (any(is.na(row_idx)))
  stop(sprintf("CNEE not found: %s", paste(cnee_names[is.na(row_idx)], collapse=", ")))

# Pre-compute log transform offset from FULL paths (needs 500-row sample)
message("Pre-computing log transform offset...")
set.seed(1)
iis <- sample(nrow(treesObj$paths), 500)
tmpdat <- treesObj$paths[iis, ]
mval <- min(tmpdat[tmpdat > 0], na.rm = TRUE)
message(sprintf("  Log offset (mval): %.6f", mval))

message("Subsetting treesObj to 40 CNEEs...")
treesObj_sub <- treesObj
treesObj_sub$trees    <- treesObj$trees[row_idx]
treesObj_sub$numTrees <- length(row_idx)
treesObj_sub$report   <- treesObj$report[row_idx, , drop = FALSE]
treesObj_sub$paths    <- treesObj$paths[row_idx, , drop = FALSE]
# Apply log transform manually (avoid transformPaths which samples 500 rows)
treesObj_sub$paths <- log(treesObj_sub$paths + mval)
treesObj_sub$transform <- "log"
treesObj_sub$transformFunc <- function(x) log(x + mval)
treesObj_sub$transformInv <- function(x) exp(x) - mval
treesObj_sub$logOffset <- mval
treesObj_sub$pathsImputed <- NULL
treesObj_sub$weights <- treesObj_sub$paths
treesObj_sub$weights[] <- 1

rm(treesObj)
gc()

# ── Compute RER residuals ───────────────────────────────────────────────────
message("Computing RER residuals (log transform, 40 CNEEs, no imputation)...")
resids <- coreGetResiduals(treesObj_sub, n.pcs = 0, use.weights = FALSE)
r <- getRMat(resids, all = FALSE, use.rows = 1:length(row_idx), norm = "scale")
message(sprintf("  RER matrix: %d CNEEs x %d branches", nrow(r), ncol(r)))

# ── Aggregate branch-level RER to species-level ─────────────────────────────
# Each column in r is named by species accession (with duplicates for multiple branches)
branch_species <- colnames(r)
unique_species <- unique(branch_species)
# Keep only species present in phenotype
unique_species <- intersect(unique_species, rownames(pheno))
message(sprintf("  Branches: %d, Unique species (in pheno): %d", length(branch_species), length(unique_species)))

# For each CNEE, average RER across branches per species
species_rer <- matrix(NA_real_, nrow = length(cnee_names), ncol = length(unique_species),
                      dimnames = list(cnee_names, unique_species))
for (i in seq_along(cnee_names)) {
  for (j in seq_along(unique_species)) {
    sp <- unique_species[j]
    vals <- r[i, branch_species == sp]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      species_rer[i, j] <- mean(vals)
    }
  }
}
message(sprintf("  Species-level RER: %d CNEEs x %d species", nrow(species_rer), ncol(species_rer)))

# ── Build long-format data frame ────────────────────────────────────────────
message("Building long-format data...")
accessions <- colnames(species_rer)
pheno_match <- pheno[accessions, ]

rows <- list()
for (i in seq_along(cnee_names)) {
  cnee <- cnee_names[i]
  info <- all40[all40$CNEE == cnee, ]
  rer_vals <- as.numeric(species_rer[i, ])
  for (j in seq_along(accessions)) {
    acc <- accessions[j]
    rows[[length(rows) + 1]] <- data.frame(
      CNEE              = cnee,
      species_accession = acc,
      scientific_name   = pheno_match[acc, "ScientificName"],
      RER               = rer_vals[j],
      TECHINNO2025_ResEff = pheno_match[acc, "TECHINNO2025_ResEff"],
      Rho               = info$Rho,
      permpval          = info$permpval,
      P                 = info$P,
      direction         = info$direction,
      geneId            = info$geneId,
      geneSymbol        = info$geneSymbol,
      distanceToTSS     = info$distanceToTSS,
      abs_distance_to_TSS = info$abs_dist,
      stringsAsFactors  = FALSE
    )
  }
}
long_df <- do.call(rbind, rows)

# Remove rows with NA RER or NA phenotype
n_before <- nrow(long_df)
long_df <- long_df[!is.na(long_df$RER) & !is.na(long_df$TECHINNO2025_ResEff), ]
n_after <- nrow(long_df)
message(sprintf("  Rows: %d total, %d after NA removal (%d removed)",
                n_before, n_after, n_before - n_after))

# ── Save 40-CNEE data ───────────────────────────────────────────────────────
out40 <- file.path(OUTPUT_DIR, "promoter_rer_data_40.csv")
write.csv(long_df, out40, row.names = FALSE)
message(sprintf("Saved: %s (%d rows)", out40, nrow(long_df)))

# ── Select 20 for the plot (10 pos + 10 neg, by |Rho| desc) ──────────────────
pos10 <- pos20[order(abs(pos20$Rho), decreasing = TRUE), ][1:10, ]
neg10 <- neg20[order(abs(neg20$Rho), decreasing = TRUE), ][1:10, ]
plot20_cnees <- c(neg10$CNEE, pos10$CNEE)

long20 <- long_df[long_df$CNEE %in% plot20_cnees, ]
facet_order <- data.frame(CNEE = plot20_cnees, facet_order = seq_along(plot20_cnees))
long20 <- merge(long20, facet_order, by = "CNEE", sort = FALSE)
long20 <- long20[order(long20$facet_order, long20$species_accession), ]

out20 <- file.path(OUTPUT_DIR, "promoter_rer_data_20.csv")
write.csv(long20, out20, row.names = FALSE)
message(sprintf("Saved: %s (%d rows)", out20, nrow(long20)))

# ── Save selection key ──────────────────────────────────────────────────────
selection <- data.frame(
  CNEE       = all40$CNEE,
  direction  = all40$direction,
  rank_by_permpval = ave(all40$permpval, all40$direction, FUN = function(x) rank(x, ties.method = "first")),
  Rho        = all40$Rho,
  permpval   = all40$permpval,
  geneSymbol = all40$geneSymbol,
  geneId     = all40$geneId,
  distanceToTSS   = all40$distanceToTSS,
  abs_distance_to_TSS = all40$abs_dist,
  within_1kb = all40$abs_dist <= 1000,
  in_plot    = all40$CNEE %in% plot20_cnees,
  stringsAsFactors = FALSE
)
selection <- selection[order(selection$direction, selection$rank_by_permpval), ]

out_sel <- file.path(OUTPUT_DIR, "promoter_cnee_selection.csv")
write.csv(selection, out_sel, row.names = FALSE)
message(sprintf("Saved: %s (%d CNEEs)", out_sel, nrow(selection)))

# ── Summary ─────────────────────────────────────────────────────────────────
message("\n=== Summary ===")
message(sprintf("Total promoter CNEEs in robust: %d", nrow(prom)))
message(sprintf("Selected 40: %d pos + %d neg", sum(all40$direction == "positive"), sum(all40$direction == "negative")))
message(sprintf("Within 1kb: pos=%d, neg=%d",
                sum(all40$direction == "positive" & all40$abs_dist <= 1000),
                sum(all40$direction == "negative" & all40$abs_dist <= 1000)))
message(sprintf("In plot (20): pos=%d, neg=%d",
                sum(selection$in_plot & selection$direction == "positive"),
                sum(selection$in_plot & selection$direction == "negative")))
message(sprintf("Species per CNEE (after NA removal): median %.0f", median(table(long_df$CNEE))))

message("\n=== Step 01 complete ===")
