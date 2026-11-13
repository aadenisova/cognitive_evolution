# phase_a_transform.R
# Phase A: Transform comparison for RERconverge continuous analysis
# Tests 4 transforms (sqrt, log, asinh, none), calibrates parametric vs permutation

suppressPackageStartupMessages({
  library(RERconverge)
  library(parallel)
})
source("/workspace/parallel_core.R")

R_LIBS_USER <- "/workspace/.Rlib"
NCORES <- 16
BATCH.SIZE <- 200
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_a_transform"

# ============================================================
# A1: Setup — load phenotype, read trees
# ============================================================

message("=== A1: Loading phenotype and reading trees ===")

# Load phenotype
pheno_df <- read.csv("/mnt/user-uploads/df_merged_for_rer_no_na.csv", row.names = 1)
tip.vals <- pheno_df$TECHINNO2025_ResEff
names(tip.vals) <- rownames(pheno_df)
message(sprintf("Phenotype loaded: %d species, range [%.4f, %.4f]",
                length(tip.vals), min(tip.vals), max(tip.vals)))

# Read trees (this takes ~7 min for 272,811 trees)
message("Reading trees (this takes ~7 min)...")
treesObj <- readTrees("/workspace/rer_trees_birds.txt", max.read = 272811)
message(sprintf("Trees loaded: %d CNEEs, %d species", treesObj$numTrees, length(treesObj$masterTree$tip.label)))

# Save treesObj RDS to workspace
saveRDS(treesObj, "/workspace/treesObj_cont.rds")
message("treesObj saved to /workspace/treesObj_cont.rds")

# Define useSpecies (all species in phenotype)
useSpecies <- names(tip.vals)
message(sprintf("useSpecies: %d species", length(useSpecies)))

# ============================================================
# A2: Compute RER + correlation for each transform
# ============================================================

transforms <- c("sqrt", "log", "asinh", "none")
res_list <- list()
rermat_list <- list()

for (trans in transforms) {
  message(sprintf("\n=== A2: Processing transform = %s ===", trans))

  # Reload treesObj from RDS (transformPaths modifies in place)
  treesObj_fresh <- readRDS("/workspace/treesObj_cont.rds")

  # Compute RER matrix using parallel version
  message(sprintf("Computing RERmat (transform=%s, %d cores)...", trans, NCORES))
  RERmat <- parallel_getAllResiduals(treesObj_fresh,
                                      transform = trans,
                                      useSpecies = useSpecies,
                                      use.weights = TRUE,
                                      norm = "scale",
                                      min.sp = 10,
                                      ncores = NCORES,
                                      batch.size = BATCH.SIZE)
  message(sprintf("RERmat: %d CNEEs x %d branches", nrow(RERmat), ncol(RERmat)))

  # Save RERmat
  saveRDS(RERmat, sprintf("/workspace/RERmat_%s.rds", trans))
  rermat_list[[trans]] <- RERmat

  # Compute correlation with phenotype
  message("Computing correlation with phenotype...")
  treesObj_fresh <- readRDS("/workspace/treesObj_cont.rds")  # reload for char2Paths
  charP <- char2Paths(tip.vals, treesObj_fresh)
  res <- correlateWithContinuousPhenotype(RERmat, charP,
                                           min.sp = 10,
                                           winsorizeRER = 3,
                                           winsorizetrait = 3)
  res$p.adj <- p.adjust(res$P, method = "BH")

  message(sprintf("Results: %d valid CNEEs, P<0.05: %d, p.adj<0.05: %d",
                  sum(!is.na(res$P)),
                  sum(res$P < 0.05, na.rm = TRUE),
                  sum(res$p.adj < 0.05, na.rm = TRUE)))

  saveRDS(res, sprintf("/workspace/res_%s.rds", trans))
  res_list[[trans]] <- res

  # Free memory
  rm(treesObj_fresh, RERmat, charP)
  gc()
}

# ============================================================
# A3: Calibration check (100 perms x top 5,000 x 4 transforms)
# ============================================================

message("\n=== A3: Calibration check ===")

NPERMS_CAL <- 100
TOPN_CAL <- 5000

calibration_results <- list()

for (trans in transforms) {
  message(sprintf("\n--- Calibration for transform = %s ---", trans))

  res <- res_list[[trans]]
  RERmat <- readRDS(sprintf("/workspace/RERmat_%s.rds", trans))

  # Select top 5,000 CNEEs by parametric P
  valid_p <- res$P[!is.na(res$P)]
  top_threshold <- sort(valid_p)[min(TOPN_CAL, length(valid_p))]
  top_cnees <- rownames(res)[which(res$P <= top_threshold)]
  message(sprintf("Top %d CNEEs selected (P <= %.6e)", length(top_cnees), top_threshold))

  # Subset RERmat to top CNEEs
  RERmat_top <- RERmat[top_cnees, , drop = FALSE]

  # Load treesObj for permulations
  treesObj_fresh <- readRDS("/workspace/treesObj_cont.rds")
  masterTree <- treesObj_fresh$masterTree

  # Observed statistics
  obs_stat <- sign(res$Rho[top_cnees]) * -log10(res$P[top_cnees])
  obs_abs <- abs(obs_stat)

  # Run permulations in parallel
  message(sprintf("Running %d permulations on %d CNEEs (%d cores)...",
                  NPERMS_CAL, length(top_cnees), NCORES))

  runPerm <- function(perm_idx, RERmat_top, tip.vals, masterTree, treesObj_fresh) {
    vec <- simpermvec(tip.vals, masterTree)
    paths <- char2Paths(na.omit(vec), treesObj_fresh)
    out <- getAllCor(RERmat_top, paths, method = "p", min.pos = 0,
                     winsorizeRER = 3, winsorizetrait = 3)
    perm_stat <- sign(out$Rho) * -log10(out$P)
    return(abs(perm_stat))
  }

  perm_results <- mclapply(1:NPERMS_CAL, function(i) {
    runPerm(i, RERmat_top, tip.vals, masterTree, treesObj_fresh)
  }, mc.cores = NCORES)

  # Compute permpval for each CNEE
  perm_matrix <- do.call(rbind, perm_results)  # NPERMS x NCNEEs
  counters <- colSums(perm_matrix >= matrix(obs_abs, nrow = NPERMS_CAL, ncol = length(obs_abs), byrow = TRUE))
  permpval <- (counters + 1) / (NPERMS_CAL + 1)

  # Calibration metrics
  spearman_cal <- cor(-log10(res$P[top_cnees]), -log10(permpval), method = "spearman", use = "complete.obs")

  # Lambda (genomic inflation factor)
  chi2 <- qchisq(res$P[top_cnees], 1, lower.tail = FALSE)
  # CAVEAT: lambda is computed on the top-5,000 CNEEs ranked by P, so it is large by construction;
  # it is NOT a genome-wide inflation factor (see 02b_genomewide_lambda.R). Same for prop_p_gt_0.5.
  lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)

  # P-value distribution
  pvals <- res$P[top_cnees]
  p_uniform <- mean(pvals > 0.5, na.rm = TRUE)  # proportion with P > 0.5 (should be ~0.5 under null)

  message(sprintf("Spearman calibration: %.4f", spearman_cal))
  message(sprintf("Lambda (inflation factor): %.4f", lambda))
  message(sprintf("Prop P > 0.5: %.4f (expected ~0.5 under null)", p_uniform))

  calibration_results[[trans]] <- data.frame(
    transform = trans,
    n_cnees = length(top_cnees),
    spearman_cal = spearman_cal,
    lambda = lambda,
    prop_p_gt_0.5 = p_uniform,
    min_permpval = min(permpval),
    n_permpval_0.05 = sum(permpval < 0.05),
    stringsAsFactors = FALSE
  )

  # Save permpval for this transform
  res$permpval_cal <- NA
  res$permpval_cal[top_cnees] <- permpval
  saveRDS(res, sprintf("/workspace/res_%s.rds", trans))

  rm(treesObj_fresh, RERmat, RERmat_top, perm_results, perm_matrix)
  gc()
}

# ============================================================
# A4: Select best transform
# ============================================================

message("\n=== A4: Selecting best transform ===")

cal_df <- do.call(rbind, calibration_results)
cal_df$lambda_dist <- abs(cal_df$lambda - 1.0)  # distance from ideal lambda=1
cal_df <- cal_df[order(-cal_df$spearman_cal, cal_df$lambda_dist), ]

print(cal_df)

best_transform <- cal_df$transform[1]
message(sprintf("\nBest transform: %s (Spearman=%.4f, Lambda=%.4f)",
                best_transform, cal_df$spearman_cal[1], cal_df$lambda[1]))

# Save comparison table
write.csv(cal_df, file.path(OUTPUT_DIR, "transform_comparison.csv"), row.names = FALSE)
message(sprintf("Saved: %s/transform_comparison.csv", OUTPUT_DIR))

# Save best transform name
writeLines(best_transform, file.path(OUTPUT_DIR, "best_transform.txt"))
message(sprintf("Saved: %s/best_transform.txt", OUTPUT_DIR))

message("\n=== Phase A complete ===")
