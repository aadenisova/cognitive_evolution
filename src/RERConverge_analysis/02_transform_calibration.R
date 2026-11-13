# phase_a_calibration.R
# Re-run just the calibration step (transforms already computed and saved)
# Fixes: (1) root tree for simpermvec, (2) use df[rows, "col"] indexing

suppressPackageStartupMessages({
  library(RERconverge)
  library(ape)
  library(parallel)
})

OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_a_transform"
ROOT_SPECIES <- "GCF_028389875"  # Rhea pennata

# Load phenotype
pheno_df <- read.csv("/mnt/user-uploads/df_merged_for_rer_no_na.csv", row.names = 1)
tip.vals <- pheno_df$TECHINNO2025_ResEff
names(tip.vals) <- rownames(pheno_df)

# Load treesObj and root the master tree for simpermvec
treesObj <- readRDS("/workspace/treesObj_cont.rds")
masterTree <- treesObj$masterTree
masterTree_rooted <- root(masterTree, outgroup = ROOT_SPECIES, resolve.root = TRUE)
message(sprintf("Tree rooted: %s, binary: %s", is.rooted(masterTree_rooted), is.binary(masterTree_rooted)))

transforms <- c("sqrt", "log", "asinh", "none")
NPERMS_CAL <- 100
TOPN_CAL <- 5000

calibration_results <- list()

for (trans in transforms) {
  message(sprintf("\n--- Calibration: transform = %s ---", trans))
  t0 <- Sys.time()

  res <- readRDS(sprintf("/workspace/res_%s.rds", trans))
  RERmat <- readRDS(sprintf("/workspace/RERmat_%s.rds", trans))

  # Select top 5,000 by parametric P (among valid P values)
  # Use positional indexing, NOT name indexing on $P (which lacks names attr)
  valid_idx <- which(!is.na(res$P))
  top_idx <- valid_idx[order(res$P[valid_idx])[1:min(TOPN_CAL, length(valid_idx))]]
  top_cnees <- rownames(res)[top_idx]
  message(sprintf("Top %d CNEEs (P range [%.6e, %.6e])", length(top_cnees),
                  min(res$P[top_idx]), max(res$P[top_idx])))

  RERmat_top <- RERmat[top_cnees, , drop = FALSE]

  # Observed statistics for top CNEEs — use positional indexing
  obs_rho <- res$Rho[top_idx]
  obs_p <- res$P[top_idx]
  obs_stat <- sign(obs_rho) * (-log10(obs_p))
  obs_abs <- abs(obs_stat)
  names(obs_abs) <- top_cnees

  message(sprintf("Running %d perms (sequential)...", NPERMS_CAL))

  perm_abs_list <- vector("list", NPERMS_CAL)

  for (i in 1:NPERMS_CAL) {
    # Use rooted tree for simpermvec
    vec <- simpermvec(tip.vals, masterTree_rooted)
    # Use original treesObj for char2Paths (works with unrooted)
    paths <- char2Paths(na.omit(vec), treesObj)
    out <- getAllCor(RERmat_top, paths, method = "p", min.pos = 0,
                     winsorizeRER = 3, winsorizetrait = 3)

    # Align perm results to top_cnees by name
    perm_stat <- rep(0, length(top_cnees))
    names(perm_stat) <- top_cnees
    common <- intersect(top_cnees, rownames(out))
    perm_stat[common] <- sign(out$Rho[common]) * (-log10(out$P[common]))
    perm_stat[is.na(perm_stat)] <- 0

    perm_abs_list[[i]] <- abs(perm_stat)

    if (i %% 20 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      message(sprintf("  perm %d/%d (%.1f sec)", i, NPERMS_CAL, elapsed))
    }
  }

  # Build perm matrix: NPERMS_CAL x length(top_cnees)
  perm_matrix <- do.call(rbind, perm_abs_list)
  message(sprintf("perm_matrix: %d x %d", nrow(perm_matrix), ncol(perm_matrix)))

  # Permpval: for each CNEE, count perms where |perm_stat| >= |obs_stat|
  obs_abs_vec <- as.numeric(obs_abs)
  counters <- integer(length(top_cnees))
  for (i in 1:NPERMS_CAL) {
    counters <- counters + as.integer(perm_matrix[i, ] >= obs_abs_vec)
  }
  permpval <- (counters + 1) / (NPERMS_CAL + 1)
  names(permpval) <- top_cnees

  # Calibration metrics
  spearman_cal <- cor(-log10(obs_p), -log10(permpval), method = "spearman", use = "complete.obs")
  chi2 <- qchisq(1 - obs_p, 1)
  lambda <- median(chi2, na.rm = TRUE) / 0.455
  p_uniform <- mean(obs_p > 0.5, na.rm = TRUE)

  message(sprintf("Spearman: %.4f, Lambda: %.4f, Prop P>0.5: %.4f (%.1f sec)",
                  spearman_cal, lambda, p_uniform,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  calibration_results[[trans]] <- data.frame(
    transform = trans, n_cnees = length(top_cnees),
    spearman_cal = spearman_cal, lambda = lambda,
    prop_p_gt_0.5 = p_uniform, min_permpval = min(permpval),
    n_permpval_0.05 = sum(permpval < 0.05), stringsAsFactors = FALSE)

  rm(RERmat, RERmat_top, perm_abs_list, perm_matrix)
  gc()
}

# Select best transform
message("\n=== Selecting best transform ===")
cal_df <- do.call(rbind, calibration_results)
cal_df$lambda_dist <- abs(cal_df$lambda - 1.0)
cal_df <- cal_df[order(-cal_df$spearman_cal, cal_df$lambda_dist), ]
print(cal_df)

best_transform <- cal_df$transform[1]
message(sprintf("\nBest transform: %s (Spearman=%.4f, Lambda=%.4f)",
                best_transform, cal_df$spearman_cal[1], cal_df$lambda[1]))

write.csv(cal_df, file.path(OUTPUT_DIR, "transform_comparison.csv"), row.names = FALSE)
writeLines(best_transform, file.path(OUTPUT_DIR, "best_transform.txt"))
message("Saved transform_comparison.csv and best_transform.txt")
message("\n=== Phase A calibration complete ===")
