#!/usr/bin/env Rscript
# 02b_genomewide_lambda.R
# Genome-wide genomic-inflation factor (lambda) of the PARAMETRIC p-values for each path transform,
# computed on ALL CNEEs with a valid P (not only the top 5,000), plus agreement of Rho between transforms.
# Input : /workspace/res_<transform>.rds   (written by 01_transform_selection.R)
# Output: transform_genomewide_lambda.csv
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_a_transform"
transforms <- c("sqrt", "log", "asinh", "none")

res <- lapply(transforms, function(t) readRDS(sprintf("/workspace/res_%s.rds", t)))
names(res) <- transforms

lambda_of <- function(P) {
  P <- P[!is.na(P) & P > 0 & P <= 1]
  median(qchisq(P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

out <- do.call(rbind, lapply(transforms, function(t) {
  r <- res[[t]]
  data.frame(
    transform         = t,
    n_valid           = sum(!is.na(r$P)),
    lambda_genomewide = lambda_of(r$P),
    n_P_lt_0.05       = sum(r$P < 0.05, na.rm = TRUE),
    spearman_rho_vs_log = suppressWarnings(cor(r$Rho, res[["log"]]$Rho, method = "spearman",
                                               use = "pairwise.complete.obs")),
    stringsAsFactors = FALSE
  )
}))
print(out)
write.csv(out, file.path(OUTPUT_DIR, "transform_genomewide_lambda.csv"), row.names = FALSE)
