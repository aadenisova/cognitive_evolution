#!/usr/bin/env Rscript
# 07_gene_symbol_map.R
# Create comprehensive geneId -> symbol mapping using org.Gg.eg.db.
# Covers ALL geneIds appearing in robust_cnees, not just attributed ones.

suppressPackageStartupMessages({
  library(org.Gg.eg.db)
  library(AnnotationDbi)
})

BASE <- "/mnt/results/rerconverge_cont_v2"
OUTPUT <- file.path(BASE, "gene_symbol_map.csv")

# ============================================================
# Collect all unique geneIds from robust_cnees
# ============================================================
message("=== Loading robust CNEEs for geneId collection ===")
robust <- read.csv(file.path(BASE, "robust_cnees_annotated_v3.csv"), stringsAsFactors = FALSE)
all_ids <- unique(robust$geneId)
all_ids <- all_ids[!is.na(all_ids)]
all_ids <- as.character(as.integer(all_ids))  # remove trailing .0
message(sprintf("Unique geneIds in robust CNEEs: %d", length(all_ids)))

# Also include geneIds from gene_level_stats
gene_stats <- read.csv(file.path(BASE, "gene_level_stats_v3.csv"), stringsAsFactors = FALSE)
gs_ids <- as.character(gene_stats$geneId)
all_ids <- unique(c(all_ids, gs_ids))
message(sprintf("Total unique geneIds to map: %d", length(all_ids)))

# ============================================================
# Map via org.Gg.eg.db
# ============================================================
message("=== Mapping via org.Gg.eg.db ===")
mapped <- AnnotationDbi::select(org.Gg.eg.db,
                                keys = all_ids,
                                columns = c("SYMBOL", "GENENAME"),
                                keytype = "ENTREZID")

colnames(mapped) <- c("geneId", "symbol", "genename")
mapped$geneId <- as.integer(mapped$geneId)

# Check unmapped
unmapped <- all_ids[!all_ids %in% mapped$geneId[!is.na(mapped$symbol)]]
message(sprintf("Mapped: %d / %d, Unmapped: %d",
                sum(!is.na(mapped$symbol)), length(all_ids), length(unmapped)))
if (length(unmapped) > 0) {
  message("Unmapped IDs: ", paste(head(unmapped, 20), collapse=", "))
}

# ============================================================
# Save
# ============================================================
write.csv(mapped, OUTPUT, row.names = FALSE)
message(sprintf("\nSaved: %s (%d rows)", OUTPUT, nrow(mapped)))

# Quick stats
n_with_sym <- sum(!is.na(mapped$symbol))
message(sprintf("Rows with symbol: %d / %d (%.1f%%)", n_with_sym, nrow(mapped), 100*n_with_sym/nrow(mapped)))

message("\n=== Step 07 complete ===")
