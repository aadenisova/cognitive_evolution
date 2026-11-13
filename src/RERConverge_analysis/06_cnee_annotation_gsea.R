# phase_d_gsea.R
# Phase D: Gene-level GSEA enrichment using full ranked list

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
  library(GenomicRanges)
  library(rtracklayer)
  library(org.Gg.eg.db)
  library(TxDb.Ggallus.UCSC.galGal6.refGene)
  library(homologene)
  library(ggplot2)
})

R_LIBS_USER <- "/workspace/.Rlib"
OUTPUT_DIR <- "/mnt/results/rerconverge_cont_v2/phase_d_gsea"

# ============================================================
# D1: Annotate ALL CNEEs with ChIPseeker
# ============================================================

message("=== D1: Annotating CNEEs with ChIPseeker ===")

cnee_bed <- read.table("/mnt/user-uploads/lifted_annotations_genomic.bed",
                        header = FALSE, stringsAsFactors = FALSE)
colnames(cnee_bed) <- c("chr", "start", "end", "name")
message(sprintf("CNEE BED: %d entries", nrow(cnee_bed)))

cnee_gr <- makeGRangesFromDataFrame(cnee_bed, keep.extra.columns = TRUE,
                                     starts.in.df.are.0based = TRUE)

txdb <- TxDb.Ggallus.UCSC.galGal6.refGene
message("Running ChIPseeker annotatePeak...")
peak_anno <- annotatePeak(cnee_gr, TxDb = txdb, level = "gene",
                           tssRegion = c(-3000, 3000), verbose = FALSE)

anno_df <- as.data.frame(peak_anno)
message(sprintf("Annotation complete: %d CNEEs", nrow(anno_df)))
message(sprintf("Feature distribution:\n%s",
                paste(capture.output(table(anno_df$annotation)), collapse="\n")))

write.csv(anno_df, file.path(OUTPUT_DIR, "cnee_annotation_cont_full.csv"), row.names = FALSE)
message(sprintf("Saved: %s/cnee_annotation_cont_full.csv", OUTPUT_DIR))

# ============================================================
# D2: Aggregate CNEE to gene level
# ============================================================

message("=== D2: Aggregating CNEE to gene level ===")

results_df <- readRDS("/workspace/res_final_cont_v2.rds")

# Merge annotation with results
merged <- merge(anno_df, results_df, by.x = "name", by.y = "CNEE", all.x = TRUE)
merged_valid <- merged[!is.na(merged$permpval), ]
message(sprintf("Merged: %d CNEEs with valid results", nrow(merged_valid)))

# For each gene: find CNEE with minimum permpval
gene_stats <- do.call(rbind, lapply(split(merged_valid, merged_valid$geneId), function(df) {
  idx <- which.min(df$permpval)
  data.frame(
    geneId = df$geneId[idx],
    geneSymbol = df$geneSymbol[idx],
    min_permpval = df$permpval[idx],
    Rho = df$Rho[idx],
    P = df$P[idx],
    n_cnees = nrow(df),
    stringsAsFactors = FALSE
  )
}))

message(sprintf("Gene-level: %d genes", nrow(gene_stats)))

# Gene-level statistic: sign(Rho) * -log10(permpval)
gene_stats$stat <- sign(gene_stats$Rho) * (-log10(gene_stats$min_permpval))
gene_stats <- gene_stats[order(gene_stats$stat, decreasing = TRUE), ]

gene_list <- gene_stats$stat
names(gene_list) <- gene_stats$geneId
gene_list <- sort(gene_list, decreasing = TRUE)

message(sprintf("Gene list: %d genes, range [%.4f, %.4f]",
                length(gene_list), min(gene_list), max(gene_list)))

# ============================================================
# D3: Run GSEA
# ============================================================

message("=== D3: Running GSEA ===")

# Chicken GO
message("Chicken GO GSEA...")
gsea_go_chicken <- tryCatch({
  gseGO(geneList = gene_list, OrgDb = org.Gg.eg.db, keyType = "ENTREZID",
        ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,
        minGSSize = 10, maxGSSize = 500, verbose = FALSE)
}, error = function(e) { message(sprintf("Error: %s", e$message)); NULL })

if (!is.null(gsea_go_chicken) && nrow(as.data.frame(gsea_go_chicken)) > 0) {
  write.csv(as.data.frame(gsea_go_chicken), file.path(OUTPUT_DIR, "gsea_chicken_go_cont.csv"), row.names = FALSE)
  message(sprintf("Chicken GO: %d terms (FDR<0.05)", nrow(as.data.frame(gsea_go_chicken))))
} else {
  write.csv(data.frame(), file.path(OUTPUT_DIR, "gsea_chicken_go_cont.csv"), row.names = FALSE)
  message("Chicken GO: no significant terms")
}

# Human ortholog mapping
message("Mapping to human orthologs...")
human_map <- homologene::homologene(gene_stats$geneId, inTax = 9031, outTax = 9606)
message(sprintf("Mapped %d chicken -> %d human genes",
                length(unique(human_map$query)), length(unique(human_map$gene_id))))

gene_stats_human <- merge(gene_stats, human_map, by.x = "geneId", by.y = "query", all.x = FALSE)
gene_stats_human <- gene_stats_human[!duplicated(gene_stats_human$gene_id), ]

human_list <- gene_stats_human$stat
names(human_list) <- gene_stats_human$gene_id
human_list <- sort(human_list, decreasing = TRUE)
message(sprintf("Human gene list: %d genes", length(human_list)))

# Human GO
message("Human GO GSEA...")
gsea_go_human <- tryCatch({
  gseGO(geneList = human_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,
        minGSSize = 10, maxGSSize = 500, verbose = FALSE)
}, error = function(e) { message(sprintf("Error: %s", e$message)); NULL })

if (!is.null(gsea_go_human) && nrow(as.data.frame(gsea_go_human)) > 0) {
  write.csv(as.data.frame(gsea_go_human), file.path(OUTPUT_DIR, "gsea_human_go_cont.csv"), row.names = FALSE)
  message(sprintf("Human GO: %d terms (FDR<0.05)", nrow(as.data.frame(gsea_go_human))))
} else {
  write.csv(data.frame(), file.path(OUTPUT_DIR, "gsea_human_go_cont.csv"), row.names = FALSE)
  message("Human GO: no significant terms")
}

# KEGG
message("KEGG GSEA...")
gsea_kegg <- tryCatch({
  gseKEGG(geneList = human_list, organism = "hsa", pAdjustMethod = "BH",
          pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, verbose = FALSE)
}, error = function(e) { message(sprintf("Error: %s", e$message)); NULL })

if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  write.csv(as.data.frame(gsea_kegg), file.path(OUTPUT_DIR, "gsea_kegg_cont.csv"), row.names = FALSE)
  message(sprintf("KEGG: %d pathways (FDR<0.05)", nrow(as.data.frame(gsea_kegg))))
} else {
  write.csv(data.frame(), file.path(OUTPUT_DIR, "gsea_kegg_cont.csv"), row.names = FALSE)
  message("KEGG: no significant pathways")
}

# ============================================================
# D4: Compare GSEA vs ORA
# ============================================================

message("=== D4: Comparing GSEA vs ORA ===")

ora_file <- "/mnt/results/enrichment_intronic_chicken_go_cont.csv"
if (file.exists(ora_file)) {
  ora_df <- read.csv(ora_file)
  gsea_terms <- if (!is.null(gsea_go_chicken)) as.data.frame(gsea_go_chicken)$Description else character(0)
  ora_terms <- ora_df$Description

  comparison <- data.frame(
    category = c("GSEA only", "ORA only", "Shared"),
    count = c(length(setdiff(gsea_terms, ora_terms)),
              length(setdiff(ora_terms, gsea_terms)),
              length(intersect(gsea_terms, ora_terms)))
  )
  write.csv(comparison, file.path(OUTPUT_DIR, "gsea_vs_ora_comparison.csv"), row.names = FALSE)
  message(sprintf("GSEA vs ORA: %d shared, %d GSEA-only, %d ORA-only",
                  comparison$count[3], comparison$count[1], comparison$count[2]))
} else {
  message("ORA results not found, skipping comparison")
}

write.csv(gene_stats, file.path(OUTPUT_DIR, "gene_level_stats.csv"), row.names = FALSE)
message(sprintf("Saved: %s/gene_level_stats.csv", OUTPUT_DIR))

message("\n=== Phase D complete ===")
