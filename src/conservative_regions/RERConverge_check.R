# conda activate r-environment


library(RERconverge)
library(tidyverse)
library(ape)

PATH_TO_DATA="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/conserved_regions/"
PATH_TO_WGS<-"/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/vgpBirds_137_way"

treeFile <- paste(PATH_TO_DATA, "Sackton_et_al/rer_trees_birds_no_na_by_final_cnees_long.txt", sep="/")
inno_rates <- read_csv(paste(PATH_TO_DATA, "df_merged_for_rer_no_na.csv", sep = "/"))

cols_of_interest <- c("TOTALINNOVATIONS2025_ResEff", "FOODINNO2025_ResEff", "TECHINNO2025_ResEff", "Relative_brain_size", "Wing_length")
col <- "FOODINNO2025_ResEff"

result_df <- read_tsv(paste(paste(PATH_TO_DATA, "Sackton_et_al/RER_results/" , sep="/"), col, "RERConverge_result_04_13.92.tsv", sep="_"))
result_df %>% filter(P <= 0.05)


ourTrees=readTrees(
    treeFile
)
