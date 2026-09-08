# conda activate r-environment
# cat TOGA_ALL_selected_polished_for_RER/* > TOGA_ALL_new_tree_reest.txt

library(RERconverge)
library(tidyverse)
library(ape)

PATH_TO_DATA="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/conserved_regions/"
PATH_TO_WGS<-"/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/vgpBirds_137_way"

treeFile <- paste(
    PATH_TO_DATA, "Sackton_et_al/rer_trees_birds_no_na_by_final_cnees_long.txt", 
    sep="/"
)
cons_regs <- read_tsv(
  treeFile,
  col_names = c("gene", "tree"),
)

inno_rates <- read_csv(paste(PATH_TO_DATA, "df_merged_for_rer_no_na.csv", sep = "/"))

cols_of_interest <- c(
    "TOTALINNOVATIONS2025_ResEff" 
#    "FOODINNO2025_ResEff", "TECHINNO2025_ResEff", "Relative_brain_size", "Wing_length"
)
col <- cols_of_interest[1]

for (col in cols_of_interest) {
    df <- read_tsv(
        paste("Sackton_et_al/RER_results/", col, "RERConverge_result_04_12.all.tsv", sep = "_"))
    
    df_filtered <- df %>% filter(P <= 0.0005)
    # df_filtered <- df %>% head(100)
    
    df_regs <- cons_regs %>% filter(gene %in% df_filtered$gene)
    write_tsv(
        df_regs,
        "Sackton_et_al/tmp/short_df.txt",
        col_names = FALSE
    )

    ourTrees = readTrees("Sackton_et_al/tmp/short_df.txt")

    sp_vect <- setNames(pull(inno_rates[col]), inno_rates$sci_name_2025)
    charpaths=char2Paths(sp_vect, ourTrees)

    res=RERconverge::getAllCor(
        mamRERw, 
        charpaths, 
        method = "p", 
        min.pos = 0, 
        winsorizeRER = 3,
    ) 

    stats=RERconverge::getStat(res)

    annots=RERconverge::read.gmt("Sackton_et_al/gmtfile.gmt")
    annotlist=list(annots)
    names(annotlist)="MSigDBpathways"
    enrichment=RERconverge::fastwilcoxGMTall(stats, annotlist, outputGeneVals=T, num.g=10)

    mt=ourTrees$masterTree
    
    # Palaeognathae	Tinamiformes
    mt=root.phylo(
        mt, 
        outgroup="GCA_047922985", 
        resolve.root=T
    )
    perms=RERconverge::getPermsContinuous(
        100, 
        sp_vect, 
        mamRERw, 
        annotlist,
        ourTrees, 
        mt,
        calculateenrich=F
    )

    corpermpvals=RERconverge::permpvalcor(res, perms)

    res$permpval=corpermpvals[match(rownames(res), rownames(corpermpvals)),]
    res$permpvaladj=p.adjust(res$permpval$permpval, method="BH")
    # count=1

    # while(count<=length(enrichment)){
    #     enrichment[[count]]$permpval=enrichpermpvals[[count]][match(rownames(enrichment[[count]]), names(enrichpermpvals[[count]]))]
    #     enrichment[[count]]$permpvaladj=p.adjust(enrichment[[count]]$permpval, method="BH")
    #     count=count+1
    # }
    
    write.csv(
        res, 
        file = paste(
            paste(PATH_TO_DATA, "Sackton_et_al/permutation_analysis/" , sep="/"), 
            col, 
            "RERConverge_result_04_20.all.sign.tsv", 
            sep="_"
        )
    )
}



# df <- read_tsv(
#     paste("Sackton_et_al/RER_results/", col, "RERConverge_result_04_12.all.tsv", sep = "_"))

# df_filtered <- df %>% filter(P <= 0.0005)
# # df_filtered <- df %>% head(100)

# df_regs <- cons_regs %>% filter(gene %in% df_filtered$gene)
# write_tsv(
#     df_regs,
#     "Sackton_et_al/tmp/short_df.txt",
#     col_names = FALSE
# )

# ourTrees = readTrees("Sackton_et_al/tmp/short_df.txt")

# sp_vect <- setNames(pull(inno_rates[col]), inno_rates$sci_name_2025)
# charpaths=char2Paths(sp_vect, ourTrees)

# res=RERconverge::getAllCor(
#     mamRERw, 
#     charpaths, 
#     method = "p", 
#     min.pos = 0, 
#     winsorizeRER = 3,
# ) 

# stats=RERconverge::getStat(res)

# annots=RERconverge::read.gmt("Sackton_et_al/gmtfile.gmt")
# annotlist=list(annots)
# names(annotlist)="MSigDBpathways"
# enrichment=RERconverge::fastwilcoxGMTall(stats, annotlist, outputGeneVals=T, num.g=10)

# mt=ourTrees$masterTree

# # Palaeognathae	Tinamiformes
# mt=root.phylo(
#     mt, 
#     outgroup="GCA_047922985", 
#     resolve.root=T
# )
# perms=RERconverge::getPermsContinuous(
#     100, 
#     sp_vect, 
#     mamRERw, 
#     annotlist,
#     ourTrees, 
#     mt,
#     calculateenrich=F
# )

# corpermpvals=RERconverge::permpvalcor(res, perms)
# # enrichpermpvals=RERconverge::permpvalenrich(enrichment, perms)

# res$permpval=corpermpvals[match(rownames(res), rownames(corpermpvals)),]
# res$permpvaladj=p.adjust(res$permpval$permpval, method="BH")
# # count=1

# # while(count<=length(enrichment)){
# #     enrichment[[count]]$permpval=enrichpermpvals[[count]][match(rownames(enrichment[[count]]), names(enrichpermpvals[[count]]))]
# #     enrichment[[count]]$permpvaladj=p.adjust(enrichment[[count]]$permpval, method="BH")
# #     count=count+1
# # }

# write.csv(
#     res, 
#     file = paste(
#         paste(PATH_TO_DATA, "Sackton_et_al/permutation_analysis/" , sep="/"), 
#         col, 
#         "RERConverge_result_04_20.all.sign.tsv", 
#         sep="_"
#     )
# )