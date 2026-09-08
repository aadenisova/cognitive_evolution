# conda activate r-environment
# cat TOGA_ALL_selected_polished_for_RER/* > TOGA_ALL_new_tree_reest.txt

library(RERconverge)
library(tidyverse)
library(ape)

PATH_TO_DATA="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/conserved_regions/"
PATH_TO_WGS<-"/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/vgpBirds_137_way"

treeFile <- paste(PATH_TO_DATA, "Sackton_et_al/rer_trees_birds_no_na_by_final_cnees_long_head.txt", sep="/")
inno_rates <- read_csv(paste(PATH_TO_DATA, "df_merged_for_rer_no_na.csv", sep = "/"))

ourTrees=readTrees(
    treeFile
)

masterTree <- ourTrees$masterTree

mamRERw=RERconverge::getAllResiduals(
    ourTrees,
    useSpecies=inno_rates$sci_name_2025,
    transform = "sqrt", 
    weighted = T,
    scale = T,
)

cols_of_interest <- c("TOTALINNOVATIONS2025_ResEff", "FOODINNO2025_ResEff", "TECHINNO2025_ResEff", "Relative_brain_size", "Wing_length")

col=cols_of_interest[1]

for (col in cols_of_interest) {
    print(col)

    sp_vect <- setNames(pull(inno_rates[col]), inno_rates$sci_name_2025)
    charpaths=char2Paths(sp_vect, ourTrees)

    res=RERconverge::getAllCor(
        mamRERw, 
        charpaths, 
        method = "p", 
        min.pos = 0, 
        winsorizeRER = 3,
    ) 

    # ordered_results = res[order(res$P),]
    # head(ordered_results)
    # ordered_results$gene <- row.names(ordered_results)

    # ordered_results %>% head()

    # write_tsv(
    #     ordered_results,
    #     paste(paste(PATH_TO_DATA, "Sackton_et_al/RER_results/" , sep="/"), col, "RERConverge_result_04_13.92.tsv", sep="_"),
    # )

    restats=RERconverge::getStat(res)
    mt=ourTrees$masterTree
    
    # Palaeognathae	Tinamiformes
    mt=root.phylo(
        mt, 
        outgroup="GCA_047922985", 
        resolve.root=T
    )
    
    perms=RERconverge::getPermsContinuous(
        1000, 
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
    
    write.csv(
        res, 
        file = paste(
            paste(PATH_TO_DATA, "Sackton_et_al/permutation_analysis/" , sep="/"), 
            col, 
            "RERConverge_result_04_20.all.head.head.tsv", 
            sep="_"
        )
    )
}
