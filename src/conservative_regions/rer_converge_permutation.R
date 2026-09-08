library(RERconverge)
library(tidyverse)

setwd("/Users/aleksandradenisova/Laboratory/innovation/Inno_2025/conserved_regions")

cols_of_interest <- c(
  # "TOTALINNOVATIONS2025_ResEff", 
  # "FOODINNO2025_ResEff", 
  "TECHINNO2025_ResEff", 
  # "Wing_length",
  "Relative_brain_size"
)

cons_regs <- read_tsv(
  "selected_cons_regs/rer_trees_birds_no_na_by_final_cnees_long.txt",
  col_names = c("gene", "tree"),
)
inno_rates <- read_csv("selected_cons_regs/df_merged_for_rer_no_na.csv")

for (col in cols_of_interest) {
  df <- read_tsv(paste("", col, "RERConverge_result_04_12.all.tsv", sep = "_"))
  
  df_filtered <- df %>% filter(P <= 0.05)
  # df_filtered <- df %>% head(100)
  
  df_regs <- cons_regs %>% filter(gene %in% df_filtered$gene)
  write_tsv(
    df_regs,
    "short_df.txt",
    col_names = FALSE
  )

  ourTrees = readTrees("short_df.txt")

  mamRERw=RERconverge::getAllResiduals(
    ourTrees,
    useSpecies=inno_rates$sci_name_2025,
    transform = "sqrt", 
  )
  
  sp_vect <- setNames(pull(inno_rates[col]), inno_rates$sci_name_2025)
  charpaths=char2Paths(sp_vect, ourTrees)
  
  res=correlateWithContinuousPhenotype(
    mamRERw, 
    charpaths, 
    min.sp = 10,
    winsorizeRER = 3, 
    winsorizetrait = 3
  )
  
  restats=RERconverge::getStat(df)

  mt=ourTrees$masterTree
  mt
  
  # Palaeognathae	Tinamiformes
  plot(mt)
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
      "permutation_analysis", 
      paste(col, "txt", sep = "."), 
      sep = "/"
    )
  )
  
}

ourTrees = readTrees("short_df.txt")



write_tsv(
  res,
  paste(
    "permutation_analysis", 
    paste(col, "txt", sep = "/"), 
    sep = "/"
  )
)
res
paste(
  "permutation_analysis", 
  paste(col, "txt", sep = "/"), 
  sep = "/"
)

res %>% nrow()

col <- "TECHINNO2025_ResEff"
df <- read_tsv(paste("", col, "RERConverge_result_04_12.all.tsv", sep = "_"))
paste("", col, "RERConverge_result_04_12.all.tsv", sep = "_")
df
cons_regs <- read_tsv(
  "selected_cons_regs/rer_trees_birds_no_na_by_final_cnees_long.txt",
  col_names = c("gene", "tree"),
)
cons_regs

df

df_filtered <- df %>% filter(P < 0.05)
df_filtered

cons_regs
df_regs <- cons_regs %>% filter(gene %in% df_filtered$gene)
df_regs 



write_tsv(
  df_regs,
  "short_df.txt",
  col_names = FALSE
)

ourTrees = readTrees("short_df.txt")
masterTree <- ourTrees$masterTree

inno_rates <- read_csv("selected_cons_regs/df_merged_for_rer_no_na.csv")

mamRERw=RERconverge::getAllResiduals(
  ourTrees,
  useSpecies=inno_rates$sci_name_2025,
  transform = "sqrt", 
)

sp_vect <- setNames(pull(inno_rates[col]), inno_rates$sci_name_2025)
charpaths=char2Paths(sp_vect, ourTrees)

res=correlateWithContinuousPhenotype(
  mamRERw, 
  charpaths, 
  min.sp = 10,
  winsorizeRER = 3, 
  winsorizetrait = 3
)
res



ordered_results = res[order(res$P),]
head(ordered_results)
ordered_results

res$gene <- res %>% rownames()
res %>% filter(gene == "mCE1684273")
df_filtered

res %>% typeof()

restats=RERconverge::getStat(df)
stats
annots=RERconverge::read.gmt("gmtfile.gmt")
annotlist=list(annots)

names(annotlist)="MSigDBpathways"
enrichment=RERconverge::fastwilcoxGMTall(restats, annotlist, outputGeneVals=T, num.g=10)

enrichment

library(ape)

mt=ourTrees$masterTree
mt

# Palaeognathae	Tinamiformes
plot(mt)
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
res[order(res$p.adj),] %>% head()
res[order(res$p.adj),] %>% filter(gene == "mCE1342012")

res[order(res$permpvaladj, decreasing=TRUE),] %>% head(50)

res$gene <- res %>% rownames()

permsnoenrich=RERconverge::getPermsContinuous(
  10, 
  sp_vect, 
  mamRERw, 
  annotlist,
  ourTrees, 
  mt, 
  calculateenrich = F
)

permswithenrich=RERconverge::getEnrichPerms(permsnoenrich, enrichment,
                                            annotlist)
corpermpvals2=RERconverge::permpvalcor(res, permswithenrich)
enrichpermpvals2=RERconverge::permpvalenrich(enrichment, permswithenrich)
plot(corpermpvals,corpermpvals2)

