suppressPackageStartupMessages({library(ape);library(phytools)})
PH <- "TECHINNO2025_ResEff"
tree <- read.tree("data/tree_clootl_overlap_25.nwk")
dat <- read.csv("data/inno_db_25_clootl_overlap.csv", stringsAsFactors=FALSE); rownames(dat) <- dat$Species
tree <- keep.tip(tree, intersect(tree$tip.label, dat$Species))
cat("tips", Ntip(tree), " NA pheno:", sum(is.na(dat[tree$tip.label, PH])), " ultrametric", is.ultrametric(tree,tol=1e-4), "\n")
x <- setNames(as.numeric(dat[tree$tip.label, PH]), tree$tip.label)
anc <- fastAnc(tree, x)
nt <- Ntip(tree); val <- c(x, anc[as.character((nt+1):(nt+tree$Nnode))])
names(val) <- NULL
e <- data.frame(parent=tree$edge[,1], child=tree$edge[,2], len=tree$edge.length, pv=val[tree$edge[,1]], cv=val[tree$edge[,2]])
write.csv(e, "tree_edges.csv", row.names=FALSE)
tips <- data.frame(idx=1:nt, Species=tree$tip.label, val=x, Order=dat[tree$tip.label,"Order"], Family=dat[tree$tip.label,"Family"])
write.csv(tips, "tree_tips.csv", row.names=FALSE)
nd <- node.depth.edgelength(tree); write.csv(data.frame(node=1:(nt+tree$Nnode), depth=nd), "tree_depth.csv", row.names=FALSE)
cat("depth", max(nd), " root", nt+1, " val range", range(val), "\n")
