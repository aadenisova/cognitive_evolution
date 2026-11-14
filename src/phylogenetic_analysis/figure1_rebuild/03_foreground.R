suppressPackageStartupMessages({library(ape);library(phytools)})
PH <- "TECHINNO2025_ResEff"
tree <- read.tree("data/tree_clootl_overlap_25.nwk")
dat <- read.csv("data/inno_db_25_clootl_overlap.csv", stringsAsFactors=FALSE); rownames(dat) <- dat$Species
sp <- read.csv("species_50.csv")$Species
t50 <- keep.tip(tree, sp)
ok <- !is.na(dat[t50$tip.label, PH]); tp <- drop.tip(t50, t50$tip.label[!ok])
x <- setNames(dat[tp$tip.label, PH], tp$tip.label)
anc <- fastAnc(tp, x); nt <- Ntip(tp)
val <- function(n) if (n <= nt) x[tp$tip.label[n]] else anc[as.character(n)]
delta <- vapply(seq_len(nrow(tp$edge)), function(i) unname(val(tp$edge[i,2]) - val(tp$edge[i,1])), 0)
thr <- mean(delta) + sd(delta)
fgc <- tp$edge[delta > thr, 2]
desc <- unique(unlist(lapply(fgc, function(n) if (n<=nt) n else c(n, getDescendants(tp,n)))))
fg <- tp$tip.label[desc[desc<=nt]]
cat("n tips:", nt, " fg branches:", length(fgc), " fg tips:", length(fg), " thr:", thr, "\n")
d <- dat[t50$tip.label,]; d$Range.Size <- suppressWarnings(as.numeric(trimws(d$Range.Size))); d <- d[!is.na(d$Range.Size) & d$Range.Size>0,]; d$fg <- d$Species %in% fg
w <- wilcox.test(log10(Range.Size) ~ fg, data=d)
cat("n bg/fg:", table(d$fg), " wilcox p:", w$p.value, "\n")
write.csv(data.frame(Species=d$Species, Order=d$Order, Family=d$Family, log10_range=log10(d$Range.Size), fg=d$fg), "boxplot_data.csv", row.names=FALSE)
write.csv(data.frame(Species=tp$tip.label, fg=tp$tip.label %in% fg), "foreground_tips.csv", row.names=FALSE)
write.csv(data.frame(p=w$p.value, W=unname(w$statistic), n_bg=sum(!d$fg), n_fg=sum(d$fg)), "wilcox.csv", row.names=FALSE)
