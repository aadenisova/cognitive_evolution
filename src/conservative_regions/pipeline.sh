# move annotations from GalGal4 to bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic
liftoff_to_gal.sh
# split maf file by chromosome 
split_maf_by_chromosome.sh
# extract conservative regions from maf_file
get_cons_regs.sh
# check number of species in the conservative regions (.maf files)


# you need to filter all alignments with huge number of missed species (e in .maf)

# reestimate branches in the conservative regions
reestimate_branches.sh
# RERConverge
run_RERConverge.sh