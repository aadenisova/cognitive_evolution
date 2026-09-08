#!/bin/bash
#SBATCH --job-name=maf_stats
#SBATCH --output=logs/maf_stats.out
#SBATCH --error=logs/maf_stats.err
#SBATCH --partition=hpc_l40s    
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Get statistics of the original MAF file

PATH_TO_WGS=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/vgpBirds_137_way
cd $PATH_TO_WGS

REF=GCF_016700215.2 # vgpBirds.137way.GCF_016700215.2

# mafStats --maf vgpBirds.137way.$REF.maf > maf.stats.$REF.txt

# Get statistics of the MAF file after filtering duplicates (not used)
# mafDuplicateFilter --maf vgpBirds.137way.$REF.maf > maf.nodup.maf
# mafStats --maf maf.nodup.maf > maf.nodup.stats.txt

# Break the MAF file into smaller files by chromosome

mkdir -p $PATH_TO_WGS/maf_by_chr_$REF
mafSplit -byTarget -useFullSequenceName custom_bed.bed maf_by_chr_$REF/ vgpBirds.137way.$REF.maf

