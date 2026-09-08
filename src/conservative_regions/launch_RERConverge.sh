#!/bin/bash

#SBATCH --time=14-00:00:00
#SBATCH --job-name=RERConverge
#SBATCH --output=logs/RERConverge_%j.out
#SBATCH --error=logs/RERConverge_%j.err
#SBATCH --nodes=1
#SBATCH --partition=vgl_b 
#SBATCH --ntasks=1
#SBATCH --mem=120G  
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=savouriess2112@gmail.com  
#SBATCH --account=jarv_condo_bank

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

# PATH_TO_FILES=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp
# Prefix="TOGA_ALL_selected"

# cd $PATH_TO_FILES

# cat ${Prefix}_polished_for_RER/* > ${Prefix}_new_tree_reest.txt

# find rer_trees_birds_no_na_by_final_cnees_long/ \
#     -type f -size +1000c -exec cat {} \; \
#     > rer_trees_birds_no_na_by_final_cnees_long.txt

# conda activate base
# python /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/whole_tree_analysis/cds/RERConverge/last_prep.py ${Prefix}

conda activate r-environment

Rscript /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/conservative_regions/rerconverge_analysis.R
# Rscript /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/conservative_regions/rerconverge_analysis_sign.R