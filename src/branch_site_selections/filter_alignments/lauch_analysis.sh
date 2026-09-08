#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=alignment_prepare
#SBATCH --output=logs/alignment_prepare_%j.out
#SBATCH --error=logs/alignment_prepare_%j.err
#SBATCH --nodes=1
#SBATCH --partition=hpc_l40_a
#SBATCH --ntasks=1
#SBATCH --mem=32G  
#SBATCH --cpus-per-task=16 
#SBATCH --mail-user=savouriess2112@gmail.com

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

conda activate base
python /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/branch_site_selections/filter_alignments/stats_table.py

# conda activate mafft

# mkdir -p "${PATH_TO_DATA}/tmp/trimmed/"
# for f in "${PATH_TO_DATA}/initial_data/cds_TOGA/*.fa"; do
#     trimal -in "$f" -out "$PATH_TO_DATA/tmp/trimmed/${f##*/}" \
#         -gt 0.8 \
#         -cons 60 \
#         -resoverlap 0.70 \
#         -seqoverlap 80 
# done



# mkdir -p "${PATH_TO_DATA}/trimmed/"

# conda activate mafft
# for f in cds_TOGA_clean_fasta/*.fa; do
#     trimal -in "$f" -out "$PATH_TO_DATA/trimmed/${f##*/}" -gt 0.8
# done