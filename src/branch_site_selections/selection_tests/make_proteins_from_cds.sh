#!/bin/bash
#SBATCH --job-name=proteins
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-100
#SBATCH --partition=hpc_l40s
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy
source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

mkdir -p logs

PREFIX="trimmed_prank_"
PATH_TO_SRC="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src"
PATH_TO_TMP="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp"
PATH_TO_INITIAL="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/initial_data"
RESULT_DIR="${PATH_TO_TMP}/${PREFIX}small_prank_fa_proteins"

mkdir -p "$RESULT_DIR"

TASK_ID=$SLURM_ARRAY_TASK_ID

gene=$(sed -n "$((TASK_ID + 1))p" ${PREFIX}busted_e_small_table_fdr05.tsv | cut -f1)

echo "Processing gene: $gene (Task ID: $TASK_ID)"

aln="${PATH_TO_TMP}/trimmed_small_prank_fa/${gene}.fa"

if [ -s "$out" ]; then
    echo "Skipping $gene — result already exists ($out)"
    continue
fi

echo ">>> Running translation on $gene"

conda activate seqkit
seqkit translate $aln \
  -o $RESULT_DIR/$gene.fasta \
  -T 1

echo "Job $SLURM_ARRAY_TASK_ID completed successfully"

