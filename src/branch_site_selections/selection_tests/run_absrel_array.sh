#!/bin/bash
#SBATCH --job-name=gene_analysis
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-100
#SBATCH --partition=hpc_l40s
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy

mkdir -p logs
PREFIX="trimmed_prank_"

PATH_TO_SRC="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src"
PATH_TO_TMP="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp"
PATH_TO_INITIAL="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/initial_data"
RESULT_DIR="${PATH_TO_TMP}/${PREFIX}absrel_small_results_inno"

#tree_file="${PATH_TO_INITIAL}/trees_labeled/pruned_tree_TOGA_small_selection_marked.newick"

# mkdir -p "$TREE_DIR"
mkdir -p "$RESULT_DIR"

# Получаем индекс задачи в массиве
TASK_ID=$SLURM_ARRAY_TASK_ID

# Получаем строку из файла (пропускаем заголовок, добавляем 1 к индексу)
gene=$(sed -n "$((TASK_ID + 1))p" ${PREFIX}busted_e_small_table_fdr05.tsv | cut -f1)

echo "Processing gene: $gene (Task ID: $TASK_ID)"

aln="${PATH_TO_TMP}/trimmed_small_phylip/${gene}.phy"
out="${RESULT_DIR}/${gene}.json"
tree_file="${PATH_TO_TMP}/trees/pruned_tree_${gene}.newick"

if [ -s "$out" ]; then
    echo "Skipping $gene — result already exists ($out)"
    continue
fi

echo ">>> Running aBSREL on $gene"
hyphy absrel \
    --alignment "$aln" \
    --tree "$tree_file" \
    --branches TEST \
    --output "$out" \
    --code Universal \
    --multiple-hits None

echo "Job $SLURM_ARRAY_TASK_ID completed successfully"