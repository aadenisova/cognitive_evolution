#!/bin/bash
#SBATCH --job-name=meme
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
RESULT_DIR="${PATH_TO_TMP}/${PREFIX}meme_prank_small_results_inno"

# tree_file="${PATH_TO_INITIAL}/trees_labeled/pruned_tree_TOGA_small_selection_marked.newick"

mkdir -p "$RESULT_DIR"

TASK_ID=$SLURM_ARRAY_TASK_ID

gene=$(sed -n "$((TASK_ID + 1))p" ${PREFIX}busted_e_small_table_fdr05.tsv | cut -f1)
gene=STAB1_ENSGALT00000079354.2

echo "Processing gene: $gene (Task ID: $TASK_ID)"

# gene=SLC16A6_rna-XM_015279952.2
aln="${PATH_TO_TMP}/trimmed_small_prank_phylip/${gene}.phy"
out="${RESULT_DIR}/${gene}.json"

# if [ -s "$out" ]; then
#     echo "Skipping $gene — result already exists ($out)"
#     continue
# fi

echo ">>> Running meme on $gene"

tree_file="${PATH_TO_TMP}/trees/pruned_tree_${gene}.newick"
# hyphy meme \
#     --alignment "$aln" \
#     --tree $tree_file \
#     --code Universal \
#     --branches TEST \
#     --output "$out" > "${RESULT_DIR}/${gene}.out"



python ${PATH_TO_SRC}/branch_site_selections/selection_tests/mark_tree_branches.py \
    $PREFIX \
    $tree_file \
    $gene

hyphy meme \
    --alignment "$aln" \
    --tree "${PATH_TO_TMP}/meme_trees/$gene.newick" \
    --code Universal \
    --branches TEST \
    --output "$out" > "${RESULT_DIR}/${gene}.out"

echo "Job $SLURM_ARRAY_TASK_ID completed successfully"



# hyphy meme \
#     --alignment tutorial_data/h3_trunk.fna \
#     --code Universal \
#     --output h3_trunk.json