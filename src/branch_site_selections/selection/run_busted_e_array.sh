#!/bin/bash
#SBATCH --job-name=gene_analysis
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-100
#SBATCH --partition=hpc_l40_a
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy

mkdir -p logs

prefix=TOGA_ALL
phenotype=food_type_innovative

PATH_TO_SRC="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src"
PATH_TO_TMP="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp"
PATH_TO_INITIAL="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/initial_data"
# TREE_DIR="${PATH_TO_TMP}/trees"

RESULT_DIR_E="${PATH_TO_TMP}/${prefix}_${phenotype}_busted_e"

# mkdir -p "$TREE_DIR"
mkdir -p "$RESULT_DIR_E"

TABLE_FILE="${PATH_TO_TMP}/${prefix}_filtered_trees.tsv"
TOTAL_ROWS=$(wc -l $TABLE_FILE | awk '{print $1}')           
MAX_JOBS=100                

SLURM_ARRAY_TASK_ID=1
BATCH_SIZE=$(( ($TOTAL_ROWS + $MAX_JOBS - 1) / $MAX_JOBS ))
START_ROW=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_ROW=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

if [ $END_ROW -gt $TOTAL_ROWS ]; then
    END_ROW=$TOTAL_ROWS
fi

echo "Job ID: $SLURM_ARRAY_TASK_ID"
echo "Processing rows: $START_ROW to $END_ROW"
echo "Batch size: $BATCH_SIZE"

BATCH_FILE="tmp_batch_${SLURM_ARRAY_TASK_ID}.tsv"

# if [ $START_ROW -eq 1 ]; then
#     # Первый батч включает заголовок
#     head -n $END_ROW "$TABLE_FILE" > "$BATCH_FILE"
# else
#     # Пропускаем заголовок для остальных батчей
#     sed -n "${START_ROW},${END_ROW}p" "$TABLE_FILE" > "$BATCH_FILE"
# fi

sed -n "${START_ROW},${END_ROW}p" "$TABLE_FILE" > "$BATCH_FILE"
LINE_NUM=0
while IFS=$'\t' read -r gene all_species || [ -n "$gene" ]; do
    ((LINE_NUM++))
    ABSOLUTE_ROW=$((START_ROW + LINE_NUM))
    
    #echo "Processing row $ABSOLUTE_ROW: gene=$gene, all_species=$all_species"
    echo "Processing row $ABSOLUTE_ROW: gene=$gene"
    aln="${PATH_TO_TMP}/${prefix}_polished_phy/${gene}.phy"
    out="${RESULT_DIR_E}/${gene}.json"
    
    tree_dir="${PATH_TO_TMP}/trees/${prefix}"
    mkdir -p $tree_dir
    python ${PATH_TO_SRC}/branch_site_selections/selection/filter_tree.py \
        $all_species \
        $gene \
        $phenotype \
        $tree_dir

    tree_file=${tree_dir}/${gene}.newick
    # tree_file="${PATH_TO_INITIAL}/trees_labeled/pruned_tree_TOGA_small_selection_marked.newick"

    hyphy busted \
        ENV=TOLERATE_NUMERICAL_ERRORS=1 \
        --error-sink Yes \
        --alignment "$aln" \
        --tree "$tree_file" \
        --branches TEST \
        --rates 3 \
        --grid-size 250 \
        --starting-points 1 \
        --output "$out" \
        --kill-zero-lengths Yes \
        --code Universal > "${RESULT_DIR_E}/${gene}.out"
   
done < "$BATCH_FILE"

rm "$BATCH_FILE"

echo "Job $SLURM_ARRAY_TASK_ID completed successfully"