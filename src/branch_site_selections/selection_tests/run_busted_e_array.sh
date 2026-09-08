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

phenotype=food_type_and_technical_innovative

PATH_TO_SRC="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src"
PATH_TO_TMP="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp"
PATH_TO_INITIAL="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/initial_data"
# TREE_DIR="${PATH_TO_TMP}/trees"
RESULT_DIR_E="${PATH_TO_TMP}/trimmed_prank_busted_e_small_results_inno"

# mkdir -p "$TREE_DIR"
mkdir -p "$RESULT_DIR_E"

# Конфигурация
TABLE_FILE="${PATH_TO_TMP}/trimmed_alignment_qc_summary_filtered.tsv"  # Ваш файл с таблицей
TOTAL_ROWS=$(wc -l $TABLE_FILE | awk '{print $1}')            # Общее количество строк
MAX_JOBS=100                # Максимальное количество джобов

# Рассчитываем размер батча
BATCH_SIZE=$(( ($TOTAL_ROWS + $MAX_JOBS - 1) / $MAX_JOBS ))

# Определяем строки для текущего джоба
START_ROW=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_ROW=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

# Убедимся, что последний батч не выходит за пределы
if [ $END_ROW -gt $TOTAL_ROWS ]; then
    END_ROW=$TOTAL_ROWS
fi

echo "Job ID: $SLURM_ARRAY_TASK_ID"
echo "Processing rows: $START_ROW to $END_ROW"
echo "Batch size: $BATCH_SIZE"

# Временный файл для текущего батча
BATCH_FILE="tmp_batch_${SLURM_ARRAY_TASK_ID}.tsv"

# Извлекаем нужные строки (предполагаем, что таблица с заголовком)
if [ $START_ROW -eq 1 ]; then
    # Первый батч включает заголовок
    head -n $END_ROW "$TABLE_FILE" > "$BATCH_FILE"
else
    # Пропускаем заголовок для остальных батчей
    sed -n "${START_ROW},${END_ROW}p" "$TABLE_FILE" > "$BATCH_FILE"
fi

# Обрабатываем каждую строку в батче
LINE_NUM=0
while IFS=$'\t' read -r gene all_species other_columns || [ -n "$gene" ]; do
    # Пропускаем заголовок если это не первый батч
    if [ $START_ROW -gt 1 ] && [ $LINE_NUM -eq 0 ]; then
        ((LINE_NUM++))
        continue
    fi
    
    ((LINE_NUM++))
    ABSOLUTE_ROW=$((START_ROW + LINE_NUM - 2))
    
    echo "Processing row $ABSOLUTE_ROW: gene=$gene, all_species=$all_species"
    aln="${PATH_TO_TMP}/trimmed_small_phylip/${gene}.phy"
    out="${RESULT_DIR_E}/${gene}.json"

    # if [ -s "$out" ]; then
    #     echo "Skipping $gene — result already exists ($out)"
    #     continue
    # fi
    
    python ${PATH_TO_SRC}/branch_site_selections/selection_tests/filter_tree.py \
        $all_species \
        $gene \
        $phenotype 
    
    tree_file="${PATH_TO_TMP}/trees/pruned_tree_${gene}.newick"
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

# Удаляем временный файл
rm "$BATCH_FILE"

echo "Job $SLURM_ARRAY_TASK_ID completed successfully"


# PATH_TO_TMP="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp"
# PATH_TO_INITIAL="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/initial_data"
# gene=CAMK2G_rna-XM_015288304.2
# tree_file="${PATH_TO_INITIAL}/trees_labeled/pruned_tree_TOGA_small_selection_marked.newick"

# aln="${PATH_TO_TMP}/cds_small_TOGA_clean_phylip/${gene}.phy"
# out="${PATH_TO_TMP}/busted_small_results_inno/${gene}.json"

# hyphy busted \
#     ENV=TOLERATE_NUMERICAL_ERRORS=1 \
#     --alignment "$aln" \
#     --tree "$tree_file" \
#     --branches TEST \
#     --rates 3 \
#     --grid-size 250 \
#     --starting-points 5 \
#     --output "$out" \
#     --kill-zero-lengths Yes \
#     --code Universal