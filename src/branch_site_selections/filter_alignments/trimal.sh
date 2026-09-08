#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=alignment_prepare
#SBATCH --output=logs/alignment_prepare_%j.out
#SBATCH --error=logs/alignment_prepare_%j.err
#SBATCH --nodes=1
#SBATCH --partition=hpc_a10_a
#SBATCH --ntasks=1
#SBATCH --mem=32G  
#SBATCH --cpus-per-task=16 
#SBATCH --mail-user=savouriess2112@gmail.com
#SBATCH --array=1-100         # 1000 задач по ~16 файлов каждая

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

conda activate mafft

mkdir -p "${PATH_TO_DATA}/tmp/trimmed/"

# Константы
FILES_PER_JOB=160  # Файлов на одну задачу
ALL_FILES=($PATH_TO_DATA/initial_data/cds_TOGA/*.fa)
TOTAL_FILES=${#ALL_FILES[@]}

# Вычисляем индексы файлов для этой задачи
START_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) * FILES_PER_JOB ))
END_IDX=$(( START_IDX + FILES_PER_JOB - 1 ))

# Обрабатываем пакет файлов
for (( i=START_IDX; i<=END_IDX && i<TOTAL_FILES; i++ )); do
    INPUT_FILE="${ALL_FILES[$i]}"
    FILENAME=$(basename "$INPUT_FILE")
    OUTPUT_FILE="${PATH_TO_DATA}/tmp/trimmed/${FILENAME}"
    
    # # Пропускаем если уже обработан
    # if [[ -f "$OUTPUT_FILE" ]]; then
    #     echo "Пропускаем (уже существует): $FILENAME"
    #     continue
    # fi
    
    echo "Обрабатываю: $FILENAME"
    
    # Обработка с использованием нескольких ядер внутри задачи
    trimal -in "$INPUT_FILE" -out "$OUTPUT_FILE" \
        -gt 0.8 \
        -cons 60 \
        -resoverlap 0.70 \
        -seqoverlap 80 &
    
    # Ограничиваем количество параллельных процессов внутри задачи
    # if (( (i - START_IDX + 1) % 50 == 0 )); then
    #     wait
    # fi
done

wait
echo "Задача $SLURM_ARRAY_TASK_ID завершена"