#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=alignment_prepare
#SBATCH --output=logs/alignment_prepare_%j.out
#SBATCH --error=logs/alignment_prepare_%j.err
#SBATCH --nodes=1
#SBATCH --partition=hpc_l40_a 
#SBATCH --ntasks=1
#SBATCH --mem=4G  
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=savouriess2112@gmail.com
#SBATCH --array=1-100        

# export PATH=$PATH:/ru-auth/local/home/adenisova/programs/prank/bin

Prefix="TOGA_ALL"

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

tree=$PATH_TO_DATA/initial_data/trees_labeled/pruned_tree_TOGA_small_selection.newick

OUTPUT_DIR_FA=${PATH_TO_DATA}/tmp/${Prefix}_prank_fa
OUTPUT_DIR_PHYLIP=${PATH_TO_DATA}/tmp/${Prefix}_prank_phylip
OUTPUT_DIR_NO_GAPS=${PATH_TO_DATA}/tmp/${Prefix}_nogaps
OUTPUT_DIR_PRANK=${PATH_TO_DATA}/tmp/${Prefix}_prank

mkdir -p "${OUTPUT_DIR_FA}" "${OUTPUT_DIR_PHYLIP}" "${OUTPUT_DIR_NO_GAPS}" "${OUTPUT_DIR_PRANK}"

# Константы
FILES_PER_JOB=160  # Файлов на одну задачу
ALL_FILES=($PATH_TO_DATA/tmp/${Prefix}_clean_fasta/*.fa)
TOTAL_FILES=${#ALL_FILES[@]}

# Вычисляем индексы файлов для этой задачи
START_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) * FILES_PER_JOB ))
END_IDX=$(( START_IDX + FILES_PER_JOB - 1 ))


# Обрабатываем пакет файлов
for (( i=START_IDX; i<=END_IDX && i<TOTAL_FILES; i++ )); do
    INPUT_FILE="${ALL_FILES[$i]}"
    RELATIVE_PATH="${INPUT_FILE#$PATH_TO_DATA/tmp/${Prefix}_clean_fasta/}"
    GENENAME="${RELATIVE_PATH%.fa}"

    echo "Обрабатываю: $GENENAME"

    conda activate seqkit
    seqkit seq -g $INPUT_FILE > ${OUTPUT_DIR_NO_GAPS}/$GENENAME.fa

    cd ${OUTPUT_DIR_PRANK}
    prank -d=${OUTPUT_DIR_NO_GAPS}/$GENENAME.fa \
        -t=$tree \
        -o=$GENENAME \
        -f=fasta -codon \
        -F -showall -once -prunetree

    conda activate mafft
    OUTPUT_FILE="${OUTPUT_DIR_FA}/${GENENAME}.fa"
    trimal -in ${OUTPUT_DIR_PRANK}/$GENENAME.best.fas -out "$OUTPUT_FILE" -gappyout

    OUTPUT_FILE="${OUTPUT_DIR_PHYLIP}/${GENENAME}.phy"
    # GENENAME=STAB1_ENSGALT00000079354.2
    trimal -in ${OUTPUT_DIR_PRANK}/$GENENAME.best.fas -out "$OUTPUT_FILE" -gappyout -phylip

done

echo "Задача $SLURM_ARRAY_TASK_ID завершена"

# GENENAME="rna-XM_004936272.3"
# file="$PATH_TO_DATA/tmp/cds_small_TOGA_clean_fasta/$GENENAME".fa
# seqkit fx2tab $file