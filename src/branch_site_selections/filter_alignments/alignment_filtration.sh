#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=alignment_prepare
#SBATCH --output=logs/alignment_prepare_%j.out
#SBATCH --error=logs/alignment_prepare_%j.err
#SBATCH --nodes=1
#SBATCH --partition=hpc_l40s
#SBATCH --ntasks=1
#SBATCH --mem=32G  
#SBATCH --cpus-per-task=16 
#SBATCH --mail-user=savouriess2112@gmail.com
#SBATCH --array=1-5         # 1000 задач по ~16 файлов каждая

export PATH=$PATH:/ru-auth/local/home/adenisova/programs/prank/bin
GUI_PATH=/ru-auth/local/home/adenisova/programs/guidance/www/Guidance

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp
cd $PATH_TO_DATA

ALL_FILES=($PATH_TO_DATA/cds_small_TOGA_clean_fasta/*.fa)

i=3

INPUT_FILE="${ALL_FILES[$i]}"

INPUT_FILE=cds_small_TOGA_clean_fasta/SLC25A12_rna-XM_015289764.2.fa
echo $INPUT_FILE

conda activate seqkit
seqkit seq -g $INPUT_FILE > cds_nogaps.fa

tree=../initial_data/trees_labeled/pruned_tree_TOGA_small_selection.newick

prank -d=cds_nogaps.fa -t=$tree -o=prank_output -f=fasta -codon -F -showall -once

# seqkit fx2tab cds_nogaps.fa | awk '{print $1, length($2)%3}'
#transeq -sequence your_alignment.fa -outseq tmp_aa.fa


# conda activate bioperl_env
# perl $GUI_PATH/guidance.pl --seqFile cds_nogaps.fa \
#     --msaProgram PRANK \
#     --seqType codon \
#     --proc_num 8 \
#     --outDir guidence_results3

# trimal -in $INPUT_FILE -out aln.trimal.fa -gappyout

# conda activate bioperl_env
# perl $GUI_PATH/guidance.pl \
#   --msaFile $INPUT_FILE \
#   --seqType codon \
#   --outDir guidance_results \
#   --program GUIDANCE

# perl $GUI_PATH/guidance.pl \
#   --seqFile $INPUT_FILE \
#   --msaProgram PRANK \
#   --msaFile $INPUT_FILE \
#   --seqType codon \
#   --outDir guidance_results \