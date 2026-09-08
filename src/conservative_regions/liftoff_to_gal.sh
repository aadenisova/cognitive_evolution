#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=liftoff_batch
#SBATCH --output=liftoff_batch_%j.out
#SBATCH --error=liftoff_batch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G  
#SBATCH --cpus-per-task=16 
#SBATCH --mail-user=savouriess2112@gmail.com

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/conserved_regions/

run_liftoff() {
    local ref_gff=$1
    local ref_fasta=$2
    local target_fasta=$3
    local output_dir=$4

    # Активируем окружение liftoff
    conda activate liftOff

    liftoff -g "$ref_gff" \
        -p 16 \
        -o "$output_dir/lifted_annotations_genomic.gff" \
        -u "$output_dir/unmapped_features.txt" \
        -infer_genes \
        -cds \
        "$target_fasta" "$ref_fasta"

    echo "=== Completed $prefix ==="
}

# For Paleognata 
# for prefix in final_cnees_long cnee_overlaps_ATAC_Seq; do
# genome_file=$PATH_TO_DATA/genome_data/GCF_000002315.6_GRCg6a_genomic.fna
# for prefix in final_cnees_long; do
#     echo "=== Processing $prefix ==="

#     gff_file=$PATH_TO_DATA/Sackton_et_al/selected_elements/$prefix.gff
#     output_dir=$PATH_TO_DATA/Sackton_et_al/liftoff_GRCg6a_$prefix

#     input_bed=$PATH_TO_DATA/Sackton_et_al/selected_elements/$prefix.bed
#     bedtools sort -i $input_bed | \
#     awk 'BEGIN{OFS="\t"} {print $1,"bedtools","gene",$2+1,$3,".",".",".","ID="$4}' \
#     > $gff_file

#     mkdir -p $output_dir
#     run_liftoff \
#         $gff_file \
#         $PATH_TO_DATA/genome_data/GCF_000002315.3_Gallus_gallus-4.0_genomic.fna \
#         $genome_file \
#         $output_dir

#     # conda activate bedops
#     # gff2bed < $output_dir/ > $output_dir/cnee_annotations_genomic.bed
# done



paper=Trigila_et_al_2025
genome_file=$PATH_TO_DATA/genome_data/GCF_000002315.6_GRCg6a_genomic.fna
# for prefix in conserved_elements_coordinates conserved_elements_coordinates_nc; do
for prefix in conserved_elements_merged accelerated_elements_merged; do
    echo "=== Processing $prefix ==="

    gff_file=$PATH_TO_DATA/$paper/selected_elements/$prefix.gff
    output_dir=$PATH_TO_DATA/$paper/liftoff_tg_$prefix
    input_bed=$PATH_TO_DATA/$paper/selected_elements/$prefix.bed

    bedtools sort -i $input_bed | \
    awk 'BEGIN{OFS="\t"} {print $1,"bedtools","gene",$2+1,$3,".",".",".","ID="$4}' \
    > $gff_file

    mkdir -p $output_dir
    run_liftoff \
        $gff_file \
        $genome_file \
        $PATH_TO_DATA/genome_data/GCA_048771995.1_bTaeGut7.mat_genomic.fna \
        $output_dir

    conda activate bedops
    gff2bed < $output_dir/lifted_annotations_genomic.gff > $output_dir/$prefix.bed
    # gff2bed < lifted_annotations_genomic.gff > $prefix.bed

done


# mkdir liftoff_final_cnees_long
# run_liftoff \
#     final_cnees_long.gff \
#     $PATH_TO_DATA/genome_data/GCF_000002315.3_Gallus_gallus-4.0_genomic.fna \
#     $PATH_TO_DATA/genome_data/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.fna \
#     liftoff_final_cnees_long

# input_bed=$PATH_TO_DATA/Sackton_et_al/selected_elements/cnee_overlaps_ATAC_Seq.bed
# bedtools sort -i $input_bed | \
# awk 'BEGIN{OFS="\t"} {print $1,"bedtools","gene",$2+1,$3,".",".",".","ID="NR}' \
# > cnee_overlaps_ATAC_Seq.gff

# mkdir liftoff_cnee_overlaps_ATAC_Seq
# run_liftoff \
#     cnee_overlaps_ATAC_Seq.gff \
#     $PATH_TO_DATA/genome_data/GCF_000002315.3_Gallus_gallus-4.0_genomic.fna \
#     $PATH_TO_DATA/genome_data/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.fna \
#     liftoff_cnee_overlaps_ATAC_Seq
