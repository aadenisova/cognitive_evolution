#!/bin/bash
#SBATCH --job-name=extract_maf
#SBATCH --output=logs/extract_maf_%A_%a.out
#SBATCH --error=logs/extract_maf_%A_%a.err
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=hpc_l40_a    
#SBATCH --time=24:00:00

PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/conserved_regions/
PATH_TO_WGS=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/vgpBirds_137_way

# SLURM_ARRAY_TASK_ID=1
# prefixes=("final_cnees_long_test")
# "cnee_overlaps_ATAC_Seq" is included in final_cnees_long

prefixes=("final_cnees_long")
for prefix in "${prefixes[@]}"; do
    BED_FILE=$PATH_TO_DATA/Sackton_et_al/liftoff_$prefix/lifted_annotations_genomic.bed

    if [ ! -f $BED_FILE ]; then
        awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4, $5, $9}' \
            $PATH_TO_DATA/Sackton_et_al/liftoff_$prefix/lifted_annotations_genomic.gff \
            > $BED_FILE
    fi

    NUM_LINES=$(wc -l < $BED_FILE)
    CHUNK_SIZE=$(( (NUM_LINES + 100 - 1)/100 ))

    START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
    END=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
    if [ $END -gt $NUM_LINES ]; then
        END=$NUM_LINES
    fi

    echo "Task $SLURM_ARRAY_TASK_ID processing lines $START to $END for prefix $prefix"

    mkdir -p $PATH_TO_DATA/Sackton_et_al/maf_by_$prefix
    mkdir -p $PATH_TO_DATA/Sackton_et_al/rer_trees_by_$prefix
    mkdir -p $PATH_TO_DATA/Sackton_et_al/phyloP_$prefix

    sed -n "${START},${END}p" $BED_FILE | while read -r LINE; do
        CHROM=$(echo $LINE | awk '{print $1}')
        START_POS=$(echo $LINE | awk '{print $2}')
        STOP_POS=$(echo $LINE | awk '{print $3}')
        ID=$(echo $LINE | awk '{print $4}'| sed -E 's/ID=([^;]+).*/\1/')

        echo "Extracting MAF for $ID: $CHROM:$START_POS-$STOP_POS"
        echo "$PATH_TO_WGS/maf_by_chr_GCF_016700215.2/2.$CHROM.maf"

        if [ -f $PATH_TO_DATA/Sackton_et_al/maf_by_$prefix/${ID}.maf ]; then
            echo "MAF for $ID already exists, skipping extraction"
            continue
        fi

        mafExtractor --maf $PATH_TO_WGS/maf_by_chr_GCF_016700215.2/2.$CHROM.maf \
            --seq GCF_016700215.2.$CHROM \
            --start $START_POS \
            --stop $STOP_POS \
            > $PATH_TO_DATA/Sackton_et_al/maf_by_$prefix/${ID}.maf

        # element_name=$ID

        # TREE_FILE=$PATH_TO_WGS/roadies_v1.1.16b_renamed.nwk
        # alignment=$PATH_TO_DATA/Sackton_et_al/maf_by_final_cnees_long/$element_name.maf

        # root=$PATH_TO_DATA/Sackton_et_al/phyloP_$prefix/$element_name
        # phyloFit \
        #     --tree $TREE_FILE \
        #     --subst-mod REV \
        #     --EM \
        #     --out-root $root \
        #     --precision HIGH \
        #     $alignment

        # echo -e "$element_name\t$(tail -n 1 $root.mod | sed 's/^TREE: //')" > $PATH_TO_DATA/Sackton_et_al/rer_trees_by_$prefix/$element_name.tree
        

        element_name=$ID

        TREE_FILE=$PATH_TO_WGS/roadies_v1.1.16b_birds_only_renamed_no_na.nwk
        alignment=$PATH_TO_DATA/Sackton_et_al/maf_by_final_cnees_long/$element_name.maf

        if [ ! -f $alignment ] ; then
            echo "Alignment file $alignment not found, skipping $element_name"
            continue
        fi

        root=$PATH_TO_DATA/Sackton_et_al/phyloP_birds_no_na_$prefix/$element_name
        phyloFit \
            --tree $TREE_FILE \
            --subst-mod REV \
            --EM \
            --out-root $root \
            --precision HIGH \
            $alignment

        echo -e "$element_name\t$(tail -n 1 $root.mod | sed 's/^TREE: //')" > $PATH_TO_DATA/Sackton_et_al/rer_trees_birds_no_na_by_$prefix/$element_name.tree
        
    done
done