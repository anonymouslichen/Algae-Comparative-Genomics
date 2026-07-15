#!/bin/bash
#SBATCH --job-name=hyphy_RELAX_paired_Tre
#SBATCH --output=hyphy_RELAX_%A_%a.out
#SBATCH --error=hyphy_RELAX_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100
#SBATCH --array=0-950

# Set directories
HYPHY_DIR="${PROJECT_DIR}/HYPHY"
ALIGN_DIR="$HYPHY_DIR/Alignments/pooled"    #CDS alignments moved from output of 7b-remove_gaps_HYPHY.py
TREE_DIR="$HYPHY_DIR/Labeled_Trees/paired_Tre"
OUT_DIR="$HYPHY_DIR/output/RELAX/paired_Tre"

mkdir -p "$OUT_DIR"

# Get alignment file based on SLURM array index
ALIGN_FILES=($ALIGN_DIR/*.fa)
SEQ_FILE=${ALIGN_FILES[$SLURM_ARRAY_TASK_ID]}
GENE_NAME=$(basename "$SEQ_FILE" _codon_alignment.fa)
TREE_FILE="$TREE_DIR/${GENE_NAME}_labeled.nwk"

# Define output directory for this gene
OUT_FILE=$OUT_DIR/RELAX_${GENE_NAME}.json

# Run HyPhy RELAX
echo 2 | conda run -n hyphy_env hyphy relax --alignment "$SEQ_FILE" --tree "$TREE_FILE" --srv Yes -test Test --reference Reference --output "$OUT_FILE"

