#!/bin/bash
#SBATCH --job-name=hyphy_tree_label_paired_Symbiochloris
#SBATCH --output=hyphy_label_tree_%A_%a.out
#SBATCH --error=hyphy_label_tree_%A_%a.err
#SBATCH --time=0:05:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100
#SBATCH --array=0-950

# Set directories
HYPHY_DIR="${PROJECT_DIR}/HYPHY"
TREE_DIR="${PROJECT_DIR}/IQ-TREE/gene_trees"
OUT_DIR="$HYPHY_DIR/Labeled_Trees/paired_Sym"

mkdir -p "$OUT_DIR"

# Taxon names — adjust these per comparison
TEST_TAXON="Sym_ret"
REF_TAXON="Sym_irr"

TREE_FILES=($TREE_DIR/*.treefile)
TREE_FILE=${TREE_FILES[$SLURM_ARRAY_TASK_ID]}
GENE_NAME=$(basename "$TREE_FILE" _aligned.treefile)
OUT_FILE=$OUT_DIR/${GENE_NAME}_labeled.nwk

# Annotate the tree
sed -e "s/${TEST_TAXON}/${TEST_TAXON}{Test}/" \
    -e "s/${REF_TAXON}/${REF_TAXON}{Reference}/" \
    "$TREE_FILE" > "$OUT_FILE"

