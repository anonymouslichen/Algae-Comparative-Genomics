#!/bin/bash

#SBATCH --job-name=iqtree_analysis
#SBATCH --output=iqtree_%j.out
#SBATCH --error=iqtree_%j.err
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-user=youemail

# Load iqtree
module load iqtree2

# Set paths
ALIGNMENT_DIR="path/to/alignmentdir"
OUTPUT_DIR="path/to/genetreedir"
TMP_DIR="path/to/tmpdir"
TRIMAL_PATH="path/to/trimal"

# Create output and temp directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

# 1. Clean FASTA headers
for file in "$ALIGNMENT_DIR"/*.fa; do
    cleaned="$TMP_DIR/$(basename "$file")"
    sed 's/|.*$//' "$file" > "$cleaned"
done

# 2. Trim alignments and build trees
for file in "$TMP_DIR"/*.fa; do
    filename=$(basename "$file" .fa)
    trimmed="$TMP_DIR/${filename}_trimmed.fa"

    echo "Processing $filename..."
    "$TRIMAL_PATH" -in "$file" -out "$trimmed" -gappyout

    iqtree2 -s "$trimmed" -st AA -m LG+F+R7 -nt AUTO -pre "$OUTPUT_DIR/${filename}"

    echo "Gene tree generated: $OUTPUT_DIR/${filename}.treefile"
done
