#!/bin/bash

#SBATCH --job-name=iqtree_analysis  # Job name
#SBATCH --output=iqtree_%j.out      # Standard output log
#SBATCH --error=iqtree_%j.err       # Standard error log
#SBATCH --time=48:00:00             # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem=16G                   # Memory per node
#SBATCH --mail-user=meye2099@umn.edu

# Load the IQ-TREE and trimAI modules
module load iqtree2

# Define input directory with alignment files and output directory for gene trees
ALIGNMENT_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/IQ-TREE/alignments #output alignmentws from 3-align_concat.py
OUTPUT_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/IQ-TREE/gene_trees
TMP_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/IQ-TREE/modified_fasta

# Create a temporary directory to store modified files
mkdir -p $TMP_DIR

# Loop through all `.fa` files in the current directory
for file in $ALIGNMENT_DIR/*.fa; do
  # Use sed to remove everything after and including "|"
  sed 's/|.*$//' "$file" > "$TMP_DIR/$(basename "$file")"
done

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each alignment file in the input directory
for file in $TMP_DIR/*.fa; do
  # Extract the filename without the directory and extension
  filename=$(basename "$file" .fa)

  # Run trimAI to clean the alignment
  /home/stan0477/meye2099/bin/trimal/source/trimal -in "$file" -out "$TMP_DIR/${filename}_trimmed.fa" -gappyout

  # Run IQ-TREE on the trimmed alignment file
  iqtree2 -s "$TMP_DIR/${filename}_trimmed.fa" -st AA -m LG+F+R7 -nt AUTO -pre "$OUTPUT_DIR/${filename}"

  echo "Protein gene tree generated: $OUTPUT_DIR/${filename}.treefile"
done
