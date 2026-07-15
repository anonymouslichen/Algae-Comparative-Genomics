#!/bin/bash
#SBATCH --job-name=iqtree_analysis  
#SBATCH --output=iqtree_%A_%a.out   
#SBATCH --error=iqtree_%A_%a.err    
#SBATCH --time=5:00:00             
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=4           
#SBATCH --mem=16G                   
#SBATCH --array=0-950

# Load the IQ-TREE and trimAI modules
module load iqtree2

# Define input directory with alignment files and output directory for gene trees
ALIGNMENT_DIR="${PROJECT_DIR}/IQ-TREE/alignments"     #MUSCLE alignments moved from step 6-align_CDS.py
OUTPUT_DIR="${PROJECT_DIR}/IQ-TREE/gene_trees"
TMP_DIR="${PROJECT_DIR}/IQ-TREE/modified_fasta"

# Create a temporary directory to store modified files
mkdir -p $TMP_DIR

# Get list of files
files=( $ALIGNMENT_DIR/*.fa )

# Get current file based on array index
file=${files[$SLURM_ARRAY_TASK_ID]}

#Get basename of file
filename=$(basename $file .fa)

# Use sed to remove everything after and including "|"
sed 's/|.*$//' "$file" > "$TMP_DIR/$filename.fa"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run trimAI to clean the alignment
trimal -in "$TMP_DIR/$filename.fa" -out "$TMP_DIR/${filename}_trimmed.fa" -gappyout

# Run IQ-TREE on the trimmed alignment file
iqtree2 -s "$TMP_DIR/${filename}_trimmed.fa"  -st AA -m LG+F+G4 -nt AUTO -pre "$OUTPUT_DIR/${filename}"

echo "Protein gene tree generated: $OUTPUT_DIR/${filename}.treefile"

