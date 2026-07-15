#!/bin/bash
#SBATCH --job-name=soft_masking
#SBATCH --output=soft_masking_%j.out
#SBATCH --error=soft_masking_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-27

# Load repeatmasker
module load trf phrap rmblast repeatmasker/4.1.1
module load python3

# Directory containing the genome files (set this to your directory)
GENOME_DIR="${PROJECT_DIR}/Genomes"

# Create an output directory for the softmasked genomes
OUTPUT_DIR="${PROJECT_DIR}/SoftMasked"
mkdir -p $OUTPUT_DIR

# Get the list of genome files
GENOME_FILES=($(ls $GENOME_DIR/*.fasta))

# Get the genome file based on SLURM_ARRAY_TASK_ID
GENOME_FILE=${GENOME_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the genome filename without the path
GENOME_NAME=$(basename $GENOME_FILE)

# Temporary file to store the processed genome with adjusted headers
PROCESSED_GENOME="$OUTPUT_DIR/processed_${GENOME_NAME}"

# Clean FASTA headers (replace spaces with underscores)
echo "Processing headers for ${GENOME_NAME}"
sed '/^>/ s/ /_/g' "$GENOME_FILE" > "$PROCESSED_GENOME"

# Run RepeatMasker on the processed genome file with softmasking enabled (-xsmall)
echo "Running RepeatMasker on ${GENOME_NAME}"
RepeatMasker -pa 4 -xsmall -dir $OUTPUT_DIR "$PROCESSED_GENOME"

echo "Finished processing ${GENOME_NAME}"
