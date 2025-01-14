#!/bin/bash
#SBATCH --job-name=soft_masking
#SBATCH --output=soft_masking_%j.out
#SBATCH --error=soft_masking_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meye2099@umn.edu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-35          #change to number of genomes - 1

# Load repeatmasker
module load trf phrap rmblast repeatmasker/4.1.1
module load python3

# Directory containing the genome files (set this to your directory)
GENOME_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/Genomes

# Create an output directory for the softmasked genomes
OUTPUT_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/SoftMasked
mkdir -p $OUTPUT_DIR

# Get the list of genome files
GENOME_FILES=($(ls $GENOME_DIR/*.fna))

# Get the genome file based on SLURM_ARRAY_TASK_ID
GENOME_FILE=${GENOME_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the genome filename without the path
GENOME_NAME=$(basename $GENOME_FILE)

# Temporary file to store the processed genome with adjusted headers
PROCESSED_GENOME="$OUTPUT_DIR/processed_${GENOME_NAME}"

# Run the Python script to process the headers (replace spaces with underscores)
echo "Processing headers for ${GENOME_NAME}"
python3 /home/stan0477/meye2099/bin/remove_whitespace.py "$GENOME_FILE" "$PROCESSED_GENOME"

# Run RepeatMasker on the processed genome file with softmasking enabled (-xsmall)
echo "Running RepeatMasker on ${GENOME_NAME}"
RepeatMasker -pa 4 -xsmall -dir $OUTPUT_DIR "$PROCESSED_GENOME"

echo "Finished processing ${GENOME_NAME}"
