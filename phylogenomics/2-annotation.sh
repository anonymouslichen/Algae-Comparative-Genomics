#!/bin/bash
#SBATCH --job-name=braker_annotation
#SBATCH --output=braker_%A_%a.log
#SBATCH --error=braker_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=0-28         


# Directory containing all the genome files in FASTA format
GENOMES_DIR="${PROJECT_DIR}/SoftMasked"
OUTPUT_DIR="${PROJECT_DIR}/Braker_min_length"

# Create output dir if it doesn't exist
mkdir -p $OUTPUT_DIR

# Set AUGUSTUS_CONFIG_PATH to your writable Augustus config directory
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH}" 

# Set BRAKER_SIF to your BRAKER Singularity image
BRAKER_SIF="${BRAKER_SIF}"

# Get a list of genome files
genomes=( $GENOMES_DIR/processed_*.fasta )

# Select the genome corresponding to the current array task
genome=${genomes[$SLURM_ARRAY_TASK_ID]}

# Get the base name of the genome file
genome_name=$(basename $genome .fasta)

# Create a directory for the current genome's output
mkdir -p $OUTPUT_DIR/$genome_name

# Run BRAKER with the current genome
singularity exec --bind ${PWD}:${PWD} $BRAKER_SIF braker.pl --genome=$genome \
          --species=$genome_name \
          --workingdir=$OUTPUT_DIR/$genome_name \
          --softmasking \
          --gff3 \
          --prot_seq=${PROJECT_DIR}/Viridiplantae.fa \
          --useexisting \
          --threads=8 \
          --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
          --augustus_args="--genemodel=complete" \
          --min_contig=10000
