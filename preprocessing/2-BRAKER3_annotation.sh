#!/bin/bash
#SBATCH --job-name=braker_annotation
#SBATCH --output=braker_%A_%a.out
#SBATCH --error=braker_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-##                  #Replace ## with number of genomes - 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail                

# Move to working directory (called upon later to mount to your Singularity container)
cd /path/to/workingdir

# Directory containing all the genome files in FASTA format
GENOMES_DIR=/path/to/genomedir
OUTPUT_DIR=/path/to/outputdir
export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config          #You will need to copy the config folder from within the singularity container to someplace writeable

# Path to BRAKER Singularity image
BRAKER_SIF=/path/to/bin/braker3.sif

# Get a list of genome files
genomes=( $GENOMES_DIR/processed_*.fna )

# Select the genome corresponding to the current array task
genome=${genomes[$SLURM_ARRAY_TASK_ID]}
genome_name=$(basename $genome .fna)

# Create a directory for the current genome's output
mkdir -p $OUTPUT_DIR/$genome_name -p

# Run BRAKER with the current genome
singularity exec --bind ${PWD}:${PWD} $BRAKER_SIF braker.pl --genome=$genome \
          --species=$genome_name \
          --workingdir=$OUTPUT_DIR/$genome_name \
          --softmasking \
          --gff3 \
          --prot_seq=Viridiplantae.fa \          #Add your desired protein file: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/
          --useexisting \
          --threads=8 \
          --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH
