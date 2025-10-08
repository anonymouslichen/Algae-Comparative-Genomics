#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=orthofinder_%A.out
#SBATCH --error=orthofinder_%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail

# Load OrthoFinder module 
module load orthofinder

# Define input directory 
INPUT_DIR=/path/to/inputdir #Directory with output from step 3

# Run OrthoFinder
orthofinder -f ${INPUT_DIR} -t ${SLURM_CPUS_PER_TASK}
