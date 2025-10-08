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

# Define working directory 
WORKDIR=/path/to/workingdir

# Create a results directory if it doesn’t exist
mkdir -p ${ORTHO_INPUT}

# Move into working directory
cd ${WORKDIR}

# Run OrthoFinder
orthofinder -f ${ORTHO_INPUT} -t ${SLURM_CPUS_PER_TASK}
