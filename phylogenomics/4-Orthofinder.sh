#!/bin/bash
#SBATCH --mem=32g             
#SBATCH --time=4:00:00        
#SBATCH --cpus-per-task=16     
#SBATCH --output=orthofinder_%A.out
#SBATCH --error=orthofinder_%A.err

module load orthofinder

cd "${PROJECT_DIR}"

orthofinder -f ./Orthofinder/ -t 16 
