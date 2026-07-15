#!/bin/bash
#SBATCH --job-name=busco_analysis
#SBATCH --output=busco_%A_%a.out
#SBATCH --error=busco_%A_%a.err
#SBATCH --array=0-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G               
#SBATCH --time=1:00:00      

# Load modules 
module load busco
module load python3/3.10.9_anaconda2023.03_libmamba

# Variables
BRAKER_DIR="${PROJECT_DIR}/Braker_min_length"                   
OUTPUT_DIR="${BRAKER_DIR}/busco_outputs_protein"       
LINEAGE_DB="chlorophyta_odb10"        
MODE="protein"                            
SAMPLE_LIST="${BRAKER_DIR}/busco_assemblies.txt"   # File listing the 7 directories to process

# Read sample directories into array
mapfile -t SAMPLES < ${SAMPLE_LIST}
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Build path to braker.aa
PROTEIN_FILE="${BRAKER_DIR}/${SAMPLE}/braker.aa"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Run BUSCO
busco -i ${PROTEIN_FILE} \
      -o ${SAMPLE}_busco \
      -l ${LINEAGE_DB} \
      -m ${MODE} \
      --out_path ${OUTPUT_DIR} \
      --cpu 16