#!/bin/bash
#SBATCH --job-name=interproscan
#SBATCH --output=interproscan_%A_%a.out
#SBATCH --error=interproscan_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meye2099@umn.edu
#SBATCH --array=0-950

# Define paths
WORK_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/Orthofinder/OrthoFinder/Results_May22/Single_Copy_Orthologue_Sequences
OUT_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/InterProScan

# Load the necessary module
module load interproscan

# Get the list of alignment files and select one based on SLURM_ARRAY_TASK_ID
SOG_FILES=($WORK_DIR/*.fa)
SOG_FILE=${SOG_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the gene name from the file name
SOG_NAME=$(basename $SOG_FILE .fa)

# make SOG_TMP tmp file name
SOG_TMP=$(mktemp /tmp/${SOG_NAME}_cleaned.XXXXXX.fa)

# Remove * from seq file, save as tmp file
sed 's/\*//g' $SOG_FILE > $SOG_TMP

# Define the output file name for each array job
OUTPUT_FILE=interproscan_${SOG_NAME}.tsv

# Run interproscan on SOGs 
interproscan.sh -f TSV -i $SOG_TMP -o $OUT_DIR/$OUTPUT_FILE -T /tmp

# Remove tmp file
rm $SOG_TMP
