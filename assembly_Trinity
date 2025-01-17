#!/bin/bash
#SBATCH --job-name=trinity_assembly
#SBATCH --output=trinity_assembly_%A_%a.out
#SBATCH --error=trinity_assembly_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-11              
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meye2099@umn.edu

# Change to working directory
cd /home/stan0477/meye2099

# Define the Singularity image for DRAP
SINGULARITY_IMAGE=bin/trinityrnaseq_latest.sif

# Define the directory containing the transcriptomes
TRANSCRIPTOME_DIR=Algae_Evolution/Trebouxiophyceae/Transcriptomes/fastq/renamed

# Define output directory
OUTPUT_DIR=/scratch.global/meye2099

# Get the list of paired-end and single-end FASTQ files
PAIRED_FILES=($(ls ${TRANSCRIPTOME_DIR}/*renamed_pass_1.fastq.gz | sort))
SINGLE_FILES=($(ls ${TRANSCRIPTOME_DIR}/*renamed_pass.fastq.gz | sort))

# Combine the lists for array job indexing
ALL_FILES=("${PAIRED_FILES[@]}" "${SINGLE_FILES[@]}")

# Get the current FASTQ file
FASTQ=${ALL_FILES[$SLURM_ARRAY_TASK_ID]}

# Determine if the file is paired-end or single-end
if [[ $FASTQ == *_pass_1.fastq.gz ]]; then
    # Paired-end processing
    FASTQ1=$FASTQ
    FASTQ2=${FASTQ1/_pass_1.fastq.gz/_pass_2.fastq.gz}
    TRANSCRIPTOME="trinity_$(basename ${FASTQ1} _renamed_pass_1.fastq.gz)"

    # Run Trinity for paired-end
    singularity exec --bind ${PWD}:${PWD},${OUTPUT_DIR}:${OUTPUT_DIR} $SINGULARITY_IMAGE Trinity \
    --seqType fq --max_memory 50G --left ${FASTQ1} --right ${FASTQ2} --output ${OUTPUT_DIR}/${TRANSCRIPTOME}


else
    # Single-end processing
    FASTQ_SINGLE=$FASTQ
    TRANSCRIPTOME="trinity_$(basename ${FASTQ_SINGLE} _renamed_pass.fastq.gz)"

    # Run DRAP for single-end
    singularity exec --bind ${PWD}:${PWD},${OUTPUT_DIR}:${OUTPUT_DIR} $SINGULARITY_IMAGE Trinity \
    --seqType fq --max_memory 50G --single ${FASTQ_SINGLE} --output ${OUTPUT_DIR}/${TRANSCRIPTOME}
fi
