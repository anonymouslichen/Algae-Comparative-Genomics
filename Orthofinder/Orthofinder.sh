#!/bin/bash                   

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g             
#SBATCH --time=1:00:00       
#SBATCH --mail-type=ALL       
#SBATCH --output=job_output.log  
#SBATCH --error=job_error.log   
#SBATCH --mail-user=

# Load Orthofinder module
module load orthofinder       

# Move to working directory
cd /home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/  

# Run OrthoFinder with specified input folder and multiple sequence alignment option
orthofinder -f ./Orthofinder/ -M msa  
