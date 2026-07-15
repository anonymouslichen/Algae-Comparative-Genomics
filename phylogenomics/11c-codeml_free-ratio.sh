#!/bin/bash
#SBATCH --job-name=codeml_job_pooled_M4
#SBATCH --output=codeml_output_%A_%a.log
#SBATCH --error=codeml_error_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --mem=100
#SBATCH --array=0-950

# Define paths
PAML_DIR="${PROJECT_DIR}/PAML"
ALIGN_DIR="${PAML_DIR}/Alignments"
TREE_DIR="${PAML_DIR}/Gene_Trees"
OUT_DIR="${PAML_DIR}/output"

# Load the necessary module
module load paml

# Get the list of alignment files and select one based on SLURM_ARRAY_TASK_ID
ALIGN_FILES=($ALIGN_DIR/*.phy)
SEQ_FILE=${ALIGN_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the gene name from the file name
GENE_NAME=$(basename $SEQ_FILE _codon_alignment.phy)

# Grab tree files
TREE_FILE=${TREE_DIR}/${GENE_NAME}_aligned_M0.txt

# Define the output file name and unique control file name for each array job
OUTPUT_FILE=Trebouxiophyceae_${GENE_NAME}_M4.out
CTL_FILE=${PAML_DIR}/Tre_${GENE_NAME}_M4


echo "Writing control file to: $CTL_FILE"
echo "Using alignment: $SEQ_FILE"
echo "Using tree: $TREE_FILE"

# Create a unique codeml control file for the current gene
cat > $CTL_FILE << EOL
      seqfile = $SEQ_FILE
      treefile = $TREE_FILE
      outfile = $OUT_DIR/$OUTPUT_FILE

      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
      CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

      ndata = 1 * number of gene alignments to be analysed
      clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis

      model = 1 * models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

      icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

      fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
      omega = .4 * initial or fixed omega, for codons or codon-based AAs

      cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

EOL

# Run codeml using the unique control file
codeml $CTL_FILE

# Check if codeml ran successfully
if [ $? -eq 0 ]; then
  echo "codeml ran successfully for $SEQ_FILE. Output is in $OUT_DIR/$OUTPUT_FILE"
else
  echo "codeml failed to run for $SEQ_FILE. Check your input files and control file."
fi

# Clean up by removing the unique control file after completion
rm -f $CTL_FILE
