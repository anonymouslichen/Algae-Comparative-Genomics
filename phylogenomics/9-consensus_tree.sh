#!/bin/bash
#SBATCH --job-name=sumtrees_consensus 
#SBATCH --output=sumtrees_%j.out   
#SBATCH --error=sumtrees_%j.err     
#SBATCH --time=4:00:00              
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=8             
#SBATCH --mem=16G                        


# Load Python / DendroPy / SumTrees
module load python3/3.10.9_anaconda2023.03_libmamba 
source activate genome_env  

# Directory of gene trees
GENE_TREE_DIR="${PROJECT_DIR}/IQ-TREE/gene_trees"

sumtrees.py \
    -s consensus \
    -f 0.5 \
    -e mean-length \
    -p \
    -d 0 \
    -M \
    --force-unrooted \
    -o "$GENE_TREE_DIR/species_tree_consensus.tre" \
    "$GENE_TREE_DIR"/*.treefile

echo "Consensus species tree written to: $GENE_TREE_DIR/species_tree_consensus.tre"