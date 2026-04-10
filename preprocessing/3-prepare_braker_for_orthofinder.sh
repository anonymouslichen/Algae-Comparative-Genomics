#!/bin/bash
#SBATCH --job-name=prepare_braker_orthofinder
#SBATCH --output=prepare_braker_%j.out
#SBATCH --error=prepare_braker_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=youremail

# Define directories
BRAKER_DIR=/path/to/brakerdir          
OUTPUT_DIR=/path/to/outputdir          

# Create output directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Define subdirectories to include
INCLUDE_SUBDIRS=(
    "processed_Myr_bis" "processed_Coc_pri" "processed_Coc_sub"
    "processed_Tre_spC0004" "processed_Tre_spC0005" "processed_Tre_spC0009"
    "processed_Tre_spC0006" "processed_Tre_spTZW2008" "processed_Tre_spC0007"
    "processed_Tre_spC0010" "processed_Coc_vir" "processed_Tre_lyn"
    "processed_Ast_eri" "processed_Ast_glo" "processed_Ast_spCNOR1"
    "processed_Sym_ret" "processed_Tre_spA1-2" "processed_Sym_irr"
)

# Copy and rename braker.aa files 
for subdir_name in "${INCLUDE_SUBDIRS[@]}"; do
    subdir="$BRAKER_DIR/$subdir_name"

    if [ -d "$subdir" ]; then
        species=$(echo "$subdir_name" | sed 's/processed_//')

        if [ -f "$subdir/braker.aa" ]; then
            new_filename="${species}_braker.fa"
            cp "$subdir/braker.aa" "$OUTPUT_DIR/$new_filename"
            echo "Copied $subdir/braker.aa → $OUTPUT_DIR/$new_filename"
        else
            echo "No braker.aa found in $subdir"
        fi
    else
        echo "Skipping $subdir_name (directory missing)"
    fi
done

# Edit FASTA headers to include genome name as prefix (e.g., "species|gene_id") 
python3 <<PYCODE
import os

output_dir = "${OUTPUT_DIR}"

def edit_braker_headers(filepath, prefix):
    """Add the genome name prefix to each FASTA header line."""
    with open(filepath, 'r') as infile:
        lines = infile.readlines()
    with open(filepath, 'w') as outfile:
        for line in lines:
            if line.startswith('>'):
                line = f">{prefix}|{line[1:].strip()}\n"
            outfile.write(line)

def edit_headers_in_directory(directory):
    for filename in os.listdir(directory):
        if filename.endswith('_braker.fa'):
            filepath = os.path.join(directory, filename)
            prefix = filename.split('_braker.fa')[0]
            edit_braker_headers(filepath, prefix)
            print(f"Edited headers in {filename}")

edit_headers_in_directory(output_dir)
PYCODE

echo "All files processed. Output saved in $OUTPUT_DIR"
