#!/usr/bin/env python3

import re
import os

# Mapping for tree labeling
# Format: "Taxon_name": (Number_label, "#1" or None)
taxa_mapping = {
    "Ast_eri": (1, "#1"),
    "Coc_sub": (2, None),
    "Coc_vir": (3, "#1"),
    "Myr_bis": (4, None),
    "Sym_irr": (5, None),
    "Sym_ret": (6, "#1"),
    "Tre_spC0010": (7, "#1")
}

# Remove branch lengths and rename taxa based on the mapping
def process_tree(tree: str) -> str:
    tree = re.sub(r":\d+\.\d+(?:[eE][-+]?\d+)?", "", tree)

    # Replace taxon names with assigned numbers and labels
    for taxa, (number, label) in taxa_mapping.items():
        replacement = f"{number} {label}" if label else f"{number}"
        tree = re.sub(fr"\b{taxa}\b", replacement, tree)

    return tree

# Input and output directories
tree_dir = "path/to/your/treedir"
processed_dir = "path/to/outdir"
os.makedirs(processed_dir, exist_ok=True)

# Process all tree files 
for filename in os.listdir(tree_dir):
    if filename.endswith(".treefile"):
        filepath = os.path.join(tree_dir, filename)
        with open(filepath, "r") as file:
            tree = file.read().strip()

        processed_tree = process_tree(tree)

        processed_filename = f"{os.path.splitext(filename)[0]}_M2.txt"  #Adjust depending on which model you're creating trees for, e.g. M0 (null model), M2 (two-ratio model)
        processed_filepath = os.path.join(processed_dir, processed_filename)

        with open(processed_filepath, "w") as file:
            file.write(processed_tree)

        print(f"Processed: {filename}")

print(f"\nAll processed trees saved to:\n{processed_dir}")
