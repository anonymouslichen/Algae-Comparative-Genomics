import re
import os
  
# Taxa -> PAML sequence numbers (must match the order in your codeml alignment)
taxa_numbers = {
    "Ast_eri": 1, "Coc_sub": 2, "Coc_vir": 3, "Myr_bis": 4,
    "Sym_irr": 5, "Sym_ret": 6, "Tre_spC0010": 7,
}

# Foreground labeling schemes: taxon -> branch label.
# "M0" is empty = no labels (used by BOTH the null and free-ratio models).
schemes = {
    "M0": {},
    "M2": {         # two-ratio: mark lichen lineages as foreground
    "Ast_eri": "#1", "Coc_vir": "#1", "Sym_ret": "#1", "Tre_spC0010": "#1",
    },
}
  
def process_tree(tree, labels):
    # Remove branch lengths
    tree = re.sub(r":\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", "", tree)
    # Rename taxa to numbers, adding a branch label where this scheme specifies one
    for taxon, number in taxa_numbers.items():
        label = labels.get(taxon)
        replacement = f"{number} {label}" if label else f"{number}"
        tree = re.sub(fr"\b{taxon}\b", replacement, tree)
    return tree

PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")
tree_dir = os.path.join(PROJECT_DIR, "IQ-TREE/gene_trees")
processed_dir = os.path.join(PROJECT_DIR, "PAML/Gene_Trees")
os.makedirs(processed_dir, exist_ok=True)

for filename in os.listdir(tree_dir):
    if not filename.endswith(".treefile"):
        continue
    with open(os.path.join(tree_dir, filename)) as f:
        tree = f.read().strip()
    base = os.path.splitext(filename)[0]
    for scheme_name, labels in schemes.items():
        out = os.path.join(processed_dir, f"{base}_{scheme_name}.txt")
        with open(out, "w") as f:
            f.write(process_tree(tree, labels))

print(f"Processed trees saved to {processed_dir}")
