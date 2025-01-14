import os
import subprocess

# Function to run MUSCLE alignment on protein sequences
def run_muscle(input_file, output_file):
    try:
        subprocess.run(["muscle", "-in", input_file, "-out", output_file], check=True)
        print(f"MUSCLE alignment created: {output_file}")
        return output_file  # Return the path to the MUSCLE alignment
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE on {input_file}: {e}")
        return None

# Function to reorder the sequences alphabetically in the alignment file
def reorder_sequences_alphabetically(alignment_file):
    with open(alignment_file, 'r') as f:
        lines = f.readlines()

    # Split into header and sequence
    sequences = []
    current_header = None
    current_sequence = []

    for line in lines:
        if line.startswith(">"):
            if current_header:  # Save previous sequence
                sequences.append((current_header, ''.join(current_sequence)))
            current_header = line.strip()
            current_sequence = []
        else:
            current_sequence.append(line.strip())

    # Append the last sequence
    if current_header:
        sequences.append((current_header, ''.join(current_sequence)))

    # Sort the sequences alphabetically by header
    sequences.sort(key=lambda x: x[0])

    # Write the sorted sequences back to the file
    with open(alignment_file, 'w') as f:
        for header, sequence in sequences:
            f.write(f"{header}\n{sequence}\n")
    print(f"Sequences reordered alphabetically in {alignment_file}")

# Function to run pal2nal
def run_pal2nal(protein_file, nucleotide_file, output_file):
    pal2nal_path = "/home/stan0477/meye2099/bin/pal2nal.v14/pal2nal.pl"  # Full path to pal2nal executable
    try:
        with open(output_file, 'w') as out_f:
            subprocess.run([pal2nal_path, protein_file, nucleotide_file, "-output", "fasta"], check=True, stdout=out_f)
        print(f"pal2nal alignment created: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running pal2nal on {protein_file} and {nucleotide_file}: {e}")

# Directories containing protein and CDS files
protein_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/Protein_dir"
cds_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/CDS_dir"
muscle_alignment_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/muscle_alignment_dir"
output_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/pal2nal_dir"

# Ensure the output directories exist
os.makedirs(muscle_alignment_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# List all protein files
protein_files = [f for f in os.listdir(protein_dir) if f.endswith(".fa") or f.endswith(".fasta")]

# Loop through each protein file
for protein_file in protein_files:
    protein_file_path = os.path.join(protein_dir, protein_file)

    # Run MUSCLE alignment on the protein file
    aligned_protein_file = os.path.join(muscle_alignment_dir, f"{os.path.splitext(protein_file)[0]}_aligned.fa")
    aligned_protein_file = run_muscle(protein_file_path, aligned_protein_file)

    if not aligned_protein_file:
        continue  # Skip to next file if MUSCLE alignment failed

    # Reorder the sequences alphabetically in the MUSCLE alignment file
    reorder_sequences_alphabetically(aligned_protein_file)

    # Construct the corresponding CDS file name
    cds_file = f"{os.path.splitext(protein_file)[0]}_cds.fa"
    cds_file_path = os.path.join(cds_dir, cds_file)

    # Check if the corresponding CDS file exists
    if not os.path.isfile(cds_file_path):
        print(f"Error: CDS file {cds_file_path} not found for {protein_file}. Skipping...")
        continue

    # Define the output file path for the codon alignment
    output_file_path = os.path.join(output_dir, f"{os.path.splitext(protein_file)[0]}_codon_alignment.fa")

    # Run pal2nal to generate the codon alignment
    run_pal2nal(aligned_protein_file, cds_file_path, output_file_path)
