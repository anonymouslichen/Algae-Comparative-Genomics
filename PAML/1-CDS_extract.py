import os

def get_species_and_gene(protein_header):
    """
    Extract the species prefix and gene identifier from the protein sequence header.
    """
    species_prefix, gene_name = protein_header.strip('>').split('|')
    return species_prefix, gene_name.strip()

def find_coding_sequence(codingseq_file, gene_name):
    """
    Locate the full coding sequence for a given gene name in a braker.codingseq file.
    """
    sequence_lines = []
    record = False

    with open(codingseq_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if record:
                    # If we were recording and hit a new header, stop recording
                    break
                # Start recording if this is the header we're looking for
                if line.strip() == f'>{gene_name}':
                    record = True
            elif record:
                # Collect sequence lines until the next header
                sequence_lines.append(line.strip())

    return ''.join(sequence_lines) if sequence_lines else None

def convert_protein_to_cds(protein_alignment_file, codingseq_base_path, output_dir):
    """
    Converts protein sequences in a given alignment file to their corresponding CDS sequences.
    """
    with open(protein_alignment_file, 'r') as protein_file:
        protein_lines = protein_file.readlines()

    output_lines = []

    for i in range(0, len(protein_lines), 2):
        # Ensure the line is a protein header and there is a subsequent sequence line
        if i+1 < len(protein_lines):
            protein_header = protein_lines[i].strip()
            protein_sequence = protein_lines[i+1].strip()

            # Ensure the line starts with '>'
            if protein_header.startswith('>'):
                # Extract species prefix and gene name from the protein header
                species_prefix, gene_name = get_species_and_gene(protein_header)

                # Define the path to the coding sequence file based on species prefix
                codingseq_file = os.path.join(codingseq_base_path, f'processed_{species_prefix}', 'braker.codingseq')

                if os.path.exists(codingseq_file):
                    # Find the corresponding coding sequence for the gene
                    cds_sequence = find_coding_sequence(codingseq_file, gene_name)

                    if cds_sequence:
                        output_lines.append(f'{protein_header}')
                        output_lines.append(cds_sequence)
                    else:
                        print(f"Warning: No matching CDS found for {protein_header}")
                else:
                    print(f"Warning: Coding sequence file not found for {protein_header} at {codingseq_file}")

    # Define output file name
    output_file_path = os.path.join(output_dir, os.path.basename(protein_alignment_file).replace('.fa', '_cds.fa'))
    with open(output_file_path, 'w') as out_file:
        out_file.write("\n".join(output_lines))

    print(f"CDS conversion for {protein_alignment_file} completed and saved to {output_file_path}")

def process_directory(input_dir, codingseq_base_path, output_dir):
    """
    Processes each protein alignment file in the specified directory.
    """
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith('.fa'):
            protein_alignment_file = os.path.join(input_dir, filename)
            convert_protein_to_cds(protein_alignment_file, codingseq_base_path, output_dir)

# Define paths
input_dir = '/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/Protein_dir'  # Directory containing protein alignment files
codingseq_base_path = '/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/Braker_min_length'
output_dir = '/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/CDS_dir'  # Directory to save CDS output files

# Run the directory processing
process_directory(input_dir, codingseq_base_path, output_dir)
