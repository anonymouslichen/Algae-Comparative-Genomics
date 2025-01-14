import os

def rewrite_fasta_single_line(input_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    header = None
    sequence = []
    output_lines = []

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                # Write the previous header and sequence
                output_lines.append(header + '\n' + ''.join(sequence) + '\n')
            # Start a new header and sequence
            header = line
            sequence = []
        else:
            # Append the sequence lines
            sequence.append(line)

    # Write the last sequence
    if header is not None:
        output_lines.append(header + '\n' + ''.join(sequence) + '\n')

    # Overwrite the original file with the new format
    with open(input_file, 'w') as outfile:
        outfile.writelines(output_lines)

def process_fasta_files_in_directory(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".fa"):
            file_path = os.path.join(directory, filename)
            rewrite_fasta_single_line(file_path)
            print(f"Processed: {filename}")

# Example usage
directory = '/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/Protein_dir'
process_fasta_files_in_directory(directory)
