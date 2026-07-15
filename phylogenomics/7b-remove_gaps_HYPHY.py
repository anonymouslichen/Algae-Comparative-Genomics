import os
import subprocess

def clean_fasta_headers(input_path, cleaned_path):
    """
    Remove everything after and including '|' in FASTA headers.
    """
    with open(input_path, "r") as infile, open(cleaned_path, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                clean_header = line.split("|")[0].strip()
                outfile.write(clean_header + "\n")
            else:
                outfile.write(line)

def run_trimal(input_dir, output_dir, trimal_path="trimal", params="-gappyout -fasta"):
    """
    Clean headers and run trimAl on all alignment files in the input directory.
    """
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith((".fa")):  
            input_file = os.path.join(input_dir, filename)
            cleaned_file = os.path.join(output_dir, f"cleaned_{filename}")
            output_file = os.path.join(output_dir, filename)

            # Clean FASTA headers
            clean_fasta_headers(input_file, cleaned_file)

            # Run trimAl
            command = [trimal_path, "-in", cleaned_file, "-out", output_file] + params.split()
            try:
                print(f"Running: {' '.join(command)}")
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while processing {input_file}: {e}")

            # Remove the cleaned file after trimAl
            os.remove(cleaned_file)

# Define your input and output directories 
PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")
input_dir = os.path.join(PROJECT_DIR, "CDS/pal2nal_dir")
output_dir = os.path.join(PROJECT_DIR, "CDS/trimmed_alignments_HYPHY")
    
# Run trimAl on the alignments
run_trimal(input_dir, output_dir)