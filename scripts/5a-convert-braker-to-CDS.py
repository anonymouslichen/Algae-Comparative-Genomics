#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO

# Find the matching CDS sequence in a BRAKER codingseq file
def extract_braker_sequence(filepath, query):
    for record in SeqIO.parse(filepath, "fasta"):
        if record.id == query:
            return str(record.seq)
    return None

# Convert BRAKER protein FASTA files to corresponding CDS sequences
def convert_protein_to_cds(input_dir, output_dir, braker_base_dir):
    os.makedirs(output_dir, exist_ok=True)

    for protein_file in os.listdir(input_dir):
        if not protein_file.endswith(".fa"):
            continue

        protein_path = os.path.join(input_dir, protein_file)
        output_file = os.path.join(output_dir, f"{os.path.splitext(protein_file)[0]}_cds.fa")

        with open(protein_path, "r") as prot_fh, open(output_file, "w") as out_fh:
            for record in SeqIO.parse(prot_fh, "fasta"):
                header = record.id
                sequence = None

                try:
                    species, gene_info = header.split("|")
                except ValueError:
                    print(f"[WARNING] Skipping malformed header: {header}")
                    continue

                sub_dir = os.path.join(braker_base_dir, f"processed_{species}")
                file_name = "braker.codingseq"
                full_path = os.path.join(sub_dir, file_name)

                if os.path.exists(full_path):
                    sequence = extract_braker_sequence(full_path, gene_info)
                else:
                    print(f"[WARNING] File not found for species {species}: {full_path}")

                if sequence:
                    out_fh.write(f">{header}\n{sequence}\n")
                else:
                    print(f"[INFO] No sequence found for header: {header}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert BRAKER protein FASTA files into CDS FASTA files."
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing input protein FASTA files")
    parser.add_argument("--output_dir", required=True, help="Directory for output CDS FASTA files")
    parser.add_argument("--braker_base_dir", required=True, help="Base directory containing processed_* BRAKER outputs")

    args = parser.parse_args()

    convert_protein_to_cds(args.input_dir, args.output_dir, args.braker_base_dir)
