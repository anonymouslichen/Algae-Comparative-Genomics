#!/usr/bin/env python3

import os
import argparse

# Rewrite one FASTA file so each sequence appears on a single line
def rewrite_fasta_single_line(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    header = None
    sequence = []
    output_lines = []

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                output_lines.append(header + '\n' + ''.join(sequence) + '\n')
            header = line
            sequence = []
        else:
            sequence.append(line)

    if header is not None:
        output_lines.append(header + '\n' + ''.join(sequence) + '\n')

    with open(output_file, 'w') as outfile:
        outfile.writelines(output_lines)


# Process all .fa files in a directory
def process_fasta_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_files = [f for f in os.listdir(input_dir) if f.endswith(".fa")]
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}")
        return

    for filename in fasta_files:
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        rewrite_fasta_single_line(input_path, output_path)
        print(f"Processed: {filename}")


def main():
    parser = argparse.ArgumentParser(description="Reformat FASTA files to single-line sequences.")
    parser.add_argument("--input_dir", required=True, help="Directory containing FASTA files to reformat")
    parser.add_argument("--output_dir", required=True, help="Directory to write reformatted FASTA files")

    args = parser.parse_args()
    process_fasta_files(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
