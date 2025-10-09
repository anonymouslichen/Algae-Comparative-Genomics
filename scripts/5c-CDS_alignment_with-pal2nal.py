#!/usr/bin/env python3

import os
import argparse
import subprocess

# Run MUSCLE to align protein sequences
def run_muscle(input_file, output_file, muscle_path="muscle"):
    try:
        subprocess.run([muscle_path, "-in", input_file, "-out", output_file], check=True)
        print(f"MUSCLE alignment created: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] MUSCLE failed on {input_file}: {e}")
        return None


# Reorder sequences in a FASTA alignment alphabetically by header
def reorder_sequences_alphabetically(alignment_file):
    with open(alignment_file, 'r') as f:
        lines = f.readlines()

    sequences = []
    current_header, current_sequence = None, []

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                sequences.append((current_header, ''.join(current_sequence)))
            current_header = line
            current_sequence = []
        else:
            current_sequence.append(line)

    if current_header:
        sequences.append((current_header, ''.join(current_sequence)))

    sequences.sort(key=lambda x: x[0])

    with open(alignment_file, 'w') as f:
        for header, sequence in sequences:
            f.write(f"{header}\n{sequence}\n")

    print(f"Sequences reordered alphabetically in {alignment_file}")


# Run PAL2NAL to generate codon alignments
def run_pal2nal(protein_file, nucleotide_file, output_file, pal2nal_path):
    try:
        with open(output_file, 'w') as out_f:
            subprocess.run(
                [pal2nal_path, protein_file, nucleotide_file, "-output", "fasta"],
                check=True, stdout=out_f
            )
        print(f"PAL2NAL codon alignment created: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] PAL2NAL failed on {protein_file}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Align protein sequences and generate codon alignments.")
    parser.add_argument("--protein_dir", required=True, help="Directory containing protein FASTA files")
    parser.add_argument("--cds_dir", required=True, help="Directory containing CDS FASTA files")
    parser.add_argument("--muscle_alignment_dir", required=True, help="Directory to store MUSCLE alignments")
    parser.add_argument("--pal2nal_output_dir", required=True, help="Directory to store PAL2NAL codon alignments")
    parser.add_argument("--muscle_path", default="muscle", help="Path to MUSCLE executable (default: muscle in PATH)")
    parser.add_argument("--pal2nal_path", required=True, help="Full path to PAL2NAL Perl script")

    args = parser.parse_args()

    # Ensure output directories exist
    os.makedirs(args.muscle_alignment_dir, exist_ok=True)
    os.makedirs(args.pal2nal_output_dir, exist_ok=True)

    # Process each protein file
    protein_files = [f for f in os.listdir(args.protein_dir) if f.endswith((".fa", ".fasta"))]

    if not protein_files:
        print(f"No protein FASTA files found in {args.protein_dir}")
        return

    for protein_file in protein_files:
        protein_path = os.path.join(args.protein_dir, protein_file)
        aligned_protein_path = os.path.join(
            args.muscle_alignment_dir,
            f"{os.path.splitext(protein_file)[0]}_aligned.fa"
        )

        # Step 1: Run MUSCLE
        aligned_file = run_muscle(protein_path, aligned_protein_path, args.muscle_path)
        if not aligned_file:
            continue

        # Step 2: Reorder alphabetically
        reorder_sequences_alphabetically(aligned_file)

        # Step 3: Find matching CDS file
        cds_file = f"{os.path.splitext(protein_file)[0]}_cds.fa"
        cds_path = os.path.join(args.cds_dir, cds_file)

        if not os.path.isfile(cds_path):
            print(f"[WARNING] No CDS file found for {protein_file}, skipping.")
            continue

        # Step 4: Run PAL2NAL
        output_codon_path = os.path.join(
            args.pal2nal_output_dir,
            f"{os.path.splitext(protein_file)[0]}_codon_alignment.fa"
        )
        run_pal2nal(aligned_file, cds_path, output_codon_path, args.pal2nal_path)


if __name__ == "__main__":
    main()
