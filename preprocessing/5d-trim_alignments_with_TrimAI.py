#!/usr/bin/env python3

import os
import subprocess
import argparse


# Remove everything after and including '|' in FASTA headers
def clean_fasta_headers(input_path, cleaned_path):
    with open(input_path, "r") as infile, open(cleaned_path, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                clean_header = line.split("|")[0].strip()
                outfile.write(clean_header + "\n")
            else:
                outfile.write(line)


# Run trimAl on all files in a directory, optionally cleaning headers first
def run_trimal(input_dir, output_dir, trimal_path, params, file_extensions=(".fa", ".fasta", ".fas"), clean_headers=True):
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if not filename.endswith(file_extensions):
            continue

        input_file = os.path.join(input_dir, filename)
        cleaned_file = os.path.join(output_dir, f"cleaned_{filename}")
        output_file = os.path.join(output_dir, filename)

        if clean_headers:
            clean_fasta_headers(input_file, cleaned_file)
            trimal_input = cleaned_file
        else:
            trimal_input = input_file

        command = [trimal_path, "-in", trimal_input, "-out", output_file] + params.split()
        print(f"Running: {' '.join(command)}")

        try:
            subprocess.run(command, check=True)
            print(f"trimAl complete: {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] trimAl failed for {input_file}: {e}")

        # Clean up temporary file
        if clean_headers and os.path.exists(cleaned_file):
            os.remove(cleaned_file)


def main():
    parser = argparse.ArgumentParser(description="Run trimAl for multiple downstream analysis types.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTA alignments.")
    parser.add_argument("--output_dir", required=True, help="Directory to save trimmed alignments.")
    parser.add_argument("--trimal_path", default="trimal", help="Path to the trimAl executable.")
    parser.add_argument("--mode", required=True, choices=["paml", "hyphy", "mk"],
                        help="Select analysis mode: paml, hyphy, or mk.")

    args = parser.parse_args()

    # Mode-specific parameters
    mode_settings = {
        "paml": {"params": "-gappyout -phylip_paml_m10", "clean_headers": False},
        "hyphy": {"params": "-gappyout -fasta", "clean_headers": True},
        "mk": {"params": "-nogaps -fasta", "clean_headers": True}
    }

    settings = mode_settings[args.mode]

    print(f"=== Running trimAl in {args.mode.upper()} mode ===")
    print(f"Input: {args.input_dir}")
    print(f"Output: {args.output_dir}")
    print(f"trimAl path: {args.trimal_path}")
    print(f"Parameters: {settings['params']}")
    print(f"Clean headers: {settings['clean_headers']}")
    print("=============================================")

    run_trimal(args.input_dir, args.output_dir, args.trimal_path,
               settings["params"], clean_headers=settings["clean_headers"])


if __name__ == "__main__":
    main()
