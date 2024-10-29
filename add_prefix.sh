#!/bin/bash

# Set the main directory containing the subdirectories
MAIN_DIR=/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/Braker_file_edits  # Replace with the path to your main directory

# Loop through each subdirectory in the main directory
for dir in $MAIN_DIR/processed_*; do
  # Extract the prefix (subdirectory name without "processed_")
  prefix=$(basename $dir | sed 's/^processed_//')

  # Loop through the target files in the current subdirectory
  for file in $dir/braker.aa $dir/braker.codingseq; do
    # Use sed to add the prefix to each sequence header line (starting with ">")
    sed -i "s/^>/>${prefix}_/" $file
  done
done

echo "Prefix added to file names and sequence headers in each braker file."
