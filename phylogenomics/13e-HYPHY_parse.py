import json
import csv
import os

# Base directory containing one RELAX output subdirectory per comparison
# PROJECT_DIR = your data root; defaults to current directory
PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")

# Base directory containing one RELAX output subdirectory per comparison (produced by 13a-13d). 
# Each subdirectory holds the per-SOG RELAX .json files.
relax_base_dir = os.path.join(PROJECT_DIR, "HYPHY/output/RELAX")

# Directory where the parsed per-comparison CSVs will be written
output_dir = os.path.join(PROJECT_DIR, "HYPHY/output")


# The four lichen-forming vs. free-living comparisons. Each entry is the name of
# the RELAX output subdirectory produced by the matching 13a-13d script.
comparisons = ['paired_Ast', 'paired_Coc', 'paired_Sym', 'paired_Tre']

# Header for each output CSV
header = ['OG', 'Branch', 'K', 'Relaxation Parameter (K)', 'p-value']


def parse_comparison(directory_path):
    """Parse all RELAX .json files in a comparison directory into rows."""
    csv_data = []

    for filename in os.listdir(directory_path):
        if not filename.endswith('.json'):
            continue

        # Extract the OG number from the file name (assuming format OG#####)
        og_number = filename.split('_')[1]  # Adjust if the filename format is different

        file_path = os.path.join(directory_path, filename)
        try:
            with open(file_path, 'r', encoding='utf-8-sig') as f:
                data = json.load(f)

            # Extract the relevant information
            test_results = data.get("test results", {})
            branch_attributes = data.get("branch attributes", {}).get("0", {})  # Access the "0" key

            # Loop through the branches and extract the relevant information
            for branch_name, attributes in branch_attributes.items():
                # Per-branch descriptive K
                k = attributes.get("k (general descriptive)", "N/A")

                # Overall relaxation/intensification parameter and p-value for the fit
                relaxation_parameter = test_results.get('relaxation or intensification parameter', 'N/A')
                p_value = test_results.get('p-value', 'N/A')

                csv_data.append([
                    og_number, branch_name, k, relaxation_parameter, p_value
                ])
        except json.JSONDecodeError as e:
            print(f"Skipping malformed file: {filename} (Error: {e})")
        except Exception as e:
            print(f"Error processing file {filename}: {e}")

    return csv_data


# Loop over each comparison and write one CSV per comparison
for comparison in comparisons:
    directory_path = os.path.join(relax_base_dir, comparison)

    if not os.path.isdir(directory_path):
        print(f"Skipping {comparison}: directory not found at {directory_path}")
        continue

    csv_data = parse_comparison(directory_path)

    output_file = os.path.join(output_dir, f'RELAX_analysis_output_{comparison}.csv')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)   # Write header
        writer.writerows(csv_data)  # Write data

    print(f"{comparison}: parsed {len(csv_data)} rows -> {output_file}")
