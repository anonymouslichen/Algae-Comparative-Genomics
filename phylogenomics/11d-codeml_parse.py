import os
import csv
from Bio.Phylo.PAML import codeml
import re


# PROJECT_DIR = your data root; defaults to current directory
PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")
  
# Directory containing the output files
directory_path = os.path.join(PROJECT_DIR, "PAML/output")
  
# Output CSV file
output_csv_path = os.path.join(PROJECT_DIR, "PAML/output/compiled_results.csv")

# Function to extract specific information from a results dictionary
def extract_info(results_dict):
    extracted_info = {
        'lnL': None,
        'branches': []
    }

    # Extract lnL from NSsites
    if 'NSsites' in results_dict:
        for site in results_dict['NSsites'].values():
            if 'lnL' in site:
                extracted_info['lnL'] = site['lnL']

            # Extract dN and dS for each branch
            if 'parameters' in site and 'branches' in site['parameters']:
                for branch, values in site['parameters']['branches'].items():
                    branch_info = {
                        'branch': branch,
                        'dN': values.get('dN', None),
                        'dS': values.get('dS', None),
                        'omega': values.get('omega', None)
                    }
                    extracted_info['branches'].append(branch_info)

    return extracted_info

# Initialize a list to store the extracted information
compiled_data = []

# Regex to extract SOG and Model from filename
filename_pattern = re.compile(r'Trebouxiophyceae_(OG\d+)_([^_]+)\.out')

# Loop through all files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith('.out'):  # Assuming the output files have .out extension
        match = filename_pattern.match(filename)
        if match:
            SOG = match.group(1)
            Model = match.group(2)

            file_path = os.path.join(directory_path, filename)
            try:
                results = codeml.read(file_path)
                extracted_info = extract_info(results)

                # Flatten the extracted_info for each branch to add to compiled_data
                lnL = extracted_info['lnL']
                for branch_info in extracted_info['branches']:
                    branch_info['lnL'] = lnL
                    #branch_info['file'] = filename  # Add the source file information
                    branch_info['SOG'] = SOG  # Add SOG information
                    branch_info['Model'] = Model  # Add Model information
                    compiled_data.append(branch_info)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

# Define the keys (column names) for the CSV file
if compiled_data:
    keys = compiled_data[0].keys()
else:
    keys = []

# Write the compiled data to a CSV file
with open(output_csv_path, 'w', newline='') as output_csv:
    dict_writer = csv.DictWriter(output_csv, fieldnames=keys)
    dict_writer.writeheader()
    dict_writer.writerows(compiled_data)

print(f"Compiled results have been written to {output_csv_path}")

