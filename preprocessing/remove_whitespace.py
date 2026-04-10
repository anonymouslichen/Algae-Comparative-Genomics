import sys  

def process_headers(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Remove white spaces and replace them with underscores
                modified_line = line.strip().replace(' ', '_')
                outfile.write(modified_line + '\n')
            else:
                # Keep all other lines as they are
                outfile.write(line)

if __name__ == "__main__":
    input_file = sys.argv[1]  # First argument = input file
    output_file = sys.argv[2]  # Second argument = output file
    process_headers(input_file, output_file)
