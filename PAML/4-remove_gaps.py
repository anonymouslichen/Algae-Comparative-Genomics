import os
import subprocess

def run_trimal(input_dir, output_dir, trimal_path="/home/stan0477/meye2099/bin/trimal/source/trimal", params="-gappyout -phylip_paml_m10"):
    """
    Run trimAl on all alignment files in the input directory.
    
    Parameters:
    - input_dir: Path to the directory containing alignment files.
    - output_dir: Path to the directory to save trimmed alignment files.
    - trimal_path: Path to the trimAl executable (default: 'trimal').
    - params: Additional parameters to pass to trimAl (default: '-gappyout -phylip_paml').
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith((".fa")):  # Modify extensions as needed
            input_file = os.path.join(input_dir, filename)
            # Change output file extension to .phy
            output_file = os.path.join(output_dir, filename.replace(".fa", ".phy"))
            
            # Run trimAl command
            command = [trimal_path, "-in", input_file, "-out", output_file] + params.split()
            try:
                print(f"Running: {' '.join(command)}")
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while processing {input_file}: {e}")

if __name__ == "__main__":
    # Define your input and output directories
    input_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/pal2nal_dir"
    output_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/trimmed_alignments"
    
    # Run trimAl on the alignments
    run_trimal(input_dir, output_dir)
