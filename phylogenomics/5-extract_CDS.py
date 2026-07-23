import os
from Bio import SeqIO

# Set the directories (PROJECT_DIR = your data root; defaults to current directory)
PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")
input_directory = os.path.join(PROJECT_DIR, "CDS/Protein_dir")
output_directory = os.path.join(PROJECT_DIR, "CDS/CDS_dir")
braker_base_directory = os.path.join(PROJECT_DIR, "Braker_min_length")

def load_cds_index(codingseq_file):
    """Index all records in a braker.codingseq file by their ID."""
    return SeqIO.to_dict(SeqIO.parse(codingseq_file, "fasta"))

def convert_protein_to_cds(input_dir, output_dir, braker_base_dir):
    os.makedirs(output_dir, exist_ok=True)
    cds_cache = {}  # cache indexes per species to avoid re-parsing

    for protein_file in os.listdir(input_dir):
        if not protein_file.endswith(".fa"):
            continue

        protein_path = os.path.join(input_dir, protein_file)
        output_file = os.path.join(
            output_dir, f"{os.path.splitext(protein_file)[0]}_cds.fa"
        )

        with open(protein_path) as prot_fh, open(output_file, "w") as out_fh:
            for record in SeqIO.parse(prot_fh, "fasta"):
                header = record.description  # full header, e.g. "species|gene"
                species, gene_info = header.split("|")
                gene_info = gene_info.strip()

                if species not in cds_cache:
                    codingseq_file = os.path.join(
                        braker_base_dir, f"processed_{species}", "braker.codingseq"
                    )
                    if not os.path.exists(codingseq_file):
                        print(f"Coding seq file not found: {codingseq_file}")
                        cds_cache[species] = None
                    else:
                        cds_cache[species] = load_cds_index(codingseq_file)

                index = cds_cache[species]
                if index and gene_info in index:
                    seq = str(index[gene_info].seq)
                    out_fh.write(f">{header}\n{seq}\n")
                else:
                    print(f"No CDS found for header: {header}")

convert_protein_to_cds(input_directory, output_directory, braker_base_directory)
