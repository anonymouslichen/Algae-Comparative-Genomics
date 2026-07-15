import os
import pandas as pd
import codonbias as cb
from Bio import SeqIO

# Set the directories (PROJECT_DIR = your data root; defaults to current directory)
PROJECT_DIR = os.environ.get("PROJECT_DIR", ".")
input_dir = os.path.join(PROJECT_DIR, "CDS/CDS_dir") #unaligned CDS directory
output_file = "codon_bias_gc_enc.csv"

# Initialize two ENC calculators
# ENC: raw effective number of codons (Wright 1990, with Sun et al. 2013 improvements)
enc_calc = cb.scores.EffectiveNumberOfCodons(bg_correction=False)
 
# ENC': background-corrected ENC (Novembre 2002), corrects for nucleotide composition
# ENC' estimates expected ENC given GC content at each codon position
enc_prime_calc = cb.scores.EffectiveNumberOfCodons(bg_correction=True)
 
# Results list
results = []
 
 
# GC content at each codon position
def compute_gc123(seq):
    gc_vals = []
    for frame in [1, 2, 3]:
        counter = cb.stats.BaseCounter(step=3, frame=frame)
        freqs = counter.count(seq).get_table(normed=True)
        gc = freqs.get('G', 0) + freqs.get('C', 0)
        gc_vals.append(gc)
    return tuple(gc_vals)
 
 
# Loop through all relevant files
for filename in os.listdir(input_dir):
    if filename.endswith(".fa"):
        file_path = os.path.join(input_dir, filename)
        gene_id = filename.split("_")[0]
 
        try:
            sequences = [str(record.seq).upper() for record in SeqIO.parse(file_path, "fasta")]
 
            # Compute ENC and ENC' for all sequences in the file
            enc_scores = enc_calc.get_score(sequences)
            enc_prime_scores = enc_prime_calc.get_score(sequences)
 
            for i, seq in enumerate(sequences):
                gc1, gc2, gc3 = compute_gc123(seq)
 
                enc = enc_scores[i] if i < len(enc_scores) else None
                enc_prime = enc_prime_scores[i] if i < len(enc_prime_scores) else None
 
 
                results.append({
                    "gene_id": gene_id,
                    "sequence_id": i,
                    "ENC": enc,
                    "ENC_prime": enc_prime,
                    "GC1": gc1,
                    "GC2": gc2,
                    "GC3": gc3,
                    "GC12": (gc1 + gc2) / 2,
                })
 
        except Exception as e:
            print(f"Error processing {filename}: {e}")
 
# Save to CSV
df = pd.DataFrame(results)
df.to_csv(output_file, index=False)
 
print(f"Results written to {output_file}")
