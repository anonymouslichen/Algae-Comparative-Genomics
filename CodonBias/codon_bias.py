import os
import pandas as pd
import codonbias as cb
from Bio import SeqIO

# Hardcoded paths
input_dir = "/home/stan0477/meye2099/Algae_Evolution/Trebouxiophyceae/CDS/trimmed_alignments"
output_file = "codon_bias_gc_enc_rscu_combined2.csv"

# Initialize ENC calculator
enc_calc = cb.scores.EffectiveNumberOfCodons()

# Initialize RSCU calculator
rscu_calc = cb.scores.RelativeSynonymousCodonUsage()

# Results list
results = []

# GC1/2/3 computation using codon parsing
# Includes ATG, TGG, and stop codons for GC1/GC2 but excludes them from GC3
def compute_gc123(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    gc1 = gc2 = gc3 = count1 = count2 = count3 = 0

    for codon in codons:
        if len(codon) != 3 or not all(base in "ATGC" for base in codon):
            continue

        # GC1 and GC2: include all codons
        gc1 += codon[0] in "GC"
        gc2 += codon[1] in "GC"
        count1 += 1
        count2 += 1

        # GC3: exclude Met (ATG), Trp (TGG), and stop codons
        if codon not in ["ATG", "TGG", "TAA", "TAG", "TGA"]:
            gc3 += codon[2] in "GC"
            count3 += 1

    return (
        gc1 / count1 if count1 else 0,
        gc2 / count2 if count2 else 0,
        gc3 / count3 if count3 else 0,
    )

# Loop through all relevant files
for filename in os.listdir(input_dir):
    if filename.endswith("_codon_alignment.phy"):
        file_path = os.path.join(input_dir, filename)
        gene_id = filename.split("_")[0]
        try:
            sequences = [str(record.seq).upper() for record in SeqIO.parse(file_path, "phylip")]
            enc_scores = enc_calc.get_score(sequences)
            rscu_weights = rscu_calc.get_weights(sequences)
            
            for i, seq in enumerate(sequences):
                gc1, gc2, gc3 = compute_gc123(seq)
                enc = enc_scores[i] if i < len(enc_scores) else None

                # Get codon order for this sequence
                codon_order = rscu_calc._get_codon_vector(seq)  # returns list of codons
                weights = rscu_weights[i] if i < len(rscu_weights) else []

                # Make a dictionary of codon -> RSCU value
                codon_rscu = dict(zip(codon_order, weights))

                # Prepare output dictionary for this sequence
                out_dict = {
                    "gene_id": gene_id,
                    "sequence_id": i,
                    "ENC": enc,
                    "GC1": gc1,
                    "GC2": gc2,
                    "GC3": gc3,
                }

                # Add RSCU codon values to output dict, skip invalid codons
                for codon, rscu_val in codon_rscu.items():
                    if (
                        not codon or 
                        not isinstance(codon, str) or 
                        len(codon) != 3 or 
                        not all(base in "ATGC" for base in codon)
                    ):
                        print(f"Invalid codon detected in {gene_id}, seq {i}: '{codon}' → skipping")
                        continue
                    out_dict[f"RSCU_{codon}"] = rscu_val

                results.append(out_dict)

        except Exception as e:
            print(f"Error processing {filename}: {e}")

# Save to CSV
df = pd.DataFrame(results)
df.to_csv(output_file, index=False)
print(f"ENC + GC1/2/3 + RSCU values written to {output_file}")
