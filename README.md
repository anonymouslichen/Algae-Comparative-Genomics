# Algae-Comparative-Genomics

### 1. Soft-masking workflow

`scripts/1-soft_masking.sh`  
Runs RepeatMasker on each genome in a directory

**Inputs:**  
- Directory of `.fna` genome files

**Outputs:**  
- Soft-masked genome files in the specified output directory

**Helper script:**  
`scripts/remove_whitespace.py` — replaces spaces in FASTA headers with underscores before masking.


### 2. Genome annotation with BRAKER3

`scripts/2-braker_annotation.sh`  
Runs BRAKER3 within a Singularity container to annotate softmasked genomes in parallel using SLURM array jobs.

**Inputs:**  
- `processed_*.fna` files from Step 1  
- `Viridiplantae.fa` (reference protein dataset)  

**Outputs:**  
- BRAKER output folders under `braker_output/`

Before running, update:
- `GENOMES_DIR`, `OUTPUT_DIR`, and `AUGUSTUS_CONFIG_PATH`
- `--array=0-##` to match the number of genomes
- `--mail-user` with your email

### 3. Prepare BRAKER protein files for OrthoFinder

`scripts/3-prepare_braker_for_orthofinder.sh`
Collects all braker.aa protein files from selected BRAKER output directories, renames them by species, and edits FASTA headers to include the species name as a prefix (e.g., >species|gene_id). The standardized files are placed together for OrthoFinder input.

**Inputs:** 
- Per-genome BRAKER output directories (e.g., processed_*, each containing a braker.aa file)

**Outputs:** 
- Renamed and header-edited FASTA files in the specified OUTPUT_DIR, ready for OrthoFinder analysis

Before running, update:
-`BRAKER_DIR` and `OUTPUT_DIR` 
- `--mail-user` with your email
- Adjust the `INCLUDE_SUBDIRS` list to match your dataset
