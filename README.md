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
- `BRAKER_DIR` and `OUTPUT_DIR` 
- `--mail-user` with your email
- Adjust the `INCLUDE_SUBDIRS` list to match your dataset

### 4. Run OrthoFinder

`scripts/4-run_orthofinder.sh`
Runs OrthoFinder on the set of protein FASTA files prepared in the previous step to infer orthogroups

**Inputs:**  
- Directory of Protein FASTA files generated from the previous step

**Outputs:**  
- OrthoFinder results directory

Before running, update:  
- `INPUT_DIR`
- `--mail-user` with your email

### 5a. Convert BRAKER protein FASTAs to CDS FASTAs

`scripts/5a-convert_braker_protein_to_cds.py`  
Converts BRAKER protein FASTA files (`*.fa`) into corresponding CDS FASTA files by locating the matching sequences in each species directory’s `braker.codingseq` file.

**Inputs:**  
- Directory of Orthogroups of BRAKER annotated protein sequences (headers must follow `species|geneID` format)  
- Directory of corresponding BRAKER outputs with CDS files

**Outputs:**  
- Direcotry of `_cds.fa` files ready for script 5b

**Before running:**  
- Install Biopython
- Confirm your BRAKER outputs are organized as processed_species/braker.codingseq.

Example command:
python scripts/5a-convert_braker_protein_to_cds.py \
    --input_dir path/to/orthogroups \
    --output_dir path/to/outputdir \
    --braker_base_dir path/to/brakerdir


### 6. Rewrite FASTA files to single-line format

`scripts/5b-rewrite_CDS_one_line.py`  
Ensures each FASTA sequence is written on a single line, improving compatibility with downstream tools that require this format.

**Inputs:**  
- Direcotry of FASTA files (`.fa`) #Output from script 5a

**Outputs:**  
- Direcotry of FASTA files with one line per sequence, ready for script 5c

Example command:
python scripts/5b-rewrite_CDS_one_line.py \
    --input_dir path/to/CDSdir \
    --output_dir path/to/outputdir
