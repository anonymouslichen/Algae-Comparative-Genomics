# Lichen symbiosis does not impose a uniform genomic syndrome on algae

Comparative-genomics pipeline and analysis code for the manuscript:

> **Meyer, A. R.** *Lichen symbiosis does not impose a uniform genomic syndrome on algae.*
> Department of Ecology, Evolution and Behavior, University of Minnesota.

## Overview

This study tests whether lichen symbiosis reduces the efficacy of selection in
lichen-forming algae, as predicted if symbiosis lowers effective population size
(*N*e). Four lichen-forming algae are compared to their closest free-living
relatives across 951 single-copy orthologous genes (SOGs) using four
complementary signatures of molecular evolution:

- **dN/dS (ω)** — strength of purifying selection (PAML `codeml`)
- **K** — relaxation/intensification of selection (HyPhy `RELAX`)
- **ENC′** — codon-usage bias, corrected for GC content
- **GC3 / GC12** — base composition at synonymous vs. non-synonymous sites

The four taxon-pair comparisons are:

| Lichen-forming | Free-living | Comparison |
|---|---|---|
| *Coccomyxa viridis* | *Coccomyxa subellipsoidea* | within-genus |
| *Symbiochloris reticulata* | *Symbiochloris irregularis* | within-genus |
| *Trebouxia* sp. C0010 | *Myrmecia bisecta* | cross-genus |
| *Asterochloris erici* | *Myrmecia bisecta* | cross-genus |

## Repository structure

```
.
├── phylogenomics/   Genome-to-alignment pipeline + selection/codon analyses (shell + Python)
├── Rscripts/        Statistical models and figure generation (R)
├── data/            Intermediate result tables consumed by the R scripts
```

## Software & dependencies

Bioinformatics tools were run on a SLURM HPC cluster (`module load` / Singularity /
conda). Versions are as reported in the manuscript Methods.

| Tool | Version | Use |
|---|---|---|
| RepeatMasker | 4.1.1 (Dfam 3.3) | Repeat soft-masking |
| BRAKER3 (GeneMark-EX, AUGUSTUS) | 3.0.8 | Structural gene annotation |
| BUSCO | 5.4.3 (`chlorophyta_odb10`) | Annotation completeness |
| OrthoFinder | 2.5.4 | Single-copy ortholog inference |
| MUSCLE | 3.8.31 | Protein alignment |
| PAL2NAL | 14 | Codon (back-translated) alignment |
| TrimAl | 1.3 | Alignment trimming |
| IQ-TREE | 2.1.2 | Gene trees (LG+F+G4) |
| DendroPy / SumTrees | 4.0.0 | Consensus species tree |
| PAML `codeml` | 4.6 | dN/dS branch models |
| `codonbias` (Python) | 0.3.5 | ENC′, GC content |
| HyPhy | 2.5.71 | Tree labeling |
| RELAX | 4.5 | Selection intensity (K) |
| InterProScan | 5.23-62.0 | Functional / GO annotation |
| R | 4.5.0 | Statistics and figures |

**R packages:** `here`, `lme4`, `lmerTest`, `emmeans`, `topGO`, `GO.db`, `ape`,
`ggtree`, `ggplot2`, `patchwork`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`,
`ggrepel`, `cowplot`, `ggResidpanel`, `broom`.

## Phylogenomics pipeline (`phylogenomics/`)

Scripts are numbered in execution order. Steps 1–9 build the SOG alignment/tree
dataset; steps 10–15 are the parallel downstream analyses that consume it.

| # | Script | Tool | Purpose |
|---|---|---|---|
| 1 | `1-softmask_repeats.sh` | RepeatMasker | Soft-mask repeats in each genome |
| 2 | `2-annotation.sh` | BRAKER3 | Structural gene annotation |
| 3 | `3-BUSCO.sh` | BUSCO | Assess annotation completeness |
| 4 | `4-Orthofinder.sh` | OrthoFinder | Identify single-copy orthologs (SOGs) |
| 5 | `5-extract_CDS.py` | Biopython | Extract CDS for each SOG protein |
| 6 | `6-align_CDS.py` | MUSCLE + PAL2NAL | Protein align → codon alignment |
| 7a | `7a-remove_gaps_PAML.py` | TrimAl | Gap-trim alignments (PHYLIP for PAML) |
| 7b | `7b-remove_gaps_HYPHY.py` | TrimAl | Gap-trim alignments (FASTA for HyPhy/RELAX) |
| 8 | `8-gene_trees.sh` | IQ-TREE | Per-gene ML trees |
| 9 | `9-consensus_tree.sh` | SumTrees | Majority-rule consensus species tree |
| 10 | `10-codeml_label_trees.py` | Python | Label foreground (lichen) branches on gene trees for the codeml models |
| 11a | `11a-codeml_null.sh` | PAML | dN/dS null (one-ratio, M0) model |
| 11b | `11b-codeml_two-ratio.sh` | PAML | Two-ratio (foreground/background, M2) model |
| 11c | `11c-codeml_free-ratio.sh` | PAML | Free-ratio (per-branch ω, M4) model |
| 11d | `11d-codeml_parse.py` | Biopython | Parse codeml output → CSV |
| 12 | `12-ENC-GC123.py` | `codonbias` | ENC/ENC′ and GC1/GC2/GC3 per gene |
| 13a–d | `13{a,b,c,d}-HYPHY_label_trees_*.sh` | HyPhy | Label test/reference branches per comparison |
| 13e | `13e-HYPHY_parse.py` | Python | Parse RELAX JSON → CSV (all four comparisons; run after 14a–d) |
| 14a–d | `14{a,b,c,d}-HYPHY_RELAX_*.sh` | RELAX | Selection-intensity test per comparison |
| 15 | `15-InterProScan.sh` | InterProScan | Functional + GO annotation of SOGs |

The single gene under "mostly intensified" selection (K > 1 in three of four
lineages; `OG0003934`) was additionally annotated with **HHpred** (Zimmermann et
al. 2018) on the MPI Bioinformatics Toolkit (https://toolkit.tuebingen.mpg.de).
The translated SOG sequence (`data/hhpred_query.fasta`) was searched against the
PDB70, Pfam, SMART, and SCOPe70 profile-HMM databases with default parameters
(accessed 2026-07-08).

**Supplementary Table 1** (BUSCO completeness per genome) is taken directly from
the BUSCO output of step 3; it is not generated by a separate script.

## Statistical analyses & figures (`Rscripts/`)

Run in numeric order. `01` and `02` build the `.RData` objects that the figure
scripts load.

| Script | Produces |
|---|---|
| `01_Data_Preparation.R` | Cleans/filters inputs → `prepared_data.RData` |
| `02_Statistical_Analyses.R` | Linear mixed models `response ~ Condition * TaxonPair + (1\|SOG)`, emmeans contrasts (Bonferroni), Wilcoxon and binomial tests → `analysis_complete.RData`; **Tables 2 & 3** |
| `03a_figure2.R` | **Figure 2** — species tree + ω, dN, dS by taxon pair |
| `03b_figure3.R` | **Figure 3** — ENC′ and GC12-vs-GC3 neutrality plot |
| `03c_figure4.R` | **Figure 4** — RELAX volcano plots (ln K vs adj. p) |
| `03d_supplementary_figure_genome_sensitivity.R` | **Figure S1** — genome-choice sensitivity |
| `03e_supplementary_figure_filtering_sensitivity.R` | **Figure S2** — dS filtering sensitivity |
| `03f_supplementary_figure_phydist_dS.R` | **Figure S3** — phylogenetic distance vs ΔdS |
| `03g_supplementary_figure_omega.R` | **Figure S4** — pooled ω, free-ratio vs two-ratio |
| `04_GO_Enrichment.R` | GO enrichment (topGO `weight01` + Fisher, BH-corrected) |
| `05_Alternative_codeml_Models_LRT.R` | Likelihood-ratio model selection (M0 vs M2 vs M4) |

## Data (`data/`)

Intermediate tables output by the pipeline and consumed by the R scripts.

| File | Description |
|---|---|
| `compiled_results_pooled.csv` | Parsed codeml branch dN, dS, ω, lnL for all SOGs (main analysis) |
| `compiled_results_*.csv` | Per-genome codeml results for the *Asterochloris*/*Trebouxia* genome-choice sensitivity analysis |
| `RELAX_analysis_output_paired_{Ast,Coc,Sym,Tre}.csv` | Parsed RELAX K and p-values per comparison |
| `codon_bias_gc_enc.csv` | ENC/ENC′ and GC1/GC2/GC3 per gene |
| `SOGs_annotated.tsv` | InterProScan functional / GO annotations per SOG |
| `consensus_tree.txt` | Majority-rule consensus species tree (Newick) |
| `gene_pairwise_distances.csv` | Pairwise phylogenetic distances between taxa |
| `phydist_vs_dS.csv` | Phylogenetic distance vs ΔdS (supplementary figure input) |

## Reproducibility notes

This repository documents the full analytical pipeline for transparency. It is
**not** a turnkey clone-and-run workflow: the shell/Python pipeline scripts read
and write under a `PROJECT_DIR` root you set yourself
(`export PROJECT_DIR=/path/to/data`), but still contain institution-specific
SLURM directives and `module load` / conda environment names, and they depend on
genome assemblies, reference databases, and HPC modules not bundled here. The R
scripts, by contrast,
run against the provided `data/` and `.RData` files and regenerate the tables and
figures. They locate all paths relative to the repository root with the `here`
package (no `setwd()` and no absolute paths), so after cloning you can run them
from RStudio or via `Rscript` from anywhere — no path editing required. Output
directories (`figures/`, `analysis_results/`) are created automatically if absent.
The intermediate `.RData` objects and the `analysis_results/` CSVs are tracked so
the figure scripts (`03*`) can be run standalone, without re-executing `01`/`02`.
Raw genome assemblies are available from the public accessions below.

## Data availability

| Species | Lifestyle | Source | Accession |
|---|---|---|---|
| *Coccomyxa viridis* | Lichen-forming | NCBI | GCA_964019345.2 |
| *Coccomyxa subellipsoidea* | Free-living | NCBI | GCF_000258705.1 |
| *Symbiochloris reticulata* | Lichen-forming | JGI | 1016105 |
| *Symbiochloris irregularis* | Free-living | NCBI | GCA_040144405.1 |
| *Trebouxia* sp. C0010 | Lichen-forming | NCBI | GCA_045269315.1 |
| *Asterochloris erici* | Lichen-forming | NCBI | GCA_019693375.1 |
| *Myrmecia bisecta* | Free-living | NCBI | GCA_040144395.1 |

## Contact

Abigail R. Meyer — Department of Ecology, Evolution and Behavior, University of
Minnesota.

---

*This README was drafted by Claude (Anthropic, Opus 4.8) from the contents of
this repository, then reviewed, edited, and verified for accuracy by a human.*
