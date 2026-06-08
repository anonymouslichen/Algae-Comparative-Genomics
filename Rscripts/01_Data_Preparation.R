################################################################################
# Script 1: Data Preparation and Filtering
################################################################################

# Load required libraries
library(dplyr) 
library(tidyr)
library(stringr)
library(readr)
library(ape)

# Set working directory (adjust as needed)
setwd("~/Desktop/Algae-Comparative-Genomics/")

################################################################################
# 1. LOAD AND PREPARE CODEML DATA
################################################################################

codeml <- read.csv("data/compiled_results_pooled.csv")

# Add degrees of freedom based on model
codeml <- codeml %>%
  mutate(DF = case_when(
    Model == "M0" ~ 13,
    Model == "M2" ~ 14,
    Model == "M4" ~ 23
  )) %>%
  separate(branch, into = c("Node", "Tip"), sep = "\\.\\.", remove = FALSE)

################################################################################
# 2. DEFINE SPECIES PAIRS AND ADD METADATA
################################################################################

# Function to add species names and condition
add_species_metadata <- function(df) {
  df %>%
    mutate(
      Taxa = case_when(
        Tip == 1 ~ "Asterochloris erici",
        Tip == 2 ~ "Coccomyxa subellipsoidea",
        Tip == 3 ~ "Coccomyxa viridis",
        Tip == 4 ~ "Myrmecia bisecta",
        Tip == 5 ~ "Symbiochloris irregularis",
        Tip == 6 ~ "Symbiochloris reticulata",
        Tip == 7 ~ "Trebouxia sp.C0010"
      ),
      Condition = case_when(
        Tip == 1 ~ "Lichen-forming",
        Tip == 2 ~ "Free-living",
        Tip == 3 ~ "Lichen-forming",
        Tip == 4 ~ "Free-living",
        Tip == 5 ~ "Free-living",
        Tip == 6 ~ "Lichen-forming",
        Tip == 7 ~ "Lichen-forming"
      ),
      TaxonPair = case_when(
        Tip %in% c(1, 4) ~ "Asterochloris_vs_Myrmecia",
        Tip %in% c(2, 3) ~ "Coccomyxa",
        Tip %in% c(5, 6) ~ "Symbiochloris",
        Tip %in% c(4, 7) ~ "Trebouxia_vs_Myrmecia"
      )
    )
}

################################################################################
# 3. FILTER M4 MODEL DATA FOR dN/dS ANALYSIS
################################################################################

# Filter for M4 model and relevant tips
M4_codeml <- codeml %>%
  filter(Model == "M4", Tip %in% 1:7) %>%
  add_species_metadata()

n_distinct(M4_codeml$SOG)

# paml_SOGs <- unique(M4_codeml$SOG)

# Apply filtering criteria (modified: dS > 3 instead of > 2)
failing_SOGs <- M4_codeml %>%
  filter(dS <= 0.01 | dS > 3 | omega >= 10) %>%
  pull(SOG) %>%
  unique()

length(failing_SOGs)

# Filter out failing genes from all taxa
M4_codeml_filter <- M4_codeml %>%
  filter(!(SOG %in% failing_SOGs))

n_distinct(M4_codeml_filter$SOG)

# Duplicate Myrmecia bisecta rows for both comparisons
myrmecia_rows <- M4_codeml_filter %>% 
  filter(Taxa == "Myrmecia bisecta")

myr_ast <- myrmecia_rows %>%
  mutate(TaxonPair = "Asterochloris_vs_Myrmecia")

myr_tre <- myrmecia_rows %>%
  mutate(TaxonPair = "Trebouxia_vs_Myrmecia")

# Remove original Myrmecia rows and add back duplicated ones
M4_codeml_filter <- M4_codeml_filter %>%
  filter(Taxa != "Myrmecia bisecta") %>%
  bind_rows(myr_ast, myr_tre)

# Set factor order for taxa (for plotting)
taxa_order <- c(
  "Symbiochloris irregularis", "Symbiochloris reticulata",
  "Coccomyxa subellipsoidea", "Coccomyxa viridis",
  "Myrmecia bisecta",
  "Asterochloris erici", "Trebouxia sp.C0010"
)

M4_codeml_filter$Taxa <- factor(M4_codeml_filter$Taxa, levels = taxa_order)
M4_codeml_filter$Condition <- factor(
  M4_codeml_filter$Condition,
  levels = c("Lichen-forming", "Free-living")
)
M4_codeml_filter$TaxonPair <- as.character(M4_codeml_filter$TaxonPair)

# Load phylogenetic tree
tree <- read.tree("data/consensus_tree.newick.txt")
print(tree$tip.label)

################################################################################
# 4. LOAD AND PREPARE CODON BIAS DATA
################################################################################

codon <- read.csv("data/codon_bias_gc_enc.csv")

# Rename first column
colnames(codon)[1] <- "SOG"

missing_SOGs <- setdiff(codon_SOGs, paml_SOGs)

# Add metadata
codon <- codon %>%
  mutate(
    Taxa = case_when(
      sequence_id == 0 ~ "Asterochloris erici",
      sequence_id == 1 ~ "Coccomyxa subellipsoidea",
      sequence_id == 2 ~ "Coccomyxa viridis",
      sequence_id == 3 ~ "Myrmecia bisecta",
      sequence_id == 4 ~ "Symbiochloris irregularis",
      sequence_id == 5 ~ "Symbiochloris reticulata",
      sequence_id == 6 ~ "Trebouxia sp.C0010"
    ),
    Condition = case_when(
      sequence_id == 0 ~ "Lichen-forming",
      sequence_id == 1 ~ "Free-living",
      sequence_id == 2 ~ "Lichen-forming",
      sequence_id == 3 ~ "Free-living",
      sequence_id == 4 ~ "Free-living",
      sequence_id == 5 ~ "Lichen-forming",
      sequence_id == 6 ~ "Lichen-forming"
    ),
    TaxonPair = case_when(
      sequence_id %in% c(0, 3) ~ "Asterochloris_vs_Myrmecia",
      sequence_id %in% c(1, 2) ~ "Coccomyxa",
      sequence_id %in% c(4, 5) ~ "Symbiochloris",
      sequence_id %in% c(3, 6) ~ "Trebouxia_vs_Myrmecia"
    )
  )

# GC12, ENC_prime, and delta_ENC are pre-computed in the CSV;
# no manual derivation needed.

# Set factor levels
codon$Taxa <- factor(codon$Taxa, levels = taxa_order)
codon$Condition <- factor(
  codon$Condition,
  levels = c("Lichen-forming", "Free-living")
)

# Duplicate Myrmecia for both comparisons
myrmecia_codon <- codon %>% filter(Taxa == "Myrmecia bisecta")

myr_ast_codon <- myrmecia_codon %>%
  mutate(TaxonPair = "Asterochloris_vs_Myrmecia")

myr_tre_codon <- myrmecia_codon %>%
  mutate(TaxonPair = "Trebouxia_vs_Myrmecia")

codon_myrmecia_dup <- codon %>%
  filter(Taxa != "Myrmecia bisecta") %>%
  bind_rows(myr_ast_codon, myr_tre_codon)

# Set TaxonPair factor order
codon_myrmecia_dup$TaxonPair <- factor(
  codon_myrmecia_dup$TaxonPair,
  levels = c(
    "Coccomyxa",
    "Symbiochloris",
    "Trebouxia_vs_Myrmecia",
    "Asterochloris_vs_Myrmecia"
  )
)

n_distinct(codon$SOG)

################################################################################
# 5. LOAD AND PREPARE RELAX DATA
################################################################################

# Function to load and process RELAX data
load_relax_results <- function(filepath, lichen_species) {
  df <- read.csv(filepath)
  
  # Clean SOG column
  colnames(df)[1] <- "SOG"
  df <- df %>%
    mutate(SOG = sub("\\.json$", "", SOG))
  
  # Add taxa labels
  df <- df %>%
    mutate(
      Taxa = case_when(
        Branch == "Ast_eri" ~ "Asterochloris erici",
        Branch == "Coc_sub" ~ "Coccomyxa subellipsoidea",
        Branch == "Coc_vir" ~ "Coccomyxa viridis",
        Branch == "Myr_bis" ~ "Myrmecia bisecta",
        Branch == "Sym_irr" ~ "Symbiochloris irregularis",
        Branch == "Sym_ret" ~ "Symbiochloris reticulata",
        Branch == "Tre_spC0010" ~ "Trebouxia sp.C0010"
      ),
      Condition = case_when(
        Branch == "Ast_eri" ~ "Lichen-forming",
        Branch == "Coc_sub" ~ "Free-living",
        Branch == "Coc_vir" ~ "Lichen-forming",
        Branch == "Myr_bis" ~ "Free-living",
        Branch == "Sym_irr" ~ "Free-living",
        Branch == "Sym_ret" ~ "Lichen-forming",
        Branch == "Tre_spC0010" ~ "Lichen-forming"
      )
    )
  
  # Filter out nodes
  df <- df %>%
    filter(!str_starts(Branch, "Node"))
  
  # FDR correction on focal taxon only
  focal <- df %>%
    filter(Taxa == lichen_species) %>%
    mutate(padj = p.adjust(p.value, method = "BH"))
  
  nonfocal <- df %>%
    filter(Taxa != lichen_species) %>%
    mutate(padj = NA)
  
  df <- bind_rows(focal, nonfocal) %>%
    mutate(
      Result = case_when(
        padj < 0.05 & Relaxation.Parameter..K. > 1 ~ "Intensified",
        padj < 0.05 & Relaxation.Parameter..K. < 1 ~ "Relaxed",
        TRUE ~ "Not Significant"
      ),
      TestedLichen = lichen_species
    )
  
  return(df)
}

# Load all 4 RELAX analyses
relax_coccomyxa <- load_relax_results(
  "data/RELAX_analysis_output_paired_Coc2.csv",
  "Coccomyxa viridis"
)

relax_symbiochloris <- load_relax_results(
  "data/RELAX_analysis_output_paired_Sym2.csv",
  "Symbiochloris reticulata"
)

relax_trebouxia <- load_relax_results(
  "data/RELAX_analysis_output_paired_Tre2.csv",
  "Trebouxia sp.C0010"
)

relax_asterochloris <- load_relax_results(
  "data/RELAX_analysis_output_paired_Ast2.csv",
  "Asterochloris erici"
)

# Combine all results
relax_all <- bind_rows(
  relax_coccomyxa,
  relax_symbiochloris,
  relax_trebouxia,
  relax_asterochloris
)

# Extract focal results (only tested lichen species)
relax_focal <- relax_all %>%
  filter(Taxa == TestedLichen) %>%
  mutate(
    Pair = case_when(
      TestedLichen == "Coccomyxa viridis" ~ "Coccomyxa",
      TestedLichen == "Symbiochloris reticulata" ~ "Symbiochloris",
      TestedLichen == "Trebouxia sp.C0010" ~ "Trebouxia",
      TestedLichen == "Asterochloris erici" ~ "Asterochloris"
    )
  )

# Apply natural log transformation to K parameter
relax_focal <- relax_focal %>%
  mutate(ln_K = log(Relaxation.Parameter..K.))

# Set factor levels
relax_focal$Pair <- factor(
  relax_focal$Pair,
  levels = c("Coccomyxa", "Symbiochloris", "Trebouxia", "Asterochloris")
)

relax_focal_filtered <- relax_focal %>%
  filter(is.finite(ln_K), abs(ln_K) < 10)

################################################################################
# 6. LOAD MCDONALD-KREITMAN DATA
################################################################################

MK_results <- read.csv("data/mk_results.csv")

################################################################################
# 7. LOAD GO ANNOTATIONS
################################################################################

annotations <- read_tsv("data/SOGs_annotated.tsv", col_names = FALSE)

colnames(annotations) <- c(
  "SOG", "Protein_accession", "Seq_MD5_digest",
  "Seq_len", "Analysis", "Signature_accession",
  "Signature_description", "Start", "Stop", "Score",
  "Status", "Date", "InterPro_accession",
  "InterPro_description", "GO_annotation"
)

# Convert Score to numeric
annotations$Score[annotations$Score == "-"] <- NA
annotations$Score <- as.numeric(annotations$Score)

# Filter by lowest evalue (Score)
annotations <- annotations %>%
  group_by(SOG) %>%
  filter(
    !is.na(GO_annotation),
    Score == min(Score, na.rm = TRUE)
  ) %>%
  ungroup()

n_distinct(annotations$SOG)

# Prepare gene-GO mapping for topGO
allGenes <- unique(annotations[, c("SOG", "GO_annotation")])
allGenes$GO_annotation <- gsub("\\|", ",", allGenes$GO_annotation)

gene2GO_list <- strsplit(allGenes$GO_annotation, ",")
names(gene2GO_list) <- allGenes$SOG


################################################################################
# 8. SAVE PREPARED DATA
################################################################################

save(
  M4_codeml,
  M4_codeml_filter,
  tree,
  codon,
  codon_myrmecia_dup,
  relax_all,
  relax_focal,
  relax_focal_filtered,
  MK_results,
  annotations,
  gene2GO_list,
  taxa_order,
  file = "Rscripts/prepared_data.RData"
)

