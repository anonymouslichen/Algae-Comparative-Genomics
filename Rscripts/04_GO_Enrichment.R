# ============================================================
# GO ENRICHMENT ANALYSIS
# Testing for functional enrichment among genes under altered selection
# Two complementary approaches:
#   1. topGO: Standard enrichment of significantly changing genes
#   2. Continuous K: Distribution differences across GO categories
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(topGO)
library(GO.db)
library(stringr)

# Load analyzed data
load("/Users/Abigail/Desktop/Claude/Rscripts:Data/analysis_complete.RData")

# ============================================================
# PREPARATION: GENE-GO MAPPING
# ============================================================

# Filter and prepare GO annotations
allGenes_annotated <- annotations %>%
  filter(!is.na(GO_annotation)) %>%
  group_by(SOG) %>%
  filter(Score == min(Score, na.rm = TRUE)) %>%  # Keep best annotation per gene
  ungroup() %>%
  dplyr::select(SOG, GO_annotation) %>%
  distinct()

# Convert pipe-separated GO terms to comma-separated
allGenes_annotated$GO_annotation <- gsub("\\|", ",", allGenes_annotated$GO_annotation)

# Create gene2GO mapping list for topGO
gene2GO_list <- strsplit(allGenes_annotated$GO_annotation, ",")
names(gene2GO_list) <- allGenes_annotated$SOG


# ============================================================
# APPROACH 1: topGO ENRICHMENT (BINARY CLASSIFICATION)
# ============================================================

# Function for lineage-specific GO enrichment
run_lineage_GO_enrichment <- function(relax_data, lineage_name, 
                                      gene2GO, ontology = "BP",
                                      min_genes = 5) {
  
  # Get genes for this lineage
  lineage_data <- relax_data %>% filter(Pair == lineage_name)
  
  # Define gene lists
  relaxed_genes <- lineage_data %>% 
    filter(Result == "Relaxed") %>% 
    pull(SOG) %>% 
    unique()
  
  strengthened_genes <- lineage_data %>% 
    filter(Result == "Strengthened") %>% 
    pull(SOG) %>% 
    unique()
  
  # All genes in this lineage
  all_lineage_genes <- unique(lineage_data$SOG)
  
  # Filter gene2GO to only genes in this lineage
  gene2GO_lineage <- gene2GO[names(gene2GO) %in% all_lineage_genes]
  
  if (length(gene2GO_lineage) < 10) {
    return(NULL)
  }
  
  results_list <- list()
  
  # Test relaxed genes
  if (length(relaxed_genes) >= min_genes) {
    
    geneList <- factor(as.integer(names(gene2GO_lineage) %in% relaxed_genes))
    names(geneList) <- names(gene2GO_lineage)
    
    tryCatch({
      GOdata <- new("topGOdata", 
                    ontology = ontology, 
                    allGenes = geneList, 
                    annot = annFUN.gene2GO, 
                    gene2GO = gene2GO_lineage,
                    nodeSize = 5)
      
      resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
      
      # Get ALL tested GO terms
      n_tested_terms <- length(usedGO(GOdata))
      
      # Extract ALL results
      results_table <- GenTable(GOdata, 
                                weight01 = resultFisher,
                                orderBy = "weight01", 
                                topNodes = n_tested_terms,  # ALL terms
                                numChar = 1000)
      
      # Convert p-values properly (handle "< 1e-05" format)
      results_table$Fisher <- sapply(results_table$weight01, function(x) {
        if (grepl("<", x)) {
          return(as.numeric(gsub("< ", "", x)))
        } else {
          return(as.numeric(x))
        }
      })
      
      # Apply FDR correction to ALL tested terms (correct practice)
      results_table$padj <- p.adjust(results_table$Fisher, method = "BH")
      results_table$Category <- "Relaxed"
      results_table$Lineage <- lineage_name
      results_table$Ontology <- ontology
      
      # Keep only top 50 for reporting (but FDR was calculated on all)
      results_table <- results_table[1:min(50, nrow(results_table)), ]
      
      results_list$relaxed <- results_table
    }, error = function(e) {
    })
  }
  
  # Test strengthened genes
  if (length(strengthened_genes) >= min_genes) {
    
    geneList <- factor(as.integer(names(gene2GO_lineage) %in% strengthened_genes))
    names(geneList) <- names(gene2GO_lineage)
    
    tryCatch({
      GOdata <- new("topGOdata", 
                    ontology = ontology, 
                    allGenes = geneList, 
                    annot = annFUN.gene2GO, 
                    gene2GO = gene2GO_lineage,
                    nodeSize = 5)
      
      resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
      
      # CORRECT: Get ALL tested GO terms, not just top 20
      n_tested_terms <- length(usedGO(GOdata))
      
      # Extract ALL results
      results_table <- GenTable(GOdata, 
                                weight01 = resultFisher,
                                orderBy = "weight01", 
                                topNodes = n_tested_terms,  # ALL terms
                                numChar = 1000)
      
      # Convert p-values properly (handle "< 1e-05" format)
      results_table$Fisher <- sapply(results_table$weight01, function(x) {
        if (grepl("<", x)) {
          return(as.numeric(gsub("< ", "", x)))
        } else {
          return(as.numeric(x))
        }
      })
      
      # Apply FDR correction to ALL tested terms (correct practice)
      results_table$padj <- p.adjust(results_table$Fisher, method = "BH")
      results_table$Category <- "Strengthened"
      results_table$Lineage <- lineage_name
      results_table$Ontology <- ontology
      
      # Keep only top 50 for reporting (but FDR was calculated on all)
      results_table <- results_table[1:min(50, nrow(results_table)), ]
      
      results_list$strengthened <- results_table
    }, error = function(e) {
    })
  }
  
  return(bind_rows(results_list))
}

# Run for all lineages and ontologies
lineages <- c("Coccomyxa", "Symbiochloris", "Trebouxia", "Asterochloris")
ontologies <- c("BP", "MF", "CC")


all_enrichment_results <- expand.grid(
  Lineage = lineages,
  Ontology = ontologies,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(results = list(
    run_lineage_GO_enrichment(
      relax_focal_filtered, 
      Lineage, 
      gene2GO_list, 
      Ontology
    )
  )) %>%
  ungroup()

# Combine all results
enrichment_combined <- all_enrichment_results %>%
  pull(results) %>%
  bind_rows()

# Filter significant results
enrichment_significant <- enrichment_combined %>%
  filter(padj < 0.05)

write.csv(enrichment_combined, "/Users/Abigail/Desktop/Claude/analysis_results/analysisGO_enrichment_topGO_all.csv", row.names = FALSE)

# ============================================================
# APPROACH 2: CONTINUOUS K ANALYSIS
# ============================================================

# Join RELAX data with GO annotations
relax_with_GO <- relax_focal_filtered %>%
  left_join(
    allGenes_annotated %>% 
      separate_rows(GO_annotation, sep = ",") %>%
      dplyr::select(SOG, GO = GO_annotation),
    by = "SOG"
  )

# Function to test K values by GO category
test_K_by_GO <- function(relax_GO_data, min_genes = 10) {
  
  GO_terms <- relax_GO_data %>%
    filter(!is.na(GO)) %>%
    count(GO) %>%
    filter(n >= min_genes) %>%
    pull(GO)
  
  results <- tibble()
  
  for (go_term in GO_terms) {
    
    # Genes with this GO term
    genes_with_GO <- relax_GO_data %>%
      filter(GO == go_term) %>%
      pull(Relaxation.Parameter..K.)
    
    # Genes without this GO term
    genes_without_GO <- relax_GO_data %>%
      filter(is.na(GO) | GO != go_term) %>%
      pull(Relaxation.Parameter..K.)
    
    # Wilcoxon test
    test_result <- wilcox.test(genes_with_GO, genes_without_GO)
    
    # Effect size (median difference)
    median_diff <- median(genes_with_GO, na.rm = TRUE) - 
      median(genes_without_GO, na.rm = TRUE)
    
    results <- bind_rows(results, tibble(
      GO_term = go_term,
      n_genes = length(genes_with_GO),
      median_K_with = median(genes_with_GO, na.rm = TRUE),
      median_K_without = median(genes_without_GO, na.rm = TRUE),
      median_diff = median_diff,
      p_value = test_result$p.value
    ))
  }
  
  results$padj <- p.adjust(results$p_value, method = "BH")
  
  return(results)
}

# Run by lineage

K_GO_results <- relax_with_GO %>%
  group_by(Pair) %>%
  group_modify(~ test_K_by_GO(.x, min_genes = 10)) %>%
  ungroup()

# Get GO term descriptions
K_GO_results$GO_description <- sapply(K_GO_results$GO_term, function(x) {
  tryCatch(Term(GOTERM[[x]]), error = function(e) NA)
})

# Get GO ontology
K_GO_results$GO_ontology <- sapply(K_GO_results$GO_term, function(x) {
  tryCatch(Ontology(GOTERM[[x]]), error = function(e) NA)
})

# Significant results
K_GO_significant <- K_GO_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)


write.csv(K_GO_results, "analysis_results/GO_continuous_K_all.csv", row.names = FALSE)
write.csv(K_GO_significant, "figures/Table_GO_continuous_K_significant.csv", row.names = FALSE)

