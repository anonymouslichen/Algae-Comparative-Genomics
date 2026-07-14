#################################################################################
# GO ENRICHMENT ANALYSIS
#################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(topGO)
library(GO.db)
library(stringr)
library(here)

# Paths resolved relative to the project root via here::here()
dir.create(here("analysis_results"), showWarnings = FALSE)

# Load prepared data
load(here("Rscripts", "prepared_data.RData"))

################################################################################
# topGO ENRICHMENT FUNCTION
################################################################################

run_topGO <- function(relax_data, lineage, condition, gene2GO,
                      ontology = "BP", node_size = 5, min_sig_genes = 5) {

  lineage_data <- relax_data %>% filter(Pair == lineage)

  focal_genes <- lineage_data %>%
    filter(Result == condition) %>%
    pull(SOG) %>%
    unique()

  all_lineage_genes <- unique(lineage_data$SOG)

  # Restrict gene universe to this lineage
  gene2GO_lineage <- gene2GO[names(gene2GO) %in% all_lineage_genes]

  if (length(focal_genes) < min_sig_genes || length(gene2GO_lineage) < 10) {
    message(sprintf("Skipping %s - %s: too few genes (%d focal, %d in universe)",
                    lineage, condition, length(focal_genes), length(gene2GO_lineage)))
    return(NULL)
  }

  # Binary gene score: 1 = focal (intensified or relaxed), 0 = background
  geneList <- factor(as.integer(names(gene2GO_lineage) %in% focal_genes))
  names(geneList) <- names(gene2GO_lineage)

  tryCatch({
    GOdata <- new("topGOdata",
                  ontology   = ontology,
                  allGenes   = geneList,
                  annot      = annFUN.gene2GO,
                  gene2GO    = gene2GO_lineage,
                  nodeSize   = node_size)

    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

    n_terms <- length(usedGO(GOdata))

    results_table <- GenTable(GOdata,
                              weight01 = resultFisher,
                              orderBy  = "weight01",
                              topNodes = n_terms,
                              numChar  = 1000)

    # Parse p-values (GenTable returns strings like "< 1e-30")
    results_table$Fisher_p <- sapply(results_table$weight01, function(x) {
      if (grepl("<", x)) as.numeric(gsub("[< ]", "", x)) else as.numeric(x)
    })

    # FDR correction across all tested terms
    results_table$padj <- p.adjust(results_table$Fisher_p, method = "BH")

    results_table$Lineage   <- lineage
    results_table$Condition <- condition
    results_table$Ontology  <- ontology
    results_table$N_focal   <- length(focal_genes)
    results_table$N_universe <- length(gene2GO_lineage)

    return(results_table)

  }, error = function(e) {
    message(sprintf("Error in %s - %s - %s: %s", lineage, condition, ontology, e$message))
    return(NULL)
  })
}

################################################################################
# RUN ALL 8 ANALYSES (4 lineages x 2 conditions) across 3 ontologies each
################################################################################

lineages   <- c("Coccomyxa", "Symbiochloris", "Trebouxia", "Asterochloris")
conditions <- c("Intensified", "Relaxed")
ontologies <- c("BP", "MF", "CC")

analyses <- expand.grid(
  Lineage   = lineages,
  Condition = conditions,
  Ontology  = ontologies,
  stringsAsFactors = FALSE
)

all_results <- vector("list", nrow(analyses))

for (i in seq_len(nrow(analyses))) {
  lin  <- analyses$Lineage[i]
  cond <- analyses$Condition[i]
  ont  <- analyses$Ontology[i]
  message(sprintf("[%d/%d] %s | %s | %s", i, nrow(analyses), lin, cond, ont))

  all_results[[i]] <- run_topGO(
    relax_data = relax_focal_filtered,
    lineage    = lin,
    condition  = cond,
    gene2GO    = gene2GO_list,
    ontology   = ont
  )
}

enrichment_all <- bind_rows(all_results)

################################################################################
# SAVE OUTPUTS
################################################################################

# Full results table
write.csv(enrichment_all,
          here("analysis_results", "GO_enrichment_topGO_all.csv"),
          row.names = FALSE)

# Significant results only (padj < 0.05)
enrichment_sig <- enrichment_all %>% filter(padj < 0.05)

write.csv(enrichment_sig,
          here("analysis_results", "GO_enrichment_topGO_significant.csv"),
          row.names = FALSE)

