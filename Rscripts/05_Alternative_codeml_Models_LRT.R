################################################################################
# Script 1: LRT with different codeml models and sensitivity analysis
################################################################################

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Set working directory (adjust as needed)
setwd("~/Desktop/Claude")

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

# Filter data frame to include only one line for each SOG model
model_lnL <- codeml %>%
  distinct(SOG, Model, .keep_all = TRUE)

# Perform Likelihood Ratio Test
perform_lr_test <- function(null_model, alt_model) {
  lr_statistic <- -2 * (null_model$lnL - alt_model$lnL)
  df_diff <- abs(null_model$DF - alt_model$DF)
  p_value <- pchisq(lr_statistic, df = df_diff, lower.tail = FALSE)
  
  return(tibble(
    gene = null_model$SOG,
    null_model = null_model$Model,
    alt_model = alt_model$Model,
    LR_statistic = lr_statistic,
    df_diff = df_diff,
    p_value = p_value
  ))
}

# Process all genes efficiently
results <- model_lnL %>%
  group_by(SOG) %>%
  group_split() %>%
  purrr::map_dfr(function(SOG_data) {
    null_model <- filter(SOG_data, Model == "M0")
    alt_models <- filter(SOG_data, Model != "M0")
    
    if (nrow(null_model) == 0 | nrow(alt_models) == 0) {
      return(tibble())  # Skip if no null or alternative model
    }
    
    purrr::map_dfr(1:nrow(alt_models), function(i) {
      perform_lr_test(null_model, alt_models[i, ])
    })
  })

# Create a function to calculate the percentage of genes with p-value < 0.01
calculate_percentage <- function(data) {
  # Filter genes with p-value < 0.01
  significant_genes <- subset(data, p_value < 0.01)
  
  # Calculate percentage
  percentage <- (nrow(significant_genes) / nrow(data)) * 100
  
  return(percentage)
}

# Calculate percentage for each model comparison
percentage_M0_M2 <- calculate_percentage(subset(results, alt_model == "M2"))
percentage_M0_M4 <- calculate_percentage(subset(results, alt_model == "M4"))


