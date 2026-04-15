################################################################################
# Sensitivity Analysis: Effect of dS filtering thresholds on lmer/emmeans results
# 
# Purpose: Test robustness of reported findings (dN, dS, omega contrasts) 
#          across different dS upper-bound filtering thresholds.
#          Fixed filters: dS <= 0.01 removed, omega >= 10 removed, dN > 2 removed
#          Varied filter: dS upper bound (2, 3, 5, Inf)
#
# Reference: Villanueva-Cañas et al. 2013 (original recommendation dS > 2)
# Reported analysis uses: dS <= 0.01 | dS > 3 | dN > 2 | omega >= 10
################################################################################

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(lme4)
library(emmeans)
library(ggplot2)
library(purrr)

# Set working directory (adjust as needed)
setwd("~/Desktop/Algae-Comparative-Genomics/")

################################################################################
# 1. LOAD DATA
################################################################################

# Load prepared data
load("Rscripts/prepared_data.RData")

# Duplicate Myrmecia bisecta rows for both comparisons
myrmecia_rows <- M4_codeml %>% 
  filter(Taxa == "Myrmecia bisecta")

myr_ast <- myrmecia_rows %>%
  mutate(TaxonPair = "Asterochloris_vs_Myrmecia")

myr_tre <- myrmecia_rows %>%
  mutate(TaxonPair = "Trebouxia_vs_Myrmecia")

# Remove original Myrmecia rows and add back duplicated ones
M4_codeml <- M4_codeml %>%
  filter(Taxa != "Myrmecia bisecta") %>%
  bind_rows(myr_ast, myr_tre)

################################################################################
# 2. DEFINE SENSITIVITY ANALYSIS PARAMETERS
################################################################################

# Fixed filters (applied at all thresholds)
OMEGA_UPPER  <- 10
DS_LOWER     <- 0.01

# dS upper bounds to test
ds_thresholds <- c(2, 3, 5, Inf)

# dN filter: apply or not (TRUE = filter dN > 2, FALSE = no dN filter)
dn_filters <- c(TRUE, FALSE)

# Response variables to model
response_vars <- c("dN", "dS", "omega")

################################################################################
# 3. FILTERING + MODEL FITTING FUNCTION
################################################################################

run_sensitivity <- function(m4_data, ds_upper, dn_filter, response_var) {
  
  # Create threshold labels
  ds_lab <- ifelse(is.infinite(ds_upper), 
                   "No dS upper filter", 
                   paste0("dS \u2264 ", ds_upper))
  dn_lab <- ifelse(dn_filter, "dN \u2264 2", "No dN filter")
  
  # Step 1: Identify SOGs where ANY branch fails quality filters
  # (gene-level removal — if one branch is bad, drop the whole gene)
  fail_condition <- m4_data %>%
    filter(dS <= DS_LOWER | dS > ds_upper | omega >= OMEGA_UPPER)
  
  # Optionally add dN filter
  if (dn_filter) {
    fail_dn <- m4_data %>% filter(dN > 2)
    fail_condition <- bind_rows(fail_condition, fail_dn)
  }
  
  failing_SOGs <- unique(fail_condition$SOG)
  
  # Step 2: Remove entire SOGs that have any failing branch
  filtered <- m4_data %>%
    filter(!(SOG %in% failing_SOGs))
  
  n_genes    <- n_distinct(filtered$SOG)
  n_branches <- nrow(filtered)
  
  # Build formula
  f <- as.formula(paste0(response_var, " ~ Condition * TaxonPair + (1 | SOG)"))
  
  # Fit lmer
  model_fit <- lmer(f, data = filtered)
  
  # Extract emmeans contrasts (Condition within each TaxonPair)
  emm <- emmeans(model_fit, ~ Condition | TaxonPair)
  contrasts_df <- as.data.frame(contrast(emm, method = "pairwise", adjust = "none"))
  
  # Append metadata
  contrasts_df %>%
    mutate(
      response     = response_var,
      ds_threshold = ds_upper,
      ds_label     = ds_lab,
      dn_filtered  = dn_filter,
      dn_label     = dn_lab,
      n_genes      = n_genes,
      n_branches   = n_branches,
    )
}

################################################################################
# 4. RUN ACROSS ALL THRESHOLDS × RESPONSE VARIABLES
################################################################################

sensitivity_results <- expand_grid(
  ds_upper     = ds_thresholds,
  dn_filter    = dn_filters,
  response_var = response_vars
) %>%
  pmap_dfr(function(ds_upper, dn_filter, response_var) {
    run_sensitivity(M4_codeml, ds_upper, dn_filter, response_var)
  })


################################################################################
# 5. SUMMARY TABLE
################################################################################

# Gene counts per threshold 
gene_counts <- sensitivity_results %>%
  distinct(ds_label, ds_threshold, dn_label, response, n_genes, n_branches) %>%
  arrange(ds_threshold, dn_label, response)


################################################################################
# 6. FOREST PLOT
################################################################################

# Create a combined filter label for the y-axis
sensitivity_results <- sensitivity_results %>%
  mutate(
    # Reported analysis: dS <= 3, no dN filter
    is_reported  = (ds_threshold == 3 & !dn_filtered),
    ds_label     = factor(ds_label, 
                          levels = c("No dS upper filter", 
                                     "dS \u2264 5", 
                                     "dS \u2264 3", 
                                     "dS \u2264 2")),
    filter_combo = paste0(ds_label, " | ", dn_label),
    # Significance flag
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

# Order the combined label: group by dS threshold, dN filter within
combo_levels <- sensitivity_results %>%
  distinct(ds_threshold, dn_filtered, filter_combo) %>%
  arrange(desc(ds_threshold), dn_filtered) %>%  # least strict at top
  pull(filter_combo)

sensitivity_results <- sensitivity_results %>%
  mutate(filter_combo = factor(filter_combo, levels = combo_levels))

# Nicer response labels for facets
response_labels <- c("dN" = "\u0394dN",
                     "dS" = "\u0394dS",
                     "omega" = "\u0394\u03c9 (dN/dS)")

# Nicer TaxonPair labels for facet rows
taxonpair_labels <- c("Asterochloris_vs_Myrmecia" = "Asterochloris",
                      "Trebouxia_vs_Myrmecia"     = "Trebouxia")

forest_plot <- ggplot(sensitivity_results,
                      aes(x = estimate,
                          y = filter_combo,
                          xmin = estimate - 1.96 * SE,
                          xmax = estimate + 1.96 * SE,
                          color = ds_label,
                          linetype = dn_label)) +
  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # CI bars
  geom_errorbarh(height = 0.3, linewidth = 0.7) +
  # Alternative thresholds: open circles
  geom_point(
    data  = filter(sensitivity_results, !is_reported),
    shape = 21, size = 2.5, stroke = 1.0,
    fill  = "white"
  ) +
  # Reported threshold: filled circle
  geom_point(
    data  = filter(sensitivity_results, is_reported),
    shape = 16, size = 3.2
  ) +
  # Significance stars
  geom_text(aes(x = estimate + 1.96 * SE, label = sig),
            hjust = -0.3, size = 3.5, show.legend = FALSE) +
  # Scales
  scale_color_manual(
    values = c("dS \u2264 2"            = "#E69F00", 
               "dS \u2264 3"            = "#D55E00",
               "dS \u2264 5"            = "#0072B2", 
               "No dS upper filter"     = "#999999"),
    name = "dS Filter"
  ) +
  scale_linetype_manual(
    values = c("dN \u2264 2" = "solid", "No dN filter" = "dashed"),
    name   = "dN Filter"
  ) +
  # Faceting: TaxonPair rows × response columns
  facet_grid(TaxonPair ~ response,
             scales = "free_x",
             labeller = labeller(response = response_labels,
                                 TaxonPair = taxonpair_labels)) +
  labs(
    title    = "Sensitivity Analysis: Effect of dS Filtering Threshold",
    x        = "Estimated difference (Lichen-forming \u2013 Free-living)",
    y        = NULL) +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey95"),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "bottom",
    legend.box         = "vertical",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 14),
  )

print(forest_plot)

# Save plot
ggsave("figures/FigureS3_sensitivity_forest_plot.pdf", forest_plot, 
       width = 13, height = 10, dpi = 300)
ggsave("figures/FigureS3_sensitivity_forest_plot.png", forest_plot, 
       width = 13, height = 10, dpi = 300)
