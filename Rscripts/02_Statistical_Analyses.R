################################################################################
# Script 2: Statistical Analyses
################################################################################

# Load required libraries
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(purrr)
library(ggResidpanel)

# Load prepared data
load("/Users/Abigail/Desktop/Claude/Rscripts/prepared_data.RData")

# Load phylogenetic tree
library(ape)

tree <- read.tree("/Users/Abigail/Desktop/Claude/data/consensus_tree.newick.txt")
print(tree$tip.label)


################################################################################
# ANALYSIS 1: ABSOLUTE RATES (dN and dS)
################################################################################

# Summary statistics for absolute rates
absolute_rates_summary <- M4_codeml_filter %>%
  group_by(TaxonPair, Condition, Taxa) %>%
  summarise(
    n_obs = n(),
    n_genes = n_distinct(SOG),
    mean_dN = mean(dN, na.rm = TRUE),
    sd_dN = sd(dN, na.rm = TRUE),
    median_dN = median(dN, na.rm = TRUE),
    mean_dS = mean(dS, na.rm = TRUE),
    sd_dS = sd(dS, na.rm = TRUE),
    median_dS = median(dS, na.rm = TRUE),
    mean_omega = mean(omega, na.rm = TRUE),
    .groups = "drop"
  )

print(absolute_rates_summary)

# Models for absolute dN
model_dN <- lmer(dN ~ Condition * TaxonPair + (1 | SOG), 
                 data = M4_codeml_filter)

resid_panel(model_dN)

print(anova(model_dN))

emm_dN <- emmeans(model_dN, ~ Condition | TaxonPair)
contrasts_dN <- contrast(emm_dN, method = "pairwise")
contrasts_dN <- rbind(contrasts_dN, adjust = "bonferroni")
print(contrasts_dN)

# Models for absolute dS
model_dS <- lmer(dS ~ Condition * TaxonPair + (1 | SOG), 
                 data = M4_codeml_filter)

resid_panel(model_dS)

print(anova(model_dS))

emm_dS <- emmeans(model_dS, ~ Condition | TaxonPair)
contrasts_dS <- contrast(emm_dS, method = "pairwise")
contrasts_dS <- rbind(contrasts_dS, adjust = "bonferroni")

print(contrasts_dS)

# Save results
absolute_rates_results <- bind_rows(
  as.data.frame(contrasts_dN) %>% mutate(Metric = "dN"),
  as.data.frame(contrasts_dS) %>% mutate(Metric = "dS")
)

write.csv(absolute_rates_results, "/Users/Abigail/Desktop/Claude/analysis_results/absolute_rates_contrasts.csv", 
          row.names = FALSE)

# Calculate effect sizes and percent changes
# First, get mean rates by pair and condition (not by individual taxa)
rates_by_pair_condition <- M4_codeml_filter %>%
  group_by(TaxonPair, Condition) %>%
  summarise(
    mean_dN = mean(dN, na.rm = TRUE),
    mean_dS = mean(dS, na.rm = TRUE),
    median_dN = median(dN, na.rm = TRUE),
    median_dS = median(dS, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

effect_sizes_absolute <- rates_by_pair_condition %>%
  dplyr::select(TaxonPair, Condition, mean_dN, mean_dS) %>%
  pivot_wider(
    names_from = Condition, 
    values_from = c(mean_dN, mean_dS),
    names_sep = "_"
  ) %>%
  mutate(
    delta_dN = `mean_dN_Lichen-forming` - `mean_dN_Free-living`,
    pct_change_dN = 100 * delta_dN / `mean_dN_Free-living`,
    delta_dS = `mean_dS_Lichen-forming` - `mean_dS_Free-living`,
    pct_change_dS = 100 * delta_dS / `mean_dS_Free-living`
  )

print(effect_sizes_absolute)
write.csv(effect_sizes_absolute, "/Users/Abigail/Desktop/Claude/analysis_results/absolute_rates_effect_sizes.csv", 
          row.names = FALSE)

################################################################################
# ANALYSIS 1B: PHYLOGENETIC DISTANCE vs dS CORRELATION
################################################################################

# Extract pairwise phylogenetic distances
phydist <- cophenetic.phylo(tree)
print(phydist)

# Create species mapping
species_map <- data.frame(
  Taxa = c("Asterochloris erici", "Coccomyxa subellipsoidea", 
           "Coccomyxa viridis", "Myrmecia bisecta",
           "Symbiochloris irregularis", "Symbiochloris reticulata", 
           "Trebouxia sp.C0010"),
  TreeLabel = c("Ast_eri", "Coc_sub", "Coc_vir", "Myr_bis",
                "Sym_irr", "Sym_ret", "Tre_spC0010")
)

# Define species pairs for comparison
pair_comparisons <- data.frame(
  Pair = c("Coccomyxa", "Symbiochloris", 
           "Trebouxia_vs_Myrmecia", "Asterochloris_vs_Myrmecia"),
  Lichen = c("Coccomyxa viridis", "Symbiochloris reticulata",
             "Trebouxia sp.C0010", "Asterochloris erici"),
  Free = c("Coccomyxa subellipsoidea", "Symbiochloris irregularis",
           "Myrmecia bisecta", "Myrmecia bisecta")
)

# Get phylogenetic distances for each pair
pair_comparisons <- pair_comparisons %>%
  left_join(species_map, by = c("Lichen" = "Taxa")) %>%
  dplyr::rename(Lichen_label = TreeLabel) %>%
  left_join(species_map, by = c("Free" = "Taxa")) %>%
  dplyr::rename(Free_label = TreeLabel) %>%
  rowwise() %>%
  mutate(
    PhyDist = phydist[Lichen_label, Free_label]
  ) %>%
  ungroup()

print(pair_comparisons)

# Get delta dS for each pair
delta_dS_by_pair <- absolute_rates_results %>%
  filter(Metric == "dS") %>%
  dplyr::select(Pair = TaxonPair, delta_dS = estimate)
  
# Join with phylogenetic distances
phydist_vs_dS <- pair_comparisons %>%
  dplyr::select(Pair, PhyDist) %>%
  left_join(delta_dS_by_pair)

print(phydist_vs_dS)

write.csv(phydist_vs_dS, "/Users/Abigail/Desktop/Claude/analysis_results/phydist_vs_dS.csv", row.names = FALSE)

################################################################################
# ANALYSIS 2: dN/dS COMPARISONS
################################################################################

# Linear mixed-effects model with SOG as random effect
model_omega <- lmer(omega ~ Condition * TaxonPair + (1 | SOG), 
                    data = M4_codeml_filter)

resid_panel(model_omega)

print(summary(model_omega))

# ANOVA to test main effects and interaction
anova_omega <- anova(model_omega)

print(anova_omega)

# Estimated marginal means by condition and pair
emm_omega <- emmeans(model_omega, ~ Condition | TaxonPair)
print(emm_omega)

# Pairwise contrasts within each pair
contrasts_omega <- contrast(emm_omega, method = "pairwise")
contrasts_omega <- rbind(contrasts_omega, adjust = "bonferroni")
print(contrasts_omega)

# Convert to data frame for saving
dnds_results <- as.data.frame(contrasts_omega) %>%
  mutate(
    Metric = "dN/dS",
    Direction = ifelse(estimate > 0, "Higher in lichen", "Lower in lichen"),
    Significant = ifelse(p.value < 0.001, "***",
                         ifelse(p.value < 0.01, "**",
                                ifelse(p.value < 0.05, "*", "ns")))
  )

write.csv(dnds_results, "/Users/Abigail/Desktop/Claude/analysis_results/dnds_contrasts.csv", row.names = FALSE)

################################################################################
# ANALYSIS 3: RELAX K PARAMETER ANALYSIS
################################################################################

# Summary by lineage
relax_summary <- relax_focal_filtered %>%
  group_by(Pair) %>%
  summarise(
    n_genes = n(),
    median_K = median(Relaxation.Parameter..K., na.rm = TRUE),
    mean_K = mean(Relaxation.Parameter..K., na.rm = TRUE),
    sd_K = sd(Relaxation.Parameter..K., na.rm = TRUE),
    n_strengthened = sum(Result == "Strengthened"),
    n_relaxed = sum(Result == "Relaxed"),
    n_ns = sum(Result == "Not Significant"),
    pct_strengthened = round(100 * n_strengthened / n_genes, 1),
    pct_relaxed = round(100 * n_relaxed / n_genes, 1),
    pct_ns = round(100 * n_ns / n_genes, 1),
    .groups = "drop"
  )

print(relax_summary)
write.csv(relax_summary, "/Users/Abigail/Desktop/Claude/analysis_results/relax_summary.csv", row.names = FALSE)

# Test if K differs from 1 for each lineage (Wilcoxon signed-rank test)
k_tests <- relax_focal_filtered %>%
  group_by(Pair) %>%
  summarise(
    n = n(),
    median_K = median(Relaxation.Parameter..K.),
    IQR_low = quantile(Relaxation.Parameter..K., 0.25),
    IQR_high = quantile(Relaxation.Parameter..K., 0.75),
    p_value = wilcox.test(Relaxation.Parameter..K., mu = 1)$p.value,
    .groups = "drop"
  )

print(k_tests)
write.csv(k_tests, "/Users/Abigail/Desktop/Claude/analysis_results/relax_k_tests.csv", row.names = FALSE)


################################################################################
# ANALYSIS 4: CODON USAGE BIAS COMPARISONS
################################################################################

# --- ENC' (ENC prime) ---
model_ENCprime <- lmer(ENC_prime ~ Condition * TaxonPair + (1 | SOG),
                       data = codon_myrmecia_dup)

resid_panel(model_ENCprime)

print(anova(model_ENCprime))

emm_ENCprime <- emmeans(model_ENCprime, ~ Condition | TaxonPair)
print(emm_ENCprime)

contrasts_ENCprime <- contrast(emm_ENCprime, method = "pairwise")
contrasts_ENCprime <- rbind(contrasts_ENCprime, adjust = "bonferroni")
print(contrasts_ENCprime)

# --- GC3 ---
model_GC3 <- lmer(GC3 ~ Condition * TaxonPair + (1 | SOG),
                  data = codon_myrmecia_dup)

resid_panel(model_GC3)

print(anova(model_GC3))

emm_GC3 <- emmeans(model_GC3, ~ Condition | TaxonPair)
print(emm_GC3)

contrasts_GC3 <- contrast(emm_GC3, method = "pairwise")
contrasts_GC3 <- rbind(contrasts_GC3, adjust = "bonferroni")
print(contrasts_GC3)


# --- GC12 ---
model_GC12 <- lmer(GC12 ~ Condition * TaxonPair + (1 | SOG),
                   data = codon_myrmecia_dup)

resid_panel(model_GC12)

print(anova(model_GC12))

emm_GC12 <- emmeans(model_GC12, ~ Condition | TaxonPair)
print(emm_GC12)

contrasts_GC12 <- contrast(emm_GC12, method = "pairwise")
contrasts_GC12 <- rbind(contrasts_GC12, adjust = "bonferroni")
print(contrasts_GC12)


# Save pairwise contrast results
enc_results <- as.data.frame(contrasts_ENCprime) %>% mutate(Metric = "ENC_prime")
gc3_results <- as.data.frame(contrasts_GC3) %>% mutate(Metric = "GC3")
gc12_results <- as.data.frame(contrasts_GC12) %>% mutate(Metric = "GC12")

codon_results <- bind_rows(enc_results, gc3_results, gc12_results)
write.csv(codon_results, "/Users/Abigail/Desktop/Claude/analysis_results/codon_bias_contrasts.csv",
          row.names = FALSE)

################################################################################
# ANALYSIS 6: MCDONALD-KREITMAN TEST
################################################################################

# Calculate genome-wide NI using Stoletzki & Eyre-Walker method
numerator <- sum((MK_results$Ds * MK_results$Pn) / 
                   (MK_results$Ps + MK_results$Ds), na.rm = TRUE)
denominator <- sum((MK_results$Ps * MK_results$Dn) / 
                     (MK_results$Ps + MK_results$Ds), na.rm = TRUE)
NI_TG <- numerator / denominator

# Wilcoxon test on gene-level NI values
wilcox_mk <- wilcox.test(MK_results$NI, mu = 1)
print(wilcox_mk)

# Bootstrap confidence interval for NI_TG
set.seed(123)
n_boot <- 1000
boot_nitg <- replicate(n_boot, {
  resampled <- MK_results[sample(nrow(MK_results), replace = TRUE), ]
  num <- sum((resampled$Ds * resampled$Pn) / (resampled$Ps + resampled$Ds), 
             na.rm = TRUE)
  den <- sum((resampled$Ps * resampled$Dn) / (resampled$Ps + resampled$Ds), 
             na.rm = TRUE)
  num / den
})

ci_nitg <- quantile(boot_nitg, c(0.025, 0.975))
print(ci_nitg)

# Calculate alpha (proportion of adaptive substitutions)
alpha <- 1 - NI_TG

# Save results
mk_summary <- data.frame(
  n_genes = nrow(MK_results),
  NI_TG = NI_TG,
  CI_lower = ci_nitg[1],
  CI_upper = ci_nitg[2],
  alpha = alpha,
  wilcox_V = wilcox_mk$statistic,
  wilcox_p = wilcox_mk$p.value
)

write.csv(mk_summary, "/Users/Abigail/Desktop/Claude/analysis_results/mk_summary.csv", row.names = FALSE)

################################################################################
# ANALYSIS 7: GENE CONSISTENCY ACROSS LINEAGES
################################################################################

# Identify genes present in all 4 lineages
genes_all_lineages <- relax_focal_filtered %>%
  group_by(SOG) %>%
  filter(n_distinct(Pair) == 4) %>%
  ungroup()


# Classify consistency patterns
gene_consistency <- genes_all_lineages %>%
  group_by(SOG) %>%
  summarise(
    n_analyses = n(),
    n_relaxed = sum(Result == "Relaxed"),
    n_strengthened = sum(Result == "Strengthened"),
    n_ns = sum(Result == "Not Significant"),
    pattern = case_when(
      n_relaxed == 4 ~ "Consistently Relaxed",
      n_strengthened == 4 ~ "Consistently Strengthened",
      n_ns == 4 ~ "Consistently NS",
      n_relaxed == 3 ~ "Mostly Relaxed",
      n_strengthened == 3 ~ "Mostly Strengthened",
      TRUE ~ "Mixed/Variable"
    ),
    .groups = "drop"
  )

# Summary of patterns
consistency_summary <- gene_consistency %>%
  group_by(pattern) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  arrange(desc(n_genes))

# Genes with consistent patterns
consistently_strengthened <- gene_consistency %>%
  filter(pattern == "Consistently Strengthened") %>%
  pull(SOG)

consistently_relaxed <- gene_consistency %>%
  filter(pattern == "Consistently Relaxed") %>%
  pull(SOG)

# ── GO annotation table for consistently/mostly strengthened genes ──────────
strengthened_focal <- gene_consistency %>%
  filter(pattern %in% c("Consistently Strengthened", "Mostly Strengthened"))

# Per-lineage K values (wide format)
k_wide <- relax_focal_filtered %>%
  filter(SOG %in% strengthened_focal$SOG) %>%
  dplyr::select(SOG, Pair, Result, Relaxation.Parameter..K.) %>%
  mutate(K_label = sprintf("%.3f (%s)", Relaxation.Parameter..K., Result)) %>%
  dplyr::select(SOG, Pair, K_label) %>%
  tidyr::pivot_wider(names_from = Pair, values_from = K_label)

# Best annotation per gene: prefer rows with InterPro description
go_annot <- annotations %>%
  filter(SOG %in% strengthened_focal$SOG) %>%
  group_by(SOG) %>%
  arrange(!is.na(InterPro_description), !is.na(GO_annotation)) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(SOG, InterPro_accession, InterPro_description, GO_annotation)

strengthened_table <- strengthened_focal %>%
  dplyr::select(SOG, pattern) %>%
  left_join(go_annot, by = "SOG") %>%
  left_join(k_wide, by = "SOG") %>%
  dplyr::rename(Gene = SOG, Pattern = pattern,
         InterPro_ID = InterPro_accession,
         Protein_function = InterPro_description,
         GO_term = GO_annotation)

print(strengthened_table)
write.csv(strengthened_table,
          "/Users/Abigail/Desktop/Claude/analysis_results/strengthened_genes_GO_table.csv",
          row.names = FALSE)

################################################################################
# SAVE ALL RESULTS AND UPDATED DATA
################################################################################


save(
  M4_codeml_filter,
  codon,
  codon_myrmecia_dup,
  relax_all,
  relax_focal,
  relax_focal_filtered,
  MK_results,
  annotations,
  gene2GO_list,
  taxa_order,
  tree,
  phydist,
  # Analysis results
  absolute_rates_summary,
  absolute_rates_results,
  effect_sizes_absolute,
  phydist_vs_dS,
  dnds_results,
  relax_summary,
  k_tests,
  codon_results,
  mk_summary,
  gene_consistency,
  consistency_summary,
  model_ENCprime,
  model_GC3,
  model_GC12,
  emm_ENCprime,
  emm_GC3,
  emm_GC12,
  file = "/Users/Abigail/Desktop/Claude/Rscripts/analysis_complete.RData"
)

