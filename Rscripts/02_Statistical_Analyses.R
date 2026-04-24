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
library(ape)

# Load prepared data
load("/Users/Abigail/Desktop/Algae-Comparative-Genomics/Rscripts/prepared_data.RData")

################################################################################
# ANALYSIS 1: ABSOLUTE RATES (dN and dS)
################################################################################

# Summary statistics for rates of molecular evolution
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

### Models for absolute dN
model_dN <- lmer(dN ~ Condition * TaxonPair + (1 | SOG), 
                 data = M4_codeml_filter)

resid_panel(model_dN)
anova(model_dN)

emm_dN <- emmeans(model_dN, ~ Condition | TaxonPair)
contrasts_dN <- contrast(emm_dN, method = "pairwise", adjust = "bonferroni")
print(contrasts_dN)

# Models for absolute dS
model_dS <- lmer(dS ~ Condition * TaxonPair + (1 | SOG), 
                 data = M4_codeml_filter)

resid_panel(model_dS)
anova(model_dS)

emm_dS <- emmeans(model_dS, ~ Condition | TaxonPair)
contrasts_dS <- contrast(emm_dS, method = "pairwise", adjust = "bonferroni")
print(contrasts_dS)

# Save results
absolute_rates_results <- bind_rows(
  as.data.frame(contrasts_dN) %>% mutate(Metric = "dN"),
  as.data.frame(contrasts_dS) %>% mutate(Metric = "dS")
)

write.csv(absolute_rates_results, "analysis_results/absolute_rates_contrasts.csv", 
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

write.csv(phydist_vs_dS, "analysis_results/phydist_vs_dS.csv", row.names = FALSE)

################################################################################
# ANALYSIS 2: dN/dS COMPARISONS
################################################################################

# -----------------------------------------------------------------------------
# Model 1: Pooled main effect of Condition across taxon pairs
# TaxonPair as random intercept to account for non-independence of SOGs
# within pairs; answers "is there an overall lifestyle effect?"
# -----------------------------------------------------------------------------
model_omega_pooled <- lmer(omega ~ Condition + (1 | SOG) + (1 | TaxonPair),
                           data = M4_codeml_filter)
summary(model_omega_pooled)
resid_panel(model_omega_pooled)

# F-test on the pooled Condition effect
anova(model_omega_pooled)

# Pooled marginal means and contrast (for Supplementary Fig. 2A-style reporting)
emm_omega_pooled <- emmeans(model_omega_pooled, ~ Condition)
print(emm_omega_pooled)
contrast_omega_pooled <- contrast(emm_omega_pooled, method = "pairwise")
print(contrast_omega_pooled)

# -----------------------------------------------------------------------------
# Model 2: Condition x TaxonPair interaction
# TaxonPair as fixed effect (only 4 levels -- too few for a random slope);
# answers "does the lifestyle effect differ across pairs?" and gives per-pair
# contrasts for Table 2 / Fig. 2D
# -----------------------------------------------------------------------------
model_omega_interaction <- lmer(omega ~ Condition * TaxonPair + (1 | SOG),
                                data = M4_codeml_filter)
summary(model_omega_interaction)
resid_panel(model_omega_interaction)

# F-tests on interaction
anova(model_omega_interaction)

# Per-pair marginal means
emm_omega_bypair <- emmeans(model_omega_interaction, ~ Condition | TaxonPair)
print(emm_omega_bypair)

# Per-pair pairwise contrasts (lichen vs. free-living within each pair)
contrasts_omega_bypair <- contrast(emm_omega_bypair, method = "pairwise",
                                   adjust = "bonferroni")
print(contrasts_omega_bypair)

# Save results for Table 2
dnds_results <- as.data.frame(contrasts_omega_bypair) %>%
  mutate(
    Metric = "dN/dS",
    Direction = ifelse(estimate > 0, "Higher in lichen", "Lower in lichen"),
    Significant = ifelse(p.value < 0.001, "***",
                         ifelse(p.value < 0.01,  "**",
                                ifelse(p.value < 0.05,  "*", "ns")))
  )

write.csv(dnds_results,
          "analysis_results/dnds_contrasts.csv",
          row.names = FALSE)

################################################################################
# ANALYSIS 3: RELAX K PARAMETER ANALYSIS
################################################################################

# Summary by lineage
relax_summary <- relax_focal_filtered %>%
  group_by(Pair) %>%
  summarise(
    n_genes = n(),
    n_intensified = sum(Result == "Intensified"),
    n_relaxed = sum(Result == "Relaxed"),
    n_ns = sum(Result == "Not Significant"),
    pct_intensified = round(100 * n_intensified / n_genes, 1),
    pct_relaxed = round(100 * n_relaxed / n_genes, 1),
    pct_ns = round(100 * n_ns / n_genes, 1),
    .groups = "drop"
  )

print(relax_summary)
write.csv(relax_summary, "analysis_results/relax_summary.csv", row.names = FALSE)

# Binomial test on intensified/relaxed contingency table
p_vals <- numeric(nrow(relax_summary))

for (i in 1:nrow(relax_summary)) {
  p_vals[i] <- binom.test(
    x = relax_summary$n_intensified[i],
    n = relax_summary$n_intensified[i] + relax_summary$n_relaxed[i],
    p = 0.5
  )$p.value
}

relax_summary$binom_p <- p_vals
relax_summary$binom_p_adj <- p.adjust(p_vals, method = "bonferroni")
relax_summary

################################################################################
# ANALYSIS 4: CODON USAGE BIAS COMPARISONS
#
# Model 1: Pooled main effect of Condition across taxon pairs
# TaxonPair as random intercept to account for non-independence of SOGs
# within pairs; answers "is there an overall lifestyle effect?"
#
# Model 2: Condition x TaxonPair interaction
# TaxonPair as fixed effect  answers "does the lifestyle effect differ across 
# pairs?" Gives per pair contrasts for Table 3
#
################################################################################

# --------- ENC Prime -------------

# Model 1: Pooled main effect of Condition across taxon pairs
model_ENC_pooled <- lmer(ENC_prime ~ Condition + (1 | SOG) + (1 | TaxonPair),
                           data = codon_myrmecia_dup)
summary(model_ENC_pooled)
resid_panel(model_ENC_pooled)

# F-test on the pooled Condition effect
anova(model_ENC_pooled)

# Pooled marginal means and contrast
emm_ENC_pooled <- emmeans(model_ENC_pooled, ~ Condition)
print(emm_ENC_pooled)
contrast_ENC_pooled <- contrast(emm_ENC_pooled, method = "pairwise")
print(contrast_ENC_pooled)

# Model 2: Condition x TaxonPair interaction
model_ENC_interaction <- lmer(ENC_prime ~ Condition * TaxonPair + (1 | SOG),
                                data = codon_myrmecia_dup)
summary(model_ENC_interaction)
resid_panel(model_ENC_interaction)

# F-tests
anova(model_ENC_interaction)

# Per-pair marginal means
emm_ENC_bypair <- emmeans(model_ENC_interaction, ~ Condition | TaxonPair)
print(emm_ENC_bypair)

# Per-pair pairwise contrasts (lichen vs. free-living within each pair)
contrasts_ENC_bypair <- contrast(emm_ENC_bypair, method = "pairwise",
                                   adjust = "bonferroni")
print(contrasts_ENC_bypair)


# --------- GC3 -------------

# Model 1: Pooled main effect of Condition across taxon pairs
model_GC3_pooled <- lmer(GC3 ~ Condition + (1 | SOG) + (1 | TaxonPair),
                         data = codon_myrmecia_dup)
summary(model_GC3_pooled)
resid_panel(model_GC3_pooled)

# F-test on the pooled Condition effect
anova(model_GC3_pooled)

# Pooled marginal means and contrast
emm_GC3_pooled <- emmeans(model_GC3_pooled, ~ Condition)
print(emm_GC3_pooled)
contrast_GC3_pooled <- contrast(emm_GC3_pooled, method = "pairwise")
print(contrast_GC3_pooled)

# Model 2: Condition x TaxonPair interaction
model_GC3_interaction <- lmer(GC3 ~ Condition * TaxonPair + (1 | SOG),
                              data = codon_myrmecia_dup)
summary(model_GC3_interaction)
resid_panel(model_GC3_interaction)

# F-tests
anova(model_GC3_interaction)

# Per-pair marginal means
emm_GC3_bypair <- emmeans(model_GC3_interaction, ~ Condition | TaxonPair)
print(emm_GC3_bypair)

# Per-pair pairwise contrasts (lichen vs. free-living within each pair)
contrasts_GC3_bypair <- contrast(emm_GC3_bypair, method = "pairwise",
                                 adjust = "bonferroni")
print(contrasts_GC3_bypair)

# --------- GC12 -------------

# Model 2: Condition x TaxonPair interaction
model_GC12 <- lmer(GC12 ~ Condition * TaxonPair + (1 | SOG),
                   data = codon_myrmecia_dup)
summary(model_GC12)
resid_panel(model_GC12)

# F-tests
anova(model_GC12)

# Per-pair marginal means
emm_GC12 <- emmeans(model_GC12, ~ Condition | TaxonPair)
print(emm_GC12)

# Per-pair pairwise contrasts (lichen vs. free-living within each pair)
contrasts_GC12_bypair <- contrast(emm_GC12, method = "pairwise",
                           adjust = "bonferroni")
print(contrasts_GC12_bypair)


# Save pairwise contrast results (Table 3)
enc_results <- as.data.frame(contrasts_ENC_bypair) %>% mutate(Metric = "ENC_prime")
gc3_results <- as.data.frame(contrasts_GC3_bypair) %>% mutate(Metric = "GC3")
gc12_results <- as.data.frame(contrasts_GC12_bypair) %>% mutate(Metric = "GC12")

codon_results <- bind_rows(enc_results, gc3_results, gc12_results)
write.csv(codon_results, "analysis_results/codon_bias_contrasts.csv",
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

write.csv(mk_summary, "analysis_results/mk_summary.csv", row.names = FALSE)

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
    n_intensified = sum(Result == "Intensified"),
    n_ns = sum(Result == "Not Significant"),
    pattern = case_when(
      n_relaxed == 4 ~ "Consistently Relaxed",
      n_intensified == 4 ~ "Consistently Intensified",
      n_ns == 4 ~ "Consistently NS",
      n_relaxed == 3 ~ "Mostly Relaxed",
      n_intensified == 3 ~ "Mostly Intensified",
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
consistently_intensified <- gene_consistency %>%
  filter(pattern == "Consistently Intensified") %>%
  pull(SOG)

consistently_relaxed <- gene_consistency %>%
  filter(pattern == "Consistently Relaxed") %>%
  pull(SOG)

# ── GO annotation table for intensified genes ──────────
intensified <- relax_focal_filtered %>%
  filter(Result == "Intensified")

# Best annotation per gene: prefer rows with InterPro description
go_annot <- annotations %>%
  filter(SOG %in% intensified$SOG) %>%
  group_by(SOG) %>%
  arrange(!is.na(InterPro_description), !is.na(GO_annotation)) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(SOG, InterPro_accession, InterPro_description, GO_annotation)

intensified_table <- intensified %>%
  left_join(go_annot, by = "SOG") %>%
  dplyr::rename(Gene = SOG,
         InterPro_ID = InterPro_accession,
         Protein_function = InterPro_description,
         GO_term = GO_annotation)

print(intensified_table)
write.csv(intensified_table,
          "analysis_results/intensified_genes_GO_table.csv",
          row.names = FALSE)

# ── GO annotation table for relaxed genes ──────────
relaxed <- relax_focal_filtered %>%
  filter(Result == "Relaxed")

# Best annotation per gene: prefer rows with InterPro description
go_annot <- annotations %>%
  filter(SOG %in% relaxed$SOG) %>%
  group_by(SOG) %>%
  arrange(!is.na(InterPro_description), !is.na(GO_annotation)) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(SOG, InterPro_accession, InterPro_description, GO_annotation)

relaxed_table <- relaxed %>%
  left_join(go_annot, by = "SOG") %>%
  dplyr::rename(Gene = SOG,
                InterPro_ID = InterPro_accession,
                Protein_function = InterPro_description,
                GO_term = GO_annotation)

print(relaxed_table)
write.csv(relaxed_table,
          "analysis_results/relaxed_genes_GO_table.csv",
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
  codon_results,
  mk_summary,
  gene_consistency,
  consistency_summary,
  model_ENC_pooled,
  model_ENC_interaction,
  model_GC3_pooled,
  model_GC3_interaction,
  model_GC12,
  emm_ENC_pooled,
  emm_ENC_bypair,
  emm_GC3_pooled,
  emm_GC3_bypair,
  emm_GC12,
  file = "Rscripts/analysis_complete.RData"
)

