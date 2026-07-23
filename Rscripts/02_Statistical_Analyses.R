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
library(here)

# Ensure the output directory exists before writing results
dir.create(here("analysis_results"), showWarnings = FALSE)

# Load prepared data
load(here("Rscripts", "prepared_data.RData"))

################################################################################
# ANALYSIS 1: dN/dS COMPARISONS
################################################################################

model_omega_interaction <- lmer(omega ~ Condition * TaxonPair + (1 | SOG),
                                data = M4_codeml_filter)
summary(model_omega_interaction)
resid_panel(model_omega_interaction)

# F-tests on main effects and interaction
anova(model_omega_interaction)

# ---- Per-pair effects (Table 2 / Fig. 2D) ----
emm_omega_bypair <- emmeans(model_omega_interaction, ~ Condition | TaxonPair)
print(emm_omega_bypair)

contrasts_omega_bypair <- contrast(emm_omega_bypair, method = "pairwise",
                                   adjust = "bonferroni")
print(contrasts_omega_bypair)

# ---- Pooled effect across pairs (Supplementary Fig. 2A) ----
# Marginal means averaging over TaxonPair levels with equal weights
emm_omega_pooled <- emmeans(model_omega_interaction, ~ Condition)
print(emm_omega_pooled)
contrast_omega_pooled <- contrast(emm_omega_pooled, method = "pairwise")
print(contrast_omega_pooled)

# ---- Wilcoxon signed-rank test: standard paired approach (vs. mixed model) ----
# Non-parametric paired comparison of lichen vs. free-living omega within each
# TaxonPair, matched by gene (SOG). Complements the lmer interaction model above.

wilcox_omega_bypair <- M4_codeml_filter %>%
  dplyr::select(SOG, TaxonPair, Condition, omega) %>%
  # average any duplicate SOG x Condition rows so pivot gives one value per cell
  group_by(TaxonPair, SOG, Condition) %>%
  pivot_wider(names_from = Condition, values_from = omega) %>%
  drop_na() %>%                       # keep only genes present in BOTH conditions
  group_by(TaxonPair) %>%
  summarise(
    n_pairs = n(),
    median_Lichen = median(`Lichen-forming`),
    median_Free   = median(`Free-living`),
    median_diff   = median(`Lichen-forming` - `Free-living`),
    test = list(wilcox.test(`Lichen-forming`, `Free-living`, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    V_statistic = map_dbl(test, "statistic"),
    p_value     = map_dbl(test, "p.value"),
    p_adj       = p.adjust(p_value, method = "bonferroni")
  ) %>%
  dplyr::select(-test)

print(wilcox_omega_bypair)

write.csv(wilcox_omega_bypair,
          here("analysis_results", "wilcoxon_omega_bypair.csv"),
          row.names = FALSE)


################################################################################
# ANALYSIS 2: ABSOLUTE RATES (dN and dS)
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
dnds_results <- bind_rows(
  as.data.frame(contrasts_dN) %>% mutate(Metric = "dN"),
  as.data.frame(contrasts_dS) %>% mutate(Metric = "dS"),
  as.data.frame(contrasts_omega_bypair) %>% mutate(Metric = "omega")
)

write.csv(dnds_results,
          here("analysis_results", "dnds_contrasts.csv"),
          row.names = FALSE)

# ---- Wilcoxon signed-rank test: standard paired approach (vs. mixed model) dS ----.

wilcox_dS_bypair <- M4_codeml_filter %>%
  dplyr::select(SOG, TaxonPair, Condition, dS) %>%
  # average any duplicate SOG x Condition rows so pivot gives one value per cell
  group_by(TaxonPair, SOG, Condition) %>%
  pivot_wider(names_from = Condition, values_from = dS) %>%
  drop_na() %>%                       # keep only genes present in BOTH conditions
  group_by(TaxonPair) %>%
  summarise(
    n_pairs = n(),
    median_Lichen = median(`Lichen-forming`),
    median_Free   = median(`Free-living`),
    median_diff   = median(`Lichen-forming` - `Free-living`),
    test = list(wilcox.test(`Lichen-forming`, `Free-living`, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    V_statistic = map_dbl(test, "statistic"),
    p_value     = map_dbl(test, "p.value"),
    p_adj       = p.adjust(p_value, method = "bonferroni")
  ) %>%
  dplyr::select(-test)

print(wilcox_dS_bypair)

write.csv(wilcox_dS_bypair,
          here("analysis_results", "wilcoxon_dS_bypair.csv"),
          row.names = FALSE)



################################################################################
# ANALYSIS 2B: PHYLOGENETIC DISTANCE vs dS CORRELATION
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
delta_dS_by_pair <- dnds_results %>%
  filter(Metric == "dS") %>%
  dplyr::select(Pair = TaxonPair, delta_dS = estimate)

# Join with phylogenetic distances
phydist_vs_dS <- pair_comparisons %>%
  dplyr::select(Pair, PhyDist) %>%
  left_join(delta_dS_by_pair)

print(phydist_vs_dS)

write.csv(phydist_vs_dS, here("analysis_results", "phydist_vs_dS.csv"), row.names = FALSE)


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
write.csv(relax_summary, here("analysis_results", "relax_summary.csv"), row.names = FALSE)

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
################################################################################

# --------- ENC Prime -------------

model_ENC <- lmer(ENC_prime ~ Condition * TaxonPair + (1 | SOG),
                              data = codon_myrmecia_dup)
summary(model_ENC)
resid_panel(model_ENC)

# F-tests on main effects and interaction
anova(model_ENC)

# Per-pair marginal means and contrasts (Table 3)
emm_ENC_bypair <- emmeans(model_ENC, ~ Condition | TaxonPair)
print(emm_ENC_bypair)

contrasts_ENC_bypair <- contrast(emm_ENC_bypair, method = "pairwise",
                                 adjust = "bonferroni")
print(contrasts_ENC_bypair)

# Pooled marginal means across pairs
emm_ENC_pooled <- emmeans(model_ENC, ~ Condition)
print(emm_ENC_pooled)
contrast_ENC_pooled <- contrast(emm_ENC_pooled, method = "pairwise")
print(contrast_ENC_pooled)


# --------- GC3 -------------

model_GC3 <- lmer(GC3 ~ Condition * TaxonPair + (1 | SOG),
                              data = codon_myrmecia_dup)
summary(model_GC3)
resid_panel(model_GC3)

# F-tests on main effects and interaction
anova(model_GC3)

# Per-pair marginal means and contrasts (Table 3)
emm_GC3_bypair <- emmeans(model_GC3, ~ Condition | TaxonPair)
print(emm_GC3_bypair)

contrasts_GC3_bypair <- contrast(emm_GC3_bypair, method = "pairwise",
                                 adjust = "bonferroni")
print(contrasts_GC3_bypair)

# Pooled marginal means across pairs
emm_GC3_pooled <- emmeans(model_GC3, ~ Condition)
print(emm_GC3_pooled)
contrast_GC3_pooled <- contrast(emm_GC3_pooled, method = "pairwise")
print(contrast_GC3_pooled)

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
write.csv(codon_results, here("analysis_results", "codon_bias_contrasts.csv"),
          row.names = FALSE)


################################################################################
# ANALYSIS 5: GENE CONSISTENCY ACROSS LINEAGES
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
          here("analysis_results", "intensified_genes_GO_table.csv"),
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
          here("analysis_results", "relaxed_genes_GO_table.csv"),
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
  annotations,
  gene2GO_list,
  taxa_order,
  tree,
  phydist,
  # Analysis results
  phydist_vs_dS,
  dnds_results,
  relax_summary,
  codon_results,
  gene_consistency,
  consistency_summary,
  model_ENC,
  model_GC3,
  model_GC12,
  emm_ENC_pooled,
  emm_ENC_bypair,
  emm_GC3_pooled,
  emm_GC3_bypair,
  emm_GC12,
  file = here("Rscripts", "analysis_complete.RData")
)

