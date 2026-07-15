library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(here)

# Paths resolved relative to the project root via here::here()
dir.create(here("figures"), showWarnings = FALSE)
load(here("Rscripts", "analysis_complete.RData"))
load(here("Rscripts", "prepared_data.RData"))

# Colors and labels
col_lichen <- "#1bbc9b"
col_free   <- "#dca81b"
color_condition <- c("Lichen-forming" = col_lichen, "Free-living" = col_free)

pair_labels <- c(
  "Coccomyxa"                  = "Coccomyxa",
  "Symbiochloris"              = "Symbiochloris",
  "Trebouxia_vs_Myrmecia"      = "Trebouxia",
  "Asterochloris_vs_Myrmecia"  = "Asterochloris"
)
pair_order <- c("Coccomyxa", "Symbiochloris",
                "Trebouxia_vs_Myrmecia", "Asterochloris_vs_Myrmecia")

# Prepare codon data with clean pair labels
codon_pairs <- codon_myrmecia_dup %>%
  filter(TaxonPair %in% pair_order) %>%
  mutate(
    PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
    Condition = factor(Condition, levels = c("Lichen-forming", "Free-living"))
  )

################################################################################
# Significance annotations from LME contrasts (figure 2 style)
################################################################################

codon_contrasts <- read.csv(here("analysis_results", "codon_bias_contrasts.csv"))

make_stats_codon <- function(metric_name) {
  codon_contrasts %>%
    filter(Metric == metric_name) %>%
    mutate(
      PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
      annotation = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE             ~ "ns"
      )
    )
}

stats_ENCprime <- make_stats_codon("ENC_prime")

################################################################################
# PANEL A: ENC' boxplots by taxon pair (figure 2 style, no jitter)
################################################################################

y_upper_enc <- quantile(codon_pairs$ENC_prime, 0.995, na.rm = TRUE) * 1.05

p3a <- ggplot(codon_pairs,
              aes(x = PairLabel, y = ENC_prime, fill = Condition)) +
  geom_boxplot(width = 0.35, outlier.shape = NA,
               position = position_dodge(0.8),
               linewidth = 0.4, alpha = 0.6, color = "gray30") +
  geom_point(alpha = 0.1, size = 0.5,
             position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.8)) +
  geom_text(data = stats_ENCprime,
            aes(x = PairLabel, y = y_upper_enc, label = annotation),
            inherit.aes = FALSE, size = 4, fontface = "bold",
            color = "gray20") +
  
  scale_fill_manual(values = color_condition, name = "Lifestyle") +
  scale_y_continuous(limits = c(NA, y_upper_enc * 1.05), expand = c(0, 0)) +
  labs(title = "Effective Number of Codons (ENC')", x = NULL, y = "ENC'") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(face = "italic", size = 16),
    axis.text.y      = element_text(size = 12),
    axis.title.x     = element_text(size = 15),
    axis.title.y     = element_text(size = 15),
    plot.title         = element_text(size = 16, face = "bold"),
    plot.margin        = margin(5, 8, 2, 5)
  )

p3a

################################################################################
# PANEL B: Centroid neutrality plot – emmeans ± 95 % CI from LME
################################################################################

emm_GC3_df  <- as.data.frame(emm_GC3_bypair)
emm_GC12_df <- as.data.frame(emm_GC12)

centroids_lme <- emm_GC3_df %>%
  dplyr::select(TaxonPair, Condition,
                GC3    = emmean,
                GC3_lo = asymp.LCL,
                GC3_hi = asymp.UCL) %>%
  left_join(
    emm_GC12_df %>%
      dplyr::select(TaxonPair, Condition,
                    GC12    = emmean,
                    GC12_lo = asymp.LCL,
                    GC12_hi = asymp.UCL),
    by = c("TaxonPair", "Condition")
  ) %>%
  filter(TaxonPair %in% pair_order) %>%
  mutate(
    PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
    Condition = factor(Condition, levels = c("Lichen-forming", "Free-living"))
  )

p3b <- ggplot(centroids_lme,
              aes(x = GC3, y = GC12, color = Condition)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "gray50", linewidth = 0.8) +
  geom_errorbar(aes(ymin = GC12_lo, ymax = GC12_hi),
                width = 0.003, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = GC3_lo, xmax = GC3_hi),
                 height = 0.003, linewidth = 0.7) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.8) +
  #scale_x_continuous(limits = c(0.58, 0.80)) +
  #scale_y_continuous(limits = c(0.48, 0.54)) +
  scale_color_manual(values = color_condition, name = "Lifestyle") +
  facet_wrap(~ PairLabel, nrow = 2, scales = "free") +
  labs(
    title = "Neutrality Plot: GC12 vs GC3",
    x     = "GC3",
    y     = "GC12"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "italic", size = 16),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),
    axis.title.x     = element_text(size = 15),
    axis.title.y     = element_text(size = 15),
    plot.title       = element_text(size = 16, face = "bold"),
    plot.margin      = margin(5, 8, 5, 5)
  )

p3b

################################################################################
# Combine into Figure 3
################################################################################

fig3 <- (p3a | p3b) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1, "cm")
  )

fig3 <- fig3 + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 25, face = "bold")))

ggsave(here("figures", "Figure3.png"), fig3,
       width = 15, height = 9, dpi = 600)
ggsave(here("figures", "Figure3.pdf"), fig3,
       width = 15, height = 9)

ggsave(here("figures", "FigureENC.png"), p3a,
       width = 9, height = 9, dpi = 600)
ggsave(here("figures", "FigureGC3.png"), p3b,
       width = 9, height = 9, dpi = 600)


