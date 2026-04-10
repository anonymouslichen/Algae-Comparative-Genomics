library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("~/Desktop/Claude")
load("Rscripts/analysis_complete.RData")
load("Rscripts/prepared_data.RData")

# Colors
color_pair <- c(
  "Coccomyxa" = "#E69F00",
  "Symbiochloris" = "#56B4E9",
  "Trebouxia" = "#009E73",
  "Asterochloris" = "#CC79A7"
)

################################################################################
# PANEL A: Translation machinery K across lineages
################################################################################

# Identify translation-related genes via GO annotations
translation_go <- c("GO:0006412")


translation_sogs <- names(gene2GO_list)[
  sapply(gene2GO_list, function(x) any(x %in% translation_go))
]

# Filter RELAX data to translation genes
relax_translation <- relax_focal_filtered %>%
  filter(SOG %in% translation_sogs)

n_trans <- n_distinct(relax_translation$SOG)

# Genome-wide medians per lineage
genome_medians <- relax_focal_filtered %>%
  group_by(Pair) %>%
  summarise(median_K = median(Relaxation.Parameter..K., na.rm = TRUE),
            .groups = "drop")

# Per-lineage test: translation K vs genome-wide K
stat_labels <- relax_translation %>%
  group_by(Pair) %>%
  summarise(
    median_trans = median(Relaxation.Parameter..K., na.rm = TRUE),
    pval = wilcox.test(
      Relaxation.Parameter..K.,
      relax_focal_filtered$Relaxation.Parameter..K.[
        relax_focal_filtered$Pair == Pair[1]]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  )

y_max <- relax_translation %>%
  group_by(Pair) %>%
  summarise(q = quantile(Relaxation.Parameter..K., 0.995, na.rm = TRUE)) %>%
  pull(q) %>% max () * 1.15

p5a <- ggplot(relax_translation,
              aes(x = Pair, y = Relaxation.Parameter..K., fill = Pair)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_boxplot(width = 0.35, outlier.shape = NA, linewidth = 0.5,
               alpha = 0.6, color = "gray30") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "black") +
  # Genome-wide median as diamond
  geom_point(data = genome_medians,
             aes(x = Pair, y = median_K),
             shape = 18, size = 4, color = "gray20", inherit.aes = FALSE) +
  # Significance labels
  geom_text(data = stat_labels,
            aes(x = Pair, y = Inf, label = sig),
            inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 1.5) +
  scale_fill_manual(values = color_pair) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  annotate("text", x = 0.55, y = 0.25, label = "K = 1",
           size = 2.8, color = "red", hjust = 0) +
  labs(
    title = paste0("Translation Machinery K"),
    x = NULL, y = "K Parameter"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 15, 5, 5)
  )

p5a

################################################################################
# PANEL B: McDonald-Kreitman NI distribution
################################################################################

mk <- MK_results %>%
  mutate(SOG = sub("_codon_alignment$", "", Gene))

# Remove infinite/NA NI values
mk_clean <- mk %>% filter(is.finite(NI) & !is.na(NI))

# Use genome-wide NI_TG and bootstrap CI from statistical analyses
mk_summary <- read.csv("analysis_results/mk_summary.csv")
NI_TG    <- mk_summary$NI_TG
ci_lower <- mk_summary$CI_lower
ci_upper <- mk_summary$CI_upper

# Trim extreme values for visualization
x_upper <- quantile(mk_clean$NI, 0.99)

p5b <- ggplot(mk_clean, aes(x = NI)) +
  # Shaded bootstrap 95% CI region
  annotate("rect", xmin = ci_lower, xmax = ci_upper,
           ymin = -Inf, ymax = Inf,
           fill = "#e74c3c", alpha = 0.12) +
  # Neutral expectation
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40",
             linewidth = 0.8) +
  # NI_TG (Stoletzki & Eyre-Walker genome-wide estimate)
  geom_vline(xintercept = NI_TG, color = "#e74c3c", linewidth = 1) +
  # Density
  geom_density(fill = "#3498db", alpha = 0.4, color = "#2c3e50",
               linewidth = 0.6) +
  # Annotations
  annotate("text", x = NI_TG - 0.02, y = Inf,
           label = sprintf("NI = %.3f", NI_TG),
           size = 3.5, fontface = "bold", color = "#e74c3c",
           hjust = 1.2, vjust = 2) +
  annotate("text", x = 1.02, y = Inf,
           label = "NI = 1\n(neutral)",
           size = 3, color = "gray40", hjust = 0, vjust = 2) +
  scale_x_continuous(limits = c(0, x_upper), expand = c(0, 0)) +
  labs(
    title = "McDonald-Kreitman Neutrality Index",
    x = "Neutrality Index (NI)",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 8, 5, 5)
  )

p5b

################################################################################
# Combine into Figure 5
################################################################################

fig5_new <- (p5a | p5b) +
  plot_annotation(tag_levels = "A")

ggsave("figures/Figure5_new.png", fig5_new,
       width = 15, height = 7, dpi = 600)
ggsave("figures/Figure5_new.pdf", fig5_new,
       width = 15, height = 7)

