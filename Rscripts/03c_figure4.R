library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggalluvial)

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

color_result <- c(
  "Strengthened" = "#e74c3c",
  "Relaxed" = "#3498db",
  "Not Significant" = "#95a5a6"
)

################################################################################
# PANEL A: K parameter distributions by lineage
################################################################################

y_max_K <- quantile(relax_focal_filtered$Relaxation.Parameter..K., 0.995) * 1.1

p4a <- ggplot(relax_focal_filtered,
              aes(x = Pair, y = Relaxation.Parameter..K., fill = Pair)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_boxplot(width = 0.35, outlier.shape = NA, linewidth = 0.5,
               alpha = 0.6, color = "gray30") +
  scale_fill_manual(values = color_pair) +
  coord_cartesian(ylim = c(0, y_max_K)) +
  annotate("text", x = 0.55, y = 0.8, label = "K = 1",
           size = 2.8, color = "red", hjust = 0) +
  labs(
    title = "K Parameter Distribution by Lineage",
    x = NULL, y = "K Parameter"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 8, 5, 5)
  )

p4a

################################################################################
# PANEL B: Alluvial diagram of gene trajectories
################################################################################

# Genes present in all 4 lineages
genes_all4 <- relax_focal_filtered %>%
  group_by(SOG) %>%
  filter(n_distinct(Pair) == 4) %>%
  ungroup() %>%
  dplyr::select(SOG, Pair, Result) %>%
  mutate(Result = factor(Result,
    levels = c("Strengthened", "Not Significant", "Relaxed")))

n_genes_all4 <- n_distinct(genes_all4$SOG)

p4b <- ggplot(genes_all4,
              aes(x = Pair, stratum = Result, alluvium = SOG, fill = Result)) +
  geom_flow(alpha = 0.4) +
  geom_stratum(alpha = 0.8, width = 0.35) +
  geom_text(stat = "stratum", aes(label = after_stat(count)),
            size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = color_result, name = "Selection Change") +
  labs(
    title = paste0("Gene Trajectories Across Lineages"),
    x = NULL, y = "Number of Genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "italic", size = 10),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 8, 5, 5)
  )

p4b

################################################################################
# Combine into Figure 4
################################################################################

fig4 <- (p4a | p4b) +
  plot_layout(widths = c(1, 1.8)) +
  plot_annotation(tag_levels = "A")

ggsave("figures/Figure4_new.png", fig4,
       width = 15, height = 8, dpi = 600)
ggsave("figures/Figure4_new.pdf", fig4,
       width = 15, height = 8)



