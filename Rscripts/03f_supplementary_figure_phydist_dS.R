################################################################################
# Supplementary Figure: Correlation between emmeans delta dS and
# phylogenetic distance across taxon pairs (species-level, n = 4)
################################################################################

library(ggplot2)
library(dplyr)
library(ggrepel)
library(here)

# Path relative to the project root via here::here()
dir.create(here("figures"), showWarnings = FALSE)

# Pair labels
pair_labels <- c(
  "Coccomyxa"                 = "Coccomyxa",
  "Symbiochloris"             = "Symbiochloris",
  "Trebouxia_vs_Myrmecia"     = "Trebouxia",
  "Asterochloris_vs_Myrmecia" = "Asterochloris"
)

################################################################################
# Load data (both computed in script 02)
################################################################################

# Phylogenetic distances (tree-derived) and emmeans delta dS per pair
phydist_data <- read.csv(here("analysis_results", "phydist_vs_dS.csv"))

# Emmeans contrasts for dS — provides SE for error bars and significance
dS_contrasts <- read.csv(here("analysis_results", "dnds_contrasts.csv")) %>%
  filter(Metric == "dS") %>%
  dplyr::select(TaxonPair, SE, p.value) %>%
  rename(Pair = TaxonPair)

# Join
plot_data <- phydist_data %>%
  left_join(dS_contrasts, by = "Pair") %>%
  mutate(
    PairLabel    = pair_labels[Pair],
    CompType     = ifelse(Pair %in% c("Coccomyxa", "Symbiochloris"),
                          "Within genus", "Across genera"),
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )


################################################################################
# Build figure
################################################################################

# Y-axis range accounting for error bars
y_lo <- min(plot_data$delta_dS - 1.96 * plot_data$SE) * 1.15
y_hi <- max(plot_data$delta_dS + 1.96 * plot_data$SE) * 1.4

fig_phydist <- ggplot(plot_data, aes(x = PhyDist, y = delta_dS, shape = CompType)) +
  # Reference line at 0
  geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5,
             linetype = "dashed") +
  # 95% CI error bars (emmeans SE × 1.96)
  geom_errorbar(
    aes(ymin = delta_dS - 1.96 * SE,
        ymax = delta_dS + 1.96 * SE),
    width = 0.02, color = "gray30", linewidth = 0.6
  ) +
  # Points — shape mapped to comparison type
  geom_point(size = 4, color = "gray15", fill = "white", stroke = 1.4) +
  # Significance labels above each point
  geom_text(
    aes(y = delta_dS + 1.96 * SE + 0.03, label = sig_label),
    size = 4.5, fontface = "bold", color = "gray20"
  ) +
  # Pair labels (ggrepel)
  geom_label_repel(
    aes(label = PairLabel),
    fontface      = "italic",
    size          = 4,
    nudge_x       = 0.03,
    box.padding   = 0.5,
    point.padding = 0.3,
    min.segment.length = 0.2,
    label.size    = NA,
    fill          = alpha("white", 0.8),
    color         = "gray15"
  ) +
  scale_shape_manual(
    values = c("Within genus" = 21, "Across genera" = 24),
    name   = "Comparison type"
  ) +
  scale_x_continuous(
    limits = c(0.45, 1.12),
    expand = expansion(mult = 0.05)
  ) +
  scale_y_continuous(
    limits = c(y_lo, y_hi),
    expand = c(0, 0)
  ) +
  labs(
    x     = "Phylogenetic distance (substitutions per site)",
    y     = expression(Delta * "dS  (lichen-forming \u2212 free-living)"),
    title = expression(bold(Delta * "dS does not scale with phylogenetic distance"))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 11),
    axis.title       = element_text(size = 12),
    plot.margin      = margin(10, 15, 5, 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10)
  )

fig_phydist

################################################################################
# Save
################################################################################

ggsave(here("figures", "FigureS3_phydist_dS.png"), fig_phydist,
       width = 7, height = 6, dpi = 600)
ggsave(here("figures", "FigureS3_phydist_dS.pdf"), fig_phydist,
       width = 7, height = 6)
