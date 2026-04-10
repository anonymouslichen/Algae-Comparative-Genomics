################################################################################
# Supplementary Figure: dN/dS (omega) between lifestyles
#   Panel A - Free-ratio model (M4): pooled overall comparison (ns)
#   Panel B - Two-ratio model (M2): global foreground vs background omega
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("~/Desktop/Claude")
dir.create("figures", showWarnings = FALSE)
load("Rscripts:Data/analysis_complete.RData")

# Color scheme (matches Figure 2)
col_lichen <- "#1bbc9b"
col_free   <- "#dca81b"
color_condition <- c("Lichen-forming" = col_lichen, "Free-living" = col_free)

################################################################################
# PANEL A: M4 free-ratio model — pooled omega comparison
################################################################################

# Use the filtered M4 data already in the RData file
plot_data_M4 <- M4_codeml_filter %>%
  mutate(Condition = factor(Condition, levels = c("Lichen-forming", "Free-living")))

# Compute gene-level median omega per condition, then paired Wilcoxon
# (each gene contributes one median per condition, removing repeated-measures issue)
M4_gene_omega <- M4_codeml_filter %>%
  group_by(SOG, Condition) %>%
  summarise(med_omega = median(omega, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = med_omega) %>%
  filter(!is.na(`Lichen-forming`), !is.na(`Free-living`))

M4_wilcox_p <- wilcox.test(
  M4_gene_omega$`Lichen-forming`,
  M4_gene_omega$`Free-living`,
  paired = TRUE
)$p.value

M4_annot <- dplyr::case_when(
  M4_wilcox_p < 0.001 ~ "***",
  M4_wilcox_p < 0.01  ~ "**",
  M4_wilcox_p < 0.05  ~ "*",
  TRUE                 ~ "ns"
)

y_upper_M4 <- quantile(plot_data_M4$omega, 0.995, na.rm = TRUE) * 1.2

p_A <- ggplot(plot_data_M4, aes(x = Condition, y = omega, fill = Condition)) +
  geom_boxplot(width = 0.35, outlier.shape = NA,
               linewidth = 0.4, alpha = 0.6, color = "gray30") +
  geom_point(alpha = 0.15, size = 0.5,
             position = position_jitter(width = 0.15, seed = 42)) +
  annotate("segment",
           x = 1, xend = 2,
           y = y_upper_M4 * 0.93, yend = y_upper_M4 * 0.93,
           color = "gray50", linewidth = 0.4) +
  annotate("text",
           x = 1.5, y = y_upper_M4,
           label = M4_annot,
           size = 5, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = color_condition, name = "Lifestyle") +
  scale_y_continuous(limits = c(0, y_upper_M4 * 1.05), expand = c(0, 0)) +
  labs(
    title = expression(bold("Free-ratio model (M4)")),
    subtitle = paste0("Paired Wilcoxon, p = ", signif(M4_wilcox_p, 3)),
    x = NULL,
    y = expression(dN/dS ~ (omega))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.margin = margin(5, 15, 2, 5)
  )

################################################################################
# PANEL B: M2 two-ratio model — global foreground vs background omega
################################################################################

# Load raw codeml data to extract M2 results
codeml_raw <- read.csv("data/compiled_results_pooled.csv")

codeml_raw <- codeml_raw %>%
  tidyr::separate(branch, into = c("Node", "Tip"), sep = "\\.\\.",
                  remove = FALSE) %>%
  mutate(Tip = as.integer(Tip))

# Filter to M2 tip rows only (Tip 1–7 = the 7 terminal taxa)
# In the two-ratio model: Lichen-forming tips (1,3,6,7) share omega_foreground;
# Free-living tips (2,4,5) share omega_background
M2_tips <- codeml_raw %>%
  filter(Model == "M2", Tip %in% 1:7) %>%
  mutate(
    Condition = case_when(
      Tip %in% c(1, 3, 6, 7) ~ "Lichen-forming",
      Tip %in% c(2, 4, 5)    ~ "Free-living"
    )
  )

# Apply same filtering thresholds as M4 (dS <= 0.01 or > 3, omega >= 10)
failing_SOGs_M2 <- M2_tips %>%
  filter(dS <= 0.01 | dS > 3 | omega >= 10) %>%
  pull(SOG) %>%
  unique()

M2_tips_filter <- M2_tips %>%
  filter(!(SOG %in% failing_SOGs_M2))

# Deduplicate to one omega per gene per condition
# (all lichen tips share one omega per gene; all free-living tips share one omega)
M2_gene_omega <- M2_tips_filter %>%
  dplyr::select(SOG, Condition, omega) %>%
  distinct(SOG, Condition, .keep_all = TRUE)

# Paired Wilcoxon: foreground (lichen) vs background (free-living) per gene
M2_omega_wide <- M2_gene_omega %>%
  pivot_wider(names_from = Condition, values_from = omega) %>%
  filter(!is.na(`Lichen-forming`), !is.na(`Free-living`))

M2_wilcox_p <- wilcox.test(
  M2_omega_wide$`Lichen-forming`,
  M2_omega_wide$`Free-living`,
  paired = TRUE
)$p.value

M2_annot <- dplyr::case_when(
  M2_wilcox_p < 0.001 ~ "***",
  M2_wilcox_p < 0.01  ~ "**",
  M2_wilcox_p < 0.05  ~ "*",
  TRUE                 ~ "ns"
)

plot_data_M2 <- M2_gene_omega %>%
  mutate(Condition = factor(Condition, levels = c("Lichen-forming", "Free-living")))

y_upper_M2 <- quantile(plot_data_M2$omega, 0.995, na.rm = TRUE) * 1.2

p_B <- ggplot(plot_data_M2, aes(x = Condition, y = omega, fill = Condition)) +
  geom_boxplot(width = 0.35, outlier.shape = NA,
               linewidth = 0.4, alpha = 0.6, color = "gray30") +
  geom_point(alpha = 0.15, size = 0.5,
             position = position_jitter(width = 0.15, seed = 42)) +
  annotate("segment",
           x = 1, xend = 2,
           y = y_upper_M2 * 0.93, yend = y_upper_M2 * 0.93,
           color = "gray50", linewidth = 0.4) +
  annotate("text",
           x = 1.5, y = y_upper_M2,
           label = M2_annot,
           size = 5, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = color_condition, name = "Lifestyle") +
  scale_y_continuous(limits = c(0, y_upper_M2 * 1.05), expand = c(0, 0)) +
  labs(
    title = expression(bold("Two-ratio model (M2)")),
    subtitle = paste0("Paired Wilcoxon, p = ", signif(M2_wilcox_p, 3)),
    x = NULL,
    y = expression(dN/dS ~ (omega))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.margin = margin(5, 15, 2, 5)
  )

################################################################################
# Combine panels and save
################################################################################

fig_supp <- (p_A | p_B) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_supp <- fig_supp +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 18, face = "bold"))
  )

fig_supp

ggsave("figures/FigureS2_omega_models.png", fig_supp,
       width = 10, height = 6, dpi = 600)
ggsave("figures/FigureS2_omega_models.pdf", fig_supp,
       width = 10, height = 6)
