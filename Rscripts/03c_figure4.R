library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(here)

# Paths resolved relative to the project root via here::here()
dir.create(here("figures"), showWarnings = FALSE)
load(here("Rscripts", "analysis_complete.RData"))
load(here("Rscripts", "prepared_data.RData"))

# Colors
color_result <- c(
  "Intensified" = "#e74c3c",
  "Relaxed" = "#3498db",
  "Not Significant" = "#95a5a6"
)

pair_labels <- c(
  "Coccomyxa"      = "Coccomyxa viridis",
  "Symbiochloris"  = "Symbiochloris reticulata",
  "Trebouxia"      = "Trebouxia sp.C0010",
  "Asterochloris"  = "Asterochloris erici"
)

################################################################################
# FIGURE 4: Volcano plots (ln(K) vs -log10(padj)) for each taxa test
################################################################################

make_volcano <- function(pair_name, df, colors) {
  dat <- df %>%
    filter(Pair == pair_name) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      Result = factor(Result, levels = c("Intensified", "Not Significant", "Relaxed"))
    )

  n_sig <- sum(dat$Result != "Not Significant", na.rm = TRUE)

  ggplot(dat, aes(x = ln_K, y = neg_log10_padj, color = Result)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "gray40", linewidth = 0.6) +
    scale_color_manual(values = colors, name = "Selection Change") +
    labs(
      title = pair_labels[pair_name],
      x = "ln(K)",
      y = expression(-log[10](p[adj]))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position  = "bottom",
      plot.title = element_text(face = "bold.italic", size = 16),
      axis.text.x      = element_text(size = 12),
      axis.text.y      = element_text(size = 12),
      axis.title.x     = element_text(size = 15),
      axis.title.y     = element_text(size = 15),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 8, 5, 5)
    )
}

p_coc <- make_volcano("Coccomyxa",     relax_focal_filtered, color_result)
p_sym <- make_volcano("Symbiochloris", relax_focal_filtered, color_result)
p_tre <- make_volcano("Trebouxia",     relax_focal_filtered, color_result)
p_ast <- make_volcano("Asterochloris", relax_focal_filtered, color_result)
p_coc

# Combine into 2x2 grid with shared legend
fig4 <- (p_coc | p_sym) / (p_tre | p_ast) +
  plot_layout(guides = "collect")  &
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1, "cm")
  )

fig4 <- fig4 + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 25, face = "bold")))

fig4

ggsave(here("figures", "Figure4.png"), fig4,
       width = 15, height =10, dpi = 600)
ggsave(here("figures", "Figure4.pdf"), fig4,
       width = 15, height = 10)
