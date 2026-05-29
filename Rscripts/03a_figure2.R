library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/Desktop/Algae-Comparative-Genomics/")
dir.create("figures", showWarnings = FALSE)
load("Rscripts/analysis_complete.RData")

# Color scheme
col_lichen <- "#1bbc9b"
col_free   <- "#dca81b"
color_condition <- c("Lichen-forming" = col_lichen, "Free-living" = col_free)

# Clean taxon pair labels
pair_labels <- c(
  "Coccomyxa" = "Coccomyxa",
  "Symbiochloris" = "Symbiochloris",
  "Trebouxia_vs_Myrmecia" = "Trebouxia",
  "Asterochloris_vs_Myrmecia" = "Asterochloris"
)
pair_order <- c("Coccomyxa", "Symbiochloris",
                "Trebouxia_vs_Myrmecia", "Asterochloris_vs_Myrmecia")

# Prepare data
plot_data <- M4_codeml_filter %>%
  filter(TaxonPair %in% pair_order) %>%
  mutate(
    PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
    Condition = factor(Condition, levels = c("Lichen-forming", "Free-living"))
  )

################################################################################
#  Compute stats for annotations
################################################################################

# Read pre-computed contrasts from linear model / emmeans
contrasts_rates <- read.csv("analysis_results/absolute_rates_contrasts.csv")
contrasts_dnds  <- read.csv("analysis_results/dnds_contrasts.csv")
contrasts_all   <- rbind(
  contrasts_rates[, c("TaxonPair", "estimate", "p.value", "Metric")],
  contrasts_dnds[, c("TaxonPair", "estimate", "p.value", "Metric")]
)

make_stats <- function(contrasts_df, metric_name, digits = 4) {
  contrasts_df %>%
    filter(Metric == metric_name) %>%
    mutate(
      PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
      sig_label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE             ~ "ns"
      ),
      annotation = sig_label
    )
}

stats_dN    <- make_stats(contrasts_all, "dN", digits = 4)
stats_dS    <- make_stats(contrasts_all, "dS", digits = 4)
stats_omega <- make_stats(contrasts_all, "omega", digits = 4)

################################################################################
#  PANEL A: Phylogeny 
################################################################################

tree_text <- readLines("data/consensus_tree.newick.txt", warn = FALSE)
tree_text_std <- gsub("):(\\d+\\.\\d+)\\[(\\d+\\.\\d+)\\]", ")\\2:\\1", tree_text)
tree <- read.tree(text = tree_text_std)

tip_info <- data.frame(
  label = c("Tre_spC0010", "Ast_eri", "Myr_bis",
            "Coc_vir", "Coc_sub", "Sym_ret", "Sym_irr"),
  Species = c("Trebouxia sp. C0010", "Asterochloris erici",
              "Myrmecia bisecta*",
              "Coccomyxa viridis", "Coccomyxa subellipsoidea",
              "Symbiochloris reticulata", "Symbiochloris irregularis"),
  Lifestyle = c("Lichen-forming", "Lichen-forming", "Free-living",
                 "Lichen-forming", "Free-living", "Lichen-forming",
                 "Free-living"),
  is_reference = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

grp <- list(
  "Lichen-forming" = tip_info$label[tip_info$Lifestyle == "Lichen-forming"],
  "Free-living"    = tip_info$label[tip_info$Lifestyle == "Free-living"]
)
tree_grouped <- groupOTU(tree, grp)

p2a <- ggtree(tree_grouped, aes(color = group), size = 1.4,
              ladderize = TRUE) %<+% tip_info +
  geom_tiplab(aes(label = Species,
                  fontface = "italic"),
              size = 6.5, offset = 0.02, align = TRUE, color = "black",
              linetype = "dashed", linesize = 0.3, linecolour = "gray70") +
  geom_nodelab(aes(label = ifelse(!is.na(label) & label != "",
                                   paste0(round(as.numeric(label) * 100), "%"), "")),
               size = 5, hjust = 1.1, vjust = -0.4, color = "gray40") +
  scale_color_manual(
    values = c("Lichen-forming" = col_lichen, "Free-living" = col_free),
    breaks = c("Lichen-forming", "Free-living"),
    name = "Lifestyle"
  ) +
  annotate("text", x = 0.80, y = 0.35, label = "* shared free-living reference",
           size = 6.25, color = "gray30", hjust = 0, fontface = "italic") +
  annotate("segment", x = 0, xend = 0.1, y = 0.25, yend = 0.25,
           color = "black", linewidth = 0.7) +
  annotate("text", x = 0.05, y = 0.28, label = "0.1",
           size = 5, vjust = -0.3) +
  xlim(NA, 2.1) +
  theme_tree() +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 0, 5, 5)
  )

p2a

# ################################################################################
#  Function to create rate distribution panel
################################################################################

make_rate_panel <- function(df, yvar, stats_df, title, ylab) {
  y_upper <- quantile(df[[yvar]], 0.995, na.rm = TRUE) * 1.2

  ggplot(df, aes(x = PairLabel, y = .data[[yvar]], fill = Condition)) +
    geom_boxplot(width = 0.35, outlier.shape = NA, position = position_dodge(0.8),
                 linewidth = 0.4, alpha = 0.6, color = "gray30") +
    geom_point(alpha = 0.3, size = 0.8,
               position = position_jitterdodge(jitter.width = 0.2,
                                               dodge.width = 0.8)) +
    geom_text(data = stats_df,
              aes(x = PairLabel, y = y_upper, label = annotation),
              inherit.aes = FALSE, size = 6, fontface = "bold",
              color = "gray10", parse = FALSE) +
    scale_fill_manual(values = color_condition, name = "Lifestyle") +
    scale_y_continuous(limits = c(0, y_upper * 1.05), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(face = "italic", size = 16),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(size = 16, face = "bold"),
      plot.margin = margin(5, 8, 2, 5)
    )
}

#################################################################################
# 3 PANELS: B, C, D 
################################################################################
p2b <- make_rate_panel(plot_data, "omega", stats_omega,
                       expression(bold("dN/dS ("*omega*")")),
                       expression(dN/dS~(omega)))

p2c <- make_rate_panel(plot_data, "dN", stats_dN,
                       "Nonsynonymous Rate (dN)", "dN")

p2d <- make_rate_panel(plot_data, "dS", stats_dS,
                       "Synonymous Rate (dS)", "dS")


################################################################################
#  Combine into Figure 1
################################################################################

# Layout: A on the left spanning full height, B/C/D stacked on the right
right_panels <- (p2b / p2c / p2d) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1, "cm")
  )

fig2 <- p2a | right_panels
fig2 <- fig2 +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(size = 25, face = "bold")))

fig2

ggsave("figures/Figure2_new.png", fig2,
       width = 15, height = 13, dpi = 300)
ggsave("figures/Figure2_new.pdf", fig2,
       width = 15, height = 11)


