library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/Desktop/Algae-Comparative-Genomics/")
load("Rscripts/prepared_data.RData")

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

codon_pairs <- codon_myrmecia_dup %>%
  filter(TaxonPair %in% pair_order) %>%
  mutate(
    PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
    Condition = factor(Condition, levels = c("Lichen-forming", "Free-living"))
  )

# Wright (1990) expected curve: ENC_exp = 2 + s + 29/(s^2 + (1-s)^2)
wright_curve <- data.frame(GC3 = seq(0, 1, by = 0.002)) %>%
  mutate(ENC_exp = 2 + GC3 + 29 / (GC3^2 + (1 - GC3)^2))

p_nc <- ggplot(codon_pairs,
               aes(x = GC3, y = ENC, color = Condition)) +
  geom_line(data = wright_curve,
            aes(x = GC3, y = ENC_exp),
            inherit.aes = FALSE,
            color = "gray40", linewidth = 0.8, linetype = "solid") +
  geom_point(size = 0.6, alpha = 0.25, shape = 16) +
  scale_color_manual(values = color_condition, name = "Lifestyle") +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(30, 61)) +
  facet_wrap(~ PairLabel, nrow = 2) +
  labs(
    x = "GC3",
    y = "ENC"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "italic", size = 10),
    plot.margin      = margin(5, 8, 5, 5)
  )

ggsave("figures/Supp_Fig_nc_plot.png", p_nc,
       width = 8, height = 7, dpi = 600)
ggsave("figures/Supp_Fig_nc_plot.pdf", p_nc,
       width = 8, height = 7)
