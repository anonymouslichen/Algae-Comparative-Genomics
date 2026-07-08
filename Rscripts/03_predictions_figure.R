library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/Desktop/Algae-Comparative-Genomics/")
dir.create("figures", showWarnings = FALSE)

# ── Color scheme (matches 03a) ────────────────────────────────────────────────
col_lichen <- "#1bbc9b"
col_free   <- "#dca81b"
color_condition <- c("Lichen-forming" = col_lichen, "Free-living" = col_free)

pair_labels <- c(
  "Coccomyxa"              = "Coccomyxa",
  "Symbiochloris"          = "Symbiochloris",
  "Trebouxia_vs_Myrmecia"  = "Trebouxia",
  "Asterochloris_vs_Myrmecia" = "Asterochloris"
)
pair_order <- c("Coccomyxa", "Symbiochloris",
                "Trebouxia_vs_Myrmecia", "Asterochloris_vs_Myrmecia")

# ── Simulate mock data showing predicted pattern ──────────────────────────────
set.seed(42)
n_genes <- 600

make_pair <- function(pair_id, n = n_genes,
                      lichen_omega_mu, free_omega_mu,
                      lichen_dN_mu,    free_dN_mu,
                      shared_dS_mu) {
  bind_rows(
    data.frame(
      TaxonPair = pair_id,
      Condition = "Lichen-forming",
      omega = rlnorm(n, log(lichen_omega_mu), 0.55),
      dN    = rlnorm(n, log(lichen_dN_mu),    0.55),
      dS    = rlnorm(n, log(shared_dS_mu),    0.40)
    ),
    data.frame(
      TaxonPair = pair_id,
      Condition = "Free-living",
      omega = rlnorm(n, log(free_omega_mu),  0.55),
      dN    = rlnorm(n, log(free_dN_mu),     0.55),
      dS    = rlnorm(n, log(shared_dS_mu),   0.40)
    )
  )
}

mock_data <- bind_rows(
  make_pair("Coccomyxa",
            lichen_omega_mu = 0.22, free_omega_mu = 0.12,
            lichen_dN_mu    = 0.09, free_dN_mu    = 0.045,
            shared_dS_mu    = 0.42),
  make_pair("Symbiochloris",
            lichen_omega_mu = 0.20, free_omega_mu = 0.11,
            lichen_dN_mu    = 0.08, free_dN_mu    = 0.040,
            shared_dS_mu    = 0.38),
  make_pair("Trebouxia_vs_Myrmecia",
            lichen_omega_mu = 0.24, free_omega_mu = 0.13,
            lichen_dN_mu    = 0.10, free_dN_mu    = 0.050,
            shared_dS_mu    = 0.55),
  make_pair("Asterochloris_vs_Myrmecia",
            lichen_omega_mu = 0.21, free_omega_mu = 0.12,
            lichen_dN_mu    = 0.085, free_dN_mu   = 0.042,
            shared_dS_mu    = 0.55)
) %>%
  # cap extreme outliers so y-axes stay clean
  group_by(TaxonPair, Condition) %>%
  mutate(
    omega = pmin(omega, quantile(omega, 0.995)),
    dN    = pmin(dN,    quantile(dN,    0.995)),
    dS    = pmin(dS,    quantile(dS,    0.995))
  ) %>%
  ungroup() %>%
  mutate(
    PairLabel = factor(pair_labels[TaxonPair], levels = pair_labels[pair_order]),
    Condition = factor(Condition, levels = c("Lichen-forming", "Free-living"))
  )

# ── Panel-building function (matches 03a) ─────────────────────────────────────
make_pred_panel <- function(df, yvar, title, ylab) {
  y_upper <- quantile(df[[yvar]], 0.995, na.rm = TRUE) * 1.2

  ggplot(df, aes(x = PairLabel, y = .data[[yvar]], fill = Condition)) +
    geom_boxplot(width = 0.35, outlier.shape = NA,
                 position = position_dodge(0.8),
                 linewidth = 0.4, alpha = 0.6, color = "gray30") +
    geom_point(alpha = 0.3, size = 0.8,
               position = position_jitterdodge(jitter.width = 0.2,
                                               dodge.width  = 0.8)) +
    scale_fill_manual(values = color_condition, name = "Lifestyle") +
    scale_y_continuous(limits = c(0, y_upper * 1.05), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x  = element_text(face = "italic", size = 16),
      axis.text.y  = element_text(size = 12),
      axis.title.y = element_text(size = 15),
      plot.title   = element_text(size = 16, face = "bold"),
      plot.margin  = margin(5, 8, 2, 5)
    )
}

# ── Build panels ──────────────────────────────────────────────────────────────
p_omega <- make_pred_panel(mock_data, "omega",
                           expression(bold("dN/dS ("*omega*")")),
                           expression(dN/dS~(omega)))

p_dN    <- make_pred_panel(mock_data, "dN",
                           "Nonsynonymous Rate (dN)", "dN")

p_dS    <- make_pred_panel(mock_data, "dS",
                           "Synonymous Rate (dS)", "dS")

# ── Combine ───────────────────────────────────────────────────────────────────
fig_pred <- (p_omega / p_dN / p_dS) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.key.size  = unit(1, "cm")
  )

fig_pred <- fig_pred +
  plot_annotation(
    title = "Predictions",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    )
  )

fig_pred

ggsave("figures/Figure_predictions.png", fig_pred,
       width = 8, height = 13, dpi = 300)
ggsave("figures/Figure_predictions.pdf", fig_pred,
       width = 8, height = 13)
