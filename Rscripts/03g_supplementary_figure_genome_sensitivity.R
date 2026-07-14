################################################################################
# Supplementary Figure: Alternative Genome Choice Sensitivity Analysis
################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(patchwork)
library(here)

# Paths resolved relative to the project root via here::here()
dir.create(here("figures"), showWarnings = FALSE)

################################################################################
# 1. GENOME METADATA LOOKUP
################################################################################

genome_meta <- tibble(
  genome_id    = c("Ast_eri",               "Ast_glo",
                   "Tre_spC0010",           "Tre_spA1-2",
                   "Tre_spC0006",           "Tre_lyn",
                   "Tre_spTZW2008"),
  display_name = c("A. erici",    "A. glomerata",
                   "T. sp. C0010","T. sp. A1-2",
                   "T. sp. C0006",          "T. lynnae",
                   "T. sp. TZW2008"),
  genus        = c("Asterochloris",         "Asterochloris",
                   "Trebouxia",             "Trebouxia",
                   "Trebouxia",             "Trebouxia",
                   "Trebouxia"),
  focal_pair   = c("Asterochloris_vs_Myrmecia", "Asterochloris_vs_Myrmecia",
                   "Trebouxia_vs_Myrmecia",      "Trebouxia_vs_Myrmecia",
                   "Trebouxia_vs_Myrmecia",      "Trebouxia_vs_Myrmecia",
                   "Trebouxia_vs_Myrmecia"),
  is_primary   = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
)

genome_files <- c(
  "Ast_eri"       = here("data", "compiled_results_Ast_eri.csv"),
  "Ast_glo"       = here("data", "compiled_results_Ast_glo.csv"),
  "Tre_spC0010"   = here("data", "compiled_results_Tre_spC0010.csv"),
  "Tre_spA1-2"    = here("data", "compiled_results_Tre_spA1-2.csv"),
  "Tre_spC0006"   = here("data", "compiled_results_Tre_spC0006.csv"),
  "Tre_lyn"       = here("data", "compiled_results_Tre_lyn.csv"),
  "Tre_spTZW2008" = here("data", "compiled_results_Tre_spTZW2008.csv")
)

################################################################################
# 2. TAXA METADATA HELPER
################################################################################

# Condition set as factor with Lichen-forming first so emmeans pairwise gives
# Lichen-forming − Free-living (positive = higher in lichen-forming).
# TaxonPair uses first-match logic in case_when; Myrmecia (Tip 4) initially
# maps to Asterochloris_vs_Myrmecia and is later duplicated for Trebouxia.

add_taxa_metadata <- function(df) {
  df %>%
    mutate(
      Taxa = case_when(
        Tip == "1" ~ "Asterochloris",
        Tip == "2" ~ "Coccomyxa subellipsoidea",
        Tip == "3" ~ "Coccomyxa viridis",
        Tip == "4" ~ "Myrmecia bisecta",
        Tip == "5" ~ "Symbiochloris irregularis",
        Tip == "6" ~ "Symbiochloris reticulata",
        Tip == "7" ~ "Trebouxia"
      ),
      Condition = factor(
        case_when(
          Tip == "1" ~ "Lichen-forming",
          Tip == "2" ~ "Free-living",
          Tip == "3" ~ "Lichen-forming",
          Tip == "4" ~ "Free-living",
          Tip == "5" ~ "Free-living",
          Tip == "6" ~ "Lichen-forming",
          Tip == "7" ~ "Lichen-forming"
        ),
        levels = c("Lichen-forming", "Free-living")
      ),
      TaxonPair = case_when(
        Tip %in% c("1", "4") ~ "Asterochloris_vs_Myrmecia",
        Tip %in% c("2", "3") ~ "Coccomyxa",
        Tip %in% c("5", "6") ~ "Symbiochloris",
        Tip %in% c("4", "7") ~ "Trebouxia_vs_Myrmecia"
      )
    )
}

################################################################################
# 3. GENOME PROCESSING PIPELINE
################################################################################

process_genome <- function(genome_id, file_path, focal_pair) {

  message("Processing: ", genome_id)

  raw <- read.csv(file_path)

  # Parse "Node..Tip" branch format into separate columns
  raw <- raw %>%
    separate(branch, into = c("Node", "Tip"), sep = "\\.\\.", remove = FALSE)

  # Filter to M4 model, tip branches 1–7 only
  m4 <- raw %>%
    filter(Model == "M4", Tip %in% as.character(1:7)) %>%
    add_taxa_metadata()

  # Gene-level quality filter: remove any SOG where at least one branch fails
  failing_SOGs <- m4 %>%
    filter(dS <= 0.01 | dS > 3 | omega >= 10) %>%
    pull(SOG) %>%
    unique()

  m4_filtered <- m4 %>%
    filter(!(SOG %in% failing_SOGs))

  n_genes <- n_distinct(m4_filtered$SOG)
  message("  SOGs retained after filtering: ", n_genes)

  # Duplicate Myrmecia bisecta for both TaxonPair comparisons
  myr <- m4_filtered %>% filter(Taxa == "Myrmecia bisecta")

  m4_model <- m4_filtered %>%
    filter(Taxa != "Myrmecia bisecta") %>%
    bind_rows(
      myr %>% mutate(TaxonPair = "Asterochloris_vs_Myrmecia"),
      myr %>% mutate(TaxonPair = "Trebouxia_vs_Myrmecia")
    )

  # Fit lmer and extract focal contrast for each response variable
  results <- lapply(c("dN", "dS", "omega"), function(resp) {

    f <- as.formula(paste0(resp, " ~ Condition * TaxonPair + (1 | SOG)"))

    model <- tryCatch(
      suppressMessages(lmer(f, data = m4_model)),
      error = function(e) {
        message("  lmer failed for ", resp, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(model)) return(NULL)

    emm   <- emmeans(model, ~ Condition | TaxonPair)
    contr <- as.data.frame(contrast(emm, method = "pairwise", adjust = "none"))

    # Retain only the relevant taxon-pair contrast
    contr %>%
      filter(TaxonPair == focal_pair) %>%
      mutate(
        response  = resp,
        genome_id = genome_id,
        n_genes   = n_genes
      )
  })

  bind_rows(results)
}

################################################################################
# 4. RUN ALL GENOMES
################################################################################

all_results <- mapply(
  FUN        = process_genome,
  genome_id  = names(genome_files),
  file_path  = unname(genome_files),
  focal_pair = genome_meta$focal_pair[match(names(genome_files), genome_meta$genome_id)],
  SIMPLIFY   = FALSE
) %>%
  bind_rows()

# Verify contrast direction (should be "Lichen-forming - Free-living")
print(unique(all_results$contrast))

# Gene retention per genome
print(
  all_results %>%
    distinct(genome_id, n_genes) %>%
    left_join(select(genome_meta, genome_id, display_name), by = "genome_id")
)

################################################################################
# 5. PREPARE PLOT DATA
################################################################################

plot_data <- all_results %>%
  left_join(genome_meta, by = "genome_id") %>%
  mutate(
    ci_lo = estimate - 1.96 * SE,
    ci_hi = estimate + 1.96 * SE,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE             ~ ""
    ),
    response_label = factor(
      response,
      levels = c("dN",    "dS",    "omega"),
      labels = c("\u0394dN", "\u0394dS", "\u0394\u03c9")
    )
  )

# Y-axis factor: levels ordered bottom-to-top within each genus group
# Primary genome sits at the top of each group
y_levels <- c(
  # Trebouxia group
  "T. sp. TZW2008", "T. lynnae", "T. sp. C0006",
  "T. sp. A1-2",    "T. sp. C0010",
  # Asterochloris group
  "A. glomerata",   "A. erici"
)

plot_data <- plot_data %>%
  mutate(
    display_name = factor(display_name, levels = y_levels),
    genus        = factor(genus, levels = c("Asterochloris", "Trebouxia"))
  )

################################################################################
# 6. BUILD FOREST PLOT
################################################################################

genus_colors <- c("Asterochloris" = "#D55E00", "Trebouxia" = "#0072B2")

forest_plot <- ggplot(
  plot_data,
  aes(x = estimate, y = display_name, color = genus)
) +
  # Reference line
  geom_vline(
    xintercept = 0, linetype = "dashed",
    color = "grey45", linewidth = 0.45
  ) +
  # 95% CI error bars
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi),
    height = 0.3, linewidth = 0.65
  ) +
  # Alternatives: open circles
  geom_point(
    data  = filter(plot_data, !is_primary),
    shape = 21, size = 3.0, stroke = 1.0,
    fill  = "white"
  ) +
  # Primary genomes: filled circles (slightly larger)
  geom_point(
    data  = filter(plot_data, is_primary),
    shape = 16, size = 3.6
  ) +
  # Significance stars (right of CI upper bound)
  geom_text(
    aes(x = ci_hi, label = sig),
    vjust = 1.5,
    size  = 3.5, color = "black",
    angle = 90,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = genus_colors,
    guide  = "none"
  ) +
  # Facet: column = response variable, row = genus group
  # space = "free_y" proportionally sizes rows to n_genomes per genus
  facet_grid(
    genus ~ response_label,
    scales = "free",
    space  = "free_y"
  ) +
  labs(
    title = "Sensitivity Analysis: Trebouxia and Asterochloris genome choice",
    x = "Estimated difference (Lichen-forming \u2212 Free-living)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey94", color = "grey70"),
    strip.text.x       = element_text(face = "bold",        size = 11),
    strip.text.y.right = element_text(face = "bold.italic", size = 9, angle = 270),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 9, face = "italic"),
    axis.title.x       = element_text(size = 10, margin = margin(t = 6)),
    plot.margin        = margin(6, 18, 6, 6),
    plot.title         = element_text(face = "bold", size = 14)
  )


print(forest_plot)

################################################################################
# 7. SAVE OUTPUT
################################################################################

ggsave(
  here("figures", "FigureS4_genome_sensitivity.pdf"),
  forest_plot,
  width  = 10,
  height = 5,
  device = cairo_pdf
)

ggsave(
  here("figures", "FigureS4_genome_sensitivity.png"),
  forest_plot,
  width  = 10,
  height = 5,
  dpi    = 300
)

