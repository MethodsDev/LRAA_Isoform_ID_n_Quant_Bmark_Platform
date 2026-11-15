#!/usr/bin/env Rscript
# Supplementary figure: SQANTI-like categories of the reads (a, b) and of the LRAA isoforms (c, d)
# for the PBMC and FTD (ALS LP00714) single-cell datasets. These were Figures 5a,b and 6a,b.
#
# Reads: each eval tree's __SQLIKE_read_analysis/*read_cats_sqlike.tsv (the aligned reads against
# the reference annotation; independent of the LRAA version).
# Isoforms: each eval tree's __SQLIKE_gtf_analysis/*_sqlike_spC.tsv (the current LRAA version's
# four runs, isoforms with the same splice pattern counted once).
#
# Writes SuppFig_categories.pdf at the manuscript text width, fonts embedded.
# Run from this directory:  Rscript SuppFig_categories.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
})
source("../../LRAA_paper-PBMCs_study/__Paper_Figures/paper_fig_style.R")

EVAL = c(PBMC = "../../LRAA_paper-PBMCs_study", FTD = "../../LRAA_paper-FTD_study")
READ_CATS = c(PBMC = "__SQLIKE_read_analysis/PBMC-read_cats_sqlike.tsv",
              FTD = "__SQLIKE_read_analysis/ALS-LP00714-read_cats_sqlike.tsv")
ISO_CATS = c(PBMC = "__SQLIKE_gtf_analysis/PBMCs-LRAA_sqlike_spC.tsv",
             FTD = "__SQLIKE_gtf_analysis/LP00714-ALS-LRAA_sqlike_spC.tsv")
TITLE = c(PBMC = "PBMC", FTD = "FTD frontal cortex")

CATS = c("FSM", "ISM", "NIC", "NNIC", "genic", "intronic", "antisense", "intergenic",
         "se_FM", "se_IM", "se_genic", "se_exonic", "se_intronic", "se_antisense", "se_intergenic")
RUNS = c("denovo-basic-spC" = "de novo, pseudobulk",
         "denovo-scg-spC" = "de novo, cluster-guided",
         "refGuided-basic-spC" = "ref-guided, pseudobulk",
         "refGuided-scg-spC" = "ref-guided, cluster-guided")
RUN_COLORS = c("#f8766d", "#7cae00", "#00bfc4", "#c77cff")   # as in the draft figures

WIDTH_IN = 6.5
HEIGHT_IN = 7

read_plot = function(ds) {
    d = read_tsv(file.path(EVAL[[ds]], READ_CATS[[ds]]), show_col_types = FALSE) %>%
        mutate(Category = factor(Category, levels = CATS))
    total = sum(d$Count)
    ggplot(d, aes(Category, Count / 1e6, fill = Category)) +
        geom_col(width = 0.75, show.legend = FALSE) +
        # share of reads on each bar (user, 2026-11-14), one decimal below 10%
        geom_text(aes(label = ifelse(Count / total < 0.1, sprintf("%.1f%%", 100 * Count / total),
                                     sprintf("%.0f%%", 100 * Count / total))),
                  angle = 90, hjust = -0.15, size = PAPER$small_size / ggplot2::.pt,
                  family = PAPER$font) +
        scale_y_continuous(name = "Reads (millions)", expand = expansion(mult = c(0, 0.18)),
                           sec.axis = sec_axis(~ . * 1e6 / total, name = "Share of reads",
                                               labels = scales::percent_format(accuracy = 1))) +
        labs(x = NULL, title = sprintf("%s reads (%.1fM)", TITLE[[ds]], total / 1e6)) +
        theme_paper() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

iso_plot = function(ds) {
    d = read_tsv(file.path(EVAL[[ds]], ISO_CATS[[ds]]), show_col_types = FALSE) %>%
        filter(type %in% names(RUNS)) %>%
        mutate(Category = factor(Category, levels = CATS),
               run = factor(RUNS[type], levels = RUNS))
    ggplot(d, aes(run, Count / 1e3, fill = run)) +
        geom_col(width = 0.85) +
        facet_grid(. ~ Category, switch = "x") +
        scale_fill_manual(values = setNames(RUN_COLORS, RUNS), name = NULL) +
        guides(fill = guide_legend(nrow = 1)) +
        scale_y_continuous(name = "Isoforms (thousands)", expand = expansion(mult = c(0, 0.05))) +
        labs(x = NULL, title = sprintf("%s isoforms (splice patterns counted once)", TITLE[[ds]])) +
        theme_paper() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              strip.placement = "outside",
              strip.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = PAPER$small_size),
              panel.spacing.x = unit(0.8, "mm"),
              panel.border = element_rect(fill = NA, linewidth = PAPER$line),
              legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
}

fig = (read_plot("PBMC") | read_plot("FTD")) / iso_plot("PBMC") / iso_plot("FTD") +
    plot_layout(heights = c(1, 1, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 11, face = "bold", family = PAPER$font),
          legend.position = "bottom")

ggsave("SuppFig_categories.pdf", fig, width = WIDTH_IN, height = HEIGHT_IN, device = grDevices::pdf, useDingbats = FALSE)
embed_fonts("SuppFig_categories.pdf")   # paper_fig_style.R: embeds Helvetica (Nimbus Sans) without changing the layout
cat("wrote SuppFig_categories.pdf\n")
