#!/usr/bin/env Rscript
# Supplementary figure: agreement of LRAA's TSS and PolyA sites with FANTOM5 CAGE peaks and
# PolyASite 2.0 clusters, PBMC and FTD (ALS LP00714) single-cell datasets.
#
# Sites: the de novo runs' splice-pattern-collapsed TSS / PolyA beds (pseudobulk and
# cluster-guided), called from the reads alone. Each eval tree's __TSS_PolyA_atlas_agreement
# (make) writes the per-site distance tables read here (<TSS|PolyA>.<basic|scg>.dist.tsv:
# distance to the nearest same-strand atlas feature, and to its representative point; the
# same sites shuffled within the LRAA transcript spans give the agreement expected by chance).
# The panels follow that dir's TSS_PolyA_atlas_agreement.Rmd.
#
# (a) share of sites inside an atlas feature, with the shuffled control;
# (b) that share by read support at the site;
# (c) signed offset from the atlas representative point, sites within 100 nt.
#
# Writes SuppFig_site_atlas.pdf at the manuscript text width, fonts embedded, and
# SuppFig_site_atlas.values.tsv (the numbers the text cites).
# Run from this directory:  Rscript SuppFig_site_atlas.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
})
source("../../LRAA_paper-PBMCs_study/__Paper_Figures/paper_fig_style.R")

EVAL = c(PBMC = "../../LRAA_paper-PBMCs_study", FTD = "../../LRAA_paper-FTD_study")
SETS = c(basic = "pseudobulk", scg = "cluster-guided")
SITE_LABEL = c(TSS = "TSS vs FANTOM5 CAGE", PolyA = "PolyA vs PolyASite 2.0")
DS_COLORS = c(PBMC = "#2a78d6", FTD = "#eb6834")

WIDTH_IN = 6.5
HEIGHT_IN = 6.5

read_dist = function(ds, kind, set) {
    read_tsv(file.path(EVAL[[ds]], "__TSS_PolyA_atlas_agreement", sprintf("%s.%s.dist.tsv", kind, set)),
             show_col_types = FALSE,
             col_types = cols_only(support = "d", match_name = "c", dist_interval = "d",
                                   point_name = "c", dist_point = "d", data = "c")) %>%
        # no same-strand atlas feature on the chromosome: bedtools reports "." and -1
        mutate(dist_interval = if_else(match_name == ".", NA_real_, dist_interval),
               dist_point = if_else(point_name == ".", NA_real_, dist_point),
               hit = !is.na(dist_interval) & dist_interval == 0,
               dataset = ds, site_type = kind, set = SETS[[set]])
}

dist = expand_grid(ds = names(EVAL), kind = names(SITE_LABEL), set = names(SETS)) %>%
    pmap_dfr(read_dist) %>%
    mutate(dataset = factor(dataset, levels = names(EVAL)),
           site_type = factor(site_type, levels = names(SITE_LABEL), labels = SITE_LABEL),
           set = factor(set, levels = SETS))

# ---- (a) share inside an atlas feature
agree = dist %>%
    group_by(dataset, site_type, set, data) %>%
    summarize(sites = n(), frac = mean(hit), .groups = "drop")
obs = filter(agree, data == "observed")
chance = filter(agree, data == "shuffled")

p_a = ggplot(obs, aes(set, frac, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = scales::percent(frac, accuracy = 1)), position = position_dodge(width = 0.75),
              vjust = -0.4, size = pt_to_size(PAPER$small_size), family = PAPER$font) +
    geom_point(data = chance, aes(set, frac, group = dataset), position = position_dodge(width = 0.75),
               shape = 4, size = 1, stroke = 0.4, show.legend = FALSE) +
    facet_wrap(~ site_type) +
    scale_fill_manual(values = DS_COLORS, name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "de novo LRAA sites", y = "Sites inside an atlas feature") +
    theme_paper()

# ---- (b) by read support
support_bins = function(x) cut(x, breaks = c(0, 9, 99, 999, Inf), labels = c("<10", "10-99", "100-999", ">=1,000"))
by_support = dist %>%
    filter(data == "observed") %>%
    mutate(support_bin = support_bins(support)) %>%
    group_by(dataset, site_type, set, support_bin) %>%
    summarize(sites = n(), frac = mean(hit), .groups = "drop")

p_b = ggplot(by_support, aes(support_bin, frac, color = dataset, linetype = set,
                             group = interaction(dataset, set))) +
    geom_line(linewidth = 0.4) + geom_point(size = 0.8) +
    facet_wrap(~ site_type) +
    scale_color_manual(values = DS_COLORS, name = NULL) +
    scale_linetype_manual(values = c(pseudobulk = "solid", `cluster-guided` = "dashed"), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Reads supporting the site", y = "Sites inside an atlas feature") +
    theme_paper()

# ---- (c) signed offset from the atlas representative point: share of the sites within 100 nt
# per 5-nt bin (square-root scale, so the shoulders beside the 0-nt peak show)
offsets = dist %>%
    filter(data == "observed", abs(dist_point) <= 100) %>%
    mutate(bin = round(dist_point / 5) * 5) %>%
    count(dataset, site_type, set, bin) %>%
    group_by(dataset, site_type, set) %>%
    mutate(share = n / sum(n)) %>%
    ungroup()

p_c = ggplot(offsets, aes(bin, share, color = set)) +
    geom_step(direction = "mid", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = PAPER$line, color = "grey50") +
    facet_grid(dataset ~ site_type) +
    scale_color_manual(values = c(pseudobulk = "#8c8b86", `cluster-guided` = "#1baf7a"), name = NULL) +
    scale_y_sqrt(labels = scales::percent, breaks = c(0.01, 0.1, 0.25, 0.5, 0.75)) +
    labs(x = "LRAA site - atlas representative point (nt; <0 upstream)", y = "Sites within 100 nt (sqrt scale)") +
    theme_paper() +
    theme(strip.text.y = element_text(angle = 0))

fig = (p_a / p_b / p_c) + plot_layout(heights = c(1, 1, 1.4)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 11, face = "bold", family = PAPER$font))

ggsave("SuppFig_site_atlas.pdf", fig, width = WIDTH_IN, height = HEIGHT_IN, device = grDevices::pdf, useDingbats = FALSE)
embed_fonts("SuppFig_site_atlas.pdf")

values = bind_rows(
    obs %>% transmute(dataset, site_type, set, measure = "frac_in_atlas", value = frac, n = sites),
    chance %>% transmute(dataset, site_type, set, measure = "frac_in_atlas_shuffled", value = frac, n = sites),
    by_support %>% transmute(dataset, site_type, set, measure = paste0("frac_in_atlas_support_", support_bin),
                             value = frac, n = sites)) %>%
    mutate(site_type = sub(" vs .*", "", site_type))
write_tsv(values, "SuppFig_site_atlas.values.tsv")
cat("wrote SuppFig_site_atlas.pdf and SuppFig_site_atlas.values.tsv\n")
