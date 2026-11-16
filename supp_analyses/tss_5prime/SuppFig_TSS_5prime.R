#!/usr/bin/env Rscript
# Supplementary figure: 5' completeness of the reads and the untemplated-G (cap) signature, FTD
# frontal cortex nuclei (ALS LP00714) vs PBMCs.
#
# (a) aligned read length by cluster (genome-wide read sample);
# (b) 5' completeness K(d) by cell-type class: the share of a gene's polyA-anchored reads
#     extending at least d nt from the 3' end, averaged over single-TSS genes;
# (c) control genes: share of reads starting at the observed TSS vs ending at the observed 3' end;
# (d) share of reads with an untemplated-G signature, for reads starting at the TSS vs inside the
#     transcript body (single-TSS genes);
# (e) the capped share of reads at each read 5' pile-up (TSS peak called from the reads alone), split by
#     whether a CAGE peak lies within 25 nt (since 2026-11-11; earlier versions used TSS switch sites).
#
# Inputs: ALS_LP00714/ALS_LRAA_eval/__TSS_Truncation_Check/{ALS,PBMC}.{lengths,reads}.tsv.gz (a, c),
# __TSS_5prime_Completeness/{ALS,PBMC}.{genes,extent,capsig} (b, d), and for (e)
# __TSS_5prime_Completeness/cap_calibration/out/{FTD,PBMC}.peaks.tsv.gz: every read 5' pile-up LRAA's
# pre-pass calls, with informative and capped read counts by LRAA's trusted-5'-end rules (G run above the
# sc Kinnex baseline of 0) and CAGE overlap within 25 nt.
# The computations follow those dirs' tss_truncation_check.Rmd and tss_5prime_completeness.Rmd.
# Signature rates use reads whose genomic base upstream of the 5' end is not G, so that a G in the
# soft clip cannot be templated.
#
# Writes SuppFig_TSS_5prime.pdf at the manuscript text width, fonts embedded, and
# SuppFig_TSS_5prime.values.tsv (the numbers the text cites).
# Run from this directory:  Rscript SuppFig_TSS_5prime.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(patchwork)
})
source("../../LRAA_paper-PBMCs_study/__Paper_Figures/paper_fig_style.R")

TRUNC = "truncation_check"
COMPL = "completeness"
PEAKS = c(ALS = file.path(COMPL, "cap_calibration/out/FTD.peaks.tsv.gz"),
          PBMC = file.path(COMPL, "cap_calibration/out/PBMC.peaks.tsv.gz"))
DS = c("ALS", "PBMC")
DS_LABEL = c(ALS = "FTD", PBMC = "PBMC")

class_map = list(
    ALS = c("0" = "Oligodendrocyte", "3" = "Oligodendrocyte", "1" = "Excitatory neuron", "8" = "Excitatory neuron",
            "4" = "Inhibitory neuron", "6" = "Inhibitory neuron", "2" = "Astrocyte", "7" = "Astrocyte",
            "5" = "Microglia", "9" = "OPC", "10" = "Vascular / meningeal", "11" = "Vascular / meningeal",
            "12" = "Vascular / meningeal", "13" = "Vascular / meningeal"),
    PBMC = c("1" = "Monocyte / DC", "11" = "Monocyte / DC", "13" = "Monocyte / DC", "6" = "B", "9" = "B",
             "3" = "CD4 T", "7" = "CD4 T", "0" = "CD4 T", "5" = "CD4 T", "12" = "CD4 T",
             "2" = "CD8 T / NK", "4" = "CD8 T / NK", "8" = "CD8 T / NK", "10" = "CD8 T / NK"))
class_order = c("Excitatory neuron", "Inhibitory neuron", "Oligodendrocyte", "OPC", "Astrocyte", "Microglia",
                "Vascular / meningeal", "Monocyte / DC", "B", "CD4 T", "CD8 T / NK")
class_colors = setNames(c("#b2182b", "#ef8a62", "#2166ac", "#67a9cf", "#1b7837", "#7fbf7b", "#8c8b86",
                          "#762a83", "#e08214", "#35978f", "#c51b7d"), class_order)
cls = function(ds, cluster) unname(class_map[[ds]][sub("Cluster_", "", cluster)])
neuron = c("Excitatory neuron", "Inhibitory neuron")

WIDTH_IN = 6.5
HEIGHT_IN = 8.6
# space between facets, so the last tick label of one facet does not run into the first of the next
FACET_GAP = theme(panel.spacing.x = unit(6, "mm"))
# a box around each facet where the panels would otherwise read as one
FACET_BOX = theme(panel.border = element_rect(fill = NA, color = "black", linewidth = PAPER$line),
                  axis.line = element_blank(), panel.spacing.x = unit(3, "mm"))
values = list()
val = function(...) values[[length(values) + 1]] <<- tibble(...)

# ---- (a) read length by cluster
lengths = rbindlist(lapply(DS, function(ds) fread(file.path(TRUNC, paste0(ds, ".lengths.tsv.gz")))[, dataset := ds]))
lengths[, aligned_len := as.numeric(aligned_len)]
lengths[, class := cls(dataset, cluster), by = dataset]
lengths[, cl := sub("Cluster_", "C", cluster)]
len_stats = lengths[, .(median = median(aligned_len), q25 = quantile(aligned_len, 0.25),
                        q75 = quantile(aligned_len, 0.75)), by = .(dataset, cl, class)]
len_stats[, `:=`(dataset = factor(DS_LABEL[dataset], levels = DS_LABEL), class = factor(class, levels = class_order))]
val(lengths[, .(median = median(aligned_len)), by = .(dataset, neuron = class %in% neuron)] %>%
        transmute(dataset, group = if_else(neuron, "neurons", if_else(dataset == "PBMC", "all", "other")), measure = "median_read_length", value = median))
val(lengths[, .(median = median(aligned_len)), by = .(dataset, cl)] %>%
        transmute(dataset, group = cl, measure = "median_read_length", value = median))

p_a = ggplot(len_stats, aes(reorder(paste(dataset, cl), median), median, color = class)) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.6) +
    geom_point(size = 1) +
    facet_grid(~ dataset, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = function(x) sub(".* ", "", x)) +
    scale_color_manual(values = class_colors, drop = TRUE, guide = "none") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = NULL, y = "Read length (nt)") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# ---- (b) 5' completeness curves by cell-type class
genes = rbindlist(lapply(DS, function(ds) fread(file.path(COMPL, paste0(ds, ".genes.tsv")))[, dataset := ds]), fill = TRUE)
use = genes[status == "ok" & cage %in% c("single", "none")]
ext = rbindlist(lapply(DS, function(ds) fread(file.path(COMPL, paste0(ds, ".extent.tsv.gz")))[, dataset := ds]))
ext = merge(ext, use[, .(dataset, gene_key, L)], by = c("dataset", "gene_key"))
ext[, class := cls(dataset, cluster), by = dataset]
hist = ext[extent_hist != "", {
    p = strsplit(strsplit(extent_hist, ";")[[1]], ":")
    .(bin = as.integer(sapply(p, `[`, 1)), n = as.integer(sapply(p, `[`, 2)))
}, by = .(dataset, gene_key, cluster)]
hist = merge(hist, ext[, .(dataset, gene_key, cluster, class)], by = c("dataset", "gene_key", "cluster"))
grid = seq(0, 4000, by = 25)
tot = ext[, .(n = sum(n_anchored), r = sum(n_reach), L = L[1]), by = .(dataset, gene_key, class)][n >= 20]
h = hist[, .(n = sum(n)), by = .(dataset, gene_key, class, bin)]
hc = h[order(-bin), .(bin, cum = cumsum(n)), by = .(dataset, gene_key, class)][order(bin)]
pts = tot[, .(d = grid[grid <= L]), by = .(dataset, gene_key, class)]
pts = merge(pts, tot, by = c("dataset", "gene_key", "class"))
pts[, bin := d]
pts[, cum := hc[pts, on = c("dataset", "gene_key", "class", "bin"), roll = -Inf, x.cum]]
pts[is.na(cum), cum := 0]
pts[, S := (r + cum) / n]
K = pts[, .(K = mean(S), genes = .N), by = .(dataset, class, d)][genes >= 20]
val(K[d == 1000] %>% transmute(dataset, group = class, measure = "K_1000nt", value = K))
K[, `:=`(dataset = factor(DS_LABEL[dataset], levels = DS_LABEL), class = factor(class, levels = class_order))]

p_b = ggplot(K, aes(d, K, color = class)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ dataset) +
    scale_color_manual(values = class_colors, drop = TRUE, name = "Cell type") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Distance from the 3' end (nt)", y = "Share of reads") +
    theme_paper() + FACET_GAP + theme(legend.key.height = unit(3, "mm")) +
    guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 0.8)))

# ---- (c) control genes: reads reaching the TSS and the 3' end
reads = rbindlist(lapply(DS, function(ds) fread(file.path(TRUNC, paste0(ds, ".reads.tsv.gz")))[, dataset := ds]))
reads[, class := cls(dataset, cluster), by = dataset]
mode10 = function(x) { x = x[!is.na(x)]; b = floor(x / 10) * 10; as.numeric(names(which.max(table(b)))) + 5 }
reads[, `:=`(obs_tss = mode10(five_prime_tx_pos[five_prime_tx_pos <= first_exon_len[1] + 30]),
             obs_end = mode10(three_prime_tx_pos)), by = .(dataset, gene)]
reads[, `:=`(reach_tss = abs(five_prime_tx_pos - obs_tss) <= 50, reach_end = abs(three_prime_tx_pos - obs_end) <= 50)]
reach = reads[, .(reads = .N, reach_TSS = mean(reach_tss, na.rm = TRUE), reach_end = mean(reach_end, na.rm = TRUE)),
              by = .(dataset, gene, class)][reads >= 30]
val(reach %>% transmute(dataset, group = paste(gene, class, sep = ":"), measure = "reach_TSS", value = reach_TSS))
reach[, `:=`(dataset = factor(DS_LABEL[dataset], levels = DS_LABEL), class = factor(class, levels = class_order))]

# one row per dataset (user, 2026-11-14: FTD and PBMC were hard to tell apart as shapes in one panel)
p_c = ggplot(reach, aes(reach_TSS, reach_end, color = class)) +
    geom_point(size = 1.2) +
    facet_grid(dataset ~ gene, switch = "y") +
    scale_color_manual(values = class_colors, drop = TRUE, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    labs(x = "Share of reads starting at the observed TSS", y = "Share ending at\nthe observed 3' end") +
    theme_paper() + FACET_GAP +
    theme(strip.placement = "outside", panel.spacing.y = unit(2, "mm"),
          strip.text.y.left = element_text(angle = 0, face = "bold"))

# ---- (d) untemplated-G signature: reads starting at the TSS vs inside the transcript body
cap = rbindlist(lapply(DS, function(ds) fread(file.path(COMPL, paste0(ds, ".capsig.tsv.gz")))[, dataset := ds]))
cap = cap[use[, .(dataset, gene_key)], on = .(dataset, gene_key), nomatch = 0]
cap[, class := cls(dataset, cluster), by = dataset]
grate = cap[upstream_G == 0 & five_class %in% c("TSS", "body"),
            .(G_rate = sum(reads[clip == "G"]) / sum(reads), reads = sum(reads)), by = .(dataset, class, five_class)]
val(cap[upstream_G == 0 & five_class %in% c("TSS", "body"),
        .(value = sum(reads[clip == "G"]) / sum(reads)), by = .(dataset, five_class)] %>%
        transmute(dataset, group = five_class, measure = "G_signature_rate", value))
grate[, `:=`(dataset = factor(DS_LABEL[dataset], levels = DS_LABEL),
             five_class = factor(five_class, levels = c("TSS", "body"), labels = c("at the TSS", "in the transcript body")),
             class = factor(class, levels = class_order))]

p_d = ggplot(grate, aes(class, G_rate, fill = five_class)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(~ dataset, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("at the TSS" = "#2a78d6", "in the transcript body" = "#8c8b86"), name = "Read starts") +
    guides(fill = guide_legend(nrow = 2)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Reads with an untemplated G") +
    theme_paper() + FACET_BOX +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# ---- (e) capped share at read 5' pile-ups (LRAA's TSS peak calling from the reads, before any TSS rule),
# split by whether a FANTOM5 CAGE peak lies within 25 nt; pile-ups with >= 20 informative reads
pk = rbindlist(lapply(DS, function(ds) fread(PEAKS[[ds]])[, dataset := ds]))[n_informative >= 20]
pk[, group := ifelse(as.logical(in_CAGE), "at CAGE peak", "no CAGE peak")]
val(pk[, .(value = median(capped_frac)), by = .(dataset, group)] %>% transmute(dataset, group, measure = "pileup_capped_frac_median", value))
val(pk[, .(value = .N), by = .(dataset, group)] %>% transmute(dataset, group, measure = "pileups", value))
val(pk[, .(value = mean(capped_frac < 0.1)), by = .(dataset, group)] %>% transmute(dataset, group, measure = "pileup_share_below_0.1", value))
pk[, dataset := factor(DS_LABEL[dataset], levels = DS_LABEL)]
p_e = ggplot(pk, aes(capped_frac, after_stat(density), fill = group)) +
    geom_histogram(binwidth = 0.05, boundary = 0, position = "identity", alpha = 0.6, color = NA) +
    scale_fill_manual(values = c("at CAGE peak" = "#1b7837", "no CAGE peak" = "#762a83"), name = NULL) +
    facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    labs(x = "Capped share of reads at a read 5\' pile-up", y = "Pile-ups (density)") +
    theme_paper() + FACET_GAP + theme(legend.position = "top")

# each panel keeps its own legend (collected legends sat far from (d) and (e)); a left margin so no
# y-axis title is cut at the page edge
fig = (p_a / p_b / p_c / (p_d | p_e) +
       plot_layout(heights = c(0.9, 1, 1.1, 1.25))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 11, face = "bold", family = PAPER$font),
          plot.margin = margin(4, 4, 4, 8))

ggsave("SuppFig_TSS_5prime.pdf", fig, width = WIDTH_IN, height = HEIGHT_IN, device = grDevices::pdf, useDingbats = FALSE)
embed_fonts("SuppFig_TSS_5prime.pdf")

values = bind_rows(values) %>%
    mutate(dataset = DS_LABEL[dataset])
write_tsv(values, "SuppFig_TSS_5prime.values.tsv")
cat("wrote SuppFig_TSS_5prime.pdf and SuppFig_TSS_5prime.values.tsv\n")
