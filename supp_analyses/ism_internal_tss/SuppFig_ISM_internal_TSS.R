#!/usr/bin/env Rscript
# Supplementary figure: internal TSSs behind ISM-shaped isoforms, and the shorter proteins they imply.
#
# (a) strict ISM TSSs with a capped read pile-up, by further evidence tier, per library;
# (b) CAGE peaks within 10 nt of capped and uncapped read pile-ups, against random positions in the same exon;
# (c) each tool's own strict ISM isoforms (bulk libraries): truncation-like share of judged starts, supported TSSs;
# (d-g) COMT (known internal promoter), the FTD example, CORO1B and RECQL4: L and S models with both CDSs (blue; GENCODE,
# or pyTD's ORF for L and its in-frame part from S's ATG), a sample of the reads compatible with L or S (2026-11-07: the
# read's intron chain is a run of the model's introns and its ends lie in the matching exons), and below, aligned to the
# same compressed x-axis, the 5' ends of those reads over the locus; reads coloured by where they start (L's TSS, the S
# site, elsewhere), both TSSs marked. The title gives the share of S-site reads with the cap signature.
# The FTD example is select_ftd_example.py's rank 1 (rule in out/examples/ftd_example_candidates.tsv; user, 2026-11-04,
# replaces MBP; rule 4 added 2026-11-08, EPB41L3 replaces BTBD8); it is example "ftd" in examples.tsv. COMT, CORO1B and RECQL4 are pinned by gene.
#
# (c) replaced 2026-11-10 (user): the pooled-target recall (tool_recall.tsv) was mostly IsoSeq's own fragments; each tool is
# now judged on its own ISM isoforms (tool_own_ism.py).
#
# Inputs: __Benchmark_extras/ISM_internal_TSS/out/ (tier_counts_dataset.tsv, atlas_background.tsv,
# tool_own_ism.pooled.tsv, examples/). COMT, CORO1B, RECQL4 chosen by gene (user, 2026-10-31) from that study's selection.
# Isoforms carry -L / -S labels; nothing here uses LRAA model IDs.
#
# Writes SuppFig_ISM_internal_TSS.pdf at the manuscript text width, fonts embedded, and
# SuppFig_ISM_internal_TSS.values.tsv (the numbers the text cites).
# Run from this directory:  Rscript SuppFig_ISM_internal_TSS.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
})
source("../../LRAA_paper-PBMCs_study/__Paper_Figures/paper_fig_style.R")

IN = "out"
DS = c("PBMC", "ALS", "BT474", "K562", "HG002", "UHRR")
DS_LABEL = c(PBMC = "PBMC", ALS = "FTD", BT474 = "BT474", K562 = "K562", HG002 = "HG002", UHRR = "UHRR")
GENES = c("COMT", "CORO1B", "RECQL4")   # + the FTD example (example id "ftd")
TOOLS = c("LRAA", "IsoQuant", "StringTie", "Bambu", "ESPRESSO", "IsoSeq", "Mandalorion")
WIDTH_IN = 6.5
HEIGHT_IN = 8.6
ds_f = function(x) factor(DS_LABEL[as.character(x)], DS_LABEL)

# ---- (a) tiers
tiers = c("all" = "All", "pile-up+cap" = "Pile-up + cap", "core: pile-up+cap+atlas" = "+ CAGE or refTSS",
          "core, not canonical TSS" = "+ not canonical TSS", "core, protein altered" = "+ protein altered",
          "core, N-truncated" = "+ N-truncated")
td = read_tsv(file.path(IN, "tier_counts_dataset.tsv"), show_col_types = FALSE) %>%
    filter(class == "strict_ISM", tier %in% names(tiers)) %>%
    mutate(label = factor(tiers[tier], rev(tiers)), ds = ds_f(dataset))
# the funnel starts at read-supported starts (user, 2026-11-10): "All" counts every strict ISM model, most of
# them without a TSS call (TSS=False); it stays in the values file for the text
# linear axis (2026-11-10): the funnel shows how much each filter removes
pa = ggplot(td %>% filter(tier != "all"), aes(peaks, label)) +
    geom_col(fill = ISOFORM_COLORS[1], width = 0.7) +
    geom_text(aes(label = scales::comma(peaks)), hjust = -0.1, size = 1.9) +
    facet_wrap(~ds, nrow = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25)), breaks = c(0, 100)) +
    labs(x = "strict ISM TSSs", y = NULL, tag = "a") + theme_paper()

# ---- (b) CAGE, capped vs uncapped pile-ups
ab = read_tsv(file.path(IN, "atlas_background.tsv"), show_col_types = FALSE) %>%
    filter(atlas == "CAGE", window == 10, stratum %in% c("pile-up, capped", "pile-up, measurable, not capped")) %>%
    mutate(ds = ds_f(dataset), group = if_else(stratum == "pile-up, capped", "Capped pile-up", "Uncapped pile-up"))
bars = bind_rows(ab %>% transmute(ds, group, kind = "Site", value = obs), ab %>% transmute(ds, group, kind = "Same-exon random", value = bg))
pb = ggplot(bars, aes(ds, 100 * value, fill = kind)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    facet_wrap(~group) +
    scale_fill_manual(values = c(Site = ISOFORM_COLORS[2], `Same-exon random` = OTHER_ISOFORM_COLOR), name = NULL) +
    labs(x = NULL, y = "% within 10 nt\nof a CAGE peak", tag = "b") + theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top",
          panel.spacing.x = unit(6, "mm"))   # capped and uncapped facets kept visibly apart

# ---- (c) each tool's OWN strict ISM isoforms (S and L from the tool, or L from GENCODE), their starts judged by reads and
# atlases (tool_own_ism.py; user, 2026-11-10). x: truncation-like share of the starts the reads can judge; y: supported
# internal TSSs; dot area: the tool's ISM isoforms. Four bulk libraries pooled.
ti = read_tsv(file.path(IN, "tool_own_ism.pooled.tsv"), show_col_types = FALSE) %>% filter(tool %in% TOOLS)
pc = ggplot(ti, aes(pct_truncation_of_judged, supported)) +
    geom_point(aes(size = ism_isoforms, colour = tool == "LRAA"), alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = tool), size = 2, min.segment.length = 0, seed = 1, box.padding = 0.4, point.padding = 0.5, max.overlaps = Inf) +
    scale_colour_manual(values = c(`TRUE` = ISOFORM_COLORS[2], `FALSE` = OTHER_ISOFORM_COLOR), guide = "none") +
    scale_size(trans = "log10", range = c(0.8, 6), breaks = c(1e3, 1e4, 1e5), labels = c("1,000", "10,000", "100,000"), name = "ISM isoforms") +
    scale_x_continuous(limits = c(0, 100)) + scale_y_log10(limits = c(5, 600)) +
    labs(x = "Truncation-like ISM starts\n(% of judged; lower is better)", y = "Supported ISM TSSs\n(higher is better)", tag = "c") +
    theme_paper() + theme(legend.position = "top", legend.title = element_text(size = 5.5), legend.text = element_text(size = 5),
                          legend.key.size = unit(2.5, "mm"), legend.margin = margin(0, 0, 0, 0))

# ---- (d-g) examples
ex_all = read_tsv(file.path(IN, "examples/examples.tsv"), show_col_types = FALSE, guess_max = 100)
ex = bind_rows(ex_all %>% filter(example != "ftd", gene %in% GENES), ex_all %>% filter(example == "ftd"))
FTD_GENE = ex$gene[ex$example == "ftd"]
fc = read_tsv(file.path(IN, "examples/ftd_example_candidates.tsv"), comment = "#", show_col_types = FALSE) %>% filter(passes) %>% slice(1)
fit_stats = list()
mod = read_tsv(file.path(IN, "examples/models.tsv"), show_col_types = FALSE)
rd = read_tsv(file.path(IN, "examples/reads.tsv.gz"), show_col_types = FALSE, col_types = cols(capped = "c"))
st = read_tsv(file.path(IN, "examples/starts.tsv"), show_col_types = FALSE)
atl = read_tsv(file.path(IN, "examples/atlas.tsv"), show_col_types = FALSE)
READ_COLS = c("Starts at L's TSS" = "#2166ac", "Starts at the S site" = "#1a9850", "Starts elsewhere" = "#9e9e9e")
L_TSS_WINDOW = 50   # nt: a read starts at L's TSS
INTRON_MIN = 30     # nt: a gap between aligned blocks this long is an intron (shorter: an alignment gap)
SPLICE_TOL = 5      # nt: read and model introns match within this
END_TOL = 50        # nt: a read may pass a model's terminal exon end by this much
read_exons = function(blocks) {
    b = do.call(rbind, lapply(strsplit(strsplit(blocks, ";")[[1]], "-"), as.numeric))
    b = b[order(b[, 1]), , drop = FALSE]; out = b[1, , drop = FALSE]
    for (i in seq_len(nrow(b))[-1]) {
        if (b[i, 1] - out[nrow(out), 2] < INTRON_MIN) out[nrow(out), 2] = max(out[nrow(out), 2], b[i, 2]) else out = rbind(out, b[i, ])
    }
    out
}
# the read fits the model: its introns are a run of the model's introns, and its first and last exons lie within the
# matching model exons (a terminal model exon may be passed by END_TOL); an unspliced read lies within one model exon
compatible = function(rx, mx) {
    ri = cbind(rx[-nrow(rx), 2], rx[-1, 1]); mi = cbind(mx[-nrow(mx), 2], mx[-1, 1]); n = nrow(mx)
    tol_lo = function(k) if (k == 1) END_TOL else 0
    tol_hi = function(k) if (k == n) END_TOL else 0
    if (nrow(ri) == 0) return(any(sapply(seq_len(n), function(k) rx[1, 1] >= mx[k, 1] - tol_lo(k) && rx[1, 2] <= mx[k, 2] + tol_hi(k))))
    if (nrow(ri) > nrow(mi)) return(FALSE)
    for (k in seq_len(nrow(mi) - nrow(ri) + 1)) {
        idx = k:(k + nrow(ri) - 1)
        if (all(abs(ri[, 1] - mi[idx, 1]) <= SPLICE_TOL & abs(ri[, 2] - mi[idx, 2]) <= SPLICE_TOL) &&
            rx[1, 1] >= mx[k, 1] - tol_lo(k) && rx[nrow(rx), 2] <= mx[k + nrow(ri), 2] + tol_hi(k + nrow(ri))) return(TRUE)
    }
    FALSE
}
compressor = function(segs, sgn, gap = 120) {
    s = segs %>% arrange(start) %>% mutate(g = cumsum(start > lag(cummax(end), default = first(start)))) %>%
        group_by(g) %>% summarise(a = min(start), b = max(end), .groups = "drop")
    ga = c(); pa = c(); pos = 0
    for (i in seq_len(nrow(s))) {
        if (i > 1) pos = pos + min(gap, s$a[i] - s$b[i - 1])
        ga = c(ga, s$a[i], s$b[i]); pa = c(pa, pos, pos + s$b[i] - s$a[i]); pos = pos + s$b[i] - s$a[i]
    }
    f = approxfun(ga, pa, rule = 2)
    if (sgn > 0) f else function(v) -f(v)
}
gene_panel = function(id, tag) {
    e = ex %>% filter(example == id); g = e$gene
    sgn = if_else(e$strand == "+", 1, -1)
    m = mod %>% filter(example == id, track %in% c("L", "S")); xm0 = m %>% filter(kind %in% c("exon", "cds"))
    fx = compressor(xm0 %>% distinct(start, end), sgn)
    mx = lapply(c(L = "L", S = "S"), function(k) { x = m %>% filter(track == k, kind == "exon") %>% arrange(start); cbind(x$start, x$end) })
    a = rd %>% filter(example == id)
    a$fits = sapply(a$blocks, function(b) { rx = read_exons(b); compatible(rx, mx$L) || compatible(rx, mx$S) })
    n_reads = nrow(a)
    a = a %>% filter(fits) %>% mutate(class = factor(case_when(
        abs(start5 - e$l_tss) <= L_TSS_WINDOW ~ "Starts at L's TSS",
        at_site ~ "Starts at the S site", TRUE ~ "Starts elsewhere"), names(READ_COLS)))
    # full-chain (FSM) reads: the read's intron chain is the whole chain of L (fsm_L) or of S (fsm_S); the sample always
    # holds at least one of each that exists (the one reaching furthest 5')
    chain_of = function(b) { rx = read_exons(b); if (nrow(rx) < 2) "" else paste(rx[-nrow(rx), 2], rx[-1, 1], collapse = ",") }
    mchain = sapply(mx, function(x) paste(x[-nrow(x), 2], x[-1, 1], collapse = ","))
    rch = sapply(a$blocks, chain_of); a$fsm_L = rch == mchain[["L"]]; a$fsm_S = rch == mchain[["S"]]
    fit_stats[[id]] <<- a %>% count(class, .drop = FALSE) %>% mutate(gene = g, reads_all = n_reads, fsm_L = sum(a$fsm_L), fsm_S = sum(a$fsm_S))
    thin = function(d, n) if (nrow(d) > n) d %>% arrange(start5) %>% slice(unique(round(seq(1, n(), length.out = n)))) else d
    r = bind_rows(thin(a %>% filter(class == "Starts at L's TSS"), 10), thin(a %>% filter(class == "Starts at the S site"), 14),
                  thin(a %>% filter(class == "Starts elsewhere"), 10))
    for (k in c("fsm_L", "fsm_S")) if (!any(r[[k]]) && any(a[[k]]))
        r = bind_rows(r, a %>% filter(.data[[k]]) %>% slice_min(sgn * start5, n = 1, with_ties = FALSE))
    r = r %>% arrange(class, sgn * start5)
    trk = tibble(track = c("CAGE", paste0(g, "-L"), paste0(g, "-S")), key = c("CAGE", "L", "S"), y = c(-0.6, -3, -6))
    bl = r %>% mutate(row = -7.5 - row_number()) %>% separate_rows(blocks, sep = ";") %>%
        separate(blocks, c("a", "b"), sep = "-", convert = TRUE) %>%
        filter(b >= min(xm0$start), a <= max(xm0$end))
    sp = bl %>% group_by(row, class) %>% summarise(a = min(a), b = max(b), .groups = "drop")
    xm = xm0 %>% left_join(trk, by = c("track" = "key"))
    xs = xm %>% group_by(y) %>% summarise(a = min(start), b = max(end), .groups = "drop")
    p1 = ggplot() +
        geom_segment(data = sp, aes(fx(a), row, xend = fx(b), yend = row, colour = class), linewidth = 0.1) +
        geom_segment(data = bl, aes(fx(a), row, xend = fx(b), yend = row, colour = class), linewidth = 0.55) +
        geom_segment(data = xs, aes(fx(a), y, xend = fx(b), yend = y), linewidth = 0.2) +
        geom_rect(data = xm %>% filter(kind == "exon"), aes(xmin = fx(start), xmax = fx(end), ymin = y - 0.8, ymax = y + 0.8), fill = "white", colour = "black", linewidth = 0.2) +
        geom_rect(data = xm %>% filter(kind == "cds"), aes(xmin = fx(start), xmax = fx(end), ymin = y - 1.1, ymax = y + 1.1), fill = ISOFORM_COLORS[1], linewidth = 0) +
        geom_vline(xintercept = fx(e$peak), linetype = 2, linewidth = 0.25, colour = "#1a9850") +
        scale_colour_manual(values = READ_COLS, drop = FALSE, name = NULL) +
        scale_y_continuous(breaks = trk$y, labels = trk$track) +
        labs(x = NULL, y = NULL, tag = tag, title = sprintf("%s, %s: %d → %d aa\nS-site reads %d%% capped", g, DS_LABEL[[e$dataset]],
                                                            e$aa_total, e$aa_remaining, round(100 * e$cap_frac))) +
        theme_paper() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                              plot.title = element_text(size = 7, face = "italic"), legend.position = "none")
    # read 5' ends over the locus (aggregate counts), on p1's compressed axis; dashed lines: L's TSS (grey) and S's site (green)
    xr = range(fx(range(xm0$start, xm0$end)))
    s2 = a %>% filter(via_s_donor) %>% mutate(x = fx(start5)) %>% filter(x >= xr[1], x <= xr[2])   # reads crossing S's first donor
    bw = diff(xr) / 120
    p2 = ggplot(s2, aes(x, fill = class)) +
        geom_histogram(binwidth = bw, boundary = xr[1]) +
        geom_vline(xintercept = fx(e$l_tss), linetype = 2, linewidth = 0.25, colour = "grey40") +
        geom_vline(xintercept = fx(e$peak), linetype = 2, linewidth = 0.25, colour = "#1a9850") +
        scale_fill_manual(values = READ_COLS, drop = FALSE, name = NULL) +
        scale_y_sqrt() + coord_cartesian(xlim = xr) +
        labs(x = NULL, y = "Read 5' ends") +
        theme_paper() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
    cg = atl %>% filter(example == .env$id, atlas == "CAGE") %>% mutate(x1 = pmin(fx(start), fx(end)), x2 = pmax(fx(start), fx(end))) %>%
        filter(x2 >= xr[1], x1 <= xr[2]) %>% mutate(m = (x1 + x2) / 2, w = pmax(x2 - x1, diff(xr) / 300))
    p1 = p1 + geom_rect(data = cg, aes(xmin = m - w / 2, xmax = m + w / 2, ymin = -1.1, ymax = -0.1), fill = "black") +
        geom_vline(xintercept = fx(e$l_tss), linetype = 2, linewidth = 0.25, colour = "grey40") + coord_cartesian(xlim = xr)
    p1 / p2 + plot_layout(heights = c(2.6, 1.1))
}
ggpubr_free_legend = function(p) { g = ggplotGrob(p); g$grobs[[which(sapply(g$grobs, function(x) grepl("guide-box", x$name)))[1]]] }
legend_src = ggplot(tibble(x = 1, y = 1, class = factor(names(READ_COLS), names(READ_COLS))), aes(x, y, colour = class)) +
    geom_point(shape = 15) + scale_colour_manual(values = READ_COLS, name = NULL) + theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 5.5), legend.key.size = unit(3, "mm"))

top = wrap_elements(pa)
mid = wrap_elements(pb) + wrap_elements(pc) + plot_layout(widths = c(1.3, 1))
exid = function(g) ex$example[ex$gene == g & ex$example != "ftd"][1]
genes = wrap_elements(gene_panel(exid("COMT"), "d")) + wrap_elements(gene_panel("ftd", "e")) + wrap_elements(gene_panel(exid("CORO1B"), "f")) + wrap_elements(gene_panel(exid("RECQL4"), "g")) + plot_layout(ncol = 2)
legend = wrap_elements(ggpubr_free_legend(legend_src))
fig = wrap_plots(list(top, mid, genes, legend), ncol = 1, heights = c(1.1, 1.6, 4.3, 0.12)) & theme(plot.margin = margin(1, 1, 1, 1, "mm"), plot.tag = element_text(size = 8, face = "bold"))
ggsave("SuppFig_ISM_internal_TSS.pdf", fig, width = WIDTH_IN, height = HEIGHT_IN, device = grDevices::pdf, useDingbats = FALSE)
embed_fonts("SuppFig_ISM_internal_TSS.pdf")

values = bind_rows(
    td %>% transmute(dataset, measure = paste0("tier_", tier), value = peaks),
    ti %>% transmute(dataset = "bulk", measure = paste0("own_ism_isoforms_", tool), value = ism_isoforms),
    ti %>% transmute(dataset = "bulk", measure = paste0("own_ism_supported_", tool), value = supported),
    ti %>% transmute(dataset = "bulk", measure = paste0("own_ism_truncation_like_", tool), value = truncation_like),
    ti %>% transmute(dataset = "bulk", measure = paste0("own_ism_pct_truncation_", tool), value = pct_truncation_of_judged),
    ab %>% transmute(dataset, measure = paste0("cage10_", group), value = obs),
    ab %>% transmute(dataset, measure = paste0("cage10_bg_", group), value = bg),
    read_tsv(file.path(IN, "qualifying_summary.tsv"), show_col_types = FALSE) %>% pivot_longer(everything(), names_to = "measure") %>% transmute(dataset = "all", measure = paste0("qualifying_", measure), value = value),
    ex %>% transmute(dataset = gene, measure = "aa_total", value = aa_total), ex %>% transmute(dataset = gene, measure = "aa_remaining", value = aa_remaining),
    ex %>% transmute(dataset = gene, measure = "n_site", value = n_site), ex %>% transmute(dataset = gene, measure = "cap_pct", value = 100 * cap_frac),
    tibble(dataset = FTD_GENE, measure = c("L_fsm_reads", "S_fsm_reads", "L_tss_reads", "neuron_pct"),
           value = c(fc$L_fsm_reads, fc$S_fsm_reads, fc$L_tss_reads, 100 * fc$fragile_share)))
write_tsv(values, "SuppFig_ISM_internal_TSS.values.tsv")
fs = bind_rows(fit_stats, .id = "example")
print(fs %>% pivot_wider(names_from = class, values_from = n), width = Inf)
cat("wrote SuppFig_ISM_internal_TSS.pdf and SuppFig_ISM_internal_TSS.values.tsv\n")
