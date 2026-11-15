# 5' completeness curves K(d), as in tss_5prime_completeness.Rmd section A: for each gene and
# group of cells, S(d) = share of its polyA-anchored reads whose 5' extent from the 3' end is
# >= d (reads reaching the TSS count as full length L); K(d) = mean S(d) over genes with
# L >= d and >= 20 anchored reads in the group. Groups: "all" (whole dataset) or "cluster".
suppressPackageStartupMessages(library(data.table))

k_curves = function(datasets, grid = seq(0, 6000, by = 25)) {
    genes = rbindlist(lapply(datasets, function(ds) fread(paste0(ds, ".genes.tsv"))[, dataset := ds]), fill = TRUE)
    use = genes[status == "ok" & cage %in% c("single", "none")]
    ext = rbindlist(lapply(datasets, function(ds) fread(paste0(ds, ".extent.tsv.gz"))[, dataset := ds]))
    ext = merge(ext, use[, .(dataset, gene_key, L)], by = c("dataset", "gene_key"))
    hist = ext[extent_hist != "", {
        p = strsplit(strsplit(extent_hist, ";")[[1]], ":")
        .(bin = as.integer(sapply(p, `[`, 1)), n = as.integer(sapply(p, `[`, 2)))
    }, by = .(dataset, gene_key, cluster)]
    ext[, all := "all"]
    hist = merge(hist, ext[, .(dataset, gene_key, cluster, all)], by = c("dataset", "gene_key", "cluster"))
    surv = function(group_cols) {
        tot = ext[, .(n = sum(n_anchored), r = sum(n_reach), L = L[1]), by = c("dataset", "gene_key", group_cols)][n >= 20]
        h = hist[, .(n = sum(n)), by = c("dataset", "gene_key", group_cols, "bin")]
        hc = h[order(-bin), .(bin, cum = cumsum(n)), by = c("dataset", "gene_key", group_cols)][order(bin)]
        pts = tot[, .(d = grid[grid <= L]), by = c("dataset", "gene_key", group_cols)]
        pts = merge(pts, tot, by = c("dataset", "gene_key", group_cols))
        pts[, bin := d]
        pts[, cum := hc[pts, on = c("dataset", "gene_key", group_cols, "bin"), roll = -Inf, x.cum]]
        pts[is.na(cum), cum := 0]
        pts[, S := (r + cum) / n]
        pts
    }
    K = function(S, group_cols) S[, .(K = mean(S), genes = .N), by = c("dataset", group_cols, "d")][genes >= 20]
    list(all = K(surv("all"), "all"), cluster = K(surv("cluster"), "cluster"))
}

# K at distances d for one dataset x group (linear interpolation, flat beyond the ends)
k_at = function(Kt, d) {
    if (!nrow(Kt)) return(rep(NA_real_, length(d)))
    approx(Kt$d, Kt$K, xout = pmax(0, d), rule = 2)$y
}
