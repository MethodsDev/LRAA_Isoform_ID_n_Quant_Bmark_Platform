#!/usr/bin/env python3
"""
Supplementary Table 1: differential usage events between cell clusters (DEXSeq on
pseudo-replicates) for the PBMC and FTD (ALS LP00714) single-cell datasets.

Reads each study's event tables (the study submodules, with their Zenodo bundles unpacked):
  __Splice_Pattern_Usage_Analysis/<prefix>.dexseq.{SplicePattern,IsoformTermini}.events.tsv
  __Site_Usage_Analysis/<prefix>.dexseq.{TSS,PolyA}.events.tsv
the per-gene tables (<...>.genes.tsv) for the summary sheet, and for the site events the splicing
class of each site pair (<prefix>.dexseq.site_pairs.splicing.tsv, classify_site_pairs_by_splicing.py)
and the LRAA isoforms carrying each site (<prefix>.sites.tsv).

The paper finds alternative termini with the site test and relates each switch to its isoforms;
the IsoformTermini test is not reported. FTD TSS site events are reported with a caveat (CAVEATS):
LRAA calls TSS sites only where the reads carry the 5' cap signature, which removes most
5'-truncation sites in the nuclear reads, but neuron-favours-downstream switches remain
over-represented.

Writes
  SuppTable1_differential_usage.xlsx: a README sheet (definitions, version, column key), a Summary
    sheet (genes per level), a Clusters sheet (each cluster's cell type and cells), and one sheet per
    dataset x level with its high-confidence events (cell types next to the cluster ids);
  differential_usage_all_events.zip: every event (one TSV per dataset x level, high-confidence
    first) plus README.txt and clusters.tsv, kept in the repo (the paper points to it with the Supplementary Code).
    All events in the workbook would be ~55 MB.
The columns linking site events to the chi-square isoform test
(isoform_*_DTU_*) are left out: the paper does not use that test.

Run from this directory:  python3 make_dexseq_supp_table.py
"""

import os
import re
import sys
import zipfile

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = "SuppTable1_differential_usage.xlsx"
OUT_ALL = "differential_usage_all_events.zip"

DATASETS = [
    # name, study submodule (its Zenodo bundle unpacked), file prefix
    ("PBMC", "LRAA_paper-PBMCs_study", "PBMCs"),
    ("FTD", "LRAA_paper-FTD_study", "LP00714"),
]
LRAA_VERSION = "v0.51.0"   # the LRAA release of the study bundles
# the LRAA gtf the site analysis reads (__Site_Usage_Analysis/Makefile GTF): a model's gene symbol
SITE_GTF = {"PBMC": "pbmcs_refGuided/pbmcs_refGuided_sc_cluster_guided/data/PBMCs_pbio_CG_refguided.withGeneSymbols.gtf",
            "FTD": "als_refGuided/als_refGuided_sc_cluster_guided/data/LP00714.refguided.withGeneSymbols.gtf"}
# cell clusters (fixed, RefQuantOnly): cells per cluster, the cluster-level cell type annotation (from CAS labels),
# and the paper's short names
# (CLUSTER_INFO in the dataset's __Paper_Figures/paper_fig_dataset.R, as in the main figures)
CLUSTER_CELLS = {"PBMC": ("pbmcs_refQuantOnly/PBMCS.refQuantOnly.cell_clusters_and_cell_types.wUMAP.tsv", "cell_type_annot"),
                 "FTD": ("als_refQuantOnly/LP00714.ALS.refQuantOnly.cell_clusters_and_cell_types.wUMAP.tsv",
                         "cell_type_simplified")}
CLUSTER_LABELS = "__Paper_Figures/paper_fig_dataset.R"
LEVELS = [
    # level, analysis dir, what a feature is
    ("SplicePattern", "__Splice_Pattern_Usage_Analysis",
     "splice pattern (isoforms with the same intron chain, EM-assigned counts summed); group = gene"),
    ("TSS", "__Site_Usage_Analysis", "LRAA TSS site (read 5' ends within 25 nt); group = gene"),
    ("PolyA", "__Site_Usage_Analysis", "LRAA PolyA site (read 3' ends within 25 nt); group = gene"),
]
DROP = ["isoform_alt_termini_DTU_same_pair", "isoform_any_DTU_same_pair", "isoform_alt_termini_DTU_gene"]
# (dataset, level) reported with a caveat: in the workbook's Summary and README, and the zip's README
CAVEATS = {("FTD", "TSS"): "interpret with care: nuclear reads are often 5'-truncated. LRAA calls TSS sites "
                           "only with the 5' cap signature, which removes most truncation sites, but switches to "
                           "a downstream TSS in neurons remain over-represented (see the paper)"}

# column key. gained_* / lost_* columns describe the gained or lost feature; *_A / *_B the two clusters.
COLUMNS = {
    "gene_key": "gene (symbol|chrom|strand)",
    "splicing_class": "site events, from the reads at the two sites: alt_terminal_usage (tandem: same terminal exon and "
                      "splicing), alt_splicing:terminal_exon (alternative first / last exon, intronic PolyA, retained "
                      "intron), alt_splicing:internal (an internal intron differs), unspliced_site, unresolved",
    "gene_symbol": "gene symbol",
    "cluster_A": "first cluster of the pair",
    "cluster_B": "second cluster; the gained feature's usage is higher in B than in A",
    "cell_type_A": "cluster A's cell type (short name as in the figures; Clusters sheet)",
    "cell_type_B": "cluster B's cell type (short name as in the figures; Clusters sheet)",
    "gained_feature": "the feature whose usage rises most significantly from A to B",
    "lost_feature": "the feature whose usage falls most from A to B",
    "gained_site": "the site whose usage rises most significantly from A to B",
    "lost_site": "the site whose usage falls most from A to B",
    "abs_delta": "|delta usage| of the gained feature between the clusters (usage = share of the group's counts)",
    "high_confidence": "reciprocal switch, >= 5 unique FSM reads for both features, no flags",
    "switch_class": "from the two features' counts per million (1.5-fold): reciprocal (gained up, lost down), "
                    "concordant (both change the same way), one site changes, neither",
    "flags": "quality flags; any flag excludes the event from high confidence (see README)",
    "n_features_tested": "features of the group tested",
    "n_features_sig": "features of the group with a significant change in this pair",
    "n_sites": "sites of the gene tested",
    "n_sites_sig": "sites of the gene with a significant change in this pair",
    "n_seeds_significant": "of the 5 random dealings into pseudo-replicates, how many give gene q < 0.05",
    "median_q": "median gene-level q (perGeneQValue) over the 5 seeds",
    "group_reads_A": "counts of the group in cluster A (tested features)",
    "group_reads_B": "counts of the group in cluster B (tested features)",
    "gene_reads_A": "read ends of the gene at tested sites in cluster A",
    "gene_reads_B": "read ends of the gene at tested sites in cluster B",
    "both_FSM": ">= 5 unique FSM reads for both features",
    "both_sites_FSM": ">= 5 unique FSM reads at both sites (the most of any isoform carrying the site)",
    "event_type": "site events: tandem_3UTR, intronic_PolyA, alt_last_exon (PolyA); tandem_TSS, alt_first_exon (TSS)",
    "separation": "distance between the gained and lost sites (bp)",
    "gained_is_downstream": "the gained site lies downstream of the lost one on the transcript strand",
    "TSS_differs": "SplicePattern: the two patterns' main TSS differ",
    "PolyA_differs": "SplicePattern: the two patterns' main PolyA sites differ",
    "exon_count_differs": "SplicePattern: the two patterns differ in exon count",
}
PREFIXED = {
    "usage": "usage (share of the group's counts) in cluster {c}",
    "reads": "counts in cluster {c}",
    "padj": "pairwise adjusted P (Benjamini-Hochberg over all pairwise tests)",
    "uniq_FSM": "unique full-length (FSM) reads of the feature's isoforms, summed over the cluster quantifications",
    "log2FC": "log2 fold change of counts per million between the clusters (the switch_class basis)",
    "n_isoforms": "LRAA isoforms carrying the feature",
    "n_seeds_confirmed": "seeds in which stageR confirms the feature",
    "transcript_ids": "LRAA isoforms making up the feature (site events: the gene's isoforms carrying the site)",
    "num_exons": "exon count",
    "main_TSS": "TSS site of the feature's most expressed isoform",
    "main_PolyA": "PolyA site of the feature's most expressed isoform",
    "TSS_site": "TSS site of the feature",
    "PolyA_site": "PolyA site of the feature",
    "PolyA_site_evidence": "PAS hexamer, PolyASite 2.0 sample fraction, A's in the 20 genomic bases downstream",
    "PolyA_site_supported": "PAS hexamer, or PolyASite 2.0 cluster used by >= 10% of samples, or <= 7 downstream A's",
    "site_evidence": "PAS hexamer, PolyASite 2.0 sample fraction, A's in the 20 genomic bases downstream",
    "site_supported": "PAS hexamer, or PolyASite 2.0 cluster used by >= 10% of samples, or <= 7 downstream A's",
    "pos": "site position",
    "source": "LRAA site catalog the site comes from",
    "n_merged": "LRAA sites merged into this one",
    "max_isoform_uniq_FSM": "most unique FSM reads of any isoform carrying the site",
    "monoexonic_only": "the site is carried only by single-exon isoforms",
    "PAS": "polyadenylation signal hexamer upstream of the site",
    "downstream_A_of_20": "A's in the 20 genomic bases downstream of the site",
}

README = [
    ("Supplementary Table 1", "Differential usage between cell clusters, DEXSeq on pseudo-replicates"),
    ("LRAA version", "{version} (reference-guided, cluster-guided isoform set)"),
    ("Datasets", "PBMC: Kinnex 10x 3' scRNA-seq, 14 clusters. FTD: frontal cortex Kinnex 10x 3' snRNA-seq, 14 clusters."),
    ("Levels", ""),
] + [(f"  {lvl}", what) for lvl, _, what in LEVELS] + [
    ("Test", "Each cluster's cells dealt at random into 3 pseudo-replicates, counts summed per feature. DEXSeq "
             "~ sample + exon + cluster:exon vs ~ sample + exon over all clusters with >= 30 cells (features as "
             "exon bins, groups as genes); gene q from perGeneQValue; stageR confirms the features carrying the change. "
             "Features tested if >= 10% of the group's counts in some cluster with >= 20 group counts there, and >= 10 "
             "counts overall; groups need >= 2 features."),
    ("Significant gene", "q < 0.05 in at least 4 of 5 random dealings (seeds)"),
    ("Event (one row)", "significant gene x cluster pair, pairwise test on the pair's pseudo-replicates with the "
                        "all-cluster dispersions: a feature at adjusted P < 0.05, |delta usage| >= 0.2, and >= 20 group "
                        "counts in each cluster. Oriented so the largest significant change is a gain in cluster_B."),
    ("High confidence", "reciprocal switch class, >= 5 unique FSM reads for both features, no flags."),
    ("Flags", "monoexonic (site only on single-exon isoforms); close_sites (< 30 bp apart); A_rich_downstream "
              "(>= 12 A's in the 20 bases past a PolyA site); polyA_site_unsupported (no PAS, no PolyASite 2.0 "
              "support, > 7 downstream A's: internal priming); downstream_TSS_no_FSM (gained TSS downstream of the "
              "lost one without full-length reads: 5' truncation)"),
    ("Clusters sheet", "per dataset: each cluster's short name (as in the figures), lineage, cell type annotation "
                       "(the Cell Annotation Service label most common among its cells; FTD clusters 12 and 13 "
                       "relabelled from marker genes, see Methods) and number of cells"),
    ("Summary sheet", "per dataset x level: genes tested, significant genes (q < 0.05 in >= 4 of 5 seeds), genes with "
                      "events and events (significant genes x cluster pairs), high-confidence genes and events"),
    ("FTD TSS", CAVEATS[("FTD", "TSS")]),
    ("Caveat", "Pseudo-replicates come from one sample: they measure cell-to-cell sampling noise, not donor "
               "variation, so P values are optimistic. Rank by abs_delta and the high-confidence filter."),
    ("Source", "LRAA util/sc/site_usage/ (site_usage_dexseq.R, annotate_site_usage_events.py, "
               "annotate_feature_usage_events.py); built by __PaperFigs/supp_tables/make_dexseq_supp_table.py"),
]


def round_floats(df):
    """3 significant digits for P/q values, 3 decimals for the rest: the workbook stays a usable size."""
    for c in df.select_dtypes("float").columns:
        if re.search(r"padj|_q$|pvalue", c):
            df[c] = df[c].map(lambda x: float(f"{x:.3g}") if pd.notna(x) else x)
        else:
            df[c] = df[c].round(3)
    return df


def describe(col):
    if col in COLUMNS:
        return COLUMNS[col]
    m = re.match(r"(gained|lost)_(.+?)(?:_([AB]))?$", col)
    if m:
        side, stem, cl = m.groups()
        for key, text in PREFIXED.items():
            if stem == key:
                return f"{side} feature: " + text.format(c=cl or "")
    return ""


def cluster_table(ds, root):
    """Cluster id (Cluster_N) -> short name, lineage, cell type annotation, cells."""
    cells_file, annot_col = CLUSTER_CELLS[ds]
    cells = pd.read_csv(os.path.join(ROOT, root, cells_file), sep="\t")
    labels = {}
    with open(os.path.join(ROOT, root, CLUSTER_LABELS)) as fh:
        block = fh.read().split("CLUSTER_INFO = tibble::tribble(", 1)[1].split("\n\n", 1)[0]
    for m in re.finditer(r'^\s*(\d+),\s*"([^"]+)",\s*"([^"]+)"', block, re.M):
        labels[int(m.group(1))] = (m.group(2), m.group(3))
    rows = []
    for cl, d in cells.groupby("seurat_clusters"):
        if cl not in labels:
            sys.exit(f"{ds}: cluster {cl} has no short name in {CLUSTER_LABELS}")
        rows.append({"dataset": ds, "cluster": f"Cluster_{cl}", "cell_type": labels[cl][0], "lineage": labels[cl][1],
                     "annotation": d[annot_col].value_counts().index[0], "cells": len(d)})
    if len(rows) != len(labels):
        sys.exit(f"{ds}: {CLUSTER_LABELS} names {len(labels)} clusters, the cell table has {len(rows)}")
    # in the figures' order (CLUSTER_INFO)
    order = list(labels)
    return pd.DataFrame(sorted(rows, key=lambda r: order.index(int(r["cluster"].split("_")[1]))))


def model_symbols(gtf):
    """LRAA model id (t:chr:strand:comp-N:iso-M) -> gene symbol, from the gtf's SYMBOL^ transcript ids."""
    sym = {}
    with open(gtf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "transcript":
                tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
                sym[tid.split("^")[-1]] = tid.split("^")[0]
    return sym


def annotate_sites(ev, root, prefix, lvl, sym):
    """Site events: the splicing class of the site pair and the gene's own LRAA isoforms carrying each
    site (a site can also be carried by another gene's models in the same LRAA component)."""
    base = os.path.join(ROOT, root, "__Site_Usage_Analysis", prefix)
    spl = pd.read_csv(f"{base}.dexseq.site_pairs.splicing.tsv", sep="\t")
    spl = spl[spl.kind == lvl]
    cls = {(g, frozenset((a, b))): c for g, a, b, c in spl[["gene_key", "gained_site", "lost_site", "splicing_class"]].itertuples(index=False)}
    ev.insert(ev.columns.get_loc("lost_site") + 1, "splicing_class",
              [cls.get((g, frozenset((a, b))), "") for g, a, b in ev[["gene_key", "gained_site", "lost_site"]].itertuples(index=False)])
    tx = pd.read_csv(f"{base}.sites.tsv", sep="\t", usecols=["site_id", "transcript_ids"]).set_index("site_id").transcript_ids
    def own(ids, gene):
        ids = str(ids).split(",")
        mine = [i for i in ids if sym.get(i) == gene]
        return ",".join(mine or ids)
    ev["gained_transcript_ids"] = [own(tx.get(s, ""), g) for s, g in zip(ev.gained_site, ev.gene_symbol)]
    ev["lost_transcript_ids"] = [own(tx.get(s, ""), g) for s, g in zip(ev.lost_site, ev.gene_symbol)]
    return ev


def main():
    version = LRAA_VERSION

    symbols = {ds: model_symbols(os.path.join(ROOT, root, SITE_GTF[ds])) for ds, root, _ in DATASETS}
    clusters = pd.concat([cluster_table(ds, root) for ds, root, _ in DATASETS], ignore_index=True)
    sheets, summary, columns = [], [], {}
    for ds, root, prefix in DATASETS:
        cell_type = clusters[clusters.dataset == ds].set_index("cluster").cell_type
        for lvl, sub, _ in LEVELS:
            base = os.path.join(ROOT, root, sub, f"{prefix}.dexseq.{lvl}")
            hcc = "high_confidence"
            ev = pd.read_csv(base + ".events.tsv", sep="\t", low_memory=False)
            genes = pd.read_csv(base + ".genes.tsv", sep="\t")
            ev = ev.drop(columns=[c for c in DROP if c in ev.columns])
            for side in "AB":
                ct = ev[f"cluster_{side}"].map(cell_type)
                assert ct.notna().all(), f"{ds} {lvl}: clusters without a cell type"
                ev.insert(ev.columns.get_loc(f"cluster_{side}") + 1, f"cell_type_{side}", ct)
            if sub == "__Site_Usage_Analysis":
                ev = annotate_sites(ev, root, prefix, lvl, symbols[ds])
            ev = ev.sort_values([hcc, "abs_delta"], ascending=[False, False])
            ev = round_floats(ev)
            sheets.append((f"{ds}_{lvl}", ev, hcc))
            hc = ev[ev[hcc].astype(str) == "True"]
            sig = set(genes.gene_key[genes.stable.astype(str).eq("True")])
            ev_sig = ev[ev.gene_key.isin(sig)]
            assert hc.gene_key.isin(sig).all()
            summary.append({"dataset": ds, "level": lvl, "caveat": CAVEATS.get((ds, lvl), ""), "genes tested": len(genes),
                            "significant genes": len(sig),
                            "genes with events": ev_sig.gene_key.nunique(), "events": len(ev_sig),
                            "high-confidence genes": hc.gene_key.nunique(), "high-confidence events": len(hc)})
            for c in ev.columns:
                columns.setdefault(c, describe(c))

    with pd.ExcelWriter(OUT, engine="xlsxwriter") as xw:
        readme = pd.DataFrame([(k, v.format(version=version)) for k, v in README], columns=["item", "description"])
        readme.to_excel(xw, sheet_name="README", index=False)
        key = pd.DataFrame(sorted(columns.items()), columns=["column", "description"])
        key.to_excel(xw, sheet_name="README", index=False, startrow=len(readme) + 2)
        xw.sheets["README"].set_column(0, 0, 32)
        xw.sheets["README"].set_column(1, 1, 120)
        pd.DataFrame(summary).to_excel(xw, sheet_name="Summary", index=False)
        clusters.to_excel(xw, sheet_name="Clusters", index=False)
        xw.sheets["Clusters"].set_column(4, 4, 60)
        for name, ev, hcc in sheets:
            ev[ev[hcc].astype(str) == "True"].to_excel(xw, sheet_name=name, index=False)
            xw.sheets[name].freeze_panes(1, 0)

    with zipfile.ZipFile(OUT_ALL, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        text = "\n".join(f"{k}\t{v.format(version=version)}" for k, v in README)
        text += "\n\nAll events, one file per dataset x level, high-confidence events first."
        text += "".join(f"\n{ds}_{lvl}: {why}" for (ds, lvl), why in CAVEATS.items())
        text += "\n\ncolumn\tdescription\n"
        text += "\n".join(f"{c}\t{d}" for c, d in sorted(columns.items())) + "\n"
        zf.writestr("README.txt", text)
        zf.writestr("clusters.tsv", clusters.to_csv(sep="\t", index=False))
        for name, ev, _ in sheets:
            zf.writestr(f"{name}.events.tsv", ev.to_csv(sep="\t", index=False))

    missing = [c for c, d in columns.items() if not d]
    if missing:
        print("WARNING: columns without a description: " + ", ".join(missing), file=sys.stderr)
    print(f"wrote {OUT} and {OUT_ALL} ({version}; {sum(len(e) for _, e, _ in sheets):,} events, "
          f"{sum(r['high-confidence events'] for r in summary):,} high-confidence)")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
