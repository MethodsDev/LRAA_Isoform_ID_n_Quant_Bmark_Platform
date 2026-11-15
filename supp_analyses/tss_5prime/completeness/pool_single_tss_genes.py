#!/usr/bin/env python3
"""Pool single-TSS gene tables (genes.tsv of tss_5prime_completeness.py) from several
datasets: one row per main transcript (first dataset wins), with the gene's TSS position
(site_pos) from that dataset's site table.

usage: pool_single_tss_genes.py DATASET:genes.tsv:sites.tsv [...] > pooled.tsv"""
import csv
import sys

w, seen = None, set()
for spec in sys.argv[1:]:
    ds, genes, sites = spec.split(":")
    pos = {r["gene_key"]: r["pos"] for r in csv.DictReader(open(sites), delimiter="\t") if r["kind"] == "TSS"}
    for r in csv.DictReader(open(genes), delimiter="\t"):
        t = r["main_transcript"]
        if not t or t in seen or r["gene_key"] not in pos:
            continue
        seen.add(t)
        r["site_pos"], r["dataset"] = pos[r["gene_key"]], ds
        if w is None:
            w = csv.DictWriter(sys.stdout, fieldnames=list(r), delimiter="\t", lineterminator="\n")
            w.writeheader()
        w.writerow(r)
