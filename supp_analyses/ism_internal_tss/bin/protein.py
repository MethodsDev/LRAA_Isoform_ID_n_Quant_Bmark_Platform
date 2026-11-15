#!/usr/bin/env python3
"""Protein consequence per S/L pair. L = GENCODE v47 chain match -> its CDS. Otherwise an ORF from
pyTD (PyTransDecoder, run_orfs.py) if --orf is given (L chain -> genomic CDS, see run_orfs.py).
usage: protein.py DATASET OUT.tsv.gz [orf_table.tsv]"""
import sys, pickle
import numpy as np, pandas as pd, pysam
from common import *
from protein_lib import *
ds, out = sys.argv[1:3]
orf = {}
if len(sys.argv) > 3 and os.path.exists(sys.argv[3]):
    o = pd.read_csv(sys.argv[3], sep="\t")
    for r in o.itertuples():
        orf[(r.chrom, r.strand, r.l_chain)] = [tuple(map(int, x.split("-"))) for x in r.cds.split(",")]
here = os.path.dirname(os.path.abspath(__file__))
gc = pickle.load(open(os.path.join(here, "..", "out", "gencode47.pkl"), "rb"))
idx = {}
for t, r in gc.items():
    if r["chain"]: idx.setdefault((r["chrom"], r["strand"], r["chain"]), []).append(t)
def pick(ts):
    cod = [t for t in ts if gc[t]["cds"] and gc[t]["start"] and "cds_start_NF" not in gc[t]["tags"] and "mRNA_start_NF" not in gc[t]["tags"]]
    pc = [t for t in cod if gc[t]["ttype"] == "protein_coding"] or cod
    if pc:
        return max(pc, key=lambda t: ("MANE_Select" in gc[t]["tags"], "Ensembl_canonical" in gc[t]["tags"], sum(b - a for a, b in gc[t]["cds"])))
    return None
fa = pysam.FastaFile(DATASETS[ds]["genome"])
p = pd.read_csv(os.path.join(here, "..", "out", "pairs", f"{ds}.pairs_j.tsv.gz"), sep="\t")
s = pd.read_csv(os.path.join(here, "..", "out", "scored", f"{ds}.sites.tsv.gz"), sep="\t", usecols=["regime", "chrom", "strand", "s_chain", "s_tss_group_lo", "peak"])
p = p.merge(s, on=["regime", "chrom", "strand", "s_chain", "s_tss_group_lo"])
cache = {}; rows = []
for r in p.itertuples():
    k = (r.chrom, r.strand, r.s_chain, r.peak, r.s_end3, r.l_chain)
    if k not in cache:
        lc = parse_chain(r.l_chain)
        ts = idx.get((r.chrom, r.strand, lc), [])
        t = pick(ts)
        d = dict(EMPTY)
        peak = int(r.peak)
        if t:
            cds = gc[t]["cds"]; route = "GENCODE"; lab = t
        elif (r.chrom, r.strand, r.l_chain) in orf:
            cds = orf[(r.chrom, r.strand, r.l_chain)]; route = "pyTD"; lab = ""
        else:
            cds = None; lab = ts[0] if ts else ""; route = "L_noncoding_or_NF" if ts else "L_unannotated"
        if cds is not None:
            d = classify(fa, r.chrom, r.strand, parse_chain(r.s_chain), peak, int(r.s_end3), cds)
            if route == "pyTD":                  # its CDS includes the stop codon
                d["aa_total"] -= 1; d["aa_remaining"] -= 1
        else:
            d["protein_class"] = "L_noncoding" if ts else "L_unannotated_no_ORF"
        d.update(cds_route=route, L_gencode_tx=lab, L_gencode_type=(gc[lab]["ttype"] if lab else ""), L_mane=bool(lab and "MANE_Select" in gc[lab]["tags"]), protein_id=(gc[lab]["pid"] if lab else ""))
        cache[k] = d
    rows.append(cache[k])
pd.concat([p.reset_index(drop=True), pd.DataFrame(rows)], axis=1).to_csv(out, sep="\t", index=False)
print(p.shape, pd.Series([x["protein_class"] for x in rows]).value_counts().to_dict())
