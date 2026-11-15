#!/usr/bin/env python3
"""HC events between halves of one cluster (NULL, all false) and between different clusters' A halves (REAL),
LRAA v0.51.0 pooled final quant, from run_null_v051.sh's work/s{1,2}/. Rows are genome-wide and restricted to
chr1/6/11/22. Writes tables/N2_v051_null_vs_real.{tsv,md}."""
import os
import pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
BASE = ["0", "1", "4", "10"]
NULLP = [(c, str(100 + int(c))) for c in BASE]
CHROMS = {"chr1", "chr6", "chr11", "chr22"}
REALP = [(a, b) for i, a in enumerate(BASE) for b in BASE[i + 1:]]
def evkey(a, b):
    x = sorted([str(a).replace("Cluster_", ""), str(b).replace("Cluster_", "")], key=int)
    return f"{x[0]}-{x[1]}"
rows = []
for seed in (1, 2):
    for K in ("SplicePattern", "IsoformTermini"):
        p = f"{H}/work/s{seed}/PBMCs.dexseq.{K}"
        if not os.path.exists(f"{p}.events.tsv"):
            continue
        g = pd.read_csv(f"{p}.genes.tsv", sep="\t")
        e = pd.read_csv(f"{p}.events.tsv", sep="\t", keep_default_na=False)
        e["pair"] = [evkey(x, y) for x, y in zip(e.cluster_A, e.cluster_B)] if len(e) else []
        e["hc"] = e.high_confidence.astype(str) == "True" if len(e) else []
        stable = int((g.stable.astype(str).str.upper() == "TRUE").sum())
        # genome-wide, and restricted to part 4's chromosomes for a like-for-like row
        for arm, chroms in (("v0.51.0 pooled (genome-wide)", None), ("v0.51.0 pooled (chr1/6/11/22)", CHROMS)):
            ea = e if chroms is None or not len(e) else e[e.gene_key.str.split("|").str[1].isin(chroms)]
            ga = g if chroms is None else g[g.gene_key.str.split("|").str[1].isin(chroms)]
            for kind, pairs in (("NULL", NULLP), ("REAL", REALP)):
                keys = {evkey(a, c) for a, c in pairs}
                ee = ea[ea.pair.isin(keys)] if len(ea) else ea
                hc = int(ee.hc.sum()) if len(ee) else 0
                rows.append(dict(seed=seed, kind=K, arm=arm, pair_kind=kind, n_pairs=len(keys),
                                 genes_tested=len(ga), stable_genes_all=int((ga.stable.astype(str).str.upper() == "TRUE").sum()),
                                 events=len(ee), HC_events=hc,
                                 HC_genes=ee.loc[ee.hc, "gene_key"].nunique() if hc else 0, HC_per_pair=hc / len(keys)))
v51 = pd.DataFrame(rows)
t = v51
os.makedirs(f"{H}/tables", exist_ok=True)
t.to_csv(f"{H}/tables/N2_v051_null_vs_real.tsv", sep="\t", index=False)
s = (t.groupby(["kind", "arm", "pair_kind"]).agg(HC_events=("HC_events", "sum"), pairs=("n_pairs", "sum"))
     .assign(HC_per_pair=lambda d: (d.HC_events / d.pairs).round(2)).reset_index()
     .pivot_table(index=["kind", "arm"], columns="pair_kind", values="HC_per_pair").reset_index())
s.columns.name = None
s = s.rename(columns={"NULL": "NULL HC/pair (false; lower is better)", "REAL": "REAL HC/pair"})
with open(f"{H}/tables/N2_v051_null_vs_real.md", "w") as fh:
    fh.write("High-confidence DTU events per cluster pair, summed over 2 seeds. NULL = halves of one cluster (every call is "
             "false); REAL = different clusters at the same half depth. LRAA v0.51.0 ref-guided sc cluster-guided models, "
             "pooled final quant.\n\n")
    fh.write("| " + " | ".join(s.columns) + " |\n|" + "---|" * len(s.columns) + "\n")
    for r in s.itertuples(index=False):
        fh.write("| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |\n")
print(open(f"{H}/tables/N2_v051_null_vs_real.md").read())
