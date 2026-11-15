#!/usr/bin/env python3
"""Do the tools recover the internal-TSS isoforms? Targets come from the tool-neutral pool (all paper tools + GENCODE), kept by
read and atlas evidence only. usage: tool_recall.py"""
import pickle
import numpy as np, pandas as pd
from common import *
here = os.path.dirname(os.path.abspath(__file__)); pool = os.path.join(here, "..", "out", "pool")
gc = pickle.load(open(os.path.join(here, "..", "out", "gencode47.pkl"), "rb"))
can = {}
for t_, r_ in gc.items():
    if r_["tags"] & {"MANE_Select", "Ensembl_canonical"} and r_["ttype"] == "protein_coding":
        can.setdefault((r_["chrom"], r_["strand"]), []).append(r_["tss"])
can = {k: np.array(sorted(v)) for k, v in can.items()}
def canon(c, s, p):
    a = can.get((c, s))
    if a is None: return False
    i = np.searchsorted(a, p)
    return any(0 <= j < len(a) and abs(int(a[j]) - p) <= 50 for j in (i - 1, i))
# the rules of assemble.py, on the sites first: the full pair x site join of the pools is too large to hold
J = ["regime", "chrom", "strand", "s_chain", "s_tss_group_lo"]
fr = []
for s_ in ["BT474", "K562", "HG002", "UHRR"]:
    x = pd.read_csv(os.path.join(here, "..", "out", "scored", f"{s_}_pool.sites.tsv.gz"), sep="\t", low_memory=False)
    fl = x.n_site / x.n_local50.replace(0, np.nan); fi = x.n_site / x.n_first_intron_reads.replace(0, np.nan); cf = x.n_cap / x.n_inf.replace(0, np.nan)
    pile = (x.n_site >= 5) & (fl >= 0.5) & (fi >= 0.10) & (x.n_inf >= 5); capped = (x.n_cap >= 2) & (cf >= 0.10); atlas = (x.cage_dist <= 10) | (x.reftss_dist <= 10)
    x["group"] = np.where(pile & capped & atlas, "real", np.where(pile & ~capped & ~atlas, "truncation", ""))   # uncapped, no atlas: the 5' truncation pile-ups
    x = x[x.group != ""].copy()
    x["cap_frac"] = x.n_cap / x.n_inf
    x = x[[not canon(c, st, int(p)) for c, st, p in zip(x.chrom, x.strand, x.peak)]]
    keys = x[J]
    ps = [c.merge(keys, on=J) for c in pd.read_csv(os.path.join(here, "..", "out", "scored", f"{s_}_pool.pairs.tsv.gz"), sep="\t", low_memory=False, chunksize=500000, usecols=lambda c: c in set(J) | {"dataset", "class", "protein_class", "l_up_junc_reads"})]
    p = pd.concat(ps); p = p[(p["class"] == "strict_ISM") & (p.l_up_junc_reads >= 5)]
    p["protein_qualifying"] = p.protein_class.isin(["N_truncated", "N_truncated_after_uORF"])
    fr.append(p.merge(x.drop(columns=[c for c in x.columns if c in ("dataset", "class")]), on=J))
    print(s_, len(fr[-1]), flush=True)
d = pd.concat(fr, ignore_index=True)
d["qual"] = d.protein_qualifying.fillna(False).astype(bool)
# One target per site + first intron, not per chain: exact S chains would favour whichever tool contributed most fragment chains
# to the pool. A tool recovers a target with any model that uses S's first intron and starts within 50 nt of the site.
d["first"] = [tuple(int(v) for v in x.split("-")) for x in d.s_first_intron] if "s_first_intron" in d else [(parse_chain(c)[0] if st == "+" else parse_chain(c)[-1]) for c, st in zip(d.s_chain, d.strand)]
K = ["dataset", "chrom", "strand", "first", "peak"]
t = d.groupby(K, as_index=False).agg(group=("group", "first"), gene=("gene", "first"), n_site=("n_site", "max"), cap_frac=("cap_frac", "max"), qual=("qual", "max"), annot=("annot_same_start_intron_nonNF", "max"))
rows, per = [], []
for ds, g in t.groupby("dataset"):
    s = ds.replace("_pool", ""); idx = pickle.load(open(f"{pool}/{s}.tool_index.pkl", "rb"))
    g = g.copy()
    for tool, m in idx.items():
        fi = {}
        for (c, st, ch), v in m.items():
            fi.setdefault((c, st, ch[0] if st == "+" else ch[-1]), []).extend(t for t, _ in v)   # index holds (5' end, 3' end)
        hit = np.array([any(abs(x - pk) <= 50 for x in fi.get((c, st, f), [])) for c, st, f, pk in zip(g.chrom, g.strand, g["first"], g.peak)])
        g[tool] = hit
        real = (g.group == "real").values
        for tier, mask in (("internal TSS", real), ("N-truncated", real & g.qual.values.astype(bool)), ("internal TSS, not annotated", real & ~g.annot.values.astype(bool)), ("uncapped pile-up", ~real)):
            rows.append(dict(sample=s, tier=tier, tool=tool, targets=int(mask.sum()), recovered=int(hit[mask].sum())))
    per.append(g)
r = pd.DataFrame(rows); r["pct"] = 100 * r.recovered / r.targets
r.to_csv(os.path.join(here, "..", "out", "tool_recall.tsv"), sep="\t", index=False)
pd.concat(per).to_csv(os.path.join(here, "..", "out", "tool_recall.targets.tsv"), sep="\t", index=False)
a = r.groupby(["tier", "tool"])[["targets", "recovered"]].sum(); a["pct"] = (100 * a.recovered / a.targets).round(1)
print(a.to_string())
print(r.pivot_table(index=["tier", "tool"], columns="sample", values="pct").round(1).to_string())
