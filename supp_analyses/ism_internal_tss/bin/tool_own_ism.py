#!/usr/bin/env python3
"""Each tool's own ISM isoforms, judged by the reads and atlases (user, 2026-11-10; replaces the pooled-target recall panel,
whose targets were mostly IsoSeq's own fragments). For each bulk library and tool: the tool's de novo models S that are strict
ISMs of a longer model L of the SAME tool or of GENCODE v47 (S's chain a proper suffix of L's, S's 5' end inside the L exon at
S's first donor, >= 5 reads over L's dropped intron, S's start > 50 nt from a canonical protein-coding TSS). Pairs and site
evidence come from the pool tables (make_pool_gtf.py, score_sites.py), which hold every tool's models; a pool pair counts for a
tool when the tool has a model with S's chain starting in S's TSS group and ending within END3_TOL nt of the 3' end of
a model with L's chain from the same tool or GENCODE (shared 3' end, user 2026-11-11).
Each start (site + S's first intron) is classed by the same rules as tool_recall.py:
  supported         read pile-up (>= 5 informative reads) with the cap signature, CAGE or refTSS peak within 10 nt
  truncation-like   pile-up without the cap signature and with no atlas peak
  other             no pile-up, too few informative reads, or mixed evidence
Writes out/tool_own_ism.tsv (per library and tool) and out/tool_own_ism.pooled.tsv (four libraries summed).
usage: tool_own_ism.py"""
import pickle
import numpy as np, pandas as pd
from common import *
here = os.path.dirname(os.path.abspath(__file__)); pool = os.path.join(here, "..", "out", "pool")
TOOLS = ["LRAA", "IsoQuant", "StringTie", "Bambu", "ESPRESSO", "IsoSeq", "Mandalorion"]
gc = pickle.load(open(os.path.join(here, "..", "out", "gencode47.pkl"), "rb"))
gends = {}   # GENCODE chain -> 3' ends
for r in gc.values():
    if r["chain"]: gends.setdefault((r["chrom"], r["strand"], r["chain"]), []).append(r["exons"][-1][1] if r["strand"] == "+" else r["exons"][0][0])
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
J = ["regime", "chrom", "strand", "s_chain", "s_tss_group_lo"]
RANK = {"supported": 2, "truncation-like": 1, "other": 0}
rows = []
for s_ in ["BT474", "K562", "HG002", "UHRR"]:
    x = pd.read_csv(os.path.join(here, "..", "out", "scored", f"{s_}_pool.sites.tsv.gz"), sep="\t", low_memory=False,
                    usecols=J + ["peak", "n_site", "n_local50", "n_first_intron_reads", "n_inf", "n_cap", "cage_dist", "reftss_dist"])
    fl = x.n_site / x.n_local50.replace(0, np.nan); fi = x.n_site / x.n_first_intron_reads.replace(0, np.nan); cf = x.n_cap / x.n_inf.replace(0, np.nan)
    pile = (x.n_site >= 5) & (fl >= 0.5) & (fi >= 0.10) & (x.n_inf >= 5); capped = (x.n_cap >= 2) & (cf >= 0.10); atlas = (x.cage_dist <= 10) | (x.reftss_dist <= 10)
    x["group"] = np.where(pile & capped & atlas, "supported", np.where(pile & ~capped & ~atlas, "truncation-like", "other"))
    x = x[J + ["group"]]
    ps = []
    for c in pd.read_csv(os.path.join(here, "..", "out", "scored", f"{s_}_pool.pairs.tsv.gz"), sep="\t", low_memory=False, chunksize=500000,
                         usecols=J + ["class", "l_up_junc_reads", "l_chain", "s_tss_group_hi", "s_first_intron", "peak"]):
        ps.append(c[(c["class"] == "strict_ISM") & (c.l_up_junc_reads >= 5)])
    p = pd.concat(ps, ignore_index=True)
    p = p[[not canon(c, st, int(k)) for c, st, k in zip(p.chrom, p.strand, p.peak)]]
    p = p.merge(x, on=J, how="left"); p["group"] = p.group.fillna("other")
    sk = [(c, st, parse_chain(ch)) for c, st, ch in zip(p.chrom, p.strand, p.s_chain)]
    lk = [(c, st, parse_chain(ch)) for c, st, ch in zip(p.chrom, p.strand, p.l_chain)]
    idx = pickle.load(open(f"{pool}/{s_}.tool_index.pkl", "rb"))
    for tool in TOOLS:
        m = idx[tool]
        def pair_ok(k, lo, hi, l):
            le = [e for _, e in m.get(l, ())] + gends.get(l, [])
            return any(lo <= t <= hi and any(abs(e - x) <= END3_TOL for x in le) for t, e in m.get(k, ()))
        ok = np.array([pair_ok(k, lo, hi, l) for k, lo, hi, l in zip(sk, p.s_tss_group_lo, p.s_tss_group_hi, lk)], dtype=bool)
        q = p[ok].copy()
        q["rank"] = q.group.map(RANK)
        st = q.groupby(["chrom", "strand", "s_first_intron", "peak"])["rank"].max()   # one start = site + S's first intron
        n_iso = q.drop_duplicates(["chrom", "strand", "s_chain", "s_tss_group_lo"]).shape[0]
        rows.append(dict(sample=s_, tool=tool, ism_isoforms=n_iso, starts=len(st), supported=int((st == 2).sum()),
                         truncation_like=int((st == 1).sum()), other=int((st == 0).sum())))
        print(rows[-1], flush=True)
r = pd.DataFrame(rows)
r.to_csv(os.path.join(here, "..", "out", "tool_own_ism.tsv"), sep="\t", index=False)
a = r.groupby("tool", as_index=False)[["ism_isoforms", "starts", "supported", "truncation_like", "other"]].sum()
a["judged"] = a.supported + a.truncation_like
a["pct_truncation_of_judged"] = 100 * a.truncation_like / a.judged
a["pct_truncation_of_starts"] = 100 * a.truncation_like / a.starts
a["pct_supported_of_starts"] = 100 * a.supported / a.starts
a.to_csv(os.path.join(here, "..", "out", "tool_own_ism.pooled.tsv"), sep="\t", index=False)
print(a.to_string())
