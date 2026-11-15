#!/usr/bin/env python3
"""candidates.tsv (one row per S/L pair) and the tier flags. usage: assemble.py OUT.tsv DATASET... """
import sys, pickle
import numpy as np, pandas as pd
from common import *
out, dss = sys.argv[1], sys.argv[2:]
here = os.path.dirname(os.path.abspath(__file__))
gc = pickle.load(open(os.path.join(here, "..", "out", "gencode47.pkl"), "rb"))
can = {}
for t, r in gc.items():
    if r["tags"] & {"MANE_Select", "Ensembl_canonical"} and r["ttype"] == "protein_coding":
        can.setdefault((r["chrom"], r["strand"]), []).append(r["tss"])
can = {k: np.array(sorted(v)) for k, v in can.items()}
def canon_dist(c, s, p):
    a = can.get((c, s))
    if a is None: return 10**9
    i = np.searchsorted(a, p)
    return min([abs(int(a[j]) - p) for j in (i - 1, i) if 0 <= j < len(a)] or [10**9])
K = ["regime", "chrom", "strand", "s_chain", "s_tss_group_lo"]
fr = []
for ds in dss:
    p = pd.read_csv(os.path.join(here, "..", "out", "scored", f"{ds}.pairs.tsv.gz"), sep="\t", low_memory=False)
    s = pd.read_csv(os.path.join(here, "..", "out", "scored", f"{ds}.sites.tsv.gz"), sep="\t", low_memory=False).drop(columns=["peak"])
    x = p.merge(s, on=K)
    fr.append(x)
d = pd.concat(fr, ignore_index=True)
d["cap_frac"] = d.n_cap / d.n_inf.replace(0, np.nan)
d["ctl_cap_frac"] = d.ctl_cap / d.ctl_inf.replace(0, np.nan)
d["site_frac_of_intron_reads"] = d.n_site / d.n_first_intron_reads.replace(0, np.nan)
d["site_frac_local50"] = d.n_site / d.n_local50.replace(0, np.nan)
d["site_share_vs_L"] = d.n_site / (d.n_site + d.l_up_junc_reads).replace(0, np.nan)
d["ev_pile"] = (d.n_site >= 5) & (d.site_frac_local50 >= 0.5) & (d.site_frac_of_intron_reads >= 0.10)
d["ev_cap_measurable"] = d.n_inf >= 5
d["ev_cap"] = d.ev_cap_measurable & (d.n_cap >= 2) & (d.cap_frac >= 0.10)
d["ev_cap_strong"] = d.ev_cap & (d.cap_frac >= 0.30)
d["ev_cage"] = d.cage_dist <= 10
d["ev_reftss"] = d.reftss_dist <= 10
d["ev_atlas"] = d.ev_cage | d.ev_reftss
d["ev_annot_tss"] = d.annot_tss_dist_nonNF.notna()
d["ev_annot_tx"] = d.annot_same_start_intron_nonNF
d["ev_du"] = d.du_events > 0
d["ev_du_HC"] = d.du_HC_events > 0
tripled = d.drop_duplicates(["chrom", "strand", "peak"])[["chrom", "strand", "peak"]].copy()
tripled["canon_tss_dist"] = [canon_dist(c, s, int(p)) for c, s, p in zip(tripled.chrom, tripled.strand, tripled.peak)]
d = d.merge(tripled, on=["chrom", "strand", "peak"])
d["s_at_canonical_tss"] = d.canon_tss_dist <= 50
gcs = {}
for t, r in gc.items():
    if r["start"]:
        gcs.setdefault((r["chrom"], r["strand"], r["chain"]), []).append(r["start"][0][0] if r["strand"] == "+" else r["start"][-1][1])
def same_start(c, s, sc, a):
    return bool(pd.notna(a)) and int(a) in gcs.get((c, s, parse_chain(sc)), [])
d["annot_S_same_ATG"] = [same_start(c, s, x, a) for c, s, x, a in zip(d.chrom, d.strand, d.s_chain, d.new_atg_gpos)]
d["protein_qualifying"] = d.protein_class.isin(["N_truncated", "N_truncated_after_uORF"])
d["protein_altered"] = d.protein_qualifying | d.protein_class.isin(["frame_changed", "N_altered_in_frame"])
# primary L of each S site: annotated coding L first (MANE, then more dropped-junction reads)
d["_o"] = d.cds_route.eq("GENCODE").astype(int) * 2 + d.cds_route.eq("pyTD").astype(int)
d = d.sort_values(["_o", "L_mane", "l_up_junc_reads"], ascending=False)
d["primary_L"] = ~d.duplicated(["dataset"] + K)
d = d.drop(columns="_o")
d["core"] = d.ev_pile & d.ev_cap & d.ev_atlas
# replication: datasets where a core site lies within 10 nt on the strand
cr = d[d.core & (d["class"] == "strict_ISM")].drop_duplicates(["dataset", "chrom", "strand", "peak"])[["dataset", "chrom", "strand", "peak"]]
rep = {}
for (c, s), g in cr.groupby(["chrom", "strand"]):
    pk = g.peak.values; dd = g.dataset.values
    for p_ in set(pk):
        rep[(c, s, int(p_))] = ";".join(sorted(set(dd[np.abs(pk - p_) <= 10])))
d["core_datasets"] = [rep.get((c, s, int(p_)), "") for c, s, p_ in zip(d.chrom, d.strand, d.peak)]
d["n_core_datasets"] = [len(x.split(";")) if x else 0 for x in d.core_datasets]
d.to_csv(out, sep="\t", index=False)
print(d.shape)
