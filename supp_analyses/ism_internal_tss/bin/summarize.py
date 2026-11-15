#!/usr/bin/env python3
"""Tables for the report, examples selection. usage: summarize.py FULL_TSV  (writes out/*.tsv, ../candidates.tsv)"""
import sys, pickle
import numpy as np, pandas as pd
from common import *
here = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(here, "..", "out")
d = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
order = ["PBMC", "ALS", "BT474", "K562", "HG002", "UHRR"]
d["n_points"] = (d.ev_cage.astype(int) + d.ev_reftss + d.ev_cap_strong + d.ev_pile + d.ev_annot_tx + d.annot_S_same_ATG + d.ev_du_HC
                 + (d.n_core_datasets >= 2).astype(int) + (d.n_core_datasets >= 3))
d["qualifying_strict"] = d.core & d.protein_qualifying & ~d.s_at_canonical_tss & (d.l_up_junc_reads >= 5) & (d["class"] == "strict_ISM")
d["altered_strict"] = d.core & d.protein_altered & ~d.s_at_canonical_tss & (d.l_up_junc_reads >= 5) & (d["class"] == "strict_ISM")
drop = ["s_chain", "l_chain", "bg_cage_d0", "bg_cage_d10", "bg_cage_d25", "bg_cage_d50", "bg_reftss_d10"]
cols = [c for c in d.columns if c not in drop]
first = ["dataset", "regime", "gene", "class", "chrom", "strand", "peak", "s_tss", "s_hash", "l_hash", "primary_L", "core", "protein_class", "qualifying_strict", "n_points"]
d[first + [c for c in cols if c not in first]].to_csv(os.path.join(here, "..", "candidates.tsv"), sep="\t", index=False)

# tier counts: pairs and unique S TSSs (S chain + TSS)
tiers = [("all", lambda x: x.dataset == x.dataset), ("pile-up", lambda x: x.ev_pile), ("cap", lambda x: x.ev_cap), ("CAGE", lambda x: x.ev_cage),
         ("refTSS", lambda x: x.ev_reftss), ("annotated start+intron", lambda x: x.ev_annot_tx), ("DU event", lambda x: x.ev_du), ("DU high confidence", lambda x: x.ev_du_HC),
         ("pile-up+cap", lambda x: x.ev_pile & x.ev_cap), ("core: pile-up+cap+atlas", lambda x: x.core),
         ("core, not canonical TSS", lambda x: x.core & ~x.s_at_canonical_tss),
         ("core, protein altered", lambda x: x.core & x.protein_altered & ~x.s_at_canonical_tss & (x.l_up_junc_reads >= 5)),
         ("core, N-truncated", lambda x: x.core & x.protein_qualifying & ~x.s_at_canonical_tss & (x.l_up_junc_reads >= 5))]
rows = []
for (ds, reg, cl), g in d.groupby(["dataset", "regime", "class"]):
    u = g[g.primary_L]
    for name, f in tiers:
        m = f(g)
        mu = f(u)
        rows.append(dict(dataset=ds, regime=reg, **{"class": cl}, tier=name, pairs=int(m.sum()), sites=int(mu.sum()), peaks=u[mu].drop_duplicates(["chrom", "strand", "peak"]).shape[0]))
pd.DataFrame(rows).to_csv(f"{O}/tier_counts.tsv", sep="\t", index=False)
# dataset-level, union of regimes
rows = []
for (ds, cl), g in d.groupby(["dataset", "class"]):
    u = g[g.primary_L]
    for name, f in tiers:
        rows.append(dict(dataset=ds, **{"class": cl}, tier=name, pairs=int(f(g).sum()), peaks=u[f(u)].drop_duplicates(["chrom", "strand", "peak"]).shape[0],
                         genes=u[f(u)].gene.nunique()))
pd.DataFrame(rows).to_csv(f"{O}/tier_counts_dataset.tsv", sep="\t", index=False)

# atlas rate vs background, strict ISM sites (bg only exists there)
s = d[d.primary_L & (d["class"] == "strict_ISM")].copy()
st = [("all S TSSs", s.dataset == s.dataset), ("pile-up, capped", s.ev_pile & s.ev_cap), ("pile-up, measurable, not capped", s.ev_pile & s.ev_cap_measurable & ~s.ev_cap),
      ("pile-up, cap not measurable", s.ev_pile & ~s.ev_cap_measurable), ("no pile-up, capped", ~s.ev_pile & s.ev_cap), ("no pile-up, not capped", ~s.ev_pile & ~s.ev_cap)]
rows = []
for ds, g in s.groupby("dataset"):
    for name, m in st:
        x = g[m.loc[g.index]].drop_duplicates(["chrom", "strand", "peak"])
        if not len(x): continue
        b = x[x.bg_cage_d10.notna()]
        for w, col in ((0, "bg_cage_d0"), (10, "bg_cage_d10"), (25, "bg_cage_d25"), (50, "bg_cage_d50")):
            rows.append(dict(dataset=ds, stratum=name, atlas="CAGE", window=w, n=len(x), obs=(b.cage_dist <= w).mean() if len(b) else np.nan, bg=b[col].mean() if len(b) else np.nan, n_bg=len(b)))
        rows.append(dict(dataset=ds, stratum=name, atlas="refTSS", window=10, n=len(x), obs=(b.reftss_dist <= 10).mean() if len(b) else np.nan, bg=b.bg_reftss_d10.mean() if len(b) else np.nan, n_bg=len(b)))
pd.DataFrame(rows).to_csv(f"{O}/atlas_background.tsv", sep="\t", index=False)
keep = ["dataset", "regime", "class", "gene", "chrom", "strand", "peak", "tss_into_exon", "tss_to_donor", "lraa_TSS_flag", "n_site", "n_inf", "n_cap", "cap_frac", "site_frac_local50", "site_frac_of_intron_reads", "cage_dist", "reftss_dist",
        "ev_pile", "ev_cap", "ev_cap_measurable", "ev_cage", "ev_atlas", "ev_annot_tx", "core", "fragile_share", "protein_class", "s_at_canonical_tss", "l_up_junc_reads", "cds_route", "aa_lost", "aa_remaining", "aa_total", "kozak"]
d[d.primary_L][keep].to_csv(f"{O}/sites.tsv.gz", sep="\t", index=False)

# off-peak capped rate: reads >50 nt from any CAGE peak, read population of the scan
import score_atlas as sa
rows = []
for ds in order:
    st = pd.read_csv(f"{O}/scan/{ds}.starts.tsv.gz", sep="\t", header=None, names="chrom strand fi pos cluster full n inf cap".split(), dtype={"cluster": str})
    far = np.zeros(len(st), bool)
    for (c, s_), g in st.groupby(["chrom", "strand"]):
        far[st.index.get_indexer(g.index)] = sa.cage.dist(c, s_, g.pos.values)[0] > 50
    rows.append(dict(dataset=ds, reads=int(st.n.sum()), inf=int(st.inf.sum()), cap=int(st.cap.sum()), off_inf=int(st.inf[far].sum()), off_cap=int(st.cap[far].sum()), on_inf=int(st.inf[~far].sum()), on_cap=int(st.cap[~far].sum())))
pd.DataFrame(rows).to_csv(f"{O}/cap_context.tsv", sep="\t", index=False)

# examples: best row per gene
q = d[d.qualifying_strict & d.primary_L].copy()
q["_route"] = q.cds_route.eq("GENCODE").astype(int)
q["_sc"] = q.dataset.isin(["PBMC", "ALS"]).astype(int)

q = q.sort_values(["n_points", "_route", "L_mane", "cap_frac", "n_site"], ascending=False).drop_duplicates(["gene"])
# draw from the dataset with most reads at that site, sc first when it has >= 20
q.head(60).to_csv(f"{O}/top_candidates.tsv", sep="\t", index=False)
print(q[["gene", "dataset", "regime", "peak", "protein_class", "cds_route", "aa_total", "aa_lost", "kozak", "n_site", "cap_frac", "n_points", "core_datasets", "L_mane", "ev_annot_tx", "annot_S_same_ATG", "du_HC_events", "fragile_share"]].head(30).to_string())
qq = d[d.qualifying_strict]
pd.DataFrame([dict(pairs=len(qq), peaks=qq.drop_duplicates(["chrom", "strand", "peak"]).shape[0], genes=qq.gene.nunique())]).to_csv(f"{O}/qualifying_summary.tsv", sep="\t", index=False)
print(d.qualifying_strict.sum(), d[d.qualifying_strict].drop_duplicates(["chrom", "strand", "peak"]).shape[0], d[d.qualifying_strict].gene.nunique())
