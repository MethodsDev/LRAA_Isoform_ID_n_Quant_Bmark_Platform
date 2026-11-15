#!/usr/bin/env python3
"""Independent evidence per S TSS (one S chain + TSS group, per regime): reads, CAGE, refTSS, background,
annotation, differential usage. usage: score_sites.py DATASET OUT.tsv.gz"""
import sys, pickle, collections
import numpy as np, pandas as pd, pysam
from common import *
ds, out = sys.argv[1:3]
cfg = DATASETS[ds]
here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(47)

from score_atlas import cage, reftss

# annotation -------------------------------------------------------------------------------
gc = pickle.load(open(os.path.join(here, "..", "out", "gencode47.pkl"), "rb"))
by_chain, by_first, by_intron, tss_idx = {}, {}, {}, collections.defaultdict(list)
for t, r in gc.items():
    if r["strand"] not in "+-": continue
    tss_idx[(r["chrom"], r["strand"])].append((r["tss"], t))
    if not r["chain"]: continue
    by_chain.setdefault((r["chrom"], r["strand"], r["chain"]), []).append(t)
    f = r["chain"][0] if r["strand"] == "+" else r["chain"][-1]
    by_first.setdefault((r["chrom"], r["strand"], f), []).append(t)
    for i in r["chain"]:
        by_intron.setdefault((r["chrom"], r["strand"], i), set()).add(r["gene"])
tss_arr = {k: (np.array([a for a, b in sorted(v)]), [b for a, b in sorted(v)]) for k, v in tss_idx.items()}
def nf(t): return "mRNA_start_NF" in gc[t]["tags"]
def near_tss(c, s, p, w=50, nonNF=True):
    a, ids = tss_arr.get((c, s), (np.array([], int), []))
    lo, hi = np.searchsorted(a, p - w), np.searchsorted(a, p + w, side="right")
    r = [(abs(int(a[i]) - p), ids[i]) for i in range(lo, hi) if not (nonNF and nf(ids[i]))]
    return min(r) if r else (None, None)
def gene_of(c, s, chain):
    g = collections.Counter(x for i in chain for x in by_intron.get((c, s, i), ()))
    return g.most_common(1)[0][0] if g else ""

# reads ------------------------------------------------------------------------------------
sc = "chrom strand fi pos cluster full n inf cap".split()
st = pd.read_csv(os.path.join(here, "..", "out", "scan", f"{ds}.starts.tsv.gz"), sep="\t", header=None, names=sc, dtype={"cluster": str})
print("start rows", len(st), "reads", st.n.sum(), flush=True)
groups = {k: v.index.values for k, v in st.groupby(["chrom", "strand", "fi"], sort=False)}
ps, ns, infs, caps, fulls, cls = st.pos.values, st.n.values, st.inf.values, st.cap.values, st.full.values, st.cluster.values

fa = pysam.FastaFile(cfg["genome"])
def sense(c, s, p, up=3, dn=3):
    if s == "+":
        return fa.fetch(c, p - 1 - up, p - 1 + dn).upper()
    x = fa.fetch(c, p - dn, p + up).upper()
    return x[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

# fragile cell types: PBMC monocytes, ALS neurons
fragile = set()
if ds == "PBMC":
    d = pd.read_csv(f"{P}/pbmcs_refGuided/data/RefQuantOnly.cell_clusters_and_cell_types.wUMAP.tsv", sep="\t", usecols=["seurat_clusters", "cell_type_annot"])
    fragile = {str(k) for k, v in d.groupby("seurat_clusters").cell_type_annot.agg(lambda x: x.mode().iat[0]).items() if "monocyte" in v.lower()}
elif ds == "ALS":
    d = pd.read_csv(f"{A}/als_refGuided/data/RefQuantOnly.LP00714.ALS.cell_clusters_and_cell_types.wUMAP.tsv", sep="\t", usecols=["seurat_clusters", "cell_type_simplified"])
    fragile = {str(k) for k, v in d.groupby("seurat_clusters").cell_type_simplified.agg(lambda x: x.mode().iat[0]).items() if "neuron" in v.lower()}
print("fragile clusters", sorted(fragile), flush=True)

# differential usage ----------------------------------------------------------------------
du = None
if cfg["site_prefix"]:
    S = pd.read_csv(cfg["site_prefix"] + ".sites.tsv", sep="\t", usecols=["site_id", "kind", "chrom", "strand", "pos", "support"])
    S = S[S.kind == "TSS"]
    E_ = pd.read_csv(cfg["site_prefix"] + ".dexseq.TSS.events.tsv", sep="\t")
    eidx = collections.defaultdict(list)
    for r in E_.itertuples():
        eidx[r.gained_site].append((r, "gained")); eidx[r.lost_site].append((r, "lost"))
    site_arr = {k: (v.pos.values, v.site_id.values) for k, v in ((k, v.sort_values("pos")) for k, v in S.groupby(["chrom", "strand"]))}
    du = (site_arr, eidx)
def du_of(c, s, p):
    if du is None: return ("", 0, 0, "")
    a = du[0].get((c, s))
    if a is None: return ("", 0, 0, "")
    lo, hi = np.searchsorted(a[0], p - 25), np.searchsorted(a[0], p + 25, side="right")
    if hi <= lo: return ("", 0, 0, "")
    sid = a[1][lo + int(np.argmin(np.abs(a[0][lo:hi] - p)))]
    ev = du[1].get(sid, [])
    if not ev: return (sid, 0, 0, "")
    hc = [e for e in ev if e[0].high_confidence]
    b = max(hc or ev, key=lambda e: e[0].abs_delta)
    r = b[0]
    return (sid, len(ev), len(hc), f"{r.cluster_A} vs {r.cluster_B}: site {b[1]}, delta {r.abs_delta:.2f}, {r.event_type}")

# sites ------------------------------------------------------------------------------------
p = pd.read_csv(os.path.join(here, "..", "out", "pairs", f"{ds}.pairs.tsv.gz"), sep="\t")
# L with most reads over the dropped junction gives the site's exon for the background
jc = pd.read_csv(os.path.join(here, "..", "out", "scan", f"{ds}.junc.tsv.gz"), sep="\t", header=None, names=["chrom", "strand", "l_up_intron", "l_up_junc_reads"])
p = p.merge(jc, on=["chrom", "strand", "l_up_intron"], how="left").fillna({"l_up_junc_reads": 0})
p["l_up_junc_reads"] = p.l_up_junc_reads.astype(int)
key = ["regime", "chrom", "strand", "s_chain", "s_tss_group_lo"]
rows = []
for k, g in p.sort_values("l_up_junc_reads", ascending=False).groupby(key, sort=False):
    regime, c, s, chain_s, glo = k
    e = g.iloc[0]
    tss, ghi = int(e.s_tss), int(e.s_tss_group_hi)
    sc_ = parse_chain(chain_s)
    idx = groups.get((c, s, e.s_first_intron), np.array([], int))
    P_, N_ = ps[idx], ns[idx]
    w = (P_ >= glo - 5) & (P_ <= ghi + 5)
    peak = tss
    if w.any():
        q = collections.Counter()
        for a, b in zip(P_[w], N_[w]): q[int(a)] += int(b)
        peak = min(q, key=lambda x: (-q[x], abs(x - tss)))
    sw = np.abs(P_ - peak) <= 5
    n_site, n_inf, n_cap = int(N_[sw].sum()), int(infs[idx][sw].sum()), int(caps[idx][sw].sum())
    away = np.abs(P_ - peak) > 50
    n_fi, n_local = int(N_.sum()), int(N_[np.abs(P_ - peak) <= 50].sum())
    n_ctl_inf, n_ctl_cap = int(infs[idx][away].sum()), int(caps[idx][away].sum())
    n_full = int(N_[sw & (fulls[idx] == 1)].sum())
    n_core = int(N_[np.abs(P_ - peak) <= 2].sum())
    cc = collections.Counter()
    for a, b in zip(cls[idx][sw], N_[sw]): cc[a] += int(b)
    frag = sum(v for a, v in cc.items() if a in fragile)
    clustered = sum(v for a, v in cc.items() if a != "NA")
    # atlases at the peak, background inside L's exon
    dc, cid, csc = cage.dist(c, s, np.array([peak])); dr, rid, _ = reftss.dist(c, s, np.array([peak]))
    dct, _, _ = cage.dist(c, s, np.array([tss]))
    bg = {}
    if e["class"] == "strict_ISM":
        lo, hi = int(e.l_exon_lo), int(e.l_exon_hi)
        pool = np.arange(lo + 10, hi - 9)
        pool = pool[np.abs(pool - peak) > 50]
        if len(pool):
            r_ = rng.choice(pool, size=10)
            d1, _, _ = cage.dist(c, s, r_); d2, _, _ = reftss.dist(c, s, r_)
            bg = dict(bg_cage_d0=np.mean(d1 == 0), bg_cage_d10=np.mean(d1 <= 10), bg_cage_d25=np.mean(d1 <= 25), bg_cage_d50=np.mean(d1 <= 50), bg_reftss_d10=np.mean(d2 <= 10))
    # annotation
    same = by_chain.get((c, s, sc_), [])
    same_nonNF = [t for t in same if not nf(t)]
    d_all, t_all = near_tss(c, s, peak, 50, nonNF=False); d_n, t_n = near_tss(c, s, peak, 50, nonNF=True)
    fi = tuple(map(int, e.s_first_intron.split("-")))
    sh = [t for t in by_first.get((c, s, fi), []) if abs(gc[t]["tss"] - peak) <= 50]
    sh_n = [t for t in sh if not nf(t)]
    gene = gene_of(c, s, sc_) or gene_of(c, s, parse_chain(e.l_chain)) or (e.lraa_symbol if isinstance(e.lraa_symbol, str) else "") or "novel"
    best_ann = (sh_n or sh or [None])[0]
    sid, n_ev, n_hc, du_txt = du_of(c, s, peak)
    rows.append(dict(zip(key, k)) | dict(class_site=e["class"], peak=peak, gene=gene, n_L=g.l_chain.nunique(),
        n_first_intron_reads=n_fi, n_site=n_site, n_site_core2=n_core, n_local50=n_local, n_full_S_chain=n_full, n_inf=n_inf, n_cap=n_cap,
        ctl_inf=n_ctl_inf, ctl_cap=n_ctl_cap, seq_m3_p3=sense(c, s, peak), fragile_share=(frag / clustered if clustered else np.nan),
        top_cluster=(max(cc, key=cc.get) if cc else ""), clusters=";".join(f"{a}:{v}" for a, v in cc.most_common(6)),
        cage_dist=int(dc[0]), cage_id=cid[0], cage_score=float(csc[0]), cage_dist_modelTSS=int(dct[0]), reftss_dist=int(dr[0]), reftss_id=rid[0],
        annot_same_chain=len(same) > 0, annot_same_chain_nonNF=len(same_nonNF) > 0,
        annot_tss_dist_any=d_all, annot_tss_dist_nonNF=d_n, annot_same_start_intron=len(sh) > 0, annot_same_start_intron_nonNF=len(sh_n) > 0,
        annot_tx=best_ann or (t_n or ""), annot_ttype=gc[best_ann or t_n]["ttype"] if (best_ann or t_n) else "",
        annot_tags=",".join(sorted(gc[best_ann or t_n]["tags"] & {"TAGENE", "MANE_Select", "Ensembl_canonical", "mRNA_start_NF", "basic"})) if (best_ann or t_n) else "",
        du_site=sid, du_events=n_ev, du_HC_events=n_hc, du_best=du_txt) | bg)
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
p.to_csv(os.path.join(here, "..", "out", "pairs", f"{ds}.pairs_j.tsv.gz"), sep="\t", index=False)
print("sites", len(rows))
