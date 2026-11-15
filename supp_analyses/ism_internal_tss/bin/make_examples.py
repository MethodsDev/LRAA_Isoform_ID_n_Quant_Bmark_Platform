#!/usr/bin/env python3
"""The 10 best S/L pairs and everything their plots need. Rule: qualifying strict ISM (core tier, start lost,
N-terminal truncation, S not at a canonical TSS, >= 5 reads over L's dropped junction), L matches a GENCODE
coding chain, >= 20 reads at the site in the dataset drawn, one per gene. Ranked by evidence points, then
single-cell datasets first, then cap fraction. Drawn from the dataset with most reads at the site. Also draws the FTD
example chosen by select_ftd_example.py (example "ftd"). CDS: L's (GENCODE or pyTD) and S's (L's from S's ATG).
usage: make_examples.py FULL_TSV"""
import sys, gzip, pickle
import numpy as np, pandas as pd
from common import *
sys.path.insert(0, os.path.join(LRAA_HOME, "pylib")); sys.path.insert(0, os.path.join(LRAA_HOME, "util/sc/diff_iso_usage"))
import pysam, Util_funcs, TSS_cap_calibration as T
from alt_termini_read_check import make_read_filter
from protein_lib import s_exons
from score_atlas import cage as CA, reftss as RT
here = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(here, "..", "out"); E = os.path.join(O, "examples"); os.makedirs(E, exist_ok=True)
d = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
d["n_points"] = (d.ev_cage.astype(int) + d.ev_reftss + d.ev_cap_strong + d.ev_pile + d.ev_annot_tx + d.annot_S_same_ATG + d.ev_du_HC + (d.n_core_datasets >= 2).astype(int) + (d.n_core_datasets >= 3))
q = d[d.core & d.protein_qualifying & ~d.s_at_canonical_tss & (d.l_up_junc_reads >= 5) & (d["class"] == "strict_ISM") & d.primary_L & (d.cds_route == "GENCODE")].copy()
q["sc"] = q.dataset.isin(["PBMC", "ALS"])
best = q.sort_values(["n_site"], ascending=False).drop_duplicates(["gene", "chrom", "strand"])   # dataset with most reads at the gene's site
rank = q.groupby("gene").agg(n_points=("n_points", "max"), sc=("sc", "max")).reset_index()
best = best.drop(columns=["n_points", "sc"]).merge(rank, on="gene")
best = best[best.n_site >= 20].sort_values(["n_points", "sc", "cap_frac"], ascending=False)
ex = best.head(10).reset_index(drop=True); ex["example"] = [f"ex{i+1:02d}" for i in range(len(ex))]
best.iloc[10:40].to_csv(f"{O}/next_candidates.tsv", sep="\t", index=False)
# the FTD example (select_ftd_example.py's rank 1, drawn as example "ftd"), any CDS route
fc = pd.read_csv(f"{E}/ftd_example_candidates.tsv", sep="\t", comment="#")
fc = fc[fc.passes].head(1)
if len(fc):
    f0 = fc.iloc[0]
    fr = d[(d.dataset == "ALS") & (d.gene == f0.gene) & (d.regime == f0.regime) & (d.arms == f0.arms) & (d.peak == f0.peak) & (d.s_hash == f0.s_hash) & (d.l_hash == f0.l_hash)].head(1).copy()
    fr["example"] = "ftd"; fr["sc"] = True
    ex = pd.concat([ex, fr], ignore_index=True)
print(ex[["example", "gene", "dataset", "regime", "peak", "n_points", "n_site", "cap_frac", "protein_class", "aa_lost", "aa_total", "kozak"]].to_string(), flush=True)
gc = pickle.load(open(os.path.join(O, "gencode47.pkl"), "rb"))
idx = {}
for t, r in gc.items():
    if r["chain"]: idx.setdefault((r["chrom"], r["strand"], r["chain"]), []).append(t)
orfs = {}
for r in pd.read_csv(f"{O}/orf/orfs.tsv", sep="\t").itertuples():
    orfs[(r.chrom, r.strand, r.l_chain)] = [tuple(map(int, x.split("-"))) for x in r.cds.split(",")]
M, C, R, H = [], [], [], []
keep_f = make_read_filter(HiFi=True, rdna_mask_bed=RDNA)
scan = {}
for e in ex.itertuples():
    cfg = DATASETS[e.dataset]
    lc, sc = parse_chain(e.l_chain), parse_chain(e.s_chain)
    sx = s_exons(sc, int(e.peak), int(e.s_end3), e.strand); lx = s_exons(lc, int(e.l_tss), int(e.s_end3), e.strand)
    for a, b in sx: M.append((e.example, "S", "S", a, b, "exon"))
    for a, b in lx: M.append((e.example, "L", "L", a, b, "exon"))
    # L's CDS (GENCODE, or the pyTD ORF; the latter includes the stop codon) and S's: L's CDS from S's ATG on
    lcds = gc[e.L_gencode_tx]["cds"] if e.cds_route == "GENCODE" else orfs.get((e.chrom, e.strand, e.l_chain), [])
    for a, b in lcds: M.append((e.example, "L", "L", a, b, "cds"))
    if lcds: M.append((e.example, "ATG", "L ATG", (min(a for a, b in lcds) if e.strand == "+" else max(b for a, b in lcds)), 0, "site"))
    if pd.notna(e.new_atg_gpos):
        g0 = int(e.new_atg_gpos); M.append((e.example, "ATG", "S ATG", g0, 0, "site"))
        for a, b in lcds:
            a2, b2 = (max(a, g0), b) if e.strand == "+" else (a, min(b, g0))
            if a2 <= b2: M.append((e.example, "S", "S", a2, b2, "cds"))
    ts = sorted(idx.get((e.chrom, e.strand, sc), []), key=lambda t: "mRNA_start_NF" in gc[t]["tags"])[:1]
    for t in ts:
        for a, b in gc[t]["exons"]: M.append((e.example, "GENCODE S chain", t, a, b, "exon"))
    a0, b0 = int(e.peak) - 1500, int(e.peak) + 1500
    for nm, at in (("CAGE", CA), ("refTSS", RT)):
        s_, e_, ids, sc_ = at.g[(e.chrom, e.strand)]
        for i in np.where((e_ >= a0) & (s_ <= b0))[0]:
            C.append((e.example, nm, s_[i] + 1, e_[i], ids[i], sc_[i]))
    # per-position starts over S's first intron (scan table)
    if e.dataset not in scan:
        scan[e.dataset] = pd.read_csv(f"{O}/scan/{e.dataset}.starts.tsv.gz", sep="\t", header=None, names="chrom strand fi pos cluster full n inf cap".split(), dtype={"cluster": str})
    t = scan[e.dataset]; t = t[(t.chrom == e.chrom) & (t.strand == e.strand) & (t.fi == e.s_first_intron)]
    h = t.groupby("pos")[["n", "inf", "cap"]].sum().reset_index(); h["example"] = e.example; H.append(h)
    # reads across S's first exon
    fa = pysam.FastaFile(cfg["genome"]); seq = fa.fetch(e.chrom); Util_funcs.register_contig_seq_for_strand_check(e.chrom, seq.upper())
    bam = pysam.AlignmentFile(cfg["bam"])
    fi = tuple(map(int, e.s_first_intron.split("-"))); up = tuple(map(int, e.l_up_intron.split("-")))
    donor = fi[0] if e.strand == "+" else fi[1]
    n = 0; seen = set()
    for r in bam.fetch(e.chrom, donor - 2, donor + 1):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not keep_f(r): continue
        if Util_funcs.transcribed_strand(r) != e.strand: continue
        x = T._five_prime_start(r, e.strand, seq)
        if x is None: continue
        p5, G, _ = x
        pos = r.reference_start; introns = []
        for o, ln in r.cigartuples:
            if o == 3: introns.append((pos + 1, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        cls = "through L intron" if up in introns else ("first intron = S's" if introns and (introns[0] if e.strand == "+" else introns[-1]) == fi else "other")
        cap = "NA" if G is None else str(G > cfg["baseline"])
        R.append((e.example, p5, cap, cls, abs(p5 - int(e.peak)) <= 5, ";".join(f"{a+1}-{b}" for a, b in r.get_blocks()), True))
        seen.add((r.query_name, r.reference_start)); n += 1
    # also the reads starting at L's TSS (within 50 nt) that end before S's first donor, so L's 5' support shows
    # (2026-11-08); via_s_donor False keeps them out of the read-end histogram, which counts reads crossing that donor
    lt = int(e.l_tss); n_l = 0
    for r in bam.fetch(e.chrom, lt - 60, lt + 60):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not keep_f(r) or (r.query_name, r.reference_start) in seen: continue
        if Util_funcs.transcribed_strand(r) != e.strand: continue
        x = T._five_prime_start(r, e.strand, seq)
        if x is None or abs(x[0] - lt) > 50: continue
        p5, G, _ = x
        cap = "NA" if G is None else str(G > cfg["baseline"])
        R.append((e.example, p5, cap, "starts at L TSS", False, ";".join(f"{a+1}-{b}" for a, b in r.get_blocks()), False))
        n_l += 1
    print(e.example, e.gene, n, "reads over S's donor,", n_l, "more starting at L's TSS", flush=True)
ex.to_csv(f"{E}/examples.tsv", sep="\t", index=False)
pd.DataFrame(M, columns="example track label start end kind".split()).to_csv(f"{E}/models.tsv", sep="\t", index=False)
pd.DataFrame(C, columns="example atlas start end id score".split()).to_csv(f"{E}/atlas.tsv", sep="\t", index=False)
pd.DataFrame(R, columns="example start5 capped rclass at_site blocks via_s_donor".split()).to_csv(f"{E}/reads.tsv.gz", sep="\t", index=False)
pd.concat(H).to_csv(f"{E}/starts.tsv", sep="\t", index=False)
# read 5' ends over the locus, per position and read class (aggregate counts; drawn under each model panel)
rr = pd.DataFrame(R, columns="example start5 capped rclass at_site blocks via_s_donor".split())
rr = rr[rr.via_s_donor]
rr["read_class"] = np.select([rr.rclass == "through L intron", rr.at_site & (rr.capped == "True"), rr.at_site & (rr.capped == "False"), rr.at_site],
                             ["L chain (through L intron)", "S site, capped", "S site, not capped", "S site, no cap information"], "starts elsewhere")
rr.groupby(["example", "start5", "read_class"]).size().rename("reads").reset_index().to_csv(f"{E}/read_starts_locus.tsv", sep="\t", index=False)
