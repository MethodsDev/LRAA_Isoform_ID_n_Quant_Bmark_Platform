#!/usr/bin/env python3
"""The FTD (ALS LP00714) example of the Supp Fig (user, 2026-11-04: replaces MBP). Rule, applied to this version's table:
  1. an FTD internal TSS in the core set (pile-up + cap signature + CAGE/refTSS site), strict ISM, not a canonical TSS,
     >= 5 reads over L's dropped intron (the study's qualifying rule);
  2. S's protein N-terminally truncated against L (protein_class N_truncated or N_truncated_after_uORF);
  3. both the long (L) and the short (S) isoform are LRAA models of the same model set, and each has full-chain (FSM)
     read support: at least one read of the FTD BAM (LRAA's HiFi read filter) whose intron chain equals the model's
     whole chain (the SQLIKE read-mode definition);
  4. (user, 2026-11-08, after BTBD8, whose L 5' end had 2 such reads) L's CDS is a GENCODE coding transcript, as for the
     other examples, and L's 5' end is read-supported: >= 10 reads start within 50 nt of L's TSS and fit L: their
     introns are L's first introns, in order, and their 3' end lies in the matching L exon (a read ending in an intron
     belongs to another transcript);
  rank by reads at the internal site (n_site), then by the smaller of the two FSM read counts. One row per gene: its
  L/S pair with the most FSM support for both (then a GENCODE L CDS first, then the reference-guided run).
Writes out/examples/ftd_example_candidates.tsv (ranked, aggregate counts only) and prints the top rows.
usage: select_ftd_example.py FULL_TSV"""
import sys
import pandas as pd
from common import *
sys.path.insert(0, os.path.join(LRAA_HOME, "pylib")); sys.path.insert(0, os.path.join(LRAA_HOME, "util/sc/diff_iso_usage"))
import pysam, Util_funcs, TSS_cap_calibration as T
from alt_termini_read_check import make_read_filter
from protein_lib import s_exons
here = os.path.dirname(os.path.abspath(__file__)); E = os.path.join(here, "..", "out", "examples"); os.makedirs(E, exist_ok=True)
d = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
q = d[(d.dataset == "ALS") & d.core & d.protein_qualifying & ~d.s_at_canonical_tss & (d.l_up_junc_reads >= 5) & (d["class"] == "strict_ISM")].copy()
keep_f = make_read_filter(HiFi=True, rdna_mask_bed=RDNA)
bam = pysam.AlignmentFile(DATASETS["ALS"]["bam"])
fsm = {}
def n_fsm(chrom, strand, chain):
    k = (chrom, strand, chain)
    if k in fsm: return fsm[k]
    c = parse_chain(chain); want = tuple(sorted(c)); lo, hi = want[0][0], want[-1][1]
    n = 0
    for r in bam.fetch(chrom, lo - 1, lo):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not keep_f(r): continue
        if r.reference_end < hi: continue
        pos = r.reference_start; introns = []
        for o, ln in r.cigartuples:
            if o == 3: introns.append((pos + 1, pos + ln))
            if o in (0, 2, 3, 7, 8): pos += ln
        if tuple(introns) == want: n += 1      # chain fixes the strand (LRAA models are spliced, strand from the GTF)
    fsm[k] = n
    return n
fa = pysam.FastaFile(DATASETS["ALS"]["genome"]); seqs = {}
def n_l_tss(chrom, strand, chain, l_tss, s_end3, win=50, intron_min=30, tol=5, end_tol=50):
    """reads starting within win nt of L's TSS that fit L: introns = L's first ones (5' side), in order, within tol nt;
    3' end inside the matching L exon (end_tol past L's last exon). Gaps < intron_min nt are alignment gaps."""
    if chrom not in seqs:
        seqs[chrom] = fa.fetch(chrom); Util_funcs.register_contig_seq_for_strand_check(chrom, seqs[chrom].upper())
    lx = sorted(s_exons(parse_chain(chain), l_tss, s_end3, strand)); Li = [(lx[i][1], lx[i + 1][0]) for i in range(len(lx) - 1)]
    n = 0
    for r in bam.fetch(chrom, l_tss - win - 10, l_tss + win + 10):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not keep_f(r): continue
        if Util_funcs.transcribed_strand(r) != strand: continue
        x = T._five_prime_start(r, strand, seqs[chrom])
        if x is None or abs(x[0] - l_tss) > win: continue
        ex = []
        for a, b in r.get_blocks():
            if ex and a + 1 - ex[-1][1] < intron_min: ex[-1][1] = max(ex[-1][1], b)
            else: ex.append([a + 1, b])
        k = len(ex) - 1
        if k == 0 or k > len(Li): continue
        ri = [(ex[i][1], ex[i + 1][0]) for i in range(k)]
        want = Li[:k] if strand == "+" else Li[-k:]
        if not all(abs(x0 - y0) <= tol and abs(x1 - y1) <= tol for (x0, x1), (y0, y1) in zip(ri, want)): continue
        if strand == "+": ok = ex[-1][1] <= lx[k][1] + (end_tol if k == len(lx) - 1 else 0)
        else: j = len(lx) - 1 - k; ok = ex[0][0] >= lx[j][0] - (end_tol if j == 0 else 0)
        n += ok
    return n
q["L_fsm_reads"] = [n_fsm(r.chrom, r.strand, r.l_chain) for r in q.itertuples()]
q["S_fsm_reads"] = [n_fsm(r.chrom, r.strand, r.s_chain) for r in q.itertuples()]
q["min_fsm_reads"] = q[["L_fsm_reads", "S_fsm_reads"]].min(axis=1)
q["L_tss_reads"] = [n_l_tss(r.chrom, r.strand, r.l_chain, int(r.l_tss), int(r.s_end3)) if r.min_fsm_reads >= 1 else 0 for r in q.itertuples()]
q["passes"] = (q.min_fsm_reads >= 1) & q.cds_route.eq("GENCODE") & (q.L_tss_reads >= 10)
q["_g"] = q.cds_route.eq("GENCODE").astype(int); q["_r"] = q.regime.eq("refGuided").astype(int)
one = q.sort_values(["passes", "min_fsm_reads", "_g", "_r"], ascending=False).drop_duplicates(["gene", "chrom", "strand", "peak"])
one = one.sort_values(["passes", "n_site", "min_fsm_reads"], ascending=False).reset_index(drop=True)
one["rank"] = range(1, len(one) + 1)
cols = ["rank", "passes", "gene", "regime", "arms", "chrom", "strand", "peak", "l_tss", "s_end3", "s_hash", "l_hash", "n_site", "cap_frac", "l_up_junc_reads",
        "L_fsm_reads", "S_fsm_reads", "min_fsm_reads", "L_tss_reads", "protein_class", "aa_total", "aa_remaining", "cds_route", "L_gencode_tx", "cage_dist", "reftss_dist",
        "fragile_share", "top_cluster", "s_chain", "l_chain"]
with open(f"{E}/ftd_example_candidates.tsv", "w") as o:
    o.write("# " + __doc__.split("usage:")[0].replace("\n", "\n# ").rstrip("# \n") + "\n")
    o.write("# columns: n_site reads at the S site; cap_frac capped share of informative site reads; l_up_junc_reads reads over L's dropped intron;\n"
            "# L_fsm_reads/S_fsm_reads reads reproducing the whole L/S intron chain; fragile_share share of site reads from neuron clusters (ALS 5' truncation artefact);\n"
            "# L_tss_reads reads starting within 50 nt of L's TSS that fit L (splicing and 3' end) (rule 4; counted for rule-3 passers only);\n"
            "# passes = rules 3 and 4 (both isoforms have >= 1 FSM read; GENCODE L CDS; L_tss_reads >= 10). Aggregate counts only, no per-read rows.\n")
    one[cols].to_csv(o, sep="\t", index=False)
print(one[cols[:-2]].to_string())
