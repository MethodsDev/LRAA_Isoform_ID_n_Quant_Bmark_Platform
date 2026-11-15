#!/usr/bin/env python3
"""Capped share of informative read starts within 50 nt of each upstream (U) and downstream (D)
site of a dataset's high-confidence TSS switch pairs (ext/<DS>_v0.45.0_TSS_event_pairs.tsv, from v0.45.0,
before LRAA's cap rules; FTD's are mostly neuron-favours-downstream).

LRAA's read rules (TSS_cap_calibration._five_prime_start). Baseline 0 G (sc Kinnex).
usage: event_site_caprates.py LRAA_HOME PAIRS_TSV BAM GENOME OUT_TSV
"""
import sys, os, csv
lraa, pairs_tsv, bam_f, fa_f, out = sys.argv[1:6]
sys.path.insert(0, os.path.join(lraa, "pylib"))
import pysam, Util_funcs, TSS_cap_calibration as T

bam = pysam.AlignmentFile(bam_f); fa = pysam.FastaFile(fa_f); held = [None, None]


def rates(site):
    _, c, p, st = site.split(":"); p = int(p)
    if held[0] != c:
        held[:] = [c, fa.fetch(c)]; Util_funcs.register_contig_seq_for_strand_check(c, held[1])
    seq = held[1]; r = {"all": [0, 0], "spliced": [0, 0]}
    for a in bam.fetch(c, max(0, p - 60), p + 60):
        if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
        s_ = Util_funcs.transcribed_strand(a)
        if s_ != st: continue
        x = T._five_prime_start(a, s_, seq)
        if x is None or abs(x[0] - p) > 50 or x[1] is None: continue
        for k in (("all", "spliced") if x[2] else ("all",)):
            r[k][0] += 1; r[k][1] += x[1] > 0
    return r


with open(out, "wt") as ofh:
    ofh.write("gene_key\trole\tsite\tn_informative\tn_capped\tn_informative_spliced\tn_capped_spliced\n")
    for e in csv.DictReader(open(pairs_tsv), delimiter="\t"):
        for role, site in (("U", e["U_site"]), ("D", e["D_site"])):
            r = rates(site)
            ofh.write("\t".join(map(str, [e["gene_key"], role, site, r["all"][0], r["all"][1], r["spliced"][0], r["spliced"][1]])) + "\n")
