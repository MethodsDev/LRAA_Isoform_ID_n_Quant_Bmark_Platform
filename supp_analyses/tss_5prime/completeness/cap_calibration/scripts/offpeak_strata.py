#!/usr/bin/env python3
"""Capped share of informative spliced read starts NOT at a peak, by distance to the nearest peak
and by how many reads start at the position: what a "non-peak background" actually contains.

usage: offpeak_strata.py LRAA_HOME BAM GENOME CONTIG BASELINE OUT_TSV
"""
import sys, os, bisect
from collections import defaultdict
lraa, bam_f, fa_f, contig, baseline, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
sys.path.insert(0, os.path.join(lraa, "pylib"))
import pysam, LRAA_Globals, Util_funcs, TSS_cap_calibration as T
from Splice_graph import aggregate_sites_within_window, filter_non_peaky_positions

cfg = LRAA_Globals.config
seq = pysam.FastaFile(fa_f).fetch(contig)
Util_funcs.register_contig_seq_for_strand_check(contig, seq)
starts = {"+": defaultdict(int), "-": defaultdict(int)}
cap = {"+": defaultdict(lambda: defaultdict(int)), "-": defaultdict(lambda: defaultdict(int))}
for r in pysam.AlignmentFile(bam_f).fetch(contig):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < cfg["min_mapping_quality"]: continue
    st = Util_funcs.transcribed_strand(r); s = T._five_prime_start(r, st, seq)
    if s is None: continue
    pos, g, spl = s; starts[st][pos] += 1
    if g is not None and spl: cap[st][pos][g] += 1
W = cfg["max_dist_between_alt_TSS_sites"]
bins = [(W + 1, 200), (201, 1000), (1001, 5000), (5001, 10**12)]
agg = defaultdict(lambda: [0, 0])
for st in "+-":
    c = starts[st]
    peaks = sorted(p for p, _ in filter_non_peaky_positions(aggregate_sites_within_window(c, W, cfg["min_alignments_define_TSS_site"]), c, contig, st))
    for pos, runs in cap[st].items():
        k = bisect.bisect_left(peaks, pos)
        d = min([abs(pos - peaks[j]) for j in (k - 1, k) if 0 <= j < len(peaks)] or [10**12])
        n = sum(runs.values()); m = sum(v for g, v in runs.items() if g > baseline)
        if d <= W: key = ("at peak (<=%d nt)" % W, "any")
        else:
            lo, hi = next(b for b in bins if b[0] <= d <= b[1])
            key = ("%d-%s nt" % (lo, hi if hi < 10**12 else "inf"), "1" if c[pos] == 1 else ("2-4" if c[pos] < 5 else ">=5"))
        agg[key][0] += n; agg[key][1] += m
with open(out, "wt") as ofh:
    ofh.write("contig\tdistance_to_peak\treads_starting_at_position\tn_informative\tn_capped\n")
    for k in sorted(agg):
        ofh.write("\t".join(map(str, [contig, k[0], k[1], agg[k][0], agg[k][1]])) + "\n")
