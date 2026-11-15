#!/usr/bin/env python3
"""Per cell cluster: how a cluster's 5' read starts split between the library's uncapped-mode
and capped-mode TSS peaks, and how often its reads carry the cap at each.

Peaks and their mode come from the library calibration (out/<lib>.peaks.tsv.gz; judged peaks,
split at the calibration's Otsu threshold t). Reads use LRAA's rules
(TSS_cap_calibration._five_prime_start): trusted 5' end, informative (no genomic G upstream),
spliced. A read counts for the nearest judged peak within 50 nt; its cluster comes from its CB
tag and the fixed cell clusters.

usage: peak_mode_by_cluster.py LRAA_HOME BAM GENOME PEAKS_TSV_GZ CALIBRATION_JSON CELL_CLUSTERS OUT_TSV [--cpu N]
"""
import sys, os, gzip, json, bisect, argparse, multiprocessing
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    for a in ("lraa_home", "bam", "genome", "peaks", "calibration", "cell_clusters", "out"):
        p.add_argument(a)
    p.add_argument("--cpu", type=int, default=8)
    return p.parse_args()


def scan_contig(task):
    lraa_home, bam_f, fa_f, contig, peaks, cells, baseline, window = task
    sys.path.insert(0, os.path.join(lraa_home, "pylib"))
    import pysam, Util_funcs, TSS_cap_calibration as T

    seq = pysam.FastaFile(fa_f).fetch(contig)
    Util_funcs.register_contig_seq_for_strand_check(contig, seq)
    pos = {s: sorted(peaks.get(s, {})) for s in "+-"}
    # cluster -> mode -> [reads, capped]
    out = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    with pysam.AlignmentFile(bam_f) as bam:
        for r in bam.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.has_tag("CB"):
                continue
            cl = cells.get(r.get_tag("CB"))
            if cl is None:
                continue
            st = Util_funcs.transcribed_strand(r)
            x = T._five_prime_start(r, st, seq)
            if x is None or x[1] is None or not x[2]:
                continue
            p = x[0]
            ps = pos[st]
            k = bisect.bisect_left(ps, p)
            near = [ps[j] for j in (k - 1, k) if 0 <= j < len(ps) and abs(ps[j] - p) <= window]
            if not near:
                continue
            site = min(near, key=lambda q: abs(q - p))
            mode = peaks[st][site]
            c = out[cl][mode]
            c[0] += 1
            c[1] += x[1] > baseline
    return {cl: dict(m) for cl, m in out.items()}


def main():
    a = parse_args()
    cal = json.load(open(a.calibration))["stats"]
    baseline, t = cal["baseline_G_run"], cal["mode_split"]
    min_reads = cal["peak_min_informative_reads"]
    peaks = defaultdict(lambda: defaultdict(dict))
    with gzip.open(a.peaks, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(h, line.rstrip("\n").split("\t")))
            if int(f["n_informative"]) < min_reads:
                continue
            mode = "capped" if float(f["capped_frac"]) >= t else "uncapped"
            peaks[f["contig"]][f["strand"]][int(f["pos"])] = mode
    cells = {}
    for line in open(a.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] not in ("", "x", "cell_barcode"):
            cells[f[0]] = f[1]
    tasks = [(a.lraa_home, a.bam, a.genome, c, {s: dict(v) for s, v in peaks[c].items()}, cells, baseline, 50)
             for c in sorted(peaks, key=lambda c: -sum(len(v) for v in peaks[c].values()))]
    total = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    with multiprocessing.get_context("spawn").Pool(min(a.cpu, len(tasks))) as pool:
        for res in pool.imap_unordered(scan_contig, tasks):
            for cl, m in res.items():
                for mode, (n, k) in m.items():
                    total[cl][mode][0] += n
                    total[cl][mode][1] += k
    with open(a.out, "wt") as ofh:
        ofh.write("cluster\tmode\treads\tcapped_reads\n")
        for cl in sorted(total, key=lambda x: int(x) if x.isdigit() else x):
            for mode in ("uncapped", "capped"):
                n, k = total[cl].get(mode, [0, 0])
                ofh.write("{}\t{}\t{}\t{}\n".format(cl, mode, n, k))
    print("clusters", len(total), "peaks", sum(len(v) for c in peaks.values() for v in c.values()), "t", t)


if __name__ == "__main__":
    main()
