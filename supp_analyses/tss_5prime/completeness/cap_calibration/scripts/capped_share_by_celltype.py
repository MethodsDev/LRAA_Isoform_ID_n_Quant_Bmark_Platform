#!/usr/bin/env python3
"""Per cell type: share of its reads whose 5' end is capped -- a per-read measure of how many
molecules reach their capped 5' end, independent of peak calling and so of read depth.

Reads use LRAA's rules (TSS_cap_calibration._five_prime_start): trusted 5' end, informative
(no genomic G upstream), spliced; capped = 5' G run above the library's baseline. Also the
aligned length of those reads (median), for comparison with 5' completeness.

usage: capped_share_by_celltype.py LRAA_HOME BAM GENOME CALIBRATION_JSON CELL_CLUSTERS DATASET OUT_TSV [--cpu N]
"""
import sys, os, json, argparse, multiprocessing
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from depth_matched_peaks import CELL_TYPE


def parse_args():
    p = argparse.ArgumentParser()
    for a in ("lraa_home", "bam", "genome", "calibration", "cell_clusters", "dataset", "out"):
        p.add_argument(a)
    p.add_argument("--cpu", type=int, default=8)
    return p.parse_args()


def scan_contig(task):
    lraa_home, bam_f, fa_f, contig, cells, baseline = task
    sys.path.insert(0, os.path.join(lraa_home, "pylib"))
    import pysam, Util_funcs, TSS_cap_calibration as T
    seq = pysam.FastaFile(fa_f).fetch(contig)
    Util_funcs.register_contig_seq_for_strand_check(contig, seq)
    n = defaultdict(lambda: [0, 0])       # cell type -> [informative spliced reads, capped]
    lengths = defaultdict(lambda: defaultdict(int))  # cell type -> aligned length bin (50 nt) -> reads
    with pysam.AlignmentFile(bam_f) as bam:
        for r in bam.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.has_tag("CB"):
                continue
            ct = cells.get(r.get_tag("CB"))
            if ct is None:
                continue
            x = T._five_prime_start(r, Util_funcs.transcribed_strand(r), seq)
            if x is None or x[1] is None or not x[2]:
                continue
            n[ct][0] += 1
            n[ct][1] += x[1] > baseline
            lengths[ct][r.query_alignment_length // 50] += 1
    return {ct: v for ct, v in n.items()}, {ct: dict(v) for ct, v in lengths.items()}


def main():
    a = parse_args()
    sys.path.insert(0, os.path.join(a.lraa_home, "pylib"))
    import pysam
    baseline = json.load(open(a.calibration))["stats"]["baseline_G_run"]
    cmap = CELL_TYPE[a.dataset]
    cells = {}
    for line in open(a.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1] in cmap:
            cells[f[0]] = cmap[f[1]]
    with pysam.AlignmentFile(a.bam) as bam:
        contigs = [s.contig for s in bam.get_index_statistics() if s.mapped > 0 and "_" not in s.contig and s.contig != "chrM"]
    with pysam.FastaFile(a.genome) as fa:
        contigs = [c for c in contigs if c in set(fa.references)]
    tot = defaultdict(lambda: [0, 0])
    lens = defaultdict(lambda: defaultdict(int))
    with multiprocessing.get_context("spawn").Pool(a.cpu) as pool:
        for n, L in pool.imap_unordered(scan_contig, [(a.lraa_home, a.bam, a.genome, c, cells, baseline) for c in contigs]):
            for ct, (i, k) in n.items():
                tot[ct][0] += i
                tot[ct][1] += k
            for ct, d in L.items():
                for b, c in d.items():
                    lens[ct][b] += c

    def median_len(d):
        total, run = sum(d.values()), 0
        for b in sorted(d):
            run += d[b]
            if run >= total / 2:
                return b * 50 + 25
    with open(a.out, "wt") as ofh:
        ofh.write("dataset\tcell_type\tinformative_spliced_reads\tcapped_reads\tcapped_share\tmedian_aligned_length\n")
        for ct in sorted(tot, key=lambda x: -tot[x][0]):
            i, k = tot[ct]
            ofh.write(f"{a.dataset}\t{ct}\t{i}\t{k}\t{k / i:.4f}\t{median_len(lens[ct])}\n")


if __name__ == "__main__":
    main()
