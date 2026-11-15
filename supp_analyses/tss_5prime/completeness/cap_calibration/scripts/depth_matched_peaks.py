#!/usr/bin/env python3
"""Is neurons' excess of uncapped TSS peaks a depth effect? Peaks called within each cell type
after thinning every cell type to the same number of reads.

Pass 1 counts each cell type's primary reads with a cell barcode in the fixed clusters. Each
cell type is then thinned to the smallest cell type's count (a read is kept when a hash of its
name falls under the cell type's keep fraction, so the choice is fixed and reproducible). Pass 2
collects, per cell type, the trusted 5' starts of the kept reads (LRAA's rules), calls peaks on
them as an LRAA graph build does (aggregate_sites_within_window, filter_non_peaky_positions),
and counts each peak's informative spliced reads by capped (G run > baseline) or not. A peak is
judged with >= MIN_READS such reads and classed uncapped below the library's Otsu split t.

Writes per cell type: barcoded reads and cells before thinning (reads per cell), reads kept,
peaks called, peaks judged, share of judged peaks uncapped, share of their reads at uncapped
peaks. Also the same without thinning, for the depth effect itself.

usage: depth_matched_peaks.py LRAA_HOME BAM GENOME CALIBRATION_JSON CELL_CLUSTERS DATASET OUT_TSV [--cpu N] [--min_reads N]
"""
import sys, os, json, zlib, argparse, multiprocessing
from collections import defaultdict

CELL_TYPE = {
    "FTD": {"0": "Oligodendrocyte", "3": "Oligodendrocyte", "1": "Excitatory neuron", "8": "Excitatory neuron",
            "4": "Inhibitory neuron", "6": "Inhibitory neuron", "2": "Astrocyte", "7": "Astrocyte", "12": "Astrocyte",
            "5": "Microglia", "9": "OPC", "10": "Vascular / meningeal", "11": "Vascular / meningeal",
            "13": "Vascular / meningeal"},
    "PBMC": {"1": "Monocyte / DC", "11": "Monocyte / DC", "13": "Monocyte / DC", "6": "B", "9": "B",
             "3": "CD4 T", "7": "CD4 T", "0": "CD4 T", "5": "CD4 T", "12": "CD4 T",
             "2": "CD8 T / NK", "4": "CD8 T / NK", "8": "CD8 T / NK", "10": "CD8 T / NK"},
}


def parse_args():
    p = argparse.ArgumentParser()
    for a in ("lraa_home", "bam", "genome", "calibration", "cell_clusters", "dataset", "out"):
        p.add_argument(a)
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--min_reads", type=int, default=20)
    return p.parse_args()


def keep(name, frac):
    return frac >= 1.0 or (zlib.crc32(name.encode()) % 1000000) / 1e6 < frac


def count_contig(task):
    bam_f, contig, cells = task
    import pysam
    n = defaultdict(int)
    with pysam.AlignmentFile(bam_f) as bam:
        for r in bam.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.has_tag("CB"):
                continue
            ct = cells.get(r.get_tag("CB"))
            if ct is not None:
                n[ct] += 1
    return dict(n)


def scan_contig(task):
    lraa_home, bam_f, fa_f, contig, cells, fracs, baseline, cfg = task
    sys.path.insert(0, os.path.join(lraa_home, "pylib"))
    import bisect, pysam, LRAA_Globals, Util_funcs, TSS_cap_calibration as T
    from Splice_graph import aggregate_sites_within_window, filter_non_peaky_positions
    LRAA_Globals.config.update(cfg)
    seq = pysam.FastaFile(fa_f).fetch(contig)
    Util_funcs.register_contig_seq_for_strand_check(contig, seq)
    W = cfg["max_dist_between_alt_TSS_sites"]
    res = {}
    for mode in ("thinned", "all"):
        starts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))   # ct -> strand -> pos -> n
        cap = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0])))  # ct -> strand -> pos -> [n, capped]
        with pysam.AlignmentFile(bam_f) as bam:
            for r in bam.fetch(contig):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.has_tag("CB"):
                    continue
                ct = cells.get(r.get_tag("CB"))
                if ct is None or (mode == "thinned" and not keep(r.query_name, fracs[ct])):
                    continue
                st = Util_funcs.transcribed_strand(r)
                x = T._five_prime_start(r, st, seq)
                if x is None:
                    continue
                p, G, spliced = x
                starts[ct][st][p] += 1
                if G is not None and spliced:
                    c = cap[ct][st][p]
                    c[0] += 1
                    c[1] += G > baseline
        out = {}
        for ct in starts:
            peaks_called, judged = 0, []
            for st in starts[ct]:
                counter = starts[ct][st]
                grouped = aggregate_sites_within_window(counter, W, cfg["min_alignments_define_TSS_site"])
                peaks = sorted(p for p, _ in filter_non_peaky_positions(grouped, counter, contig, st))
                peaks_called += len(peaks)
                positions = sorted(cap[ct][st])
                for site in peaks:
                    lo = bisect.bisect_left(positions, site - W)
                    hi = bisect.bisect_right(positions, site + W)
                    n = m = 0
                    for q in positions[lo:hi]:
                        n += cap[ct][st][q][0]
                        m += cap[ct][st][q][1]
                    judged.append((n, m))
            out[ct] = (peaks_called, judged)
        res[mode] = out
    return res


def main():
    a = parse_args()
    sys.path.insert(0, os.path.join(a.lraa_home, "pylib"))
    import pysam, LRAA_Globals
    cal = json.load(open(a.calibration))["stats"]
    baseline, t = cal["baseline_G_run"], cal["mode_split"]
    cmap = CELL_TYPE[a.dataset]
    cells, n_cells = {}, defaultdict(int)
    for line in open(a.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1] in cmap:
            cells[f[0]] = cmap[f[1]]
            n_cells[cmap[f[1]]] += 1
    with pysam.AlignmentFile(a.bam) as bam:
        contigs = [s.contig for s in bam.get_index_statistics() if s.mapped > 0 and "_" not in s.contig and s.contig != "chrM"]
    with pysam.FastaFile(a.genome) as fa:
        contigs = [c for c in contigs if c in set(fa.references)]
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(a.cpu) as pool:
        n_reads = defaultdict(int)
        for d in pool.imap_unordered(count_contig, [(a.bam, c, cells) for c in contigs]):
            for ct, n in d.items():
                n_reads[ct] += n
        target = min(n_reads.values())
        fracs = {ct: target / n for ct, n in n_reads.items()}
        cfg = dict(LRAA_Globals.config)
        totals = {m: defaultdict(lambda: [0, []]) for m in ("thinned", "all")}
        for res in pool.imap_unordered(scan_contig, [(a.lraa_home, a.bam, a.genome, c, cells, fracs, baseline, cfg) for c in contigs]):
            for mode, out in res.items():
                for ct, (pc, judged) in out.items():
                    totals[mode][ct][0] += pc
                    totals[mode][ct][1].extend(judged)
    with open(a.out, "wt") as ofh:
        ofh.write("dataset\tcell_type\tcells\tbarcoded_reads\treads_per_cell\tdepth\treads_kept\tpeaks_called\t"
                  "peaks_judged\tuncapped_peak_share\treads_at_judged_peaks\tuncapped_read_share\tcapped_frac_at_capped_peaks\n")
        for mode in ("all", "thinned"):
            for ct in sorted(n_reads, key=lambda x: -n_reads[x]):
                pc, judged = totals[mode][ct]
                j = [(n, m) for n, m in judged if n >= a.min_reads]
                unc = [(n, m) for n, m in j if m / n < t]
                capd = [(n, m) for n, m in j if m / n >= t]
                nj = sum(n for n, _ in j)
                ofh.write("\t".join(map(str, [
                    a.dataset, ct, n_cells[ct], n_reads[ct], round(n_reads[ct] / max(n_cells[ct], 1), 1),
                    "matched" if mode == "thinned" else "full",
                    target if mode == "thinned" else n_reads[ct], pc, len(j),
                    round(len(unc) / len(j), 4) if j else "NA", nj,
                    round(sum(n for n, _ in unc) / nj, 4) if nj else "NA",
                    round(sum(m for _, m in capd) / max(sum(n for n, _ in capd), 1), 4)])) + "\n")
    print("target reads per cell type", target)


if __name__ == "__main__":
    main()
