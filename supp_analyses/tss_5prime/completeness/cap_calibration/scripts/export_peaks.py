#!/usr/bin/env python3
"""Per-peak cap-signature counts for one library, with LRAA's library calibration.

Runs LRAA's own pre-pass (pylib/TSS_cap_calibration.py: scan, decide, write_report), so the
peaks, read rules and decision are exactly what LRAA uses. Writes
  <prefix>.peaks.tsv.gz                 one row per TSS peak (all peaks; the report filters)
  <prefix>.TSS_cap_calibration.{json,hist.tsv,svg}   LRAA's calibration report

usage: export_peaks.py LRAA_HOME BAM GENOME OUT_PREFIX [--cpu N] [--cage BED.gz] [--contigs a,b]
"""
import sys, os, gzip, bisect, argparse, collections


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("lraa_home"); p.add_argument("bam"); p.add_argument("genome"); p.add_argument("out_prefix")
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--cage", default=None, help="FANTOM5 CAGE peaks bed.gz (human only)")
    p.add_argument("--contigs", default=None)
    return p.parse_args()


def main():
    a = parse_args()
    sys.path.insert(0, os.path.join(a.lraa_home, "pylib"))
    import TSS_cap_calibration as T

    per_contig = T.scan(a.bam, a.genome, cpu=a.cpu, contigs=a.contigs.split(",") if a.contigs else None)
    baseline, active, stats = T.decide(per_contig)
    overrides = {}
    if stats["decision"] in ("on", "off"):
        overrides = {"TSS_cap_signature": stats["decision"], "TSS_cap_G_baseline": baseline}
    T.write_report(a.out_prefix, overrides, stats)

    cage = collections.defaultdict(list)
    if a.cage:
        for l in gzip.open(a.cage, "rt"):
            f = l.split("\t")
            cage[(f[0], f[5].strip())].append((int(f[1]), int(f[2])))
        for k in cage:
            cage[k].sort()

    def in_cage(c, s, pos, slop=25):
        iv = cage[(c, s)]
        i = bisect.bisect_right(iv, (pos, 10**12)) - 1
        return any(iv[j][0] - slop < pos <= iv[j][1] + slop for j in range(max(0, i - 5), min(len(iv), i + 2)))

    b = baseline if baseline is not None else 0
    with gzip.open(a.out_prefix + ".peaks.tsv.gz", "wt") as ofh:
        ofh.write("\t".join(["contig", "strand", "pos", "n_informative", "n_capped", "capped_frac", "in_CAGE"]
                            + ["G%d" % g for g in range(T.MAX_G_RUN + 1)]) + "\n")
        for c in sorted(per_contig):
            for s in ("+", "-"):
                for site, runs in per_contig[c][s]["peaks"]:
                    runs = {int(k): v for k, v in runs.items()}
                    n = sum(runs.values()); m = sum(v for g, v in runs.items() if g > b)
                    ofh.write("\t".join(map(str, [c, s, site, n, m, round(m / n, 4) if n else "NA",
                                                  ("TRUE" if in_cage(c, s, site) else "FALSE") if a.cage else "NA"]
                                            + [runs.get(g, 0) for g in range(T.MAX_G_RUN + 1)])) + "\n")
    print(a.out_prefix, stats["decision"], "baseline", baseline, "peaks judged", stats.get("peaks_judged"))


if __name__ == "__main__":
    main()
