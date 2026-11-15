#!/usr/bin/env python3
"""One pass of a BAM for the study's read evidence. LRAA's read filters (make_read_filter, HiFi) and
5' end rules (TSS_cap_calibration._five_prime_start; G run > baseline = capped).

For reads whose transcript-first intron is a candidate S first intron (+-2 nt, as tandem_tss_features):
  starts: strand first_intron pos cluster fullS  n n_informative n_capped  (fullS: read chain == an S chain)
For every read, counts of the L upstream introns (exact): junc.

usage: scan_reads.py DATASET CPU
"""
import sys, os, gzip, bisect, collections, multiprocessing as mp
from common import *
import pandas as pd
sys.path[:0] = [os.path.join(LRAA_HOME, "pylib"), os.path.join(LRAA_HOME, "util/sc/diff_iso_usage")]
ds, cpu = sys.argv[1], int(sys.argv[2])
cfg = DATASETS[ds]
OUT = os.path.join(os.path.dirname(__file__), "..", "out", "scan")

def load():
    p = pd.read_csv(os.path.join(OUT, "..", "pairs", f"{ds}.pairs.tsv.gz"), sep="\t", usecols=["chrom", "strand", "s_first_intron", "l_up_intron", "s_chain"])
    return p

def job(chrom):
    import pysam, Util_funcs, LRAA_Globals, TSS_cap_calibration as T
    from alt_termini_read_check import make_read_filter
    keep = make_read_filter(HiFi=True, rdna_mask_bed=RDNA)
    cells = {}
    if cfg["clusters"]:
        with open(cfg["clusters"]) as fh:
            next(fh)
            for l in fh:
                b, c = l.split()
                cells[b] = c
    p = PAIRS[PAIRS.chrom == chrom]
    first, schains, junc = {"+": set(), "-": set()}, {"+": set(), "-": set()}, {"+": set(), "-": set()}
    for s, f, u, sc in zip(p.strand, p.s_first_intron, p.l_up_intron, p.s_chain):
        first[s].add(tuple(map(int, f.split("-")))); junc[s].add(tuple(map(int, u.split("-")))); schains[s].add(parse_chain(sc))
    fsort = {s: sorted(first[s]) for s in first}
    def match(s, it):
        if it in first[s]: return it
        L = fsort[s]; i = bisect.bisect_left(L, (it[0] - 2, 0))
        while i < len(L) and L[i][0] <= it[0] + 2:
            if abs(L[i][1] - it[1]) <= 2: return L[i]
            i += 1
        return None
    fa = pysam.FastaFile(cfg["genome"]); seq = fa.fetch(chrom)
    Util_funcs.register_contig_seq_for_strand_check(chrom, seq.upper())
    bam = pysam.AlignmentFile(cfg["bam"])
    st = collections.defaultdict(lambda: [0, 0, 0]); jc = collections.Counter(); seen = 0; used = 0
    base = cfg["baseline"]
    for r in bam.fetch(chrom):
        seen += 1
        ct = r.cigartuples
        if r.is_unmapped or r.is_secondary or r.is_supplementary or not ct or not any(o == 3 for o, _ in ct): continue
        if not keep(r): continue
        strand = Util_funcs.transcribed_strand(r)
        if strand not in ("+", "-"): continue
        pos = r.reference_start; introns = []
        for o, n in ct:
            if o == 3: introns.append((pos + 1, pos + n))
            if o in (0, 2, 3, 7, 8): pos += n
        for it in introns:
            if it in junc[strand]: jc[(strand, it)] += 1
        fi = match(strand, introns[0] if strand == "+" else introns[-1])
        if fi is None: continue
        x = T._five_prime_start(r, strand, seq)
        if x is None: continue
        p5, G, _ = x
        full = tuple(introns) in schains[strand]
        cl = cells.get(r.get_tag("CB"), "NA") if (cells and r.has_tag("CB")) else ("NA" if cells else "all")
        v = st[(strand, fi, p5, cl, full)]
        v[0] += 1
        if G is not None:
            v[1] += 1
            if G > base: v[2] += 1
        used += 1
    with gzip.open(f"{OUT}/{ds}.{chrom}.starts.tsv.gz", "wt") as o:
        for (s, fi, p5, cl, full), v in st.items():
            o.write(f"{chrom}\t{s}\t{fi[0]}-{fi[1]}\t{p5}\t{cl}\t{int(full)}\t{v[0]}\t{v[1]}\t{v[2]}\n")
    with gzip.open(f"{OUT}/{ds}.{chrom}.junc.tsv.gz", "wt") as o:
        for (s, it), n in jc.items():
            o.write(f"{chrom}\t{s}\t{it[0]}-{it[1]}\t{n}\n")
    return chrom, seen, used

PAIRS = load()
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    chroms = [c for c in PRIMARY if (PAIRS.chrom == c).any()]
    with mp.Pool(cpu) as pool:
        for c, seen, used in pool.imap_unordered(job, chroms):
            print(c, seen, used, flush=True)
    os.system(f"cd {OUT} && cat {' '.join(f'{ds}.{c}.starts.tsv.gz' for c in chroms)} > {ds}.starts.tsv.gz && cat {' '.join(f'{ds}.{c}.junc.tsv.gz' for c in chroms)} > {ds}.junc.tsv.gz && rm {ds}.chr*.tsv.gz")
    print("done")
