#!/usr/bin/env python3
"""
5' reach of polyA-anchored reads at control genes, with the same sites in both datasets.

For each gene, from the reads of all datasets pooled (reads on the gene's strand, passing LRAA's read filters as
read_end_profile_by_cluster.py applies them, from barcodes with a cluster):
  - dominant PolyA site: the most common read 3' end (10-nt bins, the bin's most common position);
  - dominant TSS: the most common read 5' end, likewise;
  - dominant isoform: the most common intron chain among reads starting within WIN nt of the TSS and ending
    within WIN nt of the PolyA site (taken from the reads, so annotation artefacts such as OAZ1's 1-nt
    frameshift "intron" don't matter).
A read is counted when it ends within WIN nt of the PolyA site and follows the dominant isoform from its 5' end
on: its 5' end lies in an exon of the isoform, and its introns are exactly the isoform's introns downstream of
that point (+/- 2 bp). reaches_tss: the read starts within WIN nt of the TSS.

Outputs:
  <prefix>.sites.tsv          per gene: the two sites and the isoform's intron chain
  <prefix>.reads.tsv.gz       per counted read: dataset, gene, cluster, reaches_tss
  <prefix>.anchoring.tsv      per gene and dataset: reads following the isoform's splicing, and the share of
                              them ending at the PolyA site
"""

import argparse
import collections
import gzip
import os
import sys

import pysam

LRAA_HOME = os.environ.get("LRAA_HOME", os.path.expanduser("~/GITHUB/MDL/LongReadAlignmentAssembler"))
sys.path.insert(0, os.path.join(LRAA_HOME, "util/sc/diff_iso_usage"))
from alt_termini_read_check import (make_read_filter, parse_cell_clusters, read_geometry,  # noqa: E402
                                    transcript_strand, register_contig)

WIN = 50
INTRON_TOL = 2


def gene_spans(gtf, genes):
    want, spans = set(genes), {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or 'gene_name "' not in f[8]:
            continue
        g = f[8].split('gene_name "', 1)[1].split('"', 1)[0]
        if g in want:
            spans[g] = (f[0], int(f[3]), int(f[4]), f[6])
    missing = want - set(spans)
    if missing:
        sys.exit(f"genes not in {gtf}: {sorted(missing)}")
    return spans


def mode_pos(positions):
    """most common position within the most populated 10-nt bin"""
    bins = collections.Counter(p // 10 for p in positions)
    b = bins.most_common(1)[0][0]
    return collections.Counter(p for p in positions if p // 10 == b).most_common(1)[0][0]


def same_intron(a, b):
    return abs(a[0] - b[0]) <= INTRON_TOL and abs(a[1] - b[1]) <= INTRON_TOL


def follows(read, chain, strand):
    """read's 5' end in an exon of the isoform and its introns = the isoform's introns downstream of it"""
    five = read["five"]
    if any(i[0] <= five <= i[1] for i in chain):
        return False
    down = [i for i in chain if i[0] > five] if strand == "+" else [i for i in chain if i[1] < five]
    introns = sorted(read["introns"])
    return len(introns) == len(down) and all(same_intron(a, b) for a, b in zip(introns, sorted(down)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", action="append", nargs=3, metavar=("NAME", "BAM", "CELL_CLUSTERS"), required=True)
    ap.add_argument("--gtf", required=True, help="gene spans (gene_name attributes)")
    ap.add_argument("--genes", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--HiFi", action="store_true")
    ap.add_argument("--rdna_mask_bed", default=None)
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    keep = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
    genome = pysam.FastaFile(args.genome)
    genes = args.genes.split(",")
    spans = gene_spans(args.gtf, genes)

    reads = collections.defaultdict(list)       # gene -> read dicts, all datasets
    for name, bam_path, clusters in args.dataset:
        cell_to_cluster = parse_cell_clusters(clusters)
        bam = pysam.AlignmentFile(bam_path)
        for g in genes:
            contig, lo, hi, strand = spans[g]
            register_contig(genome, contig)
            for r in bam.fetch(contig, lo - 1, hi):
                if not r.has_tag("CB") or r.get_tag("CB") not in cell_to_cluster:
                    continue
                if not keep(r) or transcript_strand(r, genome) != strand:
                    continue
                lend, rend, introns = read_geometry(r)
                if rend < lo or lend > hi:
                    continue
                five, three = (lend, rend) if strand == "+" else (rend, lend)
                reads[g].append(dict(dataset=name, cluster=cell_to_cluster[r.get_tag("CB")], lend=lend, rend=rend,
                                     five=five, three=three, introns=[tuple(i) for i in introns]))

    with open(f"{args.output_prefix}.sites.tsv", "w") as s_out, \
            open(f"{args.output_prefix}.anchoring.tsv", "w") as a_out, \
            gzip.open(f"{args.output_prefix}.reads.tsv.gz", "wt") as r_out:
        s_out.write("gene\tcontig\tstrand\ttss\tpolyA\tisoform_introns\tfull_length_reads\tisoform_reads\n")
        a_out.write("gene\tdataset\treads_following_isoform\tshare_ending_at_polyA\n")
        r_out.write("dataset\tgene\tcluster\treaches_tss\n")
        for g in genes:
            contig, _, _, strand = spans[g]
            rs = reads[g]
            polyA = mode_pos([x["three"] for x in rs])
            tss = mode_pos([x["five"] for x in rs])
            full = [x for x in rs if abs(x["five"] - tss) <= WIN and abs(x["three"] - polyA) <= WIN]
            chains = collections.Counter(tuple(sorted(x["introns"])) for x in full)
            chain, n_chain = chains.most_common(1)[0]
            s_out.write(f"{g}\t{contig}\t{strand}\t{tss}\t{polyA}\t{';'.join(f'{a}-{b}' for a, b in chain)}\t"
                        f"{len(full)}\t{n_chain}\n")
            follow = [x for x in rs if follows(x, list(chain), strand)]
            for ds in sorted({x["dataset"] for x in rs}):
                f_ds = [x for x in follow if x["dataset"] == ds]
                at = sum(abs(x["three"] - polyA) <= WIN for x in f_ds)
                a_out.write(f"{g}\t{ds}\t{len(f_ds)}\t{at / len(f_ds) if f_ds else 'NA'}\n")
            for x in follow:
                if abs(x["three"] - polyA) <= WIN:
                    r_out.write(f"{x['dataset']}\t{g}\t{x['cluster']}\t{int(abs(x['five'] - tss) <= WIN)}\n")


if __name__ == "__main__":
    main()
