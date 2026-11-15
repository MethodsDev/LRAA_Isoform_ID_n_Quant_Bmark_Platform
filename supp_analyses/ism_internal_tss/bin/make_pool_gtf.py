#!/usr/bin/env python3
"""Union of the benchmark tools' de novo GTFs and GENCODE v47 for one bulk library, plus each tool's chain -> [(5' end, 3' end)] index.
usage: make_pool_gtf.py SAMPLE"""
import sys, gzip, pickle, csv
from common import *
s = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
pool = os.path.join(here, "..", "out", "pool"); os.makedirs(pool, exist_ok=True)
man = [r for r in csv.DictReader(open(os.path.join(REPO, "__Benchmark_extras/isoform_end_accuracy/out/manifest.tsv")), delimiter="\t")
       if r["regime"] == "denovo" and r["sample"] == f"CL_{s}_E0_human"]
idx = {}
with gzip.open(os.path.join(pool, f"{s}.pool.gtf.gz"), "wt") as o:
    def put(src, tid, r):
        t = f"{src}.{tid}"
        for a, b in r["exons"]:
            o.write(f'{r["chrom"]}\tpool\texon\t{a}\t{b}\t.\t{r["strand"]}\t.\ttranscript_id "{t}";\n')
    for m in man:
        tx = read_exons(m["gtf"]); fam = m["family"]; d = idx.setdefault(fam, {})
        for tid, r in tx.items():
            if r["chrom"] not in PRIMARY or r["strand"] not in "+-" or len(r["exons"]) < 2: continue
            tss, e3 = (r["exons"][0][0], r["exons"][-1][1]) if r["strand"] == "+" else (r["exons"][-1][1], r["exons"][0][0])
            d.setdefault((r["chrom"], r["strand"], intron_chain(r["exons"])), []).append((tss, e3))
            put(fam, tid, r)
        print(fam, len(tx), len(d), flush=True)
    for tid, r in read_exons(GENCODE47).items():
        if r["chrom"] in PRIMARY and r["strand"] in "+-" and len(r["exons"]) > 1: put("GENCODE", tid, r)
pickle.dump(idx, open(os.path.join(pool, f"{s}.tool_index.pkl"), "wb"))
