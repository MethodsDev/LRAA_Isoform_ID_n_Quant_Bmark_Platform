#!/usr/bin/env python3
"""ISM S/L pairs of one dataset, over all its model sets.

S's intron chain is a proper suffix (transcript sense) of L's: only L's 5'-most introns are
dropped. S keeps >= 1 intron, and S must end within END3_TOL nt of an L model of the same model set (shared 3' end). Models of one chain within one regime are grouped into TSS groups
(5' ends <= 25 nt apart, LRAA's site window). 5' class relative to the L exon that ends at S's
first donor: strict_ISM (inside that exon), TSS_in_L_intron (in the intron upstream of it:
a novel first exon), TSS_past_L_exon (further upstream, spanning an L exon: retained intron).

usage: enumerate_pairs.py DATASET OUT.tsv.gz
"""
import sys, gzip, collections
from common import *
ds, out = sys.argv[1:3]
cfg = DATASETS[ds]
pairs = {}
stats = []
for regime in ("denovo", "refGuided"):
    arms = [a for (r, a) in cfg["gtfs"] if r == regime]
    # (chrom, strand, chain) -> arm -> list of (tss, 3'end, flag); 
    models = collections.defaultdict(lambda: collections.defaultdict(list))
    symbol = {}
    for arm in arms:
        tx = read_exons(cfg["gtfs"][(regime, arm)])
        n = 0
        for tid, r in tx.items():
            if r["chrom"] not in PRIMARY or r["strand"] not in "+-" or len(r["exons"]) < 2: continue
            c = intron_chain(r["exons"])
            tss, end3 = (r["exons"][0][0], r["exons"][-1][1]) if r["strand"] == "+" else (r["exons"][-1][1], r["exons"][0][0])
            k = (r["chrom"], r["strand"], c)
            models[k][arm].append((tss, end3, r["tss_flag"]))
            if r["symbol"]: symbol[k] = r["symbol"]
            n += 1
        stats.append((regime, arm, len(tx), n))
    print(regime, "spliced models", stats[-len(arms):], "chains", len(models), flush=True)
    # TSS groups of each chain over the regime's arms
    keys = collections.defaultdict(list)
    for k in models:
        keys[(k[0], k[1])].append(k)
    for (chrom, strand), ks in keys.items():
        for L in ks:
            cl = L[2]
            n = len(cl)
            for i in range(1, n):        # drop i 5'-most introns
                sc = cl[i:] if strand == "+" else cl[:n - i]
                S = (chrom, strand, sc)
                if S not in models: continue
                shared_arms = set(models[S]) & set(models[L])
                if not shared_arms: continue
                up = cl[i - 1] if strand == "+" else cl[n - i]          # L's intron just 5' of S's first
                first = sc[0] if strand == "+" else sc[-1]               # S's first intron
                pts = sorted(set(m[0] for a in shared_arms for m in models[S][a]))
                groups = [[pts[0]]]
                for p in pts[1:]:
                    (groups[-1].append(p) if p - groups[-1][-1] <= 25 else groups.append([p]))
                for g in groups:
                    ms = [(a, m) for a in shared_arms for m in models[S][a] if g[0] <= m[0] <= g[-1]
                          and any(abs(m[1] - l[1]) <= END3_TOL for l in models[L][a])]
                    if not ms: continue
                    tss = sorted(m[0] for a, m in ms)[len(ms) // 2]
                    if strand == "+":
                        ex_lo, ex_hi = up[1] + 1, first[0] - 1
                        into = tss - ex_lo                    # >=0: inside the L exon (5' to 3' offset)
                        to_donor = ex_hi - tss
                        prev = cl[i - 2] if i >= 2 else None
                        cls = "strict_ISM" if tss >= ex_lo else ("TSS_in_L_intron" if tss >= up[0] else "TSS_past_L_exon")
                    else:
                        ex_lo, ex_hi = first[1] + 1, up[0] - 1
                        into = ex_hi - tss
                        to_donor = tss - ex_lo
                        cls = "strict_ISM" if tss <= ex_hi else ("TSS_in_L_intron" if tss <= up[1] else "TSS_past_L_exon")
                    l_tss = [m[0] for a in shared_arms for m in models[L][a]]
                    pk = (regime, chrom, strand, sc, g[0], cl)
                    r = pairs.setdefault(pk, dict(arms=set(), flag=False, tss=tss, grp=(g[0], g[-1]), ends=[], up=up, first=first,
                                                  exon=(ex_lo, ex_hi), into=into, to_donor=to_donor, cls=cls, dropped=i, n_s=len(sc), n_l=n,
                                                  l_tss=(min(l_tss) if strand == "+" else max(l_tss)), sym=symbol.get(S, symbol.get(L, ""))))
                    r["arms"] |= {a for a, m in ms}
                    r["flag"] |= any(m[2] for a, m in ms)
                    r["ends"] += [m[1] for a, m in ms]
print("pairs", len(pairs), flush=True)
hashes = {}
def h(chrom, strand, c):
    k = (chrom, strand, c)
    if k not in hashes: hashes[k] = splice_hash(chrom, strand, c)
    return hashes[k]
cols = "dataset regime arms chrom strand s_hash l_hash s_n_introns l_n_introns introns_dropped s_tss s_tss_group_lo s_tss_group_hi s_end3 lraa_TSS_flag class l_exon_lo l_exon_hi tss_into_exon tss_to_donor l_tss s_first_intron l_up_intron s_chain l_chain lraa_symbol".split()
with gzip.open(out, "wt") as o:
    o.write("\t".join(cols) + "\n")
    for (regime, chrom, strand, sc, g0, cl), r in pairs.items():
        ends = sorted(r["ends"]); e3 = ends[len(ends) // 2]
        o.write("\t".join(map(str, [ds, regime, ";".join(sorted(r["arms"])), chrom, strand, h(chrom, strand, sc), h(chrom, strand, cl), r["n_s"], r["n_l"], r["dropped"],
            r["tss"], r["grp"][0], r["grp"][1], e3, r["flag"], r["cls"], r["exon"][0], r["exon"][1], r["into"], r["to_donor"], r["l_tss"],
            "%d-%d" % r["first"], "%d-%d" % r["up"], chain_str(sc), chain_str(cl), r["sym"]])) + "\n")
