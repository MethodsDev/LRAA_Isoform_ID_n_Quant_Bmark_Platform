#!/usr/bin/env python3
"""Upstream (U) and downstream (D) sites of the current LRAA version's high-confidence TSS site switches,
from __Site_Usage_Analysis/<prefix>.TSS_events_vs_CAGE.events.tsv (the table the paper's TSS switch numbers
come from), one row per distinct site pair, in the format event_site_caprates.py reads.
usage: current_event_pairs.py EVENTS_TSV OUT_TSV"""
import sys, csv
ev, out = sys.argv[1:3]
seen = {}
for e in csv.DictReader(open(ev), delimiter="\t"):
    if e["high_confidence"] != "True": continue
    g, l = (e["gained_site"], e["gained_pos"]), (e["lost_site"], e["lost_pos"])
    (u, d) = (l, g) if e["gained_is_downstream"] == "True" else (g, l)
    seen.setdefault((u[0], d[0]), (e["gene_key"], u, d))
with open(out, "wt") as o:
    o.write("gene_key\tU_site\tD_site\tU_pos\tD_pos\n")
    for k, (gk, u, d) in seen.items():
        o.write(f"{gk}\t{u[0]}\t{d[0]}\t{u[1]}\t{d[1]}\n")
print(len(seen), "site pairs")
