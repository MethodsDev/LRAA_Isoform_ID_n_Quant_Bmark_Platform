#!/usr/bin/env python3
"""Cache the GENCODE v47 transcripts (exons, CDS, start codon, tags). usage: gencode_cache.py OUT.pkl"""
import sys, pickle, collections
from common import *
tx = {}
with open(GENCODE47) as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if p[2] not in ("exon", "CDS", "transcript", "start_codon", "stop_codon"): continue
        tid = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
        r = tx.setdefault(tid, dict(chrom=p[0], strand=p[6], exons=[], cds=[], start=[], stop=[]))
        if p[2] == "transcript":
            a = attrs(p[8])
            r.update(gene=a["gene_name"][0] if a["gene_name"] else a["gene_id"][0], gene_id=a["gene_id"][0], ttype=a["transcript_type"][0],
                     tags=set(a["tag"]), tsl=(a["transcript_support_level"] or ["NA"])[0], pid=(a["protein_id"] or [""])[0])
        elif p[2] == "exon": r["exons"].append((int(p[3]), int(p[4])))
        elif p[2] == "CDS": r["cds"].append((int(p[3]), int(p[4])))
        elif p[2] == "start_codon": r["start"].append((int(p[3]), int(p[4])))
        elif p[2] == "stop_codon": r["stop"].append((int(p[3]), int(p[4])))
for r in tx.values():
    r["exons"].sort(); r["cds"].sort(); r["chain"] = intron_chain(r["exons"])
    r["tss"] = r["exons"][0][0] if r["strand"] == "+" else r["exons"][-1][1]
pickle.dump(tx, open(sys.argv[1], "wb"))
print(len(tx), "transcripts", collections.Counter(r["ttype"] for r in tx.values()).most_common(6))
