#!/usr/bin/env python3
"""Domain, signal peptide and transmembrane composition of the N-truncated sites' long (L) and short (S) proteins
(user, 2026-11-11). L's protein is translated from its CDS (GENCODE or pyTD route, as in protein.py); S's is L's from residue
aa_lost + 1, the first downstream in-frame ATG. Each distinct sequence goes once to InterProScan at EBI (Pfam, SignalP
eukaryote, Phobius, TMHMM, Coils); raw results are cached in ext/interpro/<sha1>.tsv. A Pfam domain of L is lost when it
ends at or before residue aa_lost, partial when it spans it, kept otherwise.
usage: protein_domains.py EMAIL     (EBI requires a contact address per job)   writes out/protein_domains.tsv"""
import sys, os, time, hashlib, pickle, urllib.request, urllib.parse
import pandas as pd, pysam
from common import *
from protein_lib import tx_seq
here = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(here, "..", "out"); X = os.path.join(here, "..", "ext", "interpro")
os.makedirs(X, exist_ok=True)
EMAIL = sys.argv[1]
API = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
APPL = ["PfamA", "SignalP_EUK", "Phobius", "TMHMM", "Coils"]
CODON = {a + b + c: aa for (a, b, c), aa in zip(((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
         "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
def translate(s): return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))

d = pd.read_csv(os.path.join(O, "candidates.full.tsv.gz"), sep="\t", low_memory=False)
d["qualifying_strict"] = d.core & d.protein_qualifying & ~d.s_at_canonical_tss & (d.l_up_junc_reads >= 5) & (d["class"] == "strict_ISM")   # as summarize.py
d = d[d.qualifying_strict & d.primary_L].copy()
gc = pickle.load(open(os.path.join(O, "gencode47.pkl"), "rb"))
orf = {}
for r in pd.read_csv(os.path.join(O, "orf", "orfs.tsv"), sep="\t").itertuples():
    orf[(r.chrom, r.strand, r.l_chain)] = [tuple(map(int, x.split("-"))) for x in r.cds.split(",")]
fas = {}
def protein(r):
    cds = gc[r.L_gencode_tx]["cds"] if r.cds_route == "GENCODE" else orf[(r.chrom, r.strand, r.l_chain)]
    g = DATASETS[r.dataset]["genome"]; fas.setdefault(g, pysam.FastaFile(g))
    return translate(tx_seq(fas[g], r.chrom, cds, r.strand)).rstrip("*")

def ebi(seq):
    f = os.path.join(X, hashlib.sha1(seq.encode()).hexdigest() + ".tsv")
    if os.path.exists(f): return f
    body = urllib.parse.urlencode([("email", EMAIL), ("sequence", seq), ("goterms", "false"), ("pathways", "false")] + [("appl", a) for a in APPL]).encode()
    for attempt in range(6):          # EBI answers 400/5xx under load; back off and retry
        try:
            job = urllib.request.urlopen(urllib.request.Request(f"{API}/run", data=body), timeout=60).read().decode(); break
        except urllib.error.HTTPError as e:
            msg = e.read().decode()[:300]
            if attempt == 5: raise SystemExit(f"submit failed for {seq[:20]}... len {len(seq)}: {e.code} {msg}")
            time.sleep(30 * (attempt + 1))
    while True:
        time.sleep(15)
        st = urllib.request.urlopen(f"{API}/status/{job}", timeout=60).read().decode()
        if st == "FINISHED": break
        if st not in ("RUNNING", "QUEUED", "PENDING"): raise SystemExit(f"{job}: {st}")
    open(f, "w").write(urllib.request.urlopen(f"{API}/result/{job}/tsv", timeout=120).read().decode())
    return f

def hits(seq):
    f = ebi(seq)
    cols = ["acc", "md5", "len", "analysis", "sig", "desc", "start", "end", "score", "status", "date", "ipr", "ipr_desc"]
    if not os.path.getsize(f): return pd.DataFrame(columns=cols)
    t = pd.read_csv(f, sep="\t", header=None, dtype=str, engine="python", on_bad_lines="error").iloc[:, :13]   # 11-15 columns per line
    t.columns = cols[:t.shape[1]]
    t[["start", "end"]] = t[["start", "end"]].astype(int)
    return t

from concurrent.futures import ThreadPoolExecutor
prot = {r.Index: protein(r) for r in d.itertuples()}
seqs = sorted({x for i, r in zip(d.index, d.itertuples()) for x in (prot[i], prot[i][int(r.aa_lost):])})
print(len(d), "rows,", len(seqs), "distinct sequences", flush=True)
with ThreadPoolExecutor(5) as ex: list(ex.map(ebi, seqs))     # EBI allows up to 30 concurrent jobs; 5 avoids its 400s under load
rows = []
for r in d.itertuples():
    L = prot[r.Index]; a = int(r.aa_lost); S = L[a:]
    assert len(L) == int(r.aa_total) and S[:1] == "M", (r.gene, len(L), r.aa_total, S[:3])
    tl, ts = hits(L), hits(S)
    pf = tl[tl.analysis == "Pfam"]
    lost = [f"{x.desc} {x.start}-{x.end}" for x in pf.itertuples() if x.end <= a]
    part = [f"{x.desc} {x.start}-{x.end}" for x in pf.itertuples() if x.start <= a < x.end]
    kept = [f"{x.desc} {x.start}-{x.end}" for x in pf.itertuples() if x.start > a]
    def feat(t, kinds): return "; ".join(f"{x.sig} {x.start}-{x.end}" for x in t[t.analysis.isin(kinds)].itertuples())
    sp = lambda t: feat(t, ["SignalP_EUK"]) or "; ".join(f"{x.start}-{x.end}" for x in t[(t.analysis == "Phobius") & (t.sig == "SIGNAL_PEPTIDE")].itertuples())
    tm = lambda t: "; ".join(f"{x.start}-{x.end}" for x in t[t.analysis == "TMHMM"].itertuples())
    rows.append(dict(gene=r.gene, library=r.dataset.replace("ALS", "FTD"), site=r.peak, L_aa=len(L), S_aa=len(S), aa_lost=a, cds_route=r.cds_route,
                     pfam_lost="; ".join(lost), pfam_partial="; ".join(part), pfam_kept="; ".join(kept),
                     L_signal_peptide=sp(tl), S_signal_peptide=sp(ts), L_TM=tm(tl), S_TM=tm(ts)))
    print(rows[-1]["gene"], rows[-1]["library"], "lost:", rows[-1]["pfam_lost"], "| SP L/S:", rows[-1]["L_signal_peptide"], "/", rows[-1]["S_signal_peptide"],
          "| TM L/S:", rows[-1]["L_TM"], "/", rows[-1]["S_TM"], flush=True)
pd.DataFrame(rows).drop_duplicates(["gene", "library", "site"]).sort_values(["gene", "library"]).to_csv(os.path.join(O, "protein_domains.tsv"), sep="\t", index=False)
