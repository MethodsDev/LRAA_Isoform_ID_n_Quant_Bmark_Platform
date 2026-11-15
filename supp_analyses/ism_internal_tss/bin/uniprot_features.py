#!/usr/bin/env python3
"""UniProt features lost with the N-terminus, for the top examples (+ next candidates). usage: uniprot_features.py

S's protein is L's from the first downstream in-frame ATG (L residue aa_lost + 1). L's protein comes from Ensembl (protein_id,
GENCODE route only), and S is placed on the reviewed canonical UniProt sequence by its first 25 residues, so an L that is a
non-canonical isoform is still mapped (2026-11-11; before, only an L of canonical length was reported, and only Domain,
Transmembrane, Signal, Transit peptide and Propeptide features were read). A feature ending before S's start is lost; one
spanning it is partial."""
import json, urllib.request, urllib.parse, time, os
import pandas as pd
here = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(here, "..", "out"); X = os.path.join(here, "..", "ext", "uniprot")
TYPES = ("Signal", "Transit peptide", "Propeptide", "Transmembrane", "Domain", "Repeat", "Region", "Coiled coil", "Motif")
ex = pd.concat([pd.read_csv(f"{O}/examples/examples.tsv", sep="\t"), pd.read_csv(f"{O}/next_candidates.tsv", sep="\t")])

def fetch(url, path):
    for attempt in range(5):
        if os.path.exists(path): break
        try:
            b = urllib.request.urlopen(url, timeout=60).read(); open(path, "wb").write(b); time.sleep(0.2)
        except OSError:
            if attempt == 4: raise
            time.sleep(10)
    return open(path).read()

rows = []
for r in ex.itertuples():
    q = urllib.parse.quote(f"gene_exact:{r.gene} AND organism_id:9606 AND reviewed:true")
    j = json.loads(fetch(f"https://rest.uniprot.org/uniprotkb/search?query={q}&format=json&size=1", f"{X}/{r.gene}.full.json")).get("results", [])
    if not j:
        rows.append(dict(gene=r.gene, uniprot="", length=None, length_match=False, s_start_in_uniprot=None, lost_features="", partial_features=""))
        continue
    e = j[0]; can = e["sequence"]["value"]; start = None
    if isinstance(r.protein_id, str) and r.protein_id.startswith("ENSP"):
        p = r.protein_id.split(".")[0]
        L = fetch(f"https://rest.ensembl.org/sequence/id/{p}?content-type=text/plain", f"{X}/{p}.txt").strip()
        i = can.find(L[int(r.aa_lost):][:25])
        start = i + 1 if i >= 0 else None
    lost, part = [], []
    for ft in e.get("features", []):
        if ft["type"] not in TYPES or start is None: continue
        a, b = ft["location"]["start"]["value"], ft["location"]["end"]["value"]
        if a is None or b is None or a >= start: continue
        (lost if b < start else part).append(f'{ft["type"]}{(" " + ft["description"]) if ft.get("description") else ""} {a}-{b}')
    rows.append(dict(gene=r.gene, uniprot=e["primaryAccession"], length=len(can), length_match=len(can) == r.aa_total,
                     s_start_in_uniprot=start, lost_features="; ".join(lost), partial_features="; ".join(part)))
pd.DataFrame(rows).drop_duplicates().to_csv(f"{O}/uniprot_lost.tsv", sep="\t", index=False)
