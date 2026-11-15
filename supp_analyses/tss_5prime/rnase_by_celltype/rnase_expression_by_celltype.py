#!/usr/bin/env python3
"""RNase and RNA-decay gene expression by cell type (mean of the cell type's clusters' TPM).

Inputs: the per-cluster gene TPM matrix of an LRAA single-cell run (cluster_gene_TPM_matrix,
copied into ext/), and that run's withGeneSymbols GTF (gene_id "SYMBOL^lraa_gene_id") to map
LRAA gene ids to symbols; TPMs of all LRAA genes with one symbol are summed.

usage: rnase_expression_by_celltype.py DATASET TPM_TSV_GZ (WITH_SYMBOLS_GTF | GENE_SYMBOL_TSV) OUT_TSV
"""
import sys, re, gzip, collections

ds, tpm_f, gtf_f, out = sys.argv[1:5]
CT = {"FTD": {"0": "Oligodendrocyte", "3": "Oligodendrocyte", "1": "Excitatory neuron", "8": "Excitatory neuron",
              "4": "Inhibitory neuron", "6": "Inhibitory neuron", "2": "Astrocyte", "7": "Astrocyte", "12": "Astrocyte",
              "5": "Microglia", "9": "OPC", "10": "Vascular / meningeal", "11": "Vascular / meningeal", "13": "Vascular / meningeal"},
      "PBMC": {"1": "Monocyte / DC", "11": "Monocyte / DC", "13": "Monocyte / DC", "6": "B", "9": "B",
               "3": "CD4 T", "7": "CD4 T", "0": "CD4 T", "5": "CD4 T", "12": "CD4 T",
               "2": "CD8 T / NK", "4": "CD8 T / NK", "8": "CD8 T / NK", "10": "CD8 T / NK"}}[ds]
# gene -> class
GENES = {"RNASE1": "RNase A family (secreted)", "RNASE2": "RNase A family (secreted)", "RNASE3": "RNase A family (secreted)",
         "RNASE4": "RNase A family (secreted)", "ANG": "RNase A family (secreted)", "RNASE6": "RNase A family (secreted)",
         "RNASE7": "RNase A family (secreted)", "RNASET2": "lysosomal RNase", "PLD3": "lysosomal 5'-3' exonuclease",
         "PLD4": "lysosomal 5'-3' exonuclease", "RNASEL": "antiviral endoribonuclease", "ZC3H12A": "Regnase endoribonuclease",
         "ZC3H12D": "Regnase endoribonuclease", "ERN1": "IRE1 (RIDD) endoribonuclease", "XRN1": "5'-3' decay",
         "XRN2": "5'-3' decay (nuclear)", "DCP2": "decapping", "DCP1A": "decapping", "EXOSC10": "exosome (nuclear)",
         "DIS3": "exosome", "DIS3L2": "3'-5' decay", "CNOT7": "deadenylation", "PAN2": "deadenylation",
         "SMG6": "NMD endonuclease", "PNPT1": "3'-5' exoribonuclease"}

sym = {}
if gtf_f.endswith(".gtf"):
    for l in open(gtf_f):
        if "\ttranscript\t" in l:
            x = re.search(r'gene_id "([^"^]+)\^([^"]+)"', l)
            if x:
                sym[x.group(2)] = x.group(1)
else:  # lraa_gene_id <tab> gene_symbol
    for l in open(gtf_f):
        f = l.rstrip("\n").split("\t")
        if len(f) == 2 and f[0] != "lraa_gene_id":
            sym[f[0]] = f[1]
with gzip.open(tpm_f, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")[1:]
    clusters = [re.search(r"\.(\d+)\.LRAA", h).group(1) for h in hdr]
    tot = collections.defaultdict(lambda: [0.0] * len(hdr))
    for line in fh:
        f = line.rstrip("\n").split("\t")
        s = sym.get(f[0])
        if s in GENES:
            for i, v in enumerate(f[1:]):
                tot[s][i] += float(v)
with open(out, "wt") as ofh:
    ofh.write("dataset\tgene\tclass\tcell_type\tclusters\tmean_TPM\n")
    for g in GENES:
        if g not in tot:
            continue
        for ct in sorted(set(CT.values())):
            idx = [i for i, c in enumerate(clusters) if CT.get(c) == ct]
            ofh.write(f"{ds}\t{g}\t{GENES[g]}\t{ct}\t{len(idx)}\t{sum(tot[g][i] for i in idx) / len(idx):.2f}\n")
