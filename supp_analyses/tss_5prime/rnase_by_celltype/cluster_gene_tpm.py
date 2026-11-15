#!/usr/bin/env python3
"""Per-cluster gene TPM matrix from an LRAA single-cell run's final cluster quant tarball
(*.LRAA.final.cluster_quant.EXPRs.tar.gz: one quant.expr per cluster; gene TPM = sum of its isoforms' TPM).

usage: cluster_gene_tpm.py CLUSTER_QUANT_TAR_GZ OUT_TSV_GZ
"""
import sys, tarfile, gzip, os, collections

tar_f, out = sys.argv[1:3]
names, tpm = [], collections.defaultdict(dict)
with tarfile.open(tar_f) as tf:
    for m in sorted((m for m in tf.getmembers() if m.name.endswith(".quant.expr")), key=lambda m: m.name):
        col = os.path.basename(m.name)
        names.append(col)
        hdr = None
        for line in tf.extractfile(m).read().decode().splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if hdr is None:
                hdr = {c: i for i, c in enumerate(f)}
                continue
            g = f[hdr["gene_id"]]
            tpm[g][col] = tpm[g].get(col, 0.0) + float(f[hdr["TPM"]])
with gzip.open(out, "wt") as ofh:
    ofh.write("\t" + "\t".join(names) + "\n")
    for g in sorted(tpm):
        ofh.write(g + "\t" + "\t".join(f"{tpm[g].get(c, 0.0):.3f}" for c in names) + "\n")
