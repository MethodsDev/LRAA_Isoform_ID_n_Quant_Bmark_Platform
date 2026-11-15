#!/usr/bin/env python3
"""
Supplementary Table 2: internal TSSs of strict ISM isoforms that imply an N-terminally truncated protein.

Reads __Benchmark_extras/ISM_internal_TSS/out/candidates.full.tsv.gz (LRAA v0.51.0), keeps the qualifying rows
(pile-up + cap + CAGE or refTSS, not a canonical TSS, >= 5 reads over the L intron S drops, annotated start lost), one
row per dataset x site with its primary L. Isoforms are named by splice hash, never by LRAA model ID.

Writes SuppTable2_ISM_internal_TSS.xlsx and .tsv. Run from this directory:  python3 make_ISM_supp_table.py
"""
import pandas as pd
d = pd.read_csv("../supp_analyses/ism_internal_tss/out/candidates.full.tsv.gz", sep="\t", low_memory=False)
q = d[d.core & d.protein_qualifying & ~d.s_at_canonical_tss & (d.l_up_junc_reads >= 5) & (d["class"] == "strict_ISM")]
q = q.sort_values(["primary_L", "n_site"], ascending=False).drop_duplicates(["dataset", "chrom", "strand", "peak"])
q = q.assign(library=q.dataset.replace({"ALS": "FTD"}), cap_pct=(100 * q.cap_frac).round(1), annotated_start=q.ev_annot_tx.astype(bool))
cols = {"gene": "gene", "library": "library", "regime": "regime", "chrom": "chrom", "strand": "strand", "peak": "site", "n_site": "reads_at_site", "cap_pct": "capped_pct",
        "cage_dist": "CAGE_dist_nt", "reftss_dist": "refTSS_dist_nt", "annotated_start": "GENCODE_start_shares_first_intron", "protein_class": "protein_class", "aa_total": "L_aa",
        "aa_remaining": "S_aa", "kozak": "kozak", "cds_route": "L_CDS_source", "L_gencode_tx": "L_GENCODE_transcript", "s_hash": "S_splice_hash", "l_hash": "L_splice_hash",
        "l_up_junc_reads": "reads_over_dropped_L_intron", "core_datasets": "libraries_with_core_site", "du_HC_events": "high_confidence_DEXSeq_events"}
t = q[list(cols)].rename(columns=cols).sort_values(["gene", "library"])
t["libraries_with_core_site"] = t.libraries_with_core_site.str.replace("ALS", "FTD")
t.to_csv("SuppTable2_ISM_internal_TSS.tsv", sep="\t", index=False)
t.to_excel("SuppTable2_ISM_internal_TSS.xlsx", index=False)
print(len(t), t.drop_duplicates(["chrom", "strand", "site"]).shape[0], t.gene.nunique())
