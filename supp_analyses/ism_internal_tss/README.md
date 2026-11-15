# ISM isoforms with internal TSSs (LRAA v0.51.0)

Are there isoforms shaped like 5' fragments of a longer isoform (S, a chain suffix of L) whose start
is a real internal promoter, implying a shorter protein? Report: `ISM_internal_TSS.md` / `.html`,
knitted from `ISM_internal_TSS.Rmd`. Every pair: `candidates.tsv` (219 MB). Reads existing GTFs and
BAMs only; no LRAA run. `make` rebuilds it all (about 1 h, memory-capped); `make report` only knits.

Inputs: v0.51.0 PBMC and ALS GTFs (de novo, ref-guided; basic, cluster-guided), bulk PacBio cell-line
GTFs (BT474, K562, HG002, UHRR), their BAMs, GENCODE v47 (the SQLIKE reference), FANTOM5 CAGE and
refTSS v4.1 (from `../isoform_end_accuracy/data`), DEXSeq site events, UniProt (`ext/uniprot`).

| step | script | output |
|---|---|---|
| S/L pairs per dataset; since 2026-11-11 S must end within `END3_TOL` (50) nt of an L model of the same model set (`FINDINGS_ism_survival.md`) | `bin/enumerate_pairs.py` | `out/pairs/` |
| read starts, cap, junctions (one BAM pass, LRAA read rules) | `bin/scan_reads.py` | `out/scan/` |
| evidence per S TSS, CAGE background | `bin/score_sites.py` | `out/scored/*.sites.tsv.gz` |
| protein consequence per pair (`bin/protein_lib.py`) | `bin/protein.py` | `out/scored/*.pairs.tsv.gz` |
| ORF route for L without a GENCODE chain (pyTD v0.4.0, shared `orf_settings.conf`) | `bin/run_orfs.py` | `out/orf/orfs.tsv` |
| table, tiers, primary L | `bin/assemble.py`, `bin/summarize.py` | `out/candidates.full.tsv.gz`, `candidates.tsv`, `out/tier_counts*.tsv`, `out/atlas_background.tsv` |
| FTD figure example by rule (core, N-truncated, FSM reads for L and S; since 2026-11-08 also a GENCODE L CDS and >= 10 reads fitting L from its TSS) | `bin/select_ftd_example.py` | `out/examples/ftd_example_candidates.tsv` |
| tool comparison (Supp Fig c, 2026-11-10): per bulk library, a pool of every tool's de novo models + GENCODE scored by the same rules; then each tool's own strict ISM isoforms (S from the tool, L from the tool or GENCODE), starts judged supported / truncation-like / not judged. `bin/tool_recall.py` (shared pool targets) is the earlier version, kept but no longer cited | `bin/run_pool.sh` (`bin/make_pool_gtf.py`, then the per-dataset steps on each pool), `bin/tool_own_ism.py` (`make tools`) | `out/pool/`, `out/scored/*_pool.*`, `out/tool_own_ism.tsv`, `out/tool_own_ism.pooled.tsv` |
| 10 examples + the FTD example, UniProt | `bin/make_examples.py`, `bin/uniprot_features.py` | `out/examples/` (incl. `read_starts_locus.tsv`), `out/next_candidates.tsv` |

Key columns of `candidates.tsv`: `class` (strict_ISM, TSS_in_L_intron, TSS_past_L_exon); `peak` (site);
`n_site`, `n_inf`, `n_cap`, `cap_frac`; `cage_dist`, `reftss_dist`; `ev_*` (each evidence);
`core` (pile-up + cap + atlas); `primary_L` (one row per S TSS); `protein_class`, `aa_lost`,
`aa_remaining`, `kozak`; `qualifying_strict`; `fragile_share` (PBMC monocytes, ALS neurons);
`core_datasets`. Rules are in the report's Definitions.
