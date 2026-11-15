# Library-level calibration of the TSS cap signature

Started 2026-10-29. The report is `cap_calibration.md` (knitted from `cap_calibration.Rmd`).

**Question.** Does a library carry the untemplated-G cap signature, and with which baseline?
LRAA v0.46.0 answered this per splice graph, and in the genome-wide FTD run the answer failed
locally (a deep uncapped truncation pile-up switched a 10 Mb chunk strand off; small cell
clusters had too few reads). The replacement decides once per library from the distribution of
per-peak capped fractions, which is bimodal when the signature is present.

**Libraries.** FTD cortex nuclei and PBMCs (sc Kinnex 10x 3'), BT474, K562, HG002 and UHRR
(bulk PacBio Kinnex), MCF7 and A549 (SG-NEx ONT PCS109), and six without the signature: MORFs
(PacBio, ONT; 5' ends trimmed), the SIRV spike-ins of three bulk libraries, and the Arabidopsis
simulation. Inputs are listed in the `Makefile`.

**Method.** `scripts/export_peaks.py` runs LRAA's own pre-pass module
(`pylib/TSS_cap_calibration.py`: `scan`, `decide`, `write_report`), so the read rules, peak
calling and decision are LRAA's. It writes, per library, `out/<lib>.peaks.tsv.gz` (one row per
TSS peak: informative spliced reads, capped reads, G-run counts, CAGE overlap) and LRAA's
calibration report `out/<lib>.TSS_cap_calibration.{json,hist.tsv,svg}`.
`scripts/event_site_caprates.py` measures the FTD event-pair sites (`ext/`: the v0.45.0
neuron-favours-downstream TSS switches, before any cap rule). `scripts/offpeak_strata.py`
shows what a "non-peak background" contains.

```
make LRAA_HOME=<LRAA checkout with pylib/TSS_cap_calibration.py>   # tables, one library at a time
make report LRAA_HOME=...                                          # also knit the report
```

Each scan reads only 5' ends, no assembly (FTD, 48M reads: ~2.5 min on 12 workers), under a
30 GB memory cap. The pre-pass is in LRAA since v0.47.0, so the default `LRAA_HOME` (the shared
checkout) works.

**Not pursued now:** a fitted two-component beta-binomial mixture, for a model-based library
decision and per-site posteriors in place of the fixed 10% threshold (report section
"Possible extension"; LRAA `docs/TSS_cap_calibration.md`).

**Per cell type** (FTD, PBMC; fixed RefQuantOnly clusters): `scripts/peak_mode_by_cluster.py`
(reads at the library's uncapped vs capped peaks), `scripts/depth_matched_peaks.py` (peaks called
within each cell type at full depth and with every cell type thinned to the same read count) and
`scripts/capped_share_by_celltype.py` (share of each cell type's reads starting at a capped 5'
end, peak-free). Reported in `../false_tss_examples/`.
