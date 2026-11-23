# LRAA Isoform ID and Quantification Benchmark Platform

This repository is a runnable skeleton for the LRAA paper benchmark framework.
It keeps the orchestration, registries, notebook runner, and notebook templates
for the three benchmark regimes:

- `QUANT_ONLY/`: quantification-only benchmarking
- `DENOVO_ID/`: reference-free isoform identification benchmarking
- `REF_Guided/`: reference-guided isoform identification benchmarking

It also holds the paper's supplementary analyses and tables (`supp_analyses/`, `supp_tables/`), the
notebooks that generated the SIRV, MORF and cell-line data (`data_generation/`), and the two single-cell
studies as submodules (`LRAA_paper-PBMCs_study/`, `LRAA_paper-FTD_study/`). All LRAA results are from release
**v0.51.0** (tag `LRAA_v0.51.0`); in the benchmarks LRAA is the registry entry `LRAA-v0.51.0`.

External inputs and intermediate benchmark outputs are intentionally not tracked here.
Reference inputs and deposited comparator predictions are published in
[Zenodo record 20719930](https://zenodo.org/records/20719930) and can be
installed directly with:

```console
make install-inputs
```

The installers download and verify the reference and prediction archives,
unpack them under `.benchmark_inputs/`, and create the required symlinks.
Downloads are cached and partial downloads resume on the next run. Existing
non-matching files in `raw_prog_results/` are never replaced unless `--force`
is explicitly passed to `scripts/install_predictions.py`.

To install only one dataset (useful because several archives are large), run:

```console
python3 scripts/install_references.py --dataset SIRVs install
```

Use `make install-references` or `make install-predictions` when only one input
class is needed. Prediction regimes can also be selected individually, for
example `python3 scripts/install_predictions.py --regime QUANT_ONLY`.

## Layout

```text
benchmarking/                         shared Python runner, parsers, notebooks
QUANT_ONLY/                           quantification-only regime
DENOVO_ID/                            de novo isoform ID regime
REF_Guided/                           reference-guided isoform ID regime
supp_analyses/                        supplementary analyses (5' completeness, ISM internal TSS, proteoforms,
                                      site atlas, SQLIKE categories)
supp_tables/                          builders of the Supplementary Tables
data_generation/                      notebooks that generated the SIRV, MORF and cell-line data
LRAA_paper-PBMCs_study/               submodule: PBMC single-cell study
LRAA_paper-FTD_study/                 submodule: FTD single-nucleus study
```

Each regime has a `tool_registry.yaml`, a top-level `Makefile`, per-dataset
`Makefile`s, and per-sample `run*.sh` scripts. The runner discovers the nearest
`tool_registry.yaml` by walking up from the current sample directory.

## Paper benchmark figures

Figures 2 and 3 in the LRAA manuscript are composites assembled from selected
rows of the benchmark outputs below. The **report** links open the complete
knitted analyses on GitHub; the **PDF** links open the corresponding vector
plots. The PDFs retain all samples in each dataset, while the manuscript uses
the rows identified in the tables.

### Figure 2: benchmarks with defined truth sets

| Manuscript row | Quantification only | Reference-guided identification | De novo identification |
|---|---|---|---|
| Mouse, 1.6% error (`e016`) | [report](QUANT_ONLY/mouse_sim/Mouse_sim.QUANT_ONLY.summary.md) · [PDF](QUANT_ONLY/mouse_sim/mouse_sim_quant_summary_barplot.latest_tools.pdf) | [report](REF_Guided/mouse_sim/Mouse_ref_guided_eval.md) · [PDF](REF_Guided/mouse_sim/mouse_sim_summary.latest_tools.pdf) | [report](DENOVO_ID/mouse_sim/Mouse_denovo_eval.md) · [PDF](DENOVO_ID/mouse_sim/mouse_sim_summary.latest_tools.pdf) |
| Arabidopsis, 1.6% error (`e016`) | [report](QUANT_ONLY/arabidopsis_sim/Arabidopsis_sim.QUANT_ONLY.summary.md) · [PDF](QUANT_ONLY/arabidopsis_sim/arabidopsis_sim_quant_summary_barplot.latest_tools.pdf) | [report](REF_Guided/arabidopsis_sim/Arabidopsis_ref_guided_eval.md) · [PDF](REF_Guided/arabidopsis_sim/Arabidopsis_sim_summary.latest_tools.pdf) | [report](DENOVO_ID/arabidopsis_sim/Arabidopsis_denovo_eval.md) · [PDF](DENOVO_ID/arabidopsis_sim/Arabidopsis_sim_summary.latest_tools.pdf) |
| SIRVs, E1 mix | [report](QUANT_ONLY/SIRVs/SIRVs.QUANT_ONLY.summary.md) · [PDF](QUANT_ONLY/SIRVs/SIRVs_sim_quant_summary.latest_tools.pdf) | [report](REF_Guided/SIRVs/SIRVs_ref_guided_eval.md) · [PDF](REF_Guided/SIRVs/SIRVs_sim_summary.latest_tools.pdf) | [report](DENOVO_ID/SIRVs/SIRVs_denovo_eval.md) · [PDF](DENOVO_ID/SIRVs/SIRVs_sim_summary.latest_tools.pdf) |
| MORFs, PacBio | [report](QUANT_ONLY/MORFs/Morf2.QUANT_ONLY.summary.md) · [PDF](QUANT_ONLY/MORFs/combined_morf_quant_plot.latest_tools.pdf) | [report](REF_Guided/MORFs/Morfs2_ref_guided_eval.md) · [PDF](REF_Guided/MORFs/MORFs_refGuided_plot.latest_tools.pdf) | [report](DENOVO_ID/MORFs/Morf2_denovo_eval.md) · [PDF](DENOVO_ID/MORFs/combined_morf_denovo_plot.latest_tools.pdf) |

Figure 2 uses mean relative difference for quantification and F1, recall, and
precision for identification. Each linked Mouse and Arabidopsis PDF contains
all four simulated error rates; each SIRV PDF contains E0, E1, and E2; and each
MORF PDF contains both PacBio and ONT results.

### Figure 3: whole-transcriptome benchmarks

| Manuscript panels | Quantification only | Reference-guided identification | De novo identification |
|---|---|---|---|
| **3a–d**, PacBio Kinnex: BT474, K562, GM24385/HG002, and UHRR | [report](QUANT_ONLY/CellLines/examine_expr_correlations.md) · [PDF](QUANT_ONLY/CellLines/summary_CellLines_corr_plots.latest_tools.pdf) | [report](REF_Guided/CellLines/eval_cell_line_RefGuided_ID.md) · [PDF](REF_Guided/CellLines/Cell_lines_summary.latest_tools.pdf) | [report](DENOVO_ID/CellLines/eval_cell_line_denovo_ID.md) · [PDF](DENOVO_ID/CellLines/Cell_lines_summary.latest_tools.pdf) |
| **3e–f**, ONT: A549 and MCF7 | [report](QUANT_ONLY/SG-NEx/examine_expr_correlations.SGNEx.md) · [PDF](QUANT_ONLY/SG-NEx/SGNex_expr_cor_plot.latest_tools.pdf) | [report](REF_Guided/SG-NEx/eval_cell_line_RefGuided_ID.md) · [PDF](REF_Guided/SG-NEx/SG-NEx_cell_lines_stacked.latest_tools.pdf) | [report](DENOVO_ID/SG-NEx/eval_cell_line_denovo_ID.md) · [PDF](DENOVO_ID/SG-NEx/SG-NEx_cell_lines_denovo_stacked.latest_tools.pdf) |

Figure 3 uses pairwise Pearson correlations of isoform-level `log(TPM + 1)`
for quantification and proxy-truth F1, recall, and precision for identification.

### Supplementary benchmark figures

- **Supplementary Figure 2a–b** uses the remaining error-rate rows in the Mouse
  and Arabidopsis summaries linked under Figure 2. **Supplementary Figure 2c**
  uses the dedicated F1-by-error plots for [Mouse reference-guided](REF_Guided/mouse_sim/mouse_sim_F1_by_error_rate.latest_tools.pdf),
  [Mouse de novo](DENOVO_ID/mouse_sim/mouse_sim_F1_by_error_rate.latest_tools.pdf),
  [Arabidopsis reference-guided](REF_Guided/arabidopsis_sim/Arabidopsis_sim_F1_by_error_rate.latest_tools.pdf),
  and [Arabidopsis de novo](DENOVO_ID/arabidopsis_sim/Arabidopsis_sim_F1_by_error_rate.latest_tools.pdf).
- **Supplementary Figure 3** uses the E0 and E2 rows in the SIRV summaries
  linked under Figure 2.
- **Supplementary Figure 4b** uses the ONT rows in the MORF summaries linked
  under Figure 2. Supplementary Figure 4a is the MORF barcoding schematic and
  is not generated by these benchmark reports.

## Basic Usage

Install the Python dependencies and benchmark inputs:

```console
python3 -m pip install -r requirements.txt
make install-inputs
```

Run all three benchmark regimes, build their aggregate result tables, and knit
the 18 dataset summary reports as GitHub documents from the repository root:

```console
make
make summarize
```

To re-knit the reports without regenerating the aggregate tables, run
`make render-reports`.

Run a full regime:

```bash
cd QUANT_ONLY
make
make summarize
```

Run one dataset:

```bash
cd DENOVO_ID/arabidopsis_sim
make
make summarize
```

Run one sample:

```bash
cd REF_Guided/arabidopsis_sim/arabidopsis_isoseqsim_e000
./run_ref_guided_bmarking.sh
```

See `BENCHMARKING.md` for the registry schema, run-script flags, and the
benchmarking workflow.

## Supplementary analyses and tables

The supplementary analyses read the study submodules and the cross-dataset Zenodo bundle. With both in place,
each report is knit from its own directory (`Rscript -e "rmarkdown::render('<report>.Rmd')"`) and each
figure script runs there (`Rscript SuppFig_<name>.R`). The proteoform notebook runs with
`make -C supp_analyses/proteoforms notebook REF_GTF=<GENCODE v47 GTF>`; its ORFs are in the bundle
(`supp_analyses/proteoforms/orfs/`) and `run_orfs.sh` predicts them again. Supplementary Tables 1-2 are built
by `python3 make_dexseq_supp_table.py` and `python3 make_ISM_supp_table.py` in `supp_tables/`; Supplementary
Tables 3-4 (benchmark accuracy values, sequencing libraries) are provided as built.

| Paper item | Panels | What | File(s) in this repo |
|---|---|---|---|
| Supplementary Figure 6 | a-d | image (rebuilt here) | `supp_analyses/categories/SuppFig_categories.pdf` |
| Supplementary Figure 6 | a-d | generator | `supp_analyses/categories/SuppFig_categories.R` |
| Supplementary Figure 7 | all | image (rebuilt here) | `supp_analyses/proteoforms/SuppFig_proteoforms_PBMC.pdf` |
| Supplementary Figure 7 | all | generator (filtered) | `supp_analyses/proteoforms/LRAA_proteoform_analysis_PBMC_ALS_20260528.ipynb` |
| Supplementary Figure 7 | all | generator | `supp_analyses/proteoforms/orf_settings.conf` |
| Supplementary Figure 8 | all | image (rebuilt here) | `supp_analyses/proteoforms/SuppFig_proteoforms_FTD.pdf` |
| Supplementary Figure 8 | all | generator (filtered) | `supp_analyses/proteoforms/LRAA_proteoform_analysis_PBMC_ALS_20260528.ipynb` |
| Supplementary Figure 8 | all | generator | `supp_analyses/proteoforms/orf_settings.conf` |
| Supplementary Figure 9 | a-c | image (rebuilt here) | `supp_analyses/site_atlas/SuppFig_site_atlas.pdf`<br>`supp_analyses/site_atlas/SuppFig_site_atlas.values.tsv` |
| Supplementary Figure 9 | a-c | generator | `supp_analyses/site_atlas/SuppFig_site_atlas.R` |
| Supplementary Figure 10 | a,c | report (rebuilt here) | `supp_analyses/tss_5prime/truncation_check/tss_truncation_check.md`<br>+ 5 knitr figure files |
| Supplementary Figure 10 | a,c | generator | `supp_analyses/tss_5prime/truncation_check/Makefile`<br>`supp_analyses/tss_5prime/truncation_check/control_gene_tss_reach.py`<br>`supp_analyses/tss_5prime/truncation_check/tss_truncation_check.Rmd` |
| Supplementary Figure 10 | a-e | image (rebuilt here) | `supp_analyses/tss_5prime/SuppFig_TSS_5prime.pdf`<br>`supp_analyses/tss_5prime/SuppFig_TSS_5prime.values.tsv` |
| Supplementary Figure 10 | a-e | generator | `supp_analyses/tss_5prime/SuppFig_TSS_5prime.R` |
| Supplementary Figure 10 | b,d,e | report (rebuilt here) | `supp_analyses/tss_5prime/completeness/cap_calibration/cap_calibration.md`<br>`supp_analyses/tss_5prime/completeness/tss_5prime_completeness.md`<br>+ 14 knitr figure files |
| Supplementary Figure 10 | b,d,e | generator | `supp_analyses/tss_5prime/completeness/` (15 files) |
| Supplementary Figure 11 | a-g | image (rebuilt here) | `supp_analyses/ism_internal_tss/SuppFig_ISM_internal_TSS.pdf`<br>`supp_analyses/ism_internal_tss/SuppFig_ISM_internal_TSS.values.tsv` |
| Supplementary Figure 11 | a-g | generator | `supp_analyses/ism_internal_tss/SuppFig_ISM_internal_TSS.R` |
| Supplementary Figure 11 | - | report (rebuilt here) | `supp_analyses/ism_internal_tss/ISM_internal_TSS.md`<br>+ 18 knitr figure files |
| Supplementary Figure 11 | - | generator | `supp_analyses/ism_internal_tss/` (22 files) |
| Supplementary Table 1 | - | generator | `supp_tables/make_dexseq_supp_table.py` |
| Supplementary Table 1 | - | table (rebuilt here) | `supp_tables/SuppTable1_differential_usage.xlsx`<br>`supp_tables/differential_usage_all_events.zip` |
| Supplementary Table 2 | - | report (rebuilt here) | `supp_analyses/ism_internal_tss/ISM_internal_TSS.md`<br>+ 18 knitr figure files |
| Supplementary Table 2 | - | generator | `supp_analyses/ism_internal_tss/` (22 files) |
| Supplementary Table 2 | - | generator | `supp_tables/make_ISM_supp_table.py` |
| Supplementary Table 2 | - | table (rebuilt here) | `supp_tables/SuppTable2_ISM_internal_TSS.tsv`<br>`supp_tables/SuppTable2_ISM_internal_TSS.xlsx` |
| Supplementary Table 3 | - | table | `supp_tables/SuppTable3_benchmark_metrics.xlsx` |
| Supplementary Table 4 | - | table | `supp_tables/SuppTable4_libraries.xlsx` |

Cross-dataset Zenodo bundle:

| Bundle | Used by | Contents | Files |
|---|---|---|---|
| cross | Supplementary Figure 10 | tss5prime_completeness_data | `supp_analyses/tss_5prime/completeness/` (79 files) |
| cross | Supplementary Figure 10 | tss_truncation_lengths | `supp_analyses/tss_5prime/truncation_check/ALS.lengths.tsv.gz`<br>`supp_analyses/tss_5prime/truncation_check/PBMC.lengths.tsv.gz` |
| cross | Supplementary Figure 10 | tss_truncation_reads | `supp_analyses/tss_5prime/truncation_check/ALS.reads.tsv.gz`<br>`supp_analyses/tss_5prime/truncation_check/PBMC.reads.tsv.gz`<br>`supp_analyses/tss_5prime/truncation_check/control_genes.anchoring.tsv`<br>`supp_analyses/tss_5prime/truncation_check/control_genes.reads.tsv.gz`<br>`supp_analyses/tss_5prime/truncation_check/control_genes.sites.tsv` |
| cross | Supplementary Figure 11, Supplementary Table 2 | ism_data | `supp_analyses/ism_internal_tss/out/` (21 files) |
| cross | Supplementary Figure 11 | ism_example_reads | `supp_analyses/ism_internal_tss/out/examples/reads.tsv.gz` |
| cross | Supplementary Figure 7, Supplementary Figure 8 | sfig_proteoforms_orfs | `supp_analyses/proteoforms/orfs/FTD/LP00714.refGuided.scg.cdna.fasta.transdecoder.genome.gff3`<br>`supp_analyses/proteoforms/orfs/FTD/LP00714.refGuided.scg.cdna.fasta.transdecoder.pep`<br>`supp_analyses/proteoforms/orfs/FTD/ORFS_VERSION.txt`<br>`supp_analyses/proteoforms/orfs/PBMC/ORFS_VERSION.txt`<br>`supp_analyses/proteoforms/orfs/PBMC/PBMCs_pbio.LRAA.refGuided.scg.cdna.fasta.transdecoder.genome.gff3`<br>`supp_analyses/proteoforms/orfs/PBMC/PBMCs_pbio.LRAA.refGuided.scg.cdna.fasta.transdecoder.pep` |
