# LRAA Isoform ID and Quantification Benchmark Platform

This repository is a runnable skeleton for the LRAA paper benchmark framework.
It keeps the orchestration, registries, notebook runner, and notebook templates
for the three benchmark regimes:

- `QUANT_ONLY/`: quantification-only benchmarking
- `DENOVO_ID/`: reference-free isoform identification benchmarking
- `REF_Guided/`: reference-guided isoform identification benchmarking

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
QUANT_ONLY/                   quantification-only regime
DENOVO_ID/                            de novo isoform ID regime
REF_Guided/                           reference-guided isoform ID regime
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
