# LRAA Isoform ID and Quantification Benchmark Platform

This repository is a runnable skeleton for the LRAA paper benchmark framework.
It keeps the orchestration, registries, notebook runner, and notebook templates
for the three benchmark regimes:

- `QUANT_ONLY/`: quantification-only benchmarking
- `DENOVO_ID/`: reference-free isoform identification benchmarking
- `REF_Guided/`: reference-guided isoform identification benchmarking

Input files and generated benchmark outputs are intentionally not tracked here.
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
