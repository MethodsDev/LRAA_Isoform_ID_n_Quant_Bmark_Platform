# LRAA Isoform ID and Quantification Benchmark Platform

This repository is a runnable skeleton for the LRAA paper benchmark framework.
It keeps the orchestration, registries, notebook runner, and notebook templates
for the three benchmark regimes:

- `QUANT_ONLY/`: quantification-only benchmarking
- `DENOVO_ID/`: reference-free isoform identification benchmarking
- `REF_Guided/`: reference-guided isoform identification benchmarking

Input files and generated benchmark outputs are intentionally not tracked here.
Place externally hosted inputs into each dataset's `reference_data/` directory
and each sample's `raw_prog_results/` directory before running.

Reference files can be symlinked from the paper analysis inventory with:

```bash
python3 scripts/install_references.py link
```

Proxy Oarfish quantification files are also linked into the corresponding
sample `raw_prog_results/` directories.

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

Install the Python dependencies, then populate the required external inputs:

```bash
python3 -m pip install -r requirements.txt
```

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
