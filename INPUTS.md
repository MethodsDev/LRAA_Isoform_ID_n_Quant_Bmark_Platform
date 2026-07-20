# External Inputs

This skeleton does not track benchmark input data, reference annotations,
comparator outputs, or generated results.

Populate inputs as follows:

- `*/reference_data/`: reference GTFs, reduced-reference GTFs, and truth
  quantification tables used by the corresponding run scripts.
- `*/raw_prog_results/`: comparator tool outputs for a sample. Filenames must
  match the regex patterns in the regime's `tool_registry.yaml`.
- `*/processed_prog_results/`: generated parser outputs. Leave empty before
  running; the benchmark runner writes these files.

To see the exact expected filenames for a dataset, inspect its `run*.sh` files
for `--truth_gtf`, `--truth_quant`, and `--truth_reduced_gtf` arguments. To see
accepted comparator output names, inspect the regime-level `tool_registry.yaml`.

Reference data and deposited comparator predictions are downloaded from
[Zenodo record 20719930](https://zenodo.org/records/20719930), verified against
the deposited checksums, unpacked into the ignored `.benchmark_inputs/` cache,
and linked into the benchmark tree with:

```console
make install-inputs
```

The complete downloaded and unpacked cache currently occupies approximately
53 GB. Allow at least 75 GB of free space for the cache, processed benchmark
outputs, notebooks, and aggregate tables. The cache can be reused by later
runs and is excluded from Git.

Use `python3 scripts/install_references.py --dataset DATASET install` to fetch a
single dataset. Downloads are cached, and interrupted `.part` downloads resume
on the next run. Use `--dry-run` to preview link changes after installation.
Proxy Oarfish quantification files are also linked into the matching sample
`raw_prog_results/` directories unless `--skip-raw-proxy-quants` is used.

The `link` subcommand and `--inventory-dir` can still be used with an existing
local inventory without contacting Zenodo.

Use `make install-references` or `make install-predictions` to install only one
input class. Prediction archives can be selected by regime with
`python3 scripts/install_predictions.py --regime QUANT_ONLY`.
