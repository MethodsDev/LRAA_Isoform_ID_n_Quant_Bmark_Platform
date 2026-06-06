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

Reference-data symlinks can be created from the benchmark reference inventory:

```bash
python3 scripts/install_references.py link
```

Use `--dry-run` to preview changes. The `install` subcommand also unpacks local
`<dataset>.ref_data.tar.gz` inventory payloads before linking. Proxy Oarfish
quantification files are also linked into the matching sample
`raw_prog_results/` directories unless `--skip-raw-proxy-quants` is used.
