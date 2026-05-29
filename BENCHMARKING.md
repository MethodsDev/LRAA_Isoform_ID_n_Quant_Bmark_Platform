# Benchmarking framework

This document describes how the long-read RNA-seq benchmarking system in this
repo is organized after the registry-driven refactor, and walks through how to
add a new tool, a new version of an existing tool, or a new dataset.

## Three benchmarking regimes

The benchmarking is split into three regimes, each living in its own top-level
directory and using the same machinery:

| Regime | Top-level dir | What it measures | `--analysisType` |
|---|---|---|---|
| Quantification only | `QUANT_ONLY_outputs/` | Expression accuracy: tools quantify against the reference annotation | `quant_only` (or `quant_only_no_truthset`) |
| De novo ID | `DENOVO_ID/` | Isoform discovery accuracy: tools build transcripts without an annotation | `ref_free` |
| Ref-guided ID | `REF_Guided/` | Both quantification + novel discovery, with the annotation as a starting point | `ref_guided` |

The same set of comparator tools (Bambu, ESPRESSO, FLAIR, FLAMES, IsoQuant,
IsoSeq, Isosceles, LRAA, Mandalorion, Oarfish, StringTie, TALON) shows up across
the regimes, plus LRAA, with the subset that's relevant to each mode (e.g.
Oarfish is a quantifier-only tool so it doesn't appear in DENOVO_ID).

## Architecture

```
                                                     +-----------------+
                                                     |   tool_registry |
                                                     |     .yaml       |
                                                     +--------+--------+
                                                              |
                  one-shot deposit                            v
  +----------------+   ------>    +-------------------+   parser     +------------------+
  | comparator     |              |  raw_prog_results |  --------->  | processed_prog_  |
  |  outputs       |              |       /           |              |  results/        |
  | (per family)   |              |  one per sample   |              |  <name>.tsv      |
  +----------------+              +-------------------+              +--------+---------+
                                                                              |
                                                                       runner |
                                                                              v
                                                                     +------------------+
                                                                     | papermill exec   |
                                                                     |  template.<mode> |
                                                                     |    .ipynb        |
                                                                     +--------+---------+
                                                                              |
                                                                              v
                                                                     +------------------+
                                  per-dataset                        | per-sample TSVs  |
                                  +-------+                          | + bmark.ok       |
                                  | make  |  <-------- summarize --- |                  |
                                  +---+---+                          +------------------+
                                      |
                                      v
                              +------------------+
                              | aggregated TSVs  | <----------- Rmd report renders
                              +------------------+              from aggregated TSVs
```

A single declarative YAML registry per regime drives parsing, plotting, and
truth-set composition. Templates and runner code are tool-agnostic. Adding a
new tool or version means editing one YAML file.

### What lives where

The framework is packaged as one repository:

```
LRAA_Isoform_ID_n_Quant_Bmark_Platform/
+--- QUANT_ONLY_outputs/
|     +--- tool_registry.yaml                     <-- registry for this regime
|     +--- arabidopsis_sim/, CellLines/, MORFs/, SG-NEx/, SIRVs/, mouse_sim/, IGROV-1_sim/
|     |     +--- Makefile                         <-- dataset-level orchestration
|     |     +--- <sample>/
|     |     |     +--- run_quant_bmarking.sh     (or similar; the per-sample driver)
|     |     |     +--- raw_prog_results/         <-- inputs (untracked; deposited)
|     |     |     +--- processed_prog_results/   <-- parsed inputs (generated)
|     |     |     +--- bmark.ok                  <-- success sentinel (generated)
|     |     |     +--- <analysisType>.ipynb      <-- executed notebook (generated)
|     |     |     +--- *.tsv                     <-- per-sample metrics (generated)
|     |     +--- <dataset>.spearman_cor.tsv      <-- aggregated metric TSV
|     |     +--- ... etc
|     +--- __comparator_results_all_QUANT_ONLY/
|           +--- deposit_comparator_outputs.py    <-- routes new tool outputs into raw_prog_results/
+--- DENOVO_ID/                                   <-- mirrors QUANT_ONLY shape
|     +--- tool_registry.yaml
|     +--- ... (same per-dataset layout)
+--- REF_Guided/                                  <-- mirrors QUANT_ONLY shape
|     +--- tool_registry.yaml
|     +--- ... (same per-dataset layout)
+--- benchmarking/
      +--- pylib/QuantParser.py                   <-- registry-driven file dispatch
      +--- bmark_nb_runner.py                     <-- per-sample driver (loads registry, executes notebook)
      +--- template_notebooks/
      |     +--- template.quant-only.ipynb
      |     +--- template.quant-only-no-truthset.ipynb
      |     +--- template.denovo.ipynb
      |     +--- template.ref_reduced.ipynb
      +--- misc/aggregate_result_tables.py        <-- merges per-sample TSVs into one per dataset
```

The `benchmarking/` directory is regime-agnostic. The regime directories hold
their specific configs: registries, run scripts, Makefiles, deposited tool
outputs, and generated results.

## Tool registry schema

Each `tool_registry.yaml` is a list of entries, one per (tool, version)
combination. The parser walks entries in registry order and the first matching
entry wins, so list more-specific (versioned) patterns before any
less-specific fallback.

```yaml
tools:
  - name: IsoQuant-v3.13.0          # unique label; appears in plots and as
                                    # processed_prog_results/<name>.tsv
    family: IsoQuant                # used to group versions of the same tool
    quant_pattern: 'isoquant-v3\.13\.0\.IsoQuant\.counts\.tsv$'
                                    # regex matched against raw_prog_results/* basenames
    quant_id_col: 0                 # 0-indexed col holding the transcript_id
    quant_tpm_col: 1                # 0-indexed col holding the read count or TPM
                                    # (re-normalized to sum to 1M per file)
    quant_skip_rows: 0              # lines to skip before the header (default: 0)
    quant_no_header: false          # true if the file has no header (default: false)
    gtf_pattern: 'isoquant-v3\.13\.0\.IsoQuant\.gtf$'
                                    # (optional) regex for the tool's GTF
    gtf_source: own                 # "ref" -> use REF_gtf for intron derivation
                                    # "own" -> use the tool's own GTF
                                    # default: ref
    gtf_converter: null             # (optional) "flames_gff3" or "strip_blanks"
                                    # see "GTF converters" below
    color: blue                     # plot color (matplotlib named color or hex)
    display: true                   # include in plots and aggregated TSVs
                                    # (default: true)
    venn: true                      # contributes to consensus truth-set under
                                    # --Venn_mode (default: true)
```

### display vs venn

These two flags are independent and both default to `true`:

- `display` -- show this entry in plots and per-sample/aggregated TSVs.
- `venn` -- include this entry when seeding the consensus truth-set under
  `--Venn_mode` (used in DENOVO_ID and REF_Guided when no synthetic ground
  truth exists).

The convention now: when a tool family has multiple versions, only the newest
version gets `venn: true`; older versions are `venn: false` but stay
`display: true`. This avoids newer versions getting artificially rewarded for
agreeing with their own predecessors. All entries are still scored against the
truth set.

### gtf_source

- `gtf_source: ref` (default) means the entry's evaluation is keyed off
  `REF_gtf`. Correct for tools that quantify against the reference annotation
  (LRAA in quant-only mode, IsoQuant, Oarfish, FLAIR, StringTie, ESPRESSO,
  Bambu).
- `gtf_source: own` means the entry's own GTF is used for intron derivation.
  Correct for tools that emit private transcript IDs (FLAMES, IsoSeq,
  Mandalorion, Isosceles, TALON) and for everything in DENOVO_ID/REF_Guided.

### GTF converters

By default, GTF files matched by `gtf_pattern` are passed through to the
notebook unchanged. Two converters are available:

- `gtf_converter: flames_gff3` -- runs the FLAMES gff3-to-gtf converter
  (`misc/FLAMES_gff3_to_gtf_converter.py`) and
  writes the converted file to `processed_prog_results/<name>.gtf`.
- `gtf_converter: strip_blanks` -- copies the GTF through with blank lines
  removed. Required for LRAA's denovo and ref-guided GTFs (which contain
  blank lines as transcript-group separators that `gtfparse.read_gtf` trips
  on).

## Per-sample drivers

Each per-sample directory (e.g. `DENOVO_ID/CellLines/CL_K562_E0_human/`)
contains a small shell script:

```bash
#!/bin/bash
set -ex
if [ ! -f 'bmark.ok' ]; then
    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType ref_free \
        --truth_gtf ../reference_data/Mouse_isoseqsim_truth.gtf \
        --truth_quant ../reference_data/Mouse_isoseqsim_truth_quant.tsv
    touch bmark.ok
fi
```

The runner (`bmark_nb_runner.py`):
1. Walks up from `cwd` to find the regime's `tool_registry.yaml`.
2. Loads the registry; filters to entries with `display: true` (or to
   the explicit set named in `--report_set <yaml>`, if given).
3. Globs `raw_prog_results/*` and routes each file through `QuantParser`,
   producing `processed_prog_results/<name>.tsv` (and possibly `<name>.gtf`).
4. Builds a `program_files = {name: {"quant", "gtf", "color", "family", "venn"}}`
   dict and injects it into the chosen template via papermill.
5. The template emits per-sample metric TSVs and a per-sample notebook.

`bmark.ok` is the success sentinel: the run script writes it only after the
runner exits cleanly. Missing `bmark.ok` = failed sample.

### CLI flags

- `--analysisType {quant_only|quant_only_no_truthset|ref_free|ref_guided}`
  -- selects the template notebook.
- `--truth_gtf <path>` -- reference annotation GTF. For simulated
  DENOVO_ID/REF_Guided runs this is typically an expressed-only GTF, so
  unexpressed annotation isoforms never enter the identification truth set.
- `--truth_quant <path>` -- ground-truth quant (synthetic data) or proxy
  (Oarfish-byAlignment for human samples without ground truth).
- `--truth_reduced_gtf <path>` -- (ref_guided only) reduced 1-isoform-per-gene
  GTF used as the starting annotation.
- `--Venn_mode` -- enable consensus truth-set augmentation (see "Venn mode"
  below).
- `--registry <path>` -- override the registry-discovery walk-up.
- `--report_set <yaml>` -- override which registry entries are active and
  their `display`/`venn` settings, without editing the registry.

### Truth-set and FN semantics for isoform ID

For the DENOVO_ID and REF_Guided notebooks, the truth table is built by
joining `--truth_gtf` and `--truth_quant` on `transcript_id` and then
collapsing transcripts to `intronId`s. In code, this is an `inner` merge,
so only transcripts present in BOTH files enter `i_ref_df`.

This has two practical consequences:

- In simulated datasets, unexpressed annotation isoforms are normally
  excluded upstream by passing an expressed-only `--truth_gtf`
  (for example `*_annotation_expressed.gtf`) and a truth-quant file that
  contains only expressed transcripts. Those isoforms therefore never become
  false negatives in the isoform-identification metrics.
- During TP/FP/FN assignment, an intronId is only a truth positive if its
  `ref_tpm > 0`. Any truth intronId with zero expression would not count as
  a TP or FN anyway, but most current synthetic datasets exclude those rows
  before scoring.

## Per-dataset orchestration

Each per-dataset directory has a `Makefile` that runs all samples and
aggregates their outputs. After the refactor, every Makefile shares the same
shape:

```make
MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

DIRS := CL_BT474_E0_human CL_HG002_E0_human CL_K562_E0_human CL_UHRR_E0_human

.PHONY: all $(DIRS) report
all: $(DIRS) report

$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c "set -o pipefail && ./run_*.sh | tee script.log 2>&1"
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Finished $@"

report:
	@missing=""; \
	for d in $(DIRS); do [ -f $$d/bmark.ok ] || missing="$$missing $$d"; done; \
	if [ -n "$$missing" ]; then \
		echo "FAILED samples (no bmark.ok):"; \
		for d in $$missing; do echo "  $$d  (see $$d/script.log)"; done; \
		exit 1; \
	else echo "All samples completed successfully."; fi

clean:
	rm -f ./<glob>/*ref_guided* ./<glob>/*.gtf ./<glob>/bmark.ok

summarize:
	find ./ | grep <metric>.tsv > <name>.files.list
	../../benchmarking/misc/aggregate_result_tables.py \
		<name>.files.list > <dataset>.<metric>.tsv
```

Two design points worth knowing:

- `MAKEFLAGS += -j6` caps parallelism. Big-genome notebooks (CellLines,
  SG-NEx) parse GRCh38 and can OOM-kill each other if too many run at once.
  Six is a safe ceiling on a 64 GB host.

- The `-` prefix on the `$(DIRS)` recipe lets one sample fail without
  blocking the rest. Failure is detected via `bmark.ok` absence in the
  `report` target -- not via a `failed_targets.txt` file. The `bmark.ok`
  sentinel is created by the user-script, so missing-sentinel survives
  even an OOM-kill that takes out the make subshell mid-recipe.

## How to add a new tool or version

The most common case. Suppose IsoQuant releases v3.14.0 and you want to
benchmark it.

### Step 1: Deposit the outputs

The IsoQuant docker (in the submodule's `workflows/IsoQuant_docker/`)
already names outputs as
`<sample>.isoquant-v<version>.IsoQuant.{counts.tsv,gtf}`. Run the WDL
across your samples and drop the resulting files into
`<regime>/__comparator_results_all_<REGIME>/IsoQuant/` (or wherever the
regime's deposit script expects to find them).

Then run the deposit script to copy them into per-sample
`raw_prog_results/` dirs:

```
cd QUANT_ONLY_outputs/__comparator_results_all_QUANT_ONLY
python3 deposit_comparator_outputs.py
```

If the deposit script sees a file that no registry entry recognizes, it
prints a loud warning:

```
WARNING: 1 deposited file(s) matched no entry in the tool registry.
They will be ignored by bmark_nb_runner.py until you add a matching
entry to tool_registry.yaml:
  - CL_K562_E0_human.isoquant-v3.14.0.IsoQuant.counts.tsv
```

### Step 2: Add a registry entry

Edit `<regime>/tool_registry.yaml` and add an entry. Place the new
versioned entry above the older one within the same family so the parser
matches the version-bearing pattern first:

```yaml
  - name: IsoQuant-v3.14.0
    family: IsoQuant
    quant_pattern: 'isoquant-v3\.14\.0\.IsoQuant\.counts\.tsv$'
    quant_id_col: 0
    quant_tpm_col: 1
    gtf_pattern: 'isoquant-v3\.14\.0\.IsoQuant\.gtf$'
    gtf_source: own         # for DENOVO_ID/REF_Guided; use 'ref' for QUANT_ONLY
    color: navy
    display: true
    venn: true              # if this is now the newest IsoQuant version

  - name: IsoQuant-v3.13.0  # was venn:true; flip to false now that v3.14.0 supersedes it
    family: IsoQuant
    ...
    venn: false
```

(Update the older sibling's `venn` flag as part of the same edit so the
"newest only seeds the consensus" convention holds.)

### Step 3: Verify the registry catches the new files

```
cd <regime>
python3 - <<'PY'
import sys, os, re, glob
sys.path.insert(0, '../benchmarking/pylib')
import QuantParser
entries = QuantParser.load_registry('tool_registry.yaml')
unmatched = []
for d in glob.glob('*/raw_prog_results') + glob.glob('*/*/raw_prog_results'):
    for fn in os.listdir(d):
        path = os.path.join(d, fn)
        if not os.path.isfile(path): continue
        if any((e['quant_pattern'] and re.search(e['quant_pattern'], fn))
               or (e['gtf_pattern'] and re.search(e['gtf_pattern'], fn))
               for e in entries):
            continue
        unmatched.append(path)
print(f"Unmatched: {len(unmatched)}")
for u in unmatched[:10]: print(f"  {u}")
PY
```

Coverage of `raw_prog_results/` files (excluding scripts and stray backups)
should be 100%.

### Step 4: Re-run

```
cd <regime>
make clean   # wipes processed_prog_results, bmark.ok, generated notebooks
make         # rebuilds; failures surface via the report target
make summarize  # regenerates aggregated TSVs and rerenders summary plots
```

### Step 5: Re-render the Rmd reports (optional)

Per-dataset Rmd reports live alongside the aggregated TSVs:

```
cd <regime>/CellLines
Rscript -e 'rmarkdown::render("eval_cell_line_RefGuided_ID.Rmd",
                              output_format="github_document")'
```

The Rmds consume long-format aggregated TSVs (with a `Name` cell-value
column), so new tool entries appear automatically without Rmd edits. One
gotcha: if you ever switch an Rmd to gather a wide TSV's columns into a
long format, use `read.csv(..., check.names=FALSE)` -- otherwise hyphens
in tool names get munged to dots and the gather no longer joins back.

## How to add a new dataset

Less common. To add a new dataset family (e.g. a new sim regime):

1. Create `<regime>/<dataset>/` with subdirs per sample.
2. In each per-sample dir, write a `run_*.sh` invoking `bmark_nb_runner.py`
   with the right `--analysisType`, `--truth_gtf`, `--truth_quant`, etc.
   For synthetic isoform-ID datasets, prefer an expressed-only `--truth_gtf`
   so the truth set and FN universe are limited to simulated transcripts.
3. Drop a `Makefile` matching the template in "Per-dataset orchestration"
   above. Set `DIRS := <list of sample dirs>`.
4. Set up the `summarize` target to glob the per-sample metric TSVs you
   care about and run them through `aggregate_result_tables.py`.
5. Add the new dataset to the regime's top-level `Makefile`'s `DIRS`
   list so `make` at the regime root processes it.
6. Drop a per-dataset `.Rmd` that consumes the aggregated TSVs and
   produces a summary plot/PDF.
7. If the dataset has no synthetic ground truth (real biological samples),
   you'll likely need a proxy `--truth_quant`. The convention is to use
   Oarfish-v0.9.4 byAlignment quant -- symlink it from
   `QUANT_ONLY_outputs/<sample_dir>/raw_prog_results/`. Pass `--Venn_mode`
   to enable consensus truth-set augmentation.

## How to add a new analysis regime

Rare. The QUANT_ONLY -> DENOVO_ID -> REF_Guided pattern transferred wholesale
in this refactor; a fourth regime would follow the same recipe:

1. Create `<NEW_REGIME>/tool_registry.yaml` with the relevant tool subset.
2. Author a new `template.<analysisType>.ipynb` in the submodule
   (`benchmarking/template_notebooks/`) and wire its
   parameters/colors/`prog_quant_files` cells to derive from the
   `program_files` dict the runner injects.
3. Add the new analysisType to `bmark_nb_runner.py`'s
   `analysisType_to_notebook` mapping.
4. Create per-dataset Makefiles + per-sample run scripts.
5. Drop a regime-top-level Makefile + a deposit script under
   `__comparator_results_all_<NEW_REGIME>/`.

## Venn mode

For datasets without a synthetic ground truth (real biological samples in
CellLines and SG-NEx for DENOVO_ID and REF_Guided), the per-sample run
script passes `--Venn_mode`. The template's logic kicks in only when this
flag is set:

1. Among entries with `venn: true`, count how many programs predict each
   intronId.
2. IntronIds that are NOT in the reference annotation but are predicted by
   at least 2 venn-eligible programs are added to the truth set (i.e. they
   stop counting as false positives for any program).
3. The reference set is restricted to introns predicted by at least one
   venn-eligible program. Therefore, a reference isoform that no
   venn-eligible method called is removed from the truth set before TP/FN
   classification and is not counted as a false negative.
4. ALL programs (regardless of `venn`) are scored against this augmented
   truth set.

This means: only the newest version of each tool family contributes to
defining what "truth" looks like, but every version is still measured
against it.

## Common gotchas

- **OOM on parallel make**: Each per-dataset Makefile is now capped at
  `-j6`. If you launch multiple dataset families in parallel from a shell
  (`(cd a && make) & (cd b && make) &`), the global parallelism is
  6 x N_families. On a 64 GB host that's fine for ~4 families.

- **Stale Docker images**: WDL workflows pin `:latest`. After updating a
  tool's runner script, rebuild + push the docker image (each docker dir
  has `build_docker.sh` and `push_docker.sh`) before invoking the WDL
  again. Symptoms of a stale image: the WDL fallback fires, an unexpected
  output file gets picked, downstream metrics look wrong.

- **LRAA's blank-line GTF**: LRAA's denovo and ref-guided GTFs include
  blank lines as transcript-group separators. `gtfparse.read_gtf` trips
  on those, so LRAA's registry entry must use `gtf_converter: strip_blanks`.

- **read.csv hyphen munging in Rmds**: R's `read.csv` defaults to
  `check.names=TRUE`, which munges hyphens in column headers to dots.
  Long-format TSVs (with `Name` as a cell-value column) are unaffected,
  but wide-format TSVs (one column per tool) need
  `read.csv(..., check.names=FALSE)` plus `rename(methodA = `Unnamed: 0`)`.

- **bmark.ok is the source of truth**: per-sample success is signaled by
  the file `<sample>/bmark.ok`. Don't `touch` this manually; let the run
  script write it after `bmark_nb_runner.py` completes successfully.

## Source

This skeleton was extracted from the LRAA paper analysis repository after the
registry-driven benchmarking refactor. It intentionally excludes input payloads
and generated outputs.
