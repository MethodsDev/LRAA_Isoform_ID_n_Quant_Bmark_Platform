MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

# List of directories
DIRS := morf2_ont_merged_annot_compat_1isoform  morf2_pacbio_merged_annot_compat_1isoform

# Define targets for each directory
.PHONY: all $(DIRS) report

all: $(DIRS) report

# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c 'set -o pipefail; ./run_quant_bmarking.sh 2>&1 | tee -a script.log'
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Finished $@" | tee -a "$@/script.log"

# Report failures at the end
report: $(DIRS)
	@missing=""; \
	for d in $(DIRS); do [ -f $$d/bmark.ok ] || missing="$$missing $$d"; done; \
	if [ -n "$$missing" ]; then \
		echo "FAILED samples (no bmark.ok):"; \
		for d in $$missing; do echo "  $$d  (see $$d/script.log)"; done; \
		exit 1; \
	else echo "All samples completed successfully."; fi

# Clean logs and failure reports
.PHONY: clean
clean:
	rm -f ./morf2_*/*.tsv ./morf2_*/*ipynb ./morf2_*/processed_prog_results/*
	rm -rf ./*/.ipynb_checkpoints
	rm -f ./*/bmark.ok
	rm -f ./*.list ./*.tsv


summarize:
	find morf2_* | grep spearman_expr_cor.tsv$ > spearman_expr_cor.files.list
	../../benchmarking/misc/aggregate_result_tables.py spearman_expr_cor.files.list > morf2.spearman_cor.tsv

	find morf2_* | grep pearson_expr_cor.tsv$ > pearson_expr_cor.files.list
	../../benchmarking/misc/aggregate_result_tables.py pearson_expr_cor.files.list > morf2.pearson_cor.tsv

	find morf2_* | grep median_rel_diff.tsv$ > median_rel_diff.files.list
	../../benchmarking/misc/aggregate_result_tables.py median_rel_diff.files.list > morf2.median_rel_diff.tsv

	find morf2_* | grep mean_rel_diff.tsv$ > mean_rel_diff.files.list
	../../benchmarking/misc/aggregate_result_tables.py mean_rel_diff.files.list > morf2.mean_rel_diff.tsv

	find morf2_* | grep oarfish_style_metrics.tsv$ > oarfish_style_metrics.files.list
	../../benchmarking/misc/aggregate_result_tables.py oarfish_style_metrics.files.list > morf2.oarfish_style_metrics.tsv
