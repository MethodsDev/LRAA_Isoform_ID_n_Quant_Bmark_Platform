MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

# List of directories
DIRS := mouse_isoseqsim_e000  mouse_isoseqsim_e016  mouse_isoseqsim_e050  mouse_isoseqsim_e085 

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
	rm -f ./mouse_isoseqsim_*/*.tsv ./mouse_isoseqsim_*/*ipynb ./mouse_isoseqsim_*/processed_prog_results/*
	rm -rf ./*/.ipynb_checkpoints
	rm -f ./*.list ./*.tsv
	rm -f ./*/bmark.ok

summarize:
	find mouse_isoseqsim_* | grep spearman_expr_cor.tsv$ > spearman_expr_cor.files.list
	../../benchmarking/misc/aggregate_result_tables.py spearman_expr_cor.files.list > mouse_isoseqsim.spearman_cor.tsv

	find mouse_isoseqsim_* | grep pearson_expr_cor.tsv$ > pearson_expr_cor.files.list
	../../benchmarking/misc/aggregate_result_tables.py pearson_expr_cor.files.list > mouse_isoseqsim.pearson_cor.tsv

	find mouse_isoseqsim_* | grep median_rel_diff.tsv$ > median_rel_diff.files.list
	../../benchmarking/misc/aggregate_result_tables.py median_rel_diff.files.list > mouse_isoseqsim.median_rel_diff.tsv

	find mouse_isoseqsim_* | grep mean_rel_diff.tsv$ > mean_rel_diff.files.list
	../../benchmarking/misc/aggregate_result_tables.py mean_rel_diff.files.list > mouse_isoseqsim.mean_rel_diff.tsv

	find mouse_isoseqsim_* | grep oarfish_style_metrics.tsv$ > oarfish_style_metrics.files.list
	../../benchmarking/misc/aggregate_result_tables.py oarfish_style_metrics.files.list > mouse_isoseqsim.oarfish_style_metrics.tsv

