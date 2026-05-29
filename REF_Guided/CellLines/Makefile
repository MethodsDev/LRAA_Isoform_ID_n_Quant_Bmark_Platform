MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

# List of directories
DIRS := CL_BT474_E0_human CL_HG002_E0_human CL_K562_E0_human CL_UHRR_E0_human

# Define targets for each directory
.PHONY: all $(DIRS) report

all: $(DIRS) report

# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c 'set -o pipefail; ./run_ref_guided_bmarking.sh 2>&1 | tee -a script.log'
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
	rm -f ./CL_*/*.tsv ./CL_*/*ipynb ./CL_*/processed_prog_results/*
	rm -rf ./CL_*/.ipynb_checkpoints
	rm -f ./CL_*/*.gtf
	rm -f ./CL_*/bmark.ok




summarize:
	find ./ | grep known_and_novel_TPR_FDR_F1.tsv  > refguided_accuracy_summary.files.list
	../../benchmarking/misc/aggregate_result_tables.py refguided_accuracy_summary.files.list > refguided_accuracy_summary.tsv
