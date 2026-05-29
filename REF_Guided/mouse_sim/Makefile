MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

# List of directories
DIRS := mouse_isoseqsim_e000  mouse_isoseqsim_e016  mouse_isoseqsim_e050  mouse_isoseqsim_e085

# Define targets for each directory
.PHONY: all $(DIRS) report

all: $(DIRS) report

# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c 'set -o pipefail; ./run_ref_reduced_bmarking.sh 2>&1 | tee -a script.log'
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

.PHONY: summarize
summarize:
	find mouse_isoseqsim_* | grep known_and_novel_TPR_FDR_F1.tsv$ > known_and_novel_TPR_FDR_F1.files.list
	../../benchmarking/misc/aggregate_result_tables.py known_and_novel_TPR_FDR_F1.files.list > mouse_isoseqsim.known_and_novel_TPR_FDR_F1.tsv


