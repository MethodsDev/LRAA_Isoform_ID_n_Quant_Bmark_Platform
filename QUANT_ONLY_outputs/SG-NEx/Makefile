MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

DIRS=SGNex_A549 SGNex_MCF7 

.PHONY: all $(DIRS)

all: $(DIRS) report
# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c 'set -o pipefail; ./run_quant_only_comparison.sh 2>&1 | tee -a script.log'
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Finished $@" | tee -a "$@/script.log"


# Clean logs and failure reports
.PHONY: clean
clean:
	rm -f ./SG*/*denovo*
	rm -f ./SG*/*.gtf
	rm -f ./*/bmark.ok

summarize:
	find ./ | grep correlation_matrix.csv > cor_files.list
	../../benchmarking/misc/aggregate_result_tables.py cor_files.list > all_expr_correlations.tsv

# Report failures at the end (no bmark.ok = failure)
.PHONY: report
report: $(DIRS)
	@missing=""; \
	for d in $(DIRS); do [ -f $$d/bmark.ok ] || missing="$$missing $$d"; done; \
	if [ -n "$$missing" ]; then \
		echo "FAILED samples (no bmark.ok):"; \
		for d in $$missing; do echo "  $$d  (see $$d/script.log)"; done; \
		exit 1; \
	else echo "All samples completed successfully."; fi
