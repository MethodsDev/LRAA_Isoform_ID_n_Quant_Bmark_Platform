MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

DIRS=SGNex_A549 SGNex_MCF7 


all: $(DIRS) report
.PHONY: all $(DIRS)

# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c 'set -o pipefail; ./run_ref_annot_free_ID_bmarking.sh 2>&1 | tee -a script.log'
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Finished $@" | tee -a "$@/script.log"


# Clean logs and failure reports
.PHONY: clean
clean:
	rm -f ./SG*/*denovo*
	rm -f ./SG*/*.gtf
	rm -f ./SG*/bmark.ok
	rm -f ./denovo_ID_summary.files.list ./all.Denovo_ID_summary.tsv

summarize:
	find ./ | grep full_transcriptome_TPR_FDR_F1.tsv > denovo_ID_summary.files.list
	../../benchmarking/misc/aggregate_result_tables.py denovo_ID_summary.files.list > all.Denovo_ID_summary.tsv

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
