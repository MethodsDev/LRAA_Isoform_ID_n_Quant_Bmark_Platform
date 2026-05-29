MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks

# List of directories
DIRS := CL_BT474_E0_human CL_HG002_E0_human CL_K562_E0_human CL_UHRR_E0_human

# Define targets for each directory
.PHONY: all $(DIRS) report

all: $(DIRS) report

# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$$(date '+%Y-%m-%d %H:%M:%S')] Starting $@" | tee "$@/script.log"
	-@cd $@ && bash -c 'set -o pipefail; ./run_ref_annot_free_ID_bmarking.sh 2>&1 | tee -a script.log'
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
	rm -f ./CL*/*denovo*
	rm -f ./CL*/*.gtf
	rm -f ./CL*/bmark.ok
	rm -f ./denovo_ID_summary.files.list ./all.Denovo_ID_summary.tsv

summarize:
	find ./ | grep full_transcriptome_TPR_FDR_F1.tsv > denovo_ID_summary.files.list
	../../benchmarking/misc/aggregate_result_tables.py denovo_ID_summary.files.list > all.Denovo_ID_summary.tsv

