MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks


.PHONY: all report

# Each directory runs its script with logging and error handling
all:
	@echo "[$(shell date '+%Y-%m-%d %H:%M:%S')] Starting" | tee "script.log"
	-bash -c "set -o pipefail && ./run_ref_reduced_bmarking.sh | tee script.log 2>&1"
	@echo "[$(shell date '+%Y-%m-%d %H:%M:%S')] Finished" | tee "script.log"

# Report failures at the end
report:
	@if [ -f bmark.ok ]; then \
		echo "Run completed successfully."; \
	else \
		echo "FAILED: no bmark.ok (see script.log)"; \
		exit 1; \
	fi

# Clean logs and failure reports
.PHONY: clean
clean:
	rm -f ./*.tsv ./processed_prog_results/*
	rm -f ./*ipynb
	rm -rf ./.ipynb_checkpoints

.PHONY: summarize
summarize:
	echo "nothing to do here"




