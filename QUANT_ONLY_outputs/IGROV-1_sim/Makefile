MAKEFLAGS += -j6   # cap parallelism to avoid OOM on big-genome notebooks


# Define targets for each directory
.PHONY: all report

all:
	@echo "[$(shell date '+%Y-%m-%d %H:%M:%S')] Starting " | tee "script.log"
	-bash -c "set -o pipefail && ./run_quant_bmarking.sh | tee script.log 2>&1"
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
	rm -f ./*.tsv ./*ipynb ./processed_prog_results/*
	rm -rf ./.ipynb_checkpoints


summarize:
	echo "nothing to do here"
