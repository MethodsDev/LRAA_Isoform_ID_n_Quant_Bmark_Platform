# ISM S/L pairs and internal TSSs (LRAA v0.51.0), PBMC, ALS, bulk PacBio. Reads existing GTFs and BAMs only; no LRAA run.
# Whole pipeline: about 1 h. Heavy steps run in a memory-capped scope.
SHELL := /bin/bash
RUN := systemd-run --user --scope -q -p MemoryMax=20G -p MemorySwapMax=0
DS := PBMC ALS BT474 K562 HG002 UHRR
.PHONY: all tables report
all: tables report
tables:
	cd bin && $(RUN) python3 gencode_cache.py ../out/gencode47.pkl
	cd bin && for d in $(DS); do $(RUN) python3 enumerate_pairs.py $$d ../out/pairs/$$d.pairs.tsv.gz; done
	cd bin && for d in $(DS); do $(RUN) python3 scan_reads.py $$d 6; done
	cd bin && for d in $(DS); do $(RUN) python3 score_sites.py $$d ../out/scored/$$d.sites.tsv.gz && $(RUN) python3 protein.py $$d ../out/scored/$$d.pairs.tsv.gz; done
	cd bin && $(RUN) python3 assemble.py ../out/stage1.tsv $(DS) && $(RUN) python3 run_orfs.py ../out/stage1.tsv
	cd bin && for d in $(DS); do $(RUN) python3 protein.py $$d ../out/scored/$$d.pairs.tsv.gz ../out/orf/orfs.tsv; done
	cd bin && $(RUN) python3 assemble.py ../out/candidates.full.tsv $(DS) && gzip -f ../out/candidates.full.tsv
	cd bin && $(RUN) python3 summarize.py ../out/candidates.full.tsv.gz && $(RUN) python3 select_ftd_example.py ../out/candidates.full.tsv.gz && $(RUN) python3 make_examples.py ../out/candidates.full.tsv.gz && python3 uniprot_features.py
	rm -f out/stage1.tsv
report:
	rm -rf ISM_internal_TSS_files
	$(RUN) Rscript -e 'rmarkdown::render("ISM_internal_TSS.Rmd")' > logs/knit.log 2>&1
# tool comparison (bulk, tool-neutral pools; tool_own_ism.py = Supp Fig panel c): about 2 h
tools:
	cd bin && ./run_pool.sh && $(RUN) python3 tool_recall.py && $(RUN) python3 tool_own_ism.py
