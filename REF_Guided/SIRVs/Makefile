DIRS=e0_sirv  e1_sirv  e2_sirv

run:
	@for i in $(DIRS); do \
    echo "Running $$i..."; \
    (cd $$i; $(MAKE) -j) || exit $$?; done


clean:
	@for i in $(DIRS); do \
    echo "Running clean in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
	rm -f ./*.list ./*.tsv

.PHONY: summarize
summarize:
	find e* | grep known_and_novel_TPR_FDR_F1.tsv$ > known_and_novel_TPR_FDR_F1.files.list
	../../benchmarking/misc/aggregate_result_tables.py known_and_novel_TPR_FDR_F1.files.list > SIRVs.known_and_novel_TPR_FDR_F1.tsv

