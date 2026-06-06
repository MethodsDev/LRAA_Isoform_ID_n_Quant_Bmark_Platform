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

summarize:
	find e* | grep spearman_expr_cor.tsv$ > spearman_expr_cor.files.list
	../../benchmarking/misc/aggregate_result_tables.py spearman_expr_cor.files.list > SIRVs.spearman_cor.tsv

	find e* | grep pearson_expr_cor.tsv$ > pearson_expr_cor.files.list
	../../benchmarking/misc/aggregate_result_tables.py pearson_expr_cor.files.list > SIRVs.pearson_cor.tsv

	find e* | grep median_rel_diff.tsv$ > median_rel_diff.files.list
	../../benchmarking/misc/aggregate_result_tables.py median_rel_diff.files.list > SIRVs.median_rel_diff.tsv

	find e* | grep mean_rel_diff.tsv$ > mean_rel_diff.files.list
	../../benchmarking/misc/aggregate_result_tables.py mean_rel_diff.files.list > SIRVs.mean_rel_diff.tsv

	find e* | grep oarfish_style_metrics.tsv$ > oarfish_style_metrics.files.list
	../../benchmarking/misc/aggregate_result_tables.py oarfish_style_metrics.files.list > SIRVs.oarfish_style_metrics.tsv

