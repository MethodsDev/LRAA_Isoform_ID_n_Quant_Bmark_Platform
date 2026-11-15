# Proteoform analysis (the paper's two proteoform supplementary figures), for both studies at once.
#
#   make orfs       ORF prediction on the PBMC and FTD ref-guided sc cluster-guided GTFs (run_orfs.sh,
#                   settings in orf_settings.conf). Optional: the Zenodo cross bundle has the paper's ORFs
#                   in orfs/.
#   make notebook   execute the notebook (into *.executed.ipynb): per-protein tables,
#                   SuppFig_proteoforms_{PBMC,FTD}.pdf, proteoform_counts.tsv. Needs REF_GTF (the GENCODE v47 GTF).
#   make all        orfs, then notebook
#
# The notebook reads the ORFs in orfs/ and, from each study submodule, the SQLIKE category table and the
# sc cluster-guided GTF (their Zenodo bundles unpacked).
NB := LRAA_proteoform_analysis_PBMC_ALS_20260528.ipynb

.PHONY: all orfs notebook

all: orfs notebook

orfs:
	./run_orfs.sh PBMC
	./run_orfs.sh FTD

notebook:
	@test -n "$(REF_GTF)" || { echo "set REF_GTF to the GENCODE v47 GTF" >&2; exit 1; }
	rm -rf PBMC_proteoform_analysis_results ALS_proteoform_analysis_results
	REF_GTF=$(REF_GTF) jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(NB:.ipynb=.executed.ipynb) $(NB)
