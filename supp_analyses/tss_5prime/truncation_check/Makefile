# 5' truncation check for the ALS TSS switches: where reads of housekeeping genes and of the
# top ALS TSS-switch genes start and end along the gene's main annotated transcript, per
# cluster, with PBMC as a contrast (whole cells, no neuron-favours-downstream pattern).
# Also each cluster's aligned read lengths over a genome-wide sample.
SHELL     := /bin/bash
LRAA_HOME ?= $(HOME)/LongReadAlignmentAssembler
CTAT      := $(HOME)/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir
GTF       := $(CTAT)/ref_annot.gtf
GENOME    := $(CTAT)/ref_genome.fa
RDNA_MASK := ../../../LRAA_paper-PBMCs_study/shared_refs/rdna_mask.GRCh38.bed
GENES     := GAPDH,ACTB,RPL28,EEF1D,GPX4,OAZ1
SCRIPT    := $(LRAA_HOME)/util/sc/read_end_profile_by_cluster.py

ALS_BAM      := ../../LP00714.aligned.sorted.bam
ALS_CLUSTERS := ../als_refGuided/data/RefQuantOnly.LP00714.ALS.cell_clusters.tsv
PBMC_BAM      := ../../../PBMCs_kinnex/PBMCs_pbio.aligned.sorted.bam
PBMC_CLUSTERS := ../../../LRAA_paper-PBMCs_study/pbmcs_refGuided/data/RefQuantOnly.cell_clusters.tsv

.PHONY: all report
all: ALS.reads.tsv.gz PBMC.reads.tsv.gz
report: all
	Rscript -e 'rmarkdown::render("tss_truncation_check.Rmd")' > knit.log 2>&1

ALS.reads.tsv.gz: $(SCRIPT)
	python3 $(SCRIPT) --bam $(ALS_BAM) --gtf $(GTF) --genes $(GENES) --cell_clusters $(ALS_CLUSTERS) \
		--genome $(GENOME) --HiFi --rdna_mask_bed $(RDNA_MASK) --output_prefix ALS > ALS.log 2>&1
PBMC.reads.tsv.gz: $(SCRIPT)
	python3 $(SCRIPT) --bam $(PBMC_BAM) --gtf $(GTF) --genes $(GENES) --cell_clusters $(PBMC_CLUSTERS) \
		--genome $(GENOME) --HiFi --rdna_mask_bed $(RDNA_MASK) --output_prefix PBMC > PBMC.log 2>&1
