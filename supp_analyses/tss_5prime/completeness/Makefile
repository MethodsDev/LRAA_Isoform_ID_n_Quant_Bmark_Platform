# 5' completeness and the untemplated-G (cap) signature, PBMC (whole cells) and ALS (cortex
# nuclei): measurement only, for reviewing whether LRAA's TSS degradation pruning should
# take read 5' completeness into account. See README.md.
#
#   make          both datasets' tables (one at a time, $(CPU) workers each)
#   make report   also knit tss_5prime_completeness.Rmd

SHELL     := /bin/bash
LRAA_HOME ?= $(HOME)/LongReadAlignmentAssembler
CPU       ?= 12
SCRIPT    := $(LRAA_HOME)/util/sc/tss_5prime_completeness.py
CTAT      := $(HOME)/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir
GTF       := $(CTAT)/ref_annot.gtf
GENOME    := $(CTAT)/ref_genome.fa
PBMC      := ../PBMCs_kinnex/LRAA_PBMCs_eval
ALS       := ../ALS_LP00714/ALS_LRAA_eval
RDNA_MASK := ../__shared_refs/rdna_mask.GRCh38.bed
CAGE      := $(PBMC)/__Site_Usage_Analysis/ext/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
COMMON    := --gtf $(GTF) --genome $(GENOME) --HiFi --rdna_mask_bed $(RDNA_MASK) --cage_peaks $(CAGE) --CPU $(CPU)

.PHONY: all report
all: PBMC.genes.tsv ALS.genes.tsv

PBMC.genes.tsv: $(SCRIPT)
	python3 $(SCRIPT) $(COMMON) --bam ../PBMCs_kinnex/PBMCs_pbio.aligned.sorted.bam \
		--sites $(PBMC)/__Site_Usage_Analysis/PBMCs.sites.tsv --events $(PBMC)/__Site_Usage_Analysis/PBMCs.dexseq.TSS.events.tsv \
		--cell_clusters $(PBMC)/pbmcs_refGuided/data/RefQuantOnly.cell_clusters.tsv --output_prefix PBMC > PBMC.log 2>&1

# after PBMC: one dataset at a time
ALS.genes.tsv: $(SCRIPT) | PBMC.genes.tsv
	python3 $(SCRIPT) $(COMMON) --bam ../ALS_LP00714/LP00714.aligned.sorted.bam \
		--sites $(ALS)/__Site_Usage_Analysis/LP00714.sites.tsv --events $(ALS)/__Site_Usage_Analysis/LP00714.dexseq.TSS.events.tsv \
		--cell_clusters $(ALS)/als_refGuided/data/RefQuantOnly.LP00714.ALS.cell_clusters.tsv --output_prefix ALS > ALS.log 2>&1

report: all
	Rscript -e 'rmarkdown::render("tss_5prime_completeness.Rmd")' > knit.log 2>&1

# ---------------------------------------------------------------------------------------
# Tandem TSS pairs (two TSSs of one splice pattern, the downstream one inside the upstream
# one's first exon): read, geometry, sequence and external features that might separate a
# real downstream TSS from a truncation product, without the cap signature. Labelled sets:
# CAGE-supported tandem pairs vs read-start peaks in the first exons of single-TSS genes.
# MCF7 (SG-NEx ONT cDNA, LRAA v0.40.0 de novo) is the cross-platform check.
#
#   make tandem          the three feature tables (one dataset at a time)
#   make tandem_report   also knit tandem_tss_features.Rmd

TSCRIPT   := $(LRAA_HOME)/util/sc/tandem_tss_features.py
REFTSS    := ext/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz
TCOMMON   := --genome $(GENOME) --cage_peaks $(CAGE) --reftss $(REFTSS) --ref_gtf $(GTF) --CPU $(CPU)
MCF7_RUN  := ../__LRAA_local_runs/bmark_v0.40.0-5a45947
PBMC_CG   := $(PBMC)/pbmcs_refGuided/pbmcs_refGuided_sc_cluster_guided/data
ALS_CG    := $(ALS)/als_refGuided/als_refGuided_sc_cluster_guided/data

.PHONY: tandem tandem_report
tandem: tandem/PBMC.pairs.tsv tandem/ALS.pairs.tsv tandem/MCF7.pairs.tsv

$(REFTSS):
	mkdir -p ext && curl -sS -L -o $@ https://reftss.riken.jp/datafiles/current/human/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz

tandem/PBMC.pairs.tsv: $(TSCRIPT) PBMC.genes.tsv $(REFTSS)
	mkdir -p tandem
	python3 $(TSCRIPT) $(TCOMMON) --bam ../PBMCs_kinnex/PBMCs_pbio.aligned.sorted.bam --HiFi --rdna_mask_bed $(RDNA_MASK) \
		--lraa_gtf $(PBMC_CG)/PBMCs_pbio_CG_refguided.withGeneSymbols.gtf \
		--cell_clusters $(PBMC)/pbmcs_refGuided/data/RefQuantOnly.cell_clusters.tsv \
		--isoT_events $(PBMC)/__Splice_Pattern_Usage_Analysis/PBMCs.dexseq.IsoformTermini.events.tsv \
		--site_events $(PBMC)/__Site_Usage_Analysis/PBMCs.dexseq.TSS.events.tsv --sites $(PBMC)/__Site_Usage_Analysis/PBMCs.sites.tsv \
		--single_tss_genes PBMC.genes.tsv --negative_sites $(PBMC)/__Site_Usage_Analysis/PBMCs.sites.tsv \
		--output_prefix tandem/PBMC > tandem/PBMC.log 2>&1

tandem/ALS.pairs.tsv: $(TSCRIPT) ALS.genes.tsv $(REFTSS) | tandem/PBMC.pairs.tsv
	python3 $(TSCRIPT) $(TCOMMON) --bam ../ALS_LP00714/LP00714.aligned.sorted.bam --HiFi --rdna_mask_bed $(RDNA_MASK) \
		--lraa_gtf $(ALS_CG)/LP00714.refguided.withGeneSymbols.gtf \
		--cell_clusters $(ALS)/als_refGuided/data/RefQuantOnly.LP00714.ALS.cell_clusters.tsv \
		--isoT_events $(ALS)/__Splice_Pattern_Usage_Analysis/LP00714.dexseq.IsoformTermini.events.tsv \
		--site_events $(ALS)/__Site_Usage_Analysis/LP00714.dexseq.TSS.events.tsv --sites $(ALS)/__Site_Usage_Analysis/LP00714.sites.tsv \
		--single_tss_genes ALS.genes.tsv --negative_sites $(ALS)/__Site_Usage_Analysis/LP00714.sites.tsv \
		--output_prefix tandem/ALS > tandem/ALS.log 2>&1

# MCF7: one sample, no cell clusters, ONT read filters. Positives also include PBMC's
# CAGE-supported tandem pairs measured in the MCF7 reads (MCF7's own de novo models give
# only a few dozen). Negatives: PBMC's and ALS's single-TSS genes whose FANTOM5 CAGE status
# is "single" (one promoter across FANTOM5), at their LRAA TSS site.
# PBMC's and ALS's single-TSS genes (one row per main transcript, PBMC first) with each
# gene's TSS position from its own dataset's site table
tandem/single_tss_genes.PBMC_ALS.tsv: PBMC.genes.tsv ALS.genes.tsv pool_single_tss_genes.py
	mkdir -p tandem
	python3 pool_single_tss_genes.py PBMC:PBMC.genes.tsv:$(PBMC)/__Site_Usage_Analysis/PBMCs.sites.tsv \
		ALS:ALS.genes.tsv:$(ALS)/__Site_Usage_Analysis/LP00714.sites.tsv > $@

tandem/MCF7.pairs.tsv: $(TSCRIPT) tandem/single_tss_genes.PBMC_ALS.tsv $(REFTSS) tandem/PBMC.pairs.tsv | tandem/ALS.pairs.tsv
	python3 $(TSCRIPT) $(TCOMMON) --bam $(MCF7_RUN)/inputs/SGNex_MCF7_cDNAStranded_replicate3_run2.aligned.sorted.bam \
		--lraa_gtf $(MCF7_RUN)/stage/denovo/SGNex_MCF7.LRAA.ref-free.gtf \
		--ref_pairs tandem/PBMC.pairs.tsv \
		--single_tss_genes tandem/single_tss_genes.PBMC_ALS.tsv --negative_cage single \
		--output_prefix tandem/MCF7 > tandem/MCF7.log 2>&1

tandem_report: tandem
	Rscript -e 'rmarkdown::render("tandem_tss_features.Rmd")' > knit.tandem.log 2>&1
