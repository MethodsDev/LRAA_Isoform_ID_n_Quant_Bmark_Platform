# Library-level calibration of the TSS cap signature: per-peak capped fractions across
# libraries, the decision LRAA's pre-pass makes from them, and the FTD false TSS sites.
#
#   make            every table (one library at a time; each scan uses CPU workers)
#   make report     also knit cap_calibration.Rmd
#
# The scans run LRAA's own module (pylib/TSS_cap_calibration.py). LRAA_HOME defaults to
# the shared checkout; until the pre-pass is merged there, point it at the branch:
#   make LRAA_HOME=/path/to/worktree
# Each scan reads only 5' ends (no assembly): ~2.5 min for FTD (48M reads) on 12 workers.

SHELL     := /bin/bash
LRAA_HOME ?= $(HOME)/LongReadAlignmentAssembler
CPU       ?= 12
ROOT      := ../..
RUNS      := $(ROOT)/__LRAA_local_runs
CTAT      := $(HOME)/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa
GRCH38    := $(ROOT)/DENOVO_ID/CellLines/reference_data/GRCh38_no_alt.fa
SIRVFA    := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/refs/SIRV_isoforms_multi-fasta_170612a.fasta
MORFFA    := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/refs/minigenome.fa
ARAFA     := $(RUNS)/bmark_v0.35.0-3c8cfe3_quantonly/refs/Athaliana_447_TAIR10.fa
CAGE      := $(ROOT)/PBMCs_kinnex/LRAA_PBMCs_eval/__Site_Usage_Analysis/ext/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
CAP       := systemd-run --user --scope -q -p MemoryMax=30G -p MemorySwapMax=0

# library  bam  genome  cage(-)
LIB_FTD      := $(ROOT)/ALS_LP00714/LP00714.aligned.sorted.bam $(CTAT) $(CAGE)
LIB_PBMC     := $(ROOT)/PBMCs_kinnex/PBMCs_pbio.aligned.sorted.bam $(CTAT) $(CAGE)
LIB_BT474    := $(RUNS)/bmark_v0.40.0-5a45947/inputs/BT474_E0_merged_other_sorted.bam $(GRCH38) $(CAGE)
LIB_K562     := $(RUNS)/bmark_v0.35.0_big/inputs/K562_E0_merged_other_sorted.bam $(GRCH38) $(CAGE)
LIB_MCF7     := $(RUNS)/bmark_v0.40.0-5a45947/inputs/SGNex_MCF7_cDNAStranded_replicate3_run2.aligned.sorted.bam $(GRCH38) $(CAGE)
LIB_A549     := $(RUNS)/bmark_v0.40.0-5a45947_a549/inputs/SGNex_A549_cDNAStranded_replicate5_run2.aligned.sorted.bam $(GRCH38) $(CAGE)
LIB_MORF_PB  := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/inputs/morf2_pacbio_merged_annot_compat_sorted.bam $(MORFFA) -
LIB_MORF_ONT := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/inputs/morf2_ont_merged_annot_compat_sorted.bam $(MORFFA) -
LIB_HG002    := $(RUNS)/bmark_v0.35.0_big/inputs/HG002_E0_merged_other_sorted.bam $(GRCH38) $(CAGE)
LIB_UHRR     := $(RUNS)/bmark_v0.35.0_big/inputs/UHRR_E0_merged_other_sorted.bam $(GRCH38) $(CAGE)
LIB_SIRV     := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/inputs/BT474_E0_merged_sirv_sorted.bam $(SIRVFA) -
LIB_SIRV_HG002 := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/inputs/HG002_E0_merged_sirv_sorted.bam $(SIRVFA) -
LIB_SIRV_UHRR  := $(RUNS)/bmark_v0.43.0-a424fb6_sirvmorf/inputs/UHRR_E0_merged_sirv_sorted.bam $(SIRVFA) -
LIB_ARA      := $(RUNS)/bmark_v0.35.0-3c8cfe3_quantonly/inputs/arabidopsis_isoseqsim_e016.aligned.sorted.bam $(ARAFA) -

LIBS   := FTD PBMC BT474 K562 HG002 UHRR MCF7 A549 MORF_PB MORF_ONT SIRV SIRV_HG002 SIRV_UHRR ARA
PEAKS  := $(foreach l,$(LIBS),out/$(l).peaks.tsv.gz)
# what a "non-peak" background holds: one chromosome per library (baseline G run as calibrated)
STRATA := out/FTD.chr19.offpeak_strata.tsv out/PBMC.chr19.offpeak_strata.tsv out/BT474.chr19.offpeak_strata.tsv

# per cell cluster: reads at uncapped- vs capped-mode peaks (fixed RefQuantOnly clusters)
CL_FTD   := $(ROOT)/ALS_LP00714/ALS_LRAA_eval/als_refGuided/data/RefQuantOnly.LP00714.ALS.cell_clusters.tsv
CL_PBMC  := $(ROOT)/PBMCs_kinnex/LRAA_PBMCs_eval/pbmcs_refGuided/data/RefQuantOnly.cell_clusters.tsv
BYCLUSTER := out/FTD.peak_mode_by_cluster.tsv out/PBMC.peak_mode_by_cluster.tsv \
	out/FTD.depth_matched_peaks.tsv out/PBMC.depth_matched_peaks.tsv \
	out/FTD.capped_share_by_celltype.tsv out/PBMC.capped_share_by_celltype.tsv

all: $(PEAKS) out/FTD.event_site_caprates.tsv out/PBMC.event_site_caprates.tsv $(STRATA) $(BYCLUSTER)

out/%.peak_mode_by_cluster.tsv: scripts/peak_mode_by_cluster.py out/%.peaks.tsv.gz
	set -- $(LIB_$*); $(CAP) python3 scripts/peak_mode_by_cluster.py $(LRAA_HOME) $$1 $$2 \
		out/$*.peaks.tsv.gz out/$*.TSS_cap_calibration.json $(CL_$*) $@ --cpu $(CPU) > logs/$*.by_cluster.log 2>&1

.NOTPARALLEL:

out/%.peaks.tsv.gz: scripts/export_peaks.py
	@mkdir -p out logs
	set -- $(LIB_$*); cage=$$3; \
	$(CAP) python3 scripts/export_peaks.py $(LRAA_HOME) $$1 $$2 out/$* --cpu $(CPU) \
		$$( [ "$$cage" != "-" ] && echo --cage $$cage ) > logs/$*.log 2>&1

# the current version's high-confidence TSS switch pairs: Supplementary Figure TSS-5' panel e (since 2026-11-11)
EV_FTD := ../../../../LRAA_paper-FTD_study/__Site_Usage_Analysis/LP00714.TSS_events_vs_CAGE.events.tsv
EV_PBMC := ../../../../LRAA_paper-PBMCs_study/__Site_Usage_Analysis/PBMCs.TSS_events_vs_CAGE.events.tsv
ext/%_current_TSS_event_pairs.tsv: scripts/current_event_pairs.py
	python3 scripts/current_event_pairs.py $(EV_$*) $@
out/%.current_event_site_caprates.tsv: scripts/event_site_caprates.py ext/%_current_TSS_event_pairs.tsv
	set -- $(LIB_$*); $(CAP) python3 scripts/event_site_caprates.py $(LRAA_HOME) ext/$*_current_TSS_event_pairs.tsv $$1 $$2 $@ > logs/$*.current_event_sites.log 2>&1

# the v0.45.0 (pre-cap) high-confidence TSS switch pairs (earlier panel e; kept for the record)
out/%.event_site_caprates.tsv: scripts/event_site_caprates.py ext/%_v0.45.0_TSS_event_pairs.tsv
	set -- $(LIB_$*); $(CAP) python3 scripts/event_site_caprates.py $(LRAA_HOME) ext/$*_v0.45.0_TSS_event_pairs.tsv $$1 $$2 $@ > logs/$*.event_sites.log 2>&1

out/FTD.chr19.offpeak_strata.tsv: scripts/offpeak_strata.py
	set -- $(LIB_FTD); python3 scripts/offpeak_strata.py $(LRAA_HOME) $$1 $$2 chr19 0 $@
out/PBMC.chr19.offpeak_strata.tsv: scripts/offpeak_strata.py
	set -- $(LIB_PBMC); python3 scripts/offpeak_strata.py $(LRAA_HOME) $$1 $$2 chr19 0 $@
out/BT474.chr19.offpeak_strata.tsv: scripts/offpeak_strata.py
	set -- $(LIB_BT474); python3 scripts/offpeak_strata.py $(LRAA_HOME) $$1 $$2 chr19 1 $@

report: all
	rm -rf cap_calibration_files
	Rscript -e 'rmarkdown::render("cap_calibration.Rmd")' > knit.log 2>&1

.PHONY: all report

# peaks called within each cell type, every cell type thinned to the same read count
out/%.depth_matched_peaks.tsv: scripts/depth_matched_peaks.py out/%.TSS_cap_calibration.json
	set -- $(LIB_$*); $(CAP) python3 scripts/depth_matched_peaks.py $(LRAA_HOME) $$1 $$2 \
		out/$*.TSS_cap_calibration.json $(CL_$*) $* $@ --cpu $(CPU) > logs/$*.depth_matched.log 2>&1

# per read, peak-free: share of each cell type's 5' starts that are capped
out/%.capped_share_by_celltype.tsv: scripts/capped_share_by_celltype.py out/%.TSS_cap_calibration.json
	set -- $(LIB_$*); $(CAP) python3 scripts/capped_share_by_celltype.py $(LRAA_HOME) $$1 $$2 \
		out/$*.TSS_cap_calibration.json $(CL_$*) $* $@ --cpu $(CPU) > logs/$*.capped_share.log 2>&1
