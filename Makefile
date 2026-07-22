REGIMES := QUANT_ONLY DENOVO_ID REF_Guided
RSCRIPT ?= Rscript

REPORTS := \
	QUANT_ONLY/arabidopsis_sim/Arabidopsis_sim.QUANT_ONLY.summary.Rmd \
	QUANT_ONLY/CellLines/examine_expr_correlations.Rmd \
	QUANT_ONLY/MORFs/Morf2.QUANT_ONLY.summary.Rmd \
	QUANT_ONLY/SG-NEx/examine_expr_correlations.SGNEx.Rmd \
	QUANT_ONLY/SIRVs/SIRVs.QUANT_ONLY.summary.Rmd \
	QUANT_ONLY/mouse_sim/Mouse_sim.QUANT_ONLY.summary.Rmd \
	DENOVO_ID/arabidopsis_sim/Arabidopsis_denovo_eval.Rmd \
	DENOVO_ID/CellLines/eval_cell_line_denovo_ID.Rmd \
	DENOVO_ID/MORFs/Morf2_denovo_eval.Rmd \
	DENOVO_ID/SG-NEx/eval_cell_line_denovo_ID.Rmd \
	DENOVO_ID/SIRVs/SIRVs_denovo_eval.Rmd \
	DENOVO_ID/mouse_sim/Mouse_denovo_eval.Rmd \
	REF_Guided/arabidopsis_sim/Arabidopsis_ref_guided_eval.Rmd \
	REF_Guided/CellLines/eval_cell_line_RefGuided_ID.Rmd \
	REF_Guided/MORFs/Morfs2_ref_guided_eval.Rmd \
	REF_Guided/SG-NEx/eval_cell_line_RefGuided_ID.Rmd \
	REF_Guided/SIRVs/SIRVs_ref_guided_eval.Rmd \
	REF_Guided/mouse_sim/Mouse_ref_guided_eval.Rmd

.PHONY: all $(REGIMES) clean summarize render-reports install-inputs install-references install-predictions

all: $(REGIMES)

$(REGIMES):
	$(MAKE) -C $@

clean:
	@for regime in $(REGIMES); do \
		echo "Running clean in $$regime..."; \
		$(MAKE) -C $$regime clean || exit $$?; \
	done

summarize:
	@for regime in $(REGIMES); do \
		echo "Running summarize in $$regime..."; \
		$(MAKE) -C $$regime summarize || exit $$?; \
	done
	$(MAKE) render-reports

render-reports:
	@for report in $(REPORTS); do \
		echo "Rendering $$report..."; \
		report_dir=$${report%/*}; \
		report_file=$${report##*/}; \
		(cd "$$report_dir" && $(RSCRIPT) -e "rmarkdown::render('$$report_file')") || exit $$?; \
	done

install-inputs: install-references install-predictions

install-references:
	python3 scripts/install_references.py install

install-predictions:
	python3 scripts/install_predictions.py
