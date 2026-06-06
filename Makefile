REGIMES := QUANT_ONLY DENOVO_ID REF_Guided

.PHONY: all $(REGIMES) clean summarize

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

