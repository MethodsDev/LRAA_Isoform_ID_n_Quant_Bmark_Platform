# List of directories
DIRS := arabidopsis_sim  MORFs SIRVs  CellLines  mouse_sim SG-NEx #IGROV-1_sim 

# Define targets for each directory
.PHONY: all $(DIRS)

all: $(DIRS) 

# Each directory runs its script with logging and error handling
$(DIRS):
	@echo "[$(shell date '+%Y-%m-%d %H:%M:%S')] Starting $@"
	@cd $@ && make -j
	@echo "[$(shell date '+%Y-%m-%d %H:%M:%S')] Finished $@"

# Clean logs and failure reports
.PHONY: clean

clean:
	@for i in $(DIRS); do \
    echo "Running clean in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done

summarize:
	@for i in $(DIRS); do \
    echo "Running summarize in $$i..."; \
    (cd $$i; $(MAKE) summarize) || exit $$?; done




