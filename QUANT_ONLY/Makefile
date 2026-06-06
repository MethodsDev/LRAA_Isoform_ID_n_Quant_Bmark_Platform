# List of directories
DIRS := arabidopsis_sim  CellLines  MORFs SIRVs SG-NEx  mouse_sim 

# Define targets for each directory
.PHONY: all $(DIRS)

all: $(DIRS)

# Each directory runs its script with logging and error handling
$(DIRS):
	@cd $@ && make -j


# Clean logs and failure reports
.PHONY: clean
clean:
	@for i in $(DIRS); do \
    echo "Running clean in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done


.PHONY: summarize
summarize:
	@for i in $(DIRS); do \
    echo "Running summarize in $$i..."; \
    (cd $$i; $(MAKE) summarize) || exit $$?; done


