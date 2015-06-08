.PHONY: FORCE
ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/..
endif
HOMEDIR = $(wildcard ~)
MODELERC ?= $(HOMEDIR)/.modelErc
sinclude $(MODELERC)

export RUN
export RUNDECK = $(shell pwd)/$(RUN).mk
#RUNDIR=$(shell pwd)/$(RUN)
RUNDIR=$(SAVEDISK)/$(RUN)

# location of the main modelE code
MODEL_E_ROOT = $(shell pwd)/../src/
MODEL_DIR = $(MODEL_E_ROOT)/model
export MODEL_DIR

#hack to deal with incompatible "echo" (should be moved to Rules.make)
# ... and mpirun
ECHO = echo -e
MPIRUN_FLAGS =
ifeq ($(MPIDISTR),SCALI)
MPIRUN_FLAGS = -inherit_limits
endif
ifeq ($(shell uname),IRIX64)
ECHO = echo
endif
ifeq ($(shell uname),Darwin)
ECHO = echo
endif

sinclude $(RUNDECK)

all: $(RUNDECK)
	$(MAKE) -C $(GISSCLIM_DIR)/src depend
	$(MAKE) -C $(GISSCLIM_DIR)/src main

depend: $(RUNDECK)
	$(MAKE) -C $(GISSCLIM_DIR)/src depend

clean:
	$(MAKE) -C $(GISSCLIM_DIR)/src clean_all


setup:  all
	@echo
	@echo "===> creating run directory"
	-mkdir $(RUNDIR)
	if [ ! -d $(RUNDIR) ] ; then \
	  echo "can't create directory $(RUNDIR), aborting" ; exit 1 ; fi
	if [ -f $(RUNDIR)/lock ] ; then \
	  echo "lock file is present in $(RUNDIR), aborting"; exit 1 ; fi
	cp $(GISSCLIM_DIR)/src/$(RUN).bin $(RUNDIR)
	@echo "===> creating execution script"
	echo "#/bin/bash"                             > $(RUNDIR)/$(RUN)
	echo ""                                    >> $(RUNDIR)/$(RUN)
	echo "ulimit -s 512000"                    >> $(RUNDIR)/$(RUN)
	echo cd $(RUNDIR)                          >> $(RUNDIR)/$(RUN)
	perl -e 'while($$_=shift){ s/(\S*)=(\S*)/ln -s $$2 $$1/; \
	  print "$$_\n";}' $(INPUT_FILES)          >> $(RUNDIR)/$(RUN)
	echo touch lock                            >> $(RUNDIR)/$(RUN)
#ifeq ($(ESMF),YES)
#	echo "time mpirun -np 4 $(MPIRUN_FLAGS) ./$(RUN).bin > $(RUN).log 2>&1" >> $(RUNDIR)/$(RUN)
#else
	echo "time ./$(RUN).bin > $(RUN).log"      >> $(RUNDIR)/$(RUN)
#endif
	echo rm lock                               >> $(RUNDIR)/$(RUN)
	perl -e 'while($$_=shift){ s/(\S*)=(\S*)/rm $$1/; \
	  print "$$_\n";}' $(INPUT_FILES)          >> $(RUNDIR)/$(RUN)
	chmod 755 $(RUNDIR)/$(RUN)
	@echo "===> creating yibs_input file"
	$(ECHO) "$(RUN_PARAMETERS)" > $(RUNDIR)/yibs_input
	@echo "===> all done "

param:
	$(ECHO) "$(RUN_PARAMETERS)"

all_nodep: $(RUNDECK)
	$(MAKE) -C $(GISSCLIM_DIR)/src main

