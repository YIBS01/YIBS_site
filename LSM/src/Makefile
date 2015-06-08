# makefile for src root directory
.SUFFIXES:

# the following line should list all the components which are
# present in src directory
AVAILABLE_COMPONENTS = YIBS_1d model/YIBS  model/shared

.PHONY: FORCE
#ifndef GISSCLIM_DIR
#       GISSCLIM_DIR := $(shell pwd)/..
#       export GISSCLIM_DIR
#endif
#ifndef CONFIG_DIR
#       CONFIG_DIR := $(GISSCLIM_DIR)/config
#endif

export MOD_DIR := $(shell pwd)/../mod
export INCLUDE_DIR := $(shell pwd)/../include
export GISSCLIM_DIR := $(shell pwd)/../..
export CONFIG_DIR := $(GISSCLIM_DIR)/config
export SCRIPTS_DIR := $(GISSCLIM_DIR)/LSM/scripts
export COMPILER

HOMEDIR = $(wildcard ~)
MODELERC ?= $(HOMEDIR)/.modelErc
sinclude $(MODELERC)


#hack to work with current Rules.make
ESMF_Interface = .
#export MOD_DIR = $(GISSCLIM_DIR)/mod
#export INCLUDE_DIR = $(GISSCLIM_DIR)/include

sinclude $(RUNDECK)

OBJS = $(OBJ_LIST:=.o)

all: $(COMPONENTS:=_dir)
	echo $^

define PROGRAM_template
$(1)_dir: $$($(1)_OBJS) $$($(1)_LIBS:%=-l%)
	@echo
	@echo "===> building component $(1)"
	$(MAKE) -C $(1) $(OPTS_$(shell basename $(1)))
	mv -f .liblist .liblist_tmp
	echo $(1)/lib$(shell basename $(1)).a > .liblist
	cat .liblist_tmp >> .liblist
	rm -f .liblist_tmp
	@echo "===> component $(1) ok"
	@echo

$(1)_dep: $$($(1)_OBJS) $$($(1)_LIBS:%=-l%)
	$(MAKE) -C $(1) depend $(OPTS_$(shell basename $(1)))

endef

$(foreach prog,$(COMPONENTS),$(eval $(call PROGRAM_template,$(prog))))

clean:
	-rm *.o *.mod *.smod .depend ../mod/*.mod .liblist

clean_all: clean
	for i in $(AVAILABLE_COMPONENTS) ; do \
	  $(MAKE) -C $$i clean ; done

#YIBS_dir: shared_dir
#giss_LSM_dir: YIBS_dir shared_dir
#giss_LSM_standalone_dir: giss_LSM_dir YIBS_dir shared_dir

DEPENDFILE = .depend
ifneq ($(MAKECMDGOALS),$(DEPENDFILE))
sinclude $(DEPENDFILE)
endif

depend: $(COMPONENTS:=_dep)
	../scripts/comp_mkdep.pl $(COMPONENTS)

something:
	echo ok


#include $(GISSCLIM_DIR)/arch/base.mk
include $(GISSCLIM_DIR)/config/rules.mk

#COMPLIBS = $(patsubst %, %/lib.a, $(COMPONENTS))

COMPLIBS = $(shell cat .liblist)

do_main $(RUN).bin:  $(COMPONENTS:=_dir)  $(OBJS)
	@echo "===> linking"
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) $(OBJS) $(F90OBJS) $(ESMF_OBJS) \
	  $(COMPLIBS) $(LIBS) -o $(RUN).bin  $(LINK_OUTPUT)
	@echo "===> linking ok"
	@echo
main:
	-rm .liblist
	touch .liblist
	$(MAKE) do_main
