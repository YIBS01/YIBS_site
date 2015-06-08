################# NEW Makefile Format ###################################

.SUFFIXES:
.PHONY: FORCE

AVAILABLE_COMPONENTS = YIBS shared

export MOD_DIR := $(shell pwd)/mod
export INCLUDE_DIR := $(shell pwd)/include
export GISSCLIM_DIR := $(shell pwd)/..
export CONFIG_DIR := $(GISSCLIM_DIR)/config
export SCRIPTS_DIR := $(GISSCLIM_DIR)/exec
export COMPILER

export RUN_H = $(INCLUDE_DIR)/rundeck_opts.h

# export global options from OPTS_MAIN line here
export FVCUBED
export FVCORE

#MAKE1 = $(MAKE) -f make_subdirs.mk

sinclude $(DECKS_DIR)/$(RUN).mk

all: $(COMPONENTS:=_dir)
	echo $^

define PROGRAM_template
$(1)_dir: $$($(1)_OBJS) $$($(1)_LIBS:%=-l%)
	@echo
	@echo "===> building component $(1)"
	$(MAKE) -C $(1) $(OPTS_$(1))
#	mv -f .liblist .liblist_tmp
	echo $(1)/lib$(1).a >> .liblist
#	cat .liblist_tmp >> .liblist
#	rm -f .liblist_tmp
	@echo "===> component $(1) ok"
	@echo

$(1)_dep: $$($(1)_OBJS) $$($(1)_LIBS:%=-l%) $(RUN_H)
	$(MAKE) -C $(1) depend $(OPTS_$(1))

$(1)_htmldoc:
	$(MAKE) -C $(1) htmldoc $(OPTS_$(1))

endef

$(foreach prog,$(COMPONENTS),$(eval $(call PROGRAM_template,$(prog))))


DEPENDFILE_SUBDIRS = .depend_subdirs
ifneq ($(MAKECMDGOALS),$(DEPENDFILE_SUBDIRS))
sinclude $(DEPENDFILE_SUBDIRS)
endif

depend_all: $(COMPONENTS:=_dep)
	$(SCRIPTS_DIR)/comp_mkdep.pl $(COMPONENTS)
	+$(MAKE) depend


FSRCS_TMP = $(addsuffix .f,$(strip $(OBJ_LIST)))
F90SRCS_TMP = $(addsuffix .F90,$(strip $(OBJ_LIST)))
FORTSRCS_PRESENT = $(wildcard *.f) $(wildcard *.F90)
FSRCS = $(filter $(FORTSRCS_PRESENT), $(FSRCS_TMP))
F90SRCS = $(filter $(FORTSRCS_PRESENT), $(F90SRCS_TMP))


include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk


#include $(GISSCLIM_DIR)/config/rules.mk
#COMPLIBS = $(patsubst %, %/lib.a, $(COMPONENTS))

#$(OBJS): $(COMPONENTS:=_dir)

COMPLIBS = $(shell perl -e 'print reverse <>;' < .liblist)

do_components: $(COMPONENTS:=_dir)

do_main $(RUN).bin:   $(OBJS) #  $(COMPONENTS:=_dir)
	@echo "===> linking"
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) $(OBJS) $(F90OBJS) $(ESMF_OBJS) \
	  $(COMPLIBS) $(LIBS) -o $(RUN).bin  $(LINK_OUTPUT)
	@echo "===> linking ok"
	@echo

main gcm: $(MOD_DIR)
	-rm .liblist
	touch .liblist
	+$(MAKE)  do_components
	+$(MAKE) do_main

echo_vars:
	@echo CPP_OPTIONS = $(CPP_OPTIONS)
	@echo OBJ_LIST = $(OBJ_LIST)
	@echo COMPONENTS = $(COMPONENTS)
	@echo INPUT_FILES = $(INPUT_FILES)
	@echo RUN_PARAMETERS = $(RUN_PARAMETERS)
	@echo INPUTZ = $(INPUTZ)
	+$(MAKE) main

# the following line is for compatibility with old interface only
clean_all vclean: clean

# "clean" will clean components and then will do "clean" in main directory
clean: clean_components

clean_components:
	-rm -f $(RUN_H) $(MOD_DIR)/*.mod
	for i in $(AVAILABLE_COMPONENTS) ; do \
	  $(MAKE) -C $$i clean ; done

$(RUN_H): $(DECKS_DIR)/$(RUN).mk FORCE
	perl -e '$$_="$(CPP_OPTIONS)"; s/ *\#/\n\#/g; print "$$_\n";' \
	 > rundeck_opts.tmp
	if ! cmp -s rundeck_opts.tmp $(RUN_H) ; then \
	  mv rundeck_opts.tmp $(RUN_H) ; \
	else \
	  rm rundeck_opts.tmp ; \
	fi

$(DEPENDFILE): $(RUN_H)

$(LIB):

gcmlib: $(MOD_DIR)
	-rm .liblist
	touch .liblist
	+$(MAKE) do_components
	+$(MAKE) $(LIB)
#	mv -f .liblist .liblist_tmp
	echo $(LIB) >> .liblist
#	cat .liblist_tmp >> .liblist
#	rm -f .liblist_tmp

$(MOD_DIR):
	mkdir $(MOD_DIR)

# html documentation
htmldoc_all: $(COMPONENTS:=_htmldoc)
	+$(MAKE) htmldoc


