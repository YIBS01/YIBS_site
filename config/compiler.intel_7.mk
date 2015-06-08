# just in case, options for older intel copieler (7 and earlier)
# probably not working anyway...

F90 = efc
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -h
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler INTEL-ifc-on-LINUX
FFLAGS = -fpp -Wp,-P -O2 -w95 -w90 -cm -tpp2 -common_args
F90FLAGS = -fpp -Wp,-P -O2 -FR -w95 -w90 -cm -tpp2 -common_args
LFLAGS = -O2 -w95 -w90 -tpp2 -common_args -Vaxlib
CPPFLAGS = -DCOMPILER_Intel
F90_VERSION = $(shell $(F90) -V 2>&1 | grep Build)
ifeq ($(MP),YES)
FFLAGS += -openmp
F90FLAGS += -openmp
LFLAGS += -openmp
endif
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -WB 
LFLAGS += -WB 
F90FLAGS += -WB
LFLAGSF += -WB
endif
