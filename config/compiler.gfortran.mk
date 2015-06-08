
F90 = gfortran
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_G95
FFLAGS = -cpp -fconvert=big-endian -O2 -fno-range-check
F90FLAGS = -cpp -fconvert=big-endian -O2 -fno-range-check -ffree-line-length-none
LFLAGS =

#
# Set the following to ensure that the beginning/end of records
# in sequential-access unformatted files have 4-byte markers
# (this does impose a 2 GB limit on record sizes, however)
#
#FFLAGS += -frecord-marker=4
#F90FLAGS += -frecord-marker=4

# machine-specific options
ifeq ($(MACHINE),IRIX64)
FFLAGS += -mabi=64 
F90FLAGS += -mabi=64
LFLAGS += -mabi=64
endif

# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow,denormal -fbacktrace
F90FLAGS += -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow,denormal -fbacktrace
#LFLAGS += -lefence
endif
