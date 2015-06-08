#
# this file contains rules shared by Makefiles in "model" and "aux"
#

.PHONY:
ifdef MOD_DIR
VPATH = $(MOD_DIR)
endif

#VPATH = $(MODULES_CMPLR_OPT)
#VPATH = /usr/local/unsupported/esmf/1.0.3-551/include
######  Some user customizable settings:   ########

# EXTRA_FFLAGS specifies some extra flags you want to pass 
# to Fortarn compiler, like
# -g        - include debugging information
# -listing  - create listings (.L)
EXTRA_FFLAGS = 

# EXTRA_LFLAGS specifies some extra flags you want to pass 
# to linker. Currently needed as a hack to compile hybrid MPI/OpenMP
# code to pass "-openmp" to linker
EXTRA_LFLAGS =

# hack to force compilation errors to be written to ERR files rather than scren
ifeq ($(OUTPUT_TO_FILES),YES)
  COMP_OUTPUT = > $*.ERR 2>&1 || { r=$$? ; cat $*.ERR ; exit $$r ; }
  LINK_OUTPUT = > $(RUN).ERR 2>&1 || { r=$$? ; cat $(RUN).ERR ; exit $$r ; }
else
  COMP_OUTPUT =
  LINK_OUTPUT =
endif

# if -s specified enable some extra messages
ifeq ($(findstring s,$(MFLAGS)),s)
MSG = 
else
MSG = > /dev/null
endif

#
# starting machine - specific options
#

NO_COMMAND = echo "*****  Requested target is not supported on $(UNAME)"; \
             echo "*****  or compiler is not specified properly."; \
             echo "*****  You have COMPILER=$(COMPILER)" ; exit 1;
F90 = $(NO_COMMAND)
FMAKEDEP = $(NO_COMMAND)
CMP_MOD = cmp -s
SETUP = $(SCRIPTS_DIR)/setup_e.pl
CPP = $(NO_COMMAND)
LIBS =
INCS = 
F90_VERSION = 'Unknown compiler version'
ECHO_FLAGS =
CPPFLAGS =
# the following is default fortran flag for path to include and .mod files
# it is redefined later for compilers with non-standard flags
I = I
# by default assume that fortran compiler can do cpp
EXTERNAL_CPP = NO
# by default assume that we build 64-bit code with gcc-like compiler
CFLAGS = -O2 -m64
# define default name for m4
M4 = m4
# default runlib
RANLIB = ranlib

# RFLAGS returns rundeck options for the current object (i.e. for $*)
# it has effect only of OBJ_LIST_O is defined
RFLAGS = $(shell perl \
-e '$$_="$(OBJ_LIST_O)"; m/\b$* *\|([^|]*)\|/; print " $$1";' \
)


# check consistency of compilation flugs and redefine some
# flags if necessary
ifeq ($(ESMF),YES)
MPI = YES
endif

ifeq ($(FVCORE),YES)
MPI = YES
endif

ifeq ($(FVCUBED),YES)
MPI = YES
endif

# hack to keep Intel8 name valid (only temporarily)
ifeq ($(COMPILER),Intel8)
COMPILER=intel
endif

# include machine-specific options
MACHINE = $(shell uname)
include $(CONFIG_DIR)/machine.$(MACHINE).mk

# if COMPILER is not defined try default
ifndef COMPILER
COMPILER = $(MACHINE)
endif

#include compiler-specific options
include $(CONFIG_DIR)/compiler.$(COMPILER).mk

# hack to include netcdf_stubs only when netcdf is missing
#ifdef NETCDFHOME
#NETCDF_STUBS=
#else
#NETCDF_STUBS=-lnetcdf_stubs
#endif


# include ESMF library if necessary
ifeq ($(ESMF),YES)
  include $(CONFIG_DIR)/ESMF.default.mk
endif

ifeq ($(FVCUBED),YES)
  ifndef FVCUBED_ROOT
     FVCUBED_ROOT = false
  endif
  
  # Cubed-sphere requires FVCORE and MPP enabled
  FVCORE=YES
  MPP=YES

  # Testing options: ADIABATIC
  ADIABATIC=YES

  CUBED_SPHERE=YES

  FVINC = -I$(FVCUBED_ROOT)/$(MACHINE)/include
  FVINCS = $(FVINC) $(FVINC)/MAPL_Base $(FVINC)/MAPL_cfio $(FVINC)/FVdycoreCubed_GridComp
  INCS += $(FVINCS)
  FVINCx = $(FVCUBED_ROOT)/$(MACHINE)/include
  FVINCSx = $(FVINCx):$(FVINCx)/MAPL_Base:$(FVINCx)/FVdycoreCubed_GridComp
  ifdef SYSTEM_MOD_DIRS
    SYSTEM_MOD_DIRS = $(SYSTEM_MOD_DIRS):$(FVINCSx)
  else
    SYSTEM_MOD_DIRS = $(FVINCSx)
  endif
  LIBS += -L$(FVCUBED_ROOT)/$(MACHINE)/lib -lMAPL_cfio -lMAPL_Base -lFVdycoreCubed_GridComp -lfvdycore
  # this extra -lesmf would not be needed if the ESMF stuff came after this section
  LIBS += -lesmf
endif

ifeq ($(FVCORE),YES)
  ifndef FVCORE_ROOT
     FVCORE_ROOT = false
  endif
  CPPFLAGS += -DUSE_FVCORE 
  ifneq ($(FVCUBED),YES)
    FVINC = -I$(FVCORE_ROOT)/$(MACHINE)/include
    CPPFLAGS += -DFVCUBED_SKIPPED_THIS -DCREATE_FV_RESTART 
    INCS += $(FVINC) $(FVINC)/GEOS_Base $(FVINC)/GEOS_Shared $(FVINC)/GMAO_gfio_r8 $(FVINC)/GMAO_cfio_r8 $(FVINC)/GMAO_pilgrim $(FVINC)/FVdycore_GridComp  -I$(BASELIBDIR)/include
    LIBS += -L$(FVCORE_ROOT)/$(MACHINE)/lib  -lFVdycore_GridComp  -lGMAO_pilgrim -lGMAO_gfio_r8 -lGMAO_cfio_r8 -lGEOS_Shared -lGEOS_Base -L$(BASELIBDIR)/lib
    LIBS += -L${BASELIBDIR}/lib -lesmf
  endif
endif

ifeq ($(SKIP_FV),YES)
  CPPFLAGS+=-DSKIP_FV
endif

ifeq ($(CUBED_SPHERE),YES)
  CPPFLAGS += -DCUBED_SPHERE
  INCS += -I$(FFTW_ROOT)/include
  LIBS += -L$(FFTW_ROOT)/lib -lfftw3
endif

ifeq ($(MPP),YES)  
CPPFLAGS += -DUSE_MPP
FFLAGS += -I$(MPPDIR)/include
F90FLAGS += -I$(MPPDIR)/include
LIBS += -L$(MPPDIR)/lib -lfms_mpp_shared 
LIBS += -lfmpi -lmpi
endif

ifeq ($(USE_ENT),YES)
  CPPFLAGS += -DUSE_ENT
  FFLAGS += -$(I)Ent
  F90FLAGS += -$(I)Ent
endif

ifeq ($(ADIABATIC),YES)
  CPPFLAGS += -DADIABATIC
endif

ifdef PNETCDFHOME
  LIBS += -L$(PNETCDFHOME)/lib -lpnetcdf
  FFLAGS += -I$(PNETCDFHOME)/include
  INCS += -I$(PNETCDFHOME)/include
endif


ifeq ($(MPI),YES)
  ifeq ($(MPIDISTR),)
    # unknown distribution, just trying to add a mpi librarry ...
    LIBS += -lmpi
  else
    include $(CONFIG_DIR)/mpi.$(MPIDISTR).mk
  endif
endif

#
# Check for extra options specified in modelErc
#

ifdef NETCDFHOME
ifeq ($(MACHINE),IRIX64)
  LIBS += -L$(NETCDFHOME)/lib64 -lnetcdf
else
  LIBS += -L$(NETCDFHOME)/lib -lnetcdf
endif
  FFLAGS += -I$(NETCDFHOME)/include
  INCS += -I$(NETCDFHOME)/include
endif

#ifdef PNETCDFHOME
#  LIBS += -L$(PNETCDFHOME)/lib -lpnetcdf
#  FFLAGS += -I$(PNETCDFHOME)/include
#  INCS += -I$(PNETCDFHOME)/include
#endif

# access new interfaces in sub-directory.
ifdef ESMF_Interface
FFLAGS += -$(I)$(ESMF_Interface)
F90FLAGS += -$(I)$(ESMF_Interface)
endif
CPPFLAGS += $(INCS)

# path to the modules dir if present
ifdef MOD_DIR
  FFLAGS += -$(I)$(MOD_DIR)
  F90FLAGS += -$(I)$(MOD_DIR)
endif

ifdef INCLUDE_DIR
  CPPFLAGS += -I$(INCLUDE_DIR)
endif

# add path to MPI includes
ifdef MPIDIR
  CPPFLAGS += -I$(MPIDIR)/include
endif



ifeq ($(COMPARE_MODULES_HACK),NO)
CMP_MOD = cmp -s
endif

# add fartran flags into a single string
ifeq ($(EXTERNAL_CPP),YES)
FFLAGS_ALL =  $(FFLAGS) $(EXTRA_FFLAGS)
F90FLAGS_ALL = $(F90FLAGS) $(EXTRA_FFLAGS)
else
FFLAGS_ALL =  $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS)
F90FLAGS_ALL = $(F90FLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS)
endif

ifdef SYSTEM_MOD_DIRS
VPATH += $(SYSTEM_MOD_DIRS)
endif

#
# Pattern  rules
#

*.mod: .current_options

%.mod:
	@echo checking $@
	@if [ "$<empty" = "empty" ]; then \
	echo "No dependency for $@ : assuming it is a system module";\
	else \
	cmp $< $@ ; \
	if ! cmp -s $< $@ ; then \
	  echo "will copy $< to $@" ; \
	  cp $< $@ ; \
	fi ; \
	fi

%.smod:
	@if [ ! -f $@ ] ; then rm -f $< ; $(MAKE) $< RUN=$(RUN); fi


# Standard fortran
ifeq ($(EXTERNAL_CPP),YES)
%.o: %.f.cpp.f
else
%.o: %.f
endif
	@echo $(ECHO_FLAGS)  compiling $< ... $(MSG) \\c
	$(F90) -c -o $@ $(FFLAGS_ALL) $(RFLAGS) $< $(COMP_OUTPUT)
	@if [ -s $(DEPENDFILE) ] ; then \
	for i in \
	`perl -e 'while(<>){ if(/(\S+)\.mod: *(\w+\@$*\.smod)/){print " $$1";} }' $(DEPENDFILE)` ; \
	do \
	    cp -f $$i.mod $$i\@$*\.smod ; cp $$i\@$*\.smod $$i.mod ; \
	done ; \
	fi
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif

ifeq ($(EXTERNAL_CPP),YES)
%.o: %.F90.cpp.F90
else
%.o: %.F90
endif
	@echo $(ECHO_FLAGS)  compiling $< ... $(MSG) \\c
	$(F90) -c -o $@ $(F90FLAGS_ALL) $(RFLAGS) $< $(COMP_OUTPUT)
	@if [ -s $(DEPENDFILE) ] ; then \
	for i in \
	`perl -e 'while(<>){ if(/(\S+)\.mod: *(\w+\@$*\.smod)/){print " $$1";} }' $(DEPENDFILE)` ; \
	do \
	    cp -f $$i.mod $$i\@$*\.smod ; touch $$i.mod ; \
	done ; \
	fi
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif

# cpp preprocessing

%.f.cpp.f: %.f
	@echo $(ECHO_FLAGS)  preprocessing $< ... $(MSG) \\c
	$(CPP) $(CPPFLAGS) $*.f | sed -n '/^#pragma/!p' > $@

%.F90.cpp.F90: %.F90
	@echo $(ECHO_FLAGS)  preprocessing $< ... $(MSG) \\c
	$(CPP) $(CPPFLAGS) $*.F90 | sed -n '/^#pragma/!p' > $@

%.f.cpp: %.f
	@echo preprocessing $<  $(MSG)
	$(CPP) $(CPPFLAGS) $*.f > $*.f.cpp

%.F90.cpp: %.F90
	 @echo preprocessing $<  $(MSG)
	 $(CPP) $(CPPFLAGS) $*.F90 > $*.F90.cpp

%.o: %.c
	cc -c -O2 -m64 $*.c

ifneq ($(MACHINE),IRIX64)

%.f: %.m4f
	m4 $*.m4f > $*.f

%.F90: %.m4F90
	m4 $*.m4F90 > $*.F90
endif



# end of Pattern  rules

