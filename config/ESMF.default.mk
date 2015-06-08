# these are the options needed to compile the code with ESMF library
# (MPI support should be included separately)


CPPFLAGS += -DUSE_ESMF

# hack to get compatible compiler name for ESMF directory
ESMF_COMPILER ?= $(COMPILER)
ifeq ($(COMPILER),IRIX64)
ESMF_COMPILER = default
endif

# the following variables specify location of ESMF library and includes
# they can be overwritten in ~/.modelErc file if necessary
ifeq ($(ESMF_DIR),)
ESMFINCLUDEDIR ?= ${BASELIBDIR}/include/esmf
ESMFLIBDIR ?= ${BASELIBDIR}/lib
else
ESMFINCLUDEDIR ?= ${ESMF_DIR}/mod/mod${ESMF_BOPT}/$(MACHINE).$(ESMF_COMPILER).64.default
ESMFLIBDIR ?= ${ESMF_DIR}/lib/lib${ESMF_BOPT}/$(MACHINE).$(ESMF_COMPILER).64.default
endif

# the following tells make where to look for system mod files
VPATH += ${ESMFINCLUDEDIR}

FFLAGS += -I${ESMFINCLUDEDIR}
F90FLAGS += -I${ESMFINCLUDEDIR}
LIBS += -L${ESMFLIBDIR} -lesmf

# if we don't have netcdf library add netcdf_stubs
ifndef NETCDFHOME
LIBS +=  -lnetcdf_stubs
endif
