
ifneq (${MPIDIR},)
FFLAGS += -I${MPIDIR}/include
F90FLAGS += -I${MPIDIR}/include
LIBS += -L${MPIDIR}/lib
endif

# try to work around memory leak
CPPFLAGS += -DMPITYPE_LOOKUP_HACK

#LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
#-lfmpi -lmpi -lstdc++ -threads

##LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
##-lmpi_f90 -lmpi_f77 -lmpi -lstdc++ -threads

LIBS += -lmpi_f77 -lmpi -lmpi_cxx -lstdc++
ifneq ($(shell uname),Darwin)
  LIBS += -lrt
endif

#LIBS +=  -lmpi -lmpigc3 -lmpigc4 -lmpigf -lmpigi -lmpiic4 -lmpiic -lmpiif -threads

