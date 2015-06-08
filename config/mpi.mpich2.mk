
ifneq (${MPIDIR},)
LIBS += -L${MPIDIR}/lib
endif

LIBS += -lm -lrt -ldl -lmpich -lfmpich -lstdc++

