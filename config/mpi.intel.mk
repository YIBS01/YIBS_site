
#LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
#-lfmpi -lmpi -lstdc++ -threads

ifeq ($(IFORT_RELEASE),9.1)
LIBS += -lcprts -limf -lm -lcxa -lunwind -lrt -ldl \
-lmpiif -lmpi -lstdc++ -threads
else
LIBS += -limf -lm -lrt -ldl \
-lmpiif -lmpi_mt -lstdc++ -threads
endif

#-lmpigf -lmpi -lstdc++ -threads


#LIBS +=  -lmpi -lmpigc3 -lmpigc4 -lmpigf -lmpigi -lmpiic4 -lmpiic -lmpiif -threads

