
F90 = f95
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_NAG
FFLAGS = -fpp -O0 -maxcontin=100 -kind=byte -dusty
F90FLAGS = -fpp -O2 -free -kind=byte
LFLAGS =
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS +=
LFLAGS += -lefence
endif
