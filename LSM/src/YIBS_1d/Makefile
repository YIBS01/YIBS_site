# makefile for component YIBS_1d
.SUFFIXES:

SUPPORTED_OPTIONS = PFT_MODEL MIXED_CANOPY_OPT SOILCARB_SITE_OPT

FSRCS = yibs_prog.f  lenc.f \
        handle_err.f read_site_forcing.f calc_solarzen.f \
        read_veg_forcing.f write_site_flux.f 

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/../config
endif

include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

CPPFLAGS += -DSOILCARB_SITE  -DYIBS_1D_DIAG -DFLUXNETINIT

ifeq ($(PFT_MODEL),YIBS)
CPPFLAGS += -DPFT_MODEL_YIBS
endif

ifeq ($(MIXED_CANOPY_OPT),YES)
  CPPFLAGS += -DMIXED_CANOPY
endif

ifeq ($(SOILCARB_SITE_OPT),YES)
  CPPFLAGS += -DSOILCARB_SITE
endif

ifeq ($(O3DEP_UPTAKE),YES)
  CPPFLAGS += -DO3DEP_UPTAKE
endif

ifeq ($(O3DEP_UPTAKE_OFFLINE),YES)
  CPPFLAGS += -DO3DEP_UPTAKE_OFFLINE
endif

ifeq ($(ACTIVE_GROWTH),YES)
  CPPFLAGS += -DACTIVE_GROWTH
endif

ifeq ($(OUTPUT_FORCING),YES)
  CPPFLAGS += -DOUTPUT_FORCING
endif

ifeq ($(HEIGHT_SPINUP),YES)
  CPPFLAGS += -DHEIGHT_SPINUP
endif



