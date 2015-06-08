#makefile for the component YIBS
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF PFT_MODEL PS_MODEL RAD_MODEL MIXED_CANOPY_OPT 

FSRCS = 
F90SRCS = 

FSRCS += yibs_prescribed_drv.f 

FSRCS += yibs_mod.f yibs.f \
cohorts.f patches.f yibs_cells.f \
physutil.f canopyspitters.f \
allometryfn.f phenology.f respauto_physio.f \
soilbgc.f yibs_const.f \
yibs_types.f yibs_prescr_veg.f yibs_prescribed_updates.f \
yibs_prognostic.f yibs_pfts16_FLUXNET.f \
FBBphotosynthesis.f FBBpfts16_FLUXNET.f

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
#dependencies (should be done later automatically):
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

ifeq ($(PFT_MODEL),YIBS)
CPPFLAGS += -DPFT_MODEL_YIBS
endif

ifeq ($(SITE), YES)
CPPFLAGS += -DSOILCARB_SITE
endif

ifeq ($(YIBS_1D_DIAG),YES)
CPPFLAGS += -DYIBS_1D_DIAG
endif

ifeq ($(MIXED_CANOPY_OPT),YES)
  CPPFLAGS += -DMIXED_CANOPY
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

ifeq ($(PS_BVOC),YES)
  CPPFLAGS += -DPS_BVOC
endif

ifeq ($(FAST_SPINUP),YES)
  CPPFLAGS += -DFAST_SPINUP
endif

ifeq ($(RESTART_CPOOLS),YES)
  CPPFLAGS += -DRESTART_CPOOLS
endif

