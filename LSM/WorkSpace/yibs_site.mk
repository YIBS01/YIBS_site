# rundeck for compiling stand-alone version of GISS LSM
# coupled to YIBS with FBB photosynthesis

OBJ_LIST=main_yibs

# components to use:
COMPONENTS = YIBS_1d model/YIBS model/shared

# specific options for components:
OPTS_YIBS = PFT_MODEL=YIBS ACTIVE_GROWTH=YES RESTART_CPOOLS=YES \
            O3DEP_UPTAKE_OFFLINE=YES PS_BVOC=YES
OPTS_YIBS_1d = PFT_MODEL=YIBS O3DEP_UPTAKE_OFFLINE=YES ACTIVE_GROWTH=YES
#OPTS_giss_LSM = USE_YIBS=YES
OPTS_shared =
#OPTS_giss_LSM_standalone = USE_YIBS=YES

DATADIR=/home/YIBS_site/Input/DRIVER_DATA

INPUT_FILES = \
soil_textures=$(DATADIR)/soil_textures_top30cm_2x2.5 \
SOILCARB_global=$(DATADIR)/soilcarb_top30cm_nmaps_2x2.5bin.dat \
SOILCARB_restart=$(DATADIR)/soil_carbon_wfdei_1980_1x1 \
SOILCARB_site=$(DATADIR)/SOILCARB_site_MMSF.csv \
VEG16=$(DATADIR)/V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ij \
VHT16=$(DATADIR)/V144x90_EntMM16_height_trimmed_scaled_nocrops.ij \
VCROPS=$(DATADIR)/CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp \
VLAIX16=$(DATADIR)/V144x90_EntMM16_lai_max_trimmed_scaled_nocrops.ij \
LAIM16_01=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Jan.ij \
LAIM16_02=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Feb.ij \
LAIM16_03=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Mar.ij \
LAIM16_04=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Apr.ij \
LAIM16_05=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_May.ij \
LAIM16_06=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Jun.ij \
LAIM16_07=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Jul.ij \
LAIM16_08=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Aug.ij \
LAIM16_09=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Sep.ij \
LAIM16_10=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Oct.ij \
LAIM16_11=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Nov.ij \
LAIM16_12=$(DATADIR)/V144x90_EntMM16_lai_trimmed_scaled_nocrops_Dec.ij \
CROPS_CAL=$(DATADIR)/crop_dates.nc


RUN_PARAMETERS = \
&input_parameters \
year=1999,jday=1,year2=2004,jday2=1,dt=3600.0, \
lon=-68.747, lat=45.209, site='US-HO2', id_veg=4 \
/
