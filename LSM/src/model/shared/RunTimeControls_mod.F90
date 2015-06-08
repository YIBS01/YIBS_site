! This file uses the unix m4 macro expansion utility to manipulate
! tokens in a manner beyond the capabilities of CPP.   CPP is of course
! also used in a second pass of preprocessing.











module RunTimeControls_mod
!@sum This module contains a set of Fortran logicals that correspondend
!@+ to CPP conditionals in a 1-1 fashion.   Fortran logicals are preferable,
!@+ since they do not require recompilation, and this module facilitates
!@+ gradual replacement of CPP conditionals with run-time conditionals.
!@+ Initially the Fortran logicals will be parameters, but when a given CPP
!@+ conditional is ready to be completely replaced, true dynamism can be 
!@+ introduced.
!@auth T. Clune
  implicit none
  public






#if defined(ACCMIP_LIKE_DIAGS)
  logical, parameter :: accmip_like_diags = .true.
#else
  logical, parameter :: accmip_like_diags = .false.
#endif

#if defined(ADIABATIC)
  logical, parameter :: adiabatic = .true.
#else
  logical, parameter :: adiabatic = .false.
#endif

#if defined(AG2OG_OCEANS_BUNDLE)
  logical, parameter :: ag2og_oceans_bundle = .true.
#else
  logical, parameter :: ag2og_oceans_bundle = .false.
#endif

#if defined(AG2OG_PRECIP_BUNDLE)
  logical, parameter :: ag2og_precip_bundle = .true.
#else
  logical, parameter :: ag2og_precip_bundle = .false.
#endif

#if defined(ALTER_RADF_BY_LAT)
  logical, parameter :: alter_radf_by_lat = .true.
#else
  logical, parameter :: alter_radf_by_lat = .false.
#endif

#if defined(ALT_CASE5_WINDS)
  logical, parameter :: alt_case5_winds = .true.
#else
  logical, parameter :: alt_case5_winds = .false.
#endif

#if defined(ALT_CLDMIX_UV)
  logical, parameter :: alt_cldmix_uv = .true.
#else
  logical, parameter :: alt_cldmix_uv = .false.
#endif

#if defined(ALT_INTERP)
  logical, parameter :: alt_interp = .true.
#else
  logical, parameter :: alt_interp = .false.
#endif

#if defined(ASCII_REGIONS)
  logical, parameter :: ascii_regions = .true.
#else
  logical, parameter :: ascii_regions = .false.
#endif

#if defined(ATM2x2h_HYCOM1deg)
  logical, parameter :: atm2x2h_hycom1deg = .true.
#else
  logical, parameter :: atm2x2h_hycom1deg = .false.
#endif

#if defined(ATM2x2h_HYCOM2deg)
  logical, parameter :: atm2x2h_hycom2deg = .true.
#else
  logical, parameter :: atm2x2h_hycom2deg = .false.
#endif

#if defined(ATM4x5_HYCOM2deg)
  logical, parameter :: atm4x5_hycom2deg = .true.
#else
  logical, parameter :: atm4x5_hycom2deg = .false.
#endif

#if defined(BC_ALB)
  logical, parameter :: bc_alb = .true.
#else
  logical, parameter :: bc_alb = .false.
#endif

#if defined(BIN_OLSON)
  logical, parameter :: bin_olson = .true.
#else
  logical, parameter :: bin_olson = .false.
#endif

#if defined(BIOGENIC_EMISSIONS)
  logical, parameter :: biogenic_emissions = .true.
#else
  logical, parameter :: biogenic_emissions = .false.
#endif

#if defined(BLK_2M)
  logical, parameter :: blk_2m = .true.
#else
  logical, parameter :: blk_2m = .false.
#endif

#if defined(BLK_2MOM)
  logical, parameter :: blk_2mom = .true.
#else
  logical, parameter :: blk_2mom = .false.
#endif

#if defined(BOT_MONO)
  logical, parameter :: bot_mono = .true.
#else
  logical, parameter :: bot_mono = .false.
#endif

#if defined(BUNDLE_INTERP)
  logical, parameter :: bundle_interp = .true.
#else
  logical, parameter :: bundle_interp = .false.
#endif

#if defined(CALCULATE_FLAMMABILITY)
  logical, parameter :: calculate_flammability = .true.
#else
  logical, parameter :: calculate_flammability = .false.
#endif

#if defined(CALCULATE_LIGHTNING)
  logical, parameter :: calculate_lightning = .true.
#else
  logical, parameter :: calculate_lightning = .false.
#endif

#if defined(CALC_GWDRAG)
  logical, parameter :: calc_gwdrag = .true.
#else
  logical, parameter :: calc_gwdrag = .false.
#endif

#if defined(CHECK_GRID)
  logical, parameter :: check_grid = .true.
#else
  logical, parameter :: check_grid = .false.
#endif

#if defined(CHECK_VAN2)
  logical, parameter :: check_van2 = .true.
#else
  logical, parameter :: check_van2 = .false.
#endif

#if defined(CHL0)
  logical, parameter :: chl0 = .true.
#else
  logical, parameter :: chl0 = .false.
#endif

#if defined(CHL_from_OBIO)
  logical, parameter :: chl_from_obio = .true.
#else
  logical, parameter :: chl_from_obio = .false.
#endif

#if defined(CHL_from_SeaWIFs)
  logical, parameter :: chl_from_seawifs = .true.
#else
  logical, parameter :: chl_from_seawifs = .false.
#endif

#if defined(CLD_AER_CDNC)
  logical, parameter :: cld_aer_cdnc = .true.
#else
  logical, parameter :: cld_aer_cdnc = .false.
#endif

#if defined(COLUMN_TRACER)
  logical, parameter :: column_tracer = .true.
#else
  logical, parameter :: column_tracer = .false.
#endif

#if defined(COMPILER_ABSOFT)
  logical, parameter :: compiler_absoft = .true.
#else
  logical, parameter :: compiler_absoft = .false.
#endif

#if defined(COMPILER_G95)
  logical, parameter :: compiler_g95 = .true.
#else
  logical, parameter :: compiler_g95 = .false.
#endif

#if defined(COMPILER_Intel8)
  logical, parameter :: compiler_intel8 = .true.
#else
  logical, parameter :: compiler_intel8 = .false.
#endif

#if defined(COMPILER_NAG)
  logical, parameter :: compiler_nag = .true.
#else
  logical, parameter :: compiler_nag = .false.
#endif

#if defined(COMPILER_XLF)
  logical, parameter :: compiler_xlf = .true.
#else
  logical, parameter :: compiler_xlf = .false.
#endif

#if defined(CONVERT_BIGENDIAN)
  logical, parameter :: convert_bigendian = .true.
#else
  logical, parameter :: convert_bigendian = .false.
#endif

#if defined(CUBED_SPHERE)
  logical, parameter :: cubed_sphere = .true.
#else
  logical, parameter :: cubed_sphere = .false.
#endif

#if defined(DEBUG)
  logical, parameter :: debug = .true.
#else
  logical, parameter :: debug = .false.
#endif

#if defined(DEBUG1)
  logical, parameter :: debug1 = .true.
#else
  logical, parameter :: debug1 = .false.
#endif

#if defined(DEBUG_DECOMP)
  logical, parameter :: debug_decomp = .true.
#else
  logical, parameter :: debug_decomp = .false.
#endif

#if defined(DEBUG_SNOW)
  logical, parameter :: debug_snow = .true.
#else
  logical, parameter :: debug_snow = .false.
#endif

#if defined(DEBUG_TR_WITH_WATER)
  logical, parameter :: debug_tr_with_water = .true.
#else
  logical, parameter :: debug_tr_with_water = .false.
#endif

#if defined(DEEP_ATM)
  logical, parameter :: deep_atm = .true.
#else
  logical, parameter :: deep_atm = .false.
#endif

#if defined(DOMAIN_DECOMP_ATM_IS_1D)
  logical, parameter :: domain_decomp_atm_is_1d = .true.
#else
  logical, parameter :: domain_decomp_atm_is_1d = .false.
#endif

#if defined(DO_EXPLIC_0)
  logical, parameter :: do_explic_0 = .true.
#else
  logical, parameter :: do_explic_0 = .false.
#endif

#if defined(DYNAMICS_ZS)
  logical, parameter :: dynamics_zs = .true.
#else
  logical, parameter :: dynamics_zs = .false.
#endif

#if defined(DYNAMIC_BIOMASS_BURNING)
  logical, parameter :: dynamic_biomass_burning = .true.
#else
  logical, parameter :: dynamic_biomass_burning = .false.
#endif

#if defined(ECOSYSTEM_SCALE)
  logical, parameter :: ecosystem_scale = .true.
#else
  logical, parameter :: ecosystem_scale = .false.
#endif

#if defined(EDGAR_1995)
  logical, parameter :: edgar_1995 = .true.
#else
  logical, parameter :: edgar_1995 = .false.
#endif

#if defined(EDGAR_HYDE_SOURCES)
  logical, parameter :: edgar_hyde_sources = .true.
#else
  logical, parameter :: edgar_hyde_sources = .false.
#endif

#if defined(EFLUX_OUT)
  logical, parameter :: eflux_out = .true.
#else
  logical, parameter :: eflux_out = .false.
#endif

#if defined(EIGHT_BYTE)
  logical, parameter :: eight_byte = .true.
#else
  logical, parameter :: eight_byte = .false.
#endif

#if defined(YIBS_1D_DIAG)
  logical, parameter :: yibs_1d_diag = .true.
#else
  logical, parameter :: yibs_1d_diag = .false.
#endif

#if defined(EVAP_VEG_GROUND)
  logical, parameter :: evap_veg_ground = .true.
#else
  logical, parameter :: evap_veg_ground = .false.
#endif

#if defined(EVAP_VEG_GROUND_NEW)
  logical, parameter :: evap_veg_ground_new = .true.
#else
  logical, parameter :: evap_veg_ground_new = .false.
#endif

#if defined(EXTEND_VG)
  logical, parameter :: extend_vg = .true.
#else
  logical, parameter :: extend_vg = .false.
#endif

#if defined(FIVE_AVG)
  logical, parameter :: five_avg = .true.
#else
  logical, parameter :: five_avg = .false.
#endif

#if defined(FVCUBED_SKIPPED_THIS)
  logical, parameter :: fvcubed_skipped_this = .true.
#else
  logical, parameter :: fvcubed_skipped_this = .false.
#endif

#if defined(FV_LAND)
  logical, parameter :: fv_land = .true.
#else
  logical, parameter :: fv_land = .false.
#endif

#if defined(GFDL_NUDGE)
  logical, parameter :: gfdl_nudge = .true.
#else
  logical, parameter :: gfdl_nudge = .false.
#endif

#if defined(GFED_3D_BIOMASS)
  logical, parameter :: gfed_3d_biomass = .true.
#else
  logical, parameter :: gfed_3d_biomass = .false.
#endif

#if defined(GHY_FD_1_HACK)
  logical, parameter :: ghy_fd_1_hack = .true.
#else
  logical, parameter :: ghy_fd_1_hack = .false.
#endif

#if defined(GHY_USE_LARGESCALE_PRECIP)
  logical, parameter :: ghy_use_largescale_precip = .true.
#else
  logical, parameter :: ghy_use_largescale_precip = .false.
#endif

#if defined(GLOBAL_TRIG)
  logical, parameter :: global_trig = .true.
#else
  logical, parameter :: global_trig = .false.
#endif

#if defined(HTAP_LIKE_DIAGS)
  logical, parameter :: htap_like_diags = .true.
#else
  logical, parameter :: htap_like_diags = .false.
#endif

#if defined(HYCOM_UNFINISHED)
  logical, parameter :: hycom_unfinished = .true.
#else
  logical, parameter :: hycom_unfinished = .false.
#endif

#if defined(INITIAL_GHG_SETUP)
  logical, parameter :: initial_ghg_setup = .true.
#else
  logical, parameter :: initial_ghg_setup = .false.
#endif

#if defined(INIT_4BYTE)
  logical, parameter :: init_4byte = .true.
#else
  logical, parameter :: init_4byte = .false.
#endif

#if defined(INTEL_OPT)
  logical, parameter :: intel_opt = .true.
#else
  logical, parameter :: intel_opt = .false.
#endif

#if defined(INTERACTIVE_WETLANDS_CH4)
  logical, parameter :: interactive_wetlands_ch4 = .true.
#else
  logical, parameter :: interactive_wetlands_ch4 = .false.
#endif

#if defined(INTERCEPT_TEMPORAL)
  logical, parameter :: intercept_temporal = .true.
#else
  logical, parameter :: intercept_temporal = .false.
#endif

#if defined(IRIX64)
  logical, parameter :: irix64 = .true.
#else
  logical, parameter :: irix64 = .false.
#endif

#if defined(IRRIGATION_ON)
  logical, parameter :: irrigation_on = .true.
#else
  logical, parameter :: irrigation_on = .false.
#endif

#if defined(Jprod_based_on_cocc)
  logical, parameter :: jprod_based_on_cocc = .true.
#else
  logical, parameter :: jprod_based_on_cocc = .false.
#endif

#if defined(Jprod_based_on_pp)
  logical, parameter :: jprod_based_on_pp = .true.
#else
  logical, parameter :: jprod_based_on_pp = .false.
#endif

#if defined(LARGE_SCALE_PRECIP_INTERCEPT)
  logical, parameter :: large_scale_precip_intercept = .true.
#else
  logical, parameter :: large_scale_precip_intercept = .false.
#endif

#if defined(LATLON_CORE)
  logical, parameter :: latlon_core = .true.
#else
  logical, parameter :: latlon_core = .false.
#endif

#if defined(LUS_VERT_ADV)
  logical, parameter :: lus_vert_adv = .true.
#else
  logical, parameter :: lus_vert_adv = .false.
#endif

#if defined(MACHINE_DEC)
  logical, parameter :: machine_dec = .true.
#else
  logical, parameter :: machine_dec = .false.
#endif

#if defined(MACHINE_Linux)
  logical, parameter :: machine_linux = .true.
#else
  logical, parameter :: machine_linux = .false.
#endif

#if defined(MACHINE_MAC)
  logical, parameter :: machine_mac = .true.
#else
  logical, parameter :: machine_mac = .false.
#endif

#if defined(MACHINE_SGI)
  logical, parameter :: machine_sgi = .true.
#else
  logical, parameter :: machine_sgi = .false.
#endif

#if defined(MAPL_MODE)
  logical, parameter :: mapl_mode = .true.
#else
  logical, parameter :: mapl_mode = .false.
#endif

#if defined(MARS_GCM)
  logical, parameter :: mars_gcm = .true.
#else
  logical, parameter :: mars_gcm = .false.
#endif

#if defined(MIRROR_V)
  logical, parameter :: mirror_v = .true.
#else
  logical, parameter :: mirror_v = .false.
#endif

#if defined(MPITYPE_LOOKUP_HACK)
  logical, parameter :: mpitype_lookup_hack = .true.
#else
  logical, parameter :: mpitype_lookup_hack = .false.
#endif

#if defined(MPI_DEFS_HACK)
  logical, parameter :: mpi_defs_hack = .true.
#else
  logical, parameter :: mpi_defs_hack = .false.
#endif

#if defined(NCASA2)
  logical, parameter :: ncasa2 = .true.
#else
  logical, parameter :: ncasa2 = .false.
#endif

#if defined(NEWDIAG)
  logical, parameter :: newdiag = .true.
#else
  logical, parameter :: newdiag = .false.
#endif

#if defined(NEW_IO)
  logical, parameter :: new_io = .true.
#else
  logical, parameter :: new_io = .false.
#endif

#if defined(NEW_IO_4STRAT)
  logical, parameter :: new_io_4strat = .true.
#else
  logical, parameter :: new_io_4strat = .false.
#endif

#if defined(NEW_VECT)
  logical, parameter :: new_vect = .true.
#else
  logical, parameter :: new_vect = .false.
#endif

#if defined(NON_CONSV_Q)
  logical, parameter :: non_consv_q = .true.
#else
  logical, parameter :: non_consv_q = .false.
#endif

#if defined(NO_FORCING)
  logical, parameter :: no_forcing = .true.
#else
  logical, parameter :: no_forcing = .false.
#endif

#if defined(NO_MASS_FLUX)
  logical, parameter :: no_mass_flux = .true.
#else
  logical, parameter :: no_mass_flux = .false.
#endif

#if defined(NO_WASHOUT_IN_CLOUDS)
  logical, parameter :: no_washout_in_clouds = .true.
#else
  logical, parameter :: no_washout_in_clouds = .false.
#endif

#if defined(NO_WIND)
  logical, parameter :: no_wind = .true.
#else
  logical, parameter :: no_wind = .false.
#endif

#if defined(NUDGE_ON)
  logical, parameter :: nudge_on = .true.
#else
  logical, parameter :: nudge_on = .false.
#endif

#if defined(O18_KINETIC_FRAC)
  logical, parameter :: o18_kinetic_frac = .true.
#else
  logical, parameter :: o18_kinetic_frac = .false.
#endif

#if defined(OBIO_ON_GARYocean)
  logical, parameter :: obio_on_garyocean = .true.
#else
  logical, parameter :: obio_on_garyocean = .false.
#endif

#if defined(OBIO_RAD_coupling)
  logical, parameter :: obio_rad_coupling = .true.
#else
  logical, parameter :: obio_rad_coupling = .false.
#endif

#if defined(OFFLINE)
  logical, parameter :: offline = .true.
#else
  logical, parameter :: offline = .false.
#endif

#if defined(OFFLINE_RUN)
  logical, parameter :: offline_run = .true.
#else
  logical, parameter :: offline_run = .false.
#endif

#if defined(OG2AG_OCEANS_BUNDLE)
  logical, parameter :: og2ag_oceans_bundle = .true.
#else
  logical, parameter :: og2ag_oceans_bundle = .false.
#endif

#if defined(OG2AG_TOC2SST_BUNDLE)
  logical, parameter :: og2ag_toc2sst_bundle = .true.
#else
  logical, parameter :: og2ag_toc2sst_bundle = .false.
#endif

#if defined(OLD_A2D)
  logical, parameter :: old_a2d = .true.
#else
  logical, parameter :: old_a2d = .false.
#endif

#if defined(OLD_RAYF)
  logical, parameter :: old_rayf = .true.
#else
  logical, parameter :: old_rayf = .false.
#endif

#if defined(PBL_E1)
  logical, parameter :: pbl_e1 = .true.
#else
  logical, parameter :: pbl_e1 = .false.
#endif

#if defined(PFT_MODEL_YIBS)
  logical, parameter :: pft_model_yibs = .true.
#else
  logical, parameter :: pft_model_yibs = .false.
#endif

#if defined(PHENOLOGY_DIAG)
  logical, parameter :: phenology_diag = .true.
#else
  logical, parameter :: phenology_diag = .false.
#endif

#if defined(PRINT_GHY_VARS)
  logical, parameter :: print_ghy_vars = .true.
#else
  logical, parameter :: print_ghy_vars = .false.
#endif

#if defined(PRINT_GRID)
  logical, parameter :: print_grid = .true.
#else
  logical, parameter :: print_grid = .false.
#endif

#if defined(PS_BVOC)
  logical, parameter :: ps_bvoc = .true.
#else
  logical, parameter :: ps_bvoc = .false.
#endif

#if defined(QS_TEST)
  logical, parameter :: qs_test = .true.
#else
  logical, parameter :: qs_test = .false.
#endif

#if defined(RAD_O3_GCM_HRES)
  logical, parameter :: rad_o3_gcm_hres = .true.
#else
  logical, parameter :: rad_o3_gcm_hres = .false.
#endif

#if defined(RAD_VEG_GROUND)
  logical, parameter :: rad_veg_ground = .true.
#else
  logical, parameter :: rad_veg_ground = .false.
#endif

#if defined(RESTRICT_LITTER_FLUX)
  logical, parameter :: restrict_litter_flux = .true.
#else
  logical, parameter :: restrict_litter_flux = .false.
#endif

#if defined(RIGHT_HAND)
  logical, parameter :: right_hand = .true.
#else
  logical, parameter :: right_hand = .false.
#endif

#if defined(RMUMAX_allcocco)
  logical, parameter :: rmumax_allcocco = .true.
#else
  logical, parameter :: rmumax_allcocco = .false.
#endif

#if defined(RUNTIME_NTM)
  logical, parameter :: runtime_ntm = .true.
#else
  logical, parameter :: runtime_ntm = .false.
#endif

#if defined(SCM)
  logical, parameter :: scm = .true.
#else
  logical, parameter :: scm = .false.
#endif

#if defined(SERIAL_MODE)
  logical, parameter :: serial_mode = .true.
#else
  logical, parameter :: serial_mode = .false.
#endif

#if defined(SET_FLAG)
  logical, parameter :: set_flag = .true.
#else
  logical, parameter :: set_flag = .false.
#endif

#if defined(SET_SOILCARBON_GLOBAL_TO_ZERO)
  logical, parameter :: set_soilcarbon_global_to_zero = .true.
#else
  logical, parameter :: set_soilcarbon_global_to_zero = .false.
#endif

#if defined(SHIFT_WEST)
  logical, parameter :: shift_west = .true.
#else
  logical, parameter :: shift_west = .false.
#endif

#if defined(SHINDELL_STRAT_CHEM)
  logical, parameter :: shindell_strat_chem = .true.
#else
  logical, parameter :: shindell_strat_chem = .false.
#endif

#if defined(SHINDELL_STRAT_EXTRA)
  logical, parameter :: shindell_strat_extra = .true.
#else
  logical, parameter :: shindell_strat_extra = .false.
#endif

#if defined(SKIP_TRACERS_RAD)
  logical, parameter :: skip_tracers_rad = .true.
#else
  logical, parameter :: skip_tracers_rad = .false.
#endif

#if defined(SKIP_TRACER_DIAGS)
  logical, parameter :: skip_tracer_diags = .true.
#else
  logical, parameter :: skip_tracer_diags = .false.
#endif

#if defined(SOA_DIAGS)
  logical, parameter :: soa_diags = .true.
#else
  logical, parameter :: soa_diags = .false.
#endif

#if defined(SOILCARB_SITE)
  logical, parameter :: soilcarb_site = .true.
#else
  logical, parameter :: soilcarb_site = .false.
#endif

#if defined(SPMD)
  logical, parameter :: spmd = .true.
#else
  logical, parameter :: spmd = .false.
#endif

#if defined(SULF_ONLY_AEROSOLS)
  logical, parameter :: sulf_only_aerosols = .true.
#else
  logical, parameter :: sulf_only_aerosols = .false.
#endif

#if defined(SUMROOTSCELL)
  logical, parameter :: sumrootscell = .true.
#else
  logical, parameter :: sumrootscell = .false.
#endif

#if defined(SW_DYNAMICS)
  logical, parameter :: sw_dynamics = .true.
#else
  logical, parameter :: sw_dynamics = .false.
#endif

#if defined(TAF_DOES_NOT_LIKE)
  logical, parameter :: taf_does_not_like = .true.
#else
  logical, parameter :: taf_does_not_like = .false.
#endif

#if defined(TEST2)
  logical, parameter :: test2 = .true.
#else
  logical, parameter :: test2 = .false.
#endif

#if defined(TEST_DB)
  logical, parameter :: test_db = .true.
#else
  logical, parameter :: test_db = .false.
#endif

#if defined(TEST_GWAVES)
  logical, parameter :: test_gwaves = .true.
#else
  logical, parameter :: test_gwaves = .false.
#endif

#if defined(TEST_MONO)
  logical, parameter :: test_mono = .true.
#else
  logical, parameter :: test_mono = .false.
#endif

#if defined(TEST_TRACER)
  logical, parameter :: test_tracer = .true.
#else
  logical, parameter :: test_tracer = .false.
#endif

#if defined(TEST_VAND2)
  logical, parameter :: test_vand2 = .true.
#else
  logical, parameter :: test_vand2 = .false.
#endif

#if defined(TES_LIKE_DIAGS)
  logical, parameter :: tes_like_diags = .true.
#else
  logical, parameter :: tes_like_diags = .false.
#endif

#if defined(THIS_PART_IS_NOT_READY)
  logical, parameter :: this_part_is_not_ready = .true.
#else
  logical, parameter :: this_part_is_not_ready = .false.
#endif

#if defined(TRACERS_AEROSOLS_Koch)
  logical, parameter :: tracers_aerosols_koch = .true.
#else
  logical, parameter :: tracers_aerosols_koch = .false.
#endif

#if defined(TRACERS_AEROSOLS_SOA)
  logical, parameter :: tracers_aerosols_soa = .true.
#else
  logical, parameter :: tracers_aerosols_soa = .false.
#endif

#if defined(TRACERS_AGE_OCEAN)
  logical, parameter :: tracers_age_ocean = .true.
#else
  logical, parameter :: tracers_age_ocean = .false.
#endif

#if defined(TRACERS_AMP)
  logical, parameter :: tracers_amp = .true.
#else
  logical, parameter :: tracers_amp = .false.
#endif

#if defined(TRACERS_AMP_M1)
  logical, parameter :: tracers_amp_m1 = .true.
#else
  logical, parameter :: tracers_amp_m1 = .false.
#endif

#if defined(TRACERS_AMP_M2)
  logical, parameter :: tracers_amp_m2 = .true.
#else
  logical, parameter :: tracers_amp_m2 = .false.
#endif

#if defined(TRACERS_AMP_M3)
  logical, parameter :: tracers_amp_m3 = .true.
#else
  logical, parameter :: tracers_amp_m3 = .false.
#endif

#if defined(TRACERS_AMP_M4)
  logical, parameter :: tracers_amp_m4 = .true.
#else
  logical, parameter :: tracers_amp_m4 = .false.
#endif

#if defined(TRACERS_AMP_M5)
  logical, parameter :: tracers_amp_m5 = .true.
#else
  logical, parameter :: tracers_amp_m5 = .false.
#endif

#if defined(TRACERS_AMP_M6)
  logical, parameter :: tracers_amp_m6 = .true.
#else
  logical, parameter :: tracers_amp_m6 = .false.
#endif

#if defined(TRACERS_AMP_M7)
  logical, parameter :: tracers_amp_m7 = .true.
#else
  logical, parameter :: tracers_amp_m7 = .false.
#endif

#if defined(TRACERS_AMP_M8)
  logical, parameter :: tracers_amp_m8 = .true.
#else
  logical, parameter :: tracers_amp_m8 = .false.
#endif

#if defined(TRACERS_ATM_ONLY)
  logical, parameter :: tracers_atm_only = .true.
#else
  logical, parameter :: tracers_atm_only = .false.
#endif

#if defined(TRACERS_Alkalinity)
  logical, parameter :: tracers_alkalinity = .true.
#else
  logical, parameter :: tracers_alkalinity = .false.
#endif

#if defined(TRACERS_COSMO)
  logical, parameter :: tracers_cosmo = .true.
#else
  logical, parameter :: tracers_cosmo = .false.
#endif

#if defined(TRACERS_DRYDEP)
  logical, parameter :: tracers_drydep = .true.
#else
  logical, parameter :: tracers_drydep = .false.
#endif

#if defined(TRACERS_DUST)
  logical, parameter :: tracers_dust = .true.
#else
  logical, parameter :: tracers_dust = .false.
#endif

#if defined(TRACERS_DUST_Silt4)
  logical, parameter :: tracers_dust_silt4 = .true.
#else
  logical, parameter :: tracers_dust_silt4 = .false.
#endif

#if defined(TRACERS_GASEXCH_land)
  logical, parameter :: tracers_gasexch_land = .true.
#else
  logical, parameter :: tracers_gasexch_land = .false.
#endif

#if defined(TRACERS_GASEXCH_land_CO2)
  logical, parameter :: tracers_gasexch_land_co2 = .true.
#else
  logical, parameter :: tracers_gasexch_land_co2 = .false.
#endif

#if defined(TRACERS_GASEXCH_ocean)
  logical, parameter :: tracers_gasexch_ocean = .true.
#else
  logical, parameter :: tracers_gasexch_ocean = .false.
#endif

#if defined(TRACERS_GASEXCH_ocean_CFC)
  logical, parameter :: tracers_gasexch_ocean_cfc = .true.
#else
  logical, parameter :: tracers_gasexch_ocean_cfc = .false.
#endif

#if defined(TRACERS_GASEXCH_ocean_CO2)
  logical, parameter :: tracers_gasexch_ocean_co2 = .true.
#else
  logical, parameter :: tracers_gasexch_ocean_co2 = .false.
#endif

#if defined(TRACERS_HETCHEM)
  logical, parameter :: tracers_hetchem = .true.
#else
  logical, parameter :: tracers_hetchem = .false.
#endif

#if defined(TRACERS_HYCOM_Ventilation)
  logical, parameter :: tracers_hycom_ventilation = .true.
#else
  logical, parameter :: tracers_hycom_ventilation = .false.
#endif

#if defined(TRACERS_MINERALS)
  logical, parameter :: tracers_minerals = .true.
#else
  logical, parameter :: tracers_minerals = .false.
#endif

#if defined(TRACERS_NITRATE)
  logical, parameter :: tracers_nitrate = .true.
#else
  logical, parameter :: tracers_nitrate = .false.
#endif

#if defined(TRACERS_OCEAN)
  logical, parameter :: tracers_ocean = .true.
#else
  logical, parameter :: tracers_ocean = .false.
#endif

#if defined(TRACERS_OCEAN_INDEP)
  logical, parameter :: tracers_ocean_indep = .true.
#else
  logical, parameter :: tracers_ocean_indep = .false.
#endif

#if defined(TRACERS_OCEAN_WATER_MASSES)
  logical, parameter :: tracers_ocean_water_masses = .true.
#else
  logical, parameter :: tracers_ocean_water_masses = .false.
#endif

#if defined(TRACERS_OM_SP)
  logical, parameter :: tracers_om_sp = .true.
#else
  logical, parameter :: tracers_om_sp = .false.
#endif

#if defined(TRACERS_ON)
  logical, parameter :: tracers_on = .true.
#else
  logical, parameter :: tracers_on = .false.
#endif

#if defined(TRACERS_OceanBiology)
  logical, parameter :: tracers_oceanbiology = .true.
#else
  logical, parameter :: tracers_oceanbiology = .false.
#endif

#if defined(TRACERS_QUARZHEM)
  logical, parameter :: tracers_quarzhem = .true.
#else
  logical, parameter :: tracers_quarzhem = .false.
#endif

#if defined(TRACERS_RADON)
  logical, parameter :: tracers_radon = .true.
#else
  logical, parameter :: tracers_radon = .false.
#endif

#if defined(TRACERS_SPECIAL_Lerner)
  logical, parameter :: tracers_special_lerner = .true.
#else
  logical, parameter :: tracers_special_lerner = .false.
#endif

#if defined(TRACERS_SPECIAL_O18)
  logical, parameter :: tracers_special_o18 = .true.
#else
  logical, parameter :: tracers_special_o18 = .false.
#endif

#if defined(TRACERS_SPECIAL_Shindell)
  logical, parameter :: tracers_special_shindell = .true.
#else
  logical, parameter :: tracers_special_shindell = .false.
#endif

#if defined(TRACERS_TERP)
  logical, parameter :: tracers_terp = .true.
#else
  logical, parameter :: tracers_terp = .false.
#endif

#if defined(TRACERS_WATER)
  logical, parameter :: tracers_water = .true.
#else
  logical, parameter :: tracers_water = .false.
#endif

#if defined(TRACERS_WATER_OLD)
  logical, parameter :: tracers_water_old = .true.
#else
  logical, parameter :: tracers_water_old = .false.
#endif

#if defined(TRAC_ADV_CPU)
  logical, parameter :: trac_adv_cpu = .true.
#else
  logical, parameter :: trac_adv_cpu = .false.
#endif

#if defined(TRDIAG_WETDEPO)
  logical, parameter :: trdiag_wetdepo = .true.
#else
  logical, parameter :: trdiag_wetdepo = .false.
#endif

#if defined(MELT_FRESH_SNOW_ON_WARM_GROUND)
  logical, parameter :: melt_fresh_snow_on_warm_ground = .true.
#else
  logical, parameter :: melt_fresh_snow_on_warm_ground = .false.
#endif

#if defined(UNFINISHED_CROPS_CODE)
  logical, parameter :: unfinished_crops_code = .true.
#else
  logical, parameter :: unfinished_crops_code = .false.
#endif

#if defined(UPWIND_HALOS)
  logical, parameter :: upwind_halos = .true.
#else
  logical, parameter :: upwind_halos = .false.
#endif

#if defined(USE_2D)
  logical, parameter :: use_2d = .true.
#else
  logical, parameter :: use_2d = .false.
#endif

#if defined(USE_CONST_ZINT)
  logical, parameter :: use_const_zint = .true.
#else
  logical, parameter :: use_const_zint = .false.
#endif

#if defined(USE_DATA_ZS)
  logical, parameter :: use_data_zs = .true.
#else
  logical, parameter :: use_data_zs = .false.
#endif

#if defined(USE_DD2D_UTILS)
  logical, parameter :: use_dd2d_utils = .true.
#else
  logical, parameter :: use_dd2d_utils = .false.
#endif

#if defined(USE_YIBS)
  logical, parameter :: use_yibs = .true.
#else
  logical, parameter :: use_yibs = .false.
#endif

#if defined(USE_ESMF)
  logical, parameter :: use_esmf = .true.
#else
  logical, parameter :: use_esmf = .false.
#endif

#if defined(USE_EXTEND_CUBE)
  logical, parameter :: use_extend_cube = .true.
#else
  logical, parameter :: use_extend_cube = .false.
#endif

#if defined(USE_FFTW)
  logical, parameter :: use_fftw = .true.
#else
  logical, parameter :: use_fftw = .false.
#endif

#if defined(USE_FVCORE)
  logical, parameter :: use_fvcore = .true.
#else
  logical, parameter :: use_fvcore = .false.
#endif

#if defined(USE_FV_Q)
  logical, parameter :: use_fv_q = .true.
#else
  logical, parameter :: use_fv_q = .false.
#endif

#if defined(USE_GISS)
  logical, parameter :: use_giss = .true.
#else
  logical, parameter :: use_giss = .false.
#endif

#if defined(USE_MPI)
  logical, parameter :: use_mpi = .true.
#else
  logical, parameter :: use_mpi = .false.
#endif

#if defined(USE_MPP)
  logical, parameter :: use_mpp = .true.
#else
  logical, parameter :: use_mpp = .false.
#endif

#if defined(USE_NORM_VECT)
  logical, parameter :: use_norm_vect = .true.
#else
  logical, parameter :: use_norm_vect = .false.
#endif

#if defined(USE_PBL_E1)
  logical, parameter :: use_pbl_e1 = .true.
#else
  logical, parameter :: use_pbl_e1 = .false.
#endif

#if defined(USE_PFUNIT)
  logical, parameter :: use_pfunit = .true.
#else
  logical, parameter :: use_pfunit = .false.
#endif

#if defined(USE_RADIATION_E1)
  logical, parameter :: use_radiation_e1 = .true.
#else
  logical, parameter :: use_radiation_e1 = .false.
#endif

#if defined(USE_SYSUSAGE)
  logical, parameter :: use_sysusage = .true.
#else
  logical, parameter :: use_sysusage = .false.
#endif

#if defined(VORT_ON)
  logical, parameter :: vort_on = .true.
#else
  logical, parameter :: vort_on = .false.
#endif

#if defined(WATER_MISC_GRND_CH4_SRC)
  logical, parameter :: water_misc_grnd_ch4_src = .true.
#else
  logical, parameter :: water_misc_grnd_ch4_src = .false.
#endif

#if defined(WATER_PROPORTIONAL)
  logical, parameter :: water_proportional = .true.
#else
  logical, parameter :: water_proportional = .false.
#endif

#if defined(WAVE_FORM)
  logical, parameter :: wave_form = .true.
#else
  logical, parameter :: wave_form = .false.
#endif

#if defined(WET_DEPO_Ina)
  logical, parameter :: wet_depo_ina = .true.
#else
  logical, parameter :: wet_depo_ina = .false.
#endif

#if defined(constCO2)
  logical, parameter :: constco2 = .true.
#else
  logical, parameter :: constco2 = .false.
#endif

#if defined(DO_TOPMODEL_RUNOFF)
  logical, parameter :: do_topmodel_runoff = .true.
#else
  logical, parameter :: do_topmodel_runoff = .false.
#endif

#if defined(pCO2_ONLINE)
  logical, parameter :: pco2_online = .true.
#else
  logical, parameter :: pco2_online = .false.
#endif



  contains

    subroutine initialize()


      
    end subroutine initialize

end module RunTimeControls_mod
