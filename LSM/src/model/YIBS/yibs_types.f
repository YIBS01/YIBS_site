      module yibs_types
      use yibs_const

      implicit none

! the following two parameters are basically needed for integrity 
! checks (to avoid infinite loops when processing linked lists)
!@var MAX_PATCHES maximal number of patches per cell
!@var MAX_COHORTS maximal number of cohorts per patch
        integer, parameter :: MAX_PATCHES=32, MAX_COHORTS=64


!****************************************************************************
!*       TYPE DECLARATIONS
!****************************************************************************
      type timestruct
         !Structure for holding time variables
         integer :: year
         integer :: month
         integer :: jday  !day of year
         integer :: day   !day of month, not of year
         integer :: hour  !0-23
         integer :: minute 
         integer :: seconds
      end type timestruct
!****************************************************************************
      type pftype
         !Constant parameters specific to vegetation type
         integer :: pst ! Photosynth type 1=C3, 2=C4
         logical :: woody !Woody, FALSE=NO, TRUE=YES.
         integer :: leaftype !1=broadleaf, 2=needleleaf, 3=monocot (not crops)
         real*8 :: hwilt  !Wilting point matric potential (m)
         real*8 :: sstar  !Rel. soil moist at stress onset (Rodriguez-Iturbe)
         real*8 :: swilt  !Normalized soil water at wilting point (dim'less)
         real*8 :: nf !Canopy nitrogen factor (dimensionless) (Kull and Kruijt)
         !CASA parameters, CLM parameters
         real*8 :: sla !Specific leaf area (m^2 leaf area/kg C)
         real*8 :: r      !CLM respiration parameter (gC/gN)
         real*8 :: lrage !CASA Turnover time of leaves and roots (years)
         real*8 :: woodage !CASA Turnover time of stems (years)
         real*8 :: lit_C2N !CASA litcn_casa (C:N ratio) IS THIS FOLIAGE&ROOTS?
         real*8 :: lignin  !CASA lignin (UNITS?  lignin content of ??)
         real*8 :: croot_ratio !Coarse roots:Stem mass woody ratio

         !Phenology parameter - KIM
         !* Parameter for phenology 
         !phenotype - phenological type
         !            evergreen (1), 
         !            cold deciduous (2), 
         !            drought deciduous (3), 
         !            cold/drought deciduous (4),
         !            annual  (5)
         integer :: phenotype !phenological types
         !* Parameters for plant allomteries
         real*8 :: b1Cf !para 1 for allometric relation btw DBH & foliage C 
         real*8 :: b2Cf !para 2 for allometric relation btw DBH & foliage C
         real*8 :: b1Cd !para 1 for allometric relation btw DBH & dead C
         real*8 :: b2Cd !para 2 for allometric relation btw DBH & dead C
         real*8 :: b1Ht !para 1 for allometric relation btw DBH & height
         real*8 :: b2Ht !para 2 for allometric relation btw DBH & height

! PFT-specific Isoprene Paramaters

         real*8 :: Y_eps    ! PFT-specific fraction of electrons to make isoprene
         real*8 :: Y_eps_m  ! PFT-specific monoterpene basal emissions

      end type pftype
!****************************************************************************

!****************************************************************************
      type cohort
         integer :: pft           !* PFT number
         real*8 :: n              ! Density of individuals in cohort (#/m^2)
         type(ycelltype),pointer :: cellptr => null() !Pointer to yibs grid cell
         type(patch),pointer :: pptr    => null() !Pointer to patch
         type(cohort),pointer :: taller => null() !Pointer to next tallest cohort
         type(cohort),pointer :: shorter => null() !Pointer to next shortest cohort
         type(cohort),pointer :: csptaller => null() !Pointer to next taller conspecific
         type(cohort),pointer :: cspshorter => null() !Pointer to next shorter conspecfic

         !* PFT PARAMETERS
         ! Only need to index array of pftypes.

         !* NITROGEN status */
         !@var LAI Total cohort leaf area index (m2[leaf]/m2[ground])
         real*8 :: LAI               !* 
                                                           
         !* ALL QUANTITIES BELOW ARE FOR AN INDIVIDUAL *!
                                                           
         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* Height (m)
         real*8,pointer :: fracroot(:) => null() ! Fraction of roots in soil layer

         real*8,pointer :: height_dz(:) => null() ! height levels, equal dz
         real*8,pointer :: fp_dz(:) => null()     ! foliage profile, equal dz (could be empty)
         real*8,pointer :: height(:) => null()    ! height levels, same as patch
         real*8,pointer :: fp(:) => null()        ! Cumulative foliage profile at height (could be zero); 0 at top of canopy.
                                                           
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) !!Cohort level
         real*8 :: gcanopy        ! Conductance of water vapor/cohort (m/s)
         real*8 :: GPP            ! GPP flux/cohort/area cover (kg-C/m2/s)
         real*8 :: GPP0           ! GPP flux/cohort/area cover (kg-C/m2/s)
         real*8 :: IPP            ! Isoprene emission flux (kg-C/m2/s)
         real*8 :: MTP            ! Monoterpene emission flux (kg-C/m2/s)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         real*8 :: FO3            ! Ozone flux to stomata (nmol O3 m-2 s-1)
         real*8 :: dFO3           ! Excess ozone flux to stomata (nmol O3 m-2 s-1)
#endif
Cxyue
         real*8 :: ht_p           ! progostic vegetation height (m)
         real*8 :: ht_o           ! old progostic vegetation height (m)
         real*8 :: lai_p          ! progostic LAI (m2/m2)
         real*8 :: nstep          ! accumulative time step (days)
         real*8 :: phen_old       ! old phenology [0-1]
         real*8 :: g_leaf_ac      ! Acuumulative leaf turnover rate
         real*8 :: g_leafp_ac     ! Accumulative leaf turnover including leaf phenology
         real*8 :: npp_ac         ! Accumulative NPP (kg C /m2)
         real*8 :: resp_p         ! Plant respiration rate (kg C /m2/s)
         real*8 :: resp_l         ! Leaf respiration rate (kg C /m2/s)
         real*8 :: resp_r         ! Root respiration rate (kg C /m2/s)
         real*8 :: resp_w         ! Wood respiration rate (kg C /m2/s)
         real*8 :: resp_w_ac      ! Accumulative wood respiration (kg C /m2)
         real*8 :: C_leaf         ! Total Leaf carbon (kg C /m2)
         real*8 :: C_root         ! Total Root carbon (kg C /m2)
         real*8 :: C_wood         ! Total Wood carbon (kg C /m2)
         real*8 :: C_leaf_lit     ! Turnover of Leaf carbon (kg C /m2 /day)
         real*8 :: C_root_lit     ! Turnover of Root carbon (kg C /m2 /day)
         real*8 :: C_wood_lit     ! Turnover of Wood carbon (kg C /m2 /day)
         real*8 :: dC_leaf        ! Changes in Total Leaf carbon (kg C /m2 /day)
         real*8 :: dC_root        ! Changes in Total Root carbon (kg C /m2 /day)
         real*8 :: dC_wood        ! Changes in Total Wood carbon (kg C /m2 /day)
!         real*8 :: resp_p_m       ! Plant maintenance respiration (kg C /m2/s)
!         real*8 :: resp_p_g       ! Plant growth respiration (kg C /m2/s)
Cxyue

         real*8 :: NPP            ! NPP flux/cohort/area cover (kg-C/m2/s)  
         real*8 :: R_auto         ! Autotrophic respiration/cohort/area (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: R_root         ! Root respiration/cohort/area  (kg-C/m2/s) -PK 5/15/07
         !* PHENOLOGY - KIM             
         real*8 :: phenofactor   !phenofactor_c * phenofactor_d
         real*8 :: phenofactor_c !Cold deciduousness
         real*8 :: phenofactor_d !Drought deciduousness
         real*8 :: phenostatus          
         real*8 :: phenostatus_c
         real*8 :: phenostatus_d        
         real*8 :: betad_10d ! 10-day running average of betad
         real*8 :: CB_d !daily carbon balance
         real*8 :: turnover_amp         
         real*8 :: llspan               
         real*8 :: Sacclim ! state of acclimation/frost hardiness [deg C]
                                                           
         !* PHYSIOLOGICAL STATUS *!  !NYK
         real*8 :: stressH2O !* fraction stress factor, 0=stressed, 1=no stress
         real*8,pointer :: stressH2Ol(:) => null() !Water stress in layers.
         !* Additional C accounting     
         real*8 :: C_growth  !* Daily tissue growth respiration (kg-C/m2-cohort/day)
                             !*  Save C_growth to restart to distribute flux over the day.
                             ! - this is remainig respiration carbon 
         real*8 :: C_growth_flux ! respiration flux due to C_growth
         real*8 :: C_total   !* Hack to check for C balance.
      end type cohort

!****************************************************************************
      type patch
         real*8 :: age                !*Patch age (years)
         real*8 :: area               !*Patch area (fraction of ycell)
         type(ycelltype),pointer:: cellptr => null() !Pointer to grid cell
         type(patch),pointer :: older => null() !Pointer to next older patch
         type(patch),pointer :: younger=> null()  !Pointer to next younger patch
         type(cohort),pointer :: tallest => null() !Pointer to tallest cohort
         type(cohort),pointer :: shortest => null() !Pointer to shortest cohort

         !*- - - - - - - Cohorts summary variables - - - - - - - - - - -*!
         !  Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals or LAI.
         !  Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* DIAGNOSTICS - NITROGEN and LEAF status */
         !@var nm   Mean leaf nitrogen (g/m2[leaf])
         real*8 LAI               !*
         !@var LAIpft LAI by cover type.
         real*8,pointer :: LAIpft(:) => null() !(N_COVERTYPES)
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* mean Height (m)
         real*8,pointer :: HTpft(:) => null() !(N_COVERTYPES)
         real*8,pointer :: fracroot(:) => null() ! Fraction of roots in soil layer
         !* EXPORT - FLUXES (whole patch)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) 
         real*8 :: GCANOPY        ! Conductance of water vapor (m/s)
         !* DIAGNOSTICS - FLUXES
         real*8 :: GPP            ! GPP flux (kg-C/m2/s)
         real*8 :: GPP0           ! GPP flux (kg-C/m2/s)
         real*8 :: IPP            ! Isoprene flux (kg-C/m2/s)
         real*8 :: MTP            ! Monoterpene flux (kg-C/m2/s)
         real*8,pointer :: GPPpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: IPPpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: MTPpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: phenfpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: Cdead(:)  => null()!(N_COVERTYPES)
         real*8,pointer :: Clive(:)  => null()!(N_COVERTYPES)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         real*8 :: FO3            ! Ozone flux to stomata (nmol O3 m-2 s-1)
         real*8 :: dFO3           ! Excess ozone flux to stomata (nmol O3 m-2 s-1)
#endif
Cxyue
         logical :: l_area        !*true: allow area change
         real*8 :: area_p
         real*8 :: ht_p           ! progostic vegetation height (m)
         real*8 :: lai_p          ! progostic LAI (m2/m2)
Cxyue
         real*8 :: NPP            ! NPP flux (kg-C/m2/s)
         real*8 :: resp_p         ! Plant respiration rate (kg C /m2/s)
         real*8 :: resp_l         ! Leaf respiration rate (kg C /m2/s)
         real*8 :: resp_r         ! Root respiration rate (kg C /m2/s)
         real*8 :: resp_w         ! Wood respiration rate (kg C /m2/s)
         real*8 :: R_auto         ! Autotrophic respiration (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: R_root         ! Root respiration (kg-C/m2/s) -PK 5/15/07
         !- - - - - - - end of cohort summary variables - - - - - - - - - - - - -

         !- - - - - - - Patch total - - - - - - - - - - - - - - - - - - - - - - -

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED
         real*8 :: z0              !Roughness length, average over patch
         !* EXPORT 
         real*8 :: albedo(N_BANDS)   !Spectral albedo, average over patch
         real*8 :: betad             !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) => null() !Water stress in layers.
         real*8 :: TRANS_SW          !Transmittance of shortwave radiation (400-2500 nm) to the ground (fraction)
         real*8 :: CO2flux           !Net CO2 flux up (kg-C/m2-gnd/s)
         !* DIAGNOSTICS - soil
         real*8 :: Soil_resp         !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS) :: Tpool !(g-C/m^2, CASA Tpools, single cell) !added dim N_CASA_LAYERS -PK
         real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS) :: Tpool_o !(g-C/m^2)

         !* IMPORT - Variables calculated by GCM/EWB - downscaled from grid cell
         
      !use soil moisture (and temperature) for 2 CASA layers:    -PK 
      !0-30 cm, 30-100 cm (second might change to 30-200 cm)
      !**might change this and soiltemp to dynamically allocated arrays** -PK 7/07 
!         real*8 :: Soilmoist(N_CASA_LAYERS) !Soil moisture (volumetric fraction)   
         real*8 :: Soilmoist(N_DEPTH) !Soil moisture (rel. sat. vol. fraction) 

!         real*8 :: N_deposit    !N deposition (kgN/m2)

         !* Variables for biophysics and biogeochemistry
         
         !* Reproductive pools *!
         real*8,pointer :: Reproduction(:) => null() !Reproductive/seed pools array of length N_PFT (kgC/m2-patch) 

         !* Disturbance values
         real*8 :: fuel
         real*8 :: ignition_rate
         real*8 :: lambda1(T_SUB) !Site-averaged fire dist. rate during year
         real*8 :: disturbance_rate(N_DIST_TYPES)

         !* Soil data (needed for albedo computation)
         integer soil_type      ! 1 - sand (bright) ; 2 - dirt (dark)

         real*8 :: C_total
         real*8 :: C_growth
         real*8 :: carb_tot
         real*8 :: carb_soil
      end type patch


!****************************************************************************
      type ycelltype
Cxyue
         integer :: jday
         real*8 :: lat
         real*8 :: plantdate
         real*8 :: harvestdate
Cxyue
         real*8 :: area         !Area km^2
         type(patch), pointer:: youngest => null() 
         type(patch), pointer:: oldest => null() 

         !*- - - - - - - Cohorts summary variables - - - - - - - - - - -*!
         !* Per vegetated ground area of ycell ** excludes bare soil area.
         !  Intensive properties (e.g. geometry, LMA) are averages weighted by
         ! total number of individuals or leaf area.
         !  Extensive properties (e.g. biomass, Ntot) are totals per m2 ground

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED - NITROGEN and LEAF status */
         !@var nm   Mean leaf nitrogen (g/m2[leaf]).
         real*8 :: LAI               
         !@var LAIpft LAI by cover type.
         real*8,pointer :: LAIpft(:) => null()!(N_COVERTYPES)
         !real*8 :: LA            ! Leaf area (m2[leaf]/individual)

         !* IMPORT-PRESCRIBED, EXPORT-SIMULATED - GEOMETRY - trees:  GORT ellipsoids, grasses:leaf only
         real*8 :: h              !* mean Height (m)
         real*8,pointer :: HTpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: fracroot(:) => null() ! Fraction of roots in soil layer

         !* IMPORT/EXPORT PUBLIC - FLUXES)
         real*8 :: Ci             !*Internal foliage CO2 (mol/m3) 
         real*8 :: GCANOPY        ! Conductance of water vapor (m/s)
         !* EXPORT - FLUXES 
         real*8 :: GPP            ! GPP flux (kg-C/m2/s)
         real*8 :: GPP0           ! GPP flux (kg-C/m2/s)
         real*8 :: IPP            ! Isoprene flux (kg-C/m2/s)
         real*8 :: MTP            ! Monoterpene flux (kg-C/m2/s)
         real*8,pointer :: GPPpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: IPPpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: MTPpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: phenfpft(:) => null()!(N_COVERTYPES)
         real*8,pointer :: Cdead(:)  => null()!(N_COVERTYPES)
         real*8,pointer :: Clive(:)  => null()!(N_COVERTYPES)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         real*8 :: FO3            ! Ozone flux to stomata (nmol O3 m-2 s-1)
         real*8 :: dFO3           ! Excess ozone flux to stomata (nmol O3 m-2 s-1)
#endif
Cxyue
         real*8 :: ht_p           ! prognostic vegetation height (m)
         real*8 :: lai_p          ! prognostic LAI (m2/m2)
Cxyue
         real*8 :: NPP            ! NPP flux (kg-C/m2/s)
         real*8 :: R_auto         ! Autotrophic respiration (kg-C/m2/s)
                                  ! = growth(Acan) + maint(fol,sapwood,root)
         real*8 :: R_root         ! Root respiration (kg-C/m2/s) -PK 5/15/07
         real*8 :: resp_p         ! Plant respiration rate (kg C /m2/s)
         real*8 :: resp_l         ! Leaf respiration rate (kg C /m2/s)
         real*8 :: resp_r         ! Root respiration rate (kg C /m2/s)
         real*8 :: resp_w         ! Wood respiration rate (kg C /m2/s)
         !- - - - - - - end of cohort summary variables - - - - - - - - - - - - - - - - - - -

         !- - - - - -  Patch-level summary values - PHYSICAL ------------------
         !* EXPORT - from radiative transfer
         real*8 :: z0              !Roughness length, average over patch
         real*8 :: albedo(N_BANDS) !Albedo may be in bands or hyperspectral
         real*8 :: betad             !Water stress  # CALC FROM Soilmoist & SSTAR by PFT
         real*8,pointer :: betadl(:) => null() !Water stress in layers.
         real*8 :: TRANS_SW  !Transmittance of shortwave radiation to the ground (fraction)
         real*8 :: CO2flux           !Net CO2 flux up (kg-C/m2-gnd/s)
         !* DIAGNOSTICS - soil
         real*8 :: Soil_resp         !soil resp flux (kg-C/m2/s) -PK 6/14/06,changed umol to kg-NK 07/28/06
         real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS) :: Tpool !(g-C/m^2, CASA Tpools, single cell) !added dim N_CASA_LAYERS -PK

         !- - - - - - Ycell-level variables - - - - - - - - - - - - - - - -
         !* IMPORT/EXPORT - vegetation
         real*8 :: fv            ! vegetation fraction of ycell
         real*8 :: heat_capacity ! total veg. heat capacity
         real*8 :: fwet_canopy   ! fraction of canopy that is wet
         !* IMPORT - SOIL
         real*8 :: soil_Phi      !Soil porosity (m3/m3)
         real*8 :: soil_dry     !Soil wetness "when evapotranspiration stops"
         real*8 :: soildepth    !Soil depth (m)
         real*8 :: theta_max    !Saturated soil water volume (m/m)
         real*8 :: k_sat        !Saturated hydraulic conductivity
         real*8 :: root_Phi     !Infiltration factor promoted by roots (units?)

         !SOIL - CONSTANTS
         !Soil textures for CASA -PK
         real*8 :: soil_texture(N_SOIL_TEXTURES) ! fractions of soil textures
!         real*8 clayfrac  !fractional clay content (passed from GHY.f)
!         real*8 sandfrac  !fractional sand content (also from GHY.f)


         !IMPORT - METEOROLOGICAL STATE VARIABLES
         !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
         real*8 :: TairC ! Air temperature (Clesius) !KIM-to drive phenology
         real*8 :: TcanopyC     !Canopy temperatue (Celsius)
         real*8 :: Qf           !*Foliage surface vapor mixing ratio (kg/kg)
         real*8 :: P_mbar       !Atmospheric pressure (mb)
         real*8 :: Ca           !@Atmos CO2 conc at surface height (mol/m3).
         real*8 :: O3s          !Surf. O3 (nmol m-3)

         !next two now explicitly depth-structured (see above) -PK
!         real*8 :: Soilmoist(N_CASA_LAYERS) !Soil moisture (volumetric fraction)
!         real*8 :: Soiltemp(N_CASA_LAYERS)  !Soil temperature (Celsius)
         real*8 :: Soilmoist(N_DEPTH) !Soil moisture (volumetric fraction)
         real*8 :: Soiltemp(N_DEPTH)  !Soil temperature (Celsius)
         real*8,pointer :: Soilmp(:) => null() !Soil matric potential (m)
         real*8,pointer :: fice(:) => null() !Fraction of soil layer that is ice
         real*8 :: Ch           !Ground to surface heat transfer coefficient 
         real*8 :: U            !Surface layer wind speed (m s-1)


         !Radiation - IMPORT STATE VARIABLES
         !may later be broken down into hyperspectral increments.
         ! in an array
!         real*8 :: Ivis          !Incident visible  (W m-2)
!         real*8 :: Idir          !Incident direct visible  (W m-2)
!         real*8 :: IPAR         !Incident PAR 400-700 nm (W m-2)
         real*8 :: IPARdir        !Incident direct PAR (W m-2)
         real*8 :: IPARdif        !Incident diffuse PAR (W m-2)
         real*8 :: CosZen         !cos of solar zenith angle

         !PHENOLOGY - KIM
         real*8 :: soiltemp_10d
         real*8 :: airtemp_10d
         real*8 :: paw_10d          
         real*8 :: par_10d
         real*8 :: gdd
         real*8 :: sgdd
         real*8 :: fdd
         real*8 :: sfdd
         real*8 :: ncd
         real*8 :: nmd
         real*8 :: nmsd
         real*8 :: daylength(2) !previous & present day 
         integer :: fall
         integer :: fall_old
         real*8  :: soiltemp_ltm
         integer :: jday_old
         integer :: jday_tot
         integer :: drought_state

         ! diags and hacks
         real*8 :: C_total
         real*8 :: C_growth
         real*8 :: carb_tot
         real*8 :: carb_soil

!!! hack !!! - just to try master phenology
         real*8 :: ld
         real*8 :: light
      end type ycelltype


!****************************************************************************
!      type yibsdatatype
!        real longmin, longmax, latmin, latmax
!        integer longi,latj
!        type(timestruct),pointer :: tt      !Greenwich Mean Time
!        type (ycelltype),pointer :: grid(:,:)
!      end type yibsdatatype

!****************************************************************************
      type yibs_config
      ! this type should contain all parameters that describe the run
      ! i.e. flags, array dimensions etc. They assumed to be constant
      ! during the run but may change from run to run
        logical do_soilresp       ! do soil respiration
        logical do_phenology_activegrowth
        logical do_structuralgrowth
        logical do_frost_hardiness
        logical do_patchdynamics
        logical do_init_geo
Cxyue
        logical do_init_activegrowth
Cxyue
!        logical mixed_veg
      end type yibs_config


      end module yibs_types
