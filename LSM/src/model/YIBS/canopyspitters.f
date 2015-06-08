#ifdef YIBS_1D_DIAG
#define  DEBUG  1
#endif

      module biophysics !canopyspitters
!@sum Spitters (1986) canopy radiative transfer (sunlit/shaded leaves), and
!@+   Simpson's Rule (Price, Numerical Recipes) for canopy layering.
!@+   Call photosynthesis/conductance routines from other module to scale
!@+   up leaf-level fluxes to the canopy.

      use yibs_types
      use yibs_const
      use yibs_pfts
      use photcondmod, only : pscondleaf, ciMIN
      use FarquharBBpspar

      implicit none
      
      public photosynth_cond !This is interface for YIBS.
      !public canopyfluxes

!      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
!      real*8,parameter :: EPS=1.d-8   !Small, to prevent div zero.
      real*8,parameter :: IPARMINMOL=LOW_PAR_LIMIT  !umol m-2 s-1
      real*8,parameter :: O2frac=.20900 !fraction of Pa, partial pressure.

      type :: canraddrv
         !Canopy radiative transfer
         !@var sigma Leaf scattering coefficient (?unitless).
         real*8 :: sigma        !=0.2D0
         real*8  sqrtexpr  !This just calculated from sigma.
         !@var kdf Canopy extinction coeff. for diffuse radiation (unitless).
         real*8 :: kdf          !=0.71D0
         !@var rhor Canopy reflectivity (?unitless).
         real*8 rhor
         !@var kbl Canopy extinction coeff. for black leaves (unitless).
         real*8 kbl

         !Vegetation geometry, biology
         !UPDATE WITH CANOPY RADIATIVE TRANSFER SCHEME
         integer :: pft
         !real*8 :: leafalbedo
         real*8 :: canalbedo  
         real*8 :: LAI

         !Radiation
!         real*8 :: solarzen !Solar zenith angle (radians)
         real*8 :: CosZen !cos(solarzen) = sin(solarelev)
         real*8 :: I0df, I0dr   !Direct and diffuse incident PAR at top of canopy (umol m-2 s-1)
      end type canraddrv

      contains
!################## MAIN SUBROUTINE #########################################
      subroutine photosynth_cond(dtsec, pp)
!@sum photosynth_cond  Main routine to set up drivers and calculate 
!@sum canopy scale fluxes.
!@+   Version that calls Farquhar-Ball-Berry leaf biophysics.
!@+   Calculates photosynthesis, autotrophic respiration, H2O conductance,
!@+   looping through cohorts.
!@+   Inputs:  met drivers, radiation, Ca, Tcanopy
!@+   Outputs:  GPP, NPP, autotrophic respiration components
      use yibs_const
      use yibs_types
      use FarquharBBpspar !pspartype, psdrvtype
      use photcondmod, only : biophysdrv_setup, calc_Pspar, pspar
      use respauto_physio, only : Rdark, water_stress3
      use physutil, only:  QSAT
      use yibs_prognostic, only: leaf_fall, calc_npp_resp, secs_360

      implicit none

      real*8, intent(in) :: dtsec
      type(patch),pointer :: pp
      !----Local----------------!
      type(cohort),pointer :: cop
      type(psdrvtype) :: psdrvpar !Met biophysics drivers, except for radiation.
      real*8 :: ci_umol !umol/mol, Leaf internal CO2 
      real*8 :: ca_umol !umol/mol, Ambient air CO2
      real*8 :: TsurfK, TcanK, TsoilK, Pa !,rh
      real*8 :: CosZen !,betad
      real*8 :: IPAR            !Incident PAR 400-700 nm (W m-2)
      real*8 :: fdir            !Fraction of IPAR that is direct
      real*8 :: Gb !Leaf boundary layer conductance of water vapor(mol m-2 s-1)
      real*8 :: fdry_pft_eff ! pft-specific effective dry canoopy fraction   
      real*8 :: Anet,Atot,Atot0,Rd    !umol m-2 s-1
      real*8 :: Iemis ! Nadine's isoprene emission umol m-2 s-1
      real*8 :: Memis ! Nadine's monoterpene emission umol m-2 s-1
      real*8 :: GCANOPY,TRANS_SW ! Ci,NPP !,R_auto
      real*8 :: GCANOPYsum, Ciavg, GPPsum, NPPsum, R_autosum,
     &          R_rootsum, GPPsum0  !PK 5/15/07
      real*8 :: IPPsum,MTPsum, Resprsum, Resplsum, Respwsum, Respgsum
      real*8 :: molconc_to_umol
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8 :: o3_nmol !nmol/m3, O3 surf conc.
      real*8 :: mma
      real*8 :: mml
      real*8 :: FO3sum, dFO3sum, Focan, dFocan
#endif
      real*8 :: fsmc
      real*8 :: g_leaf


      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
        pp%TRANS_SW = 1.d0
        return
      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *,"photosynth_cond: wrong pft = ", pp%tallest%pft
        call stop_model("photosynth_cond: wrong pft",255)
      endif

      !* ZERO SOME OUTPUT VARIABLES AT PATCH LEVEL
      pp%TRANS_SW = 1.d0 !Case of zero LAI.
      !* Time-stepped outputs:  CNC, Ci, Qf.

      !* INITIALIZE SUMMARY OUTPUT VARIABLES *!
      GCANOPYsum = 0.d0
      Ciavg = 0.d0
      GPPsum = 0.d0
      NPPsum = 0.d0
      IPPsum = 0.d0
      MTPsum = 0.d0
      Resprsum = 0.d0
      Resplsum = 0.d0
      Respwsum = 0.d0
      Respgsum = 0.d0
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      FO3sum = 0.d0
      dFO3sum = 0.d0
#endif
      R_autosum = 0.d0
      R_rootsum = 0.d0

      !* SET UP DRIVERS *!
      !* Patch-level water stress only needed for Friend&Kiang conductance.
      !* Cohort-level water stress is used for cohort-level photosynthesis.
!      pp%betad = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
!     i     ,pp%fracroot(:)
!     i     ,pp%cellptr%fice(:), pfpar(pp%tallest%pft)%hwilt
!     o     , pp%betadl(:))

      !* Radiation drivers *!
      IPAR = pp%cellptr%IPARdir + pp%cellptr%IPARdif
      if (pp%cellptr%IPARdir.eq.0.d0) then
        fdir = 0.d0
      else
        fdir = pp%cellptr%IPARdir / IPAR
      endif
      CosZen = pp%cellptr%CosZen

      Pa = pp%cellptr%P_mbar * 100.d0

      !* Other photosynthesis drivers *!
      !Set up psdrvpar - pack in met drivers.
      Gb = pp%cellptr%Ch*pp%cellptr%U* Pa/
     &     (gasc*(pp%cellptr%TairC+KELVIN)) !m/s * N/m2 * mol K/J * 1/K = mol/m2/s
      molconc_to_umol = gasc * (pp%cellptr%TcanopyC + KELVIN)/Pa * 1d6
      ca_umol = pp%cellptr%Ca * molconc_to_umol  !Convert to umol/mol or ppm.
      ci_umol = 0.7d0*ca_umol !pp%cellptr%Ci * molconc_to_umol  !This is solved for is pscubic in FBBphotosynthesis.f.  Replace with dummy initialization.
      TcanK = pp%cellptr%TcanopyC + KELVIN
      TsurfK = pp%cellptr%TairC + KELVIN
      TsoilK = pp%cellptr%Soiltemp(1) + KELVIN
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      o3_nmol = pp%cellptr%O3s
      mma = (gasc*TsurfK)/Pa
      mml = (gasc*TcanK)/Pa
#endif
      call biophysdrv_setup(ca_umol,ci_umol,
     &     pp%cellptr%TcanopyC,Pa,
     &     min(1.d0,  !RH
     &     max(pp%cellptr%Qf,0.d0)/Qsat(TcanK,
     &     2500800.d0 - 2360.d0*(TsurfK-KELVIN),Pa/100.d0)),
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &     o3_nmol, mma, mml,
#endif
     &     psdrvpar)  !Equation for latent heat of vaporizaiton comes from ..?
!----------------------

      !* LOOP THROUGH COHORTS *!
      cop => pp%tallest
      do while (ASSOCIATED(cop))

        !* Assign vegpar

        if (cop%LAI.gt.0.d0) then
!SOILMOIST_OLD
          !KIM - water_stress3 uses Soilmoist as a satured fraction
          cop%stressH2O = water_stress3(cop%pft, N_DEPTH,  
     i          pp%cellptr%Soilmoist(:), 
     &          cop%fracroot, pp%cellptr%fice(:), cop%stressH2Ol(:))

          call calc_Pspar(dtsec,cop%pft,psdrvpar%Pa,psdrvpar%Tc
     i         ,O2frac*psdrvpar%Pa
     i         ,cop%stressH2O,cop%Sacclim,cop%llspan)

          call canopyfluxes(dtsec, cop%pft
     &         ,pp%albedo(1)
     &         ,pp%LAI
     &         ,cop%LAI,IPAR*4.55d0 !4.05 see canopyfluxes comments.
     &         ,CosZen,fdir
     &         ,Gb
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &         ,Focan,dFocan
#endif
     &         ,psdrvpar
     &         ,GCANOPY,Anet,Atot,Atot0,Rd !NOTE: Ci should be cohort level
     &         ,Iemis
     &         ,Memis
     &         ,TRANS_SW)       !NOTE:  Should include stressH2O.
!     &       ,if_ci)  

          if (pfpar(cop%pft)%leaftype.eq.BROADLEAF) then
            ! stomata on underside of leaves so max stomatal blocking = 0
            fdry_pft_eff = 1.d0
          elseif(pfpar(cop%pft)%leaftype.eq.NEEDLELEAF) then
            ! Bosveld & Bouten (2003) max stomatal blocking = 1/3
            fdry_pft_eff = 1.d0 - min(pp%cellptr%fwet_canopy, 0.333d0)
          else
             fdry_pft_eff = 1.d0
          endif

         !* Assign outputs to cohort *!
         !* Account for wet leaves with pp%cellptr%fwet_canopy.
          cop%GCANOPY = GCANOPY*fdry_pft_eff*(gasc*TsurfK)/Pa  !Convert mol-H2O m-2 s-1 to m/s
          cop%Ci = psdrvpar%ci * !ci is in mole fraction
     &         psdrvpar%Pa/(gasc * (psdrvpar%Tc+KELVIN)) !mol m-3
          cop%GPP = Atot * fdry_pft_eff * 0.012d-6 !umol m-2 s-1 to kg-C/m2-ground/s
          cop%GPP0 = Atot0 * fdry_pft_eff * 0.012d-6 !umol m-2 s-1 to kg-C/m2-ground/s
          cop%IPP = Iemis * 0.0600d-6 !umol m-2 s-1 to kg-C-isoprene/m2-ground/s
          cop%MTP = Memis * 0.1200d-6  ! umolm-2s-1 to KgC-monoterpene/m2-ground/s
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
          cop%FO3 = Focan  !  nmol O3 m-2 s-1 
          cop%dFO3 = dFocan  !  nmol O3 m-2 s-1 
#endif
          Anet = Atot - Rd 

          call leaf_fall(cop%pft, pp%cellptr%TairC, 
     &                   pp%cellptr%Soilmoist(:), fsmc, g_leaf)
          cop%g_leaf_ac = cop%g_leaf_ac + g_leaf*dtsec/secs_360
          
#ifdef ACTIVE_GROWTH
          call calc_npp_resp(cop%pft, cop%ht_p, cop%lai_p, psdrvpar%Tc,
     &                       fsmc, cop%gpp*cop%lai_p/(cop%lai+1.0d-6),
     &                       rd*0.012d-6*cop%lai_p/(cop%lai+1.0d-6),
     &                       cop%npp, cop%resp_w,
     &                       cop%resp_r, cop%resp_l, cop%resp_p)
#else
          call calc_npp_resp(cop%pft, cop%h, cop%lai, psdrvpar%Tc,
     &                       fsmc, cop%gpp, rd*0.012d-6,
     &                       cop%npp, cop%resp_w,
     &                       cop%resp_r, cop%resp_l, cop%resp_p)
#endif

          cop%npp_ac = cop%npp_ac + cop%npp*dtsec
          cop%resp_w_ac = cop%resp_w_ac + cop%resp_w*dtsec

       else                     !Zero LAI or no light
          cop%GCANOPY=0.d0 !May want minimum conductance for stems.
          cop%Ci = EPS
          cop%NPP = 0.d0
          cop%GPP = 0.d0
          cop%GPP0 = 0.d0
          cop%IPP = 0.d0
          cop%MTP = 0.d0
          cop%resp_r = 0.d0
          cop%resp_l = 0.d0
          cop%resp_w = 0.d0
          cop%resp_p = 0.d0
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
          cop%FO3 = 0.0d0  !  nmol O3 m-2 s-1 
          cop%dFO3 = 0.0d0  !  nmol O3 m-2 s-1 
#endif
          if (cop%LAI.eq.0.d0) then 
             TRANS_SW = 1.d0
          else !(IPAR*4.05.lt.LOW_LIGHT_LIMIT)) then
             TRANS_SW = 0.d0
          endif
          Rd = Rdark(pspar%Vcmax)*cop%LAI
       endif
        !* Update cohort respiration components, NPP
       cop%R_auto = cop%resp_p
       cop%R_root = cop%resp_r

!        call Allocate_NPP_to_labile(dtsec, cop)

        ! update total carbon
        cop%C_total = cop%C_total + cop%NPP*dtsec

        !* pp cohort flux summaries
        GCANOPYsum = GCANOPYsum + cop%GCANOPY
        Ciavg = Ciavg + cop%Ci*cop%LAI
        GPPsum = GPPsum + cop%GPP
        GPPsum0= GPPsum0+ cop%GPP0
        IPPsum = IPPsum + cop%IPP
        MTPsum = MTPsum + cop%MTP
        NPPsum = NPPsum + cop%NPP
        Resprsum = Resprsum + cop%resp_r
        Resplsum = Resplsum + cop%resp_l
        Respwsum = Respwsum + cop%resp_w
        Respgsum = Respgsum + cop%resp_p
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
        FO3sum  = FO3sum + cop%FO3
        dFO3sum = dFO3sum + cop%dFO3
#endif
        R_autosum = R_autosum + cop%R_auto
        R_rootsum = R_rootsum + cop%R_root  !PK 5/15/07

        !set values for debugging
        cop => cop%shorter
      end do

      !* Patch-level OUTPUTS *!
      pp%GCANOPY = GCANOPYsum
      if ( pp%LAI > 0.d0 ) then
        pp%Ci = Ciavg/pp%LAI
      else
        pp%Ci = 0.d0
      endif
      pp%GPP = GPPsum
      pp%GPP0= GPPsum0
      pp%IPP = IPPsum
      pp%MTP = MTPsum
      pp%NPP = NPPsum
      pp%resp_r = Resprsum
      pp%resp_l = Resplsum
      pp%resp_w = Respwsum
      pp%resp_p = Respgsum
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      pp%FO3 = FO3sum
      pp%dFO3 = dFO3sum
#endif
      pp%R_auto = R_autosum
      pp%R_root = R_rootsum
      pp%TRANS_SW = TRANS_SW ! looks like a hack which will not work for
                             ! multiple cohorts...

      end subroutine photosynth_cond


!---------------------------------------------------------------------------
      subroutine canopyfluxes(dt, pft
     i     ,canalbedo,LAIcanopy,LAIcohort,IPAR,CosZen,fdir
     i     ,Gb
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &     ,Focan,dFocan
#endif
     i     ,psdrvpar
     o     ,Gs,Anet,Atot,Atot0,Rd,Iemis,Memis,TRANS_SW)
!     i     ,if_ci)
!@sum canopyfluxes Calculates cohort photosynthesis and conductance with
!@sum Farqhuar et al. (1980) photosynthesis, Ball-Berry stomatal conductance,
!@um  and Spitters (1986, 1987) canopy radiation (sunlit, shaded leaves).
!@sum Integrates vertically over the canopy with Simpson's Rule.
!@sum Ci is updated at the canopy level using the canopy boundary layer
!@sum conductance as in Friend and Kiang (2005). 
!@+
!@+   Should be renamed cohortfluxes.
!@+   If PAR is not directly available, the following conversions may be used:
      !  From total shortwave (W m-2) to PAR (umol m-2 s-1) (Monteith & Unsworth):
      !          PAR(umol m-2 s-1) = 2.3(umol/J)*SW(W m-2)
      !  From PAR (W m-2) to PAR (umol m-2 s-1) (U.Maryland, Dept. of Met., PAR Project),
      !  suggest nominal 485 nm for conversion, which gives:
      !          PAR(umol m-2 s-1) = 4.05(umol/J) * PAR(W m-2)
      !  Dye, D.G. (2004) suggests a slightly different conversion:
      !          PAR(umol m-2 s-1) = 4.56(umol/J) * PAR(W m-2)
      
      implicit none
      real*8,intent(in) :: dt   !time step (seconds)
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: canalbedo !Whole-canopy albedo
      real*8,intent(in) :: LAIcanopy  !Whole-canopy LAI
      real*8,intent(in) :: LAIcohort  !cohort LAI
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
!      real*8,intent(in) :: solarzen !Solar zenith angle
      real*8,intent(in) :: CosZen !Cosine ofSolar zenith angle
      real*8,intent(in) :: fdir !Fraction of PAR that is direct
      real*8,intent(in) :: Gb   !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
!      real*8,intent(in) :: ca   !CO2 mole fraction at reference height (umol mol-1)
!      real*8,intent(in) :: Tc   !Canopy (foliage) temperature (Celsius)
!      real*8,intent(in) :: rh   !Relative humidity (fraction)
!      real*8,intent(in) :: Pa   !Atmospheric pressure (Pa)
      type(psdrvtype) :: psdrvpar !Photosynthesis drivers, except for radiation.
!      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction (umol mol-1)
      real*8,intent(inout) :: Gs !Canopy stomatal conductance of water vapor (mol m-2 s-1)
      real*8,intent(out) :: Anet !Leaf net photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Atot0 !Leaf gross photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Leaf respiration (umol/m2/s)
      real*8,intent(out) :: TRANS_SW !Transmittance of shortwave to ground surface.
      real*8,intent(out) :: Iemis ! Leaf isoprene emission Nadine (micromol m-2 s-1)
      real*8,intent(out) :: Memis ! Leaf monoterpene emission Nadine (micromol m-2 s-1
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8,intent(out) :: Focan ! Canopy O3 uptake flux (nmol m-2 s-1)
      real*8,intent(out) :: dFocan ! Canopy excess O3 uptake (nmol m-2 s-1)
#endif
!      integer,intent(in) :: if_ci !FLAG 0-don't calculate ci, 1-calculate ci
      
      !-----Local---------------------
      real*8 :: Gsint           !Gs canopy conductance from qsimp integration (mol/m2/s)
      real*8 :: Rdint           !Rd canopy foliage respiration form qsimp integration (umol/m2/s)
      real*8 :: Iint            !Iint Nadine's isoprene emission
      real*8 :: Mint            !Mint Nadine's monoterpene emission
      type(canraddrv) :: cradpar
      !real*8 :: Ciconc          !Leaf internal CO2 mole fraction calculate (umol mol-1)


      !Set up cradpar and calculate diffuse and direct PAR.
      call canopy_rad_setup(pft,CosZen,fdir,IPAR,
     &     LAIcanopy,canalbedo,cradpar) 

      !Calculate net photosynthesis, integrate vertically with Simpson's Rule.
      !call qsimp(cradpar%LAI,cradpar,ci,Tc,Pa,rh,Anet,Gsint) 
      !### LAIcanopy for radiation and LAIcohort for photosynthesis need to be distinguished.
      call qsimp(cradpar%LAI,cradpar,psdrvpar,Gb,Atot,Atot0,
     & Gsint,Rdint,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     & Focan,dFocan,
#endif
     & Iint,Mint) 

      Rd = Rdint
      Gs = Gsint
      Anet = Atot - Rd
      Iemis = Iint
      Memis = Mint

      call canopy_transmittance(TRANS_SW,CosZen,fdir,cradpar)
      !Return:  Ci, Gs, Anet
!#ifdef DEBUG        
!        write(92,*) CosZen,IPAR,cradpar,psdrvpar
!     &       ,Gb,Gsint,Gs,Atot,Anet,Rd,TRANS_SW
!#endif
      end subroutine canopyfluxes


!################## PHOTOSYNTHESIS #########################################

      subroutine photosynth_sunshd(
!@sum photosynth_sunshd  Calculates sunlit and shaded leaf fluxes.
!@auth N.Y.Kiang
      !Spitters parameters
     i     Lcum                 !Cumulative LAI from top of canopy (m2/m2)
     i     ,crp                !Canopy radiation parameters 
     i     ,psd                 !Photosynthesis met drivers
     i     ,Gb                  !Leaf boundary layer conductance (mol/m2/s)
     o     ,Aleaf              !Leaf Net assimilation of CO2 in layer (umol m-2 s-1)
     o     ,Aleaf0             !Leaf Net assimilation of CO2 in layer (umol m-2 s-1)
     o     ,gsleaf            !Leaf Conductance of water vapor in layer (mol m-2 s-1)
     o     ,Rdleaf             !Leaf respiration (umol m-2 s-1)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     o     ,Fo3leaf
     o     ,DFo3leaf
#endif
     o     ,Ileaf           ! Leaf isoprene emission (umol m-2 s-1)
     o     ,Mleaf)           ! Leaf monoterpene emission (umol m-2 s-1)  
      implicit none
      real*8,intent(in) :: Lcum
      type(canraddrv) :: crp
      type(psdrvtype) :: psd
      real*8,intent(in) :: Gb
      !real*8,intent(in):: cs,Tl,Pa,rh
      !real*8,intent(inout) :: ci
      real*8,intent(out) :: Aleaf !Flux for single leaf
      real*8,intent(out) :: Aleaf0 !Flux for single leaf
      real*8,intent(out) :: gsleaf !Conductance for single leaf
      real*8,intent(out) :: Rdleaf !Respiration for single leaf
      real*8,intent(out) :: Ileaf ! Isop emis for single leaf
      real*8,intent(out) :: Mleaf ! MTP emis for single leaf
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8,intent(out) :: Fo3leaf
      real*8,intent(out) :: DFo3leaf
#endif
      !------Local---------
      real*8 fsl  !Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 Isl !PAR absorbed by sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 Ish !PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 Asl !Anet from sunlit leaves (umol/m2/s)
      real*8 Ash !Anet from shaded leaves (umol/m2/s)
      real*8 gssl,gssh !Leaf conductance, sunlit,shaded (mol-H2O/m2/s)
      real*8 Rdsl, Rdsh
      real*8 Iemisl  ! Nadine's isop emis from sunlit leaves (umol/m2/s)
      real*8 Iemiss  ! Nadine's isop emis from shaded leaves (umol/m2/s)
      real*8 Memisl  ! Nadine's mtp emis from sunlit leaves (umol/m2/s)
      real*8 Memiss  ! Nadine's mtp emis from shaded leaves (umol/m2/s)
      integer :: sunlitshaded !1-sunlit, 2-shaded
      real*8 Asl0 !Anet from sunlit leaves (umol/m2/s)
      real*8 Ash0 !Anet from shaded leaves (umol/m2/s)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8 Fo3l, dFo3l
      real*8 Fo3s, dFo3s
#endif


      call canopy_rad(Lcum,crp,Isl,Ish,fsl)
!      write(997,*) 'Lcum,sigma, sqrtexpr,kdf,rhor,kbl,pft,canalbedo,
!     & LAI,solarzen,I0df,I0dr,Isl,Ish,fsl',Lcum,crp,Isl,Ish,fsl

      !Calculate photosynthesis and stomatal conductance.
!      write(991,*) 'sunlit'

      sunlitshaded = 1
      call pscondleaf(crp%pft,Isl,psd,Gb,gssl,Asl,Rdsl,sunlitshaded,
     & Iemisl,Memisl,Asl0)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      Fo3l  = psd%Fo3
      dFo3l = psd%dFo3
#endif
!      write(992,*) 'shaded'
      sunlitshaded = 2
      call pscondleaf(crp%pft,Ish,psd,Gb,gssh,Ash,Rdsh,sunlitshaded,
     & Iemiss,Memiss,Ash0)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      Fo3s  = psd%Fo3
      dFo3s = psd%dFo3
#endif
                   
      Aleaf  = fsl*Asl + (1.0d0 - fsl)*Ash
      Aleaf0 = fsl*Asl0 + (1.0d0 - fsl)*Ash0
      gsleaf = fsl*gssl + (1.0d0 - fsl)*gssh
      Rdleaf = fsl*Rdsl + (1.0d0 - fsl)*Rdsh
      Ileaf = fsl*Iemisl + (1.0d0 - fsl)*Iemiss
      Mleaf = fsl*Memisl + (1.0d0 - fsl)*Memiss
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      Fo3leaf  = fsl*Fo3l + (1.0d0 - fsl)*Fo3s
      DFo3leaf = fsl*dFo3l + (1.0d0 - fsl)*dFo3s
#endif

      end subroutine photosynth_sunshd

!################# CANOPY CONDUCTANCE ######################################

      subroutine Gs_bound(dt, LAI, Gsnew, Gsinout)  
!@sum Gs_bound Limit rate of change of Gs (umol m-2 s-1)
!@+   Call to this was commented out - Igor?
!@auth A.Friend
      ! Required change in canopy conductance to reach equilibrium (m/s).
      implicit none
      real*8,intent(in) :: dt, LAI
      real*8,intent(in) :: Gsnew !Canopy conductance of curr time step (mol m-2 s-1)
      real*8,intent(inout) :: Gsinout !Bounded canopy conductance (mol m-2 s-1)
      !---Local----!
      real*8 :: Gsold           !Canopy conductance at prev time step (mol m-2 s-1)
      real*8, parameter :: rhoH2O = 998.2 !Density of water (1000 kg/m^3)
      real*8, parameter :: MW_H2O = 18.015 !Molecular weight of water (g/mol)
      !real*8, parameter :: ghi = 0.006d0*rhoh2o*1000./MW_H2O !Upper limit of gs leaf (mol m-2 s-1)
      !real*8, parameter :: glo = 0.000001d0*rhoH2O*1000./MW_H2O !Lower limit of gs leaf (mol m-2 s-1), See Ball and Berry paper.
      real*8, parameter :: ghi = 333.0 !Conversion from 6 mm s-1 upper limit.(mol m-2 s-1)
      real*8, parameter :: glo = .015 !Temperate grassland. Korner (1994) (mol m-2 s-1)

      real*8 :: dGs, dGs_max 

      dGs=Gsnew-Gsinout
      Gsold = Gsinout
      Gsinout = Gsnew
      !nu Limit Gs change over timestep because of guard cell mechanics (m/s)
      dGs_max=dt*LAI*(ghi-glo)/1800.0D0
      if( dGs.gt.dGs_max) Gsinout = Gsold + dGs_max
      if(-dGs.gt.dGs_max) Gsinout = Gsold - dGs_max
      ! Biological limits of absolute Gs (m/s).
      if(Gsinout.gt.ghi*LAI) Gsinout=ghi*LAI
      if(Gsinout.lt.glo*LAI) Gsinout=glo*LAI

      end subroutine Gs_bound


!################# RADIATIVE TRANSFER ######################################
      subroutine canopy_rad_setup(
     i     pft, CosZen, fdir,IPAR,LAI,canalbedo,
     o     crp)
!@sum Calculate diffuse and direct PAR, canopy albedo, if unknown, from given
!@sum solar zenith angle, and other canopy radiation parameters.

      implicit none

! NOTES FOR OBTAINING PAR
! If only incident shortwave is available, convert 
!  total shortwave (W/m2) to PAR (umol/m2/s).
!  2.3 umol/J = shortwave to fraction that is PAR (Monteith & Unsworth).
!      PAR=2.3d0*max(srht,zero)/(1.0D0-0.08D0)
! Current: Replaced back-calculation with actual incident PAR at surface.

! GENERAL PAR CONVERSIONS.
! nyk  For 400-700 nm, energy is 3.3-6 umol photons/J.
!      Top-of-atmosphere flux-weighted average of energy is 4.54 umol/J.
!      Univ. Maryland, Dept. of Met., PAR Project, suggests nominal
!      485 nm for conversion, which gives 4.05 umol/J.

      !Input parameters
      integer,intent(in) :: pft
!      real*8,intent(in) :: solarzen !Solar zenith angle
      real*8,intent(in) :: CosZen !Cosine of Solar zenith angle
      real*8,intent(in) :: fdir   !Fraction of shortwave radiation that is direct
      real*8,intent(in) :: IPAR    !PAR incident at top of canopy (umol/m2/s)
      real*8,intent(in) :: LAI    !Leaf area index of whole canopy
      real*8,intent(in) :: canalbedo !Canopy albedo, if known

      !Output parameters
      type(canraddrv) :: crp
      !real*8,intent(out) :: I0df, I0dr !Diffuse and direct PAR incident at top of canopy (umol m-2 s-1)
      !----Local var-------
      real*8 :: sbeta  !Sine of solar zenith angle
      ! CONSTANTS
      real*8,parameter :: sigma=0.2d0  !Leaf scattering coefficient
      real*8,parameter :: kdf=0.71d0  !Canopy extinction coefficient
      real*8,parameter :: EPS=1.d-8   !Small, to prevent div zero.

      sbeta = CosZen ! cos(solarzen) = sin(solarelev)

!      crp%solarzen = acos(CosZen)
      crp%CosZen = CosZen
! Direct beam PAR incident on canopy (umol/m2/s).
      crp%I0dr=fdir*IPAR
! Diffuse PAR incident on canopy (umol/m2/s).
      crp%I0df=(1.0D0-fdir)*IPAR
! Canopy reflectivity depends on sun angle. Or take prescribed albedos.
      crp%sigma = sigma
      crp%kdf = kdf
      crp%sqrtexpr=sqrt(1.0D0-crp%sigma)
! Canopy extinction coefficient for direct beam radiation depends on
! sun angle (unitless).
      crp%kbl=0.5D0*crp%kdf/(0.8D0*crp%sqrtexpr*sbeta+EPS)

      !Hack canopy reflectivity.  This assumes a closed canopy!
      if (canalbedo.eq.0.0) then !Use if radiation gets initialized late.
        crp%rhor=((1.0D0-crp%sqrtexpr)/(1.0D0+crp%sqrtexpr))*
     &       (2.0D0/(1.0D0+1.6D0*sbeta))
        crp%canalbedo = crp%rhor
      else
       !Prescribed veg albedos until have canopy scheme.
        crp%rhor = canalbedo
        crp%canalbedo = crp%rhor
      end if

      crp%pft = pft
      crp%LAI = LAI

      end subroutine canopy_rad_setup

!-----------------------------------------------------------------------------
      subroutine canopy_rad(
!@sum Calculate PAR on sunlit and shaded foliage in layer Lc      
     i     Lc,                  !Cumulative LAI at layer from top of canopy
     i     crp,                !veg parameters
     o     Isla,                !PAR on shaded foliage in layer Lc
     o     Isha,                 !PAR on sunlit foliage in layer Lic
     o     fsl)                 !Fraction of sunlit foliage at Lc (unitless).

      implicit none
! Get incident radiation on layer, diffuse/direct, sunlit/shaded.
      real*8,intent(in) :: Lc
      type(canraddrv) :: crp
!@var Isla PAR absorbed by sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8,intent(out) :: Isla
!@var Isha PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s).
      real*8,intent(out) :: Isha
!@var fsl Fraction of sunlit foliage in layer at Lc (unitless).
      real*8,intent(out) :: fsl

      !-------Local vars--------------------------------------------------
      real*8 :: sbeta           !cos(solarzen) = sin(solarelev)
      real*8 :: I0df, I0dr !Diffuse and direct IPAR (umol m-2 s-1)
!@var Idfa Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idfa
!@var Idra Direct PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idra
!@var Idrdra Direct PAR at Lc that remains direct (umol/m[ground]2/s).
      real*8 Idrdra
      real*8,parameter :: zero = 0.d0
      real*8 :: Ish, Isl !Incident rather than absorbed PAR (umol/m2/s).

      sbeta = crp%CosZen !cos(solarzen) = sin(solarelev)
      I0df = crp%I0df
      I0dr = crp%I0dr
! When sun is above the horizon.
      if(sbeta.gt.EPS)then
!? Diffuse PAR in canopy at Lc (umol/m[ground]2/s). ?Old comment.
! Diffuse PAR that is absorbed at Lc (umol/m2-gnd/s). Eq. 10.
        Idfa=(1.0D0-crp%rhor)*I0df*crp%kdf *exp(-crp%kdf*Lc)
!? Direct PAR in canopy at Lc (umol/m[ground]2/s). ?Old comments.
! Direct PAR that is absorbed at Lc (umol/m2-gnd/s). Eq. 11.
        Idra=(1.0D0-crp%rhor)*I0dr* crp%sqrtexpr * crp%kbl *
     $        exp(-crp%sqrtexpr*crp%kbl*Lc)
!? Direct PAR at Lc that remains direct (umol/m[ground]2/s).?Old comment
! Non-scattered part of direct flux that is absorbed (umol/m2-gnd/s). Eq. 12.
        Idrdra=(1.0D0-crp%sigma)*I0dr * crp%kbl *
     $        exp(-crp%sqrtexpr*crp%kbl*Lc)
! PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s). Diffuse + Direct portion that is scattered (becomes diffuse)
        Isha=Idfa+(Idra-Idrdra)
! PAR aborbed by sunlit foliage at Lc (umol/m[foliage]2/s)._
        Isla=Isha+(1.0D0-crp%sigma)* crp%kbL * I0dr

!######## HACK - above are absorbed radiation, not incident ###########
!######## Replace with Spitters' equations for incident radiation #####
!        Isl = (1.0D0-crp%rhor)*I0dr*exp(-crp%kbL*crp%sqrtexpr*Lc)
!        Ish = (1.0D0-crp%rhor)*I0df*exp(-crp%kdf*Lc)
!        Isla = I0dr + Ish  !## HACK - this is to pass out Isl, using Isla var.
!        Isha = Ish !## HACK - this is to pass out Ish, using Isha var.
! Fraction of sunlit foliage in layer at Lc (unitless).
        fsl=exp(-crp%kbL*Lc)
        if(Isha.le.0.0D0)Isha=0.0D0
        if(Isla.le.0.0D0)Isla=0.0D0
        if(fsl.lt.zero)fsl=zero
      else
        Isha=0.0D0
        Isla=0.0D0
        fsl=zero
      endif

      end subroutine canopy_rad

      subroutine canopy_transmittance(TRANS_SW, sbeta,fdir,crp)
!@sum Calculate the transmittance (fraction) of radiation through the canopy
!@+   to the soil surface. Verbatim from Spitters et al. (1986,1987)     
!     Transmission of shortwave radiation through canopy. This has errors.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Parameters -------------------------------------------------------------
!@var trans_sw Total canopy transmittance of shortwave (fraction)
      real*8, intent(out) :: TRANS_SW  
!@var sbeta Cos of solar zenith angle
      real*8, intent(in) :: sbeta
!@var fdir Fraction of surface visible radiation that is direct
      real*8, intent(in) :: fdir
! Canopy radiation parameters
      type(canraddrv) :: crp

!Local variables --------------------------------------------------------
! Get assigned from crp passed in.
      real*8 sigma,sqrtexpr,rhor,kbl,kdf,alai
!@var absdf Absorbance of diffuse SW radiation (fraction)
      real*8 absdf
!@var absdr Absorbance of direct SW radiation (fraction)
      real*8 absdr
!@var absdrdr Absorbance of direct SW that remains direct (fraction)
      real*8 absdrdr
!@var abssh Absorbance of SW through shaded foliage (fraction)
      real*8 abssh
!@var abssl Absorbance of SW through sunlit foliage (fraction)
      real*8 abssl
!@var fracsl Fraction of leaves in canopy that are sunlit (fraction)
      real*8 fracsl
!------------------------------------------------------------------------
      sigma=crp%sigma
      sqrtexpr=crp%sqrtexpr
      rhor=crp%rhor
      kbl=crp%kbl
      kdf=crp%kdf
      alai=crp%LAI

! Diffuse shortwave in canopy (umol/m[ground]2/s).
      if(sbeta.gt.zero)then
        absdf=(1-fdir)*(1.0D0-rhor)*kdf*exp(-kdf*alai)
! Direct shortwave in canopy (umol/m[ground]2/s).
        absdr=fdir*(1.0D0-rhor)*sqrtexpr*kbl*exp(-sqrtexpr*kbl*alai)
! Direct shortwave that remains direct (umol/m[ground]2/s).
        absdrdr=absdr*(1.0D0-sigma)*kbl*exp(-sqrtexpr*kbl*alai)
! Shortwave penetrating shaded foliage (umol/m[foliage]2/s).
        abssh=absdf + (absdr-absdrdr) 
! Shortwave penetrating sunlit foliage (umol/m[foliage]2/s).
        abssl=abssh+(1.0D0-sigma)*kbL*fdir
        if(abssh.lt.0.0D0) abssh=0.0D0
        if(abssl.lt.0.0D0) abssl=0.0D0
        if(abssh.gt.1.0D0) abssh=1.0D0
        if(abssl.gt.0.0D0) abssl=1.0D0
        fracsl=exp(-kbl*alai)
      else
        abssh=0.001D0
        abssl=0.0D0
        fracsl=0
      endif

!      TRANS_SW = 1.d0 - ((1.d0-fracsl)*abssh + fracsl*abssl)
      TRANS_SW = ((1.d0-fracsl)*abssh + fracsl*abssl)

!      write(110,*) TRANS_SW,I0dr,I0df,abssh, abssl, sbeta, sigma, kbl


      end subroutine canopy_transmittance

!#############################################################################
!###################### INTEGRATION ##########################################

      subroutine qsimp(Xlim,crp,psp,Gb,S,S0,Sg,Sr,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     *   Sf, Sd,
#endif
     *   Si, Sm)
!----------------------------------------------------------------------!
!@sum qsimp Numerical routine to calculate canopy photosynthesis by 
!@+   increasing the number of
!@+   layers in the integration until the result (S) changes by less than
!@+    0.1 umol[CO2]/m2/s.
!@auth A.D.Friend
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Passed parameters
      real*8,intent(in) :: Xlim
      type(canraddrv) :: crp
      type(psdrvtype) :: psp
      real*8,intent(in) :: Gb
      !real*8,intent(in) :: ca, ci,T,P,RH
      real*8,intent(out) :: S,S0,Sg,Sr,Si,Sm
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8,intent(out) :: Sf, Sd
#endif
!----------------------------------------------------------------------!
!Local variables
!nu   real*8, parameter :: EPS=1.D-3
      integer, parameter :: MAXIT=6
      real*8,parameter :: ERRLIM=0.1d0
      integer IT, layers
      real*8 :: A, B, OST,OS,ST, ERROR
      real*8 :: Ac,Bc
      real*8 :: OST0, OS0, ST0 !NK
      real*8 :: OSTg, OSg, STg !NK
      real*8 :: OSTr, OSr, STr !NK
      real*8 :: OSTi, OSi, STi !NK + NU
      real*8 :: OSTm, OSm, STm !NU
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8 :: OSTf, OSf, STf
      real*8 :: OSTd, OSd, STd
#endif
!----------------------------------------------------------------------!
! Calculate canopy radiative transfer and photosynthesis in increasing 
!   number of layers.
! NOTE:  ERROR is only calculated on photosynthesis value (ST).
!        Corresponding conductance value is then returned also (STg).
!        and isoprene emission value (STi)
      A=0.0D0
      B=Xlim      !Canopy LAI.  Ideally, qsimp should not see inside crp.
      Ac=0.0d0
      Bc=crp%LAI  !Cohort LAI, assumes spans height. Later:make different heights.
      OST=-1.D30
      OS= -1.D30
      OST0=-1.D30
      OS0= -1.D30
      OSTg=-1.D30 !NK
      OSg= -1.D30 !NK
      OSTr=-1.D30 !NK
      OSr= -1.D30 !NK
      OSTi= -1.D30 !NK + NU
      OSi= -1.D30 !NK + NU
      OSTm= -1.D30 !NU
      OSm= -1.D30 !NU
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      OSTf= -1.D30 !NK + NU
      OSf= -1.D30 !NK + NU
      OSTd= -1.D30 !NK + NU
      OSd= -1.D30 !NK + NU
#endif
      layers=1
      do 11 IT=1,MAXIT
         CALL TRAPZD(A,B,Ac,Bc,ST,ST0,STg,STr,STi,STm,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     *     STf, STd,
#endif
     *     IT,layers,crp,psp,Gb)
         S=(4.D0*ST-OST)/3.D0
         S0=(4.D0*ST0-OST0)/3.D0
         Sg=(4.D0*STg-OSTg)/3.D0 !NK
         Sr=(4.D0*STr-OSTr)/3.D0 !NK
         Si=(4.D0*STi-OSTi)/3.D0 !NK+NU
         Sm=(4.D0*STm-OSTm)/3.D0 !NU
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         Sf=(4.D0*STf-OSTf)/3.D0 !NK+NU
         Sd=(4.D0*STd-OSTd)/3.D0 !NK+NU
#endif
         ERROR=ABS(S-OS)
         IF (ERROR.lt.ERRLIM) RETURN
         OS=S
         OST=ST
         OS0=S0
         OST0=ST0
         OSg=Sg !NK
         OSTg=STg !NK
         OSr=Sr !NK
         OSTr=STr !NK
         OSi=Si !NK + NU
         OSTi=STi !NK + NU
         OSm=Sm !NU
         OSTm=STm !NU
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         OSf=Sf !NK + NU
         OSTf=STf !NK + NU
         OSd=Sd !NK + NU
         OSTd=STd !NK + NU
#endif

         if(IT.gt.1) layers=layers*2
   11 enddo
      return
      end subroutine qsimp

!======================================================================
      subroutine trapzd(L1,L2,L1c,L2c,S,S0,Sg,Sr,Si,Sm,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     *     Sf, Sd,
#endif
     *     N,layers,crp,psp,Gb)
!----------------------------------------------------------------------!
!@sum Integrates canopy photosynthesis over canopy layers using Simpson's
!@+   Rule (Press et al., 19??).
!----------------------------------------------------------------------!

      implicit none
!----------------------------------------------------------------------!
      real*8 L1,L2 !LAI of whole canopy
      real*8 L1c,L2c,S,S0,Sg,Sr,Si,Sm !LAI and sums of cohort
      integer N,layers
      type(canraddrv) :: crp 
      type(psdrvtype) :: psp
      real*8,intent(in) :: Gb
      !-----Local------------------
      integer L
      real*8 DEL,X,SUM,SUM0,SUMg,SUMr,SUMi,SUMm,RCL
      real*8 DELc
      real*8 Aleaf1 !Mean net photosynthesis at layer bottom (umol[CO2]/m2/s).
      real*8 Aleaf2 !Mean net photosynthesis at layer top (umol[CO2]/m2/s).
      real*8 Dleaf1 !Mean net photosynthesis at layer bottom (umol[CO2]/m2/s).
      real*8 Dleaf2 !Mean net photosynthesis at layer top (umol[CO2]/m2/s).
      real*8 gleaf1, gleaf2 !Mean conductance at L1, L2 (umol[H2O]/m2/s).
      real*8 Rdleaf1,Rdleaf2 !Mean leaf respiration at L1, L2 (umol[CO2]/m2/s)
      real*8 Ileaf1,Ileaf2 ! mean isop emis at L1,L2 (umol[C]/m2/s)
      real*8 Mleaf1,Mleaf2 ! mean monoterpene emis at L1,L2 (umol[C]/m2/s)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8 Fo3leaf1, Fo3leaf2, Sf, SUMf
      real*8 dFo3leaf1, dFo3leaf2, Sd, SUMd
#endif

      if(N.eq.1)then
         call photosynth_sunshd(L1,crp,psp,Gb,Aleaf1,Dleaf1,
     * gleaf1,Rdleaf1,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     * Fo3leaf1, DFo3leaf1,
#endif
     * Ileaf1,Mleaf1)

         call photosynth_sunshd(L2,crp,psp,Gb,Aleaf2,Dleaf2,
     * gleaf2,Rdleaf2,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     * Fo3leaf2, DFo3leaf2,
#endif
     * Ileaf2,Mleaf2)
!         S=0.5D0*(L2-L1)*(Aleaf1+Aleaf2)
!         Sg=0.5d0*(L2-L1)*(gleaf1+gleaf2)
!         Sr=0.5d0*(L2-L1)*(Rdleaf1+Rdleaf2)
         S=0.5D0*(L2c-L1c)*(Aleaf1+Aleaf2) !Cohort LAI for leaf processes.
         S0=0.5D0*(L2c-L1c)*(Dleaf1+Dleaf2) !Cohort LAI for leaf processes.
         Sg=0.5d0*(L2c-L1c)*(gleaf1+gleaf2)
         Sr=0.5d0*(L2c-L1c)*(Rdleaf1+Rdleaf2)
         Si=0.5d0*(L2c-L1c)*(Ileaf1+Ileaf2)
         Sm=0.5d0*(L2c-L1c)*(Mleaf1+Mleaf2)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         Sf=0.5d0*(L2c-L1c)*(Fo3leaf1+Fo3leaf2)
         Sd=0.5d0*(L2c-L1c)*(dFo3leaf1+dFo3leaf2)
#endif
      else
         RCL=layers  ! convert to real*8
         DEL=(L2-L1)/RCL !Canopy LAI for radiative transfer.
         DELc=(L2c-L1c)/RCL !Cohort LAI for radiative transfer.
         X=L1+0.5D0*DEL
         SUM=0.D0
         SUM0=0.D0
         SUMg=0.d0
         SUMr=0.d0
         SUMi=0.d0
         SUMm=0.d0
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         SUMf=0.d0
         SUMd=0.d0
#endif
         do 11 L=1,layers
           call photosynth_sunshd(X,crp,psp,Gb,Aleaf1,Dleaf1,
     *  gleaf1,Rdleaf1,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     *  Fo3leaf1, DFo3leaf1,
#endif
     *  Ileaf1,Mleaf1)
           SUM = SUM + Aleaf1
           SUM0 = SUM0 + Dleaf1
           SUMg = SUMg + gleaf1
           SUMr = SUMr + Rdleaf1
           SUMi = SUMi + Ileaf1
           SUMm = SUMm + Mleaf1
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
           SUMf = SUMf + Fo3leaf1
           SUMd = SUMd + DFo3leaf1
#endif
           X=X+DEL
   11    continue
!         S=0.5D0*(S+(L2c-L1c)*SUM/RCL)
!         Sg=0.5D0*(Sg+(L2c-L1c)*SUMg/RCL)
!         Sr=0.5D0*(Sr+(L2c-L1c)*SUMr/RCL)
         S=0.5D0*(S+(L2-L1)*SUM/RCL)
         S0=0.5D0*(S0+(L2-L1)*SUM0/RCL)
         Sg=0.5D0*(Sg+(L2-L1)*SUMg/RCL)
         Sr=0.5D0*(Sr+(L2-L1)*SUMr/RCL)
         Si=0.5D0*(Si+(L2-L1)*SUMi/RCL)
         Sm=0.5D0*(Sm+(L2-L1)*SUMm/RCL)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
         Sf=0.5D0*(Sf+(L2-L1)*SUMf/RCL)
         Sd=0.5D0*(Sd+(L2-L1)*SUMd/RCL)
#endif
      endif
      return
      end subroutine trapzd



!============================================================================

      end module biophysics !canopyspitters
