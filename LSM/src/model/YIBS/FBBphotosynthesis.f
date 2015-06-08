#include "rundeck_opts.h"
!#define DEBUG 1
!#define USE_NR_SOLVER_FOR_FBB
      module  photcondmod
!@sum Photosynthesis and stomatal conductance at the leaf level at the
!@+   physical time step.
!@+   Photosynthesis from Farquhar and von Caemmerer (1982) and
!@+   Stomatal conductance from Ball and Berry (1985, 1987).
!@+   Cubic solution to coupled equations by I. Aleinov, this module.
!@+
!@+   Main routine pscondleaf is called by the canopy biophysics module
!@+   which scales up fluxes from the leaf to the canopy level.
!@+
!@+   Physiological status routines, frost_hardiness and par-phenology, are
!@+   necessary for seasonal variation in leaf photosynthetic capacity in
!@+   biophysics-only runs, and are also used by the phenology module when
!@+   prognostic phenology is run.
      
      use FarquharBBpspar
      use yibs_const

      implicit none

      save

      public init_ci, pscondleaf, biophysdrv_setup,calc_Pspar,ciMIN
      !,fbb_night
      public frost_hardiness, par_phenology
      public photosynthpar, pspar

      !=====CONSTANTS=====!
      real*8,parameter :: ciMIN = 1.d-8  !Small error
!      real*8,parameter :: Kelvin = 273.15d0
      real*8,parameter :: Rgas = gasc !8.314510d0 !gas constant (8.314510 J/mol K)
      real*8,parameter :: Kc = 30.0       !Michaelis-Menten coeff.constant for CO2 (Pa), Collatz (1991)
      real*8,parameter :: Ko = 3.d4       !Michaelis-Menten coeff.constant for O2 (Pa), Collatz (1991)
      real*8,parameter :: KcQ10 = 2.1d0   !Kc Q10 exponent, Collatz (1991)
      real*8,parameter :: KoQ10 = 1.2d0   !Ko Q10 exponent, Collatz (1991)

      !=====DECLARED TYPES======!
      type photosynthpar       !Calculated values from pft-dependent pspartypes
      integer :: pft           !Plant functional type.  1-C3 grassland
      real*8 :: PARabsorb      !Leaf PAR absorptance (fraction)
      real*8 :: Vcmax          !Maximum photosynthetic capacity (umol m-2 s-1)
      real*8 :: RVcmax         !Maximum photosynthetic capacity for Respiration (umol m-2 s-1)
      real*8 :: Kc              !Michaelis-Menten parameter for CO2 (Pa)
      real*8 :: Ko              !Michaelis-Menten parameter for O2 (Pa)
      real*8 :: Gammastar       !CO2 compensation point (Pa)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
      real*8 :: Nleaf           !g-N/m^2[leaf] - May want to take this from Tpool instead.
      real*8 :: stressH2O       !Water stress factor (fraction, 1=no stress)
      logical :: first_call     !For optimizing run
      real*8 :: Ac              !Save Ac for calculation only once per timestep.
      real*8 :: As              !Save As for calculation only once per timestep.
      logical :: reset_ci_cubic1
      end type photosynthpar

      private
      !=====GLOBAL VARIABLES (MODULE ONLY)=====!
      type(photosynthpar) :: pspar !Photosynthesis parameters.
!-----------------------------------------------------------------------------

      contains

      subroutine init_ci(ca, ci)
!@sum init_ci  Initialize leaf internal CO2 concentration.
      implicit none
      real*8,intent(in) :: ca   !Ambient air CO2 concentration (umol mol-1)
      real*8,intent(inout) :: ci !Leaf internal CO2 mole fraction  (umol mol-1)

      !ci should be initialized. For rule of thumb, initialize to typically
      !observed ratio:
      ci = 0.7d0*ca

      end subroutine init_ci
      
!-----------------------------------------------------------------------------
      subroutine biophysdrv_setup(ca,ci,Tc,Pa,rh,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &  o3,mol_to_ma, mol_to_ml,
#endif
     &  psdrvpar)

!@sum Set up met drivers for photosynthesis at every physical time step.
      implicit none
      real*8,intent(in) :: ca, ci, Tc, Pa, rh
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8,intent(in) :: o3
      real*8,intent(in) :: mol_to_ma
      real*8,intent(in) :: mol_to_ml
#endif
      type(psdrvtype),intent(out) :: psdrvpar

      psdrvpar%ca = ca
      psdrvpar%ci = ci
      psdrvpar%Tc = Tc
      psdrvpar%Pa = Pa
      psdrvpar%rh = rh

#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      psdrvpar%o3 = o3
      psdrvpar%mol_to_ma = mol_to_ma
      psdrvpar%mol_to_ml = mol_to_ml
      psdrvpar%FO3  = 0.0d0
      psdrvpar%dFO3 = 0.0d0
#endif

      end subroutine biophysdrv_setup
!-----------------------------------------------------------------------------
      subroutine pscondleaf(pft,IPAR,psd,Gb,gsout,Aout,Rdout
     &     ,sunlitshaded,ISPout,MTPout,Atot0)
!@sum pscondleaf  Main routine to obtain leaf-level photosynthesis, 
!@+   stomatal conductance, and dark respiration.  
!@+   See Photosynth_analyticsoln for units.
      implicit none
      integer,intent(in) :: pft
      real*8,intent(in) :: IPAR !umol m-2 s-1. Absorbed PAR. Should APAR.
      type(psdrvtype) :: psd
      real*8,intent(in) :: Gb !mol m-2 s-1
      real*8,intent(out) :: gsout, Aout, Rdout, Atot0 !ci recorded in psd
      real*8,intent(out) :: ISPout,MTPout
      integer,intent(in) :: sunlitshaded
      !---Local---
      real*8 :: ci!, cs
      real*8 :: Fo3   ! ozone flux to stomata (nmol O3 m-2 s-1)
      real*8 :: dFo3  ! excessive ozone flux to stomata (nmol O3 m-2 s-1)
      real*8,parameter :: LOW_LIGHT_LIMIT = 2.5d0 !umol m-2 s-1.  Nobel 1999, lower light limit for green plants is 0.7 W m-2 ~ 3 umol m-2 s-1.
      
c      if (IPAR.lt.LOW_LIGHT_LIMIT) then
c        Rdout = Respveg(pftpar(pft)%Nleaf,psd%Tc)  !Should be only leaf respiration!
c        Aout = 0.d0
c        cs = ca - (Aout-Rdout)*1.37d0/Gb
c        gsout = pftpar(pft)%b
c        psd%ci = ca             !Dummy assignment, no need to solve for ci 
c      else
cddd      print *,"called Photosynth_analyticsoln",
cddd     &     pft,IPAR,psd%ca,ci,
cddd     &     psd%Tc,psd%Pa,psd%rh,Gb,gsout,Aout,Rdout,sunlitshaded

        call Photosynth_analyticsoln(pft,IPAR,psd%ca,ci,
     &     psd%Tc,psd%Pa,psd%rh,Gb,gsout,Aout,Atot0,Rdout,sunlitshaded,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     o     psd%o3, psd%mol_to_ma, psd%mol_to_ml, fo3, dfo3,
#endif
     &     ISPout,MTPout)
        psd%ci = ci             !Ball-Berry:  ci is analytically solved.  F-K: ci saved between time steps.
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
        psd%fo3 = fo3
        psd%dFo3 = dfo3
#endif


c      endif
        
cddd      !Biological limits for gs - cuticular conductance?
cddd      if(gsout.lt.(0.00006d0*psd%Pa/(gasc*(psd%Tc+KELVIN)))) then
cddd        gsout=0.00006d0*psd%Pa/(gasc*(psd%Tc+KELVIN))
cddd      endif

      end subroutine pscondleaf

!-----------------------------------------------------------------------------

!      subroutine fbb_night(Atot,Gs,Rd,Iso)
!      real*8, intent(out) :: Atot,Gs,Rd,Iso
!
!      Atot = 0.d0
!      Gs = BallBerry(0.d0, 1.d0, 1.d0, pspar)
!      Rd = 0.015d0 * pspar%Vcmax
!      Iso = 0.d0
!
!      end subroutine fbb_night

      subroutine Photosynth_analyticsoln(pft,IPAR,ca,ci,Tl,Pa,rh,gb,
     o     gs,Atot,Atot0,Rd,sunlitshaded
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     o    ,o3, mol_to_ma, mol_to_ml, Fo3, dFo3
#endif
     o    ,isp,mtp)
!@sum Photosynth_analyticsoln  Selects the correct root for the 
!@+   solution to the cubic coupled equations of  Farquhar photosynthesis and 
!@+   Ball-Berry conductance, including dark respiration.  
!@+   ci is solved for analytically for each of the limiting cases.
!@+   Outputs gs, Atot, Rd. May output also other VOC fluxes.
!@auth  N.Y.Kiang, I.Aleinov
      use yibs_pfts, only : pfpar, GRASSC4
      use respauto_physio, only : Rdark
      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: IPAR !Absorbed PAR.  WRONG OLD COMMENT:Incident PAR (umol m-2 s-1) 
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)      
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      integer,intent(in) :: sunlitshaded !For diagnostic outputs only.
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8,intent(in) :: o3
      real*8,intent(in) :: mol_to_ma
      real*8,intent(in) :: mol_to_ml
      real*8,intent(out) :: Fo3  ! ozone flux to stomata (nmol O3 m-2 s-1)
      real*8,intent(out) :: dFo3  ! excessive ozone flux to stomata (nmol O3 m-2 s-1)
#endif
      real*8,intent(out) :: ci   !Leaf internal CO2 concentration (umol mol-1)      
      real*8,intent(out) :: gs  !Leaf stomatal conductance (mol-H2O m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (CO2 uptake, micromol m-2 s-1)
      real*8,intent(out) :: Atot0 !Leaf gross photosynthesis (CO2 uptake, micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Dark = above-ground growth + maintenance respiration (umol m-2 s-1)
      real*8,intent(out) :: isp ! Isoprene emission (umol C m-2 s-1)
      real*8,intent(out) :: mtp ! Monoterpene emission (umol C m-2 s-1)
        !---Local----
!      type(photosynthpar) :: pspar !Moved to global to module.
      real*8,parameter :: O2pres=20900.d0 !O2 partial pressure in leaf (Pa) Not exactly .209*101325.
      real*8 :: cie, cic, cis   !Leaf internal CO2 (umol mol-1)
!      real*8 :: Je1, Jc1, Js1   !Assimilation of CO2, 3 limiting cases
      real*8 :: Anet            !Net assimilation of CO2 = Atot - aboveground respir (umol m-2 s-1)
      real*8 :: Aiso            ! Rate of photosynthesis for isoprene emissions (umol m-2 s-1)
      real*8 :: cs   !CO2 mole fraction at the leaf surface (umol mol-1)
      real*8 :: Ae, Ac, As      !* These are Anet!
      real*8, save :: a1c=1.d30, f1c=-1.d30
      real*8 :: a1e, f1e
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8 :: F        ! O3-induced reduction fraction
      real*8 :: Ra       ! boundary resistance (s m-1)
      real*8 :: gp       ! leaf conductance (m s-1)
#endif

      real*8 :: alpha !Intrinsic quantum efficiency for CO2 uptake
#ifdef PS_BVOC
      logical, parameter :: need_isoprene = .true.
#else
      logical, parameter :: need_isoprene = .false.
#endif
      integer, save :: counter = 0
      counter = counter + 1

      !write(888,*) "counter=", counter

      Rd = Rdark(pspar%Vcmax)

      alpha = 0.08d0
      If (pft.eq.GRASSC4) alpha = 0.04d0

      if ( IPAR < .000001d0 ) then
        Atot = 0.d0
        Atot0= 0.d0
        Anet = - Rd
        cs = ca - Anet*1.37d0/gb
        gs = BallBerry(Anet, rh, cs, pspar)
        ci = cs - Anet/(gs/1.65d0) 
        isp = 0.d0
        mtp = 0.d0
        return
      endif


      !* Photosynthetic rate limited by RuBP saturation
      !* Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
!      call Ci_Jc(ca,gb,rh,IPAR,Pa,pspar, Rd,O2pres, cic, Jc1)
      ! Jc_RuBP = pspar%Vcmax*(Cip - pspar%Gammastar)/
      !           (Cip + pspar%Kc*(1 + O2/pspar%Ko))

      !Assimilation is of the form a1*(Ci - Gammastar)/(e1*Ci + f)
      if ( pspar%first_call ) then
         if (pfpar(pspar%pft)%pst.eq.C3) then
            a1c = pspar%Vcmax
            f1c = pspar%Kc*(1.d0 + O2pres/pspar%Ko) * 1.d06/Pa !umol/mol

            !call ci_cubic (ca,rh,gb,Pa,Rd,a1c,f1c,pspar,Axxx)
            call ci_cubic(ca,rh,gb,Pa,Rd,a1c,f1c,pspar,Ac)
            !if ( Ac >= -Rd ) write(578,*) Axxx, Ac, Ac - Axxx
            !write(888,*) "Ac", ca,rh,gb,Pa,Rd,a1,f1,pspar,Ac
         else !C4 photosynthesis
            Ac = pspar%Vcmax - Rd
         endif
         pspar%Ac = Ac
         !pspar%first_call = .false. !Probably bug-prone, but reset after As
      else 
         Ac = pspar%Ac
      endif


!      call Ci_Je(ca,gb,rh,IPAR,Pa, pspar, Rd, cie, Je1)
      ! Photosynthetic rate limited by light electron transport (umol m-2 s-1)
      ! Je_light = (pspar%PARabsorb*IPAR)*alpha*(Cip-pspar%Gammastar)/
      !            (Cip+2*pspar%Gammastar)

      !Assimilation is of the form a1*(ci - Gammastar.umol)/(e1*ci + f1)

      if (pfpar(pspar%pft)%pst.eq.C3) then
         !a1 = pspar%PARabsorb*IPAR*alpha
         a1e = IPAR*alpha       !### HACK:  IPAR from canopyspitters.f is APAR.  When we switch to Wenze's canopyrad, then leaf PARabsorb will be used -NK ###
         f1e = 2*pspar%Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

         if ( a1e < a1c .or. 
     &        f1e > f1c .or.
     &        need_isoprene ) then
            !call ci_cubic (ca,rh,gb,Pa,Rd,a1e,f1e,pspar,Axxx)
            call ci_cubic(ca,rh,gb,Pa,Rd,a1e,f1e,pspar,Ae)
            !write(888,*) "Ae", ca,rh,gb,Pa,Rd,a1,f1,pspar,Ae 
cddd        call ci_cubic1(ca,rh,gb,Pa,Rd,a1,f1,pspar,Axxx)
cddd        write(579,*) Ae, Axxx
cddd        if ( Ae > 0.d0 ) write(578,*) Axxx - Ae
        !if ( Ae >= -Rd ) write(578,*) Axxx, Ae, Ae - Axxx
         else
            Ae = 1.d30
         endif
      else !C4 photosynthesis
         Ae = IPAR*alpha - Rd
      endif


      !* Photosynthetic rate limited by utilization of photosynthetic products:
      !* (umol m-2 s-1)  Triosphosphate (TPU limitation for C3,
      !*                 PEP carboxylase limitation for C4.
      if (pfpar(pspar%pft)%pst.eq.C3) then
           !call Ci_Js(ca,gb,rh,IPAR,Pa,pspar,Rd, cis, Js1)
           !Js_sucrose = pspar%Vcmax/2.d0
            As = pspar%Vcmax/2.d0 - Rd  !Anet
            pspar%As = As
           !write(888,*) "As", As
      else  !C4 photosynthesis
            !As = 4000.d0*pspar%Vcmax*ci - Rd
         if (pspar%first_call) then
            call Asnet_C4(ca,rh,gb,Rd,pspar,As) !This is Anet
            pspar%As = As
         else
            As = pspar%As
            pspar%first_call = .false.
         endif
      endif

      Anet = min(Ae, Ac, As)

#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      cs = ca - Anet*1.37d0/gb
      gs = BallBerry(Anet, rh, cs, pspar)
      gp = gs*mol_to_ml
      Ra = 1.d0/(gb*mol_to_ma)
      call  o3dep_uptake(Ra, gp, pspar%pft, o3, .true., F, Fo3, dFo3)
      if (F .gt. 1.0d0 .or. F .lt. 0)
     & call stop_model('F>1 or F<0 in o3dep_damage',255)
#endif

#ifdef O3DEP_UPTAKE
      Anet = Anet*F
      Rd   = Rd*F
#endif

      Atot = Anet + Rd
      Aiso = Ae + Rd

      if (Atot.lt.0.d0) then
      ! can only happen if ca < Gammastar . Does it make sense? -Yes-NK
#ifdef OFFLINE
         write(997,*) "Error, Atot<0.0:",Atot,Ae,Ac,As,ca,gb,rh,IPAR
     &        ,Pa,pspar,sunlitshaded, pspar%Gammastar * 1.d06/Pa
#endif
         Atot = 0.d0
         Anet = - Rd
!!       ci = pspar%Gammastar * 1.d06/Pa  
!!       gs = 0. ! MK: setting to 0 to avoid erratic results

!!    else
      endif

      cs = ca - Anet*1.37d0/gb
      gs = BallBerry(Anet, rh, cs, pspar)
      ci = cs - Anet/(gs/1.65d0)

      Atot0 = Atot
#ifdef O3DEP_UPTAKE_OFFLINE
      Atot0 = Atot*F
#endif


!!!        endif

#ifdef PS_BVOC
         call Voccalc(pft,pa,ca,ci,Tl,pspar%Gammastar,
     & isp,mtp,Aiso)

#else
       isp=0.0d0
       mtp=0.0d0
#endif

      !write(888,*) "gs,ci,cs", gs,ci,cs

      end subroutine Photosynth_analyticsoln

!-----------------------------------------------------------------------------
      subroutine Voccalc(pft,pa,ca,ci,Tl,Gammastar,isp,mtp,Aiso)
!@sum Voccalc  Isoprene emissions coupled to photosynthesis
!@auth Nadine Unger
      use yibs_pfts

      implicit none
      integer,intent(in) :: pft !Plant functional type, 1-C3 grassland
      real*8,intent(in) :: ca   !Ambient air CO2 mole fraction (umol mol-1)      
      real*8,intent(in) :: Pa   !Pressure (Pa)
      real*8,intent(in) :: Tl   !Leaf temperature (Celsius)
      real*8,intent(in) :: Gammastar   
      real*8,intent(in) :: Aiso   !(umol m-2 s-1)
      real*8,intent(in) :: ci
      real*8,intent(out) :: isp ! isoprene emission (umol C m-2 s-1)
      real*8,intent(out) :: mtp ! monoterpene emission (umol C m-2 s-1)
      type(photosynthpar) :: pspar
        !---Local----
      integer, parameter :: numpft = 8
      real*8 :: gammamol,fact
      real*8 :: IBASER, Y_alpha, kapco2
      real*8 :: tauiso,taumtp
C April 2009 values
c      real*8, parameter, dimension(numpft) :: Y_eps = !
c     & (/0.0d0,2.10d-02,8.24d-02,6.48d-02,1.08d-01,
c     & 4.44d-02,1.38d-01,0.0d0/)
C New values June 2009
C      real*8, parameter, dimension(numpft) :: Y_eps = !
C     & (/0.0d0,1.91d-02,7.19d-02,5.13d-02,8.79d-02,
C     & 2.18d-02,8.35d-02,0.0d0/)

      gammamol =  Gammastar * 1.d06/Pa !Convert from Pa to umol/mol

C Y_alpha, Y_eps unitless
  
      Y_alpha=(ci-gammamol)/(6.0*((4.67*ci)+(9.33*gammamol)))

      isp = pfpar(pft)%Y_eps*Aiso*Y_alpha

      mtp = pfpar(pft)%Y_eps_m

C Include CO2 effects

      kapco2 = (0.7*370.0)/ci

C Include temperature effects

      tauiso = exp(0.1*(Tl-30.0))
      taumtp = exp(0.09*(Tl-30.0))

C Include seasonal effects? Add later.
C Note can switch on and off kapco2

      isp = isp*kapco2*tauiso
      mtp = mtp*taumtp*kapco2

      end subroutine Voccalc

#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)

      subroutine o3dep_uptake(Ra, gp, ipft, o3conc, high, F, Fo3, dFo3)

      use yibs_const,  only: io3d, fo3c, ah, al, ko3

      implicit none

      integer, intent(in)  :: ipft           ! PFT types
      real*8,  intent(in)  :: Ra             ! aerodynamic resist. between leaf and ref. level (s m-1)
      real*8,  intent(in)  :: gp             ! leaf conductance of water vapor (m/s)
      real*8,  intent(in)  :: o3conc         ! [O3] concentrations (nmol m-3)
      logical, intent(in)  :: high           ! ture: high sens. false (default): low sens.
      real*8,  intent(out) :: F              ! fraction of GPP by o3 to the orignial GPP
      real*8,  intent(out) :: Fo3            ! ozone flux to stomata (nmol O3 m-2 s-1)
      real*8,  intent(out) :: dFo3           ! excessive ozone flux to stomata (nmol O3 m-2 s-1)
      integer iflag, id
      real*8  aa, bb, cc, f1, f2
      real*8  a0, b2, bx, fc

      F = 1.0d0
      iflag = 0
      id  = io3d(ipft)
      if (o3conc .lt. 1.0d-8 .or. gp .lt. 1.0d-8
     &   .or. Ra .lt. 1.0d-8 .or. id .lt. 0) return

      a0 = al(id)                            ! for a low sensitivity
      if (high) a0 = ah(id)                  ! for a high sensitivity
      aa = gp*Ra
      bb = ko3 - gp*Ra*(1.0d0+a0*fo3c(id)) + a0*o3conc*gp
      cc = -ko3*(1.0d0+a0*fo3c(id))
      b2 = bb*bb - 4.0d0*aa*cc
      bx = sqrt(b2)
      f1 = (-bb+bx)/(2.0d0*aa)
      f2 = (-bb-bx)/(2.0d0*aa)

      fc = max(ko3*fo3c(id)/(o3conc-fo3c(id)*Ra)/gp, 0.0d0)

      if (f1 .gt. fc .and. f1 .le. 1.0d0) then
         F = f1
         iflag = 1
      endif
      if (f2 .gt. fc .and. f2 .le. 1.0d0) then
         F = f2
         if (iflag .gt. 0)
     &      call stop_model('double F in o3dep_damage',255)
      endif

      Fo3 = o3conc/(Ra+ko3/(gp*F))
      dFO3 = Max(Fo3-fo3c(id), 0.0d0)

      return
      end subroutine o3dep_uptake

#endif


      function calc_CO2compp(O2,Kc,Ko,Tl) Result(Gammastar)
!@sum calc_CO2compp CO2 compensation point in absence of dark respiration (Pa)

      implicit none
      real*8,intent(in) :: O2 !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: Kc   !Michaelis-Menten parameter for CO2 (Pa)
      real*8,intent(in) :: Ko   !Michaelis-Menten parameter for O2 (Pa)
      real*8,intent(in) :: Tl !Leaf temperature (Celsius)
      real*8 :: Gammastar  !CO2 compensation point (Pa)
      !----Local-----
      real*8,parameter :: tau=2600.d0  !CO2/O2-specificity ratio

!      Gammastar = O2/(2.d0*tau*Q10fn(0.57d0,Tl)) !Collatz (A3)
!      Gammastar = O2*Q10fn(1.75,Tl)/(2.d0*tau) !Collatz (A3) Same as above, !KcQ10/KoQ10 = 2.1/1.2 = 1.75 = 1/.57 
!      Gammastar = 0.5d0*(Kc/Ko)*0.21*O2 !CLM Tech Note. Gives smaller Gammastar than Collatz.

C Nadine - use previous T-dep version
#ifdef PS_BVOC
       Gammastar = O2/(2.d0*tau*Q10fn(0.57d0,Tl)) !Collatz (A3) 
#else
       Gammastar = 0.5d0*(Kc/Ko)*0.21*O2 !CLM Tech Note. Gives smaller Gammastar than Collatz.
#endif

      end function calc_CO2compp
!-----------------------------------------------------------------------------

      
      function arrhenius(Tcelsius,c1,c2) Result(arrh)
!@sum arrhenius Arrhenius response to temperature for biological kinetics.
!@+   From David Medvigy's lphys.f90
      implicit none
      real*8 :: Tcelsius
      real*8 :: c1,c2 !Process-specific temperature response parameters.
      real*8 :: arrh

      arrh = c1*exp(c2*(1.d0/288.15d0-1.d0/(Tcelsius+Kelvin)))
      return
      end function arrhenius
!=================================================
      function Q10fn(Q10par,Tcelsius) Result(Q10factor)
!@sum Q10fn  Q10 function, biological response to temperature.
!@+   From Collatz, et al. (1991)
      implicit none
      real*8 :: Q10par, Tcelsius !parameter, temperature 
      real*8 :: Q10factor

      Q10factor = Q10par**((Tcelsius-25.d0)/10.d0)

      end function Q10fn
!=================================================

      function Tresponse(c,deltaH,Tcelsius) Result(Tfactor)
!@sum Tresponse  Arrhenius temperature response function that accounts for
!@+   activation energy.
!@+   From Bernacchi, et al. (2001).
      implicit none
      real*8,intent(in) :: c !Scaling factor
      real*8,intent(in) :: deltaH !Activation energy 
      real*8,intent(in) :: Tcelsius !Temperature (Celsius)
      real*8 :: Tfactor

      Tfactor = exp(c - deltaH/(Rgas * (Tcelsius + Kelvin)))

      end function Tresponse
!=================================================
      subroutine  Asnet_C4(ca,rh,gb,Rd,pspar,Asnet)
!@sum Asnet_C4 PEP carboxlase-limited carbon assimilation for C4 photosynthesis
!@+   After Collatz, and CLM's correction of the coefficient
!@+   Returns Asnet = Astot - Rd
!@+   Solving for Asnet via the equations:
!@+   1) Asnet = Astot - Rd
!@+           = (ca - ci)/[(1.37*rb + 1.65*rs)] 
!@+           = (ca - cs)/(1.37*rb) 
!@+           = (cs - ci)/(1.65*rs)
!@+   2) 1/rs = gs = m*A*rh/cs + b     
!@+   3) Astot = 4000.d0*pspar%Vcmax*(ci*1e-06)  
!@+        !4000 is CLM, Collatz had 1800. Convert ci from umol/mol to mol/mol
!@+   Do subsitutions to eliminate cs and ci and solve for Asnet.


      real*8,intent(in) :: ca  !Surface air CO2 concentration (umol/mol)
      real*8,intent(in) :: rh   !Relative humidity
      real*8,intent(in) :: gb   !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: Rd   !Leaf mitochondrial respiration (umol m-2 s-1)
      type(photosynthpar) :: pspar
      real*8,intent(out) :: Asnet  !Net assimilation of carbon (umol m-2 s-1)
      !---Local----
      real*8 :: ci !Leaf internal CO2 concentration (umol/mol)
      real*8 :: K1, K2, K3, K4
      real*8 :: X, Y, Z
      real*8 :: b0, a0, c0
      real*8 :: Aspos, Asneg

      K1 = pspar%m * rh
      K2 = 4000.d0*pspar%Vcmax * 1.d-06
      K3 = 1.37d0/gb
      K4 = 1.65d0*ca
      
      !Anet^2*X + A*Y + Z = 0
      X = (K1 - pspar%b*K3)*(1/K2 + K3) - 1.65d0*K3
C      Y = ca*(pspar%b/K2 - K1 + 1.65d0) + Rd/K2*(K1 - pspar%b*K3)
C      Z = -pspar%b * ca * ( Rd/K2 - ca )
Cxyue
      Y = ca*(pspar%b/K2 - K1 + 2*pspar%b*K3 + 1.65d0)
     &  + Rd/K2*(K1 - pspar%b*K3)
      Z = pspar%b * ca * ( Rd/K2 - ca )
Cxyue


      !This section is correct to solve for Atot, solves to Anet+Rd.
      !Easier just to add Rd to Asnet.
      !a0*Atot^2 + b0*Atot + c0 = 0
      !a0 = X
      !b0 = -2.d0*Rd*X + Y
      !c0 = (Rd**2.d0)*X - Rd*Y + Z
      !Aspos = (-b0 + sqrt(b0**2.d0 - 4.d0*a0*c0))/(2*a0)
      !Asneg = (-b0 - sqrt(b0**2.d0 - 4.d0*a0*c0))/(2*a0)
      !Astot = max(Aspos, Asneg)

      Asnet = max (-Rd
     &     ,(-Y + sqrt(Y**2.d0 - 4.d0*X*Z))/(2.d0*X)) !Positive root is max.

      end subroutine Asnet_C4

!=================================================

#ifndef USE_NR_SOLVER_FOR_FBB
      subroutine ci_cubic(ca,rh,gb,Pa,Rd,a1,f1,pspar,A)
!@sum ci_cubic Analytical solution for cubic equation of coupled
!@+   Ball-Berry/Farquhar stomatal conductance/photosynthesis.
!@+   Version that uses analytical equation solution.
!@+   Solves for Anet.
!@auth I.Aleinov
      !@sum ci (umol/mol)
      !@sum For the case of assimilation being of the form:
      !@sum         A = Anet = a*(Cip - Gammastar)/(e*Cip + f) - Rd
      !@sum Numerator and denominator are converted from (Pa/Pa) to (umol mol-1)/(umol mol-1)
      !@sum         A = Anet = a1*(ci - gammamol) /(e1*ci + fmol) - Rd
      !@sum where gammamol = Gammastar*1d06/Pa, fmol = f1 = f*1d06/Pa

      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: rh              !Relative humidity
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      real*8 :: Rd              !Leaf mitochondrial respiration (umol m-2 s-1)
      real*8 :: a1              !Coefficient in linear Farquhar equ.
      real*8 :: f1              !Coefficient in linear Farquhar equ.
      real*8, intent(out) :: A
      type(photosynthpar) :: pspar
      !----Local----
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 :: Ra, b, K, gamol, A_d_asymp
      real*8 :: X, Y, Z, Y1  ! tmp vars
      real*8 :: c3, c2, c1, c   !Coefficients of the cubic of ci (c3*ci^3 + c2*ci^2 + c1*ci + c)
      real*8 :: cixx(3) ! solutions of cubic
!      real*8 :: cs, Rs ! needed to compute ci
      integer :: nroots, i
!      real*8 ci

      Ra = 1/gb * S_ATM
      b = pspar%b / S_STOM
      K = pspar%m * rh / S_STOM
      gamol = pspar%Gammastar * 1.d06/Pa !Convert Pa to umol/mol
      A_d_asymp = - b*Ca / (K - b*Ra) ! asymptotic val of A from diffusion eq.

      ! first check some special cases
      if ( A_d_asymp >= 0.d0 ) then
        ! this can happen only for very low humidity
        ! probably should never happen in the real world, but if it does,
        ! this case should be considered separately
        !!print *,"!!! A_d_asymp >= 0.d0 !!!", A_d_asymp
        A_d_asymp = -1.d30 !!! hack
        !!print *,"K<b*Ra: m,rh,b,Ra:",pspar%m,rh,b,Ra
        !!call stop_model("ci_cubic: rh too small ?",255)
      endif

      ! dependence on e1 if needed
cddd      Y= f1/e1
cddd      X= -a1/e1 * (gamol+f1/e1)
cddd      Z= a1/e1 -Rd
      Y= f1
      X= -a1 * (gamol+f1)
      Z= a1 -Rd

      if ( Z > 0.d0 ) then
        ! Farquhar curve is above zero. May have solution A > 0
        c = -(b*Ca*(X + (Ca + Y)*Z))
        c1 = Ca*Z - K*(X + Ca*Z + Y*Z) + 
     &       b*(Ca**2 + Ca*(Y + 2*Ra*Z) + Ra*(X + Y*Z))
        c2 = Ca*(-1 + K - 2*b*Ra) + K*(Y + Ra*Z) - Ra*(b*Y + Z + b*Ra*Z)
        c3 = Ra*(1 - K + b*Ra)

        call cubicroot(c3, c2, c1, c, cixx, nroots)

        !!print *,"roots= ", cixx(1:nroots)

        ! find minimal root above the asymptotic value
        A = 1.d30
        do i=1,nroots
          if ( cixx(i) < A .and. cixx(i) > A_d_asymp ) A = cixx(i)
        enddo
        if ( A == 1.d30 )  then
          print *," m,rh,b,Ra:",pspar%m,rh,b,Ra
          print *,"ca,gb,Pa:",ca,gb,Pa
          print *,"pspar:",pspar
          print *," A_d_asymp,K,gamol,f1,a1,Rd",
     &         A_d_asymp,K,gamol,f1,a1,Rd
          print *,"c3,c2,c1,c", c3,c2,c1,c
          print *,"nroots,cixx",nroots,cixx(1:nroots)
          call stop_model("ci_cubic: no solution",255)
        endif

        if ( A >= 0 ) then
cddd          cs = ca - A*Ra
cddd          Rs = 1.d0 / (K*A/cs + b)
cddd          ci = cs - A*Rs
cddd          ! just in case, check consistency
cddd          if ( ci < 0.d0 ) call stop_model("ci_cubic: ci<0",255)
cddd          if ( cs < 0.d0 ) call stop_model("ci_cubic: cs<0",255)
cddd          !!print *,'QQQQ ',A,ci
          return
        endif

      endif

      ! if we got here then A<0 : have to solve quaratic equation

      Y1 = Y + ca
      c2 = Ra + 1.d0/b
      c1 = - (Y1 + c2*Z)
      c  = X + Y1*Z

      ! just in case,
      if (  c1*c1 - 4.d0*c2*c < 0.d0 )
     &     call stop_model("ci_cubic: no solution to quadratic",255)
      A = ( - c1 - sqrt( c1*c1 - 4.d0*c2*c ) ) / ( 2.d0 * c2 )
cddd      cs = ca - A*Ra
cddd      Rs = 1.d0 / ( b)
cddd      ci = cs - A*Rs
cddd      !!print *,"q ", ci, cs, A, Rs, Ra
cddd      ! just in case, check consistency
cddd      if ( ci < 0.d0 ) call stop_model("ci_cubic: q: ci<0",255)
cddd      if ( cs < 0.d0 ) call stop_model("ci_cubic: q: cs<0",255)
cddd      !!print *,'QQQQ ',A,ci

      end subroutine ci_cubic


!=================================================
      subroutine cubicroot(a,b,c,d,x,n) 
!@sum cubicroot  Solve cubic equation: a x^3 + b x^2 + c x + d = 0 *!
!@auth I.Aleinov
      !* Written by Igor Aleinov from solution by Cardano in
      !* Korn, Korn, Mathematical Handbook.
      implicit none
      real*8,intent(in) :: a,b,c,d  ! coefficients of cubic
      real*8, intent(out) :: x(:)   ! results ( 0-3 roots )
      integer, intent(out) :: n     ! number of roots
      real*8 :: x0,x1,x2
      real*8 :: a0,a1,a2,Q1,R1,D1
      real*8, parameter :: EPS0 = 1.d-8 ! 1.d-15
      real*8, parameter :: one3rd = 1.d0/3.d0
      real*8 :: arg, S, T
      complex*16 :: ST

      !print *,"cubicroot:",a,b,c,d

      if (abs(a) < (abs(b)+abs(c)+abs(d))*EPS0 ) then
        if (abs(b) < (abs(c)+abs(d))*EPS0) then
          if (abs(c) < abs(d)*EPS0) then
            write(*,*) "Internal Error in Cardano: no solution."
            stop
          endif
          x0 = -d/c
          x(1) = x0
          !write(*,*) "Cardano: returning",x0
          n = 1
        else
          !write(*,*) "What's this?"
          D1 = c*c - 4.d0*b*d
          
          if (D1 > 0.d0) then
            Q1 = sqrt(D1)
            x0 = (-c + Q1) / (2.d0 * b)
            x1 = (-c - Q1) / (2.d0 * b)
            !return
            n = 2
          else if (D1.eq.0.) then
            x0 = -c / (2.d0 * b)
            x1 = x0
            n = 1
          else 
            x0 = -c /(2.d0 *b)
            x1 = sqrt(-D1) / (2.d0* b)
            n = 0
          end if
        end if
        !print *,"CX1",x0,x1
        !x = max(x0,x1)
        x(1) = x0
        x(2) = x1
      else
        a2 = b/a
        a1 = c/a
        a0 = d/a
        Q1 = (3.d0 * a1 - a2*a2 ) / 9.d0
        R1 = (9.d0 * a2 * a1 - 27.d0 * a0 - 2.d0 * a2*a2*a2) /54.d0
        D1 = Q1*Q1*Q1 + R1*R1
        !write(*,*) "abcda2a1a0Q1R1D1",a,b,c,d,a2,a1,a0,Q1,R1,D1
        if (D1 > 0.d0) then       !* only one real root *!
          !write(*,*) "One real root."
          arg = R1 + sqrt(D1)
          S = sign(1.d0, arg) * (abs(arg)**one3rd)
          arg = R1 - sqrt(D1)
          T = sign(1.d0, arg) * (abs(arg)**one3rd)
          x0 = -a2/3.d0 + S + T
          x1 = -a2/3.d0 - (S+T)*0.5d0
          x2 = sqrt(3.d0) * (S-T)*0.5d0
          !print *,"CX2",x0,x1,x2
          n = 1
        else if (D1.eq.0.) then !* two roots coincide * *!
          !write(*,*) "Two roots coincide."
          S = sign(1.d0, R1) * (abs(R1)**one3rd)
          x0 = -a2/3.d0 + 2.d0*S
          x1 = -a2/3.d0 - S
          x2 = x1
          !print *,"CX3",x0,x1,x2
          n =2
        else                    !* three different real roots *!
          !call CRtCube( R1, sqrt(-D1), S, T)
          !write(*,*) "Three different real roots. a2R1D1ST",a2,R1,D1,S,T
          ST = ( cmplx(R1, sqrt(-D1),kind(1.d0)) )**one3rd
          S = real (ST)
          T = aimag(ST)
          x0 = -a2/3.d0 + 2.d0*S
          x1 = -a2/3.d0 - S + sqrt(3.d0)*T
          x2 = -a2/3.d0 - S - sqrt(3.d0)*T
          !print *,"CX4",x0,x1,x2
          n = 3
        end if
        !x = max(x0,x1,x2)
        !x = x2
        x(1) = x0
        x(2) = x1
        x(3) = x2
      end if
      end subroutine cubicroot
#endif

!=================================================

      subroutine calc_Pspar(dtsec,pft,Pa,Tl,O2pres,stressH2O,
     &                      Sacclim,llspan)
!@sum calc_Pspar Collatz photosynthesis parameters in data structure pspar.
!@    pspar is GLOBAL TO MODULE.
      !Later need to replace these with von Caemmerer book Arrhenius
      !function sensitivities (her Table 2.3)
      implicit none
      integer,intent(in) :: pft   !Plant functional type, 1=C3 grassland
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: Pa     !Atmospheric pressure (Pa)
      real*8,intent(in) :: Tl     !Leaf temperature (Celsius)
      real*8,intent(in) :: O2pres !O2 partial pressure in leaf (Pa)
      real*8,intent(in) :: stressH2O
      real*8,intent(in) :: Sacclim !state of acclimation/frost hardiness
      real*8,intent(in) :: llspan !mean leaf life span
!      type(photosynthpar),intent(inout) :: pspar !Moved to global to module.
      integer :: p

      !----Local-----
      real*8 :: facclim ! acclimation/forst hardiness factor [-]
      !Below parameters are declared at top of module, though only used here.
!      real*8,parameter :: Kc        !Michaelis-Menten constant for CO2 (Pa)
!      real*8,parameter :: Ko        !Michaelis-Menten constant for O2 (Pa)
!      real*8,parameter :: KcQ10           !Kc Q10 exponent
!      real*8,parameter :: KoQ10           !Ko Q10 exponent
      real*8 :: fparlimit !light(i.e.,PAR) control
      integer, save :: counter = 0
      counter = counter + 1

      facclim = frost_hardiness(Sacclim)
     
      fparlimit = par_phenology(pft,llspan)
!      fparlimit = 1.d0

!!! this var is not reproducible on restart, please figure out why
!      fparlimit = 1.d0 ! seems to be ok now

      !write(877,*) "counter", counter
      !write(877,*) "facclim", facclim
      !write(877,*) "fparlimit", fparlimit

      p = pft
      pspar%pft = pft
      pspar%PARabsorb = pftpar(p)%PARabsorb !Collatz et al. (1991)
!      pspar%Vcmax = pftpar(p)%Vcmax/(1 + exp((-220.e03+703.*(Tl+Kelvin))
!     &     /(Rgas*(Tl+Kelvin))))
      pspar%Vcmax = pftpar(p)%Vcmax * Q10fn(2.21d0, Tl)
     &            * facclim
      pspar%RVcmax = pftpar(p)%Vcmax * Q10fn(3.09d0-0.043d0*Tl, Tl)
     &            * facclim
cxyue   &            * facclim * fparlimit
      pspar%Kc = Kc*Q10fn(KcQ10,Tl) !(Collatz, eq. A12)
      pspar%Ko = Ko*Q10fn(KoQ10,Tl) !(Collatz, eq. A12)
      pspar%Gammastar = calc_CO2compp(O2pres,pspar%Kc,pspar%Ko,Tl) !(Pa) (Collatz)
      pspar%m = stressH2O*pftpar(p)%m     !Slope of Ball-Berry equation (Collatz)
!      pspar%m = pftpar(p)%m     !Slope of Ball-Berry equation (Collatz)
      pspar%b = pftpar(p)%b     !Intercept of Ball-Berry equation (mol m-2 s-1) (Collatz)
      pspar%Nleaf = pftpar(p)%Nleaf !g-N/m^2[leaf] Needed for foliar respiration.
      !pspar%Nleaf = pftpar(p)%Nleaf * phenology factor !Here can adjust Nleaf according
                                !to foliage N pools or phenology factor.
      pspar%stressH2O = stressH2O

      pspar%first_call = .true.
!      pspar%As = 0.d0 !Unnecessary but zero anyway
!      pspar%Ac = 0.d0 !Unnecessary but zero anyway
      pspar%reset_ci_cubic1 = .true.

      end subroutine calc_Pspar

!-----------------------------------------------------------------------------
      real*8 function par_phenology(pft,llspan) Result(fparlimit)  
!@sum par_phenology  Physiological status variable for Vcmax of tropical
!@+   broadleaf evergreen trees.
!@auth Y.Kim
      integer, intent(in) :: pft
      real*8, intent(in) :: llspan
      real*8, parameter :: vc_tran =12.d0 !9.d0 ! 7.2     !transition
      real*8, parameter :: vc_slop = 15.d0 !10.d0 !16.9    !slope
      real*8, parameter :: vc_amp = 15.d0 !30.d0 !29.8    !amplitude
      real*8, parameter :: vc_min = 10.d0 !25.d0 !7.7     !minimum

      if (llspan > 0.d0) then
         fparlimit = (vc_amp/(1.d0+(llspan/vc_tran)**vc_slop)+vc_min)
     &               /pftpar(pft)%Vcmax
      else
         fparlimit = 1.d0
      endif

      end function par_phenology        
!-----------------------------------------------------------------------------

      real*8 function frost_hardiness(Sacclim) Result(facclim)
!@sum frost_hardiness.  Calculate factor for adjusting photosynthetic capacity
!@+   due to frost hardiness phenology.
!@+   Based on Repo et al (1990), Hanninen & Kramer (2007),
!@+   and Makela et al (2006)
!@auth M.Puma
      real*8,intent(in) :: Sacclim 
!      real*8 :: facclim ! acclimation/frost hardiness factor [0. to 1.]
      !----Local-----
      real*8,parameter :: Tacclim=-5.93d0 ! threshold temperature for photosynthesis [deg C]
      !real*8,parameter :: Tacclim=-3.d0 ! Best tune for Hyytiala
                        ! Site specific thres. temp.: state of photosyn.acclim
                        ! Hyytiala Scots Pine, -5.93 deg C Makela et al (2006)
      !real*8,parameter :: a_const=0.0595 ! factor to convert from Sacclim [degC] to facclim [-]
                        ! Site specific; conversion (1/Sacclim_max)=1/16.8115
                        ! estimated by using the max S from Hyytiala 1998
      real*8, parameter :: a_const = 0.1d0 !Closer tune for Hyytiala

      if (Sacclim > Tacclim) then ! photosynthesis occurs 
         facclim = a_const * (Sacclim-Tacclim) 
         if (facclim > 1.d0) facclim = 1.d0
!      elseif (Sacclim < -1E10)then !UNDEFINED
      elseif (Sacclim.eq.UNDEF)then !UNDEFINED
         facclim = 1.d0         ! no acclimation for this pft and/or simualtion
      else
         facclim = 0.01d0       ! arbitrary min value so that photosyn /= zero
      endif

      end function frost_hardiness

!-----------------------------------------------------------------------------

      function BallBerry(Anet, rh, cs, pspar) Result (gsw)
!@sum Ball-Berry (1987) model of leaf stomatal conductance of 
!@    water vapor, gsw (mol m-2 s-1)      
!@auth N.Y.Kiang
      implicit none
      real*8,intent(in) :: Anet !Net assimilation of CO2 (umol m-2 s-1)
      real*8,intent(in) :: rh   !Relative humidity (fractional ratio)
      real*8,intent(in) :: cs   !Leaf surface CO2 mole fraction (umol mol-1)
      type(photosynthpar) :: pspar
      real*8 :: gsw !Leaf conductance of water vapor (mol m-2 s-1)
      !----Local-----
      
      ! just in case check cs (remove after debugging ?)
      if ( cs <= 0.d0 ) call stop_model("BallBerry: cs <= 0", 255)
      gsw = pspar%m*Anet*rh/cs + pspar%b
      if (gsw < pspar%b) gsw = pspar%b

      end function BallBerry

!-----------------------------------------------------------------------------

#ifdef USE_NR_SOLVER_FOR_FBB
      subroutine ci_cubic(ca,rh,gb,Pa,Rd,a1,f1,pspar,A)
!@sum ci_cubic Numerical solution for cubic equation of coupled
!@+   Ball-Berry/Farquhar stomatal conductance/photosynthesis.
!@+   Solves for Atot.
!@+   Version that uses Newton-Raphson solver
!@auth I.Aleinov
      implicit none
      real*8 :: ca              !Ambient air CO2 concentration (umol mol-1)
      real*8 :: rh              !Relative humidity
      real*8 :: gb              !Leaf boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: Pa              !Pressure (Pa)
      real*8 :: Rd              !Leaf mitochondrial respiration (umol m-2 s-1)
      real*8 :: a1              !Coefficient in linear Farquhar equ.
      real*8 :: f1              !Coefficient in linear Farquhar equ.
      real*8, intent(out) :: A
      type(photosynthpar) :: pspar
      !----Local----
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 :: Ra, b, K, gamol, A_d_asymp
      real*8 :: x1, x2, xacc, x2tmp, x2save
      integer :: numit, counter=0
      save

      if ( pspar%reset_ci_cubic1 ) then
        pspar%reset_ci_cubic1 = .false.
cddd      if ( pspar%reset_ci_cubic1 == .false. ) then
cddd        write(777,*) counter, Ra, b, K, gamol, A_d_asymp,
cddd     &       x1, x2, xacc, x2tmp
cddd      endif
      Ra = 1/gb * S_ATM
      b = pspar%b / S_STOM
      K = pspar%m * rh / S_STOM
      gamol = pspar%Gammastar * 1.d06/Pa !Convert Pa to umol/mol
      A_d_asymp = - b*Ca / (K - b*Ra) ! asymptotic val of A from diffusion eq.

cddd      ! first check some special cases
cddd      if ( A_d_asymp >= 0.d0 ) then
cddd        ! this can happen only for very low humidity
cddd        ! probably should never happen in the real world, but if it does,
cddd        ! this case should be considered separately
cddd        !!print *,"!!! A_d_asymp >= 0.d0 !!!", A_d_asymp
cddd      !!!  A_d_asymp = -1.d30 !!! hack
cddd        !!print *,"K<b*Ra: m,rh,b,Ra:",pspar%m,rh,b,Ra
cddd        !!call stop_model("ci_cubic: rh too small ?",255)
cddd      endif

      !x1 = 0.d0
      x1 = -Rd
      x2save = ca/Ra
      x2tmp =  b*ca / (1.d0 - K + b*Ra)
      if( x2tmp > 0.d0 ) x2save = min( x2save, x2tmp )
      x2tmp = A_d_asymp
      if( x2tmp > 0.d0 ) x2save = min( x2save, x2tmp )
      x2save = x2save - .0000001d0
      x2 = min( x2, a1 - Rd)
      xacc = .0001d0
      !xacc = .01d0
cddd      if ( pspar%reset_ci_cubic1 == .false. ) then
cddd        write(778,*) counter, Ra, b, K, gamol, A_d_asymp,
cddd     &       x1, x2, xacc, x2tmp
cddd      endif
      endif
      x2 = min( x2save, a1 - Rd)
      A = rtsafe(A_eqn, x1,x2,xacc,  Ra, b, K, gamol,  ca, a1, f1, Rd
     &     , numit)
      !write(577,*) numit

      end subroutine ci_cubic


      subroutine A_eqn(A, f, df,  Ra, b, K1, gamol,  ca, a1, f1, Rd )
!@sum Calculates the f coefficient in the coupled equ of photosynth/cond.
!@+   Igor, what is this solving for?? For f and df?
!@+   I.Aleinov
      real*8 A, f, df
      real*8 Ra, b, K1, gamol,  ca, a1, f1, Rd
      !---
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 cs, ci, dci
      real*8 byAKbcs, bycif1, K

      if ( A > 0.d0 ) then
        K = K1
      else
        K = 0.d0
      endif
      !write(579,*) "start A_eqn", A
      cs = ca - A*Ra
      byAKbcs = 1.d0/(A*K + b*cs)
      ci = cs * ( 1.d0 - A*byAKbcs )

      bycif1 = 1.d0/(ci+f1)
      f = A - ( a1*(ci-gamol)*bycif1 -Rd)

      dci = -Ra*( 1.d0 - A*byAKbcs )
     &     + cs*(- byAKbcs + A*byAKbcs*byAKbcs*(K-b*Ra) )
      df = 1 - a1*(f1+gamol)*bycif1*bycif1 * dci
      
      !write(579,*) "stop A_eqn", f, df

      end subroutine A_eqn

      subroutine A_eqn_0(A, f, Ra, b, K1, gamol,  ca, a1, f1, Rd )
!@sum Calculates coefficients in equation for coupled photosynth/cond.
!@auth I.Aleinov
      real*8 A, f
      real*8 Ra, b, K1, gamol,  ca, a1, f1, Rd
      !---
      real*8, parameter :: S_ATM=1.37d0  ! diffusivity ratio H2O/CO2 (atmosph.)
      real*8, parameter :: S_STOM=1.65d0 ! diffusivity ratio H2O/CO2 (stomatal)
      real*8 cs, ci, dci
      real*8 byAKbcs, bycif1, K

      if ( A > 0.d0 ) then
        K = K1
      else
        K = 0.d0
      endif
      !write(579,*) "start A_eqn", A
      cs = ca - A*Ra
      byAKbcs = 1.d0/(A*K + b*cs)
      ci = cs * ( 1.d0 - A*byAKbcs )

      bycif1 = 1.d0/(ci+f1)
      f = A - ( a1*(ci-gamol)*bycif1 -Rd)

      !write(579,*) "stop A_eqn", f, df

      end subroutine A_eqn_0



      FUNCTION rtsafe(funcd,x1,x2,xacc,  Ra, b, K, gamol,ca, a1, f1, Rd
     &     , numit )
!@sum Newton-Raphson solver (Numerical Recipes)
!@auth   I.Aleinov
      INTEGER MAXIT
      REAL*8 rtsafe,x1,x2,xacc
      real*8 Ra, b, K, gamol,  ca, a1, f1, Rd
      integer numit
      EXTERNAL funcd
      PARAMETER (MAXIT=100)
      INTEGER j
      REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl

cddd      ! for check
cddd      real*8 xxx
cddd      xxx = (x1+x2)/2.d0
cddd      call funcd(xxx,fl,df,  Ra, b, K, gamol,  ca, a1, f1, Rd)
cddd      xxx = xxx + .001d0
cddd      call funcd(xxx,fh,df,  Ra, b, K, gamol,  ca, a1, f1, Rd)
cddd      write(579,*) "deriv: ", (fh-fl)/.001d0, df
      
      numit = 0

      call A_eqn_0(x1,fl,  Ra, b, K, gamol,  ca, a1, f1, Rd)
      call A_eqn_0(x2,fh,  Ra, b, K, gamol,  ca, a1, f1, Rd)
      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
        rtsafe = -1.d30
        return ! for now return 0
        !call stop_model('root must be bracketed in rtsafe',255)
      endif
      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      !rtsafe=.5*(x1+x2)
      rtsafe=x1
      dxold=abs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df,  Ra, b, K, gamol,  ca, a1, f1, Rd)
      do 11 j=1,MAXIT
        numit = j
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or. abs(2.*
     *f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return
        call funcd(rtsafe,f,df,  Ra, b, K, gamol,  ca, a1, f1, Rd)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      call stop_model('rtsafe exceeding maximum iterations',255)
      return
      END FUNCTION rtsafe
#endif


      end module photcondmod

