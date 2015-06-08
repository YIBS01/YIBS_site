      MODULE CONSTANT
!@sum  CONSTANT definitions for physical constants and useful numbers
!@auth G. Schmidt
!@ver  1.0
      IMPLICIT NONE
      SAVE
C**** Conventions: 'by' implies reciprocal, 'rt' implies square root

C**** Numerical constants

      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: twopi = 2d0*pi           !@param twopi 2*pi
      real*8,parameter :: radian = pi/180d0        !@param radian pi/180
!@param zero,one 0 and 1 for occasional use as arguments
      real*8,parameter :: zero = 0d0, one=1d0
!@param rt2,byrt2   sqrt(2), 1/sqrt(2)
      real*8,parameter :: rt2 = 1.4142135623730950d0
      real*8,parameter :: byrt2 = 1./rt2
!@param rt3,byrt3   sqrt(3), 1/sqrt(3)
      real*8,parameter :: rt3 = 1.7320508075688772d0
      real*8,parameter :: byrt3 = 1./rt3
!@param rt12,byrt12   sqrt(12), 1/sqrt(12)
      real*8,parameter :: rt12 = 3.4641016151377546d0
      real*8,parameter :: byrt12 = 1./rt12
      real*8,parameter :: by3 =1./3d0  !@param by3  1/3
      real*8,parameter :: by6 =1./6d0  !@param by6  1/6
      real*8,parameter :: by9 =1./9d0  !@param by9  1/9
      real*8,parameter :: by12=1./12d0 !@param by12 1/12
!@param undef Missing value
      real*8,parameter :: undef=-1.d30
!@param teeny  small positive value used in num/(den+teeny) to avoid 0/0
      real*8,parameter :: teeny=1.d-30
      integer*8,parameter :: intNaN=-1  ! i.e. = Z'FFFFFFFFFFFFFFFF'
!@param NaN NaN
#if (defined COMPILER_PGI || defined COMPILER_NAG)
      real*8,parameter :: NaN=1d30
#else
      real*8,parameter :: NaN=transfer(intNaN,1.d0)
#endif

C**** Physical constants

!@param stbo Stefan-Boltzmann constant (W/m^2 K^4)
      real*8,parameter :: stbo =5.67051d-8 !current best estimate

c**** Latent heats:
c**** Note that for energy conservation the efective latent heat at any
c**** temperature must follow these formulae (assuming a reference
c**** temperature of 0 Celcius, and constant specific heats).
c**** If specific heats vary as a function of temperature, the extra
c**** term becomes an integral
c**** lhe(T) = lhe(0) + (shv-shw) T (in C)
c**** lhm(T) = lhm(0) + (shw-shi) T (in C)
c**** lhs(T) = lhs(0) + (shv-shi) T (in C)
!@param lhe   latent heat of evap at 0 C (2.5008d6 J/kg)
      real*8,parameter :: lhe = 2.5d6
!@param lhm   latent heat of melt at 0 C (334590 J/kg)
      real*8,parameter :: lhm = 3.34d5
!@param bylhm  1/lhm
      real*8,parameter :: bylhm = 1./lhm
!@param lhs  latent heat of sublimation at 0 C (J/kg)
      real*8,parameter :: lhs = lhe+lhm

!@param rhow density of pure water (1000 kg/m^3)
      real*8,parameter :: rhow = 1d3
!@param rhows density of average sea water (1030 kg/m^3)
      real*8,parameter :: rhows = 1030d0
!@param byrhows recip. density of average sea water (1/1030 m^3/kg)
      real*8,parameter :: byrhows = 1d0/rhows
!@param rhoi density of pure ice (916.6 kg/m^3)
      real*8,parameter :: rhoi = 916.6d0
!@param byrhoi 1/rhoi (m^3/kg)
      real*8,parameter :: byrhoi = 1d0/rhoi

!@param tf freezing point of water at 1 atm (273.16 K)
      real*8,parameter :: tf = 273.16d0
!@param bytf 1/tf (K^-1)
      real*8,parameter :: bytf = 1d0/tf

!@param shw heat capacity of water (at 20 C) (4185 J/kg C)
      real*8,parameter :: shw  = 4185.
!@param byshw 1/shw
      real*8,parameter :: byshw = 1d0/shw

!@param shi heat capacity of pure ice (at 0 C) (2060 J/kg C)
      real*8,parameter :: shi  = 2060.
!@param byshi 1/shi
      real*8,parameter :: byshi = 1d0/shi

c**** RGAS = R/M_A = 1000* 8.314510 J/mol K /28.9655 g/mol
c**** For values of CO2 much larger than present day (> 4x conc)
c**** the molecular weight of dry air M_A could change.
c**** Assume that M_O2 = 31.9988 and M_CO2 = 44.00995
c**** and current percentages 20.946% and 0.0350% (US Stand. Atm.)
c**** Assuming CO2 displaces other gases equally M_A=28.9602 + n*0.00527
c**** where n is multiple of present day CO2 conc (350 ppm)
c**** For 4xCO2  M_A = 28.9813  => rgas = 286.89
c**** For 10xCO2 M_A = 29.0129  => rgas = 286.58
!@param gasc  gas constant (8.314510 J/mol K)
      real*8,parameter :: gasc = 8.314510d0
!@param bygasc  1/gasc
      real*8,parameter :: bygasc = 1./gasc
!@param mair molecular weight of dry air (28.9655 g/mol)
      real*8,parameter :: mair = 28.9655d0
!@param rgas gas constant (287.05 J/K kg)
      real*8,parameter :: rgas = 1d3 * gasc / mair ! = 287.05...

!@param mwat molecular weight of water vapour
      real*8,parameter :: mwat = 18.015d0
!@param rvap  gas constant for water vapour (461.5 J/K kg)
c**** defined as R/M_W = 1000* 8.314510 J/mol K /18.015 g/mol
      real*8,parameter :: rvap = 1d3 * gasc / mwat ! = 461.5...

!@param mrat  mass ratio of air to water vapour (0.62197)
      real*8,parameter :: mrat = mwat/mair    ! = 0.62197....
!@param bymrat 1/mrat (1.6078)
      real*8,parameter :: bymrat = 1./mrat    ! = 1.6078....
!@param deltx coeff. of humidity in virtual temperature defn. (0.6078)
      real*8,parameter :: deltx = bymrat-1.   ! = 0.6078....

!@param srat ratio of specific heats at const. press. and vol. (=1.401)
      real*8,parameter :: srat = 1.401d0
!@param kapa ideal gas law exponent for dry air (.2862)
c**** kapa = (g-1)/g where g=1.401 = c_p/c_v
      real*8,parameter :: kapa = (srat - 1.)/srat  ! =.2862....
!@param bykapa,bykapap1,bykapap2 various useful reciprocals of kapa
      real*8,parameter :: bykapa = 1./kapa
      real*8,parameter :: bykapap1 = 1./(kapa+1.)
      real*8,parameter :: bykapap2 = 1./(kapa+2.)

!@param sha specific heat of dry air (const. pres.) (rgas/kapa J/kg C)
      real*8,parameter :: sha = rgas/kapa
!@param bysha 1/sha
      real*8,parameter :: bysha = 1./sha

!@param shv specific heat of water vapour (const. pres.) (J/kg C)
c**** shv is currently assumed to be zero to aid energy conservation in
c**** the atmosphere. Once the heat content associated with water
c**** vapour is included, this can be set to the standard value
c**** Literature values are 1911 (Arakawa), 1952 (Wallace and Hobbs)
c**** Smithsonian Met Tables = 4*rvap + delta = 1858--1869 ????
c     real*8,parameter :: shv = 4.*rvap  ????
      real*8,parameter :: shv = 0.

C**** air viscosity - temperature independent
!@var visc_air0 dynamic viscosity of air (kg/m s)
      real*8,parameter :: visc_air0 = 1.7d-5

!@var visc_air_kin0 kinematic viscosity of air (1 bar 15 deg C) (m^2/s)
      real*8,parameter :: visc_air_kin0 = 1.46d-5

!@var visc_wtr_kin kinematic viscosity of water (35 psu, 20 deg C) (m^2/s)
      real*8,parameter :: visc_wtr_kin = 1.05d-6

!@var avog Avogadro's constant (atmos/mole)
      real*8,parameter :: avog=6.023d23

C**** Astronomical constants

!@param sday  sec per day (s)
      real*8,parameter :: sday = 86400.
!@param syr  sec per year (s)
      real*8,parameter :: syr = sday*365.

!@param hrday  hours in a day (hrs)
      real*8,parameter :: hrday = sday/3600.

!@param omega earth's rotation rate (7.29 s^-1)
c      real*8,parameter :: omega = 7.2921151467d-5 ! NOVAS value
      real*8,parameter :: EDPERD = 1.
      real*8,parameter :: EDPERY = 365.
      real*8,parameter :: omega = TWOPI*(EDPERD+EDPERY)/
     *                            (EDPERD*EDPERY*SDAY)
!@param omega2 2*omega
      real*8,parameter :: omega2 = 2.*omega

!@param radius radius of the earth (6371000 m, IUGG)
      real*8,parameter :: radius = 6371000.
!@param areag surface area of the earth (m^2)
      real*8,parameter :: areag = 4.*pi*radius*radius

!@param grav gravitaional accelaration (9.80665 m/s^2)
c**** SI reference gravity (at 45 deg) = 9.80665
      real*8,parameter :: grav = 9.80665d0
!@param bygrav 1/grav
      real*8,parameter :: bygrav = 1d0/grav

C**** lapse rate related variables
!@param GAMD dry adiabatic lapse rate (=0.0098 K/m)
      real*8, parameter :: gamd = grav*kapa/rgas
!@param BMOIST moist adiabatic lapse rate (K/m)
      real*8, parameter :: bmoist = 0.0065d0
!@param BBYG moist adiabatic lapse rate divided by grav
      real*8, parameter :: bbyg = bmoist*bygrav
!@param GBYRB grav divided by rgas and bmoist
      real*8, parameter :: gbyrb = grav/(rgas*bmoist)

C**** Useful conversion factors

!@param kg2mb,mb2kg conversion from milli-bars to kg/m^2
      real*8,parameter :: kg2mb = 1d-2*grav, mb2kg = 1d2*bygrav
!@param kgpa2mm,mm2kgpa conversion from kg/m^2 water to mm
      real*8,parameter :: kgpa2mm = 1d0, mm2kgpa = 1d0

      CONTAINS

      real*8 function visc_air(T)
!@sum visc_air dynamic viscosity of air (function of T) (kg/m s)
!@auth Sutherland formula
      real*8, intent(in) :: T  ! temperature (K)
      real*8, parameter :: n0=1.827d-5, T0=291.15d0, C=120d0

      visc_air = n0*sqrt((T/T0)**3)*(T0+C)/(T+C)

      return
      end function

      real*8 function visc_air_kin(T)
!@sum visc_air_kin kinematic viscosity of air (function of T) (m2/s)
!@auth COARE formula - Andreas (1989) CRREL Rep. 89-11
      real*8, intent(in) :: T  ! temperature (K)
      real*8, parameter :: nu0=1.326d-5, a0=6.542d-3, b0=8.301d-6,
     *     c0=4.84d-9
      real*8 :: Tc  ! temperature in deg C

      Tc=T-tf
      visc_air_kin = nu0*(1.+Tc*(a0+Tc*(b0-c0*Tc)))   !m2/s

      return
      end function

      END MODULE CONSTANT
