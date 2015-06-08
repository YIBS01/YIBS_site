module JulianCalendar_mod
!@sum Specifies the parameters for the Julian calendar used in modelE.
!@auth T. Clune
  implicit none


! During the transition period, this module will continue to support the legacy
! names for constants which tend towards the terse end of the spectrum.  


  !         Legacy :  New
  public :: SDAY,     SECONDS_PER_DAY
  public :: SYR,      SECONDS_PER_YEAR
  public ::           SECONDS_PER_HOUR
  public :: HRDAY,    HOURS_PER_DAY

  public :: JDPERY,   DAYS_PER_YEAR
  public :: JMPERY,   MONTHS_PER_YEAR
  public :: JDendOfM, LAST_JULIAN_DAY_IN_MONTH
  public :: JDmidOfM, MID_JULIAN_DAY_IN_MONTH

  public :: AMONTH

  ! Not entirely clear if these constants belong here, since they actually
  ! should vary for other planets and/or eras (e.g. paleo has shorter day)
  public :: EDPERD,   EARTH_DAYS_PER_DAY
  public :: EDPERY,   EARTH_DAYS_PER_YEAR

  public :: Month_type, JULIAN_MONTHS
  public :: JANUARY,   FEBRUARY, MARCH,    APRIL
  public :: MAY,       JUNE,     JULY,     AUGUST
  public :: SEPTEMBER, OCTOBER,  NOVEMBER, DECEMBER

  
!@var DAYS_PER_YEAR (JDPERY)    number of days per year
  integer, parameter :: DAYS_PER_YEAR = 365, JDPERY = DAYS_PER_YEAR
!@var MONTHS_PER_YEAR (JMperY)  number of months per year
  integer, parameter :: MONTHS_PER_YEAR = 12, JMPERY = MONTHS_PER_YEAR

  real*8, parameter :: SECONDS_PER_DAY = 86400.,   SDAY = SECONDS_PER_DAY
  real*8, parameter :: EARTH_DAYS_PER_DAY = 1.,    EDPERD = EARTH_DAYS_PER_DAY
  real*8, parameter :: EARTH_DAYS_PER_YEAR = 365., EDPERY = EARTH_DAYS_PER_YEAR
  real*8, parameter :: SECONDS_PER_YEAR = SECONDS_PER_DAY * DAYS_PER_YEAR
  real*8, parameter :: SYR = SECONDS_PER_YEAR

  real*8, parameter :: SECONDS_PER_HOUR = 3600.
  real*8, parameter :: HOURS_PER_DAY = SECONDS_PER_DAY / SECONDS_PER_HOUR
  real*8, parameter :: HRDAY = HOURS_PER_DAY


!@var LAST_JULIAN_DAY_IN_MONTH (JDendOfM, ) last Julian day in month
  integer, parameter :: LAST_JULIAN_DAY_IN_MONTH(0:MONTHS_PER_YEAR) = (/ &
       & 0,31,59,90,120,151,181,212,243,273,304,334,365 &
       & /)
  integer, parameter :: JDendOfM(0:MONTHS_PER_YEAR) = LAST_JULIAN_DAY_IN_MONTH
!@var MID_JULIAN_DAY_IN_MONTH(0:13) (JDmidOfM(0:13)) middle Julian day in month
  integer, parameter :: MID_JULIAN_DAY_IN_MONTH(0:MONTHS_PER_YEAR+1) = (/ &
       & -15,16,45,75,106,136,167,197,228,259,289,320,350,381 &
       & /)
  integer, parameter :: JDmidOfM(0:MONTHS_PER_YEAR+1) = MID_JULIAN_DAY_IN_MONTH

  ! Months
  type Month_type
    character(len=4) :: shortName
    character(len=20) :: longName
    integer :: numDays
    integer :: lastJulianDay
    integer :: middleJulianDay
  end type Month_type

  type (Month_type), parameter :: JANUARY   = Month_type('JAN ', 'January   ', 31,  31,  16)
  type (Month_type), parameter :: FEBRUARY  = Month_type('FEB ', 'February  ', 28,  59,  45)
  type (Month_type), parameter :: MARCH     = Month_type('MAR ', 'March     ', 31,  90,  75)
  type (Month_type), parameter :: APRIL     = Month_type('APR ', 'April     ', 30, 120,  106)
  type (Month_type), parameter :: MAY       = Month_type('MAY ', 'May       ', 31, 151,  136)
  type (Month_type), parameter :: JUNE      = Month_type('JUNE', 'June      ', 30, 181,  167)
  type (Month_type), parameter :: JULY      = Month_type('JULY', 'July      ', 31, 212,  197)
  type (Month_type), parameter :: AUGUST    = Month_type('AUG ', 'August    ', 31, 243,  228)
  type (Month_type), parameter :: SEPTEMBER = Month_type('SEP ', 'September ', 30, 273,  259)
  type (Month_type), parameter :: OCTOBER   = Month_type('OCT ', 'October   ', 31, 304,  289)
  type (Month_type), parameter :: NOVEMBER  = Month_type('NOV ', 'November  ', 30, 334,  320)
  type (Month_type), parameter :: DECEMBER  = Month_type('DEC ', 'December  ', 31, 365,  350)

  type (Month_type), parameter :: JULIAN_MONTHS(MONTHS_PER_YEAR) = (/ &
       & JANUARY,  FEBRUARY, MARCH,     &
       & APRIL,    MAY,      JUNE,      &
       & JULY,     AUGUST,   SEPTEMBER, &
       & OCTOBER,  NOVEMBER, DECEMBER   &
       & /)

!@var AMONTH(0:12)  (3-4 letter) names for all months
! AMONTH(0) = 'IC' (i.e. initial conditions) only used early in a
! model run.  Should find a way to eliminate it.
#ifdef COMPILER_PGI
  ! hack to work around PGI internal error
  character*4, parameter :: AMONTH(0:12) = (/ 'IC  ','JAN ','FEB ','MAR ', &
   'APR ','MAY ','JUN ','JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/)
#else
  character*4, parameter :: AMONTH(0:12) = (/ 'IC  ', JULIAN_MONTHS%shortName /)
#endif



end module JulianCalendar_mod
