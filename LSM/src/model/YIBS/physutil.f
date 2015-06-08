      module physutil
!@sum Physics utility functions.
!@auth N.Y.Kiang

      use yibs_const
      implicit none

      public QSAT

      contains
!============================================================================
      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0 (I think this is at least version 2.0)
!      USE CONSTANT, only : mrat,rvap,tf
      IMPLICIT NONE
!@var Physical constants from GISS GCM CONST.f
      real*8, parameter :: MWAT = 18.015d0 !molecular weight of water vapour (g/mol)
      real*8, parameter :: MAIR = 28.9655d0 !molecular weight of dry air (28.9655 g/mol)
      real*8, parameter :: MRAT = MWAT/MAIR 
      real*8, parameter :: RVAP = 1d3 * GASC/MWAT !gas constant for water vapour (461.5 J/K kg)
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TFRZ)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
C**** LH = 0.5*(LH(0)+LH(t)), where LH(0)=
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(LH*(B-C/max(130.d0,TM)))/PR
      RETURN
      END FUNCTION QSAT
!============================================================================

      end module physutil
