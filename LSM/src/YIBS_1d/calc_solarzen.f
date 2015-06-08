      subroutine calc_solarzen(td,latdegrees,sbeta1)
      !* Calculate solar zenith angle **in radians**
      !* From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
      implicit none
      real*8,intent(in) :: td             ! day(to minute fraction)
      real*8,intent(in) :: latdegrees     ! latitude in degrees
      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: rad = pi/180.d0 ! Conversion from degrees to radians.
      real*8 :: hour,latrad
      real*8 :: delta                     ! declination angle
      real*8 :: td0
      real*8,intent(out) :: sbeta1        ! sbeta1=cos(zen angle)=sin(elev angle)
!      real*8,intent(out) :: solarelev    ! solar elevation angle (rad)
!      real*8,intent(out) :: solarzen     ! solar zenith angle (rad)
      
      td0 = td
      If (td0 .lt. 0.d0) td0 = td0 + 365.0d0
      If (td0 .gt. 365.0d0) td0 = td0 - 365.0d0
      hour = (td0-floor(td0))*24.d0
      latrad = latdegrees*rad
      delta = asin(-sin(rad*23.45d0)*cos(2.d0*pi*(td0+10.d0)/365.d0))
      sbeta1 = sin(latrad)*sin(delta)+
     &     cos(latrad)*cos(delta)*cos(rad* 15.d0*(hour-12.d0))
c      if (sbeta1 < 0.d0) sbeta1 = 0.d0  !**GCM does this too** 
 
      end subroutine calc_solarzen

