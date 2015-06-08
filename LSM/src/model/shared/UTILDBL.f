#include "rundeck_opts.h"
!@sum  UTILDBL Model Independent Utilities
!@auth Original Development Team
!@ver  1.0
!@cont THBAR,QSAT,DQSATDT,READT

#if ( defined USE_ESMF )  || ( defined USE_MPP )
#define USE_MPI
#endif

      FUNCTION THBAR (X,Y)
!@sum  THBAR calculates mean temperature used in vertical differencing
!@auth Gary Russell, Jean Lerner, Arakawa
!@ver  1.0
C****
C**** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
C****              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
C****      g(x) is replaced by a rational function
C****           (a+bx+cxx+dxxx+xxxx)/(e+fx+gxx)
C****      approx.error <1.E-6 for x between .9 and 1.7
C****
      IMPLICIT NONE
!@var A,B,C,D,E,F,G   expansion coefficients for THBAR
      REAL*8, PARAMETER :: A=113.4977618974100d0
      REAL*8, PARAMETER :: B=438.5012518098521d0
      REAL*8, PARAMETER :: C=88.49964112645850d0
      REAL*8, PARAMETER :: D=-11.50111432385882d0
      REAL*8, PARAMETER :: E=30.00033943846368d0
      REAL*8, PARAMETER :: F=299.9975118132485d0
      REAL*8, PARAMETER :: G=299.9994728900967d0
      REAL*8 :: Q,AL                 !@var Q,AL   working variables
      REAL*8, INTENT(IN) :: X,Y      !@var X,Y    input temperatures
      REAL*8 :: THBAR                !@var THBAR  averaged temperature
      Q=X/Y
      AL=(A+Q*(B+Q*(C+Q*(D+Q))))/(E+Q*(F+G*Q))
      THBAR=X*AL
      RETURN
      END

      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : mrat,rvap,tf
      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
C**** LH = 0.5*(LH(0)+LH(t))
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(LH*(B-C/max(130.d0,TM)))/PR
      RETURN
      END

      FUNCTION DQSATDT (TM,LH)
!@sum  DQSATDT calculates change of sat. vapour mixing ratio with temp.
!@auth Gary Russell
!@ver  1.0
C**** Note that d(qsat)/dt = qsat * lh * c / T*T
C**** Only the factor of qsat is given here
      USE CONSTANT, only : rvap
      IMPLICIT NONE
!@var C coefficient for QSAT
      REAL*8, PARAMETER :: C = 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in DQSATDT is the actual LH at TM i.e. LH=LH(TM)
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: DQSATDT         !@var DQSATDT d(qsat)/dT factor only.
      DQSATDT = LH*C/(TM*TM)    ! * QSAT(TM,LH,PR)
      RETURN
      END

      FUNCTION SLP(PS,TAS,ZS)
!@sum SLP estimates sea level pressure in the presence of topography
!@+   for better match to reanalyses.
      USE CONSTANT, only: bmoist, grav, rgas, by3
      IMPLICIT NONE
!@var PS surface pressure (mb)
!@var TAS surface temperature (K)
!@var ZS surface elevation (m)
      REAL*8, INTENT(IN) :: PS, TAS, ZS
      REAL*8 :: SLP, TSL, BETA, BZBYT, GBYRB, TASn

      IF (ZS.ne.0.) THEN
        TSL= TAS+BMOIST*ZS
        TASn=TAS
        BETA=BMOIST
        IF (TAS < 290.5 .and. TSL > 290.5) BETA= (290.5d0 - TAS)/ZS
        IF (TAS > 290.5 .and. TSL > 290.5) TASn = 0.5*(290.5d0 + TAS)
        IF (TAS < 255) TASn = 0.5*(255d0 + TAS)
        BZBYT=BETA*ZS/TASn
        GBYRB=GRAV/(RGAS*BETA)
        IF (BETA > 1d-6 ) THEN
          SLP=PS*(1.+BZBYT)**GBYRB
        ELSE
          SLP=PS*EXP((1.-0.5*BZBYT+BZBYT**(2.*by3))*GBYRB*BZBYT)
        END IF
      ELSE
        SLP=PS
      END IF
      RETURN
      END FUNCTION SLP


      SUBROUTINE READT (IUNIT,NSKIP,LENGTH,AOUT,IPOS)
!@sum   READT  read in title and real*4 array and convert to real*8
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT        !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: NSKIP    !@var  NSKIP  no. of R*4's to skip
      INTEGER, INTENT(IN) :: LENGTH       !@var  LENGTH size of array
      INTEGER, INTENT(IN) :: IPOS  !@var  IPOS   no. of recs. to advance
      !REAL*4, INTENT(OUT) :: AIN(LENGTH)  !@var  AIN    real*4 array
      REAL*8, INTENT(OUT) :: AOUT(LENGTH) !@var  AOUT   real*8 array
      REAL*4 :: X               !@var  X      dummy variable
      INTEGER :: N              !@var  N      loop variable
      CHARACTER*80 TITLE        !@var  TITLE  title of file record
      real*4, allocatable :: buf(:)

      allocate( buf(LENGTH) )

      DO N=1,IPOS-1
        READ (IUNIT,END=920)
      END DO
      READ (IUNIT,ERR=910,END=920) TITLE,(X,N=1,NSKIP),buf
C**** do transfer backwards in case AOUT and AIN are same workspace
!!! THIS DOESN''T WORK IN F90+  !!! stop using such hacks !
!      DO N=LENGTH,1,-1
!        AOUT(N)=AIN(N)
!      END DO
      AOUT(:) = buf(:)
      deallocate( buf )
      WRITE(6,*) "Read from file ",TRIM(NAME(IUNIT)),": ",TRIM(TITLE)
      RETURN
  910 WRITE(6,*) 'READ ERROR ON FILE ',NAME(IUNIT)
      call stop_model('tREAD: READ ERROR',255)
  920 WRITE(6,*) 'END OF FILE ENCOUNTERED ON FILE ',NAME(IUNIT)
      call stop_model('tREAD: No data found',255)
      END

      subroutine WRITEI (iunit,it,aout,len4)
!@sum   WRITEI  writes array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT       !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT          !@var  IT time, 1st & last word
      INTEGER, INTENT(IN) :: LEN4        !@var  LENGTH size of array
      REAL*4,  INTENT(IN) :: AOUT(LEN4)  !@var  AOUT   real*4 array

      write (iunit) it,aout,it
      call sys_flush(iunit)
      write (6,*) "Wrote to file ",TRIM(NAME(IUNIT)),", time=",it
      return
      END subroutine WRITEI

      subroutine READI (iunit,it,ain,it1,len4,iok)
!@sum  READI reads array surrounded by IT (for post processing)
!@auth  Original Development Team
!@ver   1.0
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT,LEN4
      INTEGER, INTENT(out) :: IT,it1,iok ! iok: 0='ok',1='not ok'
      real*4  AIN(LEN4)
      iok = 0
      read(iunit,end=555) it,ain,it1
      return
  555 iok=1
      return
      end subroutine readi

      subroutine WRITEI8 (iunit,it,aout,len8)
!@sum   WRITEI8 writes real*8 array surrounded by IT and secures it
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT       !@var  IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT          !@var  IT time, 1st & last word
      INTEGER, INTENT(IN) :: LEN8        !@var  LENGTH size of array
      REAL*8,  INTENT(IN) :: AOUT(LEN8)  !@var  AOUT   real*8 array

      write (iunit) it,aout,it
      call sys_flush(iunit)
      write (6,*) "Wrote to file ",TRIM(NAME(IUNIT)),", time=",it
      return
      END subroutine WRITEI8

      subroutine READI8 (iunit,it,ain,it1,len8,iok)
!@sum  READI reads array surrounded by IT (for post processing)
!@auth  Original Development Team
!@ver   1.0
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT,LEN8
      INTEGER, INTENT(out) :: IT,it1,iok ! iok: 0='ok',1='not ok'
      real*8  AIN(LEN8)
      iok = 0
      read(iunit,end=555) it,ain,it1
      return
  555 iok=1
      return
      end subroutine readi8

      subroutine io_POS (iunit,it,len4,itdif)
!@sum   io_POS  positions a seq. output file for the next write operat'n
!@auth  Original Development Team
!@ver   1.0
!@var NAME name of record being read
      USE FILEMANAGER, only : NAME=>nameunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIT !@var IUNIT  file unit number
      INTEGER, INTENT(IN) :: IT,ITdif !@var IT,ITdif current time,dt
      INTEGER, INTENT(IN) :: LEN4 !@var LENGTH of array in words
      INTEGER :: IT1,IT2   !@var  time_tags at start,end of each record
      INTEGER :: N              !@var  N      loop variable

      read (iunit,end=10,err=50) it1,(it2,n=1,len4+1)
      if(it1 .le. it) go to 30
   10 write(6,*) "Starting a new file ",TRIM(NAME(IUNIT)),", time=",it
      rewind iunit
      return

   20 read (iunit,end=35,err=50) it1,(it2,n=1,len4+1)
   30 if (it2 .ne. it1) then
        write(6,*) 'file ',TRIM(NAME(IUNIT)),' damaged: it/it1/it2=',
     *    it,it1,it2
        call stop_model('io_POS: damaged file',255)
      end if
      if (it1 .le. it) go to 20
      it1=it1-itdif
   35 backspace iunit
      if (it1+itdif .le. it) go to 40
      write (6,*) "positioned ",TRIM(NAME(IUNIT)),", it1/itime=",it1,it
      return
   40 write (6,*) "file ",TRIM(NAME(IUNIT))," too short, it1/it=",it1,it
      call stop_model('io_POS: file too short',255)
   50 write (6,*) "Read error on: ",TRIM(NAME(IUNIT)),", it1/it=",it1,it
      call stop_model('io_POS: read error',255)
      END subroutine io_POS

      SUBROUTINE CHECK3(A,IN,JN,LN,SUBR,FIELD)
!@sum  CHECK3 Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: IN,JN,LN
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(IN,JN,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK3 = .FALSE.
      INTEGER I,J,L !@var I,J,L loop variables

!$OMP PARALLEL DO PRIVATE (L,J,I) SHARED (QCHECK3)
      DO L=1,LN
      DO J=1,JN
      DO I=1,IN
        IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.) .or.
     *       ABS(A(I,J,L)) .gt.HUGE(A(I,J,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,L,A(I,J,L),'after ',SUBR
          IF (J.LT.JN.AND.J.GT.1) QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) call stop_model('CHECK3',255)
      RETURN
      END SUBROUTINE CHECK3

      SUBROUTINE CHECK3B(A,I1,I2,J1,J2,NJPOL,LN,SUBR,FIELD)
!@sum  CHECK3B Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: I1,I2,J1,J2,NJPOL,LN
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(I1:I2,J1:J2,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK3 = .FALSE.
      INTEGER I,J,L !@var I,J,L loop variables

!$OMP PARALLEL DO PRIVATE (L,J,I) SHARED (QCHECK3)
      DO L=1,LN
      DO J=J1+NJPOL,J2-NJPOL
      DO I=I1,I2
        IF (.NOT.(A(I,J,L).GT.0..OR.A(I,J,L).LE.0.) .or.
     *       ABS(A(I,J,L)) .gt.HUGE(A(I,J,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,L,A(I,J,L),'after ',SUBR
          QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) call stop_model('CHECK3',255)
      RETURN
      END SUBROUTINE CHECK3B

      SUBROUTINE CHECK3C(A,LN,I1,I2,J1,J2,NJPOL,SUBR,FIELD)
!@sum  CHECK3B Checks for NaN/INF in real 3-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,LN size of 3-D array
      INTEGER, INTENT(IN) :: LN,I1,I2,J1,J2,NJPOL
!@var SUBR identifies where CHECK3 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(LN,I1:I2,J1:J2),INTENT(IN) :: A
      LOGICAL :: QCHECK3 = .FALSE.
      INTEGER I,J,L !@var I,J,L loop variables

!$OMP PARALLEL DO PRIVATE (L,J,I) SHARED (QCHECK3)
      DO J=J1+NJPOL,J2-NJPOL
      DO I=I1,I2
      DO L=1,LN
        IF (.NOT.(A(L,I,J).GT.0..OR.A(L,I,J).LE.0.) .or.
     *       ABS(A(L,I,J)) .gt.HUGE(A(L,I,J)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',L,I,J,A(L,I,J),'after ',SUBR
          QCHECK3 = .TRUE.
        END IF
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK3) call stop_model('CHECK3',255)
      RETURN
      END SUBROUTINE CHECK3C

      SUBROUTINE CHECK4(A,IN,JN,KN,LN,SUBR,FIELD)
!@sum  CHECK4 Checks for NaN/INF in real 4-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,KN,LN size of 4-D array
      INTEGER, INTENT(IN) :: IN,JN,KN,LN
!@var SUBR identifies where CHECK4 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(IN,JN,KN,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK4 = .FALSE.
      INTEGER I,J,K,L !@var I,J,K,L loop variables

!$OMP PARALLEL DO PRIVATE (L,K,J,I) SHARED (QCHECK4)
      DO L=1,LN
      DO K=1,KN
      DO J=1,JN
      DO I=1,IN
        IF (.NOT.(A(I,J,K,L).GT.0..OR.A(I,J,K,L).LE.0.) .or.
     *       ABS(A(I,J,K,L)) .gt.HUGE(A(I,J,K,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,K,L,A(I,J,K,L),'after ',SUBR
          IF (J.LT.JN.AND.J.GT.1) QCHECK4 = .TRUE.
        END IF
      END DO
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK4) call stop_model('CHECK4',255)
      RETURN
      END SUBROUTINE CHECK4

      SUBROUTINE CHECK4B(A,I1,I2,J1,J2,NJPOL,KN,LN,SUBR,FIELD)
!@sum  CHECK4 Checks for NaN/INF in real 4-D arrays
!@auth Original development team
!@ver  1.0
      IMPLICIT NONE

!@var IN,JN,KN,LN size of 4-D array
      INTEGER, INTENT(IN) :: I1,I2,J1,J2,NJPOL,KN,LN
!@var SUBR identifies where CHECK4 was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var FIELD identifies the field being tested
      CHARACTER*6, INTENT(IN) :: FIELD
!@var A array being tested
      REAL*8, DIMENSION(I1:I2,J1:J2,KN,LN),INTENT(IN) :: A
      LOGICAL :: QCHECK4 = .FALSE.
      INTEGER I,J,K,L !@var I,J,K,L loop variables

!$OMP PARALLEL DO PRIVATE (L,K,J,I) SHARED (QCHECK4)
      DO L=1,LN
      DO K=1,KN
      DO J=J1+NJPOL,J2-NJPOL
      DO I=I1,I2
        IF (.NOT.(A(I,J,K,L).GT.0..OR.A(I,J,K,L).LE.0.) .or.
     *       ABS(A(I,J,K,L)) .gt.HUGE(A(I,J,K,L)) ) THEN
          WRITE (6,*) TRIM(FIELD),': ',I,J,K,L,A(I,J,K,L),'after ',SUBR
          QCHECK4 = .TRUE.
        END IF
      END DO
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL SYS_FLUSH(6)
      IF (QCHECK4) call stop_model('CHECK4',255)
      RETURN
      END SUBROUTINE CHECK4B

      function unit_string (pow10,ending)
!@sum Construct a units string with nice properties (no embedded blanks)
!@auth G. Schmidt, J. Lerner
C**** If a trailing ')' is supplied, it is assumed that a leading
C****      '(' is required, so it is inserted
      implicit none
      character*(*) ending,unit_string
      character*10 tpow
      integer pow10

      tpow = ' '
      if(pow10.ne.0) then
        write(tpow,'(i3)') pow10
        if (index(ending,')') .ne.0) then
          tpow='(10^'//trim(adjustl(tpow))
        else
          tpow= '10^'//trim(adjustl(tpow))
        end if
      endif
      unit_string = adjustl(trim(tpow)//" "//trim(ending))
      return
      end function unit_string



      subroutine write_run_status( message, retcode )
      implicit none
      character*(*), intent (in) :: message
      integer, intent(in) :: retcode
      integer, parameter :: iu_status = 9
      character*10 :: form_str
      integer num_digits

      ! construct format string in such a way that retcode is printed
      ! at the beginning of the line with no extra spaces
      if ( retcode .ne. 0 ) then
        num_digits = log10( real( abs(retcode), kind(1.d0) ) ) + 1
      else
        num_digits = 1
      endif
      if ( retcode < 0 ) num_digits = num_digits + 1

      write(form_str,"('(I',I1,')')") num_digits

      open( iu_status, file='run_status', form='FORMATTED',
     &     status='UNKNOWN', ERR=10 )
      write( iu_status, form_str, ERR=10 ) retcode
      write( iu_status, '(A)', ERR=10 ) message
      close( iu_status )

      return
 10   continue
      write( 0, * ) "ERROR: Can't write to the run_status file"
      write( 0, * ) "STATUS:", message
      end subroutine write_run_status

      module precision_mod
!@sum  The reduce_precision routines truncate the number of
!@+    significant digits in a real*8 number x to an approximate
!@+    precision of relacc (1d-16 < relacc << 1).  Fortran functions
!@+    can be used to define x = fraction(x) * 2**exponent(x),
!@+    where fraction(x) = O(1).  The part of fraction(x) smaller than
!@+    relacc is discarded.
!@auth M. Kelley
!@ver  1.0
      implicit none
      public :: reduce_precision
      interface reduce_precision
        module procedure reduce_precision_0d
        module procedure reduce_precision_1d
        module procedure reduce_precision_2d
        module procedure reduce_precision_3d
        module procedure reduce_precision_4d
      end interface reduce_precision
      contains
      subroutine reduce_precision_0d(x,relacc)
      real*8, intent(inout) :: x
      real*8, intent(in) :: relacc
      x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
      end subroutine reduce_precision_0d
      subroutine reduce_precision_1d(x,relacc)
      real*8, intent(inout) :: x(:)
      real*8, intent(in) :: relacc
      x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
      end subroutine reduce_precision_1d
      subroutine reduce_precision_2d(x,relacc)
      real*8, intent(inout) :: x(:,:)
      real*8, intent(in) :: relacc
      x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
      end subroutine reduce_precision_2d
      subroutine reduce_precision_3d(x,relacc)
      real*8, intent(inout) :: x(:,:,:)
      real*8, intent(in) :: relacc
      x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
      end subroutine reduce_precision_3d
      subroutine reduce_precision_4d(x,relacc)
      real*8, intent(inout) :: x(:,:,:,:)
      real*8, intent(in) :: relacc
      x = nint(fraction(x)/relacc,kind=8)*relacc*2d0**exponent(x)
      end subroutine reduce_precision_4d

      end module precision_mod

      function clean_str(string)
!@sum clean_str utility to clean strings of netcdf-unfriendly characters
      implicit none
      character(len=*), intent(in) :: string
      character(len=40) :: clean_str
      integer :: k

      clean_str=trim(string)
      do k=1,len_trim(string)
         if (clean_str(k:k).eq." " .or. clean_str(k:k).eq."+")
     $        clean_str(k:k)="_"
      end do
      return
      end function clean_str

