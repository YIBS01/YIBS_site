C     *****************
C     *****************
      SUBROUTINE HANDLE_ERR(NSTAT)
C     *****************
C     *****************

      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER NSTAT
      IF (NSTAT .NE. NF_NOERR) THEN
        PRINT*, NF_STRERROR(NSTAT)
      STOP
      ENDIF

      RETURN
      END

