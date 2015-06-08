MODULE FILEMANAGER
!@sum  FILEMANAGER keeps data concerning the files and unit numbers
!@ver  2.0
  implicit none
  save
  private

  public openunit, closeunit, print_open_units,findunit
  public nameunit

  interface openunit
    module procedure openunit_0d
    module procedure openunit_1d
  end interface

  interface closeunit
    module procedure closeunit_0d
    module procedure closeunit_1d
  end interface

!@param MINUNIT, MAXUNIT - min and max unit number allowed
  integer, parameter :: MINUNIT = 50, MAXUNIT = 528
  integer, parameter :: TOTALUNITS = MAXUNIT-MINUNIT+1

  type UnitStr
    logical in_use                     ! is the unit in use or not
    character*80 filename              ! the name on the file
  end type UnitStr

  type (UnitStr) :: Units( MINUNIT:MAXUNIT ) = UnitStr(.false.," ")

contains

  function nameunit( unit )
!@sum nameunit returns the name of the file corresponding to unit <unit>
    implicit none
    character*80 nameunit
    integer, intent(in) :: unit

    if ( unit>MAXUNIT .or. unit<MINUNIT .or. .not. Units(unit)%in_use ) then
      write(6,*) "FILEMANAGER: asked name of a wrong unit: ",unit
      call stop_model("FILEMANAGER: asked name of a wrong unit",255)
    endif
    nameunit = Units(unit)%filename
  end function nameunit


  subroutine findunit( unit )
!@sum findunit finds available unit
    implicit none
    integer, intent(out) :: unit
    logical :: opened

    !ia next line is needed for Absoft compiler on Linux (compiler bug?)
    unit = 0
    do unit=MINUNIT,MAXUNIT
      if ( unit > 98 .and. unit < 103 ) cycle
      !if ( .not. Units(unit)%in_use ) return
      inquire( unit=unit, opened=opened )
      if ( .not. opened ) return
    enddo
    write(6,*) "FILEMANAGER: Maximum file number reached"
    call print_open_units
    call stop_model("FILEMANAGER: Maximum file number reached",255)
  end subroutine findunit


  subroutine openunit_0d( filename, iunit, qbin, qold )
!@sum openunit opens the file <filename> and returns its unit in <unit>
    implicit none
!@var unit - unit of opened file
    integer, intent(out) :: iunit
!@var filename - name of the file to open
    character*(*), intent(in) :: filename
!@var qbin - .true. if file is binary (.false. by default)
!@var qold - .true. if file is old (.false. by default)
    logical, optional, intent(in) :: qbin, qold
    character*11 form, status
    integer name_len

!**** set default options
    form = "FORMATTED"
    status = "UNKNOWN"
!**** parse options
    if( present(qbin) ) then
      if( qbin ) form = "UNFORMATTED"
    end if
    if( present(qold) ) then
      if( qold ) status = "OLD"
    end if

    call findunit( iunit )

    !dbg  write(6,*) "FILEMANAGER: Before Opening file ",trim(filename) !RKF debug
    if ( form == "FORMATTED" ) then
      open( iunit, FILE=filename, FORM=form, STATUS=status, &
#ifdef CONVERT_BIGENDIAN
      &       CONVERT='BIG_ENDIAN', &
#endif
      &       RECL=65534, ERR=10 )
    else
      open( iunit, FILE=filename, FORM=form, STATUS=status, &
#ifdef CONVERT_BIGENDIAN
      &       CONVERT='BIG_ENDIAN', &
#endif
      &       ERR=10 )
    endif

    !dbg  write(6,*) "FILEMANAGER: Opened file ",trim(filename) !RKF debug

    Units(iunit)%in_use = .true.
    name_len = len_trim(filename)
    Units(iunit)%filename = filename( max(1,name_len-79) : name_len )
    return

10  write(6,*) "FILEMANAGER: Error opening file ",trim(filename)
    call stop_model('FILEMANAGER: FILE OPENING ERROR',255)
  end subroutine openunit_0d

  subroutine closeunit_0d( unit )
!@sum closeunit closes the file the file corresponding to unit <unit>
    implicit none
!@var unit - unit of the file to close
    integer, intent(in) :: unit

    if ( unit>MAXUNIT .or. unit<MINUNIT .or. .not. Units(unit)%in_use ) then
      write(6,*) "FILEMANAGER: attempt to close wrong unit: ",unit
      call stop_model("FILEMANAGER: attempt to close wrong unit",255)
    endif

    close( unit )
    Units(unit)%in_use = .false.
    Units(unit)%filename = " "
  end subroutine closeunit_0d


  subroutine openunit_1d( filename, iunit, qbin, qold )
!@sum  openunit_1d sets unit number for requested files and opens them
!@auth Gavin Schmidt, simplified by I. Aleinov
!@ver  1.0
!@var unit - unit of opened file
    integer, intent(out), dimension(:) :: iunit
!@var filename - name of the file to open
    character*(*), intent(in), dimension(:) :: filename
!@var qbin - .true. if file is binary (.false. by default)
!@var qold - .true. if file is old (.false. by default)
    logical, optional, intent(in) :: qbin, qold
    integer i !@var i loop variable

    do i=1,size(filename)
      call openunit( filename(i), iunit(i), qbin, qold )
    end do
  end subroutine openunit_1d


  subroutine closeunit_1d( iunit )
!@sum  closeunit_1d closes files corresponding to units iunit
!@auth Gavin Schmidt, I. Aleinov
!@ver  1.0
    implicit none
!@var iunit unit numbers for files in current request
    integer, intent(in), dimension(:) :: iunit
!@var nreq number of file unit numbers
    integer i

    do i=1,size(iunit)
      call closeunit( iunit(i) )
    end do
  end subroutine closeunit_1d


  subroutine print_open_units
!@sum print_open_units prints info on open units (for debugging)
    implicit none
    integer unit

    write(6,*) "FILEMANAGER: Open Units:"
    do unit=MINUNIT,MAXUNIT
      if ( Units(unit)%in_use ) then
        write(6,*) "unit = ", unit, "file = ", Units(unit)%filename
      endif
    enddo
  end subroutine print_open_units

END MODULE FILEMANAGER
