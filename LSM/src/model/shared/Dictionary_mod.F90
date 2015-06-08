!TODO delete alloc_param
module Dictionary_mod
!@sum  Provides interfaces to manipulate sets of key-value pairs.
!@+  This type of data structure is also knows as an associative array
!@+  and a map.
!@auth I. Aleinov & T.Clune
!@ver 1.1
!@usage This module has the following public subroutines. Most of them
!@+ are Fortramn 90 interfaces, i.e. they recognize the type of the
!@+ variables implicitly and call corresponding subroutine.
!@+
!@+ Simple copy routines to copy parameters to/from database:
!@+     set_param( name, value, dim, opt ) - put a parameter with the
!@+       name <name> and the value <value> to the database
!@+     get_param( name, value, dim ) - copy the value of the parameter
!@+       <name> from the database to the variable <value>
!@+
!@+ Query logical function:
!@+     is_set_param( name ) - returns .true. if parameter <name> is
!@+       present in the database, .false. otherwise
!@+     sync_param( name, value, dim ) - puts parameter <name> into the
!@+       database if it is not yet there, otherwise gets it
!@+
!@+ Reading/writing subroutines:
!@+     read_param( kunit, ovrwrt ) - reads the parameter database from
!@+       the unit <kunit>
!@+     write_param( kunit ) - writes the parameter database to the unit
!@+       <kunit>
!@+
!@+ Other useful subroutines:
!@+     print_param( kunit ) - does formatted output to the unit <kunit>
!@+       in a way similar to namelist
!@+     query_param( n, name, dim, ptype ) - returns information about
!@+       the parameter by its number in the database <n>. It returns
!@+       'EMPTY' in the <name> if no parameter with such an <n> exists
!@+       (i.e. if <n> is bigger than the number of parameters)
!@+
!@+ The formal arguments in the subroutines are:
!@+     name - character*(*) - the name of the parameter which is a
!@+       character string no longer than 32 bytes
!@+     value - a scalar variable or a linear array of type: integer,
!@+       real*8,character*1 to character*128
!@+     dim - integer - dimension of an array; omit 'dim' for scalars
!@+     opt - character*1 - an optional "option" (opt='o' means
!@+       "overwrite")
!@+     kunit - integer - unit number for reading/writing
!@+     ptype - character*1, intent(out) - returns the type of the
!@+       parameter: 'i' for integer, 'r' for real*8, 'c' for character
!@+     ovrwrt - logical - if .true. then reading overwrites those
!@+       parameters that are already in the database. If .false. then
!@+       those parameters which are already in the database are left
!@+       unchanged and only new parameters are added.
!@+
!@+   Read FAQ's for the full description.
!@+
!@+               CHANGE LOG:
!@+ 04/18/02 added 3 bytes to ParamStr so that its size is
!@+ divisible by 4 (needed for portability SGI,LINUX <-> IBM,COMPAQ).
!@+ Header renamed to "PARAM02 "
  use StringUtilities_mod, only: toLowerCase
  use KeyValuePair_mod
  use GenericType_mod
  implicit none
  save
  private

  ! Dictionary
  public :: Dictionary_type ! data type
  public :: Dictionary      ! constructor
  public :: clean           ! destructor
  public :: insert
  public :: getNumEntries
  public :: lookup
  public :: merge
  public :: hasKey
  public :: getKeys
  public :: readUnformatted
  public :: writeUnformatted
  public :: operator(==)


  ! Accessors

  
  public set_param, get_param, read_param, write_param
  public is_set_param, sync_param, print_param
  public query_param, print_unused_param
  public :: reset

  ! params
  public :: DICT_TYPE
  integer, parameter :: DICT_TYPE = 5

  public :: MAX_LEN_LINE
  integer, parameter :: MAX_LEN_LINE  = 256

  integer, parameter :: NOT_FOUND = -1

  type Dictionary_type
    private
    type (KeyValuePair_type), pointer :: pairs(:) => null()
  end type Dictionary_type

  integer, parameter :: MAX_PARAMS = 345
  integer, parameter :: MAX_RPARAMS = 288
  integer, parameter :: MAX_IPARAMS = 555
  integer, parameter :: MAX_CPARAMS = 205
  integer, parameter :: MAX_NAME_LEN = 32
  integer, parameter :: MAX_CHAR_LEN = 128

  character*80 :: MODULE_HEADER='PARAM02 '

  type ParamStr
    character(MAX_NAME_LEN) :: name = 'EMPTY'  ! parameter name
    integer :: indx = 0                 ! storage for its value
    integer :: dim = 0                  ! number of elements
    character(1) :: attrib = 'u'        ! type: real ('r') or int ('i')
    character(1) :: is_accessed = 'n'   ! was it accessed with get_*
    character(1) :: source = 'u'        ! where it came from (r=rundeck)
    character(1) :: reserved(1) = 'u'
  end type ParamStr

  type (ParamStr), target :: Params(MAX_PARAMS)
  real*8, target :: Rdata(MAX_RPARAMS)
  integer, target :: Idata(MAX_IPARAMS)
  character*(MAX_CHAR_LEN), target :: Cdata(MAX_CPARAMS)
  integer :: num_param = 0
  integer :: num_rparam = 0
  integer :: num_iparam = 0
  integer :: num_cparam = 0

  interface set_param
    module procedure set_iparam, set_rparam, set_cparam
    module procedure  set_aiparam, set_arparam, set_acparam
  end interface

  interface get_param
    module procedure get_iparam, get_rparam, get_cparam
    module procedure get_aiparam, get_arparam, get_acparam
  end interface

  interface sync_param
    module procedure sync_iparam, sync_rparam, sync_cparam
    module procedure sync_aiparam, sync_arparam, sync_acparam
  end interface

  ! Constructors
  interface Dictionary
    module procedure Dictionary_empty
    module procedure Dictionary_copy
  end interface

  interface clean
    module procedure cleanDictionary
  end interface

  interface insert
    module procedure insert_pair
    module procedure insert_integer
    module procedure insert_real64
    module procedure insert_logical
    module procedure insert_string
    module procedure insert_integerArray
    module procedure insert_real64Array
    module procedure insert_logicalArray
    module procedure insert_stringArray
  end interface

  interface merge
    module procedure merge_dictionary
    module procedure merge_pair
    module procedure merge_integer
    module procedure merge_real64
    module procedure merge_logical
    module procedure merge_string
    module procedure merge_integerArray
    module procedure merge_real64Array
    module procedure merge_logicalArray
    module procedure merge_stringArray
  end interface
  
  interface operator(==)
    module procedure equals
  end interface

  interface readUnformatted
    module procedure readUnformatted_dictionary
  end interface

  interface writeUnformatted
    module procedure writeUnformatted_dictionary
  end interface

  interface getKeys
    module procedure getKeys_dictionary
  end interface

contains

  function Dictionary_empty() 
    type (Dictionary_type) :: Dictionary_empty
    allocate(Dictionary_empty%pairs(0))
  end function Dictionary_empty

  function Dictionary_copy(original) result(copy)
    type (Dictionary_type), intent(in) :: original
    type (Dictionary_type) :: copy

    integer :: i

    allocate(copy%pairs(getNumEntries(original)))
    do i = 1, getNumEntries(original)
      copy%pairs(i) = KeyValuePair(original%pairs(i))
    end do

  end function Dictionary_copy

  subroutine reset()
    num_param = 0
    num_rparam = 0
    num_iparam = 0
    num_cparam = 0
  end subroutine reset

  function is_set_param( name_in )
    implicit none
    character*(*), intent(in) :: name_in
    logical is_set_param
    integer n

    do n=1,num_param
      if ( Params(n)%name == toLowerCase(name_in) ) exit
    enddo

    if ( n > num_param  ) then
      is_set_param = .false.
    else
      is_set_param = .true.
    endif

    return
  end function is_set_param

  subroutine set_pstr( name_in, dim, attrib, PStr, flag )
    implicit none
    character*(*), intent(in) :: name_in
    integer, intent(in) :: dim
    character*1, intent(in) ::  attrib
    type (ParamStr), pointer :: PStr
    logical, intent(in) :: flag

    if ( len(name_in) > MAX_NAME_LEN ) then
      print *, 'PARAM: parameter name too long: ', name_in
      print *, 'PARAM: maximal length allowed: ', MAX_NAME_LEN
      call stop_model('PARAM: parameter name too long: ',255)
    endif

    call get_pstr( toLowerCase(name_in), dim, attrib, PStr )

    if ( associated( PStr ) ) then
      if ( .not. flag ) then
        print *, 'PARAM: attempt to set param which is already set'
        print *, 'name: ', toLowerCase(name_in)
        call stop_model( &
             &         'PARAM: attempt to set param which is already set',255)
      else
        return  ! return PStr found by get_pstr
      endif
    endif

    if ( num_param >= MAX_PARAMS ) then
      print *, 'PARAM: Maximal number of parameters exceeded'
      print *, 'PARAM: Please recompile param with bigger MAX_PARAMS'
      print *, 'At least ',num_param+1,' is required'
      call stop_model( &
           &       'PARAM: Maximal number of parameters exceeded',255)
    endif

    num_param = num_param + 1
    PStr => Params(num_param)
    PStr%name = toLowerCase(name_in)
    PStr%attrib = attrib
    PStr%dim = dim

    select case (attrib)
    case ('i')
      if ( num_iparam+dim >= MAX_IPARAMS ) then
        print *, 'PARAM: Maximal number of int parameters exceeded'
        print *, 'PARAM: Recompile param with bigger MAX_IPARAMS'
        print *, 'At least ',num_iparam+dim+1,' is required'
        call stop_model( &
             &           'PARAM: Maximal number of int parameters exceeded',255)
      endif
      PStr%indx = num_iparam + 1
      num_iparam = num_iparam + dim
    case ('r')
      if ( num_rparam+dim >= MAX_RPARAMS ) then
        print *, 'PARAM: Need ',num_rparam+dim+1
        print *, 'PARAM: Maximal number of real parameters exceeded'
        print *, 'PARAM: Recompile param with bigger MAX_RPARAMS'
        print *, 'At least ',num_rparam+dim+1,' is required'
        call stop_model( &
             &          'PARAM: Maximal number of real parameters exceeded',255)
      endif
      PStr%indx = num_rparam + 1
      num_rparam = num_rparam + dim
    case ('c')
      if ( num_cparam+dim >= MAX_CPARAMS ) then
        print *, 'PARAM: Maximal number of char parameters exceeded'
        print *, 'PARAM: Recompile param with bigger MAX_CPARAMS'
        print *, 'At least ',num_cparam+dim+1,' is required'
        call stop_model( &
             &          'PARAM: Maximal number of char parameters exceeded',255)
      endif
      PStr%indx = num_cparam + 1
      num_cparam = num_cparam + dim
    end select

    return
  end subroutine set_pstr


  subroutine get_pstr( name_in, dim, attrib, PStr )
    implicit none
    character*(*), intent(in) :: name_in
    integer, intent(in) :: dim
    character*1, intent(in) ::  attrib
    type (ParamStr), pointer :: PStr
    integer n

    nullify( PStr )

    do n=1,num_param
      if ( Params(n)%name == toLowerCase(name_in) ) exit
    enddo

    if ( n > num_param  ) return   ! not found - return NULL

    if ( Params(n)%attrib /= attrib .or. Params(n)%dim /= dim ) then
      print *, 'PARAM: wrong type or dim of parameter: ', toLowerCase(name_in)
      print *, 'ATT: set: ', Params(n)%attrib, ' called: ', attrib
      print *, 'DIM: set: ', Params(n)%dim, ' called: ', dim
      call stop_model('PARAM: wrong type or dim of parameter',255)
    endif

    PStr => Params(n)
    return
  end subroutine get_pstr



  !***** integers ******!


  subroutine set_iparam( name, value, opt )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: value
    character(len=*), optional, intent(in) :: opt
    integer v(1)
    v(1) = value
    call set_aiparam( name, v, 1, opt )
    return
  end subroutine set_iparam


  subroutine set_aiparam( name, value, np, opt )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: np
    integer, intent(in) :: value(np)
    character(len=*), optional, intent(in) :: opt
    type (ParamStr), pointer :: PStr
    logical flag
    character*1 source

    flag = .false.
    source = 'u'
    if ( present(opt) ) then
      if ( scan(opt,'o') .ne. 0 ) flag = .true.
      if ( scan(opt,'r') .ne. 0 ) source = 'r'
    endif

    call set_pstr( name, np, 'i', PStr, flag )
    Idata( PStr%indx : PStr%indx+np-1 ) = value(1:np)
    PStr%source = source
    return
  end subroutine set_aiparam


  subroutine get_iparam( name, value )
    implicit none
    character*(*), intent(in) ::  name
    integer, intent(out) ::  value
    integer v(1)

    call get_aiparam( name, v, 1 )
    value = v(1)
    return
  end subroutine get_iparam


  subroutine get_aiparam( name, value, np )
    implicit none
    character*(*), intent(in) ::  name
    integer, intent (in) :: np
    integer, intent(out) ::  value(np)
    type (ParamStr), pointer :: PStr

    call get_pstr( name, np, 'i', PStr )
    if ( .not. associated( PStr) ) then
      print *, 'PARAM: Can''t get - not in database : ', name
      call stop_model( &
           &       'PARAM: Can''t get parameter - not in database',255)
    endif
    value(1:np) = Idata( PStr%indx : PStr%indx+np-1 )
    PStr%is_accessed = 'y'
    return
  end subroutine get_aiparam


  !***** reals ******!

  subroutine set_rparam( name, value, opt )
    implicit none
    character*(*), intent(in) :: name
    real*8, intent(in) :: value
    character(len=*), optional, intent(in) :: opt
    real*8 v(1)
    v(1) = value
    call set_arparam( name, v, 1, opt )
    return
  end subroutine set_rparam


  subroutine set_arparam( name, value, np, opt )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: np
    real*8, intent(in) :: value(np)
    character(len=*), optional, intent(in) :: opt
    type (ParamStr), pointer :: PStr
    logical flag
    character*1 source

    flag = .false.
    source = 'u'
    if ( present(opt) ) then
      if ( scan(opt,'o') .ne. 0 ) flag = .true.
      if ( scan(opt,'r') .ne. 0 ) source = 'r'
    endif

    call set_pstr( name, np, 'r', PStr, flag )
    Rdata( PStr%indx : PStr%indx+np-1 ) = value(1:np)
    PStr%source = source
    return
  end subroutine set_arparam


  subroutine get_rparam( name, value )
    implicit none
    character*(*), intent(in) ::  name
    real*8, intent(out) ::  value
    real*8 v(1)

    call get_arparam( name, v, 1 )
    value = v(1)
    return
  end subroutine get_rparam


  subroutine get_arparam( name, value, np )
    implicit none
    character*(*), intent(in) ::  name
    integer, intent (in) :: np
    real*8, intent(out) ::  value(np)
    type (ParamStr), pointer :: PStr

    call get_pstr( name, np, 'r', PStr )
    if ( .not. associated( PStr) ) then
      print *, 'PARAM: Can''t get - not in database : ', name
      call stop_model( &
           &       'PARAM: Can''t get parameter - not in database',255)
    endif
    value(1:np) = Rdata( PStr%indx : PStr%indx+np-1 )
    PStr%is_accessed = 'y'
    return
  end subroutine get_arparam



  !***** Chars ******!

  subroutine set_cparam( name, value, opt )
    implicit none
    character*(*), intent(in) :: name
    character*(*), intent(in) :: value
    character(len=*), optional, intent(in) :: opt
    character*(MAX_CHAR_LEN) v(1)
    if ( len(value) > MAX_CHAR_LEN ) then
      print *, 'PARAM: Char string too long. MAX = ', MAX_CHAR_LEN
      call stop_model('PARAM: Char string too long',255)
    endif
    v(1) = value
    call set_acparam( name, v, 1, opt )
    return
  end subroutine set_cparam


  subroutine set_acparam( name, value, np, opt )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: np
    character*(*), intent(in) :: value(np)
    character(len=*), optional, intent(in) :: opt
    type (ParamStr), pointer :: PStr
    integer n
    logical flag
    character*1 source

    flag = .false.
    source = 'u'
    if ( present(opt) ) then
      if ( scan(opt,'o') .ne. 0 ) flag = .true.
      if ( scan(opt,'r') .ne. 0 ) source = 'r'
    endif

    do n=1,np
      if ( len(value(n)) > MAX_CHAR_LEN ) then
        print *, 'PARAM: Char string too long. MAX = ', MAX_CHAR_LEN
        print *, 'You submitted LEN = ', len(value(n))
        call stop_model('PARAM: Char string too long',255)
      endif
    enddo
    call set_pstr( name, np, 'c', PStr, flag )
    Cdata( PStr%indx : PStr%indx+np-1 ) = value(1:np)
    PStr%source = source
    return
  end subroutine set_acparam


  subroutine get_cparam( name, value )
    implicit none
    character*(*), intent(in) ::  name
    character*(*), intent(out) ::  value
    character*(MAX_CHAR_LEN) v(1)

    call get_acparam( name, v, 1 )
    value = v(1)
    return
  end subroutine get_cparam


  subroutine get_acparam( name, value, np )
    implicit none
    character*(*), intent(in) ::  name
    integer, intent (in) :: np
    character*(*), intent(out) ::  value(np)
    type (ParamStr), pointer :: PStr

    call get_pstr( name, np, 'c', PStr )
    if ( .not. associated( PStr) ) then
      print *, 'PARAM: Can''t get - not in database : ', name
      call stop_model( &
           &       'PARAM: Can''t get parameter - not in database',255)
    endif
    value(1:np) = Cdata( PStr%indx : PStr%indx+np-1 )
    PStr%is_accessed = 'y'
    return
  end subroutine get_acparam


  !***** sync functions ******!

  subroutine sync_iparam( name, value )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(inout) :: value

    if ( is_set_param( name ) ) then
      call get_param( name, value )
    else
      call set_param( name, value )
    endif
  end subroutine sync_iparam


  subroutine sync_aiparam( name, value, np )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: np
    integer, intent(inout) :: value(np)

    if ( is_set_param( name ) ) then
      call get_param( name, value, np )
    else
      call set_param( name, value, np )
    endif
  end subroutine sync_aiparam


  subroutine sync_rparam( name, value )
    implicit none
    character*(*), intent(in) :: name
    real*8, intent(inout) :: value

    if ( is_set_param( name ) ) then
      call get_param( name, value )
    else
      call set_param( name, value )
    endif
  end subroutine sync_rparam


  subroutine sync_arparam( name, value, np )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: np
    real*8, intent(inout) :: value(np)

    if ( is_set_param( name ) ) then
      call get_param( name, value, np )
    else
      call set_param( name, value, np )
    endif
  end subroutine sync_arparam


  subroutine sync_cparam( name, value )
    implicit none
    character*(*), intent(in) :: name
    character*(*), intent(inout) :: value

    if ( is_set_param( name ) ) then
      call get_param( name, value )
    else
      call set_param( name, value )
    endif
  end subroutine sync_cparam


  subroutine sync_acparam( name, value, np )
    implicit none
    character*(*), intent(in) :: name
    integer, intent(in) :: np
    character*(*), intent(inout) :: value(np)

    if ( is_set_param( name ) ) then
      call get_param( name, value, np )
    else
      call set_param( name, value, np )
    endif
  end subroutine sync_acparam



  !***** input / output ******!

  subroutine read_param( kunit, ovrwrt )
    implicit none
    integer, intent(in) :: kunit
    logical, intent(in) :: ovrwrt
    integer n
    type (ParamStr), save :: LParams(MAX_PARAMS)
    real*8, save :: LRdata(MAX_RPARAMS)
    integer, save :: LIdata(MAX_IPARAMS)
    character*(MAX_CHAR_LEN), save :: LCdata(MAX_CPARAMS)
    integer lnum_param, lnum_rparam, lnum_iparam, lnum_cparam
    character*80 HEADER

    read( kunit, err=10 ) HEADER
    backspace kunit

    if (HEADER(1:8).ne.MODULE_HEADER(1:8)) then
      if (HEADER(1:8).eq.'PARAM01 ') then
        print *, 'WARNING: PARAM: Old format of parameter data.'
        call read_param_comp01( kunit, ovrwrt )
        return
      else
        print * , 'PARAM: No parameter header in input data'
        call stop_model( &
             &         'PARAM: No parameter header in input data',255)
      endif
    endif

    read( kunit, err=10 ) HEADER, &
         &     lnum_param, lnum_rparam, lnum_iparam, lnum_cparam, &
         &     ( LParams(n), n=1,min(lnum_param,MAX_PARAMS) ), &
         &     ( LRdata(n), n=1,min(lnum_rparam,MAX_RPARAMS) ), &
         &     ( LIdata(n), n=1,min(lnum_iparam,MAX_IPARAMS) ), &
         &     ( LCdata(n), n=1,min(lnum_cparam,MAX_CPARAMS) )

    if (     lnum_param  > MAX_PARAMS &
         &     .or. lnum_rparam > MAX_RPARAMS &
         &     .or. lnum_iparam > MAX_IPARAMS &
         &     .or. lnum_cparam > MAX_CPARAMS &
         &     ) then
      print *, 'PARAM: parameter list in input file too long'
      print *, 'PARAM: please recompile param with bigger MAX_?PARAMS'
      print *,'PARAM: ',lnum_param,lnum_rparam,lnum_iparam,lnum_cparam
      call stop_model( &
           &       'PARAM: parameter list in input file too long',255)
    endif

    if ( lnum_param < 1 ) return   ! no parameters in the records

#if defined(MACHINE_DEC) || defined(MACHINE_Linux)
    ! COMPAQ needs big/little endian conversion
    do n=1,lnum_param
      call swap_bytes_4( LParams(n)%indx, 1 )
      call swap_bytes_4( LParams(n)%dim,  1 )
    enddo
#endif

    ! checking big/little endian format, just in case
    if ( LParams(1)%dim > 65536 .or. LParams(1)%dim < 0 ) then
      print *, 'PARAM: wrong big/little endian format in LParams.'
      call stop_model( &
           &       'PARAM: wrong big/little endian format in LParams',255)
    endif

    ! now merge the data just read with existing database
    do n=1,lnum_param
      if ( (.not. is_set_param(LParams(n)%name)) .or. ovrwrt ) then
        select case( LParams(n)%attrib )
        case ('i')
          call set_aiparam( LParams(n)%name, LIdata(LParams(n)%indx), &
               &           LParams(n)%dim, 'o' )
        case ('r')
          call set_arparam( LParams(n)%name, LRdata(LParams(n)%indx), &
               &           LParams(n)%dim, 'o' )
        case ('c')
          call set_acparam( LParams(n)%name, LCdata(LParams(n)%indx), &
               &           LParams(n)%dim, 'o' )
        end select
      endif
    enddo
    return
10  print *, 'PARAM: Error reading, unit = ', kunit
    call stop_model('PARAM: Error reading',255)
  end subroutine read_param


  subroutine write_param( kunit )
    implicit none
    integer, intent(in) :: kunit
    integer n

    write (MODULE_HEADER(9:80),'(i10,a)') &
         &  num_param,' is the current number of parameters in database DB'

#if defined(MACHINE_DEC) || defined(MACHINE_Linux)
    ! converting it manually to big-endian for COMPAQ compiler
    do n=1,num_param
      call swap_bytes_4( Params(n)%indx, 1 )
      call swap_bytes_4( Params(n)%dim,  1 )
    enddo
#endif

    write( kunit, err=20 ) MODULE_HEADER, &
         &     num_param, num_rparam, num_iparam, num_cparam, &
         &     ( Params(n), n=1,min(num_param,MAX_PARAMS) ), &
         &     ( Rdata(n), n=1,min(num_rparam,MAX_RPARAMS) ), &
         &     ( Idata(n), n=1,min(num_iparam,MAX_IPARAMS) ), &
         &     ( Cdata(n), n=1,min(num_cparam,MAX_CPARAMS) )

#if defined(MACHINE_DEC) || defined(MACHINE_Linux)
    ! and back to little-endian ...
    do n=1,num_param
      call swap_bytes_4( Params(n)%indx, 1 )
      call swap_bytes_4( Params(n)%dim,  1 )
    enddo
#endif
    return

20  print *, 'PARAM: Error writing, unit = ', kunit
    call stop_model('PARAM: Error writing',255)
  end subroutine write_param

  subroutine print_param1( kunit )
    implicit none
    integer, intent(in) :: kunit
    integer, parameter :: nf = 7
    integer n, i

    write( kunit, * ) '&&PARAMETERS'
    do n=1, num_param
      select case( Params(n)%attrib )
      case ('i')
        write( kunit, '(1x,a16,a3,8i16)' )  &
             &         Params(n)%name, ' = ', &
             &        ( Idata(Params(n)%indx+i), i=0,min(Params(n)%dim,nf)-1 )
        if ( Params(n)%dim > nf ) &
             &         write( kunit, '(20x,8i16)' ) &
             &        ( Idata(Params(n)%indx+i), i=0,Params(n)%dim-nf-1 )
      case ('r')
        write( kunit, '(1x,a16,a3,8g16.6)' ) &
             &         Params(n)%name, ' = ', &
             &        ( Rdata(Params(n)%indx+i), i=0,min(Params(n)%dim,nf)-1 )
        if ( Params(n)%dim > nf ) &
             &         write( kunit, '(20x,8g16.6)' ) &
             &        ( Rdata(Params(n)%indx+i), i=0,Params(n)%dim-nf-1 )
      case ('c')
        write( kunit, '(1x,a16,a3,8a128)' ) &
             &         Params(n)%name, ' = ', &
             &        ( Cdata(Params(n)%indx+i), i=0,min(Params(n)%dim,nf)-1 )
        if ( Params(n)%dim > nf ) &
             &         write( kunit, '(20x,8a128)' ) &
             &        ( Cdata(Params(n)%indx+i), i=0,Params(n)%dim-nf-1 )
      end select
    enddo
    write( kunit, * ) '&&END_PARAMETERS'
  end subroutine print_param1

  subroutine print_param( kunit )
    implicit none
    integer, intent(in) :: kunit
    integer, parameter :: nf = 7
    integer n, i

    write( kunit, * ) '&&PARAMETERS'
    do n=1, num_param
      select case( Params(n)%attrib )
      case ('i')
        write( kunit, * ) &
             &         trim(Params(n)%name), ' = ', &
             &        ( Idata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
      case ('r')
        write( kunit, * ) &
             &         trim(Params(n)%name), ' = ', &
             &        ( Rdata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
      case ('c')
        write( kunit, * ) &
             &         trim(Params(n)%name), ' = ', &
             &        ( Cdata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
      end select
    enddo
    write( kunit, * ) '&&END_PARAMETERS'
  end subroutine print_param


  subroutine print_unused_param( kunit )
    implicit none
    integer, intent(in) :: kunit
    integer, parameter :: nf = 7
    integer n, i

    do n=1, num_param
      if ( Params(n)%source .ne. 'r' ) cycle
      if ( Params(n)%is_accessed == 'y' ) cycle
      select case( Params(n)%attrib )
      case ('i')
        write(kunit, *) 'WARNING: defined in rundeck but not used: ', &
             &         trim(Params(n)%name), ' = ', &
             &        ( Idata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
      case ('r')
        write(kunit, *) 'WARNING: defined in rundeck but not used: ', &
             &         trim(Params(n)%name), ' = ',&
             &        ( Rdata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
      case ('c')
        write(kunit, *) 'WARNING: defined in rundeck but not used: ', &
             &         trim(Params(n)%name), ' = ', &
             &        ( Cdata(Params(n)%indx+i), i=0,Params(n)%dim-1 )
      end select
    enddo
  end subroutine print_unused_param



  subroutine query_param( n, name, dim, ptype )
    integer, intent(in) :: n
    character*(*), intent(out) :: name
    integer, intent(out) :: dim
    character*1, intent(out) :: ptype

    if ( n>0 .and. n<=num_param ) then
      name = Params(n)%name
      dim = Params(n)%dim
      ptype = Params(n)%attrib
    else
      name = 'EMPTY'
      dim = 0
      ptype = 'U'
    endif

  end subroutine query_param

  !**** the code below is included for compatibility with older     ****
  !**** versions; it may be removed later when not needed any more  ****

  subroutine read_param_comp01( kunit, ovrwrt )
    implicit none
    type ParamStr_comp01
      character(MAX_NAME_LEN) name ! parameter name
      integer indx                 ! storage for its value
      integer dim                   ! number of elements
      character*1 attrib            ! type: real ('r') or int ('i')
    end type ParamStr_comp01
    integer, intent(in) :: kunit
    logical, intent(in) :: ovrwrt
    integer n
    type (ParamStr_comp01), save :: LParams(MAX_PARAMS)
    real*8, save :: LRdata(MAX_RPARAMS)
    integer, save :: LIdata(MAX_IPARAMS)
    character*(MAX_CHAR_LEN), save :: LCdata(MAX_CPARAMS)
    integer lnum_param, lnum_rparam, lnum_iparam, lnum_cparam
    character*80 HEADER

    read( kunit, err=10 ) HEADER
    backspace kunit

    if (HEADER(1:8).ne.'PARAM01 ') then
      print *, 'PARAM: No parameter header in input data'
      call stop_model('PARAM: No parameter header in input data',255)
    endif

    read( kunit, err=10 ) HEADER, &
         &     lnum_param, lnum_rparam, lnum_iparam, lnum_cparam, &
         &     ( LParams(n), n=1,min(lnum_param,MAX_PARAMS) ), &
         &     ( LRdata(n), n=1,min(lnum_rparam,MAX_RPARAMS) ), &
         &     ( LIdata(n), n=1,min(lnum_iparam,MAX_IPARAMS) ), &
         &     ( LCdata(n), n=1,min(lnum_cparam,MAX_CPARAMS) )

    if (     lnum_param  > MAX_PARAMS &
         &     .or. lnum_rparam > MAX_RPARAMS &
         &     .or. lnum_iparam > MAX_IPARAMS &
         &     .or. lnum_cparam > MAX_CPARAMS &
         &     ) then
      print *, 'PARAM: parameter list in input file too long'
      print *, 'PARAM: please recompile param with bigger MAX_?PARAMS'
      print *,'PARAM: ',lnum_param,lnum_rparam,lnum_iparam,lnum_cparam
      call stop_model( &
           &       'PARAM: parameter list in input file too long',255)
    endif

    ! now merge the data just read with existing database
    do n=1,lnum_param
      if ( (.not. is_set_param(LParams(n)%name)) .or. ovrwrt ) then
        select case( LParams(n)%attrib )
        case ('i')
          call set_aiparam( LParams(n)%name, LIdata(LParams(n)%indx), &
               &           LParams(n)%dim, 'o' )
        case ('r')
          call set_arparam( LParams(n)%name, LRdata(LParams(n)%indx), &
               &           LParams(n)%dim, 'o' )
        case ('c')
          call set_acparam( LParams(n)%name, LCdata(LParams(n)%indx), &
               &           LParams(n)%dim, 'o' )
        end select
      endif
    enddo
    return
10  print *, 'PARAM: Error reading, unit = ', kunit
    call stop_model('PARAM: Error reading',255)
  end subroutine read_param_comp01

  subroutine addEntry(this, key)
!@sum Creates space for new entry in dictionary.
!@+ Could be altered to grow the list by multiple entries
!@+ which would result in fewer allocations.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key

    type (KeyValuePair_type), pointer :: oldPairs(:)
    integer :: numEntries

    if (hasKey(this, key)) then
      call throwException('Dictionary: duplicate key - <'//trim(key)//'>.',14)
      return
    end if

    oldPairs => this%pairs
    numEntries = getNumEntries(this)
    allocate(this%pairs(numEntries+1))
    this%pairs(1:numEntries) = oldPairs  ! shallow copy
    deallocate(oldPairs)

  end subroutine addEntry

  ! Begin overload interface for insert()

  subroutine insert_pair(this, pair)
!@sum Insert a KeyValuePair into dictionary.
    type (Dictionary_type), intent(inout) :: this
    type (KeyValuePair_type), intent(in) :: pair

    character(len=MAX_LEN_KEY) :: key

    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(pair) ! deep copy

  end subroutine insert_pair

  subroutine insert_integer(this, key, value)
!@sum Insert an integer into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value

    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(value))

  end subroutine insert_integer

  subroutine insert_real64(this, key, value)
!@sum Insert a double into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    real*8, intent(in) :: value
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(value))

  end subroutine insert_real64

  subroutine insert_logical(this, key, value)
!@sum Insert a logical into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(value))

  end subroutine insert_logical

  subroutine insert_string(this, key, value)
!@sum Insert a string into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(value))

  end subroutine insert_string

  subroutine insert_integerArray(this, key, values)
!@sum Insert an integer array into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: values(:)
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(values))

  end subroutine insert_integerArray

  subroutine insert_real64Array(this, key, values)
!@sum Insert a real array into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    real*8, intent(in) :: values(:)
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(values))

  end subroutine insert_real64Array

  subroutine insert_logicalArray(this, key, values)
!@sum Insert a logical array into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: values(:)
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(values))

  end subroutine insert_logicalArray

  subroutine insert_stringArray(this, key, values)
!@sum Insert a string array into dictionary.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: values(:)
    
    call addEntry(this, key)
    this%pairs(getNumEntries(this)) = KeyValuePair(key, GenericType(values))

  end subroutine insert_stringArray
  ! End overload interface for insert()
  
  ! Begin overload interface for merge()
  subroutine merge_dictionary(this, other)
!@sum Merge two dictionaries. Where duplicate keys exist
!@+ result has values from the original (1st argument).
!@+ Analogous behavior to old 'sync_param'.
    type (Dictionary_type), intent(inout) :: this
    type (Dictionary_type), intent(in) :: other

    integer :: i

    do i = 1, getNumEntries(other)
      call merge(this, other%pairs(i))
    end do
  end subroutine merge_dictionary

  subroutine merge_pair(this, pair)
!@sum Merge pair into dictionary.   If key not already
!@+ present, this is equivalent to insert(), otherwise
!@+ it does nothing.
    type (Dictionary_type), intent(inout) :: this
    type (KeyValuePair_type), intent(in) :: pair

    if (hasKey(this, getKey(pair))) return
    call insert(this, pair)

  end subroutine merge_pair

  subroutine merge_integer(this, key, value)
!@sum Merge key+integer.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(inout) :: value

    if (hasKey(this, key)) then
      value = lookup(this, key)
    else
      call insert(this, key, value)
    end if

  end subroutine merge_integer

  subroutine merge_real64(this, key, value)
!@sum Merge key+double.
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    real*8, intent(inout) :: value

    if (hasKey(this, key)) then
      value = lookup(this, key)
    else
      call insert(this, key, value)
    end if

  end subroutine merge_real64

  subroutine merge_logical(this, key, value)
!@sum Merge key+logical
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(inout) :: value

    if (hasKey(this, key)) then
      value = lookup(this, key)
    else
      call insert(this, key, value)
    end if

  end subroutine merge_logical

  subroutine merge_string(this, key, value)
!@sum Merge key+string
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(inout) :: value

    if (hasKey(this, key)) then
      value = lookup(this, key)
    else
      call insert(this, key, value)
    end if

  end subroutine merge_string

  subroutine merge_integerArray(this, key, values)
!@sum Merge key+integer array
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(inout) :: values(:)

    if (hasKey(this, key)) then
      values = lookup(this, key)
    else
      call insert(this, key, values)
    end if

  end subroutine merge_integerArray

  subroutine merge_real64Array(this, key, values)
!@sum Merge key+double array
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    real*8, intent(inout) :: values(:)

    if (hasKey(this, key)) then
      values = lookup(this, key)
    else
      call insert(this, key, values)
    end if

  end subroutine merge_real64Array

  subroutine merge_logicalArray(this, key, values)
!@sum Merge key+logical array
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(inout) :: values(:)

    if (hasKey(this, key)) then
      values = lookup(this, key)
    else
      call insert(this, key, values)
    end if

  end subroutine merge_logicalArray

  subroutine merge_stringArray(this, key, values)
!@sum Merge key+stringarray
    type (Dictionary_type), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(inout) :: values(:)


    if (hasKey(this, key)) then
      values = lookup(this, key)
    else
      call insert(this, key, values)
    end if

  end subroutine merge_stringArray
  ! End of overload of merge()

  
  integer function getNumEntries(this)
!@sum Returns number of entries in dictionary.
    type (Dictionary_type), intent(in) :: this
    getNumEntries = size(this%pairs)
  end function getNumEntries

  function lookup(this, key) result(values)
!@sum Returns GenericType corresponding to given key.
!@+ note that overloading of "=" in GenericType_mod
!@+ allows the result to be given to an appropriate
!@+ intrinsic type.  
    type (Dictionary_type), intent(in) :: this
    character(len=*), intent(in) :: key
    type (GenericType_type), pointer :: values(:)
    integer :: i

    i = getIndex(this, key)
    if (i /= NOT_FOUND) then
      allocate(values(getNumValues(this%pairs(i))))
      values = getValues(this%pairs(i))
    else
      ! need to allocate something to prevent a crash
      allocate(values(0))
      call throwException('Key not found: <'//trim(key)//'>.', 14)
    end if

  end function lookup

  subroutine cleanDictionary(this)
!@sum Restore data structure to pristine state
    type (Dictionary_type), intent(inout) :: this
    integer :: i
    do i = 1, size(this%pairs)
      call clean(this%pairs(i))
    end do
    deallocate(this%pairs)
  end subroutine cleanDictionary

  logical function hasKey(this, key)
!@sum Returns true iff dictionary has given key.
    type (Dictionary_type), intent(in) :: this
    character(len=*), intent(in) :: key

    integer :: i

    hasKey = .false.
    do i = 1, getNumEntries(this)
      if (trim(key) == getKey(this%pairs(i))) then
        hasKey = .true.
        return
      end if
    end do

  end function hasKey

  function getKeys_dictionary(this) result(keys)
!@sum Returns pointer aray of keys in dictioary.
    type (Dictionary_type), intent(in) :: this
    character(len=MAX_LEN_KEY), pointer :: keys(:)
    keys => getKeys(this%pairs)
  end function getKeys_dictionary

  integer function getIndex(this, key) result(index)
!@sum Returns index in array of pairs where provided key
!@+ can be found.   Returns parameter "NOT_FOUND" if
!@+ key is not present.  
!@+ NOTE:: THIS METHOD SHOULD NOT BE MADE PUBLIC
    type (Dictionary_type), intent(in) :: this
    character(len=*), intent(in) :: key

    integer :: i
    do i = 1, getNumEntries(this)
      if (toLowerCase(key) == toLowerCase(getKey(this%pairs(i)))) then
        index = i
        return
      end if
    end do
    index = NOT_FOUND
    
  end function getIndex

  subroutine readUnformatted_dictionary(this, unit)
!@sum Populates a dicitonary from unformatted sequential file.
!@+ NOTE: No header is used, so this procedure should always
!@+ be wrapped in logic that checks for a higher level header.
    type (Dictionary_type), intent(out) :: this
    integer, intent(in) :: unit
    
    type (KeyValuePair_type) :: pair
    integer :: i, n

    this = Dictionary()

    read(unit) n
    do i = 1, n
      call readUnformatted(pair, unit)
      call insert(this, pair)
    end do

  end subroutine readUnformatted_dictionary

  subroutine writeUnformatted_dictionary(this, unit)
!@sum Stores a dicitonary to an unformatted sequential file.
!@+ NOTE: No header is used, so this procedure should always
!@+ be wrapped in logic that provides such a safe header.

    type (Dictionary_type), intent(in) :: this
    integer, intent(in) :: unit

    integer :: i, n
    n = getNumEntries(this)
    write(unit) n
    do i = 1, n
      call writeUnformatted(this%pairs(i), unit)
    end do

  end subroutine writeUnformatted_dictionary

  logical function equals(a, b)
!@sum Returns .true. if two dictionaries are identical.
!@+ I.e. if both dictionaries have the same keys, and
!@+ for each key the corresponding values are the same.
!@+ Note that no assumption about ordering of keys is made.
    type (Dictionary_type), intent(in) :: a
    type (Dictionary_type), intent(in) :: b

    integer :: i, j
    integer :: numA, numB
    character(len=MAX_LEN_KEY), pointer :: keys(:)
    character(len=MAX_LEN_KEY) :: key

    numA = getNumEntries(a)
    numB = getNumEntries(b)
    equals = (numA == numB)

    if (equals) then
      keys => getKeys(a)
      do i = 1, numA
        key = trim(keys(i))

        equals = hasKey(b, key)
        if (.not. equals) exit

        j = getIndex(b, key)
        equals = (getNumValues(a%pairs(i)) == getNumValues(b%pairs(j)))
        if (.not. equals) exit

        equals = all(getValues(a%pairs(i)) == getValues(b%pairs(j)))
        if (.not. equals) exit
      end do
    end if

  end function equals

end module Dictionary_mod


!**** this should be put somewhere else, but since it is used only ****
!**** in this module I put it here for a while ...                 ****

subroutine swap_bytes_4( c, ndim )
  integer n,ndim
  character*1 c(4,ndim),temp
!@sum  does conversion big<->little - endian for 4 byte data
!@auth I. Aleinov
!@ver 1.0
  do n=1,ndim
    temp   = c(1,n)
    c(1,n) = c(4,n)
    c(4,n) = temp
    temp   = c(2,n)
    c(2,n) = c(3,n)
    c(3,n) = temp
  end do
end subroutine swap_bytes_4


