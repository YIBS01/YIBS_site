module Parser_mod
!@sum procedures to read parameters from the rundeck into the database
!@auth I. Aleinov and T. Clune
!@ver 1.1     
  use GenericType_mod
  use KeyValuePair_mod
  use Dictionary_mod, only: MAX_LEN_LINE
  implicit none

  public :: Parser_type ! data type
  public :: parse
  public :: stripComment
  public :: skipHeader
  public :: isEndData
  public :: setCommentCharacters
  public :: setBeginData
  public :: setEndData
  public :: getValueType
  public :: splitTokens
  public :: writeFormatted

  public :: MAX_LEN_LINE
  public :: MAX_LEN_TOKEN
  integer, parameter :: MAX_COMMENT_CHARACTERS = 3
  integer, parameter :: MAX_TOKEN_SEPARATORS   = 2
  integer, parameter :: MAX_LEN_TOKEN = 32
  character(len=*), parameter :: ENTIRE_LINE = '(a257)' ! MAX_LEN_LINE + 1 char
  character(len=*), parameter :: ENTIRE_TOKEN = '(a33)' ! MAX_LEN_TOKEN + 1 char

  type Parser_type
    character(len=MAX_COMMENT_CHARACTERS) :: commentCharacters = '!#' ! legacy default
    character(len=MAX_TOKEN_SEPARATORS) :: tokenSeparators = ' =,'     ! legacy default
    character(len=MAX_LEN_LINE) :: beginData = '&&PARAMETERS'         ! legacy default
    character(len=MAX_LEN_LINE) :: endData = '&&END_PARAMETERS'     ! legacy default
  end type Parser_type

  type (Parser_type), save :: globalParser

  interface getValueType
    module procedure getValueType_single
    module procedure getValueType_multi
  end interface

  interface writeFormatted
    module procedure writeFormatted_parser
  end interface

contains

  function strip_comment( str ) result(newStr)
    ! remove comment at the end of the line. comments symbols: !#
    character*(*), intent(in) :: str
    character(len=len(str)) :: newStr
    
    newStr = stripComment(globalParser, str)

  end function strip_comment

  subroutine skip_junk( str )
    character*(*) str

    do while ( len_trim( str ) > 0 .and. scan( str, ' =,' ) == 1 )
      str = str(2:)
    enddo
    return
  end subroutine skip_junk

  subroutine sread_int( str, value )
    character*(*) str
    integer value
    integer n

    read ( str, * ) value

    ! remove chars till the next [ =,]
    n = scan( str, ' =,' )
    str = str(n+1:)

    call skip_junk( str )

    return
  end subroutine sread_int

  subroutine sread_real( str, value )
    character*(*) str
    real*8 value
    integer n

    read ( str, * ) value

    ! remove chars till the next [ =,]
    n = scan( str, ' =,' )
    str = str(n+1:)

    call skip_junk( str )

  end subroutine sread_real

  subroutine sread_char( str, value )
    character*(*) str
    character*(*) value
    !character*256 tstr
    integer n, n1

    ! replace '=' with space if not quoted
    n1 = scan( str, '''' )
    n  = scan( str, '=' )
    if ( n>0 .and. ( n1==0 .or. n<n1 ) ) str(n:n) = ' '

    read ( str, * ) value

    if ( scan( str, '''' ) == 1 ) then  ! quoted string
      str = str(2:)
      n = scan( str, '''' )
      str = str(n+1:)
    else  ! remove chars till the next [ =,]
      n = scan( str, ' =,' ) 
      str = str(n+1:)
    endif

    call skip_junk( str )

  end subroutine sread_char

  subroutine parse_params( kunit )
    use PARAM
    integer, parameter :: MAXDIM=128
    integer, intent(in) :: kunit
    character*256 bufs
    character*32 name
    character*1 type
    integer np
    integer ivars(MAXDIM)
    real*8 rvars(MAXDIM)
    character*64 cvars(MAXDIM)

    ! skip unrelated stuff
    do
      read( kunit, '(a256)', err=666, end=667 ) bufs
      if ( len_trim(bufs) < 1 ) cycle
      read( bufs, * ) name
      if ( name == '&&PARAMETERS' ) exit
    enddo

    do
      read( kunit, '(a256)', err=666, end=666 ) bufs

      if ( len_trim(bufs) < 1 ) cycle

      bufs = strip_comment( bufs )
      call skip_junk( bufs )

      if ( len_trim(bufs) < 1 ) cycle

      !read the name of the variable
      call sread_char( bufs, name )

      if ( name == '&&END_PARAMETERS' ) exit  ! end of list 

      if ( len_trim(bufs) < 1 ) then
        print *,'PARSER: no values were given to param: ', name
        call stop_model('PARSER error',255)
      endif

      ! now check the type of variables
      if ( scan( bufs, '''' ) > 0 ) then
        type = 'c'
      else if ( scan( bufs, '.' ) > 0 ) then
        type = 'r'
      else
        type = 'i'
      endif

      select case ( type )
      case ('i')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          if (np > MAXDIM) call stop_model("parse_params: increase MAXDIM",255)
          call sread_int( bufs, ivars(np) )
        end do
        call set_param( name, ivars, np, 'or' )
      case ('r')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          if (np > MAXDIM) call stop_model("parse_params: increase MAXDIM",255)
          call sread_real( bufs, rvars(np) )
        end do
        call set_param( name, rvars, np, 'or' )
      case ('c')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          if (np > MAXDIM) call stop_model("parse_params: increase MAXDIM",255)
          call sread_char( bufs, cvars(np) )
        end do
        call set_param( name, cvars, np, 'or' )
      end select

    enddo

    return
666 print *, 'PARSER: Error reading params'
    call stop_model( 'PARSER: Error reading params', 255 )
667 print *, 'PARSER: No &&PARAMETERS or &&END_PARAMETERS found'
    call stop_model( &
    &     'PARSER: No &&PARAMETERS or &&END_PARAMETERS found',255)
  end subroutine parse_params

  function parseLine(this, line) result(pair)
!@sum Attempts to create a key value pair by parsing a line of text
!@+ into separate tokens (via splitTokens).  Throws an exception if
!@+ less than two tokens are found.   (Should probably be modified
!@+ to also throw exception if first token is not an exceptable
!@+ key.)
    use Dictionary_mod
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: line
    type (KeyValuePair_type) :: pair

    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    tokens => splitTokens(this, line)
    if (size(tokens) < 2) then
      call throwException('Parser_mod: syntax error in input unit.', 14)
      return
    end if

    pair = readPair(key=tokens(1), tokens=tokens(2:))
    deallocate(tokens)

  contains

    function readPair(key, tokens) result(pair)
!@sum Helper function.  Creates a KeyValuePair from
!@+ a key and list of tokens.
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: tokens(:)
      type (KeyValuePair_type) :: pair

      type (GenericType_type), pointer :: values(:)
      integer :: numValues, i

      numValues = size(tokens)
      allocate(values(numValues))
      do i = 1, numValues
        values(i) = GenericType(tokens(i), getValueType(tokens))
      end do
      pair = KeyValuePair(key, values)
      deallocate(values)

    end function readPair

  end function parseLine

  function parse(this, unit, status) result(aDictionary)
!@sum Parses text from input unit to populate a Dictionary object.
!@+ Skips header section at top.
    use Dictionary_mod

    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: status

    type (Dictionary_type) :: aDictionary
    character(len=MAX_LEN_LINE) :: line
    type (KeyValuePair_type) :: pair

    status = 0 ! unless ...
    aDictionary = Dictionary()

    call skipHeader(this, unit, status)
    if (status /= 0) return
    
    do
      read(unit,fmt=ENTIRE_LINE,iostat=status) line
      if (status /= 0) exit
      if (isEndData(this, line)) exit
      
      line = stripComment(this, line)
      if (len_trim(line) == 0) cycle ! skip

      pair = parseLine(this, line)
      call insert(aDictionary, pair)

    end do

  end function parse
  
  subroutine setCommentCharacters(this, commentCharacters)
!@sum Set the characters which should be interpreted as comment
!@+ when processing input.  Defaults are conventional Fortran "!#",
!@+ but this routine can be used to override.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: commentCharacters

    this%commentCharacters = commentCharacters

  end subroutine setCommentCharacters

  subroutine setTokenSeparators(this, tokenSeparators)
!@sum Override default characters used to separate tokens.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: tokenSeparators
    this%tokenSeparators = tokenSeparators
  end subroutine setTokenSeparators

  function stripComment(this, str) result(newStr)
!@sum Eliminate trailing comments in line of text. 
!TO DO - ignore comments in strings
    type (Parser_type), intent(in) :: this
    character*(*), intent(in) :: str
    character(len=len(str)) :: newStr

    integer :: n

    n = scan(str, trim(this%commentCharacters))
    select case (n)
    case (0)
      newStr = trim(str)
    case (1:)
      newStr = str(:n-1)
    end select

  end function stripComment

  subroutine setBeginData(this, beginData)
!@sum Override default string signfying beginning of data,
!@+ i.e. end of header.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: beginData
    this%beginData = beginData
  end subroutine setBeginData

  logical function isEndData(this, string)
!@sum Return true if string matches the semaphore for the end of data.
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string
    
    isEndData = (trim(this%endData) == trim(string))

  end function isEndData

  subroutine setEndData(this, endData)
!@sum Override default string signfying end of data.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: endData
    this%endData = endData
  end subroutine setEndData

  subroutine skipHeader(this, unit, status)
!@sum Read from unit until semaphore for beginning of data
!@+ is located.  Returns nonzero status if semaphore is not
!@+ found before EOF.
    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: status

    character(len=MAX_LEN_LINE) :: line
    
    do
      read(unit,fmt=ENTIRE_LINE,iostat=status) line
      if (status < 0) then
        return
      end if

      if (trim(line) == this%beginData) exit
    end do

  end subroutine skipHeader

  logical function isInteger(string)
!@sum Returns true if string can be converted to an integer.
!@+ Should possibly relocate to GenericType_mod
    character(len=*), intent(in) :: string
    integer :: integerValue
    integer :: status
    
    read(string,'(i20)',iostat=status) integerValue
    isInteger = (status == 0)

  end function isInteger

  logical function isReal64(string)
!@sum Returns true if string can be converted to a double.
!@+ Should possibly relocate to GenericType_mod
    character(len=*), intent(in) :: string
    real*8 :: real64Value
    integer :: status
    
    read(string,'(g30.30)',iostat=status) real64Value
    isReal64 = (status == 0)

  end function isReal64

  logical function isLogical(string)
!@sum Returns true if string can be converted to a logical.
!@+ Allows for Fortran defaults AND other strings that are 
!@+ convenient equivalents.  (See implementation below.)
!@+ Should possibly relocate to GenericType_mod
    use StringUtilities_mod, only: toLowerCase
    character(len=*), intent(in) :: string
    logical :: logicalValue
    integer :: status
    
    select case (trim(toLowerCase(string)))
    case ('t','f','true','false','.true.','.false.')
      isLogical = .true.
    case default
      isLogical = .false.
    end select

  end function isLogical

  integer function getValueType_single(string) result(type)
!@sum Returns integer which corresponds to the _type_ that this string
!@+ can be converted to.  (INTEGER_TYPE, LOGICAL_TYPE, etc.)
!@+ Priority is given to more restrictive types, which
!@+ string being a catchall.
    use Dictionary_mod
    character(len=*), intent(in) :: string

    real*8  :: real64Value
    logical :: logicalValue
    character(len=MAX_LEN_VALUE)  :: stringValue
    integer :: status
  
    if (isInteger(string)) then
      type = INTEGER_TYPE
      return
    else if (isReal64(string)) then
      type = REAL64_TYPE
      return
    else if (isLogical(string)) then
      type = LOGICAL_TYPE
      return
    else
      type = STRING_TYPE
    end if

  end function getValueType_single

  integer function getValueType_multi(tokens) result(type)
!@sum Return the _type_ that a set of tokens can be converted to.  @+
!@+ Note that type is the lowest common denominator. E.g. mixtures of
!@+ integers and doubles will be processed as double values, and 
!@+ if any token must be treated as a string, all will be.
    use Dictionary_mod
    character(len=*), intent(in) :: tokens(:)

    integer :: numTokens
    integer :: i
    integer, allocatable :: types(:)

    numTokens = size(tokens)
    allocate(types(numTokens))
    do i = 1, numTokens
      types(i) = getValueType(tokens(i))
    end do

    if (all(types == types(1))) then ! simple case
      type = types(1)
    else ! mixed type
      if (all((types == INTEGER_TYPE) .or. (types == REAL64_TYPE))) then
        type = REAL64_TYPE ! integers treated as subset of reals
      else
        type = STRING_TYPE ! most general category
      end if
    end if

    deallocate(types)

  end function getValueType_multi

  function splitTokens(this, string) result(tokens)
!@ Attempt to split a string into a set of independent tokens
!@+ based upon separaters specified in the Parser object.
!@+ A crude attempt is made to ensure that the 1st separator is '='.
!@+ A more comprehensive treatement of this should be made.
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    character(len=len(string)) :: buffer
    integer :: numTokens
    integer :: i, idxStart
    integer :: idxNextSeparator

    numTokens = countTokens(this, string)
    allocate(tokens(numTokens))
    
    buffer = adjustl(string)
    i = 0
    do  ! pull off tokens until
      if (len_trim(buffer) == 0) exit ! done
      i = i + 1
      idxStart = skipEmbeddedSeparators(buffer)
      idxNextSeparator = scan(buffer(idxStart:), trim(this%tokenSeparators))

      if (idxNextSeparator > 0) then
        tokens(i) = trim(buffer(:idxStart+idxNextSeparator-2))
        buffer = adjustl(buffer(idxStart+idxNextSeparator:))
      else
        tokens(i) = trim(buffer) ! take the rest
        exit
      end if
    end do

    if (numTokens >= 1) then
      if (scan(trim(tokens(1)), ' ') /= 0) then
        call throwException('Parser_mod: Illegal syntax.  "=" not first separator.', 14)
      end if
    end if

  end function splitTokens

  integer function countTokens(this, string)
!@sum Sweep through string to count tokens.  This can then
!@+ be used to allocate the token array and a 2nd sweep
!@+ will fill.
!TO DO - logic with splitTokens() is duplicated. 
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string

    character(len=len(string)) :: buffer
    integer :: numTokens
    integer :: idxNextSeparator, idxStart

    buffer = adjustl(string)
    numTokens = 0
    do
      if (len_trim(buffer) == 0) exit ! done
      numTokens = numTokens + 1
      ! skip strings which might have embedded separator characters
      idxStart = skipEmbeddedSeparators(buffer)
      idxNextSeparator = scan(buffer(idxStart:), this%tokenSeparators)
      if (idxNextSeparator > 0) then
        buffer = adjustl(buffer(idxStart + idxNextSeparator:))
      else
        exit
      end if
    end do

    countTokens = numTokens

  end function countTokens

  integer function skipEmbeddedSeparators(buffer) result(idxStart)
    character(len=*), intent(in) :: buffer
!@sum Location of end of string, which might have an embedded separator.
    if (scan(buffer, '"') == 1) then
      idxStart = scan(buffer(2:),'"') + 1 + 1 ! include 1st char
    else if (scan(buffer, "'") == 1) then
      idxStart = scan(buffer(2:),"'") + 1 + 1 ! include 1st char
    else
      idxStart = 1
    end if
  end function skipEmbeddedSeparators
    
  subroutine writeFormatted_parser(this, unit, aDictionary)
!@sum Write dictionary as a text file.  Inverse of
!@+ parse().
    use Dictionary_mod
    type (Parser_type), intent(in) :: this
    type (Dictionary_type), intent(in) :: aDictionary
    integer, intent(in) :: unit

    integer :: i, j
    type (GenericType_type), pointer :: values(:)
    character(len=MAX_LEN_LINE) :: line
    character(len=MAX_LEN_KEY), pointer :: keys(:)

    keys => getKeys(aDictionary)
    write(unit,*) trim(this%beginData)
    do i = 1, getNumEntries(aDictionary)

      values =>lookup(aDictionary, keys(i))

      select case (size(values))
      case (1)
        write(line,'(2x,a32," = ",a)') keys(i), trim(toString(values(1)))
      case (2:)
        write(line,'(2x,a32," = ",a)') keys(i), trim(toString(values(1)))
        do j = 2, size(values)
          line = trim(line) // ', ' // trim(toString(values(j)))
        end do
      end select
      write(unit,*) trim(line)
    end do
    write(unit,*) trim(this%endData)

  end subroutine writeFormatted_parser

end module PARSER_MOD
