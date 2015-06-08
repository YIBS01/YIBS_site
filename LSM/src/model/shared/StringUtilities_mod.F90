module StringUtilities_mod
!@sum Container module for various procedures which manipulate strings
!@auth T.Clune
  implicit none
  private

  public :: toLowerCase

contains

  function toLowerCase(string) result(newString)
!@auth I. Aleinov
!@sum Produces a copy of the string argument, but with any
!@+ capitalized letters replaced by their lower case equivalent.
!@+ Characters which are not alphabetic letters are copied without
!@+  modification.
    character(len=*), intent(in) :: string
    character(len=len_trim(string)) :: newString

    integer n, i
    integer A, Z, shift, c

    A = iachar( 'A' )
    Z = iachar( 'Z' )
    shift = iachar( 'a' ) - iachar( 'A' )

    newString = trim(string)
    n = len(newString)
    do i=1,n
      c = iachar( newString(i:i) )
      if ( c>=A .and. c<=Z ) newString(i:i) = achar( c + shift )
    enddo

  end function toLowerCase

end module StringUtilities_mod
