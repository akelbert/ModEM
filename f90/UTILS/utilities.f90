module utilities

	use math_constants
	implicit none

	!character(80)  :: msg

Contains

  !*****************************************************************************
  subroutine errStop(msg)

    character(*), intent(in)  :: msg
    write(0,'(a9)',advance='no') 'Error: '
    write(0,*) trim(msg)
    stop

  end subroutine errStop

  !*****************************************************************************
  subroutine warning(msg)

    character(*), intent(in)  :: msg
    write(0,'(a9)',advance='no') 'Warning: '
    write(0,*) trim(msg)

  end subroutine warning

  ! **************************************************************************
  function clean(x)
    ! This is a utility routine that provides an expression used to battle
	! against machine error problems. It returns the same real or real(8)
	! as the input, but without the extra digits at the end that are often
	! a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
	! instead of x in an inequality!!!
	! LARGE_REAL is defined in the module math_constants
	! A.K.
    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec)						       :: clean

	clean = dnint(x*LARGE_REAL)/LARGE_REAL

  end function clean

  ! **************************************************************************
  function nearest_meter(x) result(clean)
    ! This is a utility routine that provides an expression used to battle
	! against machine error problems. Both input and output are values in km.
	! The function rounds the value to the nearest meter. This is useful to
	! ensure that the grid read from a file does not depend on system precision.
	! A.K.
    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec)						       :: clean

	clean = dnint(x*KM2M)/KM2M

  end function nearest_meter


  ! **************************************************************************
  function minNode(x, xNode) result(ix)
    !  This is a utility routine, used by several data functional
    !  set up routines, and for other interpolation functions
    !  Returns index ix such that  xNode(ix) <= x < xNode(ix+1)
    !  If x is out of range:
    !  x < xNode(1) returns 0; if x> xNode(nx) returns nx
    !  Assumes xNode is strictly increasing; does not check this
    !  NOTE: as presently coded, when xNode is called with center
    !  (face) node positions, this routine will return zero for
    !  the coordinates in the outer half cell nearest the boundary
    !  If evaluation over the complete model domain is to be allowed
    !  a more general interpolation rule will be required.
    !  A.K.: modified to allow input of any size, nx = size(xNode).

    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec), dimension(:), intent(in)     :: xNode

    integer                                     :: ix
    integer                                     :: i

    ix = size(xNode)
    do i = 1,size(xNode)
       if(clean(xNode(i)) .gt. clean(x)) then
          ix = i-1
          exit
       endif
    enddo

  end function minNode


  ! **************************************************************************
  function maxNode(x, xNode) result(ix)
    !  This is a utility routine, used by several data functional
    !  set up routines, and for other interpolation functions
    !  Returns index ix such that  xNode(ix) <= x < xNode(ix+1)
    !  If x is out of range:
    !  x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !  Assumes xNode is strictly decreasing; does not check this
    !  NOTE: as presently coded, when xNode is called with center
    !  (face) node positions, this routine will return zero for
    !  the coordinates in the outer half cell nearest the boundary
    !  If evaluation over the complete model domain is to be allowed
    !  a more general interpolation rule will be required.
    !  A.K.: modified to allow input of any size, nx = size(xNode).

    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec), dimension(:), intent(in)     :: xNode

    integer                                     :: ix
    integer                                     :: i

    ix = size(xNode)
    do i = 1,size(xNode)
       if(clean(xNode(i)) .lt. clean(x)) then
          ix = i-1
          exit
       endif
    enddo

  end function maxNode

! *****************************************************************************

      integer function findstr(str1,str2)
      character*(*) str1, str2
!     returns the position of str2 in str1.  Ignores case.
!     returns 0 if str2 not found in str1

      integer i, j, capdif
      logical same

      capdif= ichar('a')-ichar('A')

      do 20 i= 1, len(str1)-len(str2)+1
         do 10 j=1,len(str2)

	      same= str1(i+j-1:i+j-1) .eq. str2(j:j) .or.  &
            'A'.le.str2(j:j) .and. str2(j:j).le.'Z' .and.  &
	     ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j))+capdif .or.  &
            'a'.le.str2(j:j) .and. str2(j:j).le.'z' .and.  &
	     ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j)) - capdif

	     if( .not.same) go to 20
10       continue
         findstr=i
         return
20    continue

      findstr=0
      return
      end function findstr



! *****************************************************************************

      integer function begwrd(string,iwrd)
      integer iwrd
      character*(*) string

!     Returns the index of the first non-blank character in the iwrd'th
!     non-blank word (word are seperated by spaces, tabs or commas).
!     Returns len if iwrd'th word is not found. integer i, nword

      logical wasblk
      intrinsic len
      integer  i,nword

      wasblk=.true.

      nword= 0
      do i=1,len(string)
         if( string(i:i).eq.' ' .or.string(i:i).eq.',' .or.  &
	     string(i:i).eq.'  '    )then

!           /* current character is blank */
             wasblk=.true.
	 else
	     if (wasblk) then
		nword= nword + 1
	     endif
	     wasblk= .false.
	     if(nword.eq.iwrd)then
	        begwrd= i
		return
	     end if
	 end if
      enddo

      begwrd = len(string)
      return

      end function begwrd

! *****************************************************************************


      integer function endwrd(string,iwrd)
      integer iwrd
      character*(*) string
!     Returns the index of the last non-blank character in the iwrd'th
!     non-blank word (word are seperated by spaces, tabs or commas).
!     Returns len if iwrd'th word is not found.
      integer i, nword
      logical wasblk
      intrinsic len

      wasblk=.true.
      nword= 0
      do 100 i=1,len(string)
	if( string(i:i).eq.' ' .or.  &
	    string(i:i).eq.',' .or.  &
	    string(i:i).eq.'  '    )then

!          /* current character is blank */
           wasblk=.true.
           if(nword.eq.iwrd) RETURN

	else
           if(wasblk) nword= nword + 1
           wasblk= .false.
           if(nword.eq.iwrd) endwrd= i
	end if
100   continue

      endwrd= len(string)

      return

      end function endwrd


! *****************************************************************************

      SUBROUTINE Lenb(string,length)
      IMPLICIT NONE
      CHARACTER*(*) string
      INTEGER nstr,istr,length

      nstr = len(string)
      DO istr=nstr,1,-1
         IF (string(istr:istr).ne.' ') THEN
	    length = istr
            RETURN
	 ENDIF
      ENDDO
      length = 0


      RETURN

      END Subroutine lenb

  ! **************************************************************************
  ! Naser Meqbel included this function: apparently, it is not supported by
  ! all compilers as an intrinsic
  logical function isnan(a)

        real (kind=prec), intent(in) ::a

        if (a .ne. a) then
        	isnan = .true.
        else
        	isnan = .false.
        end if

  end function isnan

end module utilities
