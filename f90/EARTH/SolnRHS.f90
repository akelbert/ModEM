module solnrhs
!   higher level module to define EM solution and RHS objects
!   plus basic methods, linear algebra, dot products
!   EARTH version
!
! Defines: EMsoln, EMsparse, EMrhs
! Uses: EMfield

use math_constants
use utilities
use sg_vector
use sg_boundary
use sg_sparse_vector

implicit none

  type :: EMsolnMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer			:: nTx = 0
    type(cvector), pointer		:: solns(:)
    integer, pointer            :: tx(:)
    integer, pointer            :: errflag(:)
    type(grid_t), pointer       :: grid
    logical			:: allocated = .false.
  end type EMsolnMTX_t


Contains

!**********************************************************************
!           Basic EMsolnMTX methods
!**********************************************************************

   subroutine create_EMsolnMTX(nTx,eAll,grid)

      integer, intent(in)               :: nTx
      type(EMsolnMTX_t), intent(inout)  :: eAll
      type(grid_t), intent(in), target, optional :: grid

      !  local variables
      integer                           :: istat

      eAll%nTx = nTx
      allocate(eAll%solns(nTx), STAT=istat)
      allocate(eAll%tx(nTx), STAT=istat)
      allocate(eAll%errflag(nTx), STAT=istat)
      if(present(grid)) then
         eAll%grid => grid
      end if
      eAll%allocated = .true.

   end subroutine create_EMsolnMTX

   !**********************************************************************
   subroutine deall_EMsolnMTX(eAll)

      type(EMsolnMTX_t), intent(inout)  :: eAll

      !  local variables
      integer                           :: j, istat

	  do j = 1,eAll%nTx
	  	call deall_cvector(eAll%solns(j))
	  end do

      if (associated(eAll%solns)) deallocate(eAll%solns, STAT=istat)
      if (associated(eAll%tx)) deallocate(eAll%tx, STAT=istat)
      if (associated(eAll%errflag)) deallocate(eAll%errflag, STAT=istat)
      if (associated(eAll%grid)) nullify(eAll%grid)

      eAll%allocated = .false.

   end subroutine deall_EMsolnMTX

end module solnrhs
