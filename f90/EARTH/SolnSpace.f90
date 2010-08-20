module SolnSpace
!   higher level module to define EM solution and RHS objects
!   plus basic methods, linear algebra, dot products
!   EARTH version
!
! Defines: solnVector, sparseVector, rhsVector
! Uses: EMfield

use math_constants
use file_units
use utilities
use sg_vector
use sg_boundary
use sg_sparse_vector

implicit none

  type :: solnVectorMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer			:: nTx = 0
    type(cvector), pointer		:: solns(:)
    integer, pointer            :: tx(:)
    integer, pointer            :: errflag(:)
    type(grid_t), pointer       :: grid
    logical			:: allocated = .false.
  end type solnVectorMTX_t


Contains

!**********************************************************************
!           Basic solnVectorMTX methods
!**********************************************************************

   subroutine create_solnVectorMTX(nTx,eAll,grid)

      integer, intent(in)               :: nTx
      type(solnVectorMTX_t), intent(inout)  :: eAll
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

   end subroutine create_solnVectorMTX

   !**********************************************************************
   subroutine deall_solnVectorMTX(eAll)

      type(solnVectorMTX_t), intent(inout)  :: eAll

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

   end subroutine deall_solnVectorMTX

  !****************************************************************************
  ! linComb_solnVectorMTX computes linear combination of two field solutions
  ! stored for multiple transmitters; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_solnVectorMTX(c1, E1, c2, E2, E3)

    implicit none
    !   input vectors
    type (solnVectorMTX_t), intent(in)         :: E1, E2
    !  input complex scalars
    complex (kind=prec), intent(in)        :: c1, c2
    type (solnVectorMTX_t), intent(inout)      :: E3
    ! local
    integer                                :: j

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
    	call errStop('inputs not allocated yet for linComb_solnVectorMTX')
    elseif (E1%nTx .ne. E2%nTx) then
    	call errStop('inputs of different sizes for linComb_solnVectorMTX')
    end if

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
    	call create_solnVectorMTX(E1%nTx,E3,E2%grid)
    end if

    do j=1,E1%nTx
    	if (E1%tx(j) == E2%tx(j)) then
    		! initialize...
        	E3%solns(j) = E1%solns(j)
        	! form linear combination
			call linComb_cvector(c1, E1%solns(j), c2, E2%solns(j), E3%solns(j))
		else
			call errStop('inputs for different transmitters in linComb_solnVectorMTX')
		end if
    end do

    E3%tx = E1%tx
    E3%errflag = max(E1%errflag, E2%errflag)
    E3%grid => E1%grid

  end subroutine linComb_solnVectorMTX ! linComb_solnVectorMTX

  !****************************************************************************
  ! write_solnVectorMTX writes an ASCII data file containing the full solnVectorMTX
  subroutine write_solnVectorMTX(fname, E)

    implicit none
    !   input vectors
    character(*), intent(in)			   :: fname
    type (solnVectorMTX_t), intent(in)         :: E
    ! local
    integer                                :: j,ios

    if(.not.E%allocated) then
    	call errStop('input not allocated yet for write_solnVectorMTX')
    end if

	open(ioWRITE,file=fname,status='unknown',form='formatted',iostat=ios)
	write(ioWRITE,'(a36,i3,a8)') "# Full EM field solution output for ",E%nTx,"periods."
	do j = 1,E%nTx
		write(ioWRITE,'(i3)') E%tx(j)
		call write_cvector(ioWRITE,E%solns(j))
	end do
	close(ioWRITE)

  end subroutine write_solnVectorMTX ! write_solnVectorMTX

  !****************************************************************************
  ! read_solnVectorMTX reads an ASCII data file containing the full solnVectorMTX
  subroutine read_solnVectorMTX(fname, E, grid)

    implicit none
    !   input vectors
    character(*), intent(in)			   :: fname
    type (solnVectorMTX_t), intent(inout)      :: E
    type (grid_t), intent(in)			   :: grid
    ! local
    integer                                :: j,nTx,ios,istat
    character(100)						   :: comment

    if(E%allocated) then
    	call deall_solnVectorMTX(E)
    end if

	open(ioREAD,file=fname,status='unknown',form='formatted',iostat=ios)
	read(ioREAD,'(a35)',iostat=istat,advance='no') comment
	read(ioREAD,*,iostat=istat) nTx
	print *, 'Number of transmitters: ',nTx
    call create_solnVectorMTX(nTx,E,grid)
	do j = 1,nTx
		read(ioREAD,'(i3)',iostat=istat) E%tx(j)
		call read_cvector(ioREAD,E%solns(j),grid)
		E%errflag(j) = 0
	end do
	close(ioREAD)

  end subroutine read_solnVectorMTX ! read_solnVectorMTX

end module SolnSpace
