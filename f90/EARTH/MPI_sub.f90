
module MPI_sub
#ifdef MPI

use math_constants
use utilities
use ForwardSolver
use MPI_declaration






Contains

subroutine create_eAll_mpi(eAll)
use math_constants
use utilities
use ForwardSolver



    implicit none



      type(solnVectorMTX_t), intent(in)	:: eAll
      integer :: nx1,ny1,nz1,ii

     write(6,*) taskid,eAll%solns(1)%grid%nx


      nx1=eAll%solns(which_per)%grid%nx
      ny1=eAll%solns(which_per)%grid%ny
      nz1=eAll%solns(which_per)%grid%nz
     write(6,*) taskid, nx1,ny1,nz1

      typelist(0) = MPI_INTEGER
      typelist(1) = MPI_CHARACTER
      typelist(2) = MPI_DOUBLE_COMPLEX
      typelist(3) = MPI_DOUBLE_COMPLEX
      typelist(4) = MPI_DOUBLE_COMPLEX
      typelist(5) = MPI_INTEGER
      typelist(6) = MPI_LOGICAL
      typelist(7) = MPI_LOGICAL
      typelist(8) = MPI_INTEGER

      block_lengths(0) =1
      block_lengths(1) = 80
      block_lengths(2) = (nx1*(ny1+1)*(nz1+1))
      block_lengths(3) = ((nx1+1)*ny1*(nz1+1))
      block_lengths(4) = ((nx1+1)*(ny1+1)*nz1)
      block_lengths(5) = 3
      block_lengths(6) = 1
      block_lengths(7) = 1
      block_lengths(8) = 1

      call MPI_Address(eAll%solns(which_per),                          address(0), ierr)
      call MPI_Address(eAll%solns(which_per)%errflag,                     address(1), ierr)
      call MPI_Address(eAll%solns(which_per)%vec%gridType,          address(2), ierr)
      call MPI_Address(eAll%solns(which_per)%vec%x(1,1,1),          address(3), ierr)
      call MPI_Address(eAll%solns(which_per)%vec%y(1,1,1),          address(4), ierr)
      call MPI_Address(eAll%solns(which_per)%vec%z(1,1,1),          address(5), ierr)
      call MPI_Address(eAll%solns(which_per)%vec%nx,                address(6),ierr)
      call MPI_Address(eAll%solns(which_per)%vec%allocated,         address(7),ierr)
      call MPI_Address(eAll%solns(which_per)%vec%temporary,         address(8),ierr)
      call MPI_Address(eAll%solns(which_per)%tx,                       address(9),ierr)

     do ii=0,8
       displacements(ii) = address(ii+1) - address(0)
     end do

      call MPI_TYPE_STRUCT(9, block_lengths, displacements,typelist, eAll_mpi, ierr)
      call MPI_TYPE_COMMIT(eAll_mpi, ierr)
      write(6,*) taskid,eAll_mpi





end subroutine create_eAll_mpi
#endif

end module MPI_sub
