
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
      typelist(5) = MPI_DOUBLE_COMPLEX
      typelist(6) = MPI_DOUBLE_COMPLEX
      typelist(7) = MPI_DOUBLE_COMPLEX
      typelist(8) = MPI_INTEGER
      typelist(9) = MPI_LOGICAL
      typelist(10) = MPI_LOGICAL
      typelist(11) = MPI_INTEGER
       
      block_lengths(0) =1
      block_lengths(1) = 80
      block_lengths(2) = (nx1*(ny1+1)*(nz1+1))
      block_lengths(3) = ((nx1+1)*ny1*(nz1+1))
      block_lengths(4) = ((nx1+1)*(ny1+1)*nz1)
      block_lengths(5) = (nx1*(ny1+1)*(nz1+1))
      block_lengths(6) = ((nx1+1)*ny1*(nz1+1))
      block_lengths(7) = ((nx1+1)*(ny1+1)*nz1)
      block_lengths(8) = 3
      block_lengths(9) = 1
      block_lengths(10) = 1
      block_lengths(11) = 1
      
      call MPI_Address(eAll%solns(which_per),                          address(0), ierr)
      call MPI_Address(eAll%solns(which_per)%nPol,                     address(1), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%gridType,          address(2), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%x(1,1,1),          address(3), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%y(1,1,1),          address(4), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%z(1,1,1),          address(5), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(2)%x(1,1,1),          address(6), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(2)%y(1,1,1),          address(7), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(2)%z(1,1,1),          address(8), ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%nx,                address(9),ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%allocated,         address(10),ierr)
      call MPI_Address(eAll%solns(which_per)%pol(1)%temporary,         address(11),ierr)
      call MPI_Address(eAll%solns(which_per)%tx,                       address(12),ierr)
      
     do ii=0,11
       displacements(ii) = address(ii+1) - address(0)
     end do   
     
      call MPI_TYPE_STRUCT(12, block_lengths, displacements,typelist, eAll_mpi, ierr)   
      call MPI_TYPE_COMMIT(eAll_mpi, ierr)  
      write(6,*) taskid,eAll_mpi


 


end subroutine create_eAll_mpi

!*************************************************************************
!Creating an MPI version of userdef_control_MPI 
subroutine create_userdef_control_MPI(ctrl) 
    implicit none
    include 'mpif.h'
    
	type(userdef_control), intent(in)   :: ctrl
	integer :: ii,intex,CHARACTERex,sum1
	integer(MPI_ADDRESS_KIND) address1(0:20)

	
   offsets(0) = 0
   oldtypes(0) = MPI_CHARACTER
   blockcounts(0) =  80*15

   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(1) = (80*15 * extent)+offsets(0)
   oldtypes(1) = MPI_REAL8
   blockcounts(1) = 2
   		
   call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
   offsets(2) = (2* extent)+offsets(1)
   oldtypes(2) = MPI_CHARACTER
   blockcounts(2) = 80*2
   
   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(3) = (80*2* extent)+offsets(2)
   oldtypes(3) = MPI_INTEGER
   blockcounts(3) = 1
   
     call MPI_TYPE_STRUCT(4, blockcounts, offsets, oldtypes,userdef_control_MPI, ierr)
     call MPI_TYPE_COMMIT(userdef_control_MPI, ierr)

end  subroutine create_userdef_control_MPI
#endif

end module MPI_sub
