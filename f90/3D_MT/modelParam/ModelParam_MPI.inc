
#ifdef MPI
  !Create an MPI-driven dataType from modelParam_t_mpi
  subroutine create_modelParam_t_mpi(sigma_mpi)

    implicit none
    include 'mpif.h'

    type (modelParam_t), intent(in)   	        :: sigma_mpi
    integer :: nx1,ny1,nz1,ii,Nza1,NzEarth,dimension1
            integer      block_lengths(0:20)
      integer      displacements(0:20)
      integer      address(0:21)
      integer      typelist(0:21),ierr

      nx1=sigma_mpi%grid%nx
      ny1=sigma_mpi%grid%ny
      nz1=sigma_mpi%grid%Nz
      nzEarth=sigma_mpi%grid%nzEarth


    if (sigma_mpi%cellCond%gridType == CENTER) then
       dimension1=(nx1*ny1*nz1)
    else if (sigma_mpi%cellCond%gridType == CORNER) then
      dimension1=((nx1+1)*(ny1+1)*(nz1+1))
    else if(sigma_mpi%cellCond%gridType == CELL_EARTH) then
       dimension1=(nx1*ny1*nzEarth)
    end if


      typelist(0) = MPI_INTEGER
      typelist(1) = MPI_CHARACTER
      typelist(2) = MPI_DOUBLE_PRECISION
      typelist(3) = MPI_INTEGER
      typelist(4) = MPI_LOGICAL
      typelist(5) = MPI_DOUBLE_PRECISION
      typelist(6) = MPI_LOGICAL
      typelist(7) = MPI_CHARACTER
      typelist(8) = MPI_LOGICAL


      block_lengths(0) = 3
      block_lengths(1) = 80
      block_lengths(2) = dimension1
      block_lengths(3) = 3
      block_lengths(4) = 1
      block_lengths(5) = 1
      block_lengths(6) = 1
      block_lengths(7) = 80
      block_lengths(8) = 1

      call MPI_Address(sigma_mpi%Nx,                      address(0), ierr)
      call MPI_Address(sigma_mpi%cellCond,                address(1), ierr)
      call MPI_Address(sigma_mpi%cellCond%gridType,       address(2), ierr)
      call MPI_Address(sigma_mpi%cellCond%v(1,1,1),       address(3), ierr)
      call MPI_Address(sigma_mpi%cellCond%nx,             address(4), ierr)
      call MPI_Address(sigma_mpi%cellCond%allocated,      address(5), ierr)
      call MPI_Address(sigma_mpi%AirCond,                 address(6), ierr)
      call MPI_Address(sigma_mpi%allocated,               address(7), ierr)
      call MPI_Address(sigma_mpi%paramType,               address(8), ierr)
      call MPI_Address(sigma_mpi%temporary,               address(9), ierr)

     do ii=0,8
       displacements(ii) = address(ii+1) - address(0)
     end do

      call MPI_TYPE_STRUCT(9, block_lengths, displacements,typelist, modelParam_t_mpi, ierr)
      call MPI_TYPE_COMMIT(modelParam_t_mpi, ierr)

   end subroutine create_modelParam_t_mpi
#endif
