
#ifdef MPI
  !Create an MPI-driven dataType from modelParam_t_mpi
  subroutine create_modelParam_t_mpi(sigma_mpi)

    implicit none
    include 'mpif.h'
    use MPI_declaration

    type (modelParam_t), intent(in)   	        :: sigma_mpi
    integer :: ny1,nz1,ii,Nza1,NzEarth,dimension1
    integer      block_lengths(0:20)
    integer      displacements(0:20)
    integer      address(0:21)
    integer      typelist(0:21),ierr


     Nz1 = sigma_mpi%grid%Nz
     Nza1 =sigma_mpi%grid%Nza
     NzEarth = Nz-Nza
     Ny1 = sigma_mpi%grid%Ny

       dimension1=Ny1*NzEarth


      typelist(0) = MPI_INTEGER
      typelist(1) = MPI_DOUBLE_PRECISION
      typelist(2) = MPI_DOUBLE_PRECISION
      typelist(3) = MPI_LOGICAL
      typelist(4) = MPI_CHARACTER

      block_lengths(0) = 2
      block_lengths(1) = dimension1
      block_lengths(2) = 1
      block_lengths(3) = 2
      block_lengths(4) = 80


      call MPI_Address(sigma_mpi,                         address(0), ierr)
      call MPI_Address(sigma_mpi%Ny,                      address(1), ierr)
      call MPI_Address(sigma_mpi%v(1,1),                  address(2), ierr)
      call MPI_Address(sigma_mpi%AirCond,                 address(3), ierr)
      call MPI_Address(sigma_mpi%allocated,               address(4), ierr)
      call MPI_Address(sigma_mpi%paramType,               address(5), ierr)

     do ii=0,4
       displacements(ii) = address(ii+1) - address(0)
     end do

      call MPI_TYPE_STRUCT(5, block_lengths, displacements,typelist, modelParam_t_mpi, ierr)
      call MPI_TYPE_COMMIT(modelParam_t_mpi, ierr)

  end subroutine create_modelParam_t_mpi
#endif
