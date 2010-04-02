
#ifdef MPI
  !Create an MPI dataType from modelParam_t_mpi
  subroutine create_modelParam_t_mpi(sigma_mpi)

    implicit none
    include 'mpif.h'

    type (modelParam_t), intent(in)   	        :: sigma_mpi
    integer :: nx1,ny1,nz1,ii,Nza1,NzEarth,dimension1
            integer      block_lengths(0:20)
      integer      displacements(0:20)
      integer      address(0:21)
      integer      typelist(0:21),ierr,x,y,z,counter

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
      typelist(9) = MPI_LOGICAL


      block_lengths(0) = 3
      block_lengths(1) = 80
      block_lengths(2) = dimension1
      block_lengths(3) = 3
      block_lengths(4) = 1
      block_lengths(5) = 1
      block_lengths(6) = 1
      block_lengths(7) = 80
      block_lengths(8) = 1
      block_lengths(9) = 1
      
      
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
      call MPI_Address(sigma_mpi%updated,                address(10), ierr)

     do ii=0,9
       displacements(ii) = address(ii+1) - address(0)
     end do

      call MPI_TYPE_STRUCT(10, block_lengths, displacements,typelist, modelParam_t_mpi, ierr)
      call MPI_TYPE_COMMIT(modelParam_t_mpi, ierr)

   end subroutine create_modelParam_t_mpi
   
  subroutine create_model_param_place_holder(sigma_mpi)
    implicit none
    include 'mpif.h'

      type (modelParam_t), intent(in)   	        :: sigma_mpi
      integer     :: nx1,ny1,nz1,ii,Nza1,NzEarth,dimension1
      integer     ::x,y,z,counter

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
         if(allocated(model_para_vec)) then
             deallocate(model_para_vec)
         end if
       allocate(model_para_vec(dimension1))

  
    

  
  end subroutine create_model_param_place_holder 
      
   subroutine get_model_para_values(sigma_mpi)
    implicit none
    include 'mpif.h'

      type (modelParam_t), intent(in)   	        :: sigma_mpi
      integer     :: nx1,ny1,nz1,ii,Nza1,NzEarth,dimension1
      integer     ::x,y,z,counter
 
    if (sigma_mpi%cellCond%gridType == CENTER) then
      nx1=sigma_mpi%grid%nx
      ny1=sigma_mpi%grid%ny
      nz1=sigma_mpi%grid%Nz
      nzEarth=sigma_mpi%grid%nzEarth
    else if (sigma_mpi%cellCond%gridType == CORNER) then
      nx1=sigma_mpi%grid%nx+1
      ny1=sigma_mpi%grid%ny+1
      nz1=sigma_mpi%grid%Nz+1
    else if(sigma_mpi%cellCond%gridType == CELL_EARTH) then
      nx1=sigma_mpi%grid%nx
      ny1=sigma_mpi%grid%ny
      nz1=sigma_mpi%grid%nzEarth
    end if   
    
      
       counter=0
       do x=1,nx1
         do y=1,ny1
           do z=1,nz1
               counter=counter+1 
                   model_para_vec(counter)=sigma_mpi%cellCond%v(x,y,z) 
           end do
         end do
      end do
      end subroutine get_model_para_values 
            
   subroutine set_model_para_values(sigma_mpi)
    implicit none
    include 'mpif.h'

      type (modelParam_t), intent(inout)   	        :: sigma_mpi
      integer     :: nx1,ny1,nz1,ii,Nza1,NzEarth,dimension1
      integer     ::x,y,z,counter
    if (sigma_mpi%cellCond%gridType == CENTER) then
      nx1=sigma_mpi%grid%nx
      ny1=sigma_mpi%grid%ny
      nz1=sigma_mpi%grid%Nz
      nzEarth=sigma_mpi%grid%nzEarth
    else if (sigma_mpi%cellCond%gridType == CORNER) then
      nx1=sigma_mpi%grid%nx+1
      ny1=sigma_mpi%grid%ny+1
      nz1=sigma_mpi%grid%Nz+1
    else if(sigma_mpi%cellCond%gridType == CELL_EARTH) then
      nx1=sigma_mpi%grid%nx
      ny1=sigma_mpi%grid%ny
      nz1=sigma_mpi%grid%nzEarth
    end if 
          
       counter=0
       do x=1,nx1
         do y=1,ny1
           do z=1,nz1
               counter=counter+1
                   sigma_mpi%cellCond%v(x,y,z) =model_para_vec(counter)
           end do
         end do
      end do
      end subroutine set_model_para_values 
#endif
