
#ifdef MPI
  
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
         if(associated(model_para_vec)) then
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
