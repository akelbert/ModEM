
module MPI_sub
#ifdef MPI

use math_constants
use utilities
use ForwardSolver
use MPI_declaration
use userctrl



      
      

Contains
!*************************************************************************
!Packing userdef_control in a Package
  !1- Allocate a place holder
 subroutine create_userdef_control_place_holder

     implicit none
     integer Nbytes1,Nbytes2,Nbytes3,Nbytes4
 

      
       CALL MPI_PACK_SIZE(80*19, MPI_CHARACTER,        MPI_COMM_WORLD, Nbytes1,  ierr)
       CALL MPI_PACK_SIZE(3,     MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes2,  ierr)
       CALL MPI_PACK_SIZE(1,     MPI_INTEGER,          MPI_COMM_WORLD, Nbytes3,  ierr)      
        Nbytes=(Nbytes1+Nbytes2+Nbytes3)+1  

         if(associated(userdef_control_package)) then
             deallocate(userdef_control_package)
         end if
             allocate(userdef_control_package(Nbytes)) 
      
 end subroutine create_userdef_control_place_holder 
 
    !*************************************************************************
    !2- Pack ctrl into userdef_control_package
 subroutine pack_userdef_control(ctrl)
    implicit none

     	type(userdef_control), intent(in)   :: ctrl
        integer index
 

      

       index=1 

        call MPI_Pack(ctrl%job,80*19, MPI_CHARACTER, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(ctrl%lambda,3, MPI_DOUBLE_PRECISION, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(ctrl%output_level,1, MPI_INTEGER, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
 
             

end subroutine pack_userdef_control 
    !*************************************************************************
    !3- Unpack userdef_control_package into ctrl
 subroutine unpack_userdef_control(ctrl)
    implicit none

     	type(userdef_control), intent(in)   :: ctrl
        integer index
 

      

       index=1 

   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%job,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_invCtrl,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_fwdCtrl,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Grid,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Model,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Data,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_dModel,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_EMsoln,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_EMrhs,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Grid,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Model,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Data,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_dModel,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_EMsoln,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_EMrhs,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Sens,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Cov,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%search,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%test,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%lambda,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)          
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%eps,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%delta,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
             
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%output_level,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

end subroutine unpack_userdef_control


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
   blockcounts(0) =  80*16

   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(1) = (80*16 * extent)+offsets(0)
   oldtypes(1) = MPI_REAL8
   blockcounts(1) = 3
   		
   call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
   offsets(2) = (3* extent)+offsets(1)
   oldtypes(2) = MPI_CHARACTER
   blockcounts(2) = 80*3
   
   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(3) = (80*3* extent)+offsets(2)
   oldtypes(3) = MPI_INTEGER
   blockcounts(3) = 1
   
     call MPI_TYPE_STRUCT(4, blockcounts, offsets, oldtypes,userdef_control_MPI, ierr)
     call MPI_TYPE_COMMIT(userdef_control_MPI, ierr)

end  subroutine create_userdef_control_MPI
!********************************************************************
subroutine check_userdef_control_MPI (which_proc,ctrl)

	type(userdef_control), intent(in)   :: ctrl
	character(20), intent(in)           :: which_proc
	
       write(6,*)trim(which_proc),' : ctrl%wFile_Sens ',trim(ctrl%wFile_Sens) 
       write(6,*)trim(which_proc),' : ctrl%lambda ',(ctrl%lambda)      
       write(6,*)trim(which_proc),' : ctrl%eps ',(ctrl%eps)    
       write(6,*)trim(which_proc),' : ctrl%rFile_Cov ',trim(ctrl%rFile_Cov)
       write(6,*)trim(which_proc),' : ctrl%search ',trim(ctrl%search)
       write(6,*)trim(which_proc),' : ctrl%output_level ',ctrl%output_level
       write(6,*)trim(which_proc),' : ctrl%rFile_fwdCtrl ',trim(ctrl%rFile_fwdCtrl)
       write(6,*)trim(which_proc),' : ctrl%rFile_invCtrl ',trim(ctrl%rFile_invCtrl)
       

end subroutine check_userdef_control_MPI
!********************************************************************
subroutine create_eAll_param_place_holder(eAll)

     implicit none
     type(solnVectorMTX_t), intent(in)	:: eAll
     integer Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3
 

      
       Ex_size=size(eAll%solns(which_per)%pol(1)%x)
       CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
       Ey_size=size(eAll%solns(which_per)%pol(1)%y)
       CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
       Ez_size=size(eAll%solns(which_per)%pol(1)%z)
       CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)
             
         Nbytes=(2*(Nbytes1+Nbytes2+Nbytes3))+1  ! Multiple by 2 for both polarizations

         if(associated(eAll_para_vec)) then
             deallocate(eAll_para_vec)
         end if
             allocate(eAll_para_vec(Nbytes)) 

             
             
 end subroutine create_eAll_param_place_holder 
 
 
 
!********************************************************************
subroutine get_eAll_para_vec(eAll)
    implicit none

     type(solnVectorMTX_t), intent(in)	:: eAll
     integer index,Ex_size,Ey_size,Ez_size
 

      
       Ex_size=size(eAll%solns(which_per)%pol(1)%x)
       Ey_size=size(eAll%solns(which_per)%pol(1)%y)
       Ez_size=size(eAll%solns(which_per)%pol(1)%z)
       

       

 
         
index=1 

        call MPI_Pack(eAll%solns(which_per)%pol(1)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(eAll%solns(which_per)%pol(1)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(eAll%solns(which_per)%pol(1)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)

        call MPI_Pack(eAll%solns(which_per)%pol(2)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(eAll%solns(which_per)%pol(2)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(eAll%solns(which_per)%pol(2)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)  
             

end subroutine get_eAll_para_vec 

!********************************************************************
subroutine set_eAll_para_vec(eAll)
    implicit none

     type(solnVectorMTX_t), intent(inout)	:: eAll

     integer index,Ex_size,Ey_size,Ez_size
 



       Ex_size=size(eAll%solns(which_per)%pol(1)%x)
       Ey_size=size(eAll%solns(which_per)%pol(1)%y)
       Ez_size=size(eAll%solns(which_per)%pol(1)%z)
       index=1

        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%pol(1)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%pol(1)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%pol(1)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)

        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%pol(2)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%pol(2)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%pol(2)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        	   
      
    
end subroutine set_eAll_para_vec







#endif

end module MPI_sub
