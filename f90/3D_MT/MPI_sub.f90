
module MPI_sub
#ifdef MPI

use math_constants
use utilities
use SolnSpace
use UserCtrl
use MPI_declaration



Contains



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
     integer index,Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3



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
