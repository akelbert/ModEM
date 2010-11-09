
module MPI_sub
#ifdef MPI

use math_constants
use utilities
use SolnSpace
use ioTypes
use MPI_declaration

implicit none




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
   blockcounts(0) =  80*23

   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(1) = (80*23 * extent)+offsets(0)
   oldtypes(1) = MPI_REAL8
   blockcounts(1) = 2

     call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes,userdef_control_MPI, ierr)
     call MPI_TYPE_COMMIT(userdef_control_MPI, ierr)

end  subroutine create_userdef_control_MPI

!********************************************************************
subroutine check_userdef_control_MPI (which_proc,ctrl)

    type(userdef_control), intent(in)   :: ctrl
    character(20), intent(in)           :: which_proc

       write(6,*) '[',trim(which_proc),'] Parametrization: ', trim(ctrl%paramname)
       write(6,*) '[',trim(which_proc),'] Output filename: ', trim(ctrl%modelname)
       write(6,*) '[',trim(which_proc),'] Level of output: ', trim(ctrl%verbose)
       write(6,*) '[',trim(which_proc),'] What to compute: ', trim(ctrl%calculate)
       write(6,*) '[',trim(which_proc),'] R_d + mu * R_m : ', ctrl%damping
       write(6,*) '[',trim(which_proc),'] Grid specified : ', trim(ctrl%fn_grid)
       write(6,*) '[',trim(which_proc),'] Grid base model: ', trim(ctrl%fn_rho)
       write(6,*) '[',trim(which_proc),'] Thinshell distr: ', trim(ctrl%fn_shell)
       write(6,*) '[',trim(which_proc),'] Radial EM field: ', trim(ctrl%fn_field)
       write(6,*) '[',trim(which_proc),'] Period filename: ', trim(ctrl%fn_period)
       write(6,*) '[',trim(which_proc),'] Obs coordinates: ', trim(ctrl%fn_coords)
       write(6,*) '[',trim(which_proc),'] Functional info: ', trim(ctrl%fn_func)
       write(6,*) '[',trim(which_proc),'] Forward control: ', trim(ctrl%fn_ctrl)
       write(6,*) '[',trim(which_proc),'] Inverse control: ', trim(ctrl%fn_invctrl)
       write(6,*) '[',trim(which_proc),'] Radii to output: ', trim(ctrl%fn_slices)
       write(6,*) '[',trim(which_proc),'] Base model info: ', trim(ctrl%fn_param0)
       write(6,*) '[',trim(which_proc),'] Parameters info: ', trim(ctrl%fn_param)
       write(6,*) '[',trim(which_proc),'] Interior source: ', trim(ctrl%fn_source)
       write(6,*) '[',trim(which_proc),'] Field data info: ', trim(ctrl%fn_cdata)
       write(6,*) '[',trim(which_proc),'] Field data info: ', trim(ctrl%fn_ddata)
       write(6,*) '[',trim(which_proc),'] Data misfit fn : ', trim(ctrl%fn_misfit)
       write(6,*) '[',trim(which_proc),'] Gradient fname : ', trim(ctrl%fn_gradient)
       write(6,*) '[',trim(which_proc),'] Point filename : ', trim(ctrl%fn_point)


end subroutine check_userdef_control_MPI


!********************************************************************
subroutine create_eAll_param_place_holder(eAll)

     implicit none
     type(solnVectorMTX_t), intent(in)  :: eAll
     integer index,Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3



       Ex_size=size(eAll%solns(which_per)%vec%x)
       CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
       Ey_size=size(eAll%solns(which_per)%vec%y)
       CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
       Ez_size=size(eAll%solns(which_per)%vec%z)
       CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)

         Nbytes=Nbytes1+Nbytes2+Nbytes3+1

         if(associated(eAll_para_vec)) then
             deallocate(eAll_para_vec)
         end if
         allocate(eAll_para_vec(Nbytes))



 end subroutine create_eAll_param_place_holder



!********************************************************************
subroutine get_eAll_para_vec(eAll)
    implicit none

     type(solnVectorMTX_t), intent(in)  :: eAll
     integer index,Ex_size,Ey_size,Ez_size



       Ex_size=size(eAll%solns(which_per)%vec%x)
       Ey_size=size(eAll%solns(which_per)%vec%y)
       Ez_size=size(eAll%solns(which_per)%vec%z)
       index=1

        call MPI_Pack(eAll%solns(which_per)%vec%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(eAll%solns(which_per)%vec%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(eAll%solns(which_per)%vec%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine get_eAll_para_vec

!********************************************************************
subroutine set_eAll_para_vec(eAll)
    implicit none

     type(solnVectorMTX_t), intent(inout)   :: eAll

     integer index,Ex_size,Ey_size,Ez_size




       Ex_size=size(eAll%solns(which_per)%vec%x)
       Ey_size=size(eAll%solns(which_per)%vec%y)
       Ez_size=size(eAll%solns(which_per)%vec%z)
       index=1

        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%vec%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%vec%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, eAll%solns(which_per)%vec%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)


end subroutine set_eAll_para_vec



#endif

end module MPI_sub
