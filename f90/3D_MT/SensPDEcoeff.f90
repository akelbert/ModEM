! *****************************************************************************
!  Module that computes "forcings" for sensitivity calculation
!       (and adjoints).  This module is specific to the numerical
!       implementation of the solver, in this case for 3D 
!       MT finite difference modeling.   This module works with the
!       "natural" representations of conductivity: defined on edges
!        of the staggered grid.  Mappings from the potentially
!        more flexible earth conductivity parameter to these fixed, 
!        grid-specific representations are to be implemented in module 
!	 ModelParam.  This module has no dependence on the specific
!        conductivity parameterization
!
module senspdecoeff
   use math_constants
   use utilities
   use modelparameter
   use datafunc
   use sg_vector

   implicit none

   !  public routines
   public	::  Pmult, PmultT

   Contains

   !**********************************************************************
   subroutine Pmult(e0,dsigma,e)
   !   mapping from modelParam dsigma to source for forward problem
   !    (needed to calculate J*dsigma, where J is sensitivity)
   !   e0 is input background field solution;
   !    e is output ... used for forcing, created before calling 
   !    this routine

   type(EMsoln), intent(in)		:: e0
   type(modelParam_t), intent(in)		:: dsigma
   type(EMrhs), intent(inout)		:: e

   !  local variables
   complex(kind=selectedPrec)		:: minus_i_omega_mu
   type(rvector)			:: temp
   integer				:: k
   
   minus_i_omega_mu = cmplx(0.,-isign*mu*e0%omega,kind=selectedPrec)
   call create_rvector(e0%grid,temp,EDGE)

   ! map dsigma to edges, storing in array temp
   call ModelParamToEdge(dsigma,temp,e0%sigma)

   !  multiply temp by i_omeag_mu*e0, put result in e
   do k = 1,2
      call diagMult_crvector(e0%pol(k),temp,e%b(k)%s)
      call scMult_cvector(minus_i_omega_mu,e%b(k)%s,e%b(k)%s)
   enddo

   call deall_rvector(temp)

   end subroutine Pmult

   !**********************************************************************
   subroutine PmultT(e0,e,dsigmaReal,dsigmaImag)
   !   transpose of Pmult, mapping from adjoint soln e to dsigma
   !        e -> dsigma
   !   e0 is input background field solution
   !   NOTE: because the model parameter is real, while e is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...

   type(EMsoln), intent(in)			:: e0
   type(EMsoln), intent(in)			:: e
   type(modelParam_t), intent(inout)		:: dsigmaReal
   type(modelParam_t), intent(inout),optional	:: dsigmaImag

   !  local variables
   complex(kind=selectedPrec)			:: minus_i_omega_mu
   type(cvector)				:: Ctemp(2)
   type(rvector)				:: temp
   integer					:: k
   
   minus_i_omega_mu = cmplx(0.,-isign*mu*e0%omega,kind=selectedPrec)
   call create_rvector(e0%grid,temp,EDGE)
   call create_cvector(e0%grid,Ctemp(1),EDGE)
   call create_cvector(e0%grid,Ctemp(2),EDGE)

   ! multiply backward solutions (e) by minus_i_omega_mu * e0
   !   and sum over modes ... 
   do k = 1,2
      call diagMult_cvector(e0%pol(k),e%pol(k),Ctemp(k))
   enddo
   call add_cvector(Ctemp(1), Ctemp(2), Ctemp(1))
   call scMult_cvector(minus_i_omega_mu,Ctemp(1),Ctemp(1))

   ! map real/imag parts onto parameter space
   call getReal_cvector(Ctemp(1),temp)
   call EdgeToModelParam(temp,dsigmaReal,e0%sigma)

   if(present(dsigmaImag)) then
      ! also compute imaginary part
      call getImag_cvector(Ctemp(1),temp)
      call EdgeToModelParam(temp,dsigmaImag,e0%sigma)
   endif
   
   call deall_rvector(temp)
   call deall_cvector(Ctemp(1))
   call deall_cvector(Ctemp(2))

   end subroutine PmultT

   !**********************************************************************
end module senspdecoeff
