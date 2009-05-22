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
module ModelSens
   use math_constants
   use utilities
   use datafunc
   use dataspace
   use transmitters
   use datatypes

   implicit none

   !  public routines
   public	::  Pmult, PmultT
   public   ::  Qmult, QmultT, QaddT

   Contains

   !**********************************************************************
   subroutine Pmult(e0,sigma0,dsigma,e)
   !   mapping from modelParam dsigma to source for forward problem
   !    (needed to calculate J*dsigma, where J is sensitivity)
   !   e0 is input background field solution;
   !    e is output ... used for forcing, created before calling
   !    this routine

   type(EMsoln_t), intent(in)		    :: e0
   type(modelParam_t), intent(in)	:: sigma0 ! used to compute e0
   type(modelParam_t), intent(in)	:: dsigma
   type(EMrhs_t), intent(inout)		:: e

   !  local variables
   complex(kind=prec)		:: minus_i_omega_mu
   type(rvector)			        :: temp
   integer				            :: k

   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
   call create_rvector(e0%grid,temp,EDGE)

   ! map dsigma to edges, storing in array temp
   call dModelParamToEdge(dsigma,temp,sigma0)

   !  multiply temp by i_omeag_mu*e0, put result in e
   do k = 1,2
      call diagMult_crvector(e0%pol(k),temp,e%b(k)%s)
      call scMult_cvector(minus_i_omega_mu,e%b(k)%s,e%b(k)%s)
   enddo

   call deall_rvector(temp)

   end subroutine Pmult

   !**********************************************************************
   subroutine PmultT(e0,sigma0,e,dsigmaReal,dsigmaImag)
   !   transpose of Pmult, mapping from adjoint soln e to dsigma
   !        e -> dsigma
   !   e0 is input background field solution
   !   NOTE: because the model parameter is real, while e is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...

   type(EMsoln_t), intent(in)			:: e0
   type(modelParam_t), intent(in)	:: sigma0 ! used to compute e0
   type(EMsoln_t), intent(in)			:: e
   type(modelParam_t), intent(inout)		:: dsigmaReal
   type(modelParam_t), intent(inout),optional	:: dsigmaImag

   !  local variables
   complex(kind=prec)			:: minus_i_omega_mu
   type(cvector)				:: Ctemp(2)
   type(rvector)				:: temp
   integer					:: k

   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
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
   temp = real(Ctemp(1))
   call dEdgeToModelParam(temp,dsigmaReal,sigma0)

   if(present(dsigmaImag)) then
      ! also compute imaginary part
      temp = imag(Ctemp(1))
      call dEdgeToModelParam(temp,dsigmaImag,sigma0)
   endif

   call deall_rvector(temp)
   call deall_cvector(Ctemp(1))
   call deall_cvector(Ctemp(2))

   end subroutine PmultT

   !**********************************************************************
   subroutine Qmult(e0,sigma0,dsigma,d)
   ! a dummy routine at present:
   ! outputs zero data vector when Q doesn't exist

   type(EMsoln_t), intent(in)             :: e0
   type(modelParam_t), intent(in)	    :: sigma0 ! used to compute e0
   type(modelParam_t), intent(in)		:: dsigma
   type(dataVec_t), intent(inout)          	:: d
   type(dataVec_t)							:: dTemp

   dTemp = zero(d)
   d = dTemp

   end subroutine Qmult

!**********************************************************************
   subroutine QmultT(e0,sigma0,d,dsigmaReal,dsigmaImag)
   !  a dummy routine at present:
   !  outputs zero vectors when Q doesn't exist

   type(EMsoln_t), intent(in)             :: e0
   type(modelParam_t), intent(in)	    :: sigma0 ! used to compute e0
   type(dataVec_t), intent(in)             :: d
   type(modelParam_t), intent(inout)               :: dsigmaReal
   type(modelParam_t), intent(inout),optional      :: dsigmaImag

   dsigmaReal = zero_modelParam(sigma0)

   if(present(dsigmaImag)) then
      dsigmaImag = zero_modelParam(sigma0)
   endif

  end subroutine QmultT

  !****************************************************************************
  subroutine QaddT(cs,dpsiT,sigma0,dsigmaReal,dsigmaImag)

   !   QaddT (formerly EMsparseQtoModelParam)
   !   Maps the input sparse vector to the model space by multiplying it with
   !   (d\pi/dm)^T. This is different from QmultT in that the latter acts on
   !   a data vector to obtain a model vector. Thus, QmultT is used to multiply
   !   with J^T, while this routine is needed to compute the full sensitivity
   !   matrix. In practice, used to make Q_J^T = (d\pi/dm)^T * (d\psi_j/d\pi)^T,
   !   then add cs*Q_j^T to (dsigmaReal,dSigmaImag)
   !   cs is a complex constant, Q is a sparse scalar field defined on
   !     grid cells (but represented with EMsparse data object ... the
   !     xyz component indicators are ignored).  dsigmaReal/dsigmaImag
   !     are used to collect sums of real and imaginary parts; dsigmaImag
   !     is optional.
   !   In 3D MT this has not been implemented yet.
   !   This is a dummy routine that makes and adds a zero model parameter.

   complex(kind=prec),intent(in)	:: cs
   type (EMsparse_t), intent(in)                  :: dpsiT
   type (modelParam_t), intent(in)                :: sigma0
   type (modelParam_t), intent(inout)             :: dsigmaReal
   type (modelParam_t), intent(inout),optional    :: dsigmaImag

   !  local variables
   type (modelParam_t)    :: csQReal, csQImag

   csQReal = zero_modelParam(sigma0)
   call linComb_modelParam(ONE,dsigmaReal,ONE,csQReal,dsigmaReal)

   if(present(dSigmaImag)) then
      csQImag = zero_modelParam(sigma0)
      call linComb_modelParam(ONE,dsigmaImag,ONE,csQImag,dsigmaImag)
   endif

   call deall_modelParam(csQReal)
   call deall_modelParam(csQImag)

  end subroutine QaddT

   !**********************************************************************
end module ModelSens
