! *****************************************************************************
!  Module that computes "forcings" for sensitivity calculation
!       (and adjoints).  This module is specific to the numerical
!       implementation of the solver, in this case for 2D TE and TM
!       MT finite difference modeling.   This module works with the
!       "natural" representations of conductivity: defined on nodes
!        for TE, and on faces for TM; for TM data functionals conductvity
!        defined on cells is also required.  Mappings from the potentially
!        more flexible earth conductivity parameter to these fixed,
!        grid-specific representations is implemented in module CondMap,
!	 which this module uses.  This module should have no (or at least
!        minimal) dependence on the specific instance of the
!        conductivity parameterization

!
module ModelSens
   use math_constants
   use utilities
   use datafunc 	!!!! inherit: soln2d modelparameter
   use dataspace
   use transmitters
   use datatypes

   implicit none

   !  routines that are public/private
   public	::  Pmult, PmultT
   public   ::  Qmult, QmultT, QaddT

   private	::  curlB, curlE, curlE_T

   Contains
   !**********************************************************************
   subroutine curlB(b,Jy,Jz)
   ! computes curl B,  mapping from Node -> Face
   !   coded to map onto all faces, including boundaries
   !   Outputs Jy, Jz should be allocated as gridType EDGE_EARTH

   type(cvector), intent(in)	:: b
   type(cvector), intent(inout)	:: Jy,Jz

   !   local variables
   integer			:: iy, iz, Ny, Nzb,Nza
   real (kind=prec)	:: dz,dy

   if(Jy%gridType .ne. EDGE_EARTH .or.  &
		Jz%gridType .ne. EDGE_EARTH) then
      call errStop('wrong gridType for outputs Jy/Jz in curlB')
   endif

   Nza = b%grid%Nza
   Nzb = b%grid%Nz-Nza
   Ny  = b%grid%Ny

   do iy = 1,Ny+1
     do iz = 1,Nzb
       dz = b%grid%Dz(iz+Nza)
       Jy%v(iy,iz) =  (b%v(iy,iz+1) - b%v(iy,iz))/dz
     enddo
   enddo

   do iz = 1,Nzb+1
     do iy = 1,Ny
       dy = b%grid%Dy(iy)
       Jz%v(iy,iz) = -(b%v(iy+1,iz) - b%v(iy,iz))/dy
     enddo
   enddo

   end subroutine curlB

   !**********************************************************************
   subroutine curlE(Ey,Ez,b)
   ! computes curl E mapping from Edge -> Node
   !   maps only onto interior nodes
   !   Inputs Ey, Ez should be allocated as gridType EDGE_EARTH
   !   Output is of type NODE_EARTH

   type(cvector), intent(in)		:: Ey,Ez
   type(cvector), intent(inout)		:: b

   !  local variables
   integer 			:: iy, iz, Ny, Nzb, Nza
   real (kind=prec)	:: dz1,dz2,dzz, dy1, dy2, dyy

   if(Ey%gridType .ne. EDGE_EARTH .or.  &
		Ez%gridType .ne. EDGE_EARTH) then
      call errStop('wrong gridType for inputs Ey/Ez in curlE')
   endif

   Nza = Ey%grid%Nza
   Nzb = Ey%grid%Nz-Nza
   Ny  = Ey%grid%Ny

   b%v = R_ZERO
   do iy = 2,Ny
     do iz = 2,Nzb
       dz1 = Ey%grid%Dz(iz+Nza)
       dz2 = Ey%grid%Dz(iz-1+Nza)
       dzz = (dz1 + dz2)/TWO
       dy1 = Ey%grid%Dy(iy)
       dy2 = Ey%grid%Dy(iy-1)
       dyy = (dy1 + dy2)/TWO
       b%v(iy,iz) = (Ez%v(iy,iz) - Ez%v(iy-1,iz))/dyy &
                     + (Ey%v(iy,iz-1) - Ey%v(iy,iz))/dzz
     enddo
   enddo

   end subroutine curlE

   !**********************************************************************
   subroutine curlE_T(b,Ey,Ez)
   !  transpose of curlE, mapping from Node -> Face

   type(cvector), intent(in)	:: b
   type(cvector), intent(inout)	:: Ey, Ez

   ! local variables
   integer 			:: iy, iz, Ny, Nzb, Nza
   real(kind=prec)	:: dz1,dz2,dzz, dy1, dy2, dyy

   if(Ey%gridType .ne. EDGE_EARTH .or.  &
		Ez%gridType .ne. EDGE_EARTH) then
      call errStop('wrong gridType for outputs Ey/Ez in curlE_T')
   endif

   Nza = b%grid%Nza
   Nzb = b%grid%Nz-Nza
   Ny  = b%grid%Ny

   Ey%v = R_ZERO
   Ez%v = R_ZERO

   do iy = 2,Ny
     do iz = 2,Nzb
       dz1 = b%grid%Dz(iz+Nza)
       dz2 = b%grid%Dz(iz-1+Nza)
       dzz = (dz1 + dz2)/TWO
       dy1 = b%grid%Dy(iy)
       dy2 = b%grid%Dy(iy-1)
       dyy = (dy1 + dy2)/TWO
       Ez%v(iy,iz)   = Ez%v(iy,iz)   + (b%v(iy,iz))/dyy
       Ez%v(iy-1,iz) = Ez%v(iy-1,iz) - (b%v(iy,iz))/dyy
       Ey%v(iy,iz)   = Ey%v(iy,iz)   - (b%v(iy,iz))/dzz
       Ey%v(iy,iz-1) = Ey%v(iy,iz-1) + (b%v(iy,iz))/dzz
     enddo
   enddo

   end subroutine curlE_T

   !**********************************************************************
   subroutine Pmult(e0,sigma0,dsigma,e)
   !  generic Pmult (sorts out TE vs. TM from txDict):
   !   mapping from modelParam parameter dsigma to source;
   !        dsigma -> e
   !   e0 is input background field solution;
   !    (e is output used for forcing;
   !    The output is allocated before calling this routine)

   type(EMsoln_t), intent(in)             :: e0
   type(modelParam_t), intent(in)	    :: sigma0 ! used to compute e0
   type(modelParam_t), intent(in)		:: dsigma
   type(EMrhs_t), intent(inout)          	:: e
   !  local variables
   complex(kind=prec)		    :: i_omega_mu
   character*80                 	    :: gridType
   type(cvector)                 	    :: Jy,Jz,CJy,CJz

   if(e0%mode.eq.'TE') then

	   i_omega_mu = cmplx(0.,ISIGN*MU_0*e0%omega,kind=prec)
	   e%source%v = C_ZERO
	   call dModelParamToNode(dsigma,e%source,sigma0)

	   !  multiply by i * omega * mu
	   e%source%v = e%source%v*e0%vec%v*i_omega_mu

   else

	   !  allocate temporary data structures
	   gridType = EDGE_EARTH
	   call create_cvector(e0%grid,gridType,Jy)
	   call create_cvector(e0%grid,gridType,Jz)
	   call create_cvector(e0%grid,gridType,CJy)
	   call create_cvector(e0%grid,gridType,CJz)

	   call dModelParamToEdge(dsigma,CJy,CJz,sigma0)

	   call curlB(e0%vec,Jy,Jz)
	   CJy%v = CJy%v*Jy%v
	   CJz%v = CJz%v*Jz%v
	   call curlE(CJy,CJz,e%source)
	   !   If the differential operator were S = del x rho del x +- i wt
	   !    (as assumed in the derivation of P for TM mode in my notes)
	   !   then the sign of the output e%v should be reversed.
	   !   However, as implemented in WS TM forward solver the differential
	   !    operator is actually -S (this is the usual way that the
	   !       TM operator is written for 2D MT).  Hence the sign of e%v
	   !    should NOT be reversed here

	   call deall_cvector(Jy)
	   call deall_cvector(Jz)
	   call deall_cvector(CJy)
	   call deall_cvector(CJz)

   endif

   end subroutine Pmult

!**********************************************************************
   subroutine PmultT(e0,sigma0,e,dsigmaReal,dsigmaImag)
   !  generic PmultT (sorts out TE vs. TM from txDict):
   !   transpose of Pmult, mapping from adjoint soln e to sigma
   !   mapping from modelParam parameter dsigma to source;
   !        e -> dsigma
   !   e0 is input background field solution;
   !   NOTE: because the model parameter is real, while e is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...

   type(EMsoln_t), intent(in)             :: e0
   type(modelParam_t), intent(in)	    :: sigma0 ! used to compute e0
   type(EMsoln_t), intent(in)             :: e
   type(modelParam_t), intent(inout)               :: dsigmaReal
   type(modelParam_t), intent(inout),optional      :: dsigmaImag
   !  local variables
   complex(kind=prec)			:: i_omega_mu
   character*80					        :: gridType
   type(cvector)					    :: Jy,Jz,CJy,CJz
   type(cvector)					    :: temp

   if(e0%mode.eq.'TE') then

	   i_omega_mu = cmplx(0.,ISIGN*MU_0*e0%omega,kind=prec)

	   call create_cvector(e0%vec%grid,e0%vec%gridType,temp)
	   ! multiply backward solution by i_omega_mu * e0
	   ! map real/imag parts onto parameter space
	   temp%v = real(e%vec%v*e0%vec%v*i_omega_mu)

	   call dNodeToModelParam(temp,dsigmaReal,sigma0)

	   if(present(dsigmaImag)) then
	      ! also compute imaginary part
	      temp%v = imag(e%vec%v*e0%vec%v*i_omega_mu)
	      call dNodeToModelParam(temp,dsigmaImag,sigma0)
	   endif

	   call deall_cvector(temp)

   else

	   !  allocate temporary data structures
	   gridType = EDGE_EARTH
	   call create_cvector(e0%grid,gridType,Jy)
	   call create_cvector(e0%grid,gridType,Jz)
	   call create_cvector(e0%grid,gridType,CJy)
	   call create_cvector(e0%grid,gridType,CJz)

	   !  compute curls
	   call curlE_T(e%vec,CJy,CJz)
	   call curlB(e0%vec,Jy,Jz)
	   CJy%v = CJy%v*Jy%v
	   CJz%v = CJz%v*Jz%v

	   ! map from edge back to model parameter space
	   Jy%v = real(CJy%v)
	   Jz%v = real(CJz%v)
	   call dEdgeToModelParam(Jy,Jz,dsigmaReal,sigma0)

	   if(present(dsigmaImag)) then
	      Jy%v = imag(CJy%v)
	      Jz%v = imag(CJz%v)
	      call dEdgeToModelParam(Jy,Jz,dsigmaImag,sigma0)
	   endif

	   call deall_cvector(Jy)
	   call deall_cvector(Jz)
	   call deall_cvector(CJy)
	   call deall_cvector(CJz)

   endif

  end subroutine PmultT

   !**********************************************************************
   subroutine Qmult(e0,sigma0,dsigma,d)
   !  generic Qmult (sorts out TE vs. TM from txDict):
   !   mapping from modelParam parameter dsigma to the data space;
   !        dsigma -> d
   !   e0 is input background field solution

   type(EMsoln_t), intent(in)             :: e0
   type(modelParam_t), intent(in)	    :: sigma0 ! used to compute e0
   type(modelParam_t), intent(in)		:: dsigma
   type(dataVec_t), intent(inout)          	:: d
   !  local variables
   character*80                 	    :: gridType
   type(cvector)                 	    :: Jy,Jz,CJy,CJz

   if(e0%mode.eq.'TE') then

       call zero_dataVec(d)

   else

	   !  allocate temporary data structures
	   gridType = EDGE_EARTH
	   call create_cvector(e0%grid,gridType,Jy)
	   call create_cvector(e0%grid,gridType,Jz)
	   call create_cvector(e0%grid,gridType,CJy)
	   call create_cvector(e0%grid,gridType,CJz)

	   call dModelParamToEdge(dsigma,CJy,CJz,sigma0)

       ! instead, will implement
	   ! the mapping d\psi/d\pi, or multQtilde:
	   ! a linearized mapping from the edges to the data space
	   call zero_dataVec(d)

	   call deall_cvector(Jy)
	   call deall_cvector(Jz)
	   call deall_cvector(CJy)
	   call deall_cvector(CJz)

   endif

   end subroutine Qmult

!**********************************************************************
   subroutine QmultT(e0,sigma0,d,dsigmaReal,dsigmaImag)
   !  generic QmultT (sorts out TE vs. TM from txDict):
   !   transpose of Qmult, mapping from a data vector d to sigma
   !   mapping from modelParam parameter dsigma to source;
   !        d -> dsigma
   !   e0 is input background field solution;
   !   NOTE: because the model parameter is real, while e is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...
   !  outputs zero vectors when Q doesn't exist

   type(EMsoln_t), intent(in)             :: e0
   type(modelParam_t), intent(in)	    :: sigma0 ! used to compute e0
   type(dataVec_t), intent(in)             :: d
   type(modelParam_t), intent(inout)               :: dsigmaReal
   type(modelParam_t), intent(inout),optional      :: dsigmaImag
   !  local variables
   character*80					        :: gridType
   type(cvector)					    :: Jy,Jz,CJy,CJz

   if(e0%mode.eq.'TE') then

       dsigmaReal = sigma0
       call zero_modelParam(dsigmaReal)

	   if(present(dsigmaImag)) then
          dsigmaImag = sigma0
          call zero_modelParam(dsigmaImag)
	   endif

   else

	   ! after these initialization routines, will implement
	   ! the mapping (d\psi/d\pi)^T, or multQtildeT:
	   ! a linearized mapping from the data space to the edges
	   gridType = EDGE_EARTH
	   call create_cvector(e0%grid,gridType,Jy)
	   call create_cvector(e0%grid,gridType,Jz)
	   call create_cvector(e0%grid,gridType,CJy)
	   call create_cvector(e0%grid,gridType,CJz)

	   ! map from edge back to model parameter space
	   Jy%v = real(CJy%v)
	   Jz%v = real(CJz%v)
	   call dEdgeToModelParam(Jy,Jz,dsigmaReal,sigma0)

	   if(present(dsigmaImag)) then
	      Jy%v = imag(CJy%v)
	      Jz%v = imag(CJz%v)
	      call dEdgeToModelParam(Jy,Jz,dsigmaImag,sigma0)
	   endif

	   call deall_cvector(Jy)
	   call deall_cvector(Jz)
	   call deall_cvector(CJy)
	   call deall_cvector(CJz)

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

   complex(kind=prec),intent(in)	:: cs
   type (EMsparse_t), intent(in)                  :: dpsiT
   type (modelParam_t), intent(in)                :: sigma0
   type (modelParam_t), intent(inout)             :: dsigmaReal
   type (modelParam_t), intent(inout),optional    :: dsigmaImag

   !  local variables
   type (sparsevecc)	:: Ltemp
   type (modelParam_t)  :: csQReal, csQImag

   call scMult_sparsevecc(cs,dpsiT%L,Ltemp)
   call SparseCellToModelParam(Ltemp,sigma0,csQReal,csQImag)

   call linComb_modelParam(ONE,dsigmaReal,ONE,csQReal,dsigmaReal)

   if(present(dSigmaImag)) then
      call linComb_modelParam(ONE,dsigmaImag,ONE,csQImag,dsigmaImag)
   endif

   call deall_sparsevecc(Ltemp)
   call deall_modelParam(csQReal)
   call deall_modelParam(csQImag)

  end subroutine QaddT

   !**********************************************************************
end module ModelSens
