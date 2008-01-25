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
module senspdecoeff
   use math_constants
   use utilities
   use datafunc 	!!!! inherit: soln2d modelparameter

   implicit none

   !  routines that are public/private
   public	::  Pmult, PmultT
   !   routines listed below could be merged into Pmult, PmultT
   !    don't need to be separate!
   private	::  P_TE, P_TE_T, P_TM, P_TM_T
   private	::  curlB, curlE, curlE_T

   Contains
   !**********************************************************************
   subroutine curlB(b,Jy,Jz)
   ! computes curl B,  mapping from Node -> Face
   !   coded to map onto all faces, including boundaries
   !   Outputs Jy, Jz should be allocated as gridType 'FACE EARTH'

   type(cvector), intent(in)	:: b
   type(cvector), intent(inout)	:: Jy,Jz

   !   local variables
   integer			:: iy, iz, Ny, Nzb,Nza
   real (kind=selectedPrec)	:: dz,dy
   character*80			:: msg

   if(Jy%gridType .ne. 'FACE EARTH' .or.  &
		Jy%gridType .ne. 'FACE EARTH') then
      msg = 'Error: wrong gridType for outputs Jy/Jz in curlB'
      call errStop(msg)
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
   ! computes curl E mapping from Face -> Node
   !   maps only onto interior nodes
   !   Inputs Ey, Ez should be allocated as gridType 'FACE EARTH'
   !   Output is of type 'NODE EARTH'

   type(cvector), intent(in)		:: Ey,Ez
   type(cvector), intent(inout)		:: b
 
   !  local variables
   integer 			:: iy, iz, Ny, Nzb, Nza
   real (kind=selectedPrec)	:: dz1,dz2,dzz, dy1, dy2, dyy
   character*80			:: msg

   if(Ey%gridType .ne. 'FACE EARTH' .or.  &
		Ey%gridType .ne. 'FACE EARTH') then
      msg = 'Error: wrong gridType for inputs Ey/Ez in curlE'
      call errStop(msg)
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
   real(kind=selectedPrec)	:: dz1,dz2,dzz, dy1, dy2, dyy
   character*80			:: msg

   if(Ey%gridType .ne. 'FACE EARTH' .or.  &
		Ey%gridType .ne. 'FACE EARTH') then
      msg = 'Error: wrong gridType for outputs Ey/Ez in curlE_T'
      call errStop(msg)
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
   subroutine P_TE(e0,dsigma,e)
   !   mapping from EarthCond parameter dsigma to source for TE prob;
   !        dsigma -> e
   !   e0 is input background field solution;
   !    (e is output ... used for forcing, but at present we store
   !     this as type cvector, and copy into forcing later;
   !    The output is allocated before calling this routine
   !   Calls NodeToCell  ... should this be renamed???

   type(EMsoln), intent(in)		:: e0
   type(modelParam_t), intent(in)		:: dsigma
   type(EMrhs), intent(inout)		:: e

   !  local variables
   complex(kind=selectedPrec)		:: i_omega_mu
   
   i_omega_mu = cmplx(0.,isign*mu*e0%omega,kind=selectedPrec)
   e%source%v = C_ZERO
   call CellToNode(dsigma,e%source,e0%sigma)

   !  multiply by i * omega * mu
   e%source%v = e%source%v*e0%vec%v*i_omega_mu

   end subroutine P_TE

   !**********************************************************************
   subroutine P_TE_T(e0,e,dsigmaReal,dsigmaImag)
   !   transpose of P_TE, mapping from adjoint soln e to sigma
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
   complex(kind=selectedPrec)			:: i_omega_mu
   type(cvector)					:: temp
   
   i_omega_mu = cmplx(0.,isign*mu*e0%omega,kind=selectedPrec)

   call create_cvector(e0%vec%grid,e0%vec%gridType,temp)
   ! multiply backward solution by i_omega_mu * e0
   ! map real/imag parts onto parameter space
   temp%v = real(e%vec%v*e0%vec%v*i_omega_mu)

   call NodeToCell(temp,dsigmaReal,e0%sigma)

   if(present(dsigmaImag)) then
      ! also compute imaginary part
      temp%v = imag(e%vec%v*e0%vec%v*i_omega_mu)
      call NodeToCell(temp,dsigmaImag,e0%sigma)
   endif
   
   call deall_cvector(temp)

   end subroutine P_TE_T

   !**********************************************************************
   subroutine P_TM(b0,dsigma,b)
   !   mapping from EarthCond parameter sigma to source for TM prob;
   !        dsigma -> b
   !   b0 is input background field solution;
   !   Allocate output b before calling
   !   Calls param_to_face_TE, etc

   type(EMsoln), intent(in)		:: b0
   type(modelParam_t), intent(in)		:: dsigma
   type(EMrhs), intent(inout)		:: b

   !  local variables
   character*80                 	:: gridType
   type(cvector)                 	:: Jy,Jz,CJy,CJz

   !  allocate temporary data structures
   gridType = 'FACE EARTH'
   call create_cvector(b0%grid,gridType,Jy)
   call create_cvector(b0%grid,gridType,Jz)
   call create_cvector(b0%grid,gridType,CJy)
   call create_cvector(b0%grid,gridType,CJz)

   call CellToFace(dsigma,b0%sigma, CJy,CJz)

   call curlB(b0%vec,Jy,Jz)
   CJy%v = CJy%v*Jy%v
   CJz%v = CJz%v*Jz%v
   call curlE(CJy,CJz,b%source)
   !   If the differential operator were S = del x rho del x +- i wt
   !    (as assumed in the derivation of P for TM mode in my notes)
   !   then the sign of the output b%v should be reversed.
   !   However, as implemented in WS TM forward solver the differential
   !    operator is actually -S (this is the usual way that the
   !       TM operator is written for 2D MT).  Hence the sign of b%v
   !    should NOT be reversed here

   call deall_cvector(Jy)
   call deall_cvector(Jz)
   call deall_cvector(CJy)
   call deall_cvector(CJz)

   end subroutine P_TM

   !**********************************************************************
   subroutine P_TM_T(b0,b,dsigmaReal,dsigmaImag)
   !   Transpose of P_TM, mapping adjoint soln b to sigma
   !         b -> dsigma
   !   b0 is input background field solution;
   !   Calls  ...

   type(EMsoln), intent(in)			:: b0
   type(EMsoln), intent(in)			:: b
   type(modelParam_t), intent(inout)		:: dsigmaReal
   type(modelParam_t), intent(inout), optional	:: dsigmaImag

   !  local variables
   character*80					:: gridType
   type(cvector)					:: Jy,Jz,CJy,CJz

   !  allocate temporary data structures
   gridType = 'FACE EARTH'
   call create_cvector(b0%grid,gridType,Jy)
   call create_cvector(b0%grid,gridType,Jz)
   call create_cvector(b0%grid,gridType,CJy)
   call create_cvector(b0%grid,gridType,CJz)

   !  compute curls
   call curlE_T(b%vec,CJy,CJz)
   call curlB(b0%vec,Jy,Jz)
   CJy%v = CJy%v*Jy%v
   CJz%v = CJz%v*Jz%v

   ! map from face back to model parameter space
   Jy%v = real(CJy%v)
   Jz%v = real(CJz%v)
   call FaceToCell(Jy,Jz,b0%sigma,dsigmaReal)

   if(present(dsigmaImag)) then
      Jy%v = imag(CJy%v)
      Jz%v = imag(CJz%v)
      call FaceToCell(Jy,Jz,b0%sigma,dsigmaImag)
   endif

   call deall_cvector(Jy)
   call deall_cvector(Jz)
   call deall_cvector(CJy)
   call deall_cvector(CJz)

   end subroutine P_TM_T

   !**********************************************************************
   subroutine Pmult(e0,dsigma,e)
   !  generic Pmult: sorts out TE vs. TM from itx,
   !   and calls appropriate P_TE or P_TM routine
   !  (probably should be recoded to get rid of separate
   !   TE/TM routines)

   type(EMsoln), intent(in)             :: e0
   type(modelParam_t), intent(in)		:: dsigma
   type(EMrhs), intent(inout)          	:: e

   if(e0%mode.eq.'TE') then
      call P_TE(e0,dsigma,e)
   else
      call P_TM(e0,dsigma,e)
   endif
   
   end subroutine Pmult

!**********************************************************************
   subroutine PmultT(e0,e,dsigmaReal,dsigmaImag)
   !  generic PmultT: sorts our TE vs. TM from itx,
   !   and calls appropriate P_TE_T or P_TM_T routine
   !  (probably should be recoded to get rid of separate
   !   TE/TM routines)

   type(EMsoln), intent(in)             :: e0
   type(EMsoln), intent(in)             :: e
   type(modelParam_t), intent(inout)               :: dsigmaReal
   type(modelParam_t), intent(inout),optional      :: dsigmaImag


   if(present(dsigmaImag)) then
      if(e0%mode.eq.'TE') then
         call P_TE_T(e0,e,dsigmaReal,dsigmaImag)
      else
         call P_TM_T(e0,e,dsigmaReal,dsigmaImag)
      endif
   else
      if(e0%mode.eq.'TE') then
         call P_TE_T(e0,e,dsigmaReal)
      else
         call P_TM_T(e0,e,dsigmaReal)
      endif
   endif

   end subroutine PmultT

   !**********************************************************************
end module senspdecoeff
