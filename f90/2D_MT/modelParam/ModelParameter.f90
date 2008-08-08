!include "modelCov_WS.inc"

module modelparameter
!   module to define modelparam abstract data type,
!    and all methods required for creation, algebra, and inner products,
!    including covariances.  Attributes of the data type
!    are private to this module to force coding of other modules
!    to be completely independent of the specific parameterization
!    instance/implementation 
!
! Defines: modelParam
! Uses: EMfield, Grid2D, UTILS

use math_constants
use utilities
use emfield     !!!!  inherits :  grid2d
!use modelcov_ws

implicit none

  ! supported model parameter types (conductivity only)
   character(len=80), parameter		:: LOGE = 'LOGE'
   character(len=80), parameter		:: LINEAR = 'LINEAR'

type ::  modelParam_t
   !  modelParam is derived data type used to store parameters that
   !   determine conductivity/resistivity model
   !   This specific implementation is based on blocks
   !    allows for linear/log conductivity/resistivity
   private
   integer   				:: Ny = 0
   integer   				:: NzEarth = 0
   type(grid2d_t), pointer  		:: grid
   real (kind=selectedPrec), pointer, dimension(:,:) :: v
   real (kind=selectedPrec)		:: AirCond = SIGMA_AIR
   logical           			:: allocated = .false.
   !   parameter types supported:
   !   LINEAR = linear conductivity of each grid Earth cell is specified
   !   LOGE = natural log of conductivity in Earth cells specified
   character*80      			:: paramType = ''
end type modelParam_t

interface assignment (=)
   MODULE PROCEDURE copy_modelParam
end interface

interface operator (.dot.)
   MODULE PROCEDURE dotProd_modelParam
end interface

interface operator (*)
   MODULE PROCEDURE scMult_modelParam
end interface

interface modelCov
   MODULE PROCEDURE WScov
end interface

public	::	create_modelParam,deall_modelParam,dotProd_modelParam
public	::	linComb_modelParam,zero_modelParam,copy_modelParam
public  ::  scMult_modelParam, getValue_modelParam, getSize_modelParam
public  ::  maxNorm_modelParam

!   include interface for conductivity mappings
!    include CondMap2D.hd
! *****************************************************************************
!  conductivity mappings and adjoints for 2D MT modeling and inversion code 
!  routines that are public
   public	:: CellToNode, NodeToCell, CondParamToArray
   public	:: rhoC, EdgeToCell, CellToEdge, QtoModelParam


!  public          ::  WScov
!   include interface for model covariance
!   include wscovar.hd
  type, private :: mCov
    ! model covariance 
    real (kind=selectedPrec) :: ylen = R_ZERO
    real (kind=selectedPrec) :: zlen = R_ZERO
    real (kind=selectedPrec), pointer, dimension(:,:,:)   ::  YDIF
    real (kind=selectedPrec), pointer, dimension(:,:,:)   ::  ZDIF
  end type

  type (mCov), private,  save       ::  Cm
  
!   include interface for model parameter IO
!   include ModelParamIO.hd
   public writeAll_Cond2D, readAll_Cond2D

contains

!  The included file must contain subroutines create_CmSqrt, deall_CmSqrt, multBy...
   include "modelCov/WS.inc"
   
!************************************************************************
   !  allocateEarthCond allocates and initializes arrays for
   !   Earth-cell conductivity structure;
   !   Pass grid of type grid2d_t to set array sizes
   subroutine create_modelParam(grid,paramtype,cond,value,airCond)

     implicit none
     type (grid2d_t), intent(in), target   	:: grid
     character(*), intent(in)   		    :: paramtype
     type (modelParam_t), intent(inout)   	:: cond
     real (kind=selectedPrec), intent(in), optional :: value(:,:)
     real (kind=selectedPrec), intent(in), optional :: airCond
     !  local variables
     integer ::       Nz,Ny,Nza,NzEarth

	 if (cond%allocated) then
        call warning('Model parameter already allocated in create_modelParam')
        return
	 end if

     Nz = grid%Nz
     Nza = grid%Nza
     NzEarth = Nz-Nza
     Ny = grid%Ny
     cond%Ny = Ny
     cond%NzEarth = NzEarth
     cond%grid => grid
     cond%paramType = paramtype

     allocate(cond%v(Ny,NzEarth))
     if (present(value)) then
        if((size(value,1)==Ny).and.(size(value,2)==NzEarth)) then
           cond%v = value
        else
           call errStop('Wrong conductivity array size in create_modelParam')
        end if
     else
        cond%v = R_ZERO
     end if
     cond%allocated = .true.

     if (present(airCond)) then
        cond%AirCond = airCond
     end if

   end subroutine create_modelParam

  !************************************************************
   subroutine deall_modelParam(cond)
     implicit none
     type (modelParam_t), intent(inout)   :: cond

     if(cond%allocated) then
        deallocate(cond%v)
        nullify(cond%v)
        nullify(cond%grid)
        cond%allocated = .false.
        cond%paramType = ''
     endif

   end subroutine deall_modelParam

!**********************************************************************
   function dotProd_modelParam(m1,m2) result(r)

   !   dot product of two model space parameter objects
   !    here implemented for 2D blocks
     real(kind=selectedPrec)		:: r
     type(modelParam_t), intent(in)	:: m1,m2

     ! local variables
     integer	:: j,k

     if((m1%Ny .ne. m2%Ny).or. (m1%NzEarth .ne. m2%NzEarth)) then
        write(6,*) 'Ny ',m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop('size of m1, m2 incompatable in dotProd_modelParam') 
     endif

     r = R_ZERO
     do k = 1,m1%NzEarth
        do j = 1,m1%Ny
           r = r + m1%v(j,k)*m2%v(j,k) 
        enddo
     enddo
  
   end function dotProd_modelParam

!**********************************************************************
   function maxNorm_modelParam(m) result(r)

   !   Max norm = max(|m_ij|)
   !    here implemented for 2D blocks
     real(kind=selectedPrec)		:: r
     type(modelParam_t), intent(in)	:: m

     ! local variables
     integer	:: j,k

     if(.not.(m%allocated)) then
        call errStop('m not allocated in maxNorm_modelParam') 
     endif

     r = max(abs(minval(m%v)),abs(maxval(m%v)))
  
   end function maxNorm_modelParam
   
!**********************************************************************
   subroutine zero_modelParam(m)

     !  zeros a model space object
   
     type(modelParam_t), intent(inout)	:: m

     m%v = R_ZERO

   end subroutine zero_modelParam

!**********************************************************************
   subroutine linComb_modelParam(a1,m1,a2,m2,m)

     !  forms the linear combination of model parameters
     !    m = a1*m1 + a2*m2
     !  where a1 and a2 are real constants and m1 and m2
     !   are model parameters
     !   output m may overwrite m1 or m2

     real(kind=selectedPrec),intent(in)		:: a1,a2
     type(modelParam_t), intent(in)	:: m1,m2
     type(modelParam_t), intent(inout)	:: m
      
     if((m1%Ny .ne. m2%Ny).or. (m1%NzEarth .ne. m2%NzEarth)) then
        call errStop('size of m1, m2 incompatable in linComb_modelParam') 
     endif
     if(m1%paramType .ne. m2%paramType) then
        call errStop('paramType incompatable in linComb_modelParam') 
     endif

     ! make sure m is allocated, and of same size; if not same
     !   size, deallocate and allocate as correct size; otherwise
     !   do nothing (use m as input ... allowing m to overwrite m1)
     if(m%allocated) then
        if((m%Ny .ne. m1%Ny).or. (m%NzEarth .ne. m1%NzEarth)) then
           call deall_modelParam(m)
           call create_modelParam(m1%grid,m1%paramtype,m)
        endif
     else
        call create_modelParam(m1%grid,m1%paramtype,m)
     endif
     m%v = a1*m1%v + a2*m2%v
     !  we are implicitly assuming that if m1 and m2 are the same
     !   size other parameters are of the same type
     m%AirCond = m1%AirCond 

   end subroutine linComb_modelParam

  ! **********************************************************************
    function scMult_modelParam(a,mIn) result (mOut)
  !  computes mOut = a * mIn for modelParam object m and real scalar a

    real (kind=selectedPrec), intent(in)		:: a
    type(modelParam_t), intent(in)	            :: mIn
    type(modelParam_t)                            :: mOut

    ! check to see that input m is allocated
    if(.not.mIn%allocated) then
       call errStop('input not allocated on call to scDivide_DvecMTX')
    endif

	call linComb_modelParam(R_ZERO,mIn,a,mIn,mOut)
	  
  end function scMult_modelParam
  

   !**********************************************************************
   subroutine copy_modelParam(mOut,mIn)

     type(modelParam_t), intent(in)	:: mIn
     type(modelParam_t), intent(inout)	:: mOut

     ! if mOut is allocated, check to see if if is of same size as
     !   mIn; if not, deallocate and reallocate as correct size; otherwise
     !   use as input
     if(mOut%allocated) then
        if((mOut%Ny .ne. mIn%Ny).or. (mOut%NzEarth .ne. mIn%NzEarth)) then
           call deall_modelParam(mOut)
           call create_modelParam(mIn%grid,mIn%paramtype,mOut)
        endif
     else
        call create_modelParam(mIn%grid,mIn%paramtype,mOut)
     endif
     mOut%v = mIn%v
     mOut%AirCond = mIn%AirCond

   end subroutine copy_modelParam

   !************************************************************************
   !  getSize_modelParam extracts model size from a modelParam_t variable
   subroutine getSize_modelParam(cond,Ny,NzEarth)

     implicit none
     type (modelParam_t), intent(in)   	  :: cond
     integer, intent(out)                 :: Ny,NzEarth

     if (.not.cond%allocated) then
        call errStop('Model parameter not allocated in getValue_modelParam')
     end if

     Ny = cond%Ny
     NzEarth = cond%NzEarth
 
   end subroutine getSize_modelParam


   !************************************************************************
   !  getValue_modelParam extracts information for a modelParam_t variable;
   !  note that the output array must already be of correct size on input.
   !  Therefore, grid is assumed to be known before calling this subroutine.
   subroutine getValue_modelParam(cond,value,paramtype,airCond)

     implicit none
     type (modelParam_t), intent(in)   	        :: cond
     real (kind=selectedPrec), intent(inout)    :: value(:,:)
     character*80, intent(out), optional   		:: paramtype
     real (kind=selectedPrec), intent(out), optional :: airCond
     !  local variables
     integer ::       Nz,Ny,Nza,NzEarth

     if (.not.cond%allocated) then
        call errStop('Model parameter not allocated in getValue_modelParam')
     end if

     Ny = cond%Ny
     NzEarth = cond%NzEarth
 
     if((size(value,1)==Ny).and.(size(value,2)==NzEarth)) then
        value = cond%v
     else
        call errStop('Wrong conductivity array size in getValue_modelParam')
     end if
         
     if (present(paramtype)) then
        paramtype = cond%paramType
     end if
     
     if (present(airCond)) then
        airCond = cond%AirCond
     end if

   end subroutine getValue_modelParam

!  include source code for conductivity mappings
!  include CondMap2D.src

!!!!!!  >>>>>>>>>>>>>>>>>>>>   TE mode routines <<<<<<<<<<<<<<<<<<
! *****************************************************************************
      subroutine CellToNode(CellCond,NodeCond,Sigma0)
      !  LINEARIZED mapping from conductivity or log conductivity,
      !      defined on Earth cells to interior Earth node conductivities, 
      !      averaging over adjacent cells.  Overwrites inputs on 
      !      these cells, leaving  other nodes (in air, 
      !      on grid boundaries) unaltered
      !   Sigma0 is background conductivity for linearization
      !        -- this is an optional argument, 
      !           required when CellCond%paramType == 'LOGE'
      !    
      !   NOTE: CellCond array is real, while NodeCond array is
      !       complex  !!!
      !
      type (modelParam_t), intent(in)   :: CellCond
      type (modelParam_t), intent(in), optional   :: Sigma0
  
      !  OUTPUT is of type vec2D, of gridType NODE
      !   i.e., output should be defined on all grid nodes, though
      !   this routine only modifies interior nodes within the
      !   Earth.
      type (cvector), intent(inout)   :: NodeCond
      
      ! local variables
      integer iy,iz,Ny,Nz,Nza,izc, NzE
      real(kind=selectedPrec) w00,w01,w10,w11,wsum
      real(kind=selectedPrec), allocatable, dimension(:,:)	:: temp

      Ny = cellCond%grid%Ny
      Nz = cellCond%grid%Nz
      Nza = cellCond%grid%Nza
      NzE = Nz-Nza
      allocate(temp(Ny,NzE))
      if(NodeCond%gridType .NE. NODE) then 
         call errStop('NodeCond must be of gridType NODE in CellToNode')
      endif
      if(CellCond%paramType .EQ. LOGE) then 
         if(present(Sigma0)) then
            temp = CellCond%v*exp(Sigma0%v)
         else
            call errStop('Sigma0 input to CellToNode is required for case LOGE')
         endif
      else
         temp = CellCond%v 
      endif

      !  loop over interior Earth nodes 
      !        (including air/earth interface)
      do iz = Nza+1,Nz
         izc = iz-Nza
         do iy = 2,Ny
            w00 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz-1)
            w10 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1)
            w01 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz)
            w11 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz)
            wsum = w00+w10+w01+w11
            w00 = w00/wsum
            w10 = w10/wsum
            w01 = w01/wsum
            w11 = w11/wsum
            if(izc .EQ. 1) then
               NodeCond%v(iy,iz) = w01*temp(iy-1,izc) + &
			      w11*temp(iy,izc)
            else
               NodeCond%v(iy,iz) = w00*temp(iy-1,izc-1)+ &
			      w10*temp(iy,izc-1) + &
			      w01*temp(iy-1,izc) + &
			      w11*temp(iy,izc)
            endif
         enddo            
      enddo 
      deallocate(temp)           
      end subroutine CellToNode

   ! *****************************************************************************
      subroutine NodeToCell(NodeCond,CellCond,Sigma0)
      !  Maps from interior Earth nodes to Earth cells;
      !   transpose (adjoint) of EarthCellToNode
      !  INPUT is of type vec2D, of gridType NODE
      !   NodeCond is assumed defined on all nodes, but only
      !   interior Earth nodes are used 
      type (cvector), intent(in)   :: NodeCond
      

      !  OUTPUT is linear cell conductivity parameter structure
      !    (paramType should be LINEAR)
      !  CellCond is overwritten by call to this routine
      type (modelParam_t), intent(inout)   :: CellCond
      !   Optional input Sigma0 is background conductivity, required
      !        for log conductivity parameterization          
      type (modelParam_t), intent(in), optional   :: Sigma0
      
      !   NOTE: CellCond array is real, while NodeCond array is
      !       complex  !!!

      ! local variables
      integer iy,iz,izc,Ny,Nz,Nza
      real(kind=selectedPrec) w00,w01,w10,w11,wsum

      if(NodeCond%gridType .NE. NODE) then 
         call errStop('NodeCond must be of gridType NODE in NodeToCell')
      endif
      if((CellCond%paramType .eq. LOGE) .and. (.not.present(Sigma0))) then
         call errStop('Background conductivity required for paramType LOGE in NodeToCell')
      endif

      Ny = cellCond%grid%Ny
      Nz = cellCond%grid%Nz
      Nza = cellCond%grid%Nza

      !   zero cell conductivity values
      CellCond%v = R_ZERO

      !  loop over interior Earth nodes 
      !        (including air/earth interface)
      do iz = Nza+1,Nz
         izc = iz-Nza
         do iy = 2,Ny
            w00 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz-1) 
            w10 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1)
            w01 = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz)
            w11 = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz)
            wsum = w00+w10+w01+w11
            w00 = w00/wsum
            w10 = w10/wsum
            w01 = w01/wsum
            w11 = w11/wsum
            if(izc .GT. 1) then
               CellCond%v(iy-1,izc-1) = &
		 CellCond%v(iy-1,izc-1)+w00*real(NodeCond%v(iy,iz))
               CellCond%v(iy,izc-1) = &
		 CellCond%v(iy,izc-1)+w10*real(NodeCond%v(iy,iz))
            endif
            CellCond%v(iy-1,izc) = &
		 CellCond%v(iy-1,izc)+w01*real(NodeCond%v(iy,iz))
            CellCond%v(iy,izc) = &
		 CellCond%v(iy,izc)+w11*real(NodeCond%v(iy,iz))
         enddo           
      enddo
      if(CellCond%paramType .eq. LOGE) then
         CellCond%v = CellCond%v*exp(Sigma0%v)
      endif

      end subroutine NodeToCell

      !**********************************************************************
      subroutine CondParamToArray(CellCond,Ny,Nz,Cond)
      !  Copies from formal Earth cell conductivity structure
      !    to a standard real array defined for the whole grid.
      !     Sets air cells to AirCond
      !
      !  INPUT is cell conductivity parameter structure
      !    (i.e., paramType should be LINEAR)
      type (modelParam_t), intent(in)   :: CellCond
      integer, intent(in)	:: Ny,Nz
  
      real(kind=selectedPrec), intent(inout)   :: Cond(Ny,Nz)
      
      ! local variables
      integer iy,iz,izc,NyC,NzC,Nza

      NyC = cellCond%grid%Ny
      NzC = cellCond%grid%Nz
      Nza = cellCond%grid%Nza

      if((NyC .NE. Ny).or.NzC .NE. Nz) then
        call errStop('dimensions of CellCond, Cond disagree in CondParamToArray')
      endif
     
      Cond = CellCond%AirCond
      do iz = Nza+1,Nz
         izc = iz-Nza
         do iy = 1,Ny
            Cond(iy,iz) = CellCond%v(iy,izc)
         enddo            
      enddo  
      if(CellCond%paramType .eq. LOGE) then
         Cond = exp(Cond)
      endif
                
      end subroutine CondParamToArray

!!!!!!  >>>>>>>>>>>>>>>>>>>>   TM mode routines <<<<<<<<<<<<<<<<<<
      !**********************************************************************
      !   computes resistivity for cell j,k using input parameter
      !     structure sigma.  This function defines how the abstract
      !     conductivity parameter is mapped to cell resistivites needed
      !     for TM electric field interpolation functions.  The derivative
      !     of this function is required for evaluation of linearized
      !     data functionals, and for construction of the direct 
      !     parameter space part of the comb
      function rhoC(sigma,j,k) result(r)
      type (modelParam_t), intent(in)	:: sigma
      integer, intent(in)		:: j,k
      real(kind=selectedPrec)	                :: r

      ! local variables
      if(sigma%paramType .eq. LOGE) then
         r = exp(-sigma%v(j,k))
      else
         r = 1./sigma%v(j,k)
      endif

      end function rhoC

   !**********************************************************************
   subroutine CellToEdge(CellCond, sigma0, EdgeCond_y,EdgeCond_z)
   ! map from Cell to Edge on Y and Z axis

   type (modelParam_t), intent(in)    :: CellCond
   type (modelParam_t), intent(in)    :: sigma0
   type (cvector), intent(inout)   :: EdgeCond_y, EdgeCond_z

   real (kind=selectedPrec), allocatable, dimension(:,:) :: temp


   integer :: iz, iy, Nz, Nza, Ny, Nzb
   real (kind=selectedPrec)    ::  wll, wrr, wdd, wuu, wsum

   Ny  = CellCond%grid%Ny
   Nz  = CellCond%grid%Nz
   Nza = CellCond%grid%Nza
   Nzb = Nz - Nza

   allocate(temp(Ny,Nzb))
   if(CellCond%paramType .EQ. LOGE) then
      temp = - CellCond%v/exp(sigma0%v)
   else
      temp = - CellCond%v/(sigma0%v*sigma0%v)
   endif
   do iy = 2, Ny
     do iz = 1, Nzb
       wll = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz+Nza)
       wrr = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wll + wrr
       wll = wll/wsum
       wrr = wrr/wsum
       EdgeCond_y%v(iy,iz) = wll*(temp(iy-1,iz)) + &
                             wrr*(temp(iy,iz))
     enddo 
   enddo 

   do iz = 2, Nzb
     do iy = 1, Ny
       wuu = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1+Nza)
       wdd = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wuu + wdd
       wuu = wuu/wsum
       wdd = wdd/wsum
       EdgeCond_z%v(iy,iz) = wuu*(temp(iy,iz-1)) + &
                             wdd*(temp(iy,iz))
     enddo
   enddo

   deallocate(temp)

   end subroutine CellToEdge

   !**********************************************************************

   subroutine EdgeToCell(EdgeCond_y,EdgeCond_z,sigma0,CellCond)
   ! map from Edge to Cell on Y and Z axis

   type (cvector), intent(in)   :: EdgeCond_y, EdgeCond_z
   type (modelParam_t), intent(in) :: sigma0
   type (modelParam_t), intent(inout)    :: CellCond


   integer :: iz, iy, Nz, Nza, Ny, Nzb
   real (kind=selectedPrec)    ::  wll, wrr, wdd, wuu, wsum

   Ny  = CellCond%grid%Ny
   Nz  = CellCond%grid%Nz
   Nza = CellCond%grid%Nza
   Nzb = Nz - Nza

   CellCond%v = 0.
   do iy = 2, Ny
     do iz = 1, Nzb
       wll = CellCond%grid%Dy(iy-1)*CellCond%grid%Dz(iz+Nza)
       wrr = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wll + wrr
       wll = wll/wsum
       wrr = wrr/wsum
       CellCond%v(iy-1,iz) = CellCond%v(iy-1,iz) + wll*real(EdgeCond_y%v(iy,iz))
       CellCond%v(iy,iz)   = CellCond%v(iy,iz)   + wrr*real(EdgeCond_y%v(iy,iz))
     enddo 
   enddo 

   do iz = 2, Nzb
     do iy = 1, Ny
       wuu = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz-1+Nza)
       wdd = CellCond%grid%Dy(iy)*CellCond%grid%Dz(iz+Nza)
       wsum = wuu + wdd
       wuu = wuu/wsum
       wdd = wdd/wsum
       CellCond%v(iy,iz-1) = CellCond%v(iy,iz-1) + wuu*real(EdgeCond_z%v(iy,iz))
       CellCond%v(iy,iz)   = CellCond%v(iy,iz)   + wdd*real(EdgeCond_z%v(iy,iz))
     enddo
   enddo

   if (CellCond%paramType .eq. LOGE) then
      CellCond%v = - CellCond%v/exp(sigma0%v)
   else
      CellCond%v = - CellCond%v/(sigma0%v*sigma0%v)
   endif

   end subroutine EdgeToCell

   !**********************************************************************
   subroutine QtoModelParam(Qj,sigma0,dsigmaReal,dSigmaImag)

   !  given input sparse vector defined on cells,
   !   compute vector defined on parameter space
   !   add result to inputs dsigmaReal, dsigmaImag
   !    Q_j^T \partial\rho_cell/\partial\theta

   type (sparsevecc), intent(in)   :: Qj
   type (modelParam_t), intent(in) :: sigma0
   type (modelParam_t), intent(inout)    :: dsigmaReal
   type (modelParam_t), intent(inout),optional    :: dsigmaImag

   ! local variables
   integer		:: jk,j,k

   do jk = 1,Qj%nCoeff
      j = Qj%J(jk)
      k = Qj%K(jk)
      if(sigma0%paramType .eq. LOGE) then
         dsigmaReal%v(j,k) = dSigmaReal%v(j,k) &
			-real(Qj%C(jk))/exp(sigma0%v(j,k))
         if(present(dSigmaImag)) then
            dsigmaImag%v(j,k) = dSigmaImag%v(j,k) &
			-imag(Qj%C(jk))/exp(sigma0%v(j,k))
         endif
      else
         dsigmaReal%v(j,k) = dSigmaReal%v(j,k) &
			-real(Qj%C(jk))/(sigma0%v(j,k)*sigma0%v(j,k))
         if(present(dSigmaImag)) then
            dsigmaImag%v(j,k) = dSigmaImag%v(j,k) &
		-imag(Qj%C(jk))/(sigma0%v(j,k)*sigma0%v(j,k))
         endif
      endif
   enddo

   end subroutine QtoModelParam
!  include source code for covariance routines
!  include wscovar.src


!  include source code for model parameter IO
!  include ModelParamIO.src
     !******************************************************************
      subroutine writeAll_Cond2D(fid,cfile,nSigma,sigma,header)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file

      integer, intent(in)		:: fid,nSigma
      character(*), intent(in)		:: cfile, header
      type(modelParam_t), intent(in)	:: sigma(nSigma) 

      integer i
      character(80) temp

      temp = header

      open(unit=fid, file=cfile, form='unformatted')
      write(fid) temp
      write(fid) nSigma
      do i = 1,nSigma
         write(fid) sigma(i)%paramType
         write(fid) sigma(i)%Ny,sigma(i)%NzEarth
         write(fid) sigma(i)%v
         write(fid) sigma(i)%AirCond
      enddo
      close(fid)
      end subroutine writeAll_Cond2D

     !******************************************************************
      subroutine readAll_Cond2D(fid,cfile,nSigma,sigma,header)

      !  open cfile on unit fid, read nSigma objects of
      !   type modelParam , closes file
      !  sigma(nsigma) must be allocated before calling

      integer, intent(in)		:: fid
      character(*), intent(in)		:: cfile
      integer, intent(in)		:: nSigma
      character(80), intent(out)		:: header
      type(modelParam_t), intent(inout) 	:: sigma(nsigma)

      ! local variables
      integer i, nS, Ny,NzEarth

      open(unit=fid, file=cfile, form='unformatted')
      read(fid) header
      read(fid) nS
      if(nS .NE. nSigma) then
          call errStop('size of sigma does not agree with contents of file in readAll_Cond2D')
      endif

      do i = 1,nSigma
         read(fid) sigma(i)%paramType
         read(fid) Ny,NzEarth
         if((sigma(i)%Ny .NE. Ny).OR. (sigma(i)%NzEarth .NE. NzEarth)) then
            close(fid)
            call errStop('Size of cond does not agree with contents of file in readAll_Cond2D')
         else
            read(fid) sigma(i)%v
            read(fid) sigma(i)%AirCond
         endif
      enddo
      close(fid)

      end subroutine readAll_Cond2D
      
end module modelparameter
