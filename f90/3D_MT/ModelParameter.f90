!***********************************************************************
!   module to define modelparam abstract data type,
!    and all methods required for creation, algebra, and inner products,
!    including covariances.  Attributes of the data type
!    are private to this module to force coding of other modules
!    to be completely independent of the specific parameterization
!    instance/implementation
!
! Defines: modelParam
! Uses: EMfield, Grid3D, UTILS

module ModelParameter

  use grid3d
  use math_constants
  use utilities
  use sg_scalar
  use sg_vector
  use sg_sparse_vector

  implicit none

  ! supported model parameter types (conductivity only)
   character(len=80), parameter		:: LOGE = 'LOGE'
   character(len=80), parameter		:: LINEAR = 'LINEAR'
   
  type :: modelParam_t
     !  Parameters which define conductivity distribution.
     !   for the present approach, with each cell of the computational grid
     !   allowed a separate (constant) conductivity) a separate type here
     !   is hardly needed ...
     !  The idea is that the conductivity parameter should be treated as an
     !   "abstract data type", defined along with a set of routines for
     !    mapping to the internal representation of conductivity used directly
     !    by the forward solver.
     !  HERE we are starting to implement this idea ... cond_param will be
     !   public, but all attributes will be private.  Only routines in this
     !   module can use the actual attributes that define the specific realization
     !   of the model parameterization
     private
     integer			:: Nx,Ny,NzEarth
     type(rscalar)		:: cellCond
     type(grid3d_t),pointer    :: grid
     real (kind=selectedPrec)   :: AirCond
     logical			:: allocated = .false.
     character (len=80)		:: paramType = ''
     !  supported paramType at present: LINEAR and LOGE
  end type modelParam_t

!   ModelParameter module needs to contain routines of three general types
!   1)   routines which define basic model parameter methods or operations:
!       creation, destruction, zeroing, copying, linear algebra, dot products

      interface assignment (=)
         MODULE PROCEDURE copy_modelParam
      end interface

!      interface operator (.dot.)
!         MODULE PROCEDURE dotProd_modelParam
!      end interface

      public  ::      create_modelParam,deall_modelParam,dotProd_modelParam
      public  ::      linComb_modelParam,zero_modelParam,copy_modelParam

!   2) routines which define mappings between the "natural" representation
!        of conductivity/resistivity on the model grid and the formal
!        model parameter structure: in general (allowing for a non-linear
!        mapping, as will occur when log condutivity is used for inversion;
!        in the implementation here,we combine the linear and non-linear 
!        mappings in a single routine)
!        one needs the non-linear mapping, a linearized mapping, and the
!        adjoint of the linearized mapping.  If data functionals depend
!        on conductivity or resistivity additional routines are needed--
!        i.e., a real function (SigC here), with input arguments a conductivity 
!        parameter and cell indicies, which returns the conductivity of 
!        this cell (or resistivity), as well as what is essentially the
!        adjoint of the linearization of this function (sort of ... 
!         see below; this is called QtoModelParam here.)
!       Note that the functionality is always needed, but names can be
!        flexible ... these routines are called by the interpolation
!        routines in InterpEB* (for the latter two only), by the forward
!        modeling codes (the non-linear mapping only), and by SensPDEcoeff
!        (the linearized mappings and adjoint) and thus names/interfaces
!        need to be kept consistent with what is used in these routines.
!        Note that for 2D MT different mapping routines are used for
!        TE and TM modes, as the natural representations for the two cases
!        are different (conductivity vs. resistivity)
!        
      public  ::  modelParamToEdge,edgeToModelParam,SigC
      public  :: getSize_modelParam, setType_modelParam
      public  :: setValue_modelParam, getValue_modelParam


!   3) model covariance opeartors; not implemented yet, so no documenetation!

!      public  ::  modelCov,modelCovInit,modelCovDeall

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
  public          ::  WScov

!   4) read and write routines (and/or some sort of routines for setting
!       model parameters from input arrays); since attributes are private,
!       these routines must be in this module.
       public	:: writeAll_Cond3D
       
Contains
!
!**********************************************************************
!
   !  create_modelParam allocates and initializes arrays for
   !   conductivity parameter structure
   !   Pass grid of type grid3d_t to set array sizes
   ! optional arguments v and vAir are assumed to be consistent
   ! with paramType
   subroutine create_modelParam(grid,paramtype,m,v,vAir)

     implicit none
     type (grid3d_t), intent(in), target	:: grid
     character(80), intent(in)			    :: paramtype
     type (modelParam_t), intent(inout)		:: m
     type (rscalar), intent(in), optional   :: v
     real (kind=selectedPrec), intent(in), optional :: vAir

     if(m%allocated) then
        call deall_modelParam(m)
     endif

     call create_rscalar(grid,m%cellCond,CELL_EARTH)
     m%Nx = grid%Nx
     m%Ny = grid%Ny
     m%NzEarth = grid%NzEarth
     m%grid => grid
     m%paramType = paramtype
     m%allocated = .true.

     if(present(v)) then
        m%cellCond = v
     endif
     
     if(present(vAir)) then
        m%AirCond = vAir
     else
		if(paramtype .eq. LOGE) then
		   m%AirCond = log(SIGMA_AIR)
		else
		   m%AirCond = SIGMA_AIR
		endif     
     endif

   end subroutine create_modelParam

  !************************************************************
   subroutine deall_modelParam(m)
     implicit none
     type (modelParam_t), intent(inout)   :: m

     if(m%allocated) then
        call deall_rscalar(m%cellCond)
        nullify(m%grid)
        m%allocated = .false.
        m%paramType = ''
     endif

   end subroutine deall_modelParam

   !**********************************************************************
   ! Converts the input model parameter structure to paramType, by
   ! comparing paramType with m%paramType and performing the necessary
   ! computations if the two strings differ; assumes that m is allocated.
   subroutine setType_modelParam(m,paramType)

     type(modelParam_t), intent(inout)    :: m
     character(*), intent(in)		      :: paramType

     if(.not.(m%allocated)) then
        call errstop('modelParam must be allocated before calling setType_modelParam')
     endif
     
	 if(trim(paramType) .eq. trim(m%paramType)) then
	    ! we are done
	 else if((paramType == LOGE) .and. (m%paramType == LINEAR)) then
	    ! convert to log
	    m%cellCond%v = log(m%cellCond%v)
	    m%AirCond=log(m%AirCond)
	 else if((paramType == LINEAR) .and. (m%paramType == 'LOGE')) then
	    ! convert to cell
	    m%cellCond%v = exp(m%cellCond%v)
	    m%AirCond=exp(m%AirCond)
	 else
        call errstop('unknown paramType in setType_modelParam')	      
     endif

     m%paramType = paramType
     return

   end subroutine setType_modelParam

!**********************************************************************
   function dotProd_modelParam(m1,m2) result(r)

   !   dot product of two model space parameter objects
   !    here implemented for 3D numerical grid
     real(kind=selectedPrec)          :: r
     type(modelParam_t), intent(in)       :: m1,m2

     ! local variables
     integer    :: j,k

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx).or. &
		 (m1%NzEarth .ne. m2%NzEarth)) then
        write(6,*) 'Nx, Ny, Nz ',m1%Nx, m2%Nx,    &
			m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop('size of m1, m2 incompatable in dotProd_modelParam')
     endif

     r = dotProd_rscalar_f(m1%cellCond, m2%cellCond)

   end function dotProd_modelParam
   
!**********************************************************************

   subroutine zero_modelParam(m)

     !  zeros a model space object

     type(modelParam_t), intent(inout)    :: m

     call zero_rscalar(m%cellCond)

   end subroutine zero_modelParam
 
!**********************************************************************

   subroutine linComb_modelParam(a1,m1,a2,m2,m)

     !  forms the linear combination of model parameters
     !    m = a1*m1 + a2*m2
     !  where a1 and a2 are real constants and m1 and m2
     !   are model parameters
     !   output m may overwrite m1 or m2

     real(kind=selectedPrec),intent(in)       :: a1,a2
     type(modelParam_t), intent(in)		:: m1,m2
     type(modelParam_t), intent(inout)		:: m

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
		(m1%NzEarth .ne. m2%NzEarth)) then
        call errStop('size of m1, m2 incompatable in linComb_modelParam')
     endif
     if(m1%paramType .ne. m2%paramType) then
        call errStop('paramType incompatable in linComb_modelParam')
     endif

     ! make sure m is allocated, and of same size; if not same
     !   size, deallocate and allocate as correct size; otherwise
     !   do nothing (use m as input ... allowing m to overwrite m1)
     if(m%allocated) then
        if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
		(m1%NzEarth .ne. m2%NzEarth)) then
           call deall_modelParam(m)
           call create_modelParam(m1%grid,m1%paramtype,m)
        endif
     else
        call create_modelParam(m1%grid,m1%paramtype,m)
     endif
     m%cellCond%v = a1*m1%cellCond%v + a2*m2%cellCond%v
     !  we are implicitly assuming that if m1 and m2 are the same
     !   size other parameters are of the same type
     m%AirCond = m1%AirCond

   end subroutine linComb_modelParam
   
   !**********************************************************************
   ! Sets cell conductivities in a model parameter object m
   !
   ! the values of v and vAir are determined by paramType;
   ! if different from that of the model parameter, returns
   ! an error. To avoid this error, first convert m to the
   ! required paramType using setType_modelParam.   
   subroutine setValue_modelParam(m,paramType,v,vAir)

     type(modelParam_t), intent(inout)    :: m
     character(80), intent(in)		      :: paramType
     type(rscalar), intent(in)		      :: v
     real(kind=selectedPrec), intent(in), optional :: vAir

     if(.not.(m%allocated)) then
        call errstop('output modelParam must be allocated before calling setValue_modelParam')
     endif
     
     !  error checking
     if((m%Ny .ne. v%Ny).or. (m%Nx .ne. v%Nx) .or. (m%NzEarth .ne. v%Nz)) then
        call errstop('modelParam/rscalar dimensions disagree in setValue_modelParam')
     else if(paramType .ne. m%paramType) then
        call errstop('paramTypes not consistent in setValue_modelParam')
     endif
     
	 ! set values
	 m%cellCond = v
	 if(present(vAir)) then
	    m%AirCond=vAir
	 endif

   end subroutine setValue_modelParam

   !**********************************************************************
   ! Gets cell conductivities from a model parameter object m
   !
   ! Extracts the values of v and vAir;
   ! values that are extracted are converted to paramType.
   ! Different from modelParamToCellCond in that the value
   ! that gets extracted is exactly what is stored in the modelParam:
   ! it does not contain any air layers. This is needed for BC_x0_WS.
   subroutine getValue_modelParam(m,paramType,v,vAir)

     type(modelParam_t), intent(in)       :: m
     character(80), intent(inout)		  :: paramType
     type(rscalar), intent(out)		      :: v
     real(kind=selectedPrec), intent(out), optional :: vAir
     ! local variable
     type(modelParam_t)                   :: mTemp

     if(.not.(m%allocated)) then
        call errstop('input modelParam must be allocated before calling getValue_modelParam')
     endif
     
     if(trim(paramType) .eq. '') then
     	paramType = m%paramType
     endif
     
     if (v%allocated) then
        call deall_rscalar(v)
     endif
     
     ! create a temporary copy
     mTemp = m
         
     ! convert model to the required type
     call setType_modelParam(mTemp,paramType)
     
	 ! set values
	 v = mTemp%cellCond
	 if(present(vAir)) then
	    vAir = mTemp%AirCond
	 endif
	 
	 ! deallocate temporary model parameter
	 call deall_modelParam(mTemp)

   end subroutine getValue_modelParam

   !**********************************************************************
   subroutine copy_modelParam(mOut,mIn)

     type(modelParam_t), intent(in)       :: mIn
     type(modelParam_t), intent(inout)    :: mOut

     ! if mOut is allocated, check to see if if is of same size as
     !   mIn; if not, deallocate and reallocate as correct size; otherwise
     !   use as input
     if(mOut%allocated) then
        if((mOut%Ny .ne. mIn%Ny).or. (mOut%Nx .ne. mIn%Nx) .or. &
		(mOut%NzEarth .ne. mIn%NzEarth)) then
           call deall_modelParam(mOut)
           call create_modelParam(mIn%grid,mIn%paramtype,mOut)
        endif
     else
        call create_modelParam(mIn%grid,mIn%paramtype,mOut)
     endif
     mOut%cellCond%v = mIn%cellCond%v
     mOut%AirCond = mIn%AirCond

   end subroutine copy_modelParam

   !************************************************************************
   !  getSize_modelParam extracts model size from a modelParam_t variable
   subroutine getSize_modelParam(m,Nx,Ny,NzEarth)

     implicit none
     type (modelParam_t), intent(in)   	  :: m
     integer, intent(out)                 :: Nx,Ny,NzEarth

     if (.not.m%allocated) then
        call errStop('Model parameter not allocated in getValue_modelParam')
     end if

     Nx = m%Nx
     Ny = m%Ny
     NzEarth = m%NzEarth
 
   end subroutine getSize_modelParam
   
   !**********************************************************************
   subroutine modelParamToCellCond(m,cCond,paramType,grid,AirCond)

     type(modelParam_t), intent(in)        :: m
     type(rscalar), intent(inout)	       :: cCond
     character(80), intent(out), optional  :: paramType
     type(grid3d_t), intent(out), optional :: grid
     real(kind=selectedPrec), intent(out), optional :: AirCond

     if(cCond%allocated) then
        if((cCond%Ny .ne. m%Ny).or. (cCond%Nx .ne. m%Nx) .or. &
		(cCond%Nz .ne. m%grid%Nz)) then
           call deall_rscalar(cCond)
           call create_rscalar(m%grid,cCond,CENTER)
        endif
     else
        call create_rscalar(m%grid,cCond,CENTER)
     endif
     if(m%paramType .eq. LOGE) then
        cCond%v(:,:,1:m%grid%NzAir) = exp(m%AirCond)
        cCond%v(:,:,m%grid%NzAir+1:m%grid%Nz) = exp(m%cellCond%v)
     else
        cCond%v(:,:,1:m%grid%NzAir) = m%AirCond
        cCond%v(:,:,m%grid%NzAir+1:m%grid%Nz) = m%cellCond%v
     endif
     
     if (present(paramType)) then
        paramType = m%paramType
     end if

     if (present(grid)) then
        if (grid%allocated) then
           call deall_grid3d(grid)
        end if
        grid = m%grid
     end if

     if (present(AirCond)) then
        AirCond = m%AirCond
     end if
          
   end subroutine modelParamToCellCond

  !**********************************************************************
  
  subroutine ModelParamToEdge(m, Econd,m0)
  !  Maps conductivity defined on grid cells (prisms) to edge-nodes ... this
  !   can be used for both non-linear and linearized mappings when 
  !   paramtype = LOGE.  In this case, if optional argument m0 is present
  !   the linearized mapping is used; otherwise the non-linear mappring is used.
  !   If paramtype = LINEAR (linear conductivity) then the same linear mapping
  !   is done in whether m0 is present or not (and in fact, m0 is not referenced
  !   in this case)
  !  NOTE: there is a subtlty (and possible source of error) associated with
  !   treatement of air conductivity.  This is (in the present implementation)
  !   viewed as a non-adjustable parameter, but when this routine is called for
  !   the full (i.e., in general non-linear) mapping it must set the edges
  !   in the air to airCond.  On the other hand, when called to implement the
  !   linearization (i.e., to map perturbations in the model parameter) it
  !   should leave air edges set to zero.  WHen the conductivity mapping is
  !   linear, there is no way to distinguish between whether this is a linear
  !   mapping of the perturbation, or the full linear mapping.
  !   One should thus set airCond in perturbations to zero (or they will
  !   be propagated to perturbation in the air edges, only for the linear
  !   model parameter case).

    implicit none
    !  INPUTS:  model parameter 
    type(modelParam_t), intent(in)            :: m 
    !  OUTPUTS: conductivities on edge nodes, as needed for FD calc.
    !    stored as rvector; allocate Econd before calling
    type(rvector), intent(inout)         :: Econd 
    !   OPTIONAL INPUT: background model parameter for linearization
    type(modelParam_t), optional, intent(in)	:: m0

    !   LOCAL Variables
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzE,nZa
    logical		:: linearizedMapping
    ! the prefix V are the volumes S is total over surrounding cells
    !   common sides are ommited from volume calculation
    real(kind=selectedPrec)            :: Vdr,Vdl,Vur,Vul,S
    real(kind=selectedPrec)            :: Vdn,Vds,Vun,Vus
    real(kind=selectedPrec)            :: Vnr,Vnl,Vsr,Vsl,airCond
    real(kind=selectedPrec),pointer,dimension(:,:,:)	:: temp
    
    linearizedMapping = present(m0)

    if ((.not.m%allocated).or.(.not.Econd%allocated)) then
      write(0,*) 'Ccond or Econd in EdgeCond not allocated yet'
      stop
    end if

    ! Could add code to check whether the bounds are the same ...
    
    nx = m%Nx
    ny = m%Ny
    nz = m%grid%nZ
    nZE = m%nZEarth
    nZa = Nz-nZE
    
    allocate(temp(nx,ny,nzE))
    airCond = m%airCond
    temp = m%cellCond%v
   
    if(m%paramType .EQ. LOGE) then 
       if(linearizedMapping) then
          temp = temp*exp(m0%cellCond%v)
          airCond = R_ZERO
       else
          temp = exp(temp)
          airCond = exp(m%airCond)
       endif
    endif
    
    ! Now map onto edges, x, y, then z ... first  zero output
    ! Ex edges: (leave boundary edges set to 0)
    Econd%x = 0.0
    do ix = 1, nx
       do iy = 2, ny
          do iz = Nza+1,Nz
             izE = iz-Nza 
             Vdr = m%grid%dy(iy)*m%grid%dz(iz)
             Vdl = m%grid%dy(iy-1)*m%grid%dz(iz)
             Vur = m%grid%dy(iy)*m%grid%dz(iz-1)
             Vul = m%grid%dy(iy-1)*m%grid%dz(iz-1)
             S = Vdr+Vdl+Vur+Vul
             if(izE.eq.1) then
                Econd%x(ix,iy,iz) = (temp(ix,iy,izE)*Vdr + &
                                  temp(ix,iy-1,izE)*Vdl)/S
             else
                Econd%x(ix,iy,iz) = (temp(ix,iy,izE)*Vdr + &
                                  temp(ix,iy-1,izE)*Vdl + &
                                  temp(ix,iy  ,izE-1)*Vur + &
                                  temp(ix,iy-1,izE-1)*Vul)/S
             endif
          enddo
	  do iz = 1,Nza
             Econd%x(ix,iy,iz) = airCond
          enddo
       enddo
     enddo

     ! Ey edges
     Econd%y = 0.0
     do ix = 2, nx
        do iy = 1, ny
           do iz = Nza+1,Nz
              izE = iz-Nza 
              Vdn = m%grid%dx(ix)*m%grid%dz(iz)
              Vds = m%grid%dx(ix-1)*m%grid%dz(iz)
              Vun = m%grid%dx(ix)*m%grid%dz(iz-1)
              Vus = m%grid%dx(ix-1)*m%grid%dz(iz-1)
              S = Vdn+Vds+Vun+Vus
              if(izE.eq.1) then
                Econd%y(ix,iy,iz) = (temp(ix,iy,izE)*Vdn + &
                                  temp(ix-1,iy,izE)*Vds)/S
              else
                Econd%y(ix,iy,iz) = (temp(ix,iy,izE)*Vdn + &
                                  temp(ix-1,iy,izE)*Vds + &
                                  temp(ix,iy  ,izE-1)*Vun + &
                                  temp(ix-1,iy,izE-1)*Vus)/S
              endif
           enddo
	   do iz = 1,Nza
              Econd%y(ix,iy,iz) = airCond
           enddo
        enddo 
     enddo
      
     ! Ez edges
     Econd%z = 0.0
     do ix = 2, nx
        do iy = 2, ny
           do iz = Nza+1,Nz
              izE = iz-Nza 
              Vnr = m%grid%dx(ix)*m%grid%dy(iy)
              Vnl = m%grid%dx(ix)*m%grid%dy(iy-1)
              Vsr = m%grid%dx(ix-1)*m%grid%dy(iy)
              Vsl = m%grid%dx(ix-1)*m%grid%dy(iy-1)
              S =  Vnr+Vnl+Vsr+Vsl
              Econd%z(ix,iy,iz) = (temp(ix,iy,izE)*Vnr + &
                                  temp(ix,iy-1,izE)*Vnl + &
                                  temp(ix-1,iy  ,izE)*Vsr + &
                                  temp(ix-1,iy-1,izE)*Vsl)/S
           enddo
	   do iz = 1,Nza
              Econd%z(ix,iy,iz) = airCond
           enddo
        enddo
      enddo
      
      deallocate(temp)

  end subroutine ModelParamToEdge
  
   !**********************************************************************
  
  subroutine EdgeToModelParam(Econd,m,m0)
  !  Maps from a real vector (defined on edges) to modelParam;
  !    the adjoint of linear mapping implemented by ModelParamToEdge 
  !  
    implicit none
    !  INPUTS:  real vector defined on edges
    type(rvector), intent(in)            :: Econd 
    !  OUTPUTS: model parameter
    type(modelParam_t), intent(inout)         :: m 
    !  INPUT (OPTIONAL) background model parameter, 
    !         required if m%paramtype=LOGE
    type(modelParam_t), optional, intent(in)	:: m0

    !   LOCAL Variables
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzE,nZa
    ! the prefix V are the volumes S is total over surrounding cells
    !   common sides are ommited from volume calculation
    real(kind=selectedPrec)            :: Vdr,Vdl,Vur,Vul,S
    real(kind=selectedPrec)            :: Vdn,Vds,Vun,Vus
    real(kind=selectedPrec)            :: Vnr,Vnl,Vsr,Vsl
    real(kind=selectedPrec),pointer,dimension(:,:,:)	:: temp  

    
    if ((.not.m%allocated).or.(.not.Econd%allocated)) then
      call errStop('m or Econd not allocated yet in EdgeToModelParam')
      stop
    end if
    
    if((m%paramtype.eq.LOGE).and. (.not.present(m0))) then
       call errStop('Background conductivity required for paramType LOGE in EdgeToModelParam')
    endif

    ! Could add code to check whether the bounds are the same ...
    
    nx = m%Nx
    ny = m%Ny
    nz = m%grid%nZ
    nZE = m%nZEarth
    nZa = Nz-nZE
    
    allocate(temp(nx,ny,nzE))
    temp = R_ZERO 
    
    ! Now map onto edges, x, y, then z ... first  zero output
    ! Ex edges: (leave boundary edges set to 0)
    do ix = 1, nx
       do iy = 2, ny
          do iz = Nza+1,Nz
             izE = iz-Nza 
             Vdr = m%grid%dy(iy)*m%grid%dz(iz)
             Vdl = m%grid%dy(iy-1)*m%grid%dz(iz)
             Vur = m%grid%dy(iy)*m%grid%dz(iz-1)
             Vul = m%grid%dy(iy-1)*m%grid%dz(iz-1)
             S = Vdr+Vdl+Vur+Vul
             if(izE.gt.1) then
                temp(ix,iy,izE-1) = temp(ix,iy,izE-1)+ Econd%x(ix,iy,iz)*Vur/S
                temp(ix,iy-1,izE-1) = temp(ix,iy-1,izE-1)+Econd%x(ix,iy,iz)*Vul/S
             endif
             temp(ix,iy,izE) = temp(ix,iy,izE)+ Econd%x(ix,iy,iz)*Vdr/S
             temp(ix,iy-1,izE) = temp(ix,iy-1,izE)+Econd%x(ix,iy,iz)*Vdl/S
          enddo
       enddo
     enddo

     ! Ey edges
     do ix = 2, nx
        do iy = 1, ny
           do iz = Nza+1,Nz
              izE = iz-Nza 
              Vdn = m%grid%dx(ix)*m%grid%dz(iz)
              Vds = m%grid%dx(ix-1)*m%grid%dz(iz)
              Vun = m%grid%dx(ix)*m%grid%dz(iz-1)
              Vus = m%grid%dx(ix-1)*m%grid%dz(iz-1)
              S = Vdn+Vds+Vun+Vus
              if(izE.gt.1) then
	        temp(ix,iy,izE-1) = temp(ix,iy,izE-1)+ Econd%y(ix,iy,iz)*Vun/S
                temp(ix-1,iy,izE-1) = temp(ix-1,iy,izE-1)+Econd%y(ix,iy,iz)*Vus/S
              endif
              temp(ix,iy,izE) = temp(ix,iy,izE)+ Econd%y(ix,iy,iz)*Vdn/S
              temp(ix-1,iy,izE) = temp(ix-1,iy,izE)+Econd%y(ix,iy,iz)*Vds/S
           enddo
        enddo
      enddo
      
     ! Ez edges
     do ix = 2, nx
        do iy = 2, ny
           do iz = Nza+1,Nz
              izE = iz-Nza 
              Vnr = m%grid%dx(ix)*m%grid%dy(iy)
              Vnl = m%grid%dx(ix)*m%grid%dy(iy-1)
              Vsr = m%grid%dx(ix-1)*m%grid%dy(iy)
              Vsl = m%grid%dx(ix-1)*m%grid%dy(iy-1)
              S =  Vnr+Vnl+Vsr+Vsl
              temp(ix,iy,izE) = temp(ix,iy,izE)+ Econd%z(ix,iy,iz)*Vnr/S
              temp(ix,iy-1,izE) = temp(ix,iy-1,izE)+Econd%z(ix,iy,iz)*Vnl/S
              temp(ix-1,iy,izE) = temp(ix-1,iy,izE)+ Econd%z(ix,iy,iz)*Vsr/S
              temp(ix-1,iy-1,izE) = temp(ix-1,iy-1,izE)+Econd%z(ix,iy,iz)*Vsl/S
           enddo
        enddo
      enddo  
  
    if(m%paramType .EQ. LOGE) then 
       m%cellCond%v = temp*exp(m0%cellCond%v)
    else
       m%cellCond%v = temp
    endif
    
    deallocate(temp)

  end subroutine EdgeToModelParam

!**********************************************************************
   subroutine QtoModelParam(Qj,sigma0,dsigmaReal,dSigmaImag)

   !  THIS IS NOT EVEN LOOKED AT WITH REGARD TO CORRECTNESS OF
   !   FORMULAE FOR 3D!!!!! JUST COPIED FROM 2D-TM !!!!!!

   !  given input sparse vector defined on cells,
   !   compute vector defined on parameter space
   !   add result to inputs dsigmaReal, dsigmaImag
   !    Q_j^T \partial\rho_cell/\partial\theta

   type (sparseVecC), intent(in)   :: Qj
   type (modelParam_t), intent(in) :: sigma0
   type (modelParam_t), intent(inout)    :: dsigmaReal
   type (modelParam_t), intent(inout),optional    :: dsigmaImag

   ! local variables
   integer              :: jk,i,j,k

   do jk = 1,Qj%nCoeff
      i = Qj%I(jk)
      j = Qj%J(jk)
      k = Qj%K(jk)
      if(sigma0%paramType .eq. LOGE) then
         dsigmaReal%CellCond%v(i,j,k) = dSigmaReal%CellCond%v(i,j,k) &
                        -real(Qj%C(jk))/exp(sigma0%CellCond%v(i,j,k))
         if(present(dSigmaImag)) then
            dsigmaImag%CellCond%v(i,j,k) = dSigmaImag%CellCond%v(i,j,k) &
                        -imag(Qj%C(jk))/exp(sigma0%CellCond%v(i,j,k))
         endif
      else
         dsigmaReal%CellCond%v(i,j,k) = dSigmaReal%CellCond%v(i,j,k) &
           -real(Qj%C(jk))/(sigma0%CellCond%v(i,j,k)*sigma0%CellCond%v(i,j,k))
         if(present(dSigmaImag)) then
            dsigmaImag%CellCond%v(i,j,k) = dSigmaImag%CellCond%v(i,j,k) &
             -imag(Qj%C(jk))/(sigma0%CellCond%v(i,j,k)*sigma0%CellCond%v(i,j,k))
         endif
      endif
   enddo

   end subroutine QtoModelParam

  !**********************************************************************
  function sigC(m,xyz,ix,iy,iz) result(r)
  !   computes conductivity for edge xyz/i,j,k using input modelParam
  !    structure m.  This function defines how the abstract
  !     conductivity parameter is mapped to edge conductivities needed
  !     for more accurate electric field interpolation.  The derivative
  !     of this function is required for evaluation of linearized
  !     data functionals, and for construction of the direct
  !     parameter space part of the comb

    type (modelParam_t), intent(in)     :: m
    integer, intent(in)               :: xyz,ix,iy,iz
    real(kind=selectedPrec)           :: r
 
    ! local variables
    integer		:: izE,nx,ny,nz,nZE,nZa
    real(kind=selectedPrec)           :: w11,w21,w12,w22,S,temp(2,2)

    nx = m%Nx
    ny = m%Ny
    nz = m%grid%nZ
    nZE = m%nZEarth
    nZa = Nz-nZE
 
    izE = iz-Nza 
    selectcase(xyz)
    case(1)		 ! Ex edge
       w11 = m%grid%dy(iy)*m%grid%dz(iz)
       temp(1,1) = m%cellCond%v(ix,iy,izE)
       w21 = m%grid%dy(iy-1)*m%grid%dz(iz)
       temp(2,1) = m%cellCond%v(ix,iy-1,izE)
       w12 = m%grid%dy(iy)*m%grid%dz(iz-1)
       temp(1,2) = m%cellCond%v(ix,iy,izE-1)
       w22 = m%grid%dy(iy-1)*m%grid%dz(iz-1)
       temp(2,2) = m%cellCond%v(ix,iy-1,izE-1)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       if(izE.eq.1) then
          r = (temp(1,1)*w11+temp(1,2)*w21)/S
       else
          r = (temp(1,1)*w11 + temp(2,1)*w21+    &
               temp(1,2)*w12 + temp(2,2)*w22)/S
       endif

    case(2)		 ! Ey edge
       w11 = m%grid%dx(ix)*m%grid%dz(iz)
       temp(1,1) = m%cellCond%v(ix,iy,izE)
       w21 = m%grid%dx(ix-1)*m%grid%dz(iz)
       temp(2,1) = m%cellCond%v(ix-1,iy,izE)
       w12 = m%grid%dx(ix)*m%grid%dz(iz-1)
       temp(1,2) = m%cellCond%v(ix,iy,izE-1)
       w22 = m%grid%dx(ix-1)*m%grid%dz(iz-1)
       temp(2,2) = m%cellCond%v(ix-1,iy,izE-1)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       if(izE.eq.1) then
          r = (temp(1,1)*w11+temp(1,2)*w21)/S
       else
          r = (temp(1,1)*w11 + temp(2,1)*w21+ &
               temp(1,2)*w12 + temp(2,2)*w22)/S
       endif

    case(3)		 ! Ez edge
       w11 = m%grid%dx(ix)*m%grid%dy(iy)
       temp(1,1) = m%cellCond%v(ix,iy,izE)
       w21 = m%grid%dx(ix)*m%grid%dy(iy-1)
       temp(2,1) = m%cellCond%v(ix,iy-1,izE)
       w12 = m%grid%dx(ix-1)*m%grid%dy(iy)
       temp(1,2) = m%cellCond%v(ix-1,iy,izE)
       w22 = m%grid%dx(ix-1)*m%grid%dy(iy-1)
       temp(2,2) = m%cellCond%v(ix-1,iy-1,izE)
       S = w11+w21+w12+w22
       if(m%paramtype .eq. LOGE) then
          temp = exp(temp)
       endif
       r = (temp(1,1)*w11 + temp(2,1)*w21+  &
            temp(1,2)*w12 + temp(2,2)*w22)/S

    endselect
  end function sigC

     !******************************************************************
      subroutine writeAll_Cond3D(fid,cfile,nSigma,m,header)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file  

      integer, intent(in)               :: fid,nSigma
      character(*), intent(in)          :: cfile, header
      type(modelParam_t), intent(in)      :: m(nSigma)

      integer 		:: i,Nz,NzAir

      open(unit=fid, file=cfile, form='unformatted')
      write(fid) header
      write(fid) nSigma
   
      do i = 1,nSigma
         Nz = m(i)%grid%Nz
         NzAir = m(i)%grid%Nz-m(i)%NzEarth
         write(fid) m(i)%paramType
         write(fid) m(i)%Nx,m(i)%Ny,m(i)%NzEarth
         write(fid) m(i)%grid%dx
         write(fid) m(i)%grid%dy
         write(fid) m(i)%grid%dz(NzAir+1:Nz)
         write(fid) m(i)%AirCond
         write(fid) m(i)%CellCond%v
      enddo
      close(fid)
      end subroutine writeAll_Cond3D

     !******************************************************************
      subroutine readAll_Cond3D(fid,cfile,nSigma,m,header)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file  

      integer, intent(in)               :: fid
      integer, intent(in)              :: nSigma
      character(*), intent(in)          :: cfile
      character(80), intent(out)        :: header
      type(modelParam_t), intent(inout)   :: m(nSigma)

      integer 		:: i,nS,Nx,Ny,Nz,NzAir,NzEarth,istat

      open(unit=fid, file=cfile, form='unformatted')
      read(fid) header
      read(fid) nS
      
      if(nS .NE. nSigma) then
          call errStop('size of sigma does not agree with contents of file in readAll_Cond2D')
      endif

      do i = 1,nSigma
         read(fid) m(i)%paramType
         read(fid) Nx,Ny,NzEarth
         if((m(i)%Nx .NE. Nx) .OR. (m(i)%Ny .NE. Ny) .OR. (m(i)%NzEarth .NE. NzEarth)) then
            close(fid)
            call errStop('Size of cond does not agree with contents of file in readAll_Cond2D')
         else
            Nz = m(i)%grid%Nz
            NzAir = m(i)%grid%Nz-m(i)%NzEarth
            read(fid) m(i)%grid%dx
            read(fid) m(i)%grid%dy
            read(fid) m(i)%grid%dz(NzAir+1:Nz)
            read(fid) m(i)%AirCond
            read(fid) m(i)%CellCond%v
         endif
      enddo
      close(fid)
      end subroutine readAll_Cond3D
!**********************************************************************

  subroutine WScov(m)

  type (modelParam_t), intent(inout) :: m

    !call setup1DCM(m%grid)
	!call solveDiff(m)

  end subroutine WScov

end module modelParameter
