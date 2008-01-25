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
     !  supported paramType at present: 'CELL_EARTH' and 'LOG_CELL'
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
       public	:: read_Cond3D_RM, write_Cond3D, read_Cond3D, &
		writeAll_Cond3D
Contains
!
!**********************************************************************
!
   !  create_modelParam allocates and initializes arrays for
   !   conductivity parameter structure
   !   Pass grid of type grid3d_t to set array sizes
   subroutine create_modelParam(grid,paramtype,m)

     implicit none
     type (grid3d_t), intent(in), target	:: grid
     character*80, intent(in)			:: paramtype
     type (modelParam_t), intent(inout)		:: m

     call create_rscalar(grid,m%cellCond,CELL_EARTH)
     m%Nx = grid%Nx
     m%Ny = grid%Ny
     m%NzEarth = grid%NzEarth
     m%grid => grid
     m%paramType = paramtype
     m%allocated = .true.

     if(paramtype.eq.LOG_CELL) then
        m%AirCond = log(SIGMA_AIR)
     else
        m%AirCond = SIGMA_AIR
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
   function dotProd_modelParam(m1,m2) result(r)

   !   dot product of two model space parameter objects
   !    here implemented for 3D numerical grid
     real(kind=selectedPrec)          :: r
     type(modelParam_t), intent(in)       :: m1,m2

     ! local variables
     integer    :: j,k
     character*80       :: msg

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx).or. &
		 (m1%NzEarth .ne. m2%NzEarth)) then
        msg = 'Error: size of m1, m2 incompatable in dotProd_modelParam'
        write(6,*) 'Nx, Ny, Nz ',m1%Nx, m2%Nx,    &
			m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop(msg)
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

     !  local variables
     character*80                               :: msg

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
		(m1%NzEarth .ne. m2%NzEarth)) then
        msg = 'Error: size of m1, m2 incompatable in linComb_modelParam'
        call errStop(msg)
     endif
     if(m1%paramType .ne. m2%paramType) then
        msg = 'Error: paramType incompatable in linComb_modelParam'
        call errStop(msg)
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
   subroutine set_modelParam(CellCond,paramType,m,AirCond)

     ! set cell conductivities in a model parameter object m
     !   using CONDUCTIVITY values (NOT Log, NOT resistivity)
     !   from an rscalar object CellCond
     !
     !   Here paramType can be used to control log/linear 
     !      conductivity in the output modelParam object.
     !   
     !   AirCond is an optional air conductivity (NOT Log, NOT
     !    resistivity.  If ommitted default air conductivity is set
     !
     type(rscalar), intent(in)		:: CellCond
     character*80, intent(in)		:: paramType
     type(modelParam_t), intent(inout)    :: m
     real(kind=selectedPrec),optional	:: AirCond

     character*80			::msg

     !  error checking
     if(m%allocated) then
        if((m%Ny .ne. CellCond%Ny).or. (m%Nx .ne. CellCond%Nx) .or. &
		(m%NzEarth .ne. CellCond%Nz)) then
           msg = 'modelParam/rscalar dimensions disagree in set_ModelParam'
           call errstop(msg)
        else
           m%cellCond = CellCond
           if(present(AirCond)) then
              m%AirCond=AirCond
           else
              m%AirCond=SIGMA_AIR
           endif
           if(paramType(1:8).eq.LOG_CELL) then
              m%cellCond%v = log(m%cellCond%v)
              m%AirCond = log(m%AirCond)
           endif
        endif
     else
        msg = 'output modelParam must be allocated before calling set_modelParam'
        call errstop(msg)
     endif
   end subroutine set_modelParam

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

   !**********************************************************************
   subroutine modelParamToCellCond(m,cCond)

     type(modelParam_t), intent(in)       :: m
     type(rscalar), intent(inout)	:: cCond

     if(cCond%allocated) then
        if((cCond%Ny .ne. m%Ny).or. (cCond%Nx .ne. m%Nx) .or. &
		(cCond%Nz .ne. m%grid%Nz)) then
           call deall_rscalar(cCond)
           call create_rscalar(m%grid,cCond,CENTER)
        endif
     else
        call create_rscalar(m%grid,cCond,CENTER)
     endif
     if(m%paramType(1:8) .eq. LOG_CELL) then
        cCond%v(:,:,1:m%grid%NzAir) = exp(m%AirCond)
        cCond%v(:,:,m%grid%NzAir+1:m%grid%Nz) = exp(m%cellCond%v)
     else
        cCond%v(:,:,1:m%grid%NzAir) = m%AirCond
        cCond%v(:,:,m%grid%NzAir+1:m%grid%Nz) = m%cellCond%v
     endif

   end subroutine modelParamToCellCond

  !**********************************************************************
  
  subroutine ModelParamToEdge(m, Econd,m0)
  !  Maps conductivity defined on grid cells (prisms) to edge-nodes ... this
  !   can be used for both non-linear and linearized mappings when 
  !   paramtype = 'LOG_CELL'.  In this case, if optional argument m0 is present
  !   the linearized mapping is used; otherwise the non-linear mappring is used.
  !   If paramtype = 'COND' (linear conductivity) then the same linear mapping
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
   
    if(m%paramType(1:8) .EQ. LOG_CELL) then 
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
    !         required if m%paramtype=LOG_CELL
    type(modelParam_t), optional, intent(in)	:: m0

    !   LOCAL Variables
    integer		:: ix,iy,iz,ize,nx,ny,nz,nzE,nZa
    ! the prefix V are the volumes S is total over surrounding cells
    !   common sides are ommited from volume calculation
    real(kind=selectedPrec)            :: Vdr,Vdl,Vur,Vul,S
    real(kind=selectedPrec)            :: Vdn,Vds,Vun,Vus
    real(kind=selectedPrec)            :: Vnr,Vnl,Vsr,Vsl
    real(kind=selectedPrec),pointer,dimension(:,:,:)	:: temp  
    character*80 msg

    
    if ((.not.m%allocated).or.(.not.Econd%allocated)) then
      write(0,*) 'Ccond or Econd in EdgeCond not allocated yet'
      stop
    end if
    
    if((m%paramtype.eq.LOG_CELL).and. (.not.present(m0))) then
       msg = 'Background conductivity required for paramType LOG_CELL'
       call errStop(msg)
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
  
    if(m%paramType(1:8) .EQ. LOG_CELL) then 
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
      if(sigma0%paramType(1:8) .eq. LOG_CELL) then
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
       if(m%paramtype .eq. LOG_CELL) then
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
       if(m%paramtype .eq. LOG_CELL) then
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
       if(m%paramtype .eq. LOG_CELL) then
          temp = exp(temp)
       endif
       r = (temp(1,1)*w11 + temp(2,1)*w21+  &
            temp(1,2)*w12 + temp(2,2)*w22)/S

    endselect
  end function sigC

  !******************************************************************
   subroutine read_Cond3D_RM(fid,cfile,paramType,m)

  !   opens cfile on unit fid, read resistivities
  !   from Mackie format 3D model file, and returns
  !   as a paramType (loge or linear)  modelParam
  !
  !  This is a "wrapper" routine that calls a reading routine
  !   defined outside of the scope of module modelParameter to set
  !   attributes of a modelParam object
  !
  !   output modelParam is deallocated (if necessary)
  !   and (re)created based on grid defined in cfile;
  !   air conductivity is set to default (could add another
  !   argument to surboutine to change this ... airCond is
  !    not given in Mackie file)

    integer, intent(in)			:: fid
    character*80, intent(in)  		:: cfile,paramType
    type(modelParam_t), intent(inout)	:: m

    !  local variables
    type(grid3d_t)			:: tempGrid
    type(rscalar)			:: tempCond

       if(m%allocated) then
          call deall_modelParam(m)
       endif

       call readRMgridCond(fid,cfile,tempGrid,tempCond)

       if(paramType(1:8) == LOG_CELL) then
          tempCond%v = log(tempCond%v)
       endif

       call create_modelParam(tempGrid,paramtype,m)
       call copy_rscalar(m%cellCond,tempCond)

       call deall_rscalar(tempCond)
       call deall_grid3D(tempGrid)

    end subroutine read_Cond3D_RM
  !******************************************************************
   subroutine write_Cond3D(fid,cfile,m)

   !  open cfile on unit fid, writes out object of
   !   type modelParam in standard *binary* format (comparable
   !   format written by write_Cond3D), then close file

      integer, intent(in)               :: fid
      character*80, intent(in)  :: cfile
      type(modelParam_t), intent(in)              :: m

      integer		:: NzAir,Nz,Nx,Ny

      Nz = m%grid%Nz
      Nx = m%Nx
      Ny = m%Ny
      NzAir = m%grid%Nz-m%NzEarth
      open(unit=fid, file=cfile, form='unformatted')
      write(fid) m%paramType
      write(fid) Nx,Ny,m%NzEarth
      write(fid) m%grid%dx
      write(fid) m%grid%dy
      write(fid) m%grid%dz(NzAir+1:Nz)
      write(fid) m%AirCond
      write(fid) m%CellCond%v
      close(fid)
      end subroutine write_Cond3D
  !******************************************************************
   subroutine read_Cond3D(fid,cfile,m)

   !  open cfile on unit fid, writes out object of
   !   type modelParam in standard *binary* format (comparable
   !   format written by write_Cond3D), then close file

      integer, intent(in)                :: fid
      character*80, intent(in) 		 :: cfile
      type(modelParam_t), intent(inout)    :: m

      integer		:: NzAir,Nz,Nx,Ny,NzEarth
      character*80	:: paramType,msg
      real(kind=selectedPrec) 	:: AirCond

      if(m%allocated) then
          open(unit=fid, file=cfile, form='unformatted',status='old')
     
          read(fid) paramType
          read(fid) Nx,Ny,NzEarth
          if((m%Ny .NE. Ny).OR.(m%NzEarth .NE. NzEarth) &
		.or. (m%Nx .NE. Nx)) then
             close(fid)
             msg = 'Size of cond does not agree with contents of file'
             call errStop(msg)
          else
             read(fid)   !dx
             read(fid)   !dy
             read(fid)   ! dz
             read(fid) AirCond
             read(fid) m%cellCond%v
             close(fid)
          endif
       else
          msg = 'modelParam must be allocated before call to read_Cond3D'
          call errStop(msg)
       endif
      end subroutine read_Cond3D
     !******************************************************************
      subroutine writeAll_Cond3D(fid,cfile,header,nSigma,m)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file  

      integer, intent(in)               :: fid,nSigma
      character*80, intent(in)          :: cfile, header
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
  ! ***************************************************************************
  subroutine ReadRMgridCond(fidRM,inputFile,grid,Cond)
  ! this routine reads files in Mackie's 3D formats, returning the basic
  !   grid components, and optionally also conductivity
  !   If present, Cond is created during call

    implicit none

    integer,intent(in)                          :: fidRM
    character(len=80), intent(in)               :: inputFile
    type (grid3d_t), intent(inout)             :: grid
    type (rscalar), intent(inout), optional     :: Cond

    real(kind=selectedPrec)                     :: origin(3)
    real(kind=selectedPrec),pointer,dimension(:)    :: res
    integer                                     :: whichLayer
    integer                                     :: ix,iy,iz,ip,i,j
    real(kind=selectedPrec)                     :: alpha = 3.
    integer                                     :: jj,Nx,Ny,Nz,NzEarth,NzAir
    integer                                     :: status, ioerr

    character (len=80)                          :: ifValues = ''
    character (len=80)                          :: someChar = ''
    integer                                     :: jOne, jTwo
    logical                                     :: returnCond

    returnCond = present(Cond)

    ! Open file and read grid
    open(unit=fidRM,file=inputFile,status='old',ERR=9000)
    read(fidRM,*) Nx, Ny, NzEarth, nzAir, ifValues

    if (ifValues(1:6) /= 'VALUES') then
        write(0, *) 'Mapping not supported yet in:ReadGridInputRM'
        stop
    end if

    call create_Grid3D(Nx,Ny,NzAir,NzEarth,grid)

    ! In Randy Mackie's format, dx is read forwwards, as is dy and dz
    read(fidRM,*) (grid%dx(ix),ix=1,grid%nx)
    read(fidRM,*) (grid%dy(iy),iy=1,grid%ny)
    read(fidRM,*) (grid%dz(iz),iz=grid%nzAir+1,grid%nzAir+grid%nzEarth)

    !   Following is Kush's approach to setting air layers:
    ! mirror imaging the dz values in the air layer with respect to
    ! earth layer as far as we can using the following formulation
    ! air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom)
    i = grid%nzAir+1
    j = 0
    do iz = grid%nzAir, 1, -1
        j = j + 1
        grid%dz(iz) = ((alpha)**(j-1))*grid%dz(i)
        i = i + 1
    end do

    ! the topmost air layer has to be atleast 30 km
    if (grid%dz(1).lt.30000) then
        grid%dz(1) = 30000
    end if

    if(returnCond) then
       call create_rscalar(grid,Cond,CELL_EARTH)
    endif

    allocate(res(grid%Nx))
    do iz = 1,grid%nzEarth
       read(fidRM, *) whichLayer
       do iy = 1,grid%ny
          ! in Randy Mackie's format, x varies the fastest
          read(fidRM,*) res
          if(returnCond) then
              Cond%v(:,iy,iz) = 1./res
          endif
       enddo     !iy
    enddo        !iz
    deallocate(res)

    ! skip the three lines: a) WINGLINK, b) site name, and c) block numbers
    ! read WINGLINK (it can also be a blank line).
    read(fidRM, *, IOSTAT = ioerr) someChar
    ! a) WINGLINK
    if (ioerr /= 0) then
        if (someChar(1:8) == 'WINGLINK') then
           write(0, *) 'Model file created by Winglink'
        end if
    end if

    someChar = ''
    ! b) site name
    read(fidRM, *, IOSTAT = ioerr) someChar
    ! c) the block numbers
    read(fidRM, *, IOSTAT = ioerr) jOne, jTwo

    read(fidRM, *, IOSTAT = ioerr) grid%ox, grid%oy
    ! the defualt from a file read through Randy Mackie's format
    ! in Randy Mackie's format, real coordinates are in kilometers
    ! defaults in case of missing data
    if (ioerr /= 0) then
        grid%ox = 0.0
        grid%oy = 0.0
        grid%oz = 0.0
    else
       grid%ox = grid%ox*1000.0
       grid%oy = grid%oy*1000.0
       grid%oz = 0.0
    endif

    read(fidRM, *, IOSTAT = ioerr) grid%rotdeg
    if (ioerr /= 0) then
        grid%rotdeg = 0.0
    end if

    CLOSE(fidRM)

    GOTO 9999

9000 CONTINUE
    WRITE(0,*) '!!! FILE CANNOT BE FOUND !!!'
    STOP

9999 CONTINUE

  end subroutine ReadRMgridCond

!**********************************************************************

  subroutine WScov(m)

  type (modelParam_t), intent(inout) :: m

    !call setup1DCM(m%grid)
	!call solveDiff(m)

  end subroutine WScov

end module modelParameter
