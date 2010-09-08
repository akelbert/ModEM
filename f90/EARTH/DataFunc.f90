! *****************************************************************************
module DataFunc
  ! Module containing the subroutines for a-posteriori analysis of the output
  ! vector h=lH from the routine SolveMaxwells

  use model_operators
  use SolnSpace
  use dataTypes
  use transmitters
  use receivers
  use functionals
  use responses
  implicit none

  public                        :: dataResp, Lrows, Qrows

Contains

!******************************************************************************
  subroutine dataResp(H,m0,iDt,iRx,Resp)
  ! given magnetic field solution and indices into data types and receiver
  ! dictionaries, compute the single complex response function; store output
  ! in a vector or real values

  implicit none
  type (solnVector_t), intent(in)   :: H
  type (modelParam_t), intent(in)   :: m0 ! currently not used
  integer, intent(in)           :: iDt
  integer, intent(in)           :: iRx
  real(kind=prec), intent(inout)    :: Resp(:)

  !  local variables
  complex(8)            :: res
  type (functional_t), pointer   :: dataType
  type (receiver_t), pointer     :: obs
  integer               :: iComp,nComp

  dataType => TFList%info(iDt)
  obs => obsList%info(iRx)

  ! compute the complex data response (C/D ratio)

  select case (dataType%name)
     case('C')
        res = dataResp_rx(C_ratio,obs,H%vec)
     case('D')
        res = dataResp_rx(D_ratio,obs,H%vec)
     case default
        call errStop('Unknown data response specified in dataResp')
  end select

  ! save real output
  Resp(1) = dreal(res)
  Resp(2) = dimag(res)

  end subroutine dataResp

  ! ***************************************************************************
  ! * Lrows is a subroutine to output a full vector g_j defined on edges,
  ! * such that for a single frequency and a single observatory,
  ! * $\pd{psi_{\omega}^j}{veca} = g_j^* \pd{vecH}{veca}$

  subroutine Lrows(H,m0,iDt,iRx,L)

	! uses: grid
    type (solnVector_t), intent(in)                 :: H
    type (modelParam_t), intent(in)                 :: m0   ! not needed?
	integer, intent(in)					            :: iDt,iRx
	type (sparseVector_t), intent(inout)            :: L(:) ! nFunc = 1
	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz
	complex(8)										:: pd_Hx,pd_Hy,pd_Hz
	type (sparsevecc)								:: gc_sparse  ! g*
	type (sparsevecc)                  				:: g_sparse	! g
    type (functional_t), pointer                    :: dataType
    type (receiver_t), pointer                      :: obs
	real(8)											:: EARTH_R

	dataType => TFList%info(iDt)
	obs => obsList%info(iRx)

	if (.not.obs%located) then
	  write(0,*) 'Error: Observatory ',trim(obs%code),' not yet located in Lrows'
	  stop
	  !call LocateReceiver(grid,obs)
	end if

	if (.not.obs%defined) then
	  write(0,*) 'Error: Observatory ',trim(obs%code),' is not defined at Lrows'
	  stop
	end if


	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H%vec) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H%vec)
	Hz = dotProd_noConj(Lz,H%vec)

	EARTH_R = H%grid%z(H%grid%nzAir+1) * m2km

	if (dataType%name == 'C') then

!	  pd_Hx = C_ZERO
!	  pd_Hy = - km2m * (EARTH_R/2) * dtan(obs%colat*d2r) * Hz/(Hy*Hy)
!	  pd_Hz =   km2m * (EARTH_R/2) * dtan(obs%colat*d2r) * 1/Hy
	  pd_Hx = C_ZERO
	  pd_Hy = - km2m * (EARTH_R/2) * Hz/(Hy*Hy)
	  pd_Hz =   km2m * (EARTH_R/2) * 1/Hy

	  call linComb_sparsevecc(Ly,pd_Hy,Lz,pd_Hz,gc_sparse)

	else if (dataType%name == 'D') then

!	  pd_Hx =   km2m * (EARTH_R/2) * dsin(obs%colat*d2r) * 1/Hy
!	  pd_Hy = - km2m * (EARTH_R/2) * dsin(obs%colat*d2r) * Hx/(Hy*Hy)
!	  pd_Hz = C_ZERO
	  pd_Hx =   km2m * (EARTH_R/2) * 1/Hy
	  pd_Hy = - km2m * (EARTH_R/2) * Hx/(Hy*Hy)
	  pd_Hz = C_ZERO

	  call linComb_sparsevecc(Ly,pd_Hy,Lx,pd_Hx,gc_sparse)

	else

	  write(0, *) 'Error: (Lrows) unknown data functional', dataType%name
	  stop

	end if

	! This was g* (conjugated), now we want g
	!g_sparse = conjg(gc_sparse)

	! Generic output
	call create_sparseVector(H%grid,H%tx,L(1))
	!L(1)%L = g_sparse ! left for debugging
    L(1)%L = gc_sparse ! correct for Lrows

	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)
	call deall_sparsevecc(gc_sparse)
    call deall_sparsevecc(g_sparse)

  end subroutine Lrows	! Lrows

!
!****************************************************************************
  subroutine Qrows(e0,m0,iDt,iRx,Qreal,Qimag)
  !  given input background solution vector (e0) and model parameter (Sigma0)
  !  and indices into data type and receiver dictionaries
  !  compute derivative of data functional with respect to model parameters
  !  for all components of the data type ...
  !             (ZERO VECTORS FOR EARTH!!!!)

  type (solnVector_t), intent(in)       :: e0
  type (modelParam_t), intent(in)       :: m0
  integer, intent(in)                   :: iDt, iRx
  !   NOTE: Qreal and Qimag have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE: Qreal and Qimag both exist regardless of whether the data
  !     are real or complex, since Q itself is complex
  type(modelParam_t), intent(inout)     :: Qreal(:), Qimag(:)

  ! local variables
  type(functional_t), pointer :: dataType
  type(receiver_t), pointer   :: obs
  integer       :: istat,nComp,iComp
  logical       :: isComplex

  dataType => TFList%info(iDt)
  obs => obsList%info(iRx)

  ! set the rows of Q to zero
  if((size(Qreal) .ne. dataType%nComp) .or. (size(Qimag) .ne. dataType%nComp)) then
    call errStop('incorrect output size in Qrows')
  endif
  do iComp = 1, dataType%nComp
    Qreal(iComp) = m0
    call zero(Qreal(iComp))
    Qimag(iComp) = m0
    call zero(Qimag(iComp))
  enddo

  end subroutine Qrows

end module DataFunc
