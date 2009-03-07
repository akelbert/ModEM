! *****************************************************************************
module model_operators
	! This module contains the elementary operators for computations involving
	! model parametrization structures. 'Y' stands for spherical harmonics;
	! however, many of the subroutines here would be identical for a different
	! kind of layered parametrization.
	! For all parametrization <-> vector interactions, vectors are constructed
	! from those coefficients in the parametrization, whose values 'exist'.

  use math_constants
  use utilities
  use modeldef
  implicit none

  ! * BOP
  ! equal to
  interface assignment (=)
     MODULE PROCEDURE CopyParamY
	 MODULE PROCEDURE FillParamY
  end interface

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE AddParamY_f
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE SubParamY_f
  END INTERFACE

  ! multiplication
  INTERFACE OPERATOR (*)
     MODULE PROCEDURE MultParamY_f
     MODULE PROCEDURE ScMultParamY_f
  END INTERFACE

  ! division
  INTERFACE OPERATOR (/)
     MODULE PROCEDURE ScDivParamY_f
  END INTERFACE
  ! * EOP

  public			:: CreateParamY, DeleteParamY, SetUpParamY, CopyParamY
  public			:: FillParamY, FillParamValuesY, VerifyParamY, VerifyLayersY
  public			:: CreateEmptyParamY, GetParamValuesY, PrintParamInfoY
  public			:: SetLayerY, SetCoeffValueY, GetCoeffValueY, GetCoeffY
	public      :: GetCoeffArrayY
  public			:: AddParamY_f, SubParamY_f, LinCombParamY
  public			:: MultParamY_f, DotProdParamY_f, DotProdVecParamY_f
  public			:: ScMultParamY_f, ScDivParamY_f
  public			:: SmoothVParamY, SmoothHParamY
  public			:: SmoothSqrtParamY, SmoothParamY

Contains

  ! **********************************************************************
  ! * BOP
  subroutine DeleteParamY(P)

    implicit none
    type (param_info), intent(inout)               :: P
    ! * EOP

    integer                                        :: status

    if(P%allocated) then
	  deallocate(P%F,STAT=status)
	  deallocate(P%L,STAT=status)
	  deallocate(P%c,STAT=status)
    end if

	P%nF = 0
	P%nL = 0
	P%nc = 0
	P%allocated = .FALSE.

  end subroutine DeleteParamY


  ! **********************************************************************
	! * Test whether two parametrizations have the same basic definitions
	! * Returns 1 if parametrizations are compatible, zero otherwise
  ! * BOP
  function VerifyParamY(P2,P1) result (status)

    implicit none
    type (param_info), intent(in)                  :: P1
    type (param_info), intent(in)                  :: P2
		logical                                       :: status
    ! * EOP

		integer                                        :: j

	status = .FALSE.

	if(.not.(P1%allocated.and.P2%allocated)) then
		write(0,*) 'Error: (VerifyParamY) parametrization not allocated yet'
		return
	end if

	! compare functional construction
	if (P1%nL /= P2%nL) then
		write(0,*) 'Error: (VerifyParamY) the two models have a different number of layers'
		return
	end if

	! compare layer construction
	if (P1%nF /= P2%nF) then
		write(0,*) 'Error: (VerifyParamY) the two models have a different number of functionals'
		return
	end if

	! compare parameter construction
	if (P1%nc /= P2%nc) then
		write(0,*) 'Error: (VerifyParamY) the two models have a different number of coefficients'
		return
	end if

	do j = 1,P1%nL
		status = VerifyLayersY(P1%L(j),P2%L(j))
		if (status .eqv. .FALSE.) then
			write(0,*) 'Error: (VerifyParamY) incompatible layer ',j
			return
		end if
	end do

  end function VerifyParamY


  ! **********************************************************************
	! * Test whether two layer structures have the same basic definitions
	! * Returns 1 if layers are identical, zero otherwise
  ! * BOP
  function VerifyLayersY(L1,L2) result (status)

    implicit none
    type (layer_info), intent(in)      :: L1,L2
    logical                            :: status
    ! * EOP

	status = .FALSE.

	if (L1%if_log .neqv. L2%if_log) then
		write(0,*) 'Error: (VerifyLayerY) layers have different structure'
		return
	end if

	if (clean(L1%depth)  /= clean(L2%depth)) then
		write(0,*) 'Error: (VerifyLayerY) layers have different depths'
		return
	end if

	if (clean(L1%width)  /= clean(L2%width)) then
		write(0,*) 'Error: (VerifyLayerY) layers have different widths'
		return
	end if

	if (clean(L1%lbound) /= clean(L2%lbound)) then
		write(0,*) 'Error: (VerifyLayerY) layers have different lower bounds'
		return
	end if

	if (clean(L1%ubound) /= clean(L2%ubound)) then
		write(0,*) 'Error: (VerifyLayerY) layers have different upper bounds'
		return
	end if

	if (clean(L1%alpha)  /= clean(L2%alpha)) then
		write(0,*) 'Error: (VerifyLayerY) layers have different horizontal regularization'
		return
	end if

	if (clean(L1%beta)   /= clean(L2%beta)) then
		write(0,*) 'Error: (VerifyLayerY) layers have different vertical regularization'
		return
	end if

	status = .TRUE.
	return

  end function VerifyLayersY


  ! **********************************************************************
  ! * BOP
  subroutine CreateParamY(P,nL,degree)

    implicit none
    type (param_info), intent(inout)               :: P
    integer, intent(in)                            :: nL,degree
    ! * EOP

    character(80)                                  :: code
    integer                                        :: status
    integer                                        :: l,m,i,j,nV

    if (P%allocated) then
       call DeleteParamY(P)
    end if

    ! create spherical harmonic functionals
    nV=0
    do l=0,degree
      nV = nV + (2*l+1)
    end do
	allocate(P%F(nV),STAT=status)

	P%nF = nV
	P%F(:)%name='Y'

	i=1
	P%F(i)%num=i
	P%F(i)%l=0
	P%F(i)%m=0
	do l=1,degree
	  i=i+1
	  P%F(i)%num=i
	  P%F(i)%l=l
	  P%F(i)%m=0
	  do m=1,l
		i=i+1
  		P%F(i)%num=i
		P%F(i)%l=l
		P%F(i)%m=m
		i=i+1
		P%F(i)%num=i
		P%F(i)%l=l
		P%F(i)%m=-m
	  end do
	end do


	! create layer construction
	allocate(P%L(nL),STAT=status)
	P%nL = nL
	do i=1,nL
	  P%L(i)%num = i
	  P%L(i)%lbound = R_ZERO
	  P%L(i)%ubound = R_ZERO
	  P%L(i)%width  = R_ZERO
	  P%L(i)%depth  = R_ZERO
	  P%L(i)%alpha  = 0.0d0
	  P%L(i)%beta   = 1.0d0
	  P%L(i)%if_log =.FALSE.
	end do


	! create parameter construction
	allocate(P%c(nL,nV),STAT=status)

	do j=1,nL
	  do i=1,nV
		P%c(j,i)%L => P%L(j)
		P%c(j,i)%F => P%F(i)
		write(code,'(2i5)') P%F(i)%l,P%F(i)%m
		P%c(j,i)%code   = trim(code)
		P%c(j,i)%value  = R_ZERO
		P%c(j,i)%min    = R_ZERO
		P%c(j,i)%max    = R_ZERO
		P%c(j,i)%frozen =.TRUE.
		P%c(j,i)%exists =.FALSE.
		! initialize pointers
		if (j==1) then
		  nullify(P%c(j,i)%prev)
		  P%c(j,i)%next => P%c(j+1,i)
		else if (j==P%nL) then
		  P%c(j,i)%prev => P%c(j-1,i)
		  nullify(P%c(j,i)%next)
		else
		  P%c(j,i)%prev => P%c(j-1,i)
		  P%c(j,i)%next => P%c(j+1,i)
		end if
	  end do
	end do

	P%nc = 0
	P%allocated = .TRUE.

  end subroutine CreateParamY


  ! **********************************************************************
  ! * Can be used to create a new zero-valued parametrization from another
  ! * parametrization that have been already initialized
  ! * BOP
  subroutine SetUpParamY(P,L,F)

    implicit none
    type (param_info), intent(inout)               :: P
	type (layer_info), dimension(:), intent(in)    :: L
	type (func_info), dimension(:), intent(in)	   :: F
    ! * EOP

	if (P%allocated) then
	  call DeleteParamY(P)
	end if

	call CreateParamY(P,size(L),maxval(F%l))

	if (P%nF /= size(F)) then
	   write(0,*) 'Error: (SetUpParamY) spherical harmonics initialization incorrect'
	   call DeleteParamY(P)
	   stop
	end if

	P%L = L
	P%F = F

  end subroutine SetUpParamY


  ! **********************************************************************
	! * Insert new coefficients into an existing parametrization structure.
	! * Only replaces
  ! * BOP
  subroutine FillParamY(P,c)

    implicit none
    type (param_info), intent(inout)               :: P
	type (coeff_info), dimension(:), intent(in)    :: c
    ! * EOP

	integer										   :: i,j,iL,iV

    if(.not.P%allocated) then
	   write(0,*) 'Error: (FillParamY) parametrization not allocated yet'
	   return
	end if

	do i=1,size(c)
	  iL = c(i)%L%num
	  iV = c(i)%F%num
	  if (P%c(iL,iV)%code /= c(i)%code) then
		write (0,*) 'Error: (FillParamY) parametrization setup error'
		stop
	  else if(.not.P%c(iL,iV)%exists) then
			write (0,*) 'Warning: (FillParamY) unable to replace coefficient - exists=.FALSE.'
		  cycle
		end if
	  P%c(iL,iV) = c(i)
	end do

	P%nc = count(P%c%exists)


  end subroutine FillParamY


  ! **********************************************************************
  ! * BOP
  subroutine FillParamValuesY(P,v)

    implicit none
    type (param_info), intent(inout)               :: P
    real(8), dimension(:), intent(in)			   :: v
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
	   write(0,*) 'Error: (FillParamValuesY) parametrization not allocated yet'
	   return
	end if

	if(count(P%c%exists) /= size(v)) then
	   write(0,*) 'Error: (FillParamValuesY) wrong number of input coefficients'
	   return
	end if

	k=0
	do j=1,P%nL
	  do i=1,P%nF
		if(.not.P%c(j,i)%exists) then
		  cycle
		end if
		k=k+1
		P%c(j,i)%value=v(k)
	  end do
	end do

  end subroutine FillParamValuesY


  ! **********************************************************************
  ! * BOP
  subroutine GetParamValuesY(P,v)

    implicit none
    type (param_info), intent(in)				   :: P
    real(8), dimension(:), intent(out)			   :: v
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
	   write(0,*) 'Error: (GetParamValuesY) parametrization not allocated yet'
	   return
	end if

	!do j=1,P%nL
	!  do i=1,P%nF
	!	if(P%c(j,i)%exists) then
	!	  print *, j,i,P%c(j,i)%value
	!	end if
	!  end do
	!end do


	if(count(P%c%exists) /= size(v)) then
	   write(0,*) 'Error: (GetParamValuesY) wrong number of input coefficients: ',size(v)
	   return
	end if

	k=0
	do j=1,P%nL
	  do i=1,P%nF
		if(.not.P%c(j,i)%exists) then
		  cycle
		end if
		k=k+1
		v(k)=P%c(j,i)%value
	  end do
	end do

  end subroutine GetParamValuesY


  ! **********************************************************************
  ! * BOP
  subroutine SetLayerY(P,iL,upperb,lowerb,alpha,beta,if_log)

    implicit none
    type (param_info), intent(inout)               :: P
	integer, intent(in)							   :: iL
    real(8), intent(in)							   :: upperb,lowerb
    real(8), intent(in)							   :: alpha,beta
	logical, intent(in)							   :: if_log
    ! * EOP

	real(8)										   :: depth, width
    integer                                        :: status

    if(.not.P%allocated) then
	   write(0,*) 'Error: (SetLayerY) parametrization not allocated yet'
	   return
	end if

    depth = EARTH_R - lowerb
	width  = upperb - lowerb
	if (width <= R_ZERO) then
	  write(0, *) 'Error: (SetLayerY) incorrect depth of layer ',iL
	  stop
	end if

	P%L(iL)%depth  = depth
	P%L(iL)%width  = width
	P%L(iL)%lbound = lowerb
	P%L(iL)%ubound = upperb
	P%L(iL)%alpha  = alpha
	P%L(iL)%beta   = beta
	P%L(iL)%if_log = if_log


  end subroutine SetLayerY


  ! **********************************************************************
  ! * BOP
  subroutine SetCoeffValueY(P,iL,l,m,v,min,max,frozen)

    implicit none
    type (param_info), intent(inout)               :: P
	integer, intent(in)							   :: iL,l,m
    real(8), intent(in)							   :: v,min,max
	logical, intent(in)							   :: frozen
    ! * EOP

	integer										   :: iV

    if(.not.P%allocated) then
	   write(0,*) 'Error: (SetCoeffValueY) parametrization not allocated yet'
	   return
	end if

	iV = iY(l,m)
	if (P%c(iL,iV)%exists) then
	  write (0,*) 'Warning: (SetCoeffValueY) coefficient ',iL,iV,' already exists, overwrite...'
	end if

	if ((P%c(iL,iV)%F%l /= l).or.(P%c(iL,iV)%F%m /= m)) then
	  write (0,*) 'Error: (SetCoeffValueY) major error in parametrization setup'
	  stop
	end if

	P%c(iL,iV)%value = v
	P%c(iL,iV)%min = min
	P%c(iL,iV)%max = max
	P%c(iL,iV)%frozen = frozen

	P%nc = P%nc + 1

	P%c(iL,iV)%exists = .TRUE.

  end subroutine SetCoeffValueY


  ! **********************************************************************
  ! * BOP
  subroutine GetCoeffValueY(P,iL,l,m,v,min,max,frozen)

    implicit none
    type (param_info), intent(in)				   :: P
	integer, intent(in)							   :: iL,l,m
    real(8), intent(out)						   :: v,min,max
	logical, intent(out)						   :: frozen
    ! * EOP

	integer										   :: iV

    if(.not.P%allocated) then
	   write(0,*) 'Error: (GetCoeffValueY) parametrization not allocated yet'
	   return
	end if

	iV = iY(l,m)
	if (.not.P%c(iL,iV)%exists) then
	  write (0,*) 'Error: (GetCoeffValueY) required coefficient does not exist'
	  return
	end if

	v = P%c(iL,iV)%value
	min = P%c(iL,iV)%min
	max = P%c(iL,iV)%max
	frozen = P%c(iL,iV)%frozen

  end subroutine GetCoeffValueY


  ! **********************************************************************
  ! * BOP
  function GetCoeffY(P,n) result (c)

    implicit none
    type (param_info), intent(in)				   :: P
	integer, intent(in)							   :: n
	type (coeff_info)							   :: c
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
	   write(0,*) 'Error: (GetCoeffY) parametrization not allocated yet'
	   return
	end if

	k=0
	search: do j=1,P%nL
	  do i=1,P%nF
		if(P%c(j,i)%exists) then
		  k=k+1
		  if(k==n) then
			exit search
		  end if
		end if
	  end do
	end do search

	c=P%c(j,i)

  end function GetCoeffY

  ! **********************************************************************
  ! * BOP
  subroutine GetCoeffArrayY(P,c)

    implicit none
    type (param_info), intent(in)				   :: P
	type (coeff_info), dimension(:), intent(out)   :: c
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
	   write(0,*) 'Error: (GetCoeffY) parametrization not allocated yet'
	   return
	end if

	if(size(c) /= P%nc) then
	   write(0,*) 'Error: (GetCoeffY) coefficient array is of a wrong size'
	   return
	end if


	k=0
	search: do j=1,P%nL
	  do i=1,P%nF
		if(P%c(j,i)%exists) then
		  k=k+1
			c(k)=P%c(j,i)
		end if
	  end do
	end do search

  end subroutine GetCoeffArrayY

  ! **********************************************************************
  ! * BOP
  subroutine CopyParamY(P2,P1)

    implicit none
    type (param_info), intent(in)                  :: P1
    type (param_info), intent(inout)               :: P2
    ! * EOP

    integer                                        :: status
	integer										   :: i,j

	if (P2%allocated) then
	  call DeleteParamY(P2)
	end if

	! copy functional construction
	P2%nF = P1%nF
	allocate(P2%F(P2%nF),STAT=status)
	P2%F = P1%F

	! copy layer construction
	P2%nL = P1%nL
	allocate(P2%L(P2%nL),STAT=status)
	P2%L = P1%L

	! copy parameter construction
	P2%nc = P1%nc
	allocate(P2%c(P2%nL,P2%nF),STAT=status)
	P2%c = P1%c

	do j=1,P2%nL
	  do i=1,P2%nF
		P2%c(j,i)%L => P2%L(j)
		P2%c(j,i)%F => P2%F(i)
	  end do
	end do

	P2%allocated = P1%allocated

  end subroutine CopyParamY


  ! **********************************************************************
  ! * Use this function when a zero-valued copy of an existing structure
  ! * is needed
  ! * BOP
  function CreateEmptyParamY(P1) result(P)

    implicit none
    type (param_info), intent(in)		:: P1
    type (param_info)					:: P
    ! * EOP

    if(.not.(P1%allocated)) then
	   write(0,*) 'Error: (CreateEmptyParamY) parametrization not allocated yet'
	   return
	end if

	P%allocated = .FALSE.

	! Create an identical structure to P1
	call CopyParamY(P,P1)

	! Set all values to zero
	P%c%value = R_ZERO

  end function CreateEmptyParamY


  ! **********************************************************************
  ! * We can add or subtract values with exists==.FALSE. (they equal zero)
  ! * BOP
  function AddParamY_f(P1,P2) result(P)

    implicit none
    type (param_info), intent(in)		:: P1
    type (param_info), intent(in)		:: P2
    type (param_info)					:: P
    ! * EOP

    integer					:: i,j

		if (.not.VerifyParamY(P1,P2)) then
	   write(0,*) 'Error: (AddParamY) parametrization structures incompatible'
	   return
		end if

	P = CreateEmptyParamY(P1)

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = P1%c(j,i)%value + P2%c(j,i)%value
	    else
	     write(0,*) 'Error: (AddParamY) parametrization not set up correctly'
	     stop
	    end if
	  end do
	end do

  end function AddParamY_f


  ! **********************************************************************
  ! * We can add or subtract values with exists==.FALSE. (they equal zero)
  ! * BOP
  function SubParamY_f(P1,P2) result(P)

    implicit none
    type (param_info), intent(in)		:: P1
    type (param_info), intent(in)		:: P2
    type (param_info)					:: P
    ! * EOP

    integer					:: i,j

		if (.not.VerifyParamY(P1,P2)) then
	   write(0,*) 'Error: (SubParamY) parametrization structures incompatible'
	   return
		end if

	P = CreateEmptyParamY(P1)

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = P1%c(j,i)%value - P2%c(j,i)%value
	    else
	     write(0,*) 'Error: (AddParamY) parametrization not set up correctly'
	     stop
	    end if
	  end do
	end do

  end function SubParamY_f


  ! **********************************************************************
  ! * We can multiply values with exists==.FALSE. (they equal zero)
  ! * BOP
  function MultParamY_f(P1,P2) result(P)

    implicit none
    type (param_info), intent(in)		:: P1
    type (param_info), intent(in)		:: P2
    type (param_info)					:: P
    ! * EOP

    integer					:: i,j

		if (.not.VerifyParamY(P1,P2)) then
	   write(0,*) 'Error: (MultParamY) parametrization structures incompatible'
	   return
		end if

    if(.not.P%allocated) then
	  P = CreateEmptyParamY(P1)
	end if

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = P1%c(j,i)%value * P2%c(j,i)%value
	    else
		 write(0,*) 'Error: (MultParamY) parametrization not set up correctly'
	     stop
	    end if
	  end do
	end do

  end function MultParamY_f

  ! **********************************************************************
  ! * We can multiply values with exists==.FALSE. (they equal zero)
  ! * BOP
  function DotProdParamY_f(P1,P2) result(r)

    implicit none
    type (param_info), intent(in)		:: P1
    type (param_info), intent(in)		:: P2
    real(8)								:: r
    ! * EOP

    integer					:: i,j

		if (.not.VerifyParamY(P1,P2)) then
	   write(0,*) 'Error: (DotProdParamY) parametrization structures incompatible'
	   return
	end if

    r = R_ZERO

	do j=1,P2%nL
	  do i=1,P2%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  if(.not.P1%c(j,i)%frozen) then
			r = r + P1%c(j,i)%value * P2%c(j,i)%value
		  end if
	    else
		 write(0,*) 'Error: (DotProdParamY) parametrization not set up correctly'
	     stop
	    end if
	  end do
	end do

  end function DotProdParamY_f


  ! **********************************************************************
  ! * If we multiply by a vector we only use values with exists==.TRUE.
  ! * BOP
  function DotProdVecParamY_f(v,P2) result(r)

    implicit none
    type (param_info), intent(in)		:: P2
    real(8), dimension(:), intent(in)	:: v
    real(8)								:: r
    ! * EOP

    integer					:: i,j,k

    if(.not.P2%allocated) then
	   write(0,*) 'Error: (DotProdVecParamY) parametrization not allocated yet'
	   return
	end if

	if(count(P2%c%exists) /= size(v)) then
	   write(0,*) 'Error: (DocProdVecParamY) wrong number of input coefficients'
	   return
	end if

    r = R_ZERO
	k = 0

	do j=1,P2%nL
	  do i=1,P2%nF
	    if(P2%c(j,i)%exists) then
		  if(P2%c(j,i)%frozen) then
		  else
			r = r + v(k) * P2%c(j,i)%value
		  end if
		  k = k+1
	    else
		  cycle
	    end if
	  end do
	end do

  end function DotProdVecParamY_f

  ! **********************************************************************
  ! * BOP
  function ScMultParamY_f(v,P1) result(P)

    implicit none
    type (param_info), intent(in)		:: P1
    real(8), intent(in)					:: v
    type (param_info)					:: P
    ! * EOP

    if(.not.P1%allocated) then
	   write(0,*) 'Error: (ScMultParamY) parametrization not allocated yet'
	   return
	end if

    if(.not.P%allocated) then
	  P = CreateEmptyParamY(P1)
	end if

	P%c%value = v * P1%c%value


  end function ScMultParamY_f

  ! **********************************************************************
  ! * BOP
  function ScDivParamY_f(P1,v) result(P)

    implicit none
    type (param_info), intent(in)		:: P1
    real(8), intent(in)					:: v
    type (param_info)					:: P
    ! * EOP

    if(.not.P1%allocated) then
	   write(0,*) 'Error: (ScDivParamY) parametrization not allocated yet'
	   return
	end if

    if(.not.P%allocated) then
	  P = CreateEmptyParamY(P1)
	end if

	if(v==R_ZERO) then
	   write(0,*) 'Error: (ScDivParamY) division by zero'
	   return
	end if

	P%c%value = P1%c%value / v


  end function ScDivParamY_f

  ! **********************************************************************
  ! * We can add or subtract values with exists==.FALSE. (they equal zero)
  ! * BOP
  subroutine LinCombParamY(r1,P1,r2,P2,P)

    implicit none
    type (param_info), intent(in)				   :: P1
    type (param_info), intent(in)				   :: P2
    real(8), intent(in)							   :: r1
    real(8), intent(in)							   :: r2
    type (param_info), intent(out)                 :: P
    ! * EOP

	integer										   :: i,j

		if (.not.VerifyParamY(P1,P2)) then
	   write(0,*) 'Error: (LinCombParamY) parametrization structures incompatible'
	   return
	end if

	P = CreateEmptyParamY(P1)

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		 P%c(j,i)%value = r1 * P1%c(j,i)%value + r2 * P2%c(j,i)%value
	    else
		 write(0,*) 'Error: (LinCombParamY) parametrization not set up correctly'
	     stop
	    end if
	  end do
	end do

  end subroutine LinCombParamY


  ! **********************************************************************
  ! * BOP
  subroutine SmoothHParamY(P)

    implicit none
    type (param_info), intent(inout)               :: P
    ! * EOP

	integer										   :: i,j

    if(.not.P%allocated) then
	   write(0,*) 'Error: (SmoothHParamY) parametrization not allocated yet'
	   return
	end if

	do j=1,P%nL
	  do i=1,P%nF
			P%c(j,i)%value = ((P%F(i)%l+1)**(-(P%L(j)%alpha)/2)) * P%c(j,i)%value
	  end do
	end do

  end subroutine SmoothHParamY


  ! **********************************************************************
  ! * Values with exists==.FALSE. have been initialized to zero and have
  ! * a smoothing effect: neighbouring parameters will have a smaller value
  ! * BOP
  subroutine SmoothVParamY(P)

    implicit none
    type (param_info), intent(inout)               :: P
    ! * EOP

	integer										   :: i,j
	real(8),dimension(P%nL,P%nF)				   :: v,v_prev,v_next

    if(.not.P%allocated) then
	   write(0,*) 'Error: (SmoothVParamY) parametrization not allocated yet'
	   return
	end if


	do j=1,P%nL
	  do i=1,P%nF

			if(associated(P%c(j,i)%prev).and.associated(P%c(j,i)%next)) then
			v_prev(j,i) = P%c(j,i)%prev%value * (1-P%L(j-1)%beta)/2
			v(j,i)      = P%c(j,i)%value * P%L(j)%beta
			v_next(j,i) = P%c(j,i)%next%value * (1-P%L(j)%beta)/2
			else if (associated(P%c(j,i)%next)) then
			v_prev(j,i) = R_ZERO
			v(j,i)      = P%c(j,i)%value * P%L(j)%beta
			v_next(j,i) = P%c(j,i)%next%value * (1-P%L(j)%beta)/2
			else if (associated(P%c(j,i)%prev)) then
			v_prev(j,i) = P%c(j,i)%prev%value * (1-P%L(j-1)%beta)/2
			v(j,i)      = P%c(j,i)%value * P%L(j)%beta
			v_next(j,i) = R_ZERO
			else
			v_prev(j,i) = R_ZERO
			v(j,i)      = P%c(j,i)%value
			v_next(j,i) = R_ZERO
			end if

	  end do
	end do


	P%c%value = v_prev + v + v_next


  end subroutine SmoothVParamY


  ! **********************************************************************
  ! * Preconditioning operator C_p^{1/2}
  ! * BOP
  function SmoothSqrtParamY(P1) result (P)

    implicit none
    type (param_info), intent(in)				   :: P1
    type (param_info)							   :: P
    ! * EOP

    if(.not.P1%allocated) then
	   write(0,*) 'Error: (SmoothSqrtParamY) parametrization not allocated yet'
	   return
	end if

	P = P1

	! make the smoothing operator symmetric
	call SmoothVParamY(P)
	call SmoothHParamY(P)
	call SmoothVParamY(P)

  end function SmoothSqrtParamY


  ! **********************************************************************
  ! * Preconditioning operator C_p
  ! * BOP
  function SmoothParamY(P1) result (P)

    implicit none
    type (param_info), intent(in)				   :: P1
    type (param_info)							   :: P2
    type (param_info)							   :: P
    ! * EOP

    if(.not.P1%allocated) then
	   write(0,*) 'Error: (SmoothParamY) parametrization not allocated yet'
	   return
	end if

	! apply operator C_p^{1/2} twice
	P2 = SmoothSqrtParamY(P1)
	P  = SmoothSqrtParamY(P2)


  end function SmoothParamY


  ! **********************************************************************
  ! * BOP
  subroutine PrintParamInfoY(P,verbose,comment)

    implicit none
    type (param_info), intent(in)               :: P
		integer, intent(in)                   :: verbose
		character(*),intent(in),optional			:: comment
    ! * EOP

	type (coeff_info)						:: coeff
	integer										   :: i,j

    if(.not.P%allocated) then
	   write(0,*) 'Error: (PrintParamInfoY) parametrization not allocated yet'
	   return
	end if

	write(0,*)
	if(present(comment)) then
		write(0,*) comment
	end if

	if (verbose>0) then
  write(0,'(a50,i3)') 'Number of layers in script: ',P%nL
	do j=1,P%nL
    write(0,'(a46,i2,a2,i3)') 'Number of coefficients in layer ',j,': ',count(.not.P%c(j,:)%frozen)
	end do
  write(0,'(a50,i3)') 'Number of variable parameters in script: ',count(.not.P%c%frozen)
	write(0,*)
	end if

	do j=1,P%nL
    write(0,'(a46,i2,a2,g15.7)') 'Degree and order zero coefficient in layer ',j,': ',P%c(j,1)%value
	end do
	write(0,*)

	if (verbose>0) then
	do j=1,P%nL
    write(0,'(a46,i2,a2,g10.5)') 'Horizontal regularization in layer ',j,': ',P%L(j)%alpha
	end do
	write(0,*)
	end if

	if (verbose>0) then
	do j=1,P%nL
    write(0,'(a46,i2,a2,g10.5)') 'Vertical regularization in layer ',j,': ',P%L(j)%beta
	end do
	write(0,*)
	end if

	if (verbose>3) then
	do i=1,P%nc
	  coeff=GetCoeffY(P,i)
	  write (0,'(2i8,g17.9)') coeff%L%num,coeff%F%num,coeff%value
	end do
	end if

  end subroutine PrintParamInfoY


end module model_operators
