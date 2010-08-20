!*****************************************************************************
module DataSpace
  ! Routines for the generic "data space". Modules in this group are
  ! indifferent to all details of the modeling code, including the
  ! grid.

  use utilities
  use math_constants
  use datadef
  implicit none

  interface assignment (=)
     MODULE PROCEDURE copy_dataVectorMTX
  end interface

  interface create
     MODULE PROCEDURE create_dataVectorMTX
  end interface

  interface deall
     MODULE PROCEDURE deall_dataVectorMTX
  end interface

  interface compare
     MODULE PROCEDURE compare_dataVectorMTX_f
  end interface

  interface zero
     MODULE PROCEDURE zero_dataVectorMTX
  end interface

  interface linComb
     MODULE PROCEDURE linComb_dataVectorMTX
  end interface

  interface operator (*) ! interface to linComb
     MODULE PROCEDURE scMult_dataVectorMTX_f
  end interface

  interface scMultAdd ! interface to linComb
     MODULE PROCEDURE scMultAdd_dataVectorMTX
  end interface

  interface dotProd
     MODULE PROCEDURE dotProd_dataVectorMTX_f
  end interface

  interface normalizeData
     MODULE PROCEDURE normalize_dataVectorMTX
  end interface

  interface countData
     MODULE PROCEDURE count_dataVectorMTX_f
  end interface

  interface countDataVec
     MODULE PROCEDURE countVec_dataVectorMTX_f
  end interface


  ! basic operators for all dataVec types
  public			:: create_dataVectorMTX
  public            :: deall_dataVectorMTX
  public            :: compare_dataVectorMTX_f
  public            :: zero_dataVectorMTX
  public			:: copy_dataVectorMTX
  public			:: linComb_dataVectorMTX
  public            :: scMult_dataVectorMTX_f
  public            :: scMultAdd_dataVectorMTX
  public			:: dotProd_dataVectorMTX_f
  public            :: normalize_dataVectorMTX
  ! additional operators needed for the upper level code
  public			:: count_dataVectorMTX_f, countVec_dataVectorMTX_f

  ! ***************************************************************************
  ! * type response_t contains the full information about a single response at a
  ! * single frequency and observatory: value, error, anything else
  type :: response_t

	complex(8)								:: value
	real(8)									:: err
	logical									:: exists

  end type response_t

  ! ***************************************************************************
  ! * type dataValue_t contains the definition of a single data functionals used
  ! * in the expression for the penalty functional; they are specified for a
  ! * particular type of response and a frequency + observatory pair
  type :: dataValue_t

	type (transmitter_t), pointer			    :: freq
	type (receiver_t),	pointer			        :: obs
	type (functional_t),	pointer				:: func
	type (response_t)							:: resp
	! For the more general case, define dimension(freq%nMode,dataType%nComp)
	! type (response_t), pointer, dimension(:,:) :: resp

  end type dataValue_t


  ! ***************************************************************************
  ! * type dataVectorMTX_t contains the 3-D array of all data functionals we define;
  ! * it should be considered an ordered three-dimensional list of data
  ! * For some observatories data "do not exist" - either because they were
  ! * not physically present in the input data file, or since the observatory
  ! * is decided to be too close to one of the poles or the equator, so that
  ! * the responses at these locations would not be reliable. In these cases
  ! * we set v(i,j,k)%resp%exists to .FALSE.
  ! * The variable n is the number of "existing" data values for a given freq.
  ! * and functional type. It can also be computed by count(v(i,j,:)%resp%exists);
  ! * however, it is handy to have these values at hand.
  type :: dataVectorMTX_t

	type (dataValue_t), pointer, dimension(:,:,:) :: v	!nfreq,nfunc,nobs
	integer, pointer, dimension(:,:)			  :: n	!nfreq,nfunc
	! Total number of frequencies stored in this data vector
	integer                                       :: ntx
	logical                                       :: allocated = .false.
    ! needed to avoid memory leaks for temporary function outputs
    logical										  :: temporary = .false.

  end type dataVectorMTX_t

Contains

!-----------------------------!
! Basic operators: dataVectorMTX !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataVector:
  ! a vector containing data for a single transmitter and all dataType's.
  ! note that this only allocates space for the dataVec's; after running
  ! this, need to run create_dataVec to fill it in. d%allocated is only
  ! set to .true. when all of the dataVec's are also allocated

  subroutine create_dataVectorMTX(nTx, nDt, nRx, d)

    integer, intent(in)			        :: nTx, nDt, nRx
    type (dataVectorMTX_t), intent(inout)  :: d
    ! local
    integer                             :: istat

    if(d%allocated) then
        call errStop('input data vector already allocated in create_dataVectorMTX')
    endif

    d%ntx = nTx

    ! allocate space for the dataVec's
    allocate(d%v(nTx,nDt,nRx), STAT=istat)
    allocate(d%n(nTx,nDt), STAT=istat)

    d%allocated = .true.

  end subroutine create_dataVectorMTX

  !**********************************************************************
  ! deallocates all memory and reinitializes dataVectorMTX object d

  subroutine deall_dataVectorMTX(d)

    type (dataVectorMTX_t)	:: d
    ! local
    integer             :: j,istat

    if(d%allocated) then
       !  deallocate everything relevant
       deallocate(d%v, d%n, STAT=istat)
    endif

    d%allocated = .false.

  end subroutine deall_dataVectorMTX

  !**********************************************************************
  ! set the data values and error bars to zero; a function would make
  ! more sense, but d = zero(d) might create trouble, so make this a
  ! subroutine instead - much less error-prone

  subroutine zero_dataVectorMTX(d)

    type (dataVectorMTX_t), intent(inout)	:: d
    ! local
    integer                             :: j

	d%v(:,:,:)%resp%value = C_ZERO

  end subroutine zero_dataVectorMTX

  !**********************************************************************
  ! check whether the two inputs are of compatible size & nature
  ! for basic arithmetic operations

  function compare_dataVectorMTX_f(d1,d2) result (status)

    type (dataVectorMTX_t), intent(in)	    :: d1
    type (dataVectorMTX_t), intent(in)     :: d2
    logical								:: status
    ! local
    integer                             :: nTx1, nTx2
    integer                             :: nDt1, nDt2
    integer                             :: nRx1, nRx2

    status = .false.

    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
      return
    endif

    nTx1 = size(d1%v,1)
    nDt1 = size(d1%v,2)
    nRx1 = size(d1%v,3)

    nTx2 = size(d2%v,1)
    nDt2 = size(d2%v,2)
    nRx2 = size(d2%v,3)

    if ((nTx1 .eq. nTx2) .and. (nDt1 .eq. nDt2) .and. (nRx1 .eq. nRx2)) then
      status = .true.
    endif

  end function compare_dataVectorMTX_f

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVectorMTX(d2, d1)

    type (dataVectorMTX_t), intent(in)		:: d1
    type (dataVectorMTX_t), intent(inout)	:: d2
    ! local variables
    integer								:: nTx, nDt, nRx
    integer                             :: i, j, k

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVectorMTX')
    endif

    nTx = size(d1%v,1)
    nDt = size(d1%v,2)
    nRx = size(d1%v,3)

    ! check to see whether the LHS (d2) is allocated
    if (d2%allocated) then
    	call deall_dataVectorMTX(d2)
    endif

    call create_dataVectorMTX(nTx, nDt, nRx, d2)

    ! now copy the components
    d2%n = d1%n

    do i = 1, nTx
    	do j = 1, nDt
    		do k = 1, nRx
    			d2%v(i,j,k)%resp%value = d1%v(i,j,k)%resp%value
    			d2%v(i,j,k)%resp%err = d1%v(i,j,k)%resp%err
    			d2%v(i,j,k)%resp%exists = d1%v(i,j,k)%resp%exists
    			d2%v(i,j,k)%freq => d1%v(i,j,k)%freq
    			d2%v(i,j,k)%func => d1%v(i,j,k)%func
    			d2%v(i,j,k)%obs => d1%v(i,j,k)%obs
    		enddo
    	enddo
    enddo

    d2%allocated = .true.

    if (d1%temporary) then
    	call deall_dataVectorMTX(d1)
    endif

  end subroutine copy_dataVectorMTX

  ! **********************************************************************
  ! calculates linear combination of two dataVectorMTX objects a*d1 + b*d2
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_dataVec and linComb_dataVec
  !
  ! output can safely overwrite input, but output has to be allocated

  subroutine linComb_dataVectorMTX(a,d1,b,d2,dOut)

    type (dataVectorMTX_t), intent(in)			:: d1, d2
    real (kind=prec), intent(in)	        :: a, b
    type (dataVectorMTX_t), intent(inout)		:: dOut
    ! local variables
    integer									:: nTx, nDt, nRx
    integer                                 :: i, j, k

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataVectorMTX')
    endif

    ! check to see if inputs are of compatible sizes
    if (.not. compare(d1,d2)) then
       call errStop('input sizes not consistent in linComb_dataVectorMTX')
    endif

    ! create the output vector that is consistent with inputs
    if (.not. dOut%allocated) then
       call errStop('output structure has to be allocated before calling linComb_dataVectorMTX')
    end if

    ! check to see that output is of compatible size
    if (.not. compare(d1,dOut)) then
       call errStop('output size not consistent in linComb_dataVectorMTX')
    endif

    nTx = size(d1%v,1)
    nDt = size(d1%v,2)
    nRx = size(d1%v,3)

    ! finally do the linear combination ...
    do i = 1, nTx
    	do j = 1, nDt
    		do k = 1, nRx
    			dOut%v(i,j,k)%resp%value = &
    				a * d1%v(i,j,k)%resp%value + b * d2%v(i,j,k)%resp%value
    			dOut%v(i,j,k)%resp%err = &
    				a * d1%v(i,j,k)%resp%err + b * d2%v(i,j,k)%resp%err
    		enddo
    	enddo
    enddo

    dOut%allocated = d1%allocated .and. d2%allocated

  end subroutine linComb_dataVectorMTX

  ! **********************************************************************
  !  computes dOut = a * dIn for dataVectorMTX objects dIn and dOut
  !  and real scalar a
  function scMult_dataVectorMTX_f(a,dIn) result (dOut)

    real (kind=prec), intent(in)		:: a
    type (dataVectorMTX_t), intent(in)	    :: dIn
    type (dataVectorMTX_t)                 :: dOut

    dOut = dIn

	call linComb_dataVectorMTX(R_ZERO,dIn,a,dIn,dOut)

	dOut%temporary = .true.

  end function scMult_dataVectorMTX_f

  ! **********************************************************************
  !  computes d2+a*d1 for dataVectorMTX objects d1, d2 and real scalar a,
  !   overwriting d2 (is this really needed?)
  subroutine scMultAdd_dataVectorMTX(a, d1, d2)

    type (dataVectorMTX_t), intent(in)			:: d1
    real (kind=prec), intent(in)	        :: a
    type (dataVectorMTX_t), intent(inout)		:: d2

    call linComb_dataVectorMTX(a,d1,ONE,d2,d2)

  end subroutine scMultAdd_dataVectorMTX

  ! **********************************************************************
  ! dot product of two dataVectorMTX objects r  = <d1,d2>
  ! note this is the *real* dot product NOT complex

  function dotProd_dataVectorMTX_f(d1, d2) result(r)

    type (dataVectorMTX_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: nTx, nDt, nRx
    integer				                :: i, j, k
    real (kind=prec)					:: d1real, d1imag, d2real, d2imag

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVectorMTX_f')
    endif

    ! check to see if inputs are of compatible sizes
    if (.not. compare(d1,d2)) then
       call errStop('input sizes not consistent in dotProd_dataVectorMTX_f')
    endif

    nTx = size(d1%v,1)
    nDt = size(d1%v,2)
    nRx = size(d1%v,3)

    ! finally do the linear combination ...
    r = 0.0
    do i = 1, nTx
    	do j = 1, nDt
    		do k = 1, nRx
    			if (d1%v(i,j,k)%resp%exists .and. d2%v(i,j,k)%resp%exists) then
    				d1real = dreal(d1%v(i,j,k)%resp%value)
    				d2real = dreal(d2%v(i,j,k)%resp%value)
    				r = r + d1real * d2real
					d1imag = dimag(d1%v(i,j,k)%resp%value)
					d2imag = dimag(d2%v(i,j,k)%resp%value)
    				r = r + d1imag * d2imag
    			end if
    		enddo
    	enddo
    enddo

  end function dotProd_dataVectorMTX_f

  !**********************************************************************
  ! normalizes a dataVectorMTX object using error bars from dataVec's
  ! (Add attribute "normalized" to all of the dataVec objects)
  ! --> if a copy of the original vector is needed, make a copy before
  !     calling this subroutine
  ! --> the optional argument N is 1 to divide by the errors,
  !     2 to divide by the error squared; default is of course 1

  subroutine normalize_dataVectorMTX(d,N)

   type(dataVectorMTX_t),intent(inout)          :: d
   integer, intent(in), optional             :: N

   !  local variables
   integer      			:: nn, i, j, k, nTx, nDt, nRx

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVectorMTX')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

    nTx = size(d%v,1)
    nDt = size(d%v,2)
    nRx = size(d%v,3)

    ! finally do the linear combination ...
    do i = 1, nTx
    	do j = 1, nDt
    		do k = 1, nRx
        		if (abs(d%v(i,j,k)%resp%err) <= TOL6 * abs(d%v(i,j,k)%resp%value)) then
           			call errStop('data error bars too small in normalize_dataVec')
        		endif
        		d%v(i,j,k)%resp%value = d%v(i,j,k)%resp%value/(abs(d%v(i,j,k)%resp%err)**nn)
    		enddo
    	enddo
    enddo

  end subroutine normalize_dataVectorMTX

!----------------------------------!
! Additional operators: dataVectorMTX !
!----------------------------------!

  !**********************************************************************
  ! count all real data values in the full data vector dataVectorMTX

  function count_dataVectorMTX_f(d) result(ndata)

   type(dataVectorMTX_t), intent(in)		:: d
   integer				:: ndata

   ndata = 2 * sum(d%n)

  end function count_dataVectorMTX_f

  !**********************************************************************
  ! count the number of data vectors in the full data vector dataVectorMTX.
  ! there should be one for each transmitter and data type

  function countVec_dataVectorMTX_f(d) result(nvec)

   type(dataVectorMTX_t), intent(in)		:: d
   integer				:: nvec

   nvec = count(d%n > 0)

  end function countVec_dataVectorMTX_f


end module DataSpace
