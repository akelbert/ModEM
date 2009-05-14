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
     MODULE PROCEDURE copy_dataVecMTX
  end interface

  interface create
     MODULE PROCEDURE create_dataVecMTX
  end interface

  interface deall
     MODULE PROCEDURE deall_dataVecMTX
  end interface

  interface zero
     MODULE PROCEDURE zero_dataVecMTX
  end interface

  interface linComb
     MODULE PROCEDURE linComb_dataVecMTX
  end interface

  interface operator (*) ! interface to linComb
     MODULE PROCEDURE scMult_dataVecMTX_f
  end interface

  interface scMultAdd ! interface to linComb
     MODULE PROCEDURE scMultAdd_dataVecMTX
  end interface

  interface dotProd
     MODULE PROCEDURE dotProd_dataVecMTX_f
  end interface

  interface normalizeData
     MODULE PROCEDURE normalize_dataVecMTX
  end interface

  interface countData
     MODULE PROCEDURE count_dataVecMTX_f
  end interface

  interface countDataVec
     MODULE PROCEDURE countVec_dataVecMTX_f
  end interface


  ! basic operators for all dataVec types
  public			:: create_dataVecMTX
  public            :: deall_dataVecMTX
  public            :: zero_dataVecMTX
  public			:: copy_dataVecMTX
  public			:: linComb_dataVecMTX
  public            :: scMult_dataVecMTX_f
  public            :: scMultAdd_dataVecMTX
  public			:: dotProd_dataVecMTX_f
  public            :: normalize_dataVecMTX
  ! additional operators needed for the upper level code
  public			:: count_dataVecMTX_f, countVec_dataVecMTX_f

Contains

!-----------------------------!
! Basic operators: dataVecMTX !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataVecTX:
  ! a vector containing data for a single transmitter and all dataType's.
  ! note that this only allocates space for the dataVec's; after running
  ! this, need to run create_dataVec to fill it in. d%allocated is only
  ! set to .true. when all of the dataVec's are also allocated

  subroutine create_dataVecMTX(nTx, nDt, nRx, d)

    integer, intent(in)			        :: nTx, nDt, nRx
    type (dataVecMTX_t), intent(inout)  :: d
    ! local
    integer                             :: istat

    if(d%allocated) then
        call errStop('input data vector already allocated in create_dataVecMTX')
    endif

    d%ntx = nTx

    ! allocate space for the dataVec's
    allocate(d%v(nTx,nDt,nRx), STAT=istat)
    allocate(d%n(nTx,nDt), STAT=istat)

    d%allocated = .true.

  end subroutine create_dataVecMTX

  !**********************************************************************
  ! deallocates all memory and reinitializes dataVecMTX object d

  subroutine deall_dataVecMTX(d)

    type (dataVecMTX_t), intent(inout)	:: d
    ! local
    integer                             :: j,istat

    if(d%allocated) then
       !  deallocate everything relevant
       deallocate(d%v, d%n, STAT=istat)
    endif

    d%allocated = .false.

  end subroutine deall_dataVecMTX

  !**********************************************************************
  ! set the data values and error bars to zero

  function zero_dataVecMTX(d1) result (d2)

    type (dataVecMTX_t), intent(in)	    :: d1
    type (dataVecMTX_t)                 :: d2
    ! local
    integer                             :: j

	call copy_dataVecMTX(d2, d1)

	d2%v(:,:,:)%resp%value = C_ZERO

  end function zero_dataVecMTX

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVecMTX(d2, d1)

    type (dataVecMTX_t), intent(in)		:: d1
    type (dataVecMTX_t), intent(inout)	:: d2
    ! local variables
    integer								:: nTx, nDt, nRx
    integer                             :: i, j, k

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVecMTX')
    endif

    nTx = size(d1%v,1)
    nDt = size(d1%v,2)
    nRx = size(d1%v,3)

    ! check to see whether the LHS (d2) is allocated
    if (d2%allocated) then
    	call deall_dataVecMTX(d2)
    endif

    call create_dataVecMTX(nTx, nDt, nRx, d2)

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

  end subroutine copy_dataVecMTX

  ! **********************************************************************
  ! calculates linear combination of two dataVecMTX objects a*d1 + b*d2
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_dataVec and linComb_dataVec

  subroutine linComb_dataVecMTX(a,d1,b,d2,dOut)

    type (dataVecMTX_t), intent(in)			:: d1, d2
    real (kind=prec), intent(in)	        :: a, b
    type (dataVecMTX_t), intent(inout)		:: dOut
    ! local variables
    integer									:: nTx, nDt, nRx
    integer                                 :: i, j, k

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataVecMTX')
    endif

    ! check to see if inputs are compatable
    if (d1%nTx .ne. d2%nTx) then
       call errStop('input sizes not consistent in linComb_dataVecMTX')
    endif

    ! create the output vector that is consistent with inputs
    if(dOut%allocated) then
       call deall_dataVecMTX(dOut)
    end if

	dOut = zero_dataVecMTX(d1)

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

  end subroutine linComb_dataVecMTX

  ! **********************************************************************
  !  computes dOut = a * dIn for dataVecMTX objects dIn and dOut
  !  and real scalar a
  function scMult_dataVecMTX_f(a,dIn) result (dOut)

    real (kind=prec), intent(in)		:: a
    type (dataVecMTX_t), intent(in)	    :: dIn
    type (dataVecMTX_t)                 :: dOut

	call linComb_dataVecMTX(R_ZERO,dIn,a,dIn,dOut)

  end function scMult_dataVecMTX_f

  ! **********************************************************************
  !  computes d2+a*d1 for dataVecMTX objects d1, d2 and real scalar a,
  !   overwriting d2 (is this really needed?)
  subroutine scMultAdd_dataVecMTX(a, d1, d2)

    type (dataVecMTX_t), intent(in)			:: d1
    real (kind=prec), intent(in)	        :: a
    type (dataVecMTX_t), intent(inout)		:: d2

    call linComb_dataVecMTX(a,d1,ONE,d2,d2)

  end subroutine scMultAdd_dataVecMTX

  ! **********************************************************************
  ! dot product of two dataVecMTX objects r  = <d1,d2>
  ! note this is the *real* dot product NOT complex

  function dotProd_dataVecMTX_f(d1, d2) result(r)

    type (dataVecMTX_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: j
    real (kind=prec)					:: d1real, d1imag, d2real, d2imag

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVecMTX_f')
    endif

    ! check to see if inputs are compatable
    if (d1%nTx .ne. d2%nTx) then
       call errStop('input sizes not consistent in dotProd_dataVecMTX_f')
    endif

    nTx = size(d1%v,1)
    nDt = size(d1%v,2)
    nRx = size(d1%v,3)

    ! finally do the linear combination ...
    r = 0.0
    do i = 1, nTx
    	do j = 1, nDt
    		do k = 1, nRx
    			d1real = dreal(d1%v(i,j,k)%resp%value)
    			d2real = dreal(d2%v(i,j,k)%resp%value)
    			r = r + d1real * d2real
				d1imag = dimag(d1%v(i,j,k)%resp%value)
				d2imag = dimag(d2%v(i,j,k)%resp%value)
    			r = r + d1imag * d2imag
    		enddo
    	enddo
    enddo

  end function dotProd_dataVecMTX_f

  !**********************************************************************
  ! normalizes a dataVecMTX object using error bars from dataVec's
  ! (Add attribute "normalized" to all of the dataVec objects)
  ! --> if a copy of the original vector is needed, make a copy before
  !     calling this subroutine
  ! --> the optional argument N is 1 to divide by the errors,
  !     2 to divide by the error squared; default is of course 1

  subroutine normalize_dataVecMTX(d,N)

   type(dataVecMTX_t),intent(inout)          :: d
   integer, intent(in), optional             :: N

   !  local variables
   integer      			:: nn, i, j, k, nTx, nDt, nRx

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVecMTX')
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
        		if (d%v(i,j,k)%resp%err <= TOL6 * d%v(i,j,k)%resp%value) then
           			call errStop('data error bars too small in normalize_dataVec')
        		endif
        		d%v(i,j,k)%resp%value = d%v(i,j,k)%resp%value/(d%v(i,j,k)%resp%err**nn)
    		enddo
    	enddo
    enddo

  end subroutine normalize_dataVecMTX

!----------------------------------!
! Additional operators: dataVecMTX !
!----------------------------------!

  !**********************************************************************
  ! count all real data values in the full data vector dataVecMTX

  function count_dataVecMTX_f(d) result(ndata)

   type(dataVecMTX_t), intent(in)		:: d
   integer				:: ndata

   ndata = sum(d%n)

  end function count_dataVecMTX_f

  !**********************************************************************
  ! count the number of data vectors in the full data vector dataVecMTX.
  ! there should be one for each transmitter and data type

  function countVec_dataVecMTX_f(d) result(nvec)

   type(dataVecMTX_t), intent(in)		:: d
   integer				:: nvec

   nvec = count(d%n > 0)

  end function countVec_dataVecMTX_f


end module DataSpace
