!*****************************************************************************
module DataSpace
  ! Routines for the generic "data space". Modules in this group are
  ! indifferent to all details of the modeling code, including the
  ! grid.

  use utilities
  use math_constants
  implicit none

  interface assignment (=)
     MODULE PROCEDURE copy_dataVec
     MODULE PROCEDURE copy_dataVecTX
     MODULE PROCEDURE copy_dataVecMTX
  end interface

  interface create
     MODULE PROCEDURE create_dataVec
     MODULE PROCEDURE create_dataVecTX
     MODULE PROCEDURE create_dataVecMTX
  end interface

  interface deall
     MODULE PROCEDURE deall_dataVec
     MODULE PROCEDURE deall_dataVecTX
     MODULE PROCEDURE deall_dataVecMTX
  end interface

  interface zero
     MODULE PROCEDURE zero_dataVec
     MODULE PROCEDURE zero_dataVecTX
     MODULE PROCEDURE zero_dataVecMTX
  end interface

  interface linComb
     MODULE PROCEDURE linComb_dataVec
     MODULE PROCEDURE linComb_dataVecTX
     MODULE PROCEDURE linComb_dataVecMTX
  end interface

  interface operator (*) ! interface to linComb
     MODULE PROCEDURE scMult_dataVec_f
     MODULE PROCEDURE scMult_dataVecTX_f
     MODULE PROCEDURE scMult_dataVecMTX_f
  end interface

  interface scMultAdd ! interface to linComb
     MODULE PROCEDURE scMultAdd_dataVec
     MODULE PROCEDURE scMultAdd_dataVecTX
     MODULE PROCEDURE scMultAdd_dataVecMTX
  end interface

  interface dotProd
     MODULE PROCEDURE dotProd_dataVec_f
     MODULE PROCEDURE dotProd_dataVecTX_f
     MODULE PROCEDURE dotProd_dataVecMTX_f
  end interface

  interface normalizeData
     MODULE PROCEDURE normalize_dataVec
     MODULE PROCEDURE normalize_dataVecTX
     MODULE PROCEDURE normalize_dataVecMTX
  end interface

  interface countData
     MODULE PROCEDURE count_dataVecMTX_f
  end interface

  interface countDataVec
     MODULE PROCEDURE countVec_dataVecMTX_f
  end interface



  ! basic data vector containing data for a single transmitter & data type
  type :: dataVec_t

      ! nComp is number of EM components observed (e.g., 2 (3) for
      ! MT (with verticl field TFs)) times 2 to allow for real/imag;
      ! has to match the data type index dt
      integer 		:: nComp = 0

      ! nSite is number of sites where these components are observed
      integer		:: nSite = 0

      ! actual data; dimensions (nComp,nSite)
      real (kind=prec), pointer, dimension(:,:) :: value, error

      ! rx(nSite) is the array of indices into receiver dictionary
      integer, pointer, dimension(:) :: rx

      ! tx is index into transmitter dictionary
      integer		:: tx = 0

      ! dt is index into data type dictionary
      integer		:: dataType = 0

      logical       :: isComplex = .false.
      logical		:: errorBar = .false.
      logical       :: normalized = .false.
      logical		:: allocated = .false.

  end type dataVec_t


  ! collection of dataVec objects for all data types, single transmitter
  type :: dataVecTX_t
      ! the number of dataVecs (generally the number of data types) associated
      ! with this transmitter (note: not the total number of data types)
      integer       :: ndt = 0

      ! array of dataVec's, usually one for each data type (dimension ndt)
      type (dataVec_t), pointer, dimension(:)   :: data

      ! tx is the index into transmitter dictionary
      integer       :: tx = 0

      logical       :: allocated = .false.

  end type dataVecTX_t


  ! collection of dataVec objects for all transmitters
  type :: dataVecMTX_t
      ! ntx is number of transmitters, number of frequencies for MT
      integer		:: ntx = 0

      ! d is array of dataVec's for each transmitter (dimension nTX)
      type (dataVecTX_t), pointer, dimension(:)	:: d

      logical		:: allocated = .false.

  end type dataVecMTX_t

  ! basic operators for all dataVec types
  public			:: create_dataVec, create_dataVecTX, create_dataVecMTX
  public            :: deall_dataVec, deall_dataVecTX, deall_dataVecMTX
  public            :: zero_dataVec, zero_dataVecTX, zero_dataVecMTX
  public			:: copy_dataVec, copy_dataVecTX, copy_dataVecMTX
  public			:: linComb_dataVec, linComb_dataVecTX, linComb_dataVecMTX
  public            :: scMult_dataVec_f, scMult_dataVecTX_f, scMult_dataVecMTX_f
  public            :: scMultAdd_dataVec, scMultAdd_dataVecTX, scMultAdd_dataVecMTX
  public			:: dotProd_dataVec_f, dotProd_dataVecTX_f, dotProd_dataVecMTX_f
  public            :: normalize_dataVec, normalize_dataVecTX, normalize_dataVecMTX
  ! additional operators needed for the upper level code
  public			:: count_dataVecMTX_f, countVec_dataVecMTX_f

Contains

!-----------------------------!
! Basic operators: dataVec    !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataVec:
  ! a vector containing data for a single transmitter and dataType.
  ! --> nComp: number of components, nSite: number of sites
  ! --> Set d%errorBar = .true. to allocate storage for error bars also
  !     Note that keeping this attribute within dataVecMTX and dataVecTX
  !     in addition to the dataVec can easily lead to an internally
  !     inconsistent data vector. The same is true for d%normalized.
  ! --> all indices have to be set manually after the vector is created!

  subroutine create_dataVec(nComp, nSite, d, isComplex, errorBar)

    integer, intent(in)			        :: nComp, nSite
    type (dataVec_t), intent(inout)		:: d
    logical, intent(in), optional	    :: isComplex, errorBar
    ! local
    integer                             :: istat
    character(80)                       :: msg

    if(d%allocated) then

       call errStop('input data vector already allocated in create_dataVec')
    endif

    if (d%tx /= 0) then
       write(msg,*) 'please set the transmitter index ',d%tx,' after calling create_dataVec'
       call warning(msg)
    endif

    if (d%dataType /= 0) then
       write(msg,*) 'please set the data type index ',d%dataType,' after calling create_dataVec'
       call warning(msg)
    endif

    d%nComp = nComp
    d%nSite = nSite

    ! allocate and initialize the data values
    allocate(d%value(nComp, nSite), STAT=istat)
    d%value = R_ZERO

    d%tx = 0
    d%dataType = 0
    allocate(d%rx(nSite), STAT=istat)
    d%rx = 0

    ! set the isComplex parameter
    if (present(isComplex)) then
       d%isComplex = isComplex
    else
       d%isComplex = .false.
    endif

    if (d%isComplex .and. (mod(d%nComp,2) .ne. 0)) then
       call errStop('for complex data # of components must be even in create_dataVec')
    endif

    ! set the error bars
    if (present(errorBar)) then
       d%errorBar = errorBar
    else
       d%errorBar = .false.
    endif

    ! don't forget to check d%errorBar before accessing d%error in the code
    if (d%errorBar) then
       allocate(d%error(nComp, nSite), STAT=istat)
       d%error = R_ZERO
    endif

    d%normalized = .false.
    d%allocated = .true.

  end subroutine create_dataVec

  !**********************************************************************
  ! deallocates all memory and reinitialized dataVec object d
  ! NOTE: do not base your judgement on the deallocation of d%err
  ! on d%errorBar variable: it can be set to .false. somewhere in
  ! the program without deallocating d%err.

  subroutine deall_dataVec(d)

    type (dataVec_t), intent(inout)	:: d
    ! local
    integer                         :: istat

    if(d%allocated) then
       !  deallocate everything relevant
       deallocate(d%value, d%rx, STAT=istat)
       if (associated(d%error)) deallocate(d%error, STAT=istat)
    endif

    d%tx = 0
    d%dataType = 0
    d%isComplex = .false.
    d%errorBar = .false.
    d%normalized = .false.
    d%nComp = 0
    d%nSite = 0

    d%allocated = .false.

  end subroutine deall_dataVec

  !**********************************************************************
  ! set the data values and error bars to zero

  subroutine zero_dataVec(d)

    type (dataVec_t), intent(inout)	:: d

    if(d%allocated) then
       d%value = R_ZERO
       if (associated(d%error)) d%error = R_ZERO
       d%errorBar = .false.
    endif

  end subroutine zero_dataVec

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVec(d2, d1)

    type (dataVec_t), intent(in)		:: d1
    type (dataVec_t), intent(inout)		:: d2

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVec')
    endif

    ! check to see whether the LHS (d2) is allocated
    if (.not. d2%allocated) then
       call create_dataVec(d1%nComp, d1%nSite, d2, d1%isComplex, d1%errorBar)
    else
       if ((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite) .or. &
          (d1%isComplex .neqv. d2%isComplex) .or. &
          (d1%errorBar .neqv. d2%errorBar)) then
          ! deallocate d2, and reinitialize with correct parameters
          call deall_dataVec(d2)
          call create_dataVec(d1%nComp, d1%nSite, d2, d1%isComplex, d1%errorBar)
       endif
    endif

    ! now copy the components
    d2%value = d1%value
    if (d2%errorBar) then
       d2%error = d1%error
    endif
    d2%normalized = d1%normalized
    d2%rx = d1%rx
    d2%tx = d1%tx
    d2%dataType = d1%dataType

  end subroutine copy_dataVec

  ! **********************************************************************
  ! calculates linear combination of two dataVec objects a*d1 + b*d2
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_dataVec. Use the errors from one of the inputs;
  !   Error bars are not now defined when both input vectors have errors

    subroutine linComb_dataVec(a,d1,b,d2,dOut)

    type (dataVec_t), intent(in)			:: d1, d2
    real (kind=prec), intent(in)	        :: a, b
    type (dataVec_t), intent(inout)			:: dOut

    ! local variables
    logical					:: errBar = .false.

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataVec')
    endif

    ! check to see if inputs are compatable
    if ((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite) .or. &
    	(d1%isComplex .neqv. d2%isComplex)) then
       call errStop('input dataVecs not consistent in linComb_dataVec')
    endif

    ! check to see that the transmitter, dataType and receivers are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('input dataVecs correspond to different transmitters in linComb_dataVec')
    endif
    if (d1%dataType .ne. d2%dataType) then
       call errStop('input dataVecs correspond to different dataType in linComb_dataVec')
    endif
    if (maxval(abs(d1%rx - d2%rx)) > 0) then
       call errStop('input dataVecs correspond to different receiver sets in linComb_dataVec')
    endif

	! set errBar=.true. if at least one of the inputs has error bars
	errBar = (d1%errorBar .or. d2%errorBar)

    ! create the output vector that is consistent with inputs
    if(dOut%allocated) then
       call deall_dataVec(dOut)
    end if

    call create_dataVec(d1%nComp, d1%nSite, dOut, d1%isComplex, errBar)

	! set the receiver indices to those of d1
	dOut%tx = d1%tx
	dOut%dataType = d1%dataType
	dOut%rx = d1%rx

    !  finally do the linear combination ...
    dOut%value = a*d1%value + b*d2%value

	! deal with the error bars by copying them from one of the vectors;
	! currently exit if both of the input vectors have error bars defined
	! unless one of the multiplication coefficients is zero
	if (d1%errorBar .and. d2%errorBar) then
	   if ((abs(a) > R_ZERO).and.(abs(b) > R_ZERO)) then
       	  call errStop('unable to add two data vectors with error bars in linComb_dataVec')
       else
       	  dOut%error = a*d1%error + b*d2%error
       	  dOut%normalized = d1%normalized .and. d2%normalized
       end if
    else if(d1%errorBar) then
       dOut%error = a*d1%error
       dOut%normalized = d1%normalized
    else if(d2%errorBar) then
       dOut%error = b*d2%error
       dOut%normalized = d2%normalized
    end if

  end subroutine linComb_dataVec

  ! **********************************************************************
  !  computes dOut = a * dIn for dataVec objects dIn and dOut
  !  and real scalar a
  function scMult_dataVec_f(a,dIn) result (dOut)

    real (kind=prec), intent(in)		:: a
    type (dataVec_t), intent(in)	    :: dIn
    type (dataVec_t)                    :: dOut

	call linComb_dataVec(R_ZERO,dIn,a,dIn,dOut)

  end function scMult_dataVec_f

  ! **********************************************************************
  !  computes d2+a*d1 for dataVec objects d1, d2 and real scalar a,
  !   overwriting d2 (is this really needed?)
  subroutine scMultAdd_dataVec(a, d1, d2)

    type (dataVec_t), intent(in)			:: d1
    real (kind=prec), intent(in)	        :: a
    type (dataVec_t), intent(inout)		    :: d2

    call linComb_dataVec(a,d1,ONE,d2,d2)

  end subroutine scMultAdd_dataVec

  ! **********************************************************************
  ! dot product of two real dataVec objects r  = <d1,d2>

  function dotProd_dataVec_f(d1, d2) result(r)

    type (dataVec_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: j, k

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVec_f')
    endif

    ! check to see if inputs are compatable
    if ((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite) .or. &
    	(d1%isComplex .neqv. d2%isComplex)) then
       call errStop('input dataVecs not consistent in dotProd_dataVec_f')
    endif

    ! check to see that the transmitter, dataType and receivers are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('input dataVecs correspond to different transmitters in dotProd_dataVec_f')
    endif
    if (d1%dataType .ne. d2%dataType) then
       call errStop('input dataVecs correspond to different dataType in dotProd_dataVec_f')
    endif
    if (maxval(abs(d1%rx - d2%rx)) > 0) then
       call errStop('input dataVecs correspond to different receiver sets in dotProd_dataVec_f')
    endif

    r = 0.0
    do j = 1, d1%nComp
       do k = 1, d1%nSite
          r  =  r + d2%value(j,k) * d1%value(j,k)
       enddo
    enddo

  end function dotProd_dataVec_f

  !**********************************************************************
  ! normalizes a dataVec object using error bars
  ! (Add attribute "normalized" to the dataVec object)
  ! --> if a copy of the original vector is needed, make a copy before
  !     calling this subroutine
  ! --> the optional argument N is 1 to divide by the errors,
  !     2 to divide by the error squared; default is of course 1

  subroutine normalize_dataVec(d,N)

   type(dataVec_t),intent(inout)          :: d
   integer, intent(in), optional          :: N

   !  local variables
   integer      			:: nn, i, j

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVec')
   endif

   if (.not. d%errorBar) then
      call errStop('no error bars for input data in normalize_dataVec')
   endif

   if (d%normalized) then
      call errStop('data vector already normalized in normalize_dataVec')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

   do j = 1, d%nSite
     do i = 1, d%nComp
        if (d%error(i,j) <= TOL6 * d%value(i,j)) then
           call errStop('data error bars too small in normalize_dataVec')
        endif
        d%value(i,j) = d%value(i,j)/(d%error(i,j)**nn)
     enddo
   enddo

   d%normalized = .true.

  end subroutine normalize_dataVec

!-----------------------------!
! Basic operators: dataVecTX  !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataVecTX:
  ! a vector containing data for a single transmitter and all dataType's.
  ! note that this only allocates space for the dataVec's; after running
  ! this, need to run create_dataVec to fill it in. d%allocated is only
  ! set to .true. when all of the dataVec's are also allocated
  ! note that the transmitter index has to be set manually after calling
  ! this subroutine - this is done in order to be explicit about indices,
  ! so that a transmitter couldn't be specified before the vector exists.
  ! if it has been set before the vector is created, it is nullified here
  ! with a warning.

  subroutine create_dataVecTX(nDt, d)

    integer, intent(in)			        :: nDt
    type (dataVecTX_t), intent(inout)   :: d
    ! local
    integer                             :: istat
    character(80)                       :: msg

    if(d%allocated) then
        call errStop('input data vectors already allocated in create_dataVecTX')
    endif

    d%ndt = nDt

    ! allocate space for the dataVec's
    allocate(d%data(nDt), STAT=istat)

    if (d%tx /= 0) then
       write(msg,*) 'please set the transmitter index ',d%tx,' after calling create_dataVecTX'
       call warning(msg)
    endif

    d%tx = 0

    ! important - not allocated until all of the dataVec's are
    d%allocated = .false.

  end subroutine create_dataVecTX

  !**********************************************************************
  ! deallocates all memory and reinitializes dataVecTX object d

  subroutine deall_dataVecTX(d)

    type (dataVecTX_t), intent(inout)	:: d
    ! local
    integer                             :: i,istat

    if(d%allocated) then
       !  deallocate everything relevant
       do i = 1,d%nDt
          call deall_dataVec(d%data(i))
       enddo
       deallocate(d%data, STAT=istat)
    endif

    d%tx = 0
    d%allocated = .false.

  end subroutine deall_dataVecTX

  !**********************************************************************
  ! set the data values and error bars to zero

  subroutine zero_dataVecTX(d)

    type (dataVecTX_t), intent(inout)	:: d
    ! local
    integer                             :: i

    if(d%allocated) then
       do i = 1,d%nDt
          call zero_dataVec(d%data(i))
       enddo
    endif

  end subroutine zero_dataVecTX

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVecTX(d2, d1)

    type (dataVecTX_t), intent(in)		:: d1
    type (dataVecTX_t), intent(inout)		:: d2
    ! local variable
    integer                             :: i

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVecTX')
    endif

    ! check to see whether the LHS (d2) is allocated
    if (.not. d2%allocated) then
       call create_dataVecTX(d1%nDt, d2)
    else
       if (d1%nDt .ne. d2%nDt) then
          ! deallocate d2, and reinitialize with correct number of data vectors
          call deall_dataVecTX(d2)
          call create_dataVecTX(d1%nDt, d2)
       endif
    endif

    ! now copy the components
    do i = 1, d1%nDt
       d2%data(i) = d1%data(i)
    enddo
    d2%tx = d1%tx
    d2%allocated = .true.

  end subroutine copy_dataVecTX

  ! **********************************************************************
  ! calculates linear combination of two dataVecTX objects a*d1 + b*d2
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_dataVec and linComb_dataVec

  subroutine linComb_dataVecTX(a,d1,b,d2,dOut)

    type (dataVecTX_t), intent(in)			:: d1, d2
    real (kind=prec), intent(in)	        :: a, b
    type (dataVecTX_t), intent(inout)		:: dOut
    ! local
    integer                                 :: i

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataVecTX')
    endif

    ! check to see if inputs are compatable
    if (d1%nDt .ne. d2%nDt) then
       call errStop('input sizes not consistent in linComb_dataVecTX')
    endif

    ! check to see that the transmitters are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('inputs correspond to different transmitters in linComb_dataVecTX')
    endif

    ! create the output vector that is consistent with inputs
    if(dOut%allocated) then
       call deall_dataVecTX(dOut)
    end if

    call create_dataVecTX(d1%nDt, dOut)

	dOut%tx = d1%tx

    ! finally do the linear combination ...
    do i = 1, d1%nDt
       call linComb_dataVec(a,d1%data(i),b,d2%data(i),dOut%data(i))
    enddo

    dOut%allocated = d1%allocated .and. d2%allocated

  end subroutine linComb_dataVecTX

  ! **********************************************************************
  !  computes dOut = a * dIn for dataVecMTX objects dIn and dOut
  !  and real scalar a
  function scMult_dataVecTX_f(a,dIn) result (dOut)

    real (kind=prec), intent(in)		:: a
    type (dataVecTX_t), intent(in)	    :: dIn
    type (dataVecTX_t)                  :: dOut

	call linComb_dataVecTX(R_ZERO,dIn,a,dIn,dOut)

  end function scMult_dataVecTX_f

  ! **********************************************************************
  !  computes d2+a*d1 for dataVecTX objects d1, d2 and real scalar a,
  !   overwriting d2 (is this really needed?)
  subroutine scMultAdd_dataVecTX(a, d1, d2)

    type (dataVecTX_t), intent(in)			:: d1
    real (kind=prec), intent(in)	        :: a
    type (dataVecTX_t), intent(inout)		:: d2

    call linComb_dataVecTX(a,d1,ONE,d2,d2)

  end subroutine scMultAdd_dataVecTX

  ! **********************************************************************
  ! dot product of two real dataVecTX objects r  = <d1,d2>

  function dotProd_dataVecTX_f(d1, d2) result(r)

    type (dataVecTX_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: i

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVecTX_f')
    endif

    ! check to see if inputs are compatable
    if (d1%nDt .ne. d2%nDt) then
       call errStop('input sizes not consistent in dotProd_dataVecTX_f')
    endif

    ! check to see that the transmitters are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('inputs correspond to different transmitters in dotProd_dataVecTX_f')
    endif

    r = 0.0
    do i = 1, d1%nDt
       r  =  r + dotProd_dataVec_f(d1%data(i),d2%data(i))
    enddo

  end function dotProd_dataVecTX_f

  !**********************************************************************
  ! normalizes a dataVecTX object using error bars from dataVec's
  ! (Add attribute "normalized" to all of the dataVec objects)
  ! --> if a copy of the original vector is needed, make a copy before
  !     calling this subroutine
  ! --> the optional argument N is 1 to divide by the errors,
  !     2 to divide by the error squared; default is of course 1

  subroutine normalize_dataVecTX(d,N)

   type(dataVecTX_t),intent(inout)           :: d
   integer, intent(in), optional             :: N

   !  local variables
   integer      			:: nn, i

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVecMTX')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

   do i = 1, d%nDt
      call normalize_dataVec(d%data(i),nn)
   enddo

  end subroutine normalize_dataVecTX

!-----------------------------!
! Basic operators: dataVecMTX !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataVecTX:
  ! a vector containing data for a single transmitter and all dataType's.
  ! note that this only allocates space for the dataVec's; after running
  ! this, need to run create_dataVec to fill it in. d%allocated is only
  ! set to .true. when all of the dataVec's are also allocated

  subroutine create_dataVecMTX(nTx, d)

    integer, intent(in)			        :: nTx
    type (dataVecMTX_t), intent(inout)  :: d
    ! local
    integer                             :: istat

    if(d%allocated) then
        call errStop('input data vector already allocated in create_dataVecMTX')
    endif

    d%ntx = nTx

    ! allocate space for the dataVec's
    allocate(d%d(nTx), STAT=istat)

    ! not allocated until all of the dataVec's are
    d%allocated = .false.

  end subroutine create_dataVecMTX

  !**********************************************************************
  ! deallocates all memory and reinitializes dataVecMTX object d

  subroutine deall_dataVecMTX(d)

    type (dataVecMTX_t), intent(inout)	:: d
    ! local
    integer                             :: j,istat

    if(d%allocated) then
       !  deallocate everything relevant
       do j = 1,d%nTx
          call deall_dataVecTX(d%d(j))
       enddo
       deallocate(d%d, STAT=istat)
    endif

    d%allocated = .false.

  end subroutine deall_dataVecMTX

  !**********************************************************************
  ! set the data values and error bars to zero

  subroutine zero_dataVecMTX(d)

    type (dataVecMTX_t), intent(inout)	:: d
    ! local
    integer                             :: j

    if(d%allocated) then
       do j = 1,d%nTx
          call zero_dataVecTX(d%d(j))
       enddo
    endif

  end subroutine zero_dataVecMTX

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVecMTX(d2, d1)

    type (dataVecMTX_t), intent(in)		:: d1
    type (dataVecMTX_t), intent(inout)	:: d2
    ! local variable
    integer                             :: j

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVecMTX')
    endif

    ! check to see whether the LHS (d2) is allocated
    if (.not. d2%allocated) then
       call create_dataVecMTX(d1%nTx, d2)
    else
       if (d1%nTx .ne. d2%nTx) then
          ! deallocate d2, and reinitialize with correct number of data vectors
          call deall_dataVecMTX(d2)
          call create_dataVecMTX(d1%nTx, d2)
       endif
    endif

    ! now copy the components
    do j = 1, d1%nTx
       d2%d(j) = d1%d(j)
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
    ! local
    integer                                 :: j

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

    call create_dataVecMTX(d1%nTx, dOut)

    ! finally do the linear combination ...
    do j = 1, d1%nTx
       call linComb_dataVecTX(a,d1%d(j),b,d2%d(j),dOut%d(j))
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
  ! the real vs complex distinction is dealt with at the level of dataVec

  function dotProd_dataVecMTX_f(d1, d2) result(r)

    type (dataVecMTX_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: j

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVecMTX_f')
    endif

    ! check to see if inputs are compatable
    if (d1%nTx .ne. d2%nTx) then
       call errStop('input sizes not consistent in dotProd_dataVecMTX_f')
    endif

    r = 0.0
    do j = 1, d1%nTx
       r  =  r + dotProd_dataVecTX_f(d1%d(j),d2%d(j))
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
   integer      			:: nn, i, j

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVecMTX')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

   do j = 1, d%nTx
     do i = 1, d%d(j)%nDt
        call normalize_dataVec(d%d(j)%data(i),nn)
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

   ! local variables
   integer 				:: i,j
   ndata = 0
   do j = 1,d%nTx
     do i = 1,d%d(j)%nDt
       ndata = ndata + d%d(j)%data(i)%nComp * d%d(j)%data(i)%nSite
     enddo
   enddo

  end function count_dataVecMTX_f

  !**********************************************************************
  ! count the number of data vectors in the full data vector dataVecMTX.
  ! there should be one for each transmitter and data type

  function countVec_dataVecMTX_f(d) result(nvec)

   type(dataVecMTX_t), intent(in)		:: d
   integer				:: nvec

   ! local variables
   integer 				:: j
   nvec = 0
   do j = 1,d%nTx
     nvec = nvec + d%d(j)%nDt
   enddo

  end function countVec_dataVecMTX_f


end module DataSpace
