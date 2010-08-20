!*****************************************************************************
module DataSpace
  ! Routines for the generic "data space". Modules in this group are
  ! indifferent to all details of the modeling code, including the
  ! grid.

  use utilities
  use math_constants
  implicit none

  interface assignment (=)
     MODULE PROCEDURE copy_dataBlock
     MODULE PROCEDURE copy_dataVector
     MODULE PROCEDURE copy_dataVectorMTX
  end interface

  interface create
     MODULE PROCEDURE create_dataBlock
     MODULE PROCEDURE create_dataVector
     MODULE PROCEDURE create_dataVectorMTX
  end interface

  interface deall
     MODULE PROCEDURE deall_dataBlock
     MODULE PROCEDURE deall_dataVector
     MODULE PROCEDURE deall_dataVectorMTX
  end interface

  interface zero
     MODULE PROCEDURE zero_dataBlock
     MODULE PROCEDURE zero_dataVector
     MODULE PROCEDURE zero_dataVectorMTX
  end interface

  interface linComb
     MODULE PROCEDURE linComb_dataBlock
     MODULE PROCEDURE linComb_dataVector
     MODULE PROCEDURE linComb_dataVectorMTX
  end interface

  interface scMult ! operator (*) interface to linComb
     MODULE PROCEDURE scMult_dataBlock
     MODULE PROCEDURE scMult_dataVector
     MODULE PROCEDURE scMult_dataVectorMTX
  end interface

  interface scMultAdd ! interface to linComb
     MODULE PROCEDURE scMultAdd_dataBlock
     MODULE PROCEDURE scMultAdd_dataVector
     MODULE PROCEDURE scMultAdd_dataVectorMTX
  end interface

  interface dotProd
     MODULE PROCEDURE dotProd_dataBlock_f
     MODULE PROCEDURE dotProd_dataVector_f
     MODULE PROCEDURE dotProd_dataVectorMTX_f
  end interface

  interface normalizeData
     MODULE PROCEDURE normalize_dataBlock
     MODULE PROCEDURE normalize_dataVector
     MODULE PROCEDURE normalize_dataVectorMTX
  end interface

  interface countData
     MODULE PROCEDURE count_dataVectorMTX_f
  end interface

  interface countDataVec
     MODULE PROCEDURE countVec_dataVectorMTX_f
  end interface



  ! basic data vector containing data for a single transmitter & data type
  type :: dataBlock_t

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

      ! needed to avoid memory leaks for temporary function outputs
      logical		:: temporary = .false.

  end type dataBlock_t


  ! collection of dataVec objects for all data types, single transmitter
  type :: dataVector_t
      ! the number of dataVecs (generally the number of data types) associated
      ! with this transmitter (note: not the total number of data types)
      integer       :: ndt = 0

      ! array of dataVec's, usually one for each data type (dimension ndt)
      type (dataBlock_t), pointer, dimension(:)   :: data

      ! tx is the index into transmitter dictionary
      integer       :: tx = 0

      logical       :: allocated = .false.

      ! needed to avoid memory leaks for temporary function outputs
      logical		:: temporary = .false.

  end type dataVector_t


  ! collection of dataVec objects for all transmitters
  type :: dataVectorMTX_t
      ! ntx is number of transmitters, number of frequencies for MT
      integer		:: ntx = 0

      ! d is array of dataVec's for each transmitter (dimension nTX)
      type (dataVector_t), pointer, dimension(:)	:: d

      logical		:: allocated = .false.

      ! needed to avoid memory leaks for temporary function outputs
      logical		:: temporary = .false.

  end type dataVectorMTX_t

  ! basic operators for all dataVec types
  public			:: create_dataBlock, create_dataVector, create_dataVectorMTX
  public            :: deall_dataBlock, deall_dataVector, deall_dataVectorMTX
  public            :: zero_dataBlock, zero_dataVector, zero_dataVectorMTX
  public			:: copy_dataBlock, copy_dataVector, copy_dataVectorMTX
  public			:: linComb_dataBlock, linComb_dataVector, linComb_dataVectorMTX
  public            :: scMult_dataBlock, scMult_dataVector, scMult_dataVectorMTX
  public            :: scMultAdd_dataBlock, scMultAdd_dataVector, scMultAdd_dataVectorMTX
  public			:: dotProd_dataBlock_f, dotProd_dataVector_f, dotProd_dataVectorMTX_f
  public            :: normalize_dataBlock, normalize_dataVector, normalize_dataVectorMTX
  ! additional operators needed for the upper level code
  public			:: count_dataVectorMTX_f, countVec_dataVectorMTX_f

Contains

!-----------------------------!
! Basic operators: dataBlock  !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataBlock:
  ! a vector containing data for a single transmitter and dataType.
  ! --> nComp: number of components, nSite: number of sites
  ! --> Set d%errorBar = .true. to allocate storage for error bars also
  !     Note that keeping this attribute within dataVectorMTX and dataVector
  !     in addition to the dataVec can easily lead to an internally
  !     inconsistent data vector. The same is true for d%normalized.
  ! --> all indices have to be set manually after the vector is created!

  subroutine create_dataBlock(nComp, nSite, d, isComplex, errorBar)

    integer, intent(in)			        :: nComp, nSite
    type (dataBlock_t), intent(inout)		:: d
    logical, intent(in), optional	    :: isComplex, errorBar
    ! local
    integer                             :: istat
    character(80)                       :: msg

    if(d%allocated) then

       call errStop('input data vector already allocated in create_dataBlock')
    endif

    if (d%tx /= 0) then
       write(msg,*) 'please set the transmitter index ',d%tx,' after calling create_dataBlock'
       call warning(msg)
    endif

    if (d%dataType /= 0) then
       write(msg,*) 'please set the data type index ',d%dataType,' after calling create_dataBlock'
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
       call errStop('for complex data # of components must be even in create_dataBlock')
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

  end subroutine create_dataBlock

  !**********************************************************************
  ! deallocates all memory and reinitialized dataVec object d
  ! NOTE: do not base your judgement on the deallocation of d%err
  ! on d%errorBar variable: it can be set to .false. somewhere in
  ! the program without deallocating d%err.

  subroutine deall_dataBlock(d)

    type (dataBlock_t)	:: d
    ! local
    integer             :: istat

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

  end subroutine deall_dataBlock

  !**********************************************************************
  ! set the data values and error bars to zero; function would make more
  ! sense here, but d = zero(d) might create trouble, so make this a
  ! subroutine instead

  subroutine zero_dataBlock(d)

    type (dataBlock_t), intent(inout)	:: d

    if(d%allocated) then
       d%value = R_ZERO
       if (associated(d%error)) d%error = R_ZERO
       d%errorBar = .false.
    endif

  end subroutine zero_dataBlock

  ! **********************************************************************
  ! copy a data block from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataBlock(d2, d1)

    type (dataBlock_t), intent(in)		:: d1
    type (dataBlock_t), intent(inout)		:: d2

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVec')
    endif

    ! check to see whether the LHS (d2) is allocated
    if (.not. d2%allocated) then
       call create_dataBlock(d1%nComp, d1%nSite, d2, d1%isComplex, d1%errorBar)
    else
       if ((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite) .or. &
          (d1%isComplex .neqv. d2%isComplex) .or. &
          (d1%errorBar .neqv. d2%errorBar)) then
          ! deallocate d2, and reinitialize with correct parameters
          call deall_dataBlock(d2)
          call create_dataBlock(d1%nComp, d1%nSite, d2, d1%isComplex, d1%errorBar)
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

    ! if input is a temporary function output, deallocate
    if (d1%temporary) then
    	call deall_dataBlock(d1)
    endif

  end subroutine copy_dataBlock

  ! **********************************************************************
  ! calculates linear combination of two dataBlock objects a*d1 + b*d2
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_dataBlock. Use the errors from one of the inputs;
  !   Error bars are not now defined when both input vectors have errors

    subroutine linComb_dataBlock(a,d1,b,d2,dOut)

    type (dataBlock_t), intent(in)			:: d1, d2
    real (kind=prec), intent(in)	        :: a, b
    type (dataBlock_t), intent(inout)			:: dOut

    ! local variables
    logical					:: errBar = .false.
    integer					:: istat

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataBlock')
    endif

    ! check to see if inputs are compatable
    if ((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite) .or. &
    	(d1%isComplex .neqv. d2%isComplex)) then
       call errStop('input dataVecs not consistent in linComb_dataBlock')
    endif

    ! check to see that the transmitter, dataType and receivers are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('input dataVecs correspond to different transmitters in linComb_dataBlock')
    endif
    if (d1%dataType .ne. d2%dataType) then
       call errStop('input dataVecs correspond to different dataType in linComb_dataBlock')
    endif
    if (maxval(abs(d1%rx - d2%rx)) > 0) then
       call errStop('input dataVecs correspond to different receiver sets in linComb_dataBlock')
    endif

    ! create the output vector that is consistent with inputs
    if (.not. dOut%allocated) then
       call errStop('output structure has to be allocated before calling linComb_dataBlock')
    end if

    ! check to see if inputs and output are compatible
    if ((d1%nComp .ne. dOut%nComp) .or. (d1%nSite .ne. dOut%nSite) .or. &
    	(d1%isComplex .neqv. dOut%isComplex)) then
       call errStop('input and output dataVecs not consistent in linComb_dataBlock')
    endif

	! set the receiver indices to those of d1
	dOut%tx = d1%tx
	dOut%dataType = d1%dataType
	dOut%rx = d1%rx

    !  finally do the linear combination ...
    dOut%value = a*d1%value + b*d2%value

	! set errBar=.true. if at least one of the inputs has error bars
	errBar = (d1%errorBar .or. d2%errorBar)
	dOut%errorBar = errBar

	! allocate error bars, if needed
    if (errBar .and. .not. associated(dOut%error)) then
       allocate(dOut%error(dOut%nComp, dOut%nSite), STAT=istat)
    endif

	! deal with the error bars by copying them from one of the vectors;
	! currently exit if both of the input vectors have error bars defined
	! unless one of the multiplication coefficients is zero
	if (d1%errorBar .and. d2%errorBar) then
	   if ((abs(a) > R_ZERO).and.(abs(b) > R_ZERO)) then
       	  call errStop('unable to add two data vectors with error bars in linComb_dataBlock')
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

  end subroutine linComb_dataBlock

  ! **********************************************************************
  !  computes dOut = a * dIn for dataBlock objects dIn and dOut
  !  and real scalar a
  subroutine scMult_dataBlock(a,d1,d2)

    real (kind=prec), intent(in)		:: a
    type (dataBlock_t), intent(in)	    :: d1
    type (dataBlock_t), intent(inout)     :: d2

	call linComb_dataBlock(R_ZERO,d1,a,d1,d2)

  end subroutine scMult_dataBlock

  ! **********************************************************************
  !  computes d2+a*d1 for dataBlock objects d1, d2 and real scalar a,
  !   overwriting d2 (is this really needed?)
  subroutine scMultAdd_dataBlock(a, d1, d2)

    type (dataBlock_t), intent(in)			:: d1
    real (kind=prec), intent(in)	        :: a
    type (dataBlock_t), intent(inout)		    :: d2

    call linComb_dataBlock(a,d1,ONE,d2,d2)

  end subroutine scMultAdd_dataBlock

  ! **********************************************************************
  ! dot product of two real dataBlock objects r  = <d1,d2>

  function dotProd_dataBlock_f(d1, d2) result(r)

    type (dataBlock_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: j, k

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataBlock_f')
    endif

    ! check to see if inputs are compatable
    if ((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite) .or. &
    	(d1%isComplex .neqv. d2%isComplex)) then
       call errStop('input dataVecs not consistent in dotProd_dataBlock_f')
    endif

    ! check to see that the transmitter, dataType and receivers are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('input dataVecs correspond to different transmitters in dotProd_dataBlock_f')
    endif
    if (d1%dataType .ne. d2%dataType) then
       call errStop('input dataVecs correspond to different dataType in dotProd_dataBlock_f')
    endif
    if (maxval(abs(d1%rx - d2%rx)) > 0) then
       call errStop('input dataVecs correspond to different receiver sets in dotProd_dataBlock_f')
    endif

    r = 0.0
    do j = 1, d1%nComp
       do k = 1, d1%nSite
          r  =  r + d2%value(j,k) * d1%value(j,k)
       enddo
    enddo

  end function dotProd_dataBlock_f

  !**********************************************************************
  ! normalizes a dataBlock object using error bars
  ! (Add attribute "normalized" to the dataBlock object)
  ! --> if a copy of the original vector is needed, make a copy before
  !     calling this subroutine
  ! --> the optional argument N is 1 to divide by the errors,
  !     2 to divide by the error squared; default is of course 1

  subroutine normalize_dataBlock(d,N)

   type(dataBlock_t),intent(inout)          :: d
   integer, intent(in), optional          :: N

   !  local variables
   integer      			:: nn, i, j

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataBlock')
   endif

   if (.not. d%errorBar) then
      call errStop('no error bars for input data in normalize_dataBlock')
   endif

   if (d%normalized) then
      call errStop('data vector already normalized in normalize_dataBlock')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

   do j = 1, d%nSite
     do i = 1, d%nComp
        if (d%error(i,j) <= TOL6 * d%value(i,j)) then
           call errStop('data error bars too small in normalize_dataBlock')
        endif
        d%value(i,j) = d%value(i,j)/(d%error(i,j)**nn)
     enddo
   enddo

   d%normalized = .true.

  end subroutine normalize_dataBlock

!-----------------------------!
! Basic operators: dataVector !
!-----------------------------!

  ! **********************************************************************
  ! creates an object of type dataVector:
  ! a vector containing data for a single transmitter and all dataType's.
  ! note that this only allocates space for the dataVec's; after running
  ! this, need to run create_dataVec to fill it in. d%allocated is only
  ! set to .true. when all of the dataVec's are also allocated
  ! note that the transmitter index has to be set manually after calling
  ! this subroutine - this is done in order to be explicit about indices,
  ! so that a transmitter couldn't be specified before the vector exists.
  ! if it has been set before the vector is created, it is nullified here
  ! with a warning.

  subroutine create_dataVector(nDt, d)

    integer, intent(in)			        :: nDt
    type (dataVector_t), intent(inout)   :: d
    ! local
    integer                             :: istat
    character(80)                       :: msg

    if(d%allocated) then
        call errStop('input data vectors already allocated in create_dataVector')
    endif

    d%ndt = nDt

    ! allocate space for the dataVec's
    allocate(d%data(nDt), STAT=istat)

    if (d%tx /= 0) then
       write(msg,*) 'please set the transmitter index ',d%tx,' after calling create_dataVector'
       call warning(msg)
    endif

    d%tx = 0

    ! important - not allocated until all of the dataVec's are
    d%allocated = .false.

  end subroutine create_dataVector

  !**********************************************************************
  ! deallocates all memory and reinitializes dataVector object d

  subroutine deall_dataVector(d)

    type (dataVector_t)	:: d
    ! local
    integer             :: i,istat

    if(d%allocated) then
       !  deallocate everything relevant
       do i = 1,d%nDt
          call deall_dataBlock(d%data(i))
       enddo
       deallocate(d%data, STAT=istat)
    endif

    d%tx = 0
    d%allocated = .false.

  end subroutine deall_dataVector

  !**********************************************************************
  ! set the data values and error bars to zero; function would make more
  ! sense here, but d = zero(d) might create trouble, so make this a
  ! subroutine instead

  subroutine zero_dataVector(d)

    type (dataVector_t), intent(inout)	:: d
    ! local
    integer                             :: i

    if(d%allocated) then
       do i = 1,d%nDt
          call zero_dataBlock(d%data(i))
       enddo
    else
       call errStop('data vector not allocated in zero_dataVector')
    endif

  end subroutine zero_dataVector

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVector(d2, d1)

    type (dataVector_t), intent(in)		:: d1
    type (dataVector_t), intent(inout)	:: d2
    ! local variable
    integer                             :: i

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVector')
    endif

    ! check to see whether the LHS (d2) is allocated
    if (.not. d2%allocated) then
       call create_dataVector(d1%nDt, d2)
    else
       if (d1%nDt .ne. d2%nDt) then
          ! deallocate d2, and reinitialize with correct number of data vectors
          call deall_dataVector(d2)
          call create_dataVector(d1%nDt, d2)
       endif
    endif

    ! now copy the components
    do i = 1, d1%nDt
       d2%data(i) = d1%data(i)
    enddo
    d2%tx = d1%tx
    d2%allocated = .true.

    ! if input is a temporary function output, deallocate
    if (d1%temporary) then
    	call deall_dataVector(d1)
    endif

  end subroutine copy_dataVector

  ! **********************************************************************
  ! calculates linear combination of two dataVector objects a*d1 + b*d2
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_dataVec and linComb_dataVec

  subroutine linComb_dataVector(a,d1,b,d2,dOut)

    type (dataVector_t), intent(in)			:: d1, d2
    real (kind=prec), intent(in)	        :: a, b
    type (dataVector_t), intent(inout)		:: dOut
    ! local
    integer                                 :: i

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataVector')
    endif

    ! check to see if inputs are of compatible sizes
    if (d1%nDt .ne. d2%nDt) then
       call errStop('input sizes not consistent in linComb_dataVector')
    endif

    ! check to see that the transmitters are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('inputs correspond to different transmitters in linComb_dataVector')
    endif

    ! create the output vector that is consistent with inputs
    if (.not. dOut%allocated) then
       call errStop('output structure has to be allocated before calling linComb_dataVector')
    end if

    ! check to see that output is of compatible size
    if (d1%nDt .ne. dOut%nDt) then
       call errStop('output size not consistent in linComb_dataVector')
    endif

    ! check to see that the transmitters are the same
    if (d1%tx .ne. dOut%tx) then
       call errStop('input and output correspond to different transmitters in linComb_dataVector')
    endif

    ! finally do the linear combination ...
    do i = 1, d1%nDt
       call linComb_dataBlock(a,d1%data(i),b,d2%data(i),dOut%data(i))
    enddo

    dOut%allocated = d1%allocated .and. d2%allocated

  end subroutine linComb_dataVector

  ! **********************************************************************
  !  computes dOut = a * dIn for dataVectorMTX objects dIn and dOut
  !  and real scalar a
  subroutine scMult_dataVector(a,d1,d2)

    real (kind=prec), intent(in)		:: a
    type (dataVector_t), intent(in)	    :: d1
    type (dataVector_t), intent(inout)   :: d2

	call linComb_dataVector(R_ZERO,d1,a,d1,d2)

  end subroutine scMult_dataVector

  ! **********************************************************************
  !  computes d2+a*d1 for dataVector objects d1, d2 and real scalar a,
  !   overwriting d2 (is this really needed?)
  subroutine scMultAdd_dataVector(a, d1, d2)

    type (dataVector_t), intent(in)			:: d1
    real (kind=prec), intent(in)	        :: a
    type (dataVector_t), intent(inout)		:: d2

    call linComb_dataVector(a,d1,ONE,d2,d2)

  end subroutine scMultAdd_dataVector

  ! **********************************************************************
  ! dot product of two real dataVector objects r  = <d1,d2>

  function dotProd_dataVector_f(d1, d2) result(r)

    type (dataVector_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: i

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVector_f')
    endif

    ! check to see if inputs are compatable
    if (d1%nDt .ne. d2%nDt) then
       call errStop('input sizes not consistent in dotProd_dataVector_f')
    endif

    ! check to see that the transmitters are the same
    if (d1%tx .ne. d2%tx) then
       call errStop('inputs correspond to different transmitters in dotProd_dataVector_f')
    endif

    r = 0.0
    do i = 1, d1%nDt
       r  =  r + dotProd_dataBlock_f(d1%data(i),d2%data(i))
    enddo

  end function dotProd_dataVector_f

  !**********************************************************************
  ! normalizes a dataVector object using error bars from dataVec's
  ! (Add attribute "normalized" to all of the dataVec objects)
  ! --> if a copy of the original vector is needed, make a copy before
  !     calling this subroutine
  ! --> the optional argument N is 1 to divide by the errors,
  !     2 to divide by the error squared; default is of course 1

  subroutine normalize_dataVector(d,N)

   type(dataVector_t),intent(inout)           :: d
   integer, intent(in), optional             :: N

   !  local variables
   integer      			:: nn, i

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVectorMTX')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

   do i = 1, d%nDt
      call normalize_dataBlock(d%data(i),nn)
   enddo

  end subroutine normalize_dataVector

!--------------------------------!
! Basic operators: dataVectorMTX !
!--------------------------------!

  ! **********************************************************************
  ! creates an object of type dataVector:
  ! a vector containing data for a single transmitter and all dataType's.
  ! note that this only allocates space for the dataVec's; after running
  ! this, need to run create_dataVec to fill it in. d%allocated is only
  ! set to .true. when all of the dataVec's are also allocated

  subroutine create_dataVectorMTX(nTx, d)

    integer, intent(in)			        :: nTx
    type (dataVectorMTX_t), intent(inout)  :: d
    ! local
    integer                             :: istat

    if(d%allocated) then
        call errStop('input data vector already allocated in create_dataVectorMTX')
    endif

    d%ntx = nTx

    ! allocate space for the dataVec's
    allocate(d%d(nTx), STAT=istat)

    ! not allocated until all of the dataVec's are
    d%allocated = .false.

  end subroutine create_dataVectorMTX

  !**********************************************************************
  ! deallocates all memory and reinitializes dataVectorMTX object d

  subroutine deall_dataVectorMTX(d)

    type (dataVectorMTX_t)	:: d
    ! local
    integer             :: j,istat

    if(d%allocated) then
       !  deallocate everything relevant
       do j = 1,d%nTx
          call deall_dataVector(d%d(j))
       enddo
       deallocate(d%d, STAT=istat)
    endif

    d%allocated = .false.

  end subroutine deall_dataVectorMTX

  !**********************************************************************
  ! set the data values and error bars to zero; function would make more
  ! sense here, but d = zero(d) might create trouble, so make this a
  ! subroutine instead

  subroutine zero_dataVectorMTX(d)

    type (dataVectorMTX_t), intent(inout)	:: d
    ! local
    integer                             :: j

    if(d%allocated) then
       do j = 1,d%nTx
          call zero_dataVector(d%d(j))
       enddo
    else
       call errStop('data vector not allocated in zero_dataVectorMTX')
    endif

  end subroutine zero_dataVectorMTX

  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_dataVectorMTX(d2, d1)

    type (dataVectorMTX_t), intent(in)		:: d1
    type (dataVectorMTX_t), intent(inout)	:: d2
    ! local variable
    integer                             :: j

    ! check to see if RHS (d1) is allocated
    if (.not. d1%allocated) then
       call errStop('RHS not allocated yet for copy_dataVectorMTX')
    endif

    ! check to see whether the LHS (d2) is allocated
    if (.not. d2%allocated) then
       call create_dataVectorMTX(d1%nTx, d2)
    else
       if (d1%nTx .ne. d2%nTx) then
          ! deallocate d2, and reinitialize with correct number of data vectors
          call deall_dataVectorMTX(d2)
          call create_dataVectorMTX(d1%nTx, d2)
       endif
    endif

    ! now copy the components
    do j = 1, d1%nTx
       d2%d(j) = d1%d(j)
    enddo

    d2%allocated = .true.

    ! if input is a temporary function output, deallocate
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
    ! local
    integer                                 :: j

    ! check to see if inputs (d1, d2) are both allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_dataVectorMTX')
    endif

    ! check to see if inputs are of compatible sizes
    if (d1%nTx .ne. d2%nTx) then
       call errStop('input sizes not consistent in linComb_dataVectorMTX')
    endif

    ! create the output vector that is consistent with inputs
    if (.not. dOut%allocated) then
       call errStop('output structure has to be allocated before calling linComb_dataVectorMTX')
    end if

    ! check to see that output is of compatible size
    if (d1%nTx .ne. dOut%nTx) then
       call errStop('output size not consistent in linComb_dataVectorMTX')
    endif

    ! finally do the linear combination ...
    do j = 1, d1%nTx
       call linComb_dataVector(a,d1%d(j),b,d2%d(j),dOut%d(j))
    enddo

    dOut%allocated = d1%allocated .and. d2%allocated

  end subroutine linComb_dataVectorMTX

  ! **********************************************************************
  !  computes dOut = a * dIn for dataVectorMTX objects dIn and dOut
  !  and real scalar a
  subroutine scMult_dataVectorMTX(a,d1,d2)

    real (kind=prec), intent(in)		:: a
    type (dataVectorMTX_t), intent(in)	    :: d1
    type (dataVectorMTX_t), intent(inout)  :: d2

	call linComb_dataVectorMTX(R_ZERO,d1,a,d1,d2)

  end subroutine scMult_dataVectorMTX

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
  ! the real vs complex distinction is dealt with at the level of dataVec

  function dotProd_dataVectorMTX_f(d1, d2) result(r)

    type (dataVectorMTX_t), intent(in)		:: d1, d2
    real (kind=prec)		            :: r
    ! local variables
    integer				                :: j

    ! check to see if inputs (d1, d2) are allocated
    if ((.not. d1%allocated) .or. (.not. d2%allocated)) then
       call errStop('inputs not allocated yet for dotProd_dataVectorMTX_f')
    endif

    ! check to see if inputs are compatable
    if (d1%nTx .ne. d2%nTx) then
       call errStop('input sizes not consistent in dotProd_dataVectorMTX_f')
    endif

    r = 0.0
    do j = 1, d1%nTx
       r  =  r + dotProd_dataVector_f(d1%d(j),d2%d(j))
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
   integer      			:: nn, i, j

   if (.not. d%allocated) then
      call errStop('data vector not allocated in normalize_dataVectorMTX')
   endif

   if (present(N)) then
      nn = N
   else
      nn = 1
   endif

   do j = 1, d%nTx
     do i = 1, d%d(j)%nDt
        call normalize_dataBlock(d%d(j)%data(i),nn)
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

   ! local variables
   integer 				:: i,j
   ndata = 0
   do j = 1,d%nTx
     do i = 1,d%d(j)%nDt
       ndata = ndata + d%d(j)%data(i)%nComp * d%d(j)%data(i)%nSite
     enddo
   enddo

  end function count_dataVectorMTX_f

  !**********************************************************************
  ! count the number of data vectors in the full data vector dataVectorMTX.
  ! there should be one for each transmitter and data type

  function countVec_dataVectorMTX_f(d) result(nvec)

   type(dataVectorMTX_t), intent(in)		:: d
   integer				:: nvec

   ! local variables
   integer 				:: j
   nvec = 0
   do j = 1,d%nTx
     nvec = nvec + d%d(j)%nDt
   enddo

  end function countVec_dataVectorMTX_f


end module DataSpace
