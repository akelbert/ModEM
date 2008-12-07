!*****************************************************************************
module dataspace
  ! Routines for the generic "data space". Modules in this group are
  ! indifferent to all details of the modeling code, including the
  ! grid.

  !   Dec 2005: added possible storage for error bars

  use utilities
  use math_constants
  implicit none

  interface assignment (=)
     MODULE PROCEDURE copy_Dvec
     MODULE PROCEDURE copy_DvecMTX
  end interface
  
  interface dotProd
     MODULE PROCEDURE dotProd_Dvec
     MODULE PROCEDURE dotProd_DvecMTX
  end interface

  interface operator (*)
     MODULE PROCEDURE scMult_DvecMTX
  end interface

  !  vector containing data for a single "transmitter"
  type :: dvec
   SEQUENCE

      ! nComp is number of EM components observed (e.g., 2 (3) for
      ! MT (with verticl field TFs)) times 2 to allow for real/imag
      ! nSite is number of sites where these components are observed
      integer 			:: nComp = 0
      integer			:: nSite = 0

      ! data will contain actual data; 
      ! dimensions will be data(nComp,nSite)
      real (kind=selectedPrec), pointer, dimension(:,:) 	:: data,err

      ! rx will contain indices of corresponding entries in
      ! receiver definition dictionary (dimension: rx(nSite))
      integer, pointer, dimension(:)	:: rx
      !  tx is index into transmitter dictionary
      integer			:: tx = 0

      ! dataType will serve to identify rx/tx dictionaries,
      ! and instances of evaluation, comb functionals
      ! This integer will ultimately be an index into a "dataType"
      ! dictionary ...With only a single data type this is not
      ! really needed. Note that associated mode definitions for
      ! transmitter can be given in dataType dictionary.  
      integer		:: dataType = 0

      logical		:: allocated = .false.
      logical		:: errorBar = .false.
     
  end type dvec

  ! collection of dvec objects for all transmitters
  type :: dvecMTX
   SEQUENCE
      ! nTX is number of transmitters, number of frequencies for MT
      integer		:: nTX = 0
      ! Ndata is total number of REAL data (dimension of data vector)
      integer		:: Ndata = 0

      ! d is array of dvec's for each transmitter (dimension nTX)
      type(dvec), pointer, dimension(:)	:: d

      logical		:: allocated = .false.
      logical		:: errorBar = .false.
      logical		:: normalized = .false.

  end type dvecMTX

  public			:: create_Dvec, deall_Dvec, deall_DvecMTX
  public			:: copy_Dvec, copy_DvecMTX
  public			:: linComb_Dvec, linComb_DvecMTX
  public			:: dotProd_Dvec, dotProd_DvecMTX
  public			:: count_DvecMTX, SSdiff_DvecMTX
  public            :: normalize_DvecMTX, normalize2_DvecMTX
  public            :: scMult_DvecMTX

Contains

  ! **********************************************************************
  ! creates an object of type dvec. a vector containing data for a single
  ! transmitter. nComp: number of components, nSite: number of sites.  
  ! --> By default DataVec%rx is initialized to site #
  ! --> Set DataVec%errorBar = .true. to allocate storage for error bars also
  ! 
  subroutine create_Dvec(nComp, nSite, d, errBar, iTx, iDt )

    integer, intent(in)			:: nComp, nSite
    type (dvec), intent(inout)		:: d
    logical, intent(in), optional	:: errBar
    integer, intent(in), optional	:: iTx, iDt

    if(d%allocated) then
        call deall_Dvec(d)
    endif
       
    allocate(d%data(nComp, nSite)) 
    d%data = R_ZERO
    allocate(d%rx(nSite)) 
    if(present(errBar)) then
       d%errorBar = errBar
    endif
    if(present(iTx)) then
       d%tx = iTx
    endif
    if(present(iDt)) then
       d%dataType = iDt
    endif

    if(d%errorBar) then
       allocate(d%err(nComp, nSite)) 
       d%err = R_ZERO
    endif

    d%allocated = .true.
    d%nComp = nComp
    d%nSite = nSite

  end subroutine create_Dvec

  !**********************************************************************
  subroutine deall_Dvec(d)

    !  deallocates all memory and reinitialized dvec object d
    ! NOTE: do not base your judgement on the deallocation of d%err
    ! on d%errorBar variable: it can be set to .false. somewhere in
    ! the program without deallocating d%err.

    type (dvec), intent(inout)	:: d
    integer  istat

    if(d%allocated) then
       !  deallocate everything relevant
       deallocate(d%data,d%rx)
       nullify(d%data,d%rx)
       if (associated(d%err)) deallocate(d%err, STAT=istat)
       !if(d%errorBar) then
       !   deallocate(d%err)
       !   nullify(d%err)
       !endif
    endif
    d%allocated = .false.
    d%tx = 0
    d%dataType = 0
    d%ncomp = 0
    d%nSite = 0

  end subroutine deall_Dvec

  !**********************************************************************
  subroutine deall_DvecMTX(D)

    type (dvecMTX), intent(inout)	:: D

    !  local variables
    integer				:: j

    if(D%allocated) then
       do j = 1,D%nTx
          call deall_Dvec(D%d(j))
       enddo
       deallocate(D%d)
       nullify(D%d)
    endif

    D%nTx = 0
    D%Ndata = 0
    D%allocated = .false.

  end subroutine deall_DvecMTX
    
  ! **********************************************************************
  ! copy a data vector from d1 to d2 ...
  ! interface to =
  ! check for size consistency, reallocate output if needed

  subroutine copy_Dvec(d2,d1)

    type (dvec), intent(in)		:: d1
    type (dvec), intent(inout)		:: d2

    ! check to see if RHS (d1) is active (allocated)
    if(.not.d1%allocated) then
       call errStop('RHS not allocated yet for copy_Dvec')
    else
       ! check to see if everything is same
       if((d1%nComp.ne.d2%nComp).or.(d1%nSite.ne.d2%nSite) .or. &
          (d1%errorBar.neqv.d2%errorBar)) then
          !  deallocate d2 (if necessary), and reinitialize with corrct shape
          call deall_Dvec(d2)
          !  allocate SV2 as correct size ...
          call create_Dvec(d1%nComp, d1%nSite, d2,d1%errorBar)
       end if
       ! just copy the components
       d2%data = d1%data
       d2%rx = d1%rx
       d2%tx = d1%tx
       d2%dataType = d1%dataType
       if(d2%errorBar) then
          d2%err = d1%err
       endif
    end if
  end subroutine copy_Dvec

  ! **********************************************************************
  ! this will copy a collection of data vectors from D1 to D2 for all 
  ! transmitters...interface to =
  ! check for size consistency, reallocate output if needed
  subroutine copy_dvecMTX(D2,D1)

    type (dvecMTX), intent(in)		:: D1
    type (dvecMTX), intent(inout)	:: D2

    ! local variables
    integer				:: i

    ! check to see if RHS (D1) is allocated
    if(.not.D1%allocated) then
       call errStop('RHS not allocated yet for copy_dvecMTX')
    else
       ! check to see if number of TX allocated are same for D1, D2
       if(D1%nTX .ne. D2%nTX) then
          call deall_DvecMTX(D2)
          allocate(D2%d(D1%nTx))
          D2%nTX = D1%nTX
       endif
       D2%errorBar = D1%errorBar
       D2%normalized = D1%normalized
       D2%allocated = D1%allocated
       D2%Ndata = D1%Ndata
       do i =1, D1%nTX 
          Call copy_Dvec(D2%d(i),D1%d(i))
       end do
     endif
  end subroutine copy_dvecMTX
  
  ! **********************************************************************
  ! calculates linear combination of two dvec objects a*d1+ b*d2 
  !   What sort of overwriting is allowed?????
  ! Note that the error bars are allocated and the errorBar variable
  ! initialized in create_Dvec. Use the errors from one of the inputs;
  !   Error bars are not now defined when both input vectors have errors

    subroutine linComb_Dvec(a,d1,b,d2,dOut)

    type (dvec), intent(in)			:: d1, d2
    real (kind=selectedPrec), intent(in)	:: a, b
    type (dvec), intent(inout)			:: dOut

    ! local variables
    logical					:: errBar = .false.
 
    ! check to see if inputs (d1, d2) are both allocated
    if((.not.d1%allocated).or.(.not.d2%allocated)) then
       call errStop('inputs not allocated on call to linComb_Dvec')
    endif

    ! check to see if inputs are compatable
    if((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite)) then
       call errStop('input dvecs not consistent in linComb_Dvec')
    endif

	! set errBar=.true. if at least one of the inputs has error bars
	
	errBar = (d1%errorBar .or. d2%errorBar)

    ! create the output vector that is consistent with inputs

    if(dOut%allocated) then
       call deall_Dvec(dOut) 
    end if

    call create_Dvec(d1%nComp, d1%nSite, dOut, errBar, d1%tx, d1%dataType )

	! set the receiver indices to those of d1
	dOut%rx = d1%rx

    !  finally do the linear combination ...
    dOut%data = a*d1%data+b*d2%data
    
	! deal with the error bars by copying them from one of the vectors;
	!   currently exit if both of the input vectors have error bars defined
	!   unless one of the multiplication coefficients is zero
	
	if(d1%errorBar .and. d2%errorBar) then
	   if ((abs(a) > R_ZERO).and.(abs(b) > R_ZERO)) then 
       	call errStop('unable to add two data vectors with error bars in linComb_Dvec')
       else
       	dOut%err = a*d1%err+b*d2%err
       end if
    else if(d1%errorBar) then
       dOut%err = a*d1%err
    else if(d2%errorBar) then
       dOut%err = b*d2%err       
    end if	

  end subroutine linComb_Dvec

  ! **********************************************************************
  ! calculates linear combination of multi-transmitter data vector objects
  !  a*D1+ b*D2

    subroutine linComb_DvecMTX(a,D1,b,D2,Dout)

    type (dvecMTX), intent(in)			:: D1, D2
    real (kind=selectedPrec), intent(in)			:: a, b
    type (dvecMTX), intent(inout)		:: Dout

    !  local variables
    integer					:: j

    ! check to see if inputs (D1, D2) are both allocated
    if((.not.D1%allocated).or.(.not.D2%allocated)) then
       call errStop('inputs not allocated on call to linComb_DvecMTX')
    endif

    ! check to see if number of transmitters is the same
    if(D1%nTX .ne. D2%nTX) then
       call errStop('Size of inputs D1, D2 not compatbile in linComb_DvecMTX')
    endif

    !  check for allocation, sizes of Dout
    if(Dout%allocated) then
       if(Dout%nTx .ne. D1%nTx) then
          call deall_DvecMTX(Dout)
          allocate(Dout%d(D1%nTx))
          Dout%nTX = D1%nTX
       endif
    else
       allocate(Dout%d(D1%nTx))
       Dout%nTX = D1%nTX
    endif
    Dout%normalized = .false.
           
    ! Form linear combination of each dvec object (rely on error checking,
    !   possible allocation, internal to LinCombDvec)
    Dout%errorBar = .true.
    do j = 1, D1%nTX
       call linComb_Dvec(a, D1%d(j), b, D2%d(j), Dout%d(j))
       if (.not.Dout%d(j)%errorBar) Dout%errorBar = .false.
    enddo
    Dout%Ndata = D1%Ndata
    Dout%allocated = .true.
	  
  end subroutine linComb_dvecMTX
  
  ! **********************************************************************
  subroutine scMultAdd_Dvec(a,d1,d2)
  !  computes d2+a*d1 for Dvec objects d1, d2 and real scalar a,
  !   overwriting d2; ignores error field

    type (dvec), intent(in)			:: d1
    real (kind=selectedPrec), intent(in)	:: a
    type (dvec), intent(inout)			:: d2

    ! local variables
    logical					:: errBar = .false.

    ! check to see if inputs (d1, d2) are both allocated
    if((.not.d1%allocated).or.(.not.d2%allocated)) then
       call errStop('inputs not allocated on call to scMultAdd_Dvec')
    endif

    ! check to see if inputs are compatable
    if((d1%nComp .ne. d2%nComp) .or. (d1%nSite .ne. d2%nSite)) then
       call errStop('input dvecs not consistent in scMultAdd_Dvec')
    endif

    !  finally do the linear combination ...
    d2%data = a*d1%data+d2%data

  end subroutine scMultAdd_Dvec

  ! **********************************************************************
    subroutine scMultAdd_DvecMTX(a,D1,D2)
  !  computes D2+a*D1 for DvecMTX objects D1, D2 and real scalar a,
  !   overwriting D2; ignores error field

    type (dvecMTX), intent(in)		:: D1
    real (kind=selectedPrec), intent(in)		:: a
    type (dvecMTX), intent(inout)	:: D2

    !  local variables
    integer				:: j

    ! check to see if inputs (D1, D2) are both allocated
    if((.not.D1%allocated).or.(.not.D2%allocated)) then
       call errStop('inputs not allocated on call to scMultAdd_DvecMTX')
    endif

    ! check to see if number of transmitters is the same
    if(D1%nTX .ne. D2%nTX) then
       call errStop('Size of inputs D1, D2 not compatbile in scMultAddDvecMTX')
    endif

    ! Form linear combination of each dvec object (rely on error checking,
    !   possible allocation, internal to LinCombDvec)
    do j = 1, D1%nTX
       call scMultAdd_Dvec(a, D1%d(j), D2%d(j))
    enddo
	  
  end subroutine scMultAdd_DvecMTX
  
  ! **********************************************************************
    function scMult_DvecMTX(a,dIn) result (dOut)
  !  computes dOut = a * dIn for DvecMTX objects dIn and dOut 
  !  and real scalar a; ignores error field

    real (kind=selectedPrec), intent(in)		:: a
    type (dvecMTX), intent(in)	                :: dIn
    type (dvecMTX)                              :: dOut

    ! check to see that input d is allocated
    if(.not.dIn%allocated) then
       call errStop('input not allocated on call to scDivide_DvecMTX')
    endif

	call linComb_DvecMTX(R_ZERO,dIn,a,dIn,dOut)
	  
  end function scMult_DvecMTX
  
  ! **********************************************************************
  ! dot product of two real dvec objevcts r  = <d1,d2>

    function dotProd_Dvec(d1,d2) result(r)

    type (dvec), intent(in)		:: d1, d2
    real (kind=selectedPrec)		:: r

    ! local variables
    integer				:: j, k

    ! check to see if inputs (d1, d2) are allocated
    if((.not.d1%allocated).or.(.not.d2%allocated)) then
       call errStop('RHS not allocated yet for dotProd_Dvec')
    endif

    ! check to see if array sizes are the same
    if((d1%nComp.ne.d2%nComp).or.(d1%nSite.ne.d2%nSite)) then
       call errStop('size of d1, d2 not compatible in dotProd_Dvec')
    endif

    r = 0.0
    do j = 1, d1%nComp
       do k = 1, d1%nSite
          r  =  r + d2%data(j,k) * d1%data(j,k)
       end do
    end do
	  
  end function dotProd_Dvec


  ! **********************************************************************
  ! Dot product for multi-transmitter data vector object
    function dotProd_DvecMTX(D1,D2) result(r)

    type (dvecMTX), intent(in)		:: D1, D2
    real (kind=selectedPrec)		:: r

    ! local variables
    integer				:: j

    ! check to see if inputs (D1, D2) are allocated
    if((.not.D1%allocated).or.(.not.D2%allocated)) then
       call errStop('RHS not allocated yet for dotProd_DvecMTX')
    endif

    ! check to see if array sizes are the same
    if((D1%nTx.ne.D2%nTx).or.(D1%Ndata.ne.D2%Ndata)) then
       call errStop('size of D1, D2 not compatible in dotProd_DvecMTX')
    endif

    r = 0.0
    do j = 1, D1%nTX
       r =  r + dotProd_Dvec(D1%d(j), D2%d(j))
    end do
	  
  end function dotProd_DvecMTX

!**********************************************************************
   subroutine normalize_dvecMTX(d,dNorm)
   !  normalizes a dvecMTX object using error bars
   !   (Add attribute "normalized" to dvecMTX, dvec objects)
   !  dNorm is optional output; if not present d is overwritten
   !  Output dNorm must be allocated before calling 

   type(dvecMTX),intent(inout)          :: d
   type(dvecMTX),intent(inout),optional :: dNorm

   !  local variables
   integer      			:: nTx, j

   nTx = d%nTx

   if(.not.d%allocated) then
      call errStop('data vector not allocated in normalize_dvecMTX')
   endif
   if(.not.d%errorBar) then
      call errStop('no error bars for input data in normalize_dvecMTX')
   endif

   if(present(dNorm)) then
      !  make copy of input d into dNorm
      call copy_dvecMTX(dNorm,d)
      if(.not.d%normalized) then
         dNorm%normalized = .true.
         do j = 1,nTx
            dNorm%d(j)%data = d%d(j)%data/d%d(j)%err
         enddo
      endif
   else
      !  normalize input
      if(.not.d%normalized) then
         d%normalized = .true.
         do j = 1,nTx
            d%d(j)%data = d%d(j)%data/d%d(j)%err
         enddo
      endif
   endif

   end subroutine normalize_dvecMTX

!**********************************************************************
   subroutine normalize2_dvecMTX(d,dNorm)
   !  normalizes a dvecMTX object using error bars squared
   !   (Add attribute "normalized" to dvecMTX, dvec objects)
   !  dNorm is optional output; if not present d is overwritten
   !  Output dNorm must be allocated before calling 

   type(dvecMTX),intent(inout)          :: d
   type(dvecMTX),intent(inout),optional :: dNorm

   !  local variables
   integer      			:: nTx, j

   if(.not.d%allocated) then
      call errStop('data vector not allocated in normalize2_dvecMTX')
   endif
   if(.not.d%errorBar) then
      call errStop('no error bars for input data in normalize2_dvecMTX')
   endif

   nTx = d%nTx

   if(present(dNorm)) then
      !  make copy of input d into dNorm
      call copy_dvecMTX(dNorm,d)
      if(.not.d%normalized) then
         dNorm%normalized = .true.
         do j = 1,nTx
				 ! WRONG NEEDS TO BE CORRECTED: data(nComp,nSite)
            dNorm%d(j)%data = d%d(j)%data/(d%d(j)%err**2)
         enddo
      endif
   else
      !  normalize input
      if(.not.d%normalized) then
         d%normalized = .true.
         do j = 1,nTx
            d%d(j)%data = d%d(j)%data/(d%d(j)%err**2)
         enddo
      endif
   endif

   end subroutine normalize2_dvecMTX

!**********************************************************************
   subroutine setError_dvec(err,d)
   !  whether error bars exist or not, sets new error bars
	 !  according to the input parameter err, which specifies the
	 !  fractional error in the output; any errors stored in d
	 !  are overwritten

	 real (kind=selectedPrec), intent(in) :: err
   type(dvec),intent(inout)             :: d

   !  local variables
   integer      			:: nComp, nSite

   if(.not.d%allocated) then
      call errStop('data vector not allocated in setError_dvec')
   endif

   nComp = d%nComp
	 nSite = d%nSite
	 if (.not.d%errorBar) then
     allocate(d%err(nComp,nSite))
		 d%errorBar = .true.
	 end if
	 
   ! we should really be using dataType%isComplex information
	 ! create data errors for complex data, but we currently do not
	 ! have access to the dataType dictionary from here (I think);
	 ! so this is a temporary solution

	 d%err = err * abs(d%data)
   
   end subroutine setError_dvec

!**********************************************************************
   subroutine setError_dvecMTX(err,d)
   !  whether error bars exist or not, sets new error bars
	 !  according to the input parameter err, which specifies the
	 !  fractional error in the output; any errors stored in d
	 !  are overwritten

	 real (kind=selectedPrec), intent(in) :: err
   type(dvecMTX),intent(inout)          :: d

   !  local variables
   integer      			:: nTx, nData, j, k

   if(.not.d%allocated) then
      call errStop('data vector in setError_dvecMTX not allocated')
   endif

   nTx = d%nTx

   do j = 1,nTx
       call setError_dvec(err,d%d(j))
   enddo

   d%errorBar = .true.

   end subroutine setError_dvecMTX

!**********************************************************************
   function count_dvecMTX(d) result(Ndata)

   type(dvecMTX), intent(in)		:: d
   integer				:: Ndata

   ! local variables
   integer 				:: j
   Ndata = 0
   do j = 1,d%nTx
     Ndata = Ndata + d%d(j)%ncomp*d%d(j)%nSite
   enddo
   
   end function count_dvecMTX

!****************************************************************************

  subroutine SSdiff_DvecMTX(allData,allCalc,misfit)

  type (dvecMTX),intent(in)      :: allData
  type (dvecMTX),intent(in)      :: allCalc
  real (kind=selectedPrec), intent(inout)      :: misfit

  ! local variables
  integer ip
  real (kind=selectedPrec)    ::  mf2

  do ip = 1,allData%nTx
    mf2 = SUM(((allData%d(ip)%data - allCalc%d(ip)%data)/allData%d(ip)%err)**2.)
    misfit = misfit + mf2
  enddo

  end subroutine SSdiff_DvecMTX

end module dataspace
