module sensMatrix

  use math_constants
  use utilities
  use dataspace
  use modelparameter

  implicit none

  public 	:: create_sensMatrix, deall_sensMatrix

  !***********************************************************************
  ! Data type definitions for the full sensitivity matrix; a cross between
  ! a data vector and a model parameter to store the full matrix J = df/dm
  ! or the gradient of the data misfit, component-wise, i.e. (J^T)_i r_i.
  ! In this case, one column of J^T is a model parameter, so that J^T acts
  ! on a single data residual to produce the misfit gradient for that data
  ! point with respect to the model parameter; these can be summed up to
  ! obtain the full gradient. Can be useful as another measure of sensitivity.
  ! These matrices will be huge. We don't create them unless they are needed.

  ! basic data sensitivity vector for a single transmitter & data type
  type :: sensVector_t

      ! nComp is number of EM components observed (e.g., 2 (3) for
      ! MT (with verticl field TFs)) times 2 to allow for real/imag;
      ! has to match the data type index dt
      integer 		:: nComp = 0

      ! nSite is number of sites where these components are observed
      integer		:: nSite = 0

      ! sensitivities stored as model parameters; dimensions (nComp,nSite)
      type (modelParam_t), pointer, dimension(:,:) :: dm

      ! rx(nSite) is the array of indices into receiver dictionary
      integer, pointer, dimension(:) :: rx

      ! tx is index into transmitter dictionary
      integer		:: tx = 0

      ! dt is index into data type dictionary
      integer		:: dataType = 0

      logical       :: isComplex = .false.
      logical		:: allocated = .false.

      ! needed to avoid memory leaks for temporary function outputs
      logical		:: temporary = .false.

  end type sensVector_t


  ! collection of sensVector objects for all data types, single transmitter
  type :: sensMatrixTX_t
      ! the number of dataVecs (generally the number of data types) associated
      ! with this transmitter (note: not the total number of data types)
      integer       :: ndt = 0

      ! array of sensVector's, usually one for each data type (dimension ndt)
      type (sensVector_t), pointer, dimension(:)   :: v

      ! tx is the index into transmitter dictionary
      integer       :: tx = 0

      logical       :: allocated = .false.

  end type sensMatrixTX_t


  ! collection of sensVector objects for all transmitters
  type :: sensMatrix_t
      ! ntx is number of transmitters, number of frequencies for MT
      integer		:: ntx = 0

      ! d is array of dataVec's for each transmitter (dimension nTX)
      type (sensMatrixTX_t), pointer, dimension(:)	:: sens

      logical		:: allocated = .false.

  end type sensMatrix_t


Contains

  ! **********************************************************************
  ! creates an object of type sensMatrixTX:
  ! a vector containing sensitivities for a single transmitter.
  ! --> nComp: number of components, nSite: number of sites
  ! --> nDt: number of data types for transmitter tx

  subroutine create_sensMatrixTX(d, sigma0, sens)

    type (dataVecTX_t), intent(in)		 :: d
    type (modelParam_t), intent(in)		 :: sigma0
    type (sensMatrixTX_t), intent(inout) :: sens
    ! local
    integer                              :: nComp,nSite,i,j,k,istat

    if(sens%allocated) then
       call errStop('input sens matrix already allocated in create_sensMatrixTX')
    else if(.not. d%allocated) then
       call errStop('input data vector not allocated in create_sensMatrixTX')
    end if

    sens%ndt = d%ndt

    allocate(sens%v(d%ndt), STAT=istat)

    do i = 1,d%nDt

       nComp = d%data(i)%nComp
       nSite = d%data(i)%nSite

       sens%v(i)%nComp = nComp
       sens%v(i)%nSite = nSite

       allocate(sens%v(i)%dm(nComp,nSite), STAT=istat)
       do j = 1,nComp
         do k = 1,nSite
         	sens%v(i)%dm(j,k) = sigma0
         end do
       end do

       sens%v(i)%tx = d%data(i)%tx
       allocate(sens%v(i)%rx(nSite), STAT=istat)
       sens%v(i)%rx = d%data(i)%rx

       sens%v(i)%dataType = d%data(i)%dataType
       sens%v(i)%isComplex = d%data(i)%isComplex

       sens%v(i)%allocated = .true.

    end do

    sens%tx = d%tx
    sens%allocated = .true.

  end subroutine create_sensMatrixTX

  !**********************************************************************
  ! deallocates all memory and reinitialized sensMatrixTX object sens

  subroutine deall_sensMatrixTX(sens)

    type (sensMatrixTX_t)	:: sens
    ! local
    integer             	:: nComp,nSite,i,j,k,istat

    if(sens%allocated) then
       !  deallocate everything relevant
       do i = 1,sens%nDt

          sens%v(i)%tx = 0
          deallocate(sens%v(i)%rx, STAT=istat)

          sens%v(i)%dataType = 0
          sens%v(i)%isComplex = .false.

          do j = 1,nComp
            do k = 1,nSite
            	call deall_modelParam(sens%v(i)%dm(j,k))
            end do
          end do
          deallocate(sens%v(i)%dm, STAT=istat)

          sens%v(i)%nComp = 0
          sens%v(i)%nSite = 0
          sens%v(i)%allocated = .false.

       end do
       deallocate(sens%v, STAT=istat)
    end if

    sens%tx = 0
    sens%allocated = .false.

  end subroutine deall_sensMatrixTX


  ! **********************************************************************
  ! creates an object of type sensMatrix:
  ! a vector containing sensitivities for all transmitters.

  subroutine create_sensMatrix(d, sigma0, J)

    type (dataVecMTX_t), intent(in)  	:: d
    type (modelParam_t), intent(in)  	:: sigma0
    type (sensMatrix_t), intent(inout)  :: J
    ! local
    integer                             :: i,istat

    if(J%allocated) then
        call errStop('sensitivity matrix already allocated in create_sensMatrix')
    end if

    J%ntx = d%ntx

    ! allocate space for the sensitivity matrix
    allocate(J%sens(d%ntx), STAT=istat)
    do i = 1,d%ntx
		call create_sensMatrixTX(d%d(i), sigma0, J%sens(i))
	end do

    ! memory allocation complete
    J%allocated = .true.

  end subroutine create_sensMatrix

  !**********************************************************************
  ! deallocates all memory and reinitialized the full sensitivity matrix

  subroutine deall_sensMatrix(J)

    type (sensMatrix_t)	:: J
    ! local
    integer             :: i,istat

    if(J%allocated) then
       !  deallocate everything relevant
       do i = 1,J%ntx
          call deall_sensMatrixTX(J%sens(i))
       end do
       deallocate(J%sens, STAT=istat)
    endif

    J%ntx = 0
    J%allocated = .false.

  end subroutine deall_sensMatrix

end module sensMatrix
