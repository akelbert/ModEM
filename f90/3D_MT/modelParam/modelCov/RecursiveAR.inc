!----------------------------------------------------------------------!
! 3D_MT model covariance by Gary Egbert and Anna Kelbert: procedures.  !
!                                                                      !
! Recursive autoregression.
!----------------------------------------------------------------------!

! *******************************************************************
  subroutine multBy_CmSqrt(m,mhat)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m
    type (modelParam_t), intent(out)    :: mhat

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m)
    end if

    mhat = m

    call RecursiveAR(m%cellCond%v, mhat%cellCond%v)

  end subroutine multBy_CmSqrt

! *******************************************************************
  subroutine create_CmSqrt(m,cfile)

  	! Initializes CmSqrt variable stored in RecursiveAR.hd. If cfile
  	! is specified, gets this information from file.

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile
    integer                             :: istat
    logical                             :: exists
    integer                             :: Nx, Ny, NzEarth

	! initializing CmSqrt
	Nx = m%Nx
	Ny = m%Ny
	NzEarth = m%NzEarth
	allocate(CmSqrt%Sx(NzEarth), STAT=istat)
	allocate(CmSqrt%Sy(NzEarth), STAT=istat)
	CmSqrt%Sx = 0.1
	CmSqrt%Sy = 0.1
	CmSqrt%Sz = 0.1
	call create_iscalar(m%grid,CmSqrt%mask,CELL_EARTH)
	CmSqrt%mask%v = 9
    CmSqrt%allocated = .true.

    if (.not. present(cfile)) then
    	! CmSqrt already initialized
    else
		! attempt to read CmSqrt from cfile
	    inquire(FILE=cfile,EXIST=exists)
	    if (exists) then
	      call read_CmSqrt(cfile)
	    else
	      call errStop('Unable to find the input covariance file '//trim(cfile)//' in create_CmSqrt')
	    end if

		if ((CmSqrt%Nx /= m%Nx) .or. (CmSqrt%Ny /= m%Ny) .or. (CmSqrt%NzEarth /= m%NzEarth)) then
	       call errStop('Grid dimensions do not match in input covariance file '//cfile)
	    end if
    end if

  end subroutine create_CmSqrt

! *******************************************************************
  subroutine deall_CmSqrt()

    integer                             :: istat

	deallocate(CmSqrt%Sx,CmSqrt%Sy, STAT=istat)
	call deall_sparsevecc(CmSqrt%S)
    call deall_iscalar(CmSqrt%mask)
    CmSqrt%allocated = .false.

  end subroutine deall_CmSqrt

! *******************************************************************
!
!...  Copyright (C) 2008 Anna Kelbert. All rights reserved.
!
  subroutine read_CmSqrt(cfile)

	! The minimal covariance information includes the AR parameters
	! alpha(k), beta(k) for smoothing in x, y directions and gamma for
	! the vertical smoothing. Both alpha and beta could depend on the
	! vertical layer. The scaling is the identity when not specified.
	! This information is read from a file. Also, we read an integer mask
	! array that subdivides the model grid into different regions
	! (AIR, OCEAN, EARTH) and a set of rules that overrides the default
	! smoothing parameters across a particular surface between two
	! distinct regions. We use this to set up the covariance CmSqrt.
	!
	! Strictly speaking, to define the smoothing across surfaces in
	! full generality while maintaining efficiency, it has to be a sparse
	! real vector defined on FACES (sparsevecr). We only have a complex
	! sparse vector implemented (sparsevecc). We could either use that,
	! or imitate the structure.

    character(*), intent(in)                         :: cfile

    ! Exception rules
    integer, pointer, dimension(:)                   :: mask1, mask2, ii, jj, kk, xyz
    real (kind=selectedPrec), pointer, dimension(:)  :: smoothing, S

    integer       			 :: Nx, Ny, NzEarth, nrules, nS, i, j, k, n, istat
    integer                  :: fid = 30

	if (.not. CmSqrt%allocated) then
		call errStop('Model covariance must be allocated before reading from file in read_CmSqrt')
	end if

    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip the 15 lines header
    do j = 1,15
       read(fid,*)
    end do

	! read grid dimensions
	read(fid,*) Nx,Ny,NzEarth
	CmSqrt%Nx = Nx
	CmSqrt%Ny = Ny
	CmSqrt%NzEarth = NzEarth

	! read smoothing parameters
    read(fid,*) CmSqrt%Sx
    read(fid,*) CmSqrt%Sy
    read(fid,*) CmSqrt%Sz

	! read exception rules for smoothing across surfaces
    read(fid,*) nrules
    allocate(mask1(nrules),mask2(nrules),smoothing(nrules),STAT=istat)
    do n = 1,nrules
       read(fid,*) mask1(n), mask2(n), smoothing(n)
    end do

	! create and read the mask array
	call read_iscalar(fid,CmSqrt%mask)

	close(fid)

	! create a huge sparse vector to make sure we accommodate all smoothing exceptions
    call create_sparsevecc(Nx*Ny*NzEarth, CmSqrt%S, FACE)

	! now, parse the exceptions
	nS = 0
    do k = 2,NzEarth
    	do j = 2,Ny
    		do i = 2,Nx
    			do n = 1,nrules
    				! look back in the X-direction
    				if (((CmSqrt%mask%v(i-1,j,k) == mask1(n)) .and. (CmSqrt%mask%v(i,j,k) == mask2(n))) &
    					.or. ((CmSqrt%mask%v(i-1,j,k) == mask2(n)) .and. (CmSqrt%mask%v(i,j,k) == mask1(n)))) &
    					then
    					nS = nS+1
    					CmSqrt%S%i(nS) = i-1
    					CmSqrt%S%j(nS) = j
    					CmSqrt%S%k(nS) = k
    					CmSqrt%S%xyz(nS) = 1
    					CmSqrt%S%c(nS) = smoothing(n)
    				end if
    				! look back in the Y-direction
    				if (((CmSqrt%mask%v(i,j-1,k) == mask1(n)) .and. (CmSqrt%mask%v(i,j,k) == mask2(n))) &
    					.or. ((CmSqrt%mask%v(i,j-1,k) == mask2(n)) .and. (CmSqrt%mask%v(i,j,k) == mask1(n)))) &
    					then
    					nS = nS+1
    					CmSqrt%S%i(nS) = i
    					CmSqrt%S%j(nS) = j-1
    					CmSqrt%S%k(nS) = k
    					CmSqrt%S%xyz(nS) = 2
    					CmSqrt%S%c(nS) = smoothing(n)
    				end if
   					! look back in the Z-direction
    				if (((CmSqrt%mask%v(i,j,k-1) == mask1(n)) .and. (CmSqrt%mask%v(i,j,k) == mask2(n))) &
    					.or. ((CmSqrt%mask%v(i,j,k-1) == mask2(n)) .and. (CmSqrt%mask%v(i,j,k) == mask1(n)))) &
    					then
    					nS = nS+1
    					CmSqrt%S%i(nS) = i
    					CmSqrt%S%j(nS) = j
    					CmSqrt%S%k(nS) = k-1
    					CmSqrt%S%xyz(nS) = 3
    					CmSqrt%S%c(nS) = smoothing(n)
    				end if
    			end do
    		end do
    	end do
    end do
	deallocate(mask1,mask2,smoothing)

    ! now, truncate the smoothing vector to the correct number of components
    call reall_sparsevecc(nS, CmSqrt%S)

  end subroutine read_CmSqrt

! *******************************************************************
!
!...  Copyright (C) 2008 Anna Kelbert. All rights reserved.
!
  subroutine RecursiveAR(w,v)

   ! Implements the recursive autoregression algorithm for a 3D real array.
   ! In our case, the assumed-shape array would be e.g. conductivity
   ! in each cell of the Nx x Ny x NzEarth grid.

    real (kind=selectedPrec), intent(in)     :: w(:,:,:)
    real (kind=selectedPrec), intent(out)    :: v(:,:,:)
    integer                                  :: Nx, Ny, NzEarth, i, j, k

    Nx      = size(w,1)
	Ny      = size(w,2)
	NzEarth = size(w,3)

	if (maxval(abs(shape(w) - shape(v)))>0) then
		call errStop('The input arrays should be of the same shapes in RecursiveAR')
	end if

	! smooth in the X-direction (Sx)
    do k = 1,NzEarth
    	do j = 1,Ny
    		v(1,j,k) = w(1,j,k)
    		do i = 2,Nx
				v(i,j,k) = SmoothX(i-1,j,k) * v(i-1,j,k) + w(i,j,k)
    		end do
    	end do
    end do

	! smooth in the Y-direction (Sy)
    do k = 1,NzEarth
    	do i = 1,Nx
    		! v(i,1,k) = v(i,1,k)
    		do j = 2,Ny
				v(i,j,k) = SmoothY(i,j-1,k) * v(i,j-1,k) + v(i,j,k)
    		end do
    	end do
    end do

	! smooth in the Z-direction (Sz)
    do j = 1,Ny
    	do i = 1,Nx
    		! v(i,j,1) = v(i,j,1)
    		do k = 2,NzEarth
				v(i,j,k) = SmoothZ(i,j,k-1) * v(i,j,k-1) + v(i,j,k)
    		end do
    	end do
    end do

	! apply the scaling operator C
    do k = 1,NzEarth
    	do j = 1,Ny
    		do i = 1,Nx
				v(i,j,k) = Scaling(i,j,k) * v(i,j,k)
    		end do
    	end do
    end do

  end subroutine RecursiveAR

! ****************************************************************************
  function SmoothX(i,j,k) result(alpha)

  	! computes the smoothing coefficient in the x-direction based on CmSqrt

    implicit none
    integer, intent(in)                      :: i,j,k
    real (kind=selectedPrec)                 :: alpha

    ! local variables
    integer                                  :: n

    if (.not. associated(CmSqrt%Sx)) then
    	call errStop('CmSqrt%Sx has to be allocated before calling SmoothX')
    end if

    alpha = CmSqrt%Sx(k)

    if ((i < 1) .or. (i > CmSqrt%Nx)) then
    	call errStop('index i out of bounds in SmoothX(i,j,k)')
    else if ((j < 1) .or. (j > CmSqrt%Ny)) then
    	call errStop('index j out of bounds in SmoothX(i,j,k)')
    else if ((k < 1) .or. (k > CmSqrt%NzEarth)) then
    	call errStop('index k out of bounds in SmoothX(i,j,k)')
    end if

    if (CmSqrt%S%allocated .and. (CmSqrt%S%nCoeff > 0)) then
    	! scan through the special rules and possibly update the result
    	do n = 1,CmSqrt%S%nCoeff
    		if (CmSqrt%S%xyz(n) == 1) then
    			if ((CmSqrt%S%i(n) == i) .and. (CmSqrt%S%j(n) == j) .and. (CmSqrt%S%k(n) == k)) then
    				alpha = CmSqrt%S%c(n)
    				exit
    			end if
    		end if
    	end do
    end if

  end function SmoothX

! ****************************************************************************
  function SmoothY(i,j,k) result(beta)

  	! computes the smoothing coefficient in the y-direction based on CmSqrt

    implicit none
    integer, intent(in)                      :: i,j,k
    real (kind=selectedPrec)                 :: beta

    ! local variables
    integer                                  :: n

    if (.not. associated(CmSqrt%Sy)) then
    	call errStop('CmSqrt%Sy has to be allocated before calling SmoothY')
    end if

    beta = CmSqrt%Sy(k)

    if ((i < 1) .or. (i > CmSqrt%Nx)) then
    	call errStop('index i out of bounds in SmoothY(i,j,k)')
    else if ((j < 1) .or. (j > CmSqrt%Ny)) then
    	call errStop('index j out of bounds in SmoothY(i,j,k)')
    else if ((k < 1) .or. (k > CmSqrt%NzEarth)) then
    	call errStop('index k out of bounds in SmoothY(i,j,k)')
    end if

    if (CmSqrt%S%allocated .and. (CmSqrt%S%nCoeff > 0)) then
    	! scan through the special rules and possibly update the result
    	do n = 1,CmSqrt%S%nCoeff
    		if (CmSqrt%S%xyz(n) == 2) then
    			if ((CmSqrt%S%i(n) == i) .and. (CmSqrt%S%j(n) == j) .and. (CmSqrt%S%k(n) == k)) then
    				beta = CmSqrt%S%c(n)
    				exit
    			end if
    		end if
    	end do
    end if

  end function SmoothY

! ****************************************************************************
  function SmoothZ(i,j,k) result(gamma)

  	! computes the smoothing coefficient in the z-direction based on CmSqrt

    implicit none
    integer, intent(in)                      :: i,j,k
    real (kind=selectedPrec)                 :: gamma

    ! local variables
    integer                                  :: n

    gamma = CmSqrt%Sz

    if ((i < 1) .or. (i > CmSqrt%Nx)) then
    	call errStop('index i out of bounds in SmoothZ(i,j,k)')
    else if ((j < 1) .or. (j > CmSqrt%Ny)) then
    	call errStop('index j out of bounds in SmoothZ(i,j,k)')
    else if ((k < 1) .or. (k > CmSqrt%NzEarth)) then
    	call errStop('index k out of bounds in SmoothZ(i,j,k)')
    end if

    if (CmSqrt%S%allocated .and. (CmSqrt%S%nCoeff > 0)) then
    	! scan through the special rules and possibly update the result
    	do n = 1,CmSqrt%S%nCoeff
    		if (CmSqrt%S%xyz(n) == 3) then
    			if ((CmSqrt%S%i(n) == i) .and. (CmSqrt%S%j(n) == j) .and. (CmSqrt%S%k(n) == k)) then
    				gamma = CmSqrt%S%c(n)
    				exit
    			end if
    		end if
    	end do
    end if

  end function SmoothZ

! ****************************************************************************
  function Scaling(i,j,k) result(c)

  	! computes the scaling coefficient based on CmSqrt

    implicit none
    integer, intent(in)                      :: i,j,k
    real (kind=selectedPrec)                 :: c

    if ((i < 1) .or. (i > CmSqrt%Nx)) then
    	call errStop('index i out of bounds in Scaling(i,j,k)')
    else if ((j < 1) .or. (j > CmSqrt%Ny)) then
    	call errStop('index j out of bounds in Scaling(i,j,k)')
    else if ((k < 1) .or. (k > CmSqrt%NzEarth)) then
    	call errStop('index k out of bounds in Scaling(i,j,k)')
    end if

	if (CmSqrt%mask%v(i,j,k) == AIR) then
		c = 0.0
	else if (CmSqrt%mask%v(i,j,k) == OCEAN) then
		c = 0.0
	else
		c = 1.0
	end if

  end function Scaling
