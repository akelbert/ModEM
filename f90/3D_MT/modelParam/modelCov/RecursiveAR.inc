!----------------------------------------------------------------------!
! 3D_MT model covariance by Gary Egbert and Anna Kelbert: procedures.  !
!                                                                      !
! Recursive autoregression.
!----------------------------------------------------------------------!

! *******************************************************************
  subroutine multBy_CmSqrt(m1,m2)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t), intent(out)    :: m2

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if

    m2 = m1
    !call SmoothX(CmSqrt%Sx,m2)
    !call SmoothY(CmSqrt%Sy,m2)
    !call SmoothZ(CmSqrt%Sz,m2)
    !call Scaling(CmSqrt%C,m2)

  end subroutine multBy_CmSqrt

! *******************************************************************
  subroutine create_CmSqrt(m,cfile)

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile
    logical                             :: exists

    call create_rscalar(m%grid,CmSqrt%Sx,CELL_EARTH)
    CmSqrt%Sx%v = 0.1

    call create_rscalar(m%grid,CmSqrt%Sy,CELL_EARTH)
    CmSqrt%Sy%v = 0.1

    call create_rscalar(m%grid,CmSqrt%Sz,CELL_EARTH)
    CmSqrt%Sz%v = 0.1

    call create_rscalar(m%grid,CmSqrt%C,CELL_EARTH)
    CmSqrt%C%v = 1.0

    if (present(cfile)) then
      inquire(FILE=cfile,EXIST=exists)
      if (exists) then
        call read_CmSqrt(cfile,m%grid)
      else
        call errStop('Unable to find the input covariance file '//trim(cfile)//' in create_CmSqrt')
      end if
    end if

    CmSqrt%allocated = .true.

  end subroutine create_CmSqrt

! *******************************************************************
  subroutine deall_CmSqrt()

    call deall_rscalar(CmSqrt%Sx)
    call deall_rscalar(CmSqrt%Sy)
    call deall_rscalar(CmSqrt%Sz)
    call deall_rscalar(CmSqrt%C)

  end subroutine deall_CmSqrt

! *******************************************************************
!
!...  Copyright (C) 2008 Anna Kelbert. All rights reserved.
!
  subroutine read_CmSqrt(cfile,grid)

    character(*), intent(in)           :: cfile
    type (grid3d_t), intent(in)        :: grid

	! The minimal covariance information includes the AR parameters
	! alpha(k), beta(k) for smoothing in x, y directions and gamma for
	! the vertical smoothing. Both alpha and beta could depend on the
	! vertical layer. The scaling is the identity when not specified.
    real (kind=selectedPrec), pointer, dimension(:)         ::  alpha, beta
    real (kind=selectedPrec)                                ::  gamma

    ! The mask is an integer array which indicates which parameters to freeze.
    ! Example:  0 0 0 0 0 0 0 0 0 2 2 0 0 0 0 0
    !           0 0 0 0 0 0 0 0 2 2 2 2 2 2 0 0
    !           1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2
    !           1 1 1 1 1 1 2 2 2 2 2 3 3 3 2 2
    !           1 1 1 1 1 2 2 2 2 2 2 3 3 3 2 2
    !           1 1 1 1 2 2 2 2 2 3 2 3 3 3 2 2
    !           1 1 1 4 4 4 4 4 4 4 4 4 4 4 4 4
    !           4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    ! Above is a 2D model with ocean (ocean == 1) and topography (air == 0).
    ! We might wish to switch off the smoothing between the air & ocean and
    ! the rest of the model; or  between the areas 3 and 4 only; or both.
    type (iscalar)                      :: mask

    ! Exception rules
    integer, pointer, dimension(:)                         :: mask1, mask2
    real (kind=selectedPrec), pointer, dimension(:)        :: smoothing

    integer                             :: Nx, Ny, NzEarth, nrules, j, k, istat, k1, k2
    integer                             :: fid = 30

    open(unit=fid,file=cfile,form='formatted',status='old')

    ! skip the 15 lines header
    do j = 1,15
       read(fid,*)
    end do

	read(fid,*) Nx,Ny,NzEarth

	if ((Nx /= grid%Nx) .or. (Ny /= grid%Ny) .or. (NzEarth /= grid%NzEarth)) then
       call errStop('Grid dimensions do not match in input covariance file '//cfile)
    end if

	allocate(alpha(NzEarth),beta(NzEarth),STAT=istat)

    read(fid,*) alpha
    read(fid,*) beta
    read(fid,*) gamma

    read(fid,*) nrules
    allocate(mask1(nrules),mask2(nrules),smoothing(nrules),STAT=istat)
    do j = 1,nrules
       read(fid,*) mask1(j), mask2(j), smoothing(j)
    end do

	call create_iscalar(grid,mask,CELL_EARTH)
	call read_iscalar(fid,mask)

!	allocate(temp(Nx),STAT=istat)
!	do
!	    read(fid,*,iostat=istat) k1, k2
!	    if (istat /= 0) exit
!	    do j = 1,Ny
!	    	read(fid,*,iostat=istat) temp
!	    	do k = k1,k2
!	    	  	mask%v(:,j,k) = temp
!	    	end do
!	    end do
!    end do

	close(fid)

    !open(unit=30,file='temp.cov',form='formatted',status='unknown')
    !call write_iscalar(30,mask)
    !close(30)

	deallocate(mask1,mask2,smoothing)
	call deall_iscalar(mask)

  end subroutine read_CmSqrt
