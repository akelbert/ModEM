!-------------------------------------------------------------------!
! 3D_MT model covariance by Weerachai Siripunvaraporn: procedures.  !
!                                                                   !
! Converted to Fortran 90 and adapted to the modular structure by   !
! Anna Kelbert  (March 24, 2008).                                   !
!                                                                   !
! Weerachai uses the prior model (integer format) to switch off     !
! the model covariance for frozen regions in Setup1DCM.             !
! Judging from the parameters in the Setup subroutines, he has      !
! replaced the explicit time forwards scheme with the implicit time !
! backwards scheme (see solution of the diffusion equation and      !
! Crank-Nicholson scheme in the Numerical Recipes for a discussion).!
!-------------------------------------------------------------------!

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
    call solveDiff(CmSqrt,m2)

  end subroutine multBy_CmSqrt

! *******************************************************************
  subroutine create_CmSqrt(m,cfile)

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile

    if (present(cfile)) then
      ! do nothing: a placeholder for other covariance operators
    end if

    call setup1DCM(m%grid,CmSqrt)

    CmSqrt%allocated = .true.

  end subroutine create_CmSqrt

! *******************************************************************
  subroutine deall_CmSqrt()

    integer                             :: istat

    deallocate(CmSqrt%XDIF, STAT = istat)
    deallocate(CmSqrt%YDIF, STAT = istat)
    deallocate(CmSqrt%ZDIF, STAT = istat)

  end subroutine deall_CmSqrt

! *******************************************************************
!
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.
!
!----------------------------------------------------------------------
!
!     SOLVING DIFFUSION EQU. IN EXPLICIT WAY.
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

! *****************************************************************************
! Model covariance operators;
! currently only contains the WS model covariance,
!**********************************************************************
! Note that only the driver subroutine requires knowledge of modelParam.
! Without it, we would only need to include sg_scalar & grid3d.

  subroutine WScov(m)

	type (modelParam_t), intent(inout) :: m
	type (mCov_t)                      :: Cm
	type (rscalar)                     :: logcond

	call getValue_modelParam(m,LOGE,logcond)

	call setup1DCM(Cm,logcond%grid)

	dtime = 10
	call solveDiff(dtime,logcond)

	call setValue_modelParam(m,LOGE,logcond)

  end subroutine WScov

!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SOLVING DIFFUSION EQU. IN EXPLICIT WAY.
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.
!...  Adapted from Setup1DCM_NEW by Anna Kelbert

  SUBROUTINE Setup1DCM(gd,Cm)

     integer, dimension(:,:,:), intent(in)   :: PMod
     type (grid3d_t), intent(in)             :: gd
     type (mCov_t), intent(out)              :: Cm


     real (kind=prec), pointer, dimension(:,:)   ::  XDF
     real (kind=prec), pointer, dimension(:,:)   ::  YDF
     real (kind=prec), pointer, dimension(:,:)   ::  ZDF

     integer                  :: jj, ix, iy, iz, Nx, Ny, NzEarth
     real (kind=prec) :: xlen, ylen, zlen

	  Nx      = gd%Nx
	  Ny      = gd%Ny
	  NzEarth = gd%Nz - gd%Nza

	  allocate(Cm%YDIF(2,Nx,Ny,NzEarth))
	  allocate(Cm%YDIF(2,Nx,Ny,NzEarth))
	  allocate(Cm%ZDIF(2,Nx,Ny,NzEarth))
	  Cm%xlen = 0.05
	  Cm%ylen = 0.05
	  Cm%zlen = 0.05

	  XLen = Cm%xlen
	  YLen = Cm%ylen
	  ZLen = Cm%zlen

!>>>  X : Horizontal Diffusion
      DO iz = 1,NzEarth
        DO iy = 1,Ny
          CALL SetupXCM(Nx,XLen,XDF)
          DO ix = 1,Nx
            XDIF(1,ix,iy,iz) = XDF(1,ix)
            XDIF(2,ix,iy,iz) = XDF(2,ix)
          ENDDO ! ix
        ENDDO ! iy
      ENDDO ! iz


!>>>  Y : Horizontal Diffusion
      DO iz = 1,NzEarth
        DO ix = 1,Nx
          CALL SetupYCM(Ny,YLen,YDF)
          DO iy = 1,Ny
            YDIF(1,ix,iy,iz) = YDF(1,iy)
            YDIF(2,ix,iy,iz) = YDF(2,iy)
          ENDDO ! iy
        ENDDO ! ix
      ENDDO ! iz


!>>>  Vertical Diffusion
      DO ix = 1,Nx
        DO iy = 1,Ny
          CALL SetupZCM(Nzb,ZLen,ZDF)
          DO iz = 1,NzEarth
            ZDIF(1,ix,iy,iz) = ZDF(1,iz)
            ZDIF(2,ix,iy,iz) = ZDF(2,iz)
          ENDDO ! iy
        ENDDO ! iy
      ENDDO ! ix


  END SUBROUTINE Setup1DCM ! Setup1DCM

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!>>>  Horizontal Diffusion : X
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SetupXCM(Nx,xlen,XDF)

      integer Nx,info
      real (kind=prec)  XDF(2,Nx),xlen
      real (kind=prec)  cr,cl,difxc(Nx),difxr(Nx)

      integer ix

      DO ix = 1,Nx
        cr = xlen
        cl = xlen

        difxc(ix) = ONE + cr + cl
        difxr(ix) = -cr
      ENDDO ! ix

      CALL ConstantMatrixR8(XDF,2,Nx,2,Nx,R_ZERO)

      DO ix = 2,Nx
        XDF(1,ix) = difxr(ix-1)
      ENDDO
      DO ix = 1,Nx
        XDF(2,ix) = difxc(ix)
      ENDDO

      CALL DPBTRF('U',Nx,1,XDF,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : XDiff',info
        WRITE(12,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(12,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : XDiff',info
        STOP
      ENDIF

      RETURN

  END SUBROUTINE SetupXCM ! SetupXCM


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!>>>  Horizontal Diffusion : Y
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SetupYCM(Ny,ylen,YDF)

      integer Ny,info,ygd(*)
      real (kind=prec)  YDF(2,Ny),ylen
      real (kind=prec)  cr,cl,difyc(Ny),difyr(Ny)
      integer iy

      DO iy = 1,Ny
        cr = ylen
        cl = ylen

        difyc(iy) = ONE + cr + cl
        difyr(iy) = -cr
      ENDDO ! iy

      CALL ConstantMatrixR8(YDF,2,Ny,2,Ny,R_ZERO)

      DO iy = 2,Ny
        YDF(1,iy) = difyr(iy-1)
      ENDDO
      DO iy = 1,Ny
        YDF(2,iy) = difyc(iy)
      ENDDO

      CALL DPBTRF('U',Ny,1,YDF,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : YDiff',info
        WRITE(12,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(12,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : YDiff',info
        STOP
      ENDIF

      RETURN

  END SUBROUTINE SetupYCM ! SetupYCM


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!>>>  Vertical Diffusion : Z
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SetupZCM(Nz,zlen,ZDF)

      integer Nz,info,zgd(*)
      real (kind=prec)  ZDF(2,Nz),zlen
      real (kind=prec)  cr,cl,difzc(Nz),difzr(Nz)
      integer iz

      DO iz = 1,Nz
        cr = zlen
        cl = zlen

        difzc(iz) = ONE + cr + cl
        difzr(iz) = -cr
      ENDDO ! iz

      CALL ConstantMatrixR8(ZDF,2,Nz,2,Nz,R_ZERO)

      DO iz = 2,Nz
        ZDF(1,iz) = difzr(iz-1)
      ENDDO
      DO iz = 1,Nz
        ZDF(2,iz) = difzc(iz)
      ENDDO
!     DO iz = 1,Nz
!       write(6,*) 'ZDF BEFORe ',iz,ZDF(1,iz),ZDF(2,iz)
!     ENDDO

      CALL DPBTRF('U',Nz,1,ZDF,2,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
        WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : ZDiff',info
        STOP
      ENDIF

!     DO iz = 1,Nz
!       write(6,*) 'ZDF AFTER ',iz,ZDF(1,iz),ZDF(2,iz)
!     ENDDO
!     stop

      RETURN

  END SUBROUTINE SetupZCM ! SetupZCM


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE SolveDiff(DTime,Cm,u2)

      integer, intent(in)                     :: DTime
      type (mCov_t), intent(in)               :: Cm
      type (rscalar), intent(inout)           :: u2


      integer  idt

      DO idt = 1,DTime
        CALL X_Diff(Cm,u2)
        CALL Y_Diff(Cm,u2)
        CALL Z_Diff(Cm,u2)
      ENDDO ! idt

      DO idt = 1,DTime
        CALL Z_Diff(Cm,u2)
        CALL Y_Diff(Cm,u2)
        CALL X_Diff(Cm,u2)
      ENDDO ! idt

  END SUBROUTINE SolveDiff ! SolveDiff


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE X_Diff(Cm,u2)

	  type (mCov_t), intent(in)      :: Cm
	  type (rscalar), intent(inout)  :: u2


	  real (kind=prec), allocatable, dimension(:,:)   ::  XDF
	  real (kind=prec), allocatable, dimension(:)     ::  ux

	  integer            :: Nx,Ny,NzEarth,ix,iy,iz,info

      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz
      allocate(XDF(2,Nx))
      allocate(ux(Nx))

      DO iz = 1,NzEarth
        DO iy = 1,Ny
          DO ix = 1,Nx
            ux(ix) = u2%v(ix,iy,iz)
            XDF(1,ix) = Cm%XDIF(1,ix,iy,iz)
            XDF(2,ix) = Cm%XDIF(2,ix,iy,iz)
          ENDDO ! ix

          CALL DPBTRS('U',Nx,1,1,XDF,2,ux,Nx,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : XDiff',info
            STOP
          ENDIF

          DO ix = 1,Nx
            u2%v(ix,iy,iz) = ux(ix)
          ENDDO
        ENDDO ! iy
      ENDDO ! iz

      deallocate(XDF,ux)

  END SUBROUTINE X_DIFF ! X_DIFF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE Y_Diff(Cm,u2)

	  type (mCov_t), intent(in)      :: Cm
	  type (rscalar), intent(inout)  :: u2

	  real (kind=prec), allocatable, dimension(:,:)   ::  YDF
	  real (kind=prec), allocatable, dimension(:)     ::  uy

	  integer            :: Nx,Ny,NzEarth,ix,iy,iz,info

      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz
      allocate(YDF(2,Ny))
      allocate(uy(Ny))

      DO iz = 1,NzEarth
        DO ix = 1,Nx
          DO iy = 1,Ny
            uy(iy) = u2%v(ix,iy,iz)
            YDF(1,iy) = Cm%YDIF(1,ix,iy,iz)
            YDF(2,iy) = Cm%YDIF(2,ix,iy,iz)
          ENDDO ! iy

          CALL DPBTRS('U',Ny,1,1,YDF,2,uy,Ny,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : YDiff',info
            STOP
          ENDIF

          DO iy = 1,Ny
            u2%v(ix,iy,iz) = uy(iy)
          ENDDO ! iy
        ENDDO ! ix
      ENDDO ! iz

      deallocate(YDF,uy)

  END SUBROUTINE Y_Diff ! Y_DIFF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

  SUBROUTINE Z_Diff(Cm,u2)

	  type (mCov_t), intent(in)      :: Cm
	  type (rscalar), intent(inout)  :: u2

	  real (kind=prec), allocatable, dimension(:,:)   ::  ZDF
	  real (kind=prec), allocatable, dimension(:)     ::  uz

	  integer            :: Nx,Ny,NzEarth,ix,iy,iz,info


      Nx = u2%nx
      Ny = u2%ny
      NzEarth = u2%nz
      allocate(ZDF(2,NzEarth))
      allocate(uz(NzEarth))

      DO ix = 1,Nx
        DO iy = 1,Ny
          DO iz = 1,Nzb
            uz(iz) = u2%v(ix,iy,iz)
            ZDF(1,iz) = Cm%ZDIF(1,ix,iy,iz)
            ZDF(2,iz) = Cm%ZDIF(2,ix,iy,iz)
          ENDDO ! iy

          CALL DPBTRS('U',Nzb,1,1,ZDF,2,uz,ZMAX,info)
          IF (info.NE.0) THEN
            WRITE(6,*) '!!! ATTENTION, ERROR IN DECOMPOSING MATRIX'
            WRITE(6,*) 'ERROR DECOMPOSE DIFFUSION MATRIX : ZDiff',info
            STOP
          ENDIF

          DO iz = 1,Nzb
            u2%v(ix,iy,iz) = uz(iz)
          ENDDO ! iz
        ENDDO ! iy
      ENDDO ! ix

      deallocate(ZDF,uz)

  END SUBROUTINE Z_Diff ! Z_DIFF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.

end module modelCov
