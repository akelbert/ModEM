!**********************************************************************
! WScov - CmSqrt by Weerachai Siripunvaraporn

  subroutine WScov(m)

  type (modelParam_t), intent(inout) :: m

  call setup1DCM(m%grid)
  call solveDiff(m)

  end subroutine WScov

  !------------------------------------------------------------------------------

  subroutine setup1DCM(gd)

  type (grid2d_t), intent(inout) :: gd

  real (kind=selectedPrec), pointer, dimension(:,:)   ::  YDF
  real (kind=selectedPrec), pointer, dimension(:,:)   ::  ZDF

  integer  :: iy, iz, NzEarth ,Ny, istat
  real (kind=selectedPrec) :: ylen, zlen

  Ny      = gd%Ny
  NzEarth = gd%Nz - gd%Nza
  allocate(CmSqrt%YDIF(2,gd%Ny,NzEarth))
  allocate(CmSqrt%ZDIF(2,gd%Ny,NzEarth))
  CmSqrt%ylen = 0.05
  CmSqrt%zlen = 0.05

  ylen = CmSqrt%ylen
  zlen = CmSqrt%zlen

  ! horizontal diffusion
  do iz = 1,NzEarth
    call setupYCM(Ny,ylen,YDF)
    do iy = 1,Ny
      CmSqrt%YDIF(1,iy,iz) = YDF(1,iy)
      CmSqrt%YDIF(2,iy,iz) = YDF(2,iy)
    enddo
  enddo

  ! vertical diffusion
  do iy = 1,Ny
    call setupZCM(NzEarth,zlen,ZDF)
    do iz = 1,NzEarth
      CmSqrt%ZDIF(1,iy,iz) = ZDF(1,iz)
      CmSqrt%ZDIF(2,iy,iz) = ZDF(2,iz)
    enddo
  enddo

  deallocate(YDF,STAT=istat)
  deallocate(ZDF,STAT=istat)

  end subroutine setup1DCM

!------------------------------------------------------------------------------

  subroutine setupYCM(Ny,ylen,YDF)

  !use wsLAPACK
  real (kind=selectedPrec), pointer, dimension(:,:)   ::  YDF
  integer  :: Ny
  real (kind=selectedPrec) :: ylen

  integer  :: iy,info
  real (kind=selectedPrec) :: cr, cl
  real (kind=selectedPrec), allocatable, dimension(:)   ::  dfyc, dfyr


  allocate(YDF(2,Ny))
  allocate(dfyc(Ny))
  allocate(dfyr(Ny))
  do iy = 1,Ny
    dfyc(iy) = 0.
    dfyr(iy) = 0.
  enddo

  do iy = 1,Ny
    cr = ylen
    cl = ylen

    dfyc(iy) = 1. + cr + cl
    dfyr(iy) = -cr
!   dfyc(iy) = 1. - cr - cl
!   dfyr(iy) = cr
  enddo

  do iy =1,Ny
    YDF(1,iy) = 0.
    YDF(2,iy) = 0.
  enddo

  do iy =2,Ny
    YDF(1,iy) = dfyr(iy-1)
  enddo
  do iy =1,Ny
    YDF(2,iy) = dfyc(iy)
  enddo

  call dpbtrf('U',Ny,1,YDF,2,info)

  deallocate(dfyc)
  deallocate(dfyr)

  end subroutine setupYCM

!------------------------------------------------------------------------------

  subroutine setupZCM(Nzb,zlen,ZDF)

  !use wsLAPACK
  real (kind=selectedPrec), pointer,  dimension(:,:)   ::  ZDF
  integer  :: Nzb
  real (kind=selectedPrec) :: zlen

  integer  :: iz,info
  real (kind=selectedPrec) :: cr, cl
  real (kind=selectedPrec), allocatable, dimension(:)   ::  dfzc, dfzr


  allocate(ZDF(2,Nzb))
  allocate(dfzc(Nzb))
  allocate(dfzr(Nzb))
  do iz = 1,Nzb
    dfzc(iz) = 0.
    dfzr(iz) = 0.
  enddo

  do iz = 1,Nzb
    cr = zlen
    cl = zlen

    dfzc(iz) = 1. + cr + cl
    dfzr(iz) = -cr
!   dfzc(iz) = 1. - cr - cl
!   dfzr(iz) =  cr
  enddo

  do iz =1,Nzb
    ZDF(1,iz) = 0.
    ZDF(2,iz) = 0.
  enddo

  do iz =2,Nzb
    ZDF(1,iz) = dfzr(iz-1)
  enddo
  do iz =1,Nzb
    ZDF(2,iz) = dfzc(iz)
  enddo

  call dpbtrf('U',Nzb,1,ZDF,2,info)

  deallocate(dfzc)
  deallocate(dfzr)

  end subroutine setupZCM

!------------------------------------------------------------------------------

  subroutine solveDiff(CmJ)

  type (modelParam_t), intent(inout) :: CmJ

  integer            :: Ny,Nz,it, dTime

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  dTime = 10

! do iz = 1,Nz
!   do iy = 1,Ny
!     write(45,*) CmJ%v(iy,iz)
!   enddo
! enddo

  do it = 1,dTime
    call yDiff(CmJ)
    call zDiff(CmJ)
  enddo
  do it = 1,dTime
    call zDiff(CmJ)
    call yDiff(CmJ)
  enddo

! do iz = 1,Nz
!   do iy = 1,Ny
!     write(46,*) CmJ%v(iy,iz)
!   enddo
! enddo

  end subroutine solveDiff

!------------------------------------------------------------------------------

  subroutine yDiff(CmJ)

  !use wsLAPACK
  type (modelParam_t), intent(inout) :: CmJ

  real (kind=selectedPrec), allocatable, dimension(:,:)   ::  YDF
  real (kind=selectedPrec), allocatable, dimension(:)     ::  uy

  integer            :: Ny,Nz,iy,iz,info,istat

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  allocate(YDF(2,Ny))
  allocate(uy(Ny))

  do iz = 1,Nz
    do iy = 1,Ny
      uy(iy)    = CmJ%v(iy,iz)
      YDF(1,iy) = CmSqrt%YDIF(1,iy,iz)
      YDF(2,iy) = CmSqrt%YDIF(2,iy,iz)
    enddo

    call dpbtrs('U',Ny,1,1,YDF,2,uy,Ny,info)

    do iy = 1,Ny
      CmJ%v(iy,iz) = uy(iy)
    enddo

  enddo

  deallocate(YDF,STAT=istat)
  deallocate(uy,STAT=istat)

  end subroutine yDiff

!------------------------------------------------------------------------------

  subroutine zDiff(CmJ)

  !use wsLAPACK
  type (modelParam_t), intent(inout) :: CmJ

  real (kind=selectedPrec), allocatable, dimension(:,:)   ::  ZDF
  real (kind=selectedPrec), allocatable, dimension(:)     ::  uz


  integer            :: Ny,Nz,iy,iz,info,istat

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  allocate(ZDF(2,Nz))
  allocate(uz(Nz))

  do iy = 1,Ny
    do iz = 1,Nz
      uz(iz)    = CmJ%v(iy,iz)
      ZDF(1,iz) = CmSqrt%ZDIF(1,iy,iz)
      ZDF(2,iz) = CmSqrt%ZDIF(2,iy,iz)
    enddo

    call dpbtrs('U',Nz,1,1,ZDF,2,uz,Nz,info)

    do iz = 1,Nz
      CmJ%v(iy,iz) = uz(iz)
    enddo

  enddo

  deallocate(ZDF,STAT=istat)
  deallocate(uz,STAT=istat)

  end subroutine zDiff
