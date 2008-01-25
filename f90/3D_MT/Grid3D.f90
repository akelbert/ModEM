! *****************************************************************************
! initializes and does basic calculations for the grid. Computes basic derived
! grid parameters and are used repeatedly by other routines.
! Belongs to SG_Basics class: staggered cartesian grid, data
! types defined on this grid, and operations defined on these data types. Not
! specific to EM problem, no dependency on outside (from other classes) modules.
module grid3d

  use math_constants
  implicit none

  ! Initialization routines
  public                             	:: create_grid3D,deall_grid3D, &
                                GridCalcs, copy_grid3D

  ! ***************************************************************************
  ! type grid_param consists of parameters that define the basic grid geometry
  ! used for three dimensional numerical modeling
  type :: grid3d_t

     ! For possible future flexibility, store grid type as a character string
     ! Right now, the default is 'Cartesian Staggered'
     !  IS THIS EVER USED????   Makes no real sense.
     character (len=80)			:: gridType = 'Cartesian Staggered'

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is number of earth layers used in the grid modeling
     ! nzAir is number of air layers 
     ! nz is grid dimension (number of cells) in the z-direction: 
     ! nz = nzAir + nzEarth
     integer               :: nx, ny, nz, nzEarth, nzAir

     ! the origin of the model, by default set to zero
     real (kind=selectedPrec)			:: ox = 0.0, oy = 0.0, oz = 0.0
     !  the rotation angle in degrees, by defualt set to zero
     real (kind=selectedPrec)			:: rotdeg = 0.0

     ! Grid geometry:
     ! dx,delX are arrays of grid spacings in x-direction
     ! dx denotes spacings betgween cell edges: dimension: dx(nx)
     ! dxinv  = 1/ dx
     ! delX denotes spacings between cell centers: dimension: delX(nx+1)
     ! why dimensions: delX(nx+1) (delX(2) is distance between centers of cells
     ! 2 and 1)
     ! delXinv = 1/ delX
     ! dy,delY, dz, delZ are analagous arrays for other directions
     ! similarly, are dyinv and delYinv
     ! Note that the arrays are allocated dynamically
     real (kind=selectedPrec), pointer, dimension(:)        :: dx, dy, dz
     real (kind=selectedPrec), pointer, dimension(:)        :: dxinv, dyinv, dzinv
     real (kind=selectedPrec), pointer, dimension(:)        :: delX, delY, delZ
     real (kind=selectedPrec), pointer, dimension(:)        :: delXinv, delYinv, &
     delZinv

     ! Book-keeping on cumulative distances
     ! xEdge is the array for cumulative distance of the edge faces from the 
     ! coordinate axis with dimensions nx+1
     ! xCenter is the array for cumulative distance of the edge center from 
     ! the coordinate axis with dimensions nx
     ! yEdge, yCenter, zEdge, zCenter are analagous arrays for other directions
     ! Note that the arrays are allocated dynamically
     real (kind=selectedPrec), pointer, dimension(:)        :: xEdge, yEdge, zEdge
     real (kind=selectedPrec), pointer, dimension(:)        :: xCenter, yCenter, &
     zCenter
     
     ! total thickness of the air above 
     real (kind=selectedPrec)				     :: zAirThick

     ! allocated:  .true.  all the arrays have been allocated
     logical		                             :: allocated = .false.
     
  end type grid3d_t
  
Contains

  !************************************************************************
  subroutine create_Grid3D(Nx,Ny,NzAir,NzEarth,grid)
    !  creates finite differences grid3d_t structure of
    !  size  Nx x Ny Nz, allocates arrays
    !
    implicit none
    integer, intent(in)				:: Nx,Ny,NzAir,NzEarth
    type (grid3d_t) , intent(inout)		:: grid
    
    !local variables
    integer					:: Nz

    Nz = NzEarth+NzAir
    grid%NzAir = NzAir
    grid%Nx = Nx
    grid%Ny = Ny
    grid%NzEarth = NzEarth
    grid%Nz = Nz
    allocate(grid%Dx(Nx))
    allocate(grid%Dy(Ny))
    allocate(grid%Dz(Nz))

    ! dxinv  = 1/ dx and similarly for dyinv and dzinv
    allocate(grid%dxinv(Nx))
    allocate(grid%dyinv(Ny))
    allocate(grid%dzinv(Nz))

    ! delX, delY, and delZ are the distances between the electrical field
    ! defined on the center of the edges in x, y, and z axes, respectively.
    allocate(grid%delX(Nx+1))
    allocate(grid%delY(Ny+1))
    allocate(grid%delZ(Nz+1))

    ! delXinv = 1/ delX and similarly for delYinv and delZinv
    allocate(grid%delXinv(Nx+1))
    allocate(grid%delYinv(Ny+1))
    allocate(grid%delZinv(Nz+1))

    ! xEdge is the array for cumulative distance of the edge for each
    ! grid (starting from the coordinate axes) with dimensions nx+1
    ! xCenter is the array for cumulative distance of the center for each
    ! grid (starting from the coordinate axes) with dimension n
    ! yEdge, yCenter, zEdge, zCenter are analagous arrays for other directions
    allocate(grid%xEdge(Nx+1))
    allocate(grid%yEdge(Ny+1))
    allocate(grid%zEdge(Nz+1))
    allocate(grid%xCenter(Nx))
    allocate(grid%yCenter(Ny))
    allocate(grid%zCenter(Nz))

    grid%allocated = .true.

  end subroutine create_Grid3D

  ! **************************************************************************
  subroutine copy_grid3D(gridOut,gridIn)

  !  copies gridIn to gridOut; cannot overwrite, of course!

  type(grid3d_t),intent(in)		:: gridIn
  type(grid3d_t),intent(inout)		:: gridOut

     if(gridOut%allocated) then
        !  just deallocate, and start over cleanly
        call deall_grid3D(gridOut) 
     endif

     call create_Grid3D(gridIn%Nx,gridIn%Ny,gridIn%NzAir, &
             gridIn%NzEarth,gridOut)

     gridOut%Dz = gridIn%Dz
     gridOut%Dy = gridIn%Dy
     gridOut%Dx = gridIn%Dx

     call gridCalcs(gridOut)

  end subroutine copy_grid3D

  ! **************************************************************************
  subroutine deall_grid3D(grid)

    type (grid3d_t) , intent(inout)	:: grid

    deallocate(grid%Dz)
    deallocate(grid%Dy)
    deallocate(grid%Dz)

    deallocate(grid%dxinv)
    deallocate(grid%dyinv)
    deallocate(grid%dzinv)

    deallocate(grid%delX)
    deallocate(grid%delY)
    deallocate(grid%delZ)

    deallocate(grid%delXinv)
    deallocate(grid%delYinv)
    deallocate(grid%delZinv)

    deallocate(grid%xEdge)
    deallocate(grid%yEdge)
    deallocate(grid%zEdge)
    deallocate(grid%xCenter)
    deallocate(grid%yCenter)
    deallocate(grid%zCenter)

    grid%allocated = .false.
    grid%NzAir = 0
    grid%Nx = 0
    grid%Ny = 0
    grid%NzEarth = 0
    grid%Nz = 0

 end subroutine deall_grid3D

  ! **************************************************************************
  ! * GridCalcs does calculations for grid geometry, which cannot be done
  !   until dx, dy, dz, and the origin are set.
  !   Normal usage is to first call create_Grid3D to set grid dimensions
  !    and allocate arrays, read dx, dy, dz and set these elements of the grid,
  !    then call GridCalcs to do all other computations.  By including the optional
  !    origin argument, grid centers and edges are given in absolute coordinates
  !    (i.e., the origin of the grid at the Earth surface is set to the origin,
  !      and variables like xCenter, yEdge, etc. are given in the same coordinate 
  !      system).  If argument origin is not present, whatever is set already in the grid
  !      origin is used; by default this is initialized to zero.
  subroutine gridCalcs(grid, origin)

    implicit none
    type(grid3d_t), target, intent(inout)    :: grid
    real(kind=selectedPrec), intent(in), optional	  :: origin(3)
              
    integer                               :: ix,iy,iz
    integer                               :: status 
    real (kind=selectedPrec)                         :: xCum, yCum, zCum

    grid%dxinv = 1/ grid%dx
    grid%dyinv = 1/ grid%dy
    grid%dzinv = 1/ grid%dz

    grid%rotdeg = grid%rotdeg
    if (present(origin)) then
    	grid%ox = origin(1)
    	grid%oy = origin(2)
    	grid%oz = origin(3)
    end if
    
    grid%xEdge(1) = grid%ox
    grid%yEdge(1) = grid%oy
    grid%zEdge(1) = grid%oz
    xCum = 0.0
    yCum = 0.0
    zCum = 0.0
    do ix = 1, grid%nx
       xCum = xCum + grid%dx(ix)
       grid%xEdge(ix+1) = xCum + grid%ox
    enddo
    do iy = 1, grid%ny
       yCum = yCum + grid%dy(iy)
       grid%yEdge(iy+1) = yCum + grid%oy
    enddo
    !  NOTE: adjust for origin later to get airthickness, reference to origin
    !    at Earth's surface correct!
    do iz = 1, grid%nz
       zCum = zCum + grid%dz(iz)
       grid%zEdge(iz+1) = zCum
    enddo
    grid%zAirThick = grid%zEdge(grid%nzAir+1)

    ! distance between center of the grids
    grid%delX(1) = grid%dx(1)      
    DO ix = 2,grid%nx
       grid%delX(ix) = grid%dx(ix-1) + grid%dx(ix)
    ENDDO
    grid%delX(grid%nx+1) = grid%dx(grid%nx)
    grid%delX = grid%delX/2.0
    
    grid%delY(1)    = grid%dy(1)
    DO iy = 2,grid%ny
       grid%delY(iy) = grid%dy(iy-1) + grid%dy(iy)
    ENDDO
    grid%delY(grid%ny+1) = grid%dy(grid%ny)
    grid%delY = grid%delY/2.0
    
    grid%delZ(1)    = grid%dz(1)
    DO iz = 2,grid%nz
       grid%delZ(iz) = grid%dz(iz-1) + grid%dz(iz)
    ENDDO
    grid%delZ(grid%nz+1) = grid%dz(grid%nz)
    grid%delZ = grid%delZ/ 2.0

    grid%delXinv = 1/ grid%delX
    grid%delYinv = 1/ grid%delY
    grid%delZinv = 1/ grid%delZ

    xCum = 0.0
    yCum = 0.0
    zCum = 0.0
    ! cumulative distance between the centers, adjusted to model origin
    do ix = 1, grid%nx
       xCum = xCum + grid%delX(ix)
       grid%xCenter(ix) = xCum + grid%ox 
    enddo
    do iy = 1, grid%ny
       yCum = yCum + grid%delY(iy)
       grid%yCenter(iy) = yCum + grid%oy
    enddo
    do iz = 1, grid%nz
       zCum = zCum + grid%delZ(iz)
       grid%zCenter(iz) = zCum
    enddo
           
    !  need to be careful here ... grid origin is given at Earth's surface,
    !   not top of model domain!
    do iz = 1, grid%nz
       grid%zCenter(iz) = grid%zCenter(iz)-grid%zAirThick+grid%oz
       grid%zEdge(iz) = grid%zEdge(iz)-grid%zAirThick+grid%oz
    enddo
    grid%zEdge(grid%nz+1) = grid%zEdge(grid%nz+1)-grid%zAirThick+grid%oz
    

  end subroutine GridCalcs

end module grid3d
