! *****************************************************************************
! * BOP
! * module: stag_grid
! Basic data structure to define numerical simulation on a staggered grid, like
! grid structure, staggered grid vectors (grid/ cell edges, grid/ cell faces;
! full storage and sparse), scalars (grid/cell centers, grid/ cell corners),
! and tangential boundary conditions.
! Belongs to SG_Basics class: staggered cartesian grid, data
! types defined on this grid, and operations defined on these data types. Not
! specific to EM problem, no dependency on outside (from other classes)
! modules other than from SG_Basics.
! * EOP

module GridDef

  ! all the modules being used are being listed explicitly (no
  ! inheritance)
  use math_constants
  implicit none

  ! Don't forget to overload the '=' sign: depending on the compiler, might
  ! run into trouble with the default assignment, since that sometimes doesn't
  ! copy allocatable or pointer arrays
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_grid
  END INTERFACE

  ! Possible grid types for EMfield, storing the intention of use for types
  ! such as cvector, cscalar, rvector, rscalar, sparsevecc.
  character(len=80), parameter		:: FACE = 'FACE'
  character(len=80), parameter		:: EDGE = 'EDGE'
  character(len=80), parameter		:: CENTER = 'CELL'
  character(len=80), parameter		:: CORNER = 'NODE'
  character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'

  ! ***************************************************************************
  ! * BOP
  ! type grid_param consists of parameters that define the basic grid geometry
  ! used for three dimensional numerical modeling
  type :: grid_t
  ! * EOP

     ! Grid coordinate system; important - used in EMfield
     character (len=80)			:: coords = Spherical

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is number of earth layers used in the grid modeling
     ! nzAir is number of air layers
     ! nz is grid dimension (number of cells) in the z-direction:
     ! nz = nzAir + nzEarth & nz = nzMantle + nzCrust
	 ! If no crust information is provided, nzCrust=nzAir
     integer               :: nx, ny, nz, nzEarth, nzMantle, nzCrust, nzAir

     ! Grid geometry:
	 ! No grid geometry is currently defined in the grid type definition
	 ! This may (and should be) added in the future, when Spherical and
	 ! Cartesian codes are united.

     ! Book-keeping on cumulative distances
	 ! No cumulative distances are currently defined in the grid type

 	 ! storing the (spherical) grid in Randie Mackie's format
	 ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
	 real(8), pointer, dimension(:)			:: x,y,z
	 ! storing the cell node coordinates, in radians (ph,th) and km (r)
	 ! dimensions are ph(nx+1), th(ny+1), r(nz+1)
	 ! the nx+1'st longitude ph(nx+1) = 360.0 * d2r (for interpolation)
	 real(8), pointer, dimension(:)  			:: ph,th,r

     ! allocated:  .true.  all the arrays have been allocated
     logical		                             :: allocated = .false.

  end type grid_t

  public         :: create_grid, deall_grid, copy_grid

Contains

     !************************************************************************
     subroutine create_grid(nx,ny,nz,grid)
       !  creates finite differences grid_t structure of
       !  size Nz x Ny, allocates arrays
       !
       implicit none
       integer, intent(in)		        :: nx,ny,nz
       type (grid_t), intent(inout)	    :: grid
       integer                          :: istat

       grid%nx = nx
       grid%ny = ny
       grid%nz = nz
       allocate(grid%x(nx), grid%y(ny+1), grid%z(nz+1), STAT=istat)
       allocate(grid%ph(nx+1), grid%th(ny+1), grid%r(nz+1), STAT=istat)
       grid%allocated = .true.

     end subroutine create_grid

     !************************************************************************
     subroutine deall_grid(grid)
       !  deallocates finite differences grid_t structure
       !
       implicit none
       type (grid_t), intent(inout)		:: grid
       integer                          :: istat

       if (associated(grid%x)) deallocate(grid%x, STAT=istat)
       if (associated(grid%y)) deallocate(grid%y, STAT=istat)
       if (associated(grid%z)) deallocate(grid%z, STAT=istat)
       if (associated(grid%ph)) deallocate(grid%ph, STAT=istat)
       if (associated(grid%th)) deallocate(grid%th, STAT=istat)
       if (associated(grid%r)) deallocate(grid%r, STAT=istat)
       grid%allocated = .false.

     end subroutine deall_grid

     !************************************************************************
     subroutine copy_grid(gridOut,gridIn)
       !  overloads the '=' sign
       !
       implicit none
       type (grid_t), intent(in)		:: gridIn
       type (grid_t), intent(out)		:: gridOut
       integer		        			:: nx,ny,nz

       nx = gridIn%nx
       ny = gridIn%ny
       nz = gridIn%nz

       call deall_grid(gridOut)
       call create_grid(nx,ny,nz,gridOut)

       gridOut%coords = gridIn%coords
       gridOut%nzEarth = gridIn%nzEarth
       gridOut%nzMantle = gridIn%nzMantle
       gridOut%nzCrust = gridIn%nzCrust
       gridOut%nzAir = gridIn%nzAir

       gridOut%x = gridIn%x
       gridOut%y = gridIn%y
       gridOut%z = gridIn%z
       gridOut%ph = gridIn%ph
       gridOut%th = gridIn%th
       gridOut%r = gridIn%r

       gridOut%allocated = gridIn%allocated

     end subroutine copy_grid

end module GridDef
