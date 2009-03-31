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

     ! the origin is the lower-left corner at the top of the grid model
     ! the origin of the model, by default set to zero
     real (kind=prec)				:: ox = 0.0, oy = 0.0, oz = 0.0
     !  the rotation angle in degrees, by defualt set to zero
     real (kind=prec)				:: rotdeg = 0.0

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


     ! total thickness of the air above
     real (kind=prec)				     :: zAirThick = R_ZERO

     ! allocated:  .true.  all the arrays have been allocated
     logical		                             :: allocated = .false.

  end type grid_t


  ! ***************************************************************************
  ! * BOP
  ! grid_def extracts basic information from an input file and passes it on
  ! for creating grid. This is the basic information that has to be provided by
  ! an input file (any) to the forward modeling program.
  type :: grid_def
  ! * EOP

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is number of earth layers used in the grid modeling
     ! nzAir is number of air layers
     ! nz is grid dimension (number of cells) in the z-direction:
     ! nz = nzAir + nzEarth
     integer               			:: nx, ny, nz, nzEarth, nzAir

     ! the origin is the lower-left corner at the top of the grid model
     ! the origin of the model, by default set to zero
     real (kind=prec)				:: ox = 0.0, oy = 0.0, oz = 0.0
     !  the rotation angle in degrees, by defualt set to zero
     real (kind=prec)				:: rotdeg = 0.0

     ! dx, dy, and dz are arrays of grid spacing in x-, y-, and z-direction,
     ! respectively
     real (kind=prec), pointer, dimension(:)	:: dx, dy, dz

     ! allocated:  .true.  all the arrays have been allocated
     logical		                             :: allocated = .false.

  end type grid_def


end module GridDef
