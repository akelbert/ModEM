! *****************************************************************************
module modeldef
	! This module contains the type definitions for the structures that contain
	! the full information about the model parametrization

  implicit none

  ! ***************************************************************************
  ! * storing the forward solver control parameters
  type :: relaxation

      integer                                   :: ipotloopmax
	  integer									:: nrelmax
      integer                                   :: n_reldivh
      integer                                   :: ipot0
	  integer									:: ipotint
      integer                                   :: ipot_max
      real(8)                                   :: errend

  end type relaxation

  ! ***************************************************************************
  ! * type point_info contains information about a point in the Earth: coords
  ! * and relative location on the grid
  type :: point_info

	  real(8)								  :: phi,theta,r

  end type point_info


  ! ***************************************************************************
  ! * type shell_info contains the thin shell S-conductance on the uppermost
  ! * few km of the solid Earth - mainly Earth crust. A typical depth is 12650 m,
  ! * the values correspond to the given grid in GM coordinates so that
  ! * rho(i,j) = 1/(cond(i,j)/depth) = depth/cond(i,j) (where depth in in metres)
  type :: shell_info

	  real(8), dimension(:,:), pointer		  :: cond !(nx,ny) - conductance

  end type shell_info	! shell_info


  ! ***************************************************************************
  ! * type func_info contains the information about a single function that is
  ! * multiplied by the respective coefficient in the model parametrization;
  ! * add more information entries if necessary; for this, we may want to use
  ! * an analogue of "union" structure available in C, once Fortran has it. 
  type :: func_info

	  ! order number in functionals structure
	  integer								  :: num

	  ! which function is used?
	  character(80)							  :: name

	  ! degree and order, if F is a spherical harmonic term
	  integer								  :: l,m

	  ! range of values, if F is an identity function
	  ! real(8)								  :: xmin,xmax,ymin,ymax

  end type func_info	! func_info


  ! ***************************************************************************
  ! * type coeff_info contains the information about a single coefficient of the
  ! * model parametrization
  type :: coeff_info

     type (layer_info), pointer			  :: L	! which layer
     type (func_info) , pointer			  :: F	! which function it corresponds to
     type (coeff_info), pointer			  :: prev	! same coefficient prev layer
     type (coeff_info), pointer			  :: next	! same coefficient next layer
	 character(80)				  :: code ! information to identify the coefficient
     real(8)					  :: value	! value of coeff.
     real(8)                      :: min,max ! range information for inversion
	 ! if coefficient is fixed/frozen, no derivative will be provided in the output
     logical					  :: frozen	
	 ! it is possible to create a padding zero-valued coefficient that does not exist
	 logical					  :: exists 

  end type coeff_info	! coeff_info


  ! ***************************************************************************
  ! * type layer_info contains the information for a single layer in the layer
  ! * structure of the model parametrization, that is defined in modeldef.f90;
  ! * going from the first (uppermost) layer to the last that ends at CMB
  ! * Air layer is layer 0.
  type :: layer_info

	  ! order number in layer structure
	  integer											  :: num

	  ! upper layer boundary and the width of the layer in kilometers
	  real(8)											  :: lbound,ubound
	  real(8)											  :: width,depth 

	  ! layer regularization parameters
	  real(8)											  :: alpha,beta

	  ! the indicator of whether the linear combination of F_i gives log_10(rho)
	  logical											  :: if_log 

	  ! the filename of additional rho info if it exists; if so, add these
	  ! values to either rho or log_10(rho) (as specified by if_log) on the grid 
	  character(80)										  :: fn=''

  end type layer_info	! layer_info



  ! ***************************************************************************
  ! * type param_info combines together all information about the parametrization
  ! * in the most general case. When mapped onto grid, this information provides
  ! * a resistivity structure. It is also used to compute the gradients.
  ! * Since the allocatable attribute inside a derived data type definition
  ! * is not universally supported, we use pointers instead.

  type :: param_info

	  ! all information about the layer structure
	  integer											  :: nF
	  type (func_info),  dimension(:), pointer			  :: F 

	  ! all information about the layer structure
	  integer											  :: nL
	  type (layer_info), dimension(:), pointer			  :: L 

	  ! all coefficient and function information
	  !integer											  :: np
	  !type (coeff_info), dimension(:), pointer			  :: p

	  ! all coefficient and function information
	  integer											  :: nc ! number of active coeffs
	  type (coeff_info), dimension(:,:), pointer		  :: c	! structure of all coeffs

	  logical											  :: allocated

  end type param_info	! param_info


  ! ***************************************************************************
  ! * type misfit_def contains information that specifies how to calculate the
  ! * misfit; normally this is user-defined.
  type :: misfit_def

	character(80)					:: name	  ! type of data functional
	!real(8)							:: alpha  ! horizontal smoothing
	!real(8)							:: beta	  ! vertical smoothing
	real(8)							:: mu  ! damping parameter

  end type misfit_def


  ! ***************************************************************************
  ! * type misfit_info contains information about misfit and derivative
  ! * for a single frequency and functional type; normally we would create
  ! * an array type (misfit_info) misfit(nfreq,nfunc). In future this will
  ! * replace the current way of storing misfit information.
  type :: misfit_info

	real(8)							:: w  ! weight
	real(8)							:: v  ! non-averaged misfit value
	integer							:: n  ! number of data values (to average)
	type (param_info)				:: d  ! derivative

  end type misfit_info


  ! ***************************************************************************
  ! * type int_array contains an array of integer indices

  type :: int_array

	  integer, dimension(:), pointer					  :: num

  end type int_array


  ! ***************************************************************************
  ! * type real_array contains an array of real(8) indices

  type :: real_array

	  real(8), dimension(:), pointer					  :: value

  end type real_array


  ! ***************************************************************************
  ! * type grid_t contains the necessary information about user-defined grid
  ! This is now part of SG_Basics/stag_grid.f90
!  type :: grid_t
!
!	  ! storing the grid dimensions
!	  integer               					:: nx, ny, nz, nzEarth, nzAir
!	  ! storing the (spherical) grid in Randie Mackie's format
!	  ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
!	  real(8), pointer, dimension(:)			:: x,y,z
!	  ! storing the cell node coordinates, in radians (ph,th) and km (r)
!	  ! dimensions are ph(nx+1), th(ny+1), r(nz+1)
!	  real(8), pointer, dimension(:)  			:: ph,th,r
!
!  end type grid_t	! grid_t



  ! ***************************************************************************
  ! * type cell_info contains the values of phi, theta in radians and radius
  ! * in kilometers for each cell centre on the user-defined grid. Will be used
  ! * as a 3-D array of dimensions (nx,ny,nzEarth)
  type :: cell_info

	  integer									:: i,j,k
	  ! the co-ordinates of the cell centre, in radians (ph,th) and km (r)
	  real(8)									:: ph,th,r

  end type cell_info	! cell_info


Contains


  ! ***************************************************************************
  ! * These subroutines will be used as basic subroutines for various purposes


  ! Evaluate order number of Ylm in a layer

  integer function iY(l,m)

	implicit none

	integer, intent(in)	:: l,m
	integer				:: k

	iY=0
	if ((l<0).or.(m>l)) then
	  write(0,*) 'Error: (iY) wrong usage of spherical harmonics'
	  stop		
	else if (l==0) then
	  iY=1
	else
	  do k=0,l-1
		iY=iY+(2*k+1)
	  end do
	  if (m==0) then
		iY=iY+1
	  else if (m<0) then
		iY=iY+2*iabs(m)+1
	  else if (m>0) then
		iY=iY+2*m
	  end if
	end if

  end function iY


  ! Compute the value of the function at given phi and theta values

  real(8) function F_at_point(f,pt)

	use paramfunc
	implicit none

	type (func_info), intent(in)		  :: f
	type (point_info),intent(in)		  :: pt

	select case ( f%name )
	case ('Y')
	  F_at_point  = SphHarm(f%l,f%m,pt%phi,pt%theta)
	case default
	  write(0,*) 'Error: (F_at_point) sorry, only spherical harmonics are currently implemented'
	  stop		
	end select

  end function F_at_point


  ! Decide whether given radius is in the given parametrization layer

  pure logical function in_layer(rad,layer)

	type (layer_info), intent(in)					:: layer
	real(8), intent(in)								:: rad

	! Consistent with <sigma> layer definitions
	if ((rad>layer%lbound).and.(rad<=layer%ubound)) then
	  
	  in_layer = .TRUE.

	else
	
	  in_layer = .FALSE.

	end if

  end function in_layer


end module modeldef
