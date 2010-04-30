! *****************************************************************************
module iotypes
  ! This module defines the derived data type structure with all filenames

  implicit none

  ! ***************************************************************************
  ! * input_info contains the list of all essential input information currently
  ! * read in from fn_startup.
  type :: input_info

	character(80)				:: paramname  ! 'harmonic'/'mixed' parametrization
	character(80)				:: modelname  ! specify how to call the output files
	character(80)				:: verbose  ! 'none'/'model'/'responses'/'all'/'debug'
	character(80)				:: calculate  ! 'responses'/'jacobian'/'derivative'
	character(80)				:: fn_grid	! grid information
	character(80)				:: fn_shell	! GM thin shell conductance values
	character(80)				:: fn_period  ! periods or frequencies
	character(80)				:: fn_coords  ! observatory coordinates file
	character(80)				:: fn_func  ! information about data functionals
	character(80)				:: fn_ctrl	! forward solver control
	character(80)				:: fn_invctrl	! inverse solver control
	character(80)				:: fn_slices  ! grid radii at which we output the data
	character(80)				:: fn_field	! input radial field solution (if required)
	character(80)				:: fn_precond	! preconditioning parameters
	character(80)				:: fn_param0	! base model parametrization
	character(80)				:: fn_param	! information about the parametrization
	character(80)				:: fn_cdata	! name of 1st data and observatory file
	character(80)				:: fn_ddata	! name of 2nd data and observatory file
	character(80)				:: fn_misfit  ! output file for data misfit
	character(80)				:: fn_gradient	! output file for derivative
	character(80)				:: fn_point	! output file for point parametrization
	real(8)             		:: damping ! the value of damping parameter mu
	real(8)             		:: step_size ! initial step size parameter for inversion

  end type input_info ! input_info


  ! ***************************************************************************
  ! * output_info contains the list of all output file names. Files such as the
  ! * boundary condition or starting solution, or diagnostic output are not
  ! * currently required. Will be added if necessary.
  type :: output_info

	! Output files:
	character(80)				:: fn_hx, fn_hy, fn_hz
	character(80)				:: fn_jxjyjz, fn_hxhyhz
	character(80)				:: fn_err
	character(80)				:: fn_cresp, fn_dresp
	character(80)				:: fn_cdat, fn_ddat
	character(80)				:: fn_cjac, fn_djac
	character(80)				:: fn_model
	character(80)				:: fn_residuals
	character(80)				:: fn_avg_cresp, fn_avg_dresp
	! Optional input file:
    character(80)				:: fn_bv

  end type output_info  ! output_info


end module iotypes
