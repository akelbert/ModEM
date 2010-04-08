! *****************************************************************************
module Main
	! These subroutines are called from the main program only

  use modelparameter
  use dataspace ! dataVecMTX_t
  use datafunc ! to deallocate rxDict
  use emsolver ! txDict, EMsolnMTX
  use sensmatrix
  use userctrl
  use ioascii
  use dataio
  implicit none

      ! I/O units ... reuse generic read/write units if
     !   possible; for those kept open during program run,
     !   reserve a specific unit here
     integer (kind=4), save :: fidRead = 1
     integer (kind=4), save :: fidWrite = 2
     integer (kind=4), save :: fidError = 99

  ! ***************************************************************************
  ! * fwdCtrls: User-specified information about the forward solver relaxations
  !type (emsolve_control), save								:: fwdCtrls
  !type (inverse_control), save								:: invCtrls

  ! forward solver control defined in EMsolve3D
  type(EMsolve_control)  :: solverParams

  ! this is used to set up the numerical grid in SensMatrix
  type(grid_t), save	        :: grid

  ! impedance data structure
  type(dataVecMTX_t), save		:: allData

  !  storage for the "background" conductivity parameter
  type(modelParam_t), save		:: sigma0
  !  storage for a perturbation to conductivity
  type(modelParam_t), save		:: dsigma
  !  storage for the inverse solution
  type(modelParam_t), save		:: sigma1

  !  storage for the full sensitivity matrix (dimension nTx)
  type(sensMatrix_t), pointer, dimension(:), save	:: sens

  !  storage for EM solutions
  type(EMsolnMTX_t), save            :: eAll

  logical                   :: write_model, write_data, write_EMsoln



Contains

  ! ***************************************************************************
  ! * InitGlobalData is the routine to call to initialize all derived data types
  ! * and other variables defined in modules basics, modeldef, datadef and
  ! * in this module. These include:
  ! * 1) constants
  ! * 2) grid information: nx,ny,nz,x,y,z, cell centre coordinates
  ! * 3) model information: layers, coefficients, rho on the grid
  ! * 4) periods/freq info: nfreq, freq
  ! * 5) forward solver controls
  ! * 6) output file names

  subroutine initGlobalData(cUserDef)

	implicit none
	type (userdef_control), intent(in)          :: cUserDef
	integer										:: i,ios=0,istat=0
	logical                                     :: exists
	character(80)                               :: paramtype,header
	integer                                     :: nTx

	!--------------------------------------------------------------------------
	! Set global output level stored with the file units
	output_level = cUserDef%output_level

	!--------------------------------------------------------------------------
	! Check whether model parametrization file exists and read it, if exists
	inquire(FILE=cUserDef%rFile_Model,EXIST=exists)

	if (exists) then
	   ! Read background conductivity parameter and grid
       call read_modelParam(grid,sigma0,cUserDef%rFile_Model)

       ! Finish setting up the grid (if that is not done in the read subroutine)
       call setup_grid(grid)

       !  set solver control (currently using defaults)
       solverParams%UseDefaults= .true.
       call setEMsolveControl(solverParams)
	else
	  call warning('No input model parametrization')
	end if

	!--------------------------------------------------------------------------
	!  Read in data file (only a template on input--periods/sites)
	inquire(FILE=cUserDef%rFile_Data,EXIST=exists)

	if (exists) then
       !  This also sets up dictionaries
	   call read_dataVecMTX(allData,cUserDef%rFile_Data)
    else
       call warning('No input data file - unable to set up dictionaries')
    end if

	!--------------------------------------------------------------------------
	!  Initialize additional data as necessary
	select case (cUserDef%job)

     case (MULT_BY_J)
	   inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
	   if (exists) then
	      call deall_grid(grid)
	   	  call read_modelParam(grid,dsigma,cUserDef%rFile_dModel)
	   else
	      call warning('The input model perturbation file does not exist')
	   end if

     case (INVERSE)
	   inquire(FILE=cUserDef%rFile_Cov,EXIST=exists)
	   if (exists) then
          call create_CmSqrt(sigma0,cUserDef%rFile_Cov)
       else
          call create_CmSqrt(sigma0)
       end if
	   inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
	   if (exists) then
	      call deall_grid(grid)
	   	  call read_modelParam(grid,dsigma,cUserDef%rFile_dModel)
	      if (output_level > 0) then
	        write(0,*) 'Using the initial model perturbations from file ',trim(cUserDef%rFile_dModel)
	      endif
	   else
	      dsigma = sigma0
	      call zero_modelParam(dsigma)
	      if (output_level > 0) then
	        write(0,*) 'Starting search from the prior model ',trim(cUserDef%rFile_Model)
	      endif
	   end if
       select case (cUserDef%search)
       case ('NLCG')
       case ('DCG')
       case ('Hybrid')
       case default
          ! a placeholder for anything specific to a particular inversion algorithm;
          ! currently empty
       end select

     case (TEST_COV)
	   inquire(FILE=cUserDef%rFile_Cov,EXIST=exists)
	   if (exists) then
          call create_CmSqrt(sigma0,cUserDef%rFile_Cov)
       else
          call create_CmSqrt(sigma0)
       end if
       sigma1 = sigma0
       call zero(sigma1)

    end select

	!--------------------------------------------------------------------------
	!  Only write out these files if the output file names are specified
    write_model = .false.
    if (len_trim(cUserDef%wFile_Model)>1) then
       write_model = .true.
    end if
    write_data = .false.
    if (len_trim(cUserDef%wFile_Data)>1) then
       write_data = .true.
    end if
    write_EMsoln = .false.
    if (len_trim(cUserDef%wFile_EMsoln)>1) then
       write_EMsoln = .true.
    end if

	return

  end subroutine initGlobalData	! initGlobalData


  ! ***************************************************************************
  ! * DeallGlobalData deallocates all allocatable data defined globally.

  subroutine deallGlobalData()

	integer	:: i, istat

	! Deallocate global variables that have been allocated by InitGlobalData()
	if (output_level > 3) then
	   write(0,*) 'Cleaning up grid...'
	endif
	call deall_grid(grid)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up data...'
	endif
	call deall_dataVecMTX(allData)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up EM soln...'
	endif
	call deall_EMsolnMTX(eAll)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up models...'
	endif
	call deall_modelParam(sigma0)
	call deall_modelParam(dsigma)
	call deall_modelParam(sigma1)

    if (output_level > 3) then
       write(0,*) 'Cleaning up dictionaries...'
    endif
	call deall_txDict() ! 3D_MT/DICT/transmitters.f90
	call deall_rxDict() ! 3D_MT/DICT/receivers.f90
	call deall_typeDict() ! 3D_MT/DICT/dataTypes.f90

	if (associated(sens)) then
		call deall_sensMatrixMTX(sens)
	end if

	call deallEMsolveControl() ! 3D_MT/FWD/EMsolve3D.f90

	call deall_CmSqrt()

    if (output_level > 3) then
       write(0,*) 'All done.'
    endif

  end subroutine deallGlobalData	! deallGlobalData


end module Main
