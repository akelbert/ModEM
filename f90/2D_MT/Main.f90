! *****************************************************************************
module Main
	! These subroutines are called from the main program only

  use modelparameter
  use dataspace
  use emsolver
  use userctrl
  use ioascii
  implicit none
  
      ! I/O units ... reuse generic read/write units if 
     !   possible; for those kept open during program run, 
     !   reserve a specific unit here 
     integer (kind=4), save :: fidRead = 1
     integer (kind=4), save :: fidWrite = 2
     integer (kind=4), save :: fidError = 99

     integer (kind=4), save :: nPer, nSites
       
  ! ***************************************************************************
  ! * cUserDef: Vital character-based information specified by the user
  ! * fwdCtrls: User-specified information about the forward solver relaxations
  type (userdef_control), save								:: cUserDef
  !type (emsolve_control), save								:: fwdCtrls
  !type (inverse_control), save								:: invCtrls
  
  integer, save                                             :: output_level
  
  real (kind=selectedPrec), dimension(:), pointer, save	:: periods
  real (kind=selectedPrec), dimension(:,:), pointer, save	:: sites
  character*2, dimension(:), pointer, save    		:: modes

  ! grid geometry data structure
  type(grid2d_t), target, save	:: TEgrid

  ! impedance data structure
  type(dvecMTX), save		:: allData

  !  storage for the "background" conductivity parameter
  type(modelParam_t), save		:: sigma0
  !  storage for the full sensitivity matrix
  type(modelParam_t),dimension(:), pointer, save	:: sigma
  !  storage for a perturbation to conductivity
  type(modelParam_t), save		:: dsigma
  !  storage for the inverse solution
  type(modelParam_t), save		:: sigma1

  !  storage for EM solutions
  type(EMsolnMTX), save            :: eAll


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
  ! * The information about the grid, frequency and resistivity map is saved
  ! * in the module global, and the grid is also saved in module griddef.
  ! * Frequencies are listed in ascending order (by value)
  ! * Observatories are listed in ascending order (alphabetic)

  subroutine initGlobalData()

	implicit none
	integer										:: i,ios=0,istat=0
	logical                                     :: exists

	!--------------------------------------------------------------------------
	! Check whether model parametrization file exists and read it, if exists
	inquire(FILE=cUserDef%rFile_Model,EXIST=exists) 

	if (exists) then
	   ! Read background conductivity parameter and grid
       call read_Cond2D(fidRead,cUserDef%rFile_Model,sigma0,TEgrid)
     
       !  set array size parameters in WS forward code module 
       !   these stay fixed for all forward modeling with this grid
       call setWSparams(TEgrid%Ny,TEgrid%Nz,TEgrid%Nza)

       !  set grid for higher level solver routines
       call set_SolnRHS_grid(TEgrid)
	else
	  call warning('No input model parametrization')
	end if
	
	!--------------------------------------------------------------------------
	!  Read in data file (only a template on input--periods/sites)
	inquire(FILE=cUserDef%rFile_Data,EXIST=exists) 

	if (exists) then
       call read_Z(fidRead,cUserDef%rFile_Data,nPer,periods,modes,nSites,sites,allData)
       !  Using periods, sites obtained from data file
       !     set up transmitter and receiver dictionaries
       call TXdictSetUp(nPer,periods,modes) 
       call RXdictSetUp(nSites,sites)
       call TypeDictSetup()
    else
       call warning('No input data file - unable to set up dictionaries')
    end if

	!--------------------------------------------------------------------------
	!  Set the level of output based on the user input
	select case (cUserDef%verbose)
	case ('debug')
	  print *,'Output all information including debugging lines.'
	  output_level = 5
	case ('full')
	  print *,'Output full information to screen and to files.'
	  output_level = 4
	case ('regular')
	  print *,'Output information to files, and progress report to screen (default).'
	  output_level = 3
	case ('compact')
	  print *,'Output information to files, and compact summary to screen.'
	  output_level = 2
	case ('result')
	  print *,'Output information to files, and result to screen.'
	  output_level = 1
	case ('none')
	  print *,'Output nothing at all except result to screen and to files.'
	  output_level = 0
	case default
	  output_level = 3
	end select

	return

  end subroutine initGlobalData	! initGlobalData


  ! ***************************************************************************
  ! * DeallGlobalData deallocates all allocatable data defined globally.
  
  subroutine deallGlobalData()

	integer	:: istat

	! Deallocate global variables that have been allocated by InitGlobalData()
	call deall_dvecMTX(allData)
	deallocate(modes,periods,sites,STAT=istat)
	
	!if (allocated(misfitValue)) then
	!  deallocate(misfitValue)
	!end if

  end subroutine deallGlobalData	! deallGlobalData
  

end module Main
