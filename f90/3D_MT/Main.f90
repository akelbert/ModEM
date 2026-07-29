! *****************************************************************************
module Main
	! These subroutines are called from the main program only

  use ModelSpace
  use dataspace ! dataVectorMTX_t
  use datafunc ! to deallocate rxDict
  use ForwardSolver ! txDict, solnVectorMTX
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
  type(emsolve_control),save  :: solverParams

  ! this is used to set up the numerical grid in SensMatrix
  type(grid_t), save	        :: grid

  ! air layers might be set from a file, but can also use the defaults
  type(airLayers_t), save       :: airLayers

  ! impedance data structure
  type(dataVectorMTX_t), save		:: allData

  !  storage for the "background" conductivity parameter
  type(modelParam_t), save		:: sigma0
  !  storage for a perturbation to conductivity
  type(modelParam_t), save		:: dsigma
  !  storage for the inverse solution
  type(modelParam_t), save		:: sigma1
  !  currently only used for TEST_GRAD feature (otherwise, use allData)
  type(dataVectorMTX_t), save       :: predData
  !  also for TEST_GRAD feature...
  type(modelParam_t), save      :: sigmaGrad
  real(kind=prec), save         :: rms,mNorm,f1,f2,alpha
  !  storage for multi-Tx outputed from JT computation
  type(modelParam_t),pointer, dimension(:), save :: JT_multi_Tx_vec

  !  storage for the full sensitivity matrix (dimension nTx)
  type(sensMatrix_t), pointer, dimension(:), save	:: sens

  !  storage for EM solutions
  type(solnVectorMTX_t), save            :: eAll

  !  storage for EM rhs (currently only used for symmetry tests)
  !type(rhsVectorMTX_t), save            :: bAll

  logical                   :: write_model, write_data, write_EMsoln, write_EMrhs



Contains
  ! ***************************************************************************
  ! reads the "large grid" file for nested modeling; Naser Meqbel's variables
  ! Larg_Grid, eAll_larg and nTx_nPol are stored in the ForwardSolver module.
  ! This inherently assumes an equal number of polarizations for all transmitters
   subroutine read_Efield_from_file(solverControl)
    type(emsolve_control), intent(in)           :: solverControl
    ! local
    character (len=80)			                :: inFile
    character (len=20)                          :: fileVersion
    logical                                     :: exists

    nestedEM_initialized = .false.
    
    if (solverControl%read_E0_from_File) then
      inquire(FILE=solverControl%E0fileName,EXIST=exists)
      if (exists) then
        write(6,*) 'The Master reads E field from: ',trim(solverControl%E0fileName)
        inFile = trim(solverControl%E0fileName)
        call read_solnVectorMTX(Larg_Grid,eAll_larg,inFile)
        nTx_nPol=eAll_larg%nTx*eAll_larg%solns(1)%nPol
        nestedEM_initialized = .true.
      end if
     end if

    end subroutine read_Efield_from_file

  !**********************************************************************
  !   rewrite the defaults in the air layers structure
  subroutine  initAirLayers(solverControl,airLayers)
     type(emsolve_control), intent(in)    :: solverControl
     type(airLayers_t), intent(inout)     :: airLayers

     integer status

     if (.not.solverControl%AirLayersPresent) then
        ! do nothing - keep the defaults
        return
     else
        deallocate(airLayers%Dz, STAT=status)
     end if

     airLayers%method = solverControl%AirLayersMethod
     airLayers%Nz = solverControl%AirLayersNz
     allocate(airLayers%Dz(airLayers%Nz), STAT=status)
     airLayers%allocated = .true.

     if (index(airLayers%method,'mirror')>0) then
        airLayers%alpha = solverControl%AirLayersAlpha
        airLayers%MinTopDz = 1.e3*solverControl%AirLayersMinTopDz
     elseif (index(airLayers%method,'fixed height')>0) then
        airLayers%MaxHeight = 1.e3*solverControl%AirLayersMaxHeight
     elseif (index(airLayers%method,'read from file')>0) then
        airLayers%Dz = 1.e3*solverControl%AirLayersDz
     end if

  end subroutine initAirLayers

  ! ***************************************************************************
  ! * Configure the secondary field formulation (SFF) from the optional
  ! * namelist (&grid SFF=.true.), used to force SFF for -F / -I whose command
  ! * line does not carry a primary. Interprets primary_model / primary_field,
  ! * sets the ForwardSolver flags PRIMARY_MODEL_FROM_FILE / PRIMARY_E_FROM_FILE,
  ! * and copies the namelist file paths into the canonical rFile_Model1D /
  ! * rFile_EMsoln consumed by the primary-load path. Only three combinations
  ! * are supported; the two internally-computed-primary pathways are reserved
  ! * (deactivated) for now.
  subroutine configure_SFF_from_namelist(cUserDef)

    type(userdef_control), intent(inout) :: cUserDef
    logical       :: exists
    character(80) :: pmodel, pfield

    pmodel = adjustl(cUserDef%primary_model)
    pfield = adjustl(cUserDef%primary_field)

    ! --- source of the 1D background model sigma1D ---
    select case (trim(pmodel))
    case ('file')
        PRIMARY_MODEL_FROM_FILE = .true.
        cUserDef%rFile_Model1D = cUserDef%primary_model_file
        inquire(FILE=cUserDef%rFile_Model1D, EXIST=exists)
        if (.not. exists) then
            call errStop('SFF namelist: primary_model_file not found: '//trim(cUserDef%rFile_Model1D))
        end if
    case ('average')
        PRIMARY_MODEL_FROM_FILE = .false.   ! laterally average the 3D model
    case default
        call errStop("SFF namelist: primary_model must be 'file' or 'average' (got '"//trim(pmodel)//"')")
    end select

    ! --- source of the primary E-field ---
    select case (trim(pfield))
    case ('file')
        PRIMARY_E_FROM_FILE = .true.
        cUserDef%rFile_EMsoln = cUserDef%primary_field_file
        inquire(FILE=cUserDef%rFile_EMsoln, EXIST=exists)
        if (.not. exists) then
            call errStop('SFF namelist: primary_field_file not found: '//trim(cUserDef%rFile_EMsoln))
        end if
    case ('compute')
        PRIMARY_E_FROM_FILE = .false.       ! compute the 1D primary internally
    case default
        call errStop("SFF namelist: primary_field must be 'file' or 'compute' (got '"//trim(pfield)//"')")
    end select

    ! --- validate the supported combinations & guard the reserved pathways ---
    if (PRIMARY_MODEL_FROM_FILE .and. PRIMARY_E_FROM_FILE) then
        ! (1) primary_model='file', primary_field='file' -- ACTIVE
        continue
    else if (PRIMARY_MODEL_FROM_FILE .and. (.not. PRIMARY_E_FROM_FILE)) then
        ! (2) primary_model='file', primary_field='compute' -- reserved
        call errStop('SFF mode (primary_model=file, primary_field=compute) is not yet implemented')
    else if ((.not. PRIMARY_MODEL_FROM_FILE) .and. (.not. PRIMARY_E_FROM_FILE)) then
        ! (3) primary_model='average', primary_field='compute' -- reserved
        call errStop('SFF mode (primary_model=average, primary_field=compute) is not yet implemented')
    else
        ! primary_model='average', primary_field='file' -- not a supported pairing
        call errStop('SFF mode (primary_model=average, primary_field=file) is not supported')
    end if

  end subroutine configure_SFF_from_namelist

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
	! intent(inout): for a namelist-driven SFF run we copy the primary
	! model/field file paths into the canonical rFile_Model1D/rFile_EMsoln
	type (userdef_control), intent(inout)       :: cUserDef
	integer										:: i,ios=0,istat=0
	logical                                     :: exists
	character(80)                               :: paramtype,header
	integer                                     :: nTx
! Variables used to read E solution from a file.
	integer                    					:: filePer
    integer			                            :: fileMode,ioNum
    character (len=20)                          :: fileVersion
    character (len=80)			                :: inFile
    Integer                                     :: iTx,iMod
    type(solnVector_t)                          :: e_temp
    real (kind=prec)                        	:: Omega
	!--------------------------------------------------------------------------
	! Set global output level stored with the file units
	output_level = cUserDef%output_level

	!--------------------------------------------------------------------------
    !  Read forward solver control in EMsolve3D (or use defaults)
    call readEMsolveControl(solverParams,cUserDef%rFile_fwdCtrl,exists,cUserDef%eps)

    !  If solverParams contains air layers information, rewrite the defaults here
    call initAirLayers(solverParams,airLayers)

    !--------------------------------------------------------------------------
    ! Determine whether or not there is an input BC file to read
    inquire(FILE=cUserDef%rFile_EMrhs,EXIST=exists)
    if (exists) then
       BC_FROM_RHS_FILE = .true.
    end if
    ! ... or whether we're setting up BCs from a nested grid
    if (solverParams%read_E0_from_File) then
      inquire(FILE=solverParams%E0fileName,EXIST=exists)
      if (exists) then
        NESTED_BC = .true.
      end if
    end if
    ! ---------------------------------------------------------------------
    ! Decide the forward COMPUTATION mechanism (USE_SFF) and how the SFF
    ! primary is supplied -- the 1D background model sigma1D and the primary
    ! E-field. The secondary field formulation is activated either by:
    !   (a) the -E (SECONDARY_FIELD) job, whose command line carries the
    !       primary model (rFile_Model1D) and E-field (rFile_EMsoln) as
    !       positional arguments; or
    !   (b) the optional namelist SFF=.true., which lets -F / -I run the SFF
    !       for data types (MT/CSEM/TIDE/GLOBAL) whose command line does NOT
    !       carry a primary -- the primary is then described by the namelist
    !       primary_model / primary_field (and their *_file paths).
    ! Three primary-supply combinations are supported (see
    ! configure_SFF_from_namelist):
    !   (1) primary_model='file',    primary_field='file'    -- ACTIVE
    !   (2) primary_model='file',    primary_field='compute' -- reserved (stub)
    !   (3) primary_model='average', primary_field='compute' -- reserved (stub)
    ! ---------------------------------------------------------------------
    inquire(FILE=cUserDef%rFile_EMsoln,EXIST=exists)
    if (cUserDef%job==SECONDARY_FIELD) then
       ! -E job: primary 1D model and E-field come from the positional args.
       if (.not. exists) then
          call errStop('SECONDARY_FIELD (-E) requires a primary E-field solution file (rFile_EMsoln)')
       end if
       USE_SFF = .true.
       PRIMARY_MODEL_FROM_FILE = .true.
       PRIMARY_E_FROM_FILE = .true.
    elseif (cUserDef%SFF) then
       ! Namelist-driven SFF override (e.g. for -F / -I).
       USE_SFF = .true.
       call configure_SFF_from_namelist(cUserDef)
    elseif (exists) then
       ! Non-SFF job with an E-field file present: use it for BCs only.
       BC_FROM_E0_FILE = .true.
    end if

	!--------------------------------------------------------------------------
	! Check whether model parametrization file exists and read it, if exists
	inquire(FILE=cUserDef%rFile_Model,EXIST=exists)

	if (exists) then
	   ! Read background conductivity parameter and grid; complete airLayers setup
       call read_modelParam(grid,airLayers,sigma0,cUserDef%rFile_Model)
	else
	  call warning('No input model parametrization')

	  ! set up an empty grid to avoid segmentation faults in sensitivity tests
	  call create_grid(1,1,1,1,grid)
	end if


	!--------------------------------------------------------------------------
	!  Read in data file (only a template on input--periods/sites)
	inquire(FILE=cUserDef%rFile_Data,EXIST=exists)

	if (exists) then
       !  This also sets up dictionaries
	   call read_dataVectorMTX(allData,cUserDef%rFile_Data)
    else
       call warning('No input data file - unable to set up dictionaries')
    end if
    
    !--------------------------------------------------------------------------
    ! Allocate bAll in all cases except if we compute the BCs internally
    if (.not. COMPUTE_BC) then
        call create_rhsVectorMTX(allData%nTx,bAll)
    end if

	!--------------------------------------------------------------------------
	!  Initialize additional data as necessary
	select case (cUserDef%job)

     case (MULT_BY_J)
	   inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
	   if (exists) then
	      call deall_grid(grid)
	      call read_modelParam(grid,airLayers,dsigma,cUserDef%rFile_dModel)
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
	      call read_modelParam(grid,airLayers,dsigma,cUserDef%rFile_dModel)
	      if (output_level > 0) then
	        write(*,*) 'Using the initial model perturbations from file ',trim(cUserDef%rFile_dModel)
	      endif
	   else
	      dsigma = sigma0
	      call zero_modelParam(dsigma)
	      if (output_level > 0) then
	        write(*,*) 'Starting search from the prior model ',trim(cUserDef%rFile_Model)
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

     case (APPLY_COV)
       inquire(FILE=cUserDef%rFile_Cov,EXIST=exists)
       if (exists) then
          call create_CmSqrt(sigma0,cUserDef%rFile_Cov)
       else
          call create_CmSqrt(sigma0)
       end if
       dsigma = sigma0
       inquire(FILE=cUserDef%rFile_Prior,EXIST=exists)
       if (exists) then
           call deall_grid(grid)
           call read_modelParam(grid,airLayers,sigma0,cUserDef%rFile_Prior)
       else
           call zero(sigma0)
       end if
       sigma1 = sigma0
       call zero(sigma1)

     case (DATA_FROM_E) ! already reading the E field under the BC_FROM_E0_FILE flag
         inquire(FILE=cUserDef%rFile_EMsoln,EXIST=exists)
         if (.not. exists) then
            call warning('Cannot compute data from the input EM solution - file does not exist')
         end if

     case (TEST_GRAD, TEST_SENS)
         inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
         if (exists) then
             call deall_grid(grid)
             call read_modelParam(grid,airLayers,dsigma,cUserDef%rFile_dModel)
         else
             call warning('The input model perturbation file does not exist')
         end if

     case (TEST_ADJ)
       select case (cUserDef%option)
           case('J','Q','P')
               inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
               if (exists) then
                  call deall_grid(grid)
                  call read_modelParam(grid,airLayers,dsigma,cUserDef%rFile_dModel)
               else
                  call warning('The input model perturbation file does not exist')
               end if
           case default
       end select
       select case (cUserDef%option)
           case('L','P','e')
               inquire(FILE=cUserDef%rFile_EMsoln,EXIST=exists)
               if (exists) then
                  call read_solnVectorMTX(grid,eAll,cUserDef%rFile_EMsoln)
               else
                  call warning('The input EM solution file does not exist')
               end if
           case('S','b')
               inquire(FILE=cUserDef%rFile_EMrhs,EXIST=exists)
               if (exists) then
                  call read_rhsVectorMTX(grid,bAll,cUserDef%rFile_EMrhs)
               else
                  call warning('The input EM RHS file does not exist')
               end if
           case default
       end select
!       select case (cUserDef%option)
!           case('S')
!            call create_rhsVectorMTX(allData%ntx,bAll)
!            do iTx = 1,allData%ntx
!                bAll%combs(iTx)%nonzero_source = .true.
!                call create_rhsVector(grid,iTx,bAll%combs(iTx))
!            end do
!            call random_rhsVectorMTX(bAll,cUserDef%eps)
!       end select
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
    write_EMrhs = .false.
    if (len_trim(cUserDef%wFile_EMrhs)>1) then
       write_EMrhs = .true.
    end if

	return

  end subroutine initGlobalData	! initGlobalData


  subroutine ModEM_setup_IO(data_input_type, data_output_type, &
                            model_input_type, model_output_type)

    implicit none

    character(len=*), intent(in) :: data_input_type
    character(len=*), intent(in) :: data_output_type
    character(len=*), intent(in) :: model_input_type
    character(len=*), intent(in) :: model_output_type

    integer :: hdferr

    call setup_dataIO(data_input_type, data_output_type)
    call setup_modelParamIO(model_input_type, model_output_type)

  end subroutine ModEM_setup_IO


  ! ***************************************************************************
  ! * DeallGlobalData deallocates all allocatable data defined globally.

  subroutine deallGlobalData()

	integer	:: i, istat

    write(6,*) 'Cleaning up...'

	! Deallocate global variables that have been allocated by InitGlobalData()
	if (output_level > 3) then
	   write(0,*) 'Cleaning up grid...'
	endif
	call deall_grid(grid)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up data...'
	endif
	call deall_dataVectorMTX(allData)
	call deall_dataFileInfo()

	if (output_level > 3) then
	   write(0,*) 'Cleaning up EM soln...'
	endif
	call deall_solnVectorMTX(eAll)

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

	call deallEMsolveControl(solverParams) ! 3D_MT/FWD/EMsolve3D.f90

	call deall_CmSqrt()

    if (output_level > 3) then
       write(0,*) 'All done.'
    endif

  end subroutine deallGlobalData	! deallGlobalData


end module Main
