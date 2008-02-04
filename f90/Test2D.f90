! *****************************************************************************
program Test2D
!  program for linking 2D TE forward modeling/data functional
!       etc. modules to matlab inversion drivers
!     Modified from FWDtestTE
!     Idea: a self contained program called with one or
!        more command line options to execute specific mappings
!        from model parameter space to data space.  A matlab
!        script writes files, executes this program (which reads
!        input data from these files, and outputs to other files),
!        and then reads the outputs.  Generalize to a mex file!

!   All new module names ....
     use modelparameter
     use dataspace
     use solnrhs
     use wsfwd2d
     use iobinary
     use dataFunc
     use sensmatrix
     use nlcg
        
     implicit none
	
     ! I/O units ... reuse generic read/write units if 
     !   possible; for those kept open during program run, 
     !   reserve a specific unit here 
     integer (kind=4) :: fidRead = 1
     integer (kind=4) :: fidWrite = 2
     integer (kind=4) :: fidError = 99

     ! characters used to comunicate (from matlab)
     !   about which jobs is executed
     character*1		:: job = ' '
     character*1, parameter	:: FORWARD_PRED = 'F'
     character*1, parameter	:: FORWARD_SOLN = 'E'
     character*1, parameter	:: SENSITIVITIES = 'S'
     character*1, parameter	:: MULT_BY_J = 'G'
     character*1, parameter	:: MULT_BY_J_T = 'T'
     character*1, parameter	:: MULT_BY_J_MTX = 'N'
     character*1, parameter	:: MULT_BY_J_T_MTX = 'M'
     character*1, parameter	:: FORWARD_PRED_DELTA = 'D'
     character*1, parameter	:: CREATE_DATA_ERRORS = 'C'
     character*1, parameter	:: CREATE_BG_MODEL = 'B'
     character*1, parameter	:: NLCG_INVERSION = 'I'

     ! file names
     character*80  :: ioConfigFile = ''     
     character*80  :: gridFile = 'Input.grd'
     character*80  :: inputModelFile = 'Input.cpr'
     character*80  :: outputModelFile = 'Out.cpr'
     character*80  :: inputDataFile = 'Input.imp'
     character*80  :: outputDataFile = 'Out.imp'
     character*80  :: EMsolnFile = 'Out.sol'
     character*80  :: SensFile = 'Out.sns'  !'Input2.sns'   
     character*80  :: fwdConfigFile = 'forward.cfg'
     character*80  :: invConfigFile = 'inverse.cfg'     

     real (kind=selectedPrec), dimension(:), pointer	:: periods
     real (kind=selectedPrec), dimension(:,:), pointer	:: sites
     character*2, dimension(:), pointer     		:: modes

     ! grid geometry data structure
     type(grid2d_t), target 	:: TEgrid

     ! impedance data structure
     type(dvecMTX)		:: allData

     !  storage for the "background" conductivity parameter
     type(modelParam_t)		:: sigma0
     !  storage for the full sensitivity matrix
     type(modelParam_t),dimension(:),pointer	:: sigma
     !  storage for a perturbation to conductivity
     type(modelParam_t)		:: dsigma
     !  storage for the inverse solution
     type(modelParam_t)		:: sigma1

     !  storage for EM solutions
     type(EMsolnMTX)            :: eAll

     !  the damping parameter
     real(kind=selectedPrec)    :: lambda

     real(kind=selectedPrec)	:: eps
     real(kind=selectedPrec)	:: error
     
     real(kind=selectedPrec),dimension(:,:), pointer  :: bg_cond

     integer (kind=4) :: Nzb, IER, nPer, i, iy,iz,iPer,nSigma,nTx
     integer (kind=4) :: iargc,narg,k
     integer (kind=4) :: nMode=1,nComp=2,nSites
     character*80  gridType, header,arg, paramtype

     !  parse command line ...  for now just uses first argument
     !   to set job
     narg = iargc()
     if(narg .gt. 0) then
        call getarg(1,arg)
        if(arg(1:1).eq.'-') job = arg(2:2)
     else
        write(*,*) 'Usage: Test2D -[option] [config_file]'
        write(*,*)
        write(*,*) ' Available options:'
        write(*,*) 'F - calculates the predicted data'
        write(*,*) 'E - calculates the predicted data and saves the EM solution'
        write(*,*) 'S - calculates and saves the full sensitivity matrix'
        write(*,*) 'G - multiplies a model by J to create a data vector'
        write(*,*) 'T - multiplies a data vector by J^T to create a model'
        write(*,*) 'N - evaluates sum( J m_i ) over transmitters to yield a data'
        write(*,*) '    vector; reads a sensitivity matrix from file to do this'
        write(*,*) 'M - multiples d_i by J_i^T separately for each transmitter, '
        write(*,*) '    to yield a bunch of models, one for each transmitter'
        write(*,*) 'D - calculated the predicted data for a perturbed model'
        write(*,*) 'C - creates data errors in a crude way (inversion testing)'
        write(*,*) 'B - creates background model for inversion from the grid'
        write(*,*) 'I - runs an NLCG inversion to yield an inverse model'
        write(*,*)
!        write(*,*) ' If config_file parameter is present, reads the input and'
!        write(*,*) 'output file names from this file. Format of config_file: '
!        write(*,*) 'CONFIG' 
!        write(*,*) 'FWD   :'
!        write(*,*) 'INV   :'
!        write(*,*) 'INPUT '
!        write(*,*) 'GRID  :'
!        write(*,*) 'MODEL :'
!        write(*,*) 'DATA  :'
!        write(*,*) 'EMSOLN:'
!        write(*,*) 'SENS  :'
!        write(*,*) 'OUTPUT '
!        write(*,*) 'GRID  :'
!        write(*,*) 'MODEL :'
!        write(*,*) 'DATA  :'
!        write(*,*) 'EMSOLN:'
!        write(*,*) 'SENS  :'
        write(*,*) 'Config_file parameter option is not yet active. Input and'
        write(*,*) 'output file names are hard coded at the top of Test2D.f90'
        write(*,*) 'Instead, for now you can pass the input arguments: input '
        write(*,*) 'grid file name, input model, output model, input data,   '
        write(*,*) 'output data, in this order.'        
        stop
     endif
     
     if (narg > 1) then
        call getarg(2,gridFile)
     end if

     if (narg > 2) then
        call getarg(3,inputModelFile)
     end if

     if (narg > 3) then
        call getarg(4,outputModelFile)
     end if

     if (narg > 4) then
        call getarg(5,inputDataFile)
     end if

     if (narg > 5) then
        call getarg(6,outputDataFile)
     end if
        
     ! Read input grid files ....
     !  hard-wired for file names used by matlab script
     ! Read in numerical grid geometry (needed in all cases)
     call read_Grid2D(fidRead,gridFile,TEgrid)
     !  complete grid definition ... 
     call gridCalcs(TEgrid)

     !  set array size parameters in WS forward code module 
     !   these stay fixed for all forward modeling with this grid
     call setWSparams(TEgrid%Ny,TEgrid%Nz,TEgrid%Nza)

     !  set grid for higher level solver routines
     call set_SolnRHS_grid(TEgrid)

     select case (job)
     
     case (CREATE_BG_MODEL)
        write(*,*) 'Creating background model for the inversion...'
        allocate(bg_cond(TEgrid%Ny,TEgrid%Nz-TEgrid%Nza))
        bg_cond = 1e-2
        call create_modelParam(TEgrid,paramtype,sigma1,bg_cond)
        call write_Cond2D(fidWrite,outputModelFile,sigma1,TEgrid)
        deallocate(bg_cond)
        stop 
     end select
     
     ! (2) Read background conductivity parameter (allocate first,
     !    using size info obtained from grid ... this might change
     !    for different parameters!)
     call read_Cond2D(fidRead,inputModelFile,sigma0,TEgrid)
	
     ! (3) Read in data file (only a template on input--periods/sites)
     call read_Z(fidRead,inputDataFile,nPer,periods,modes,nSites,sites,allData)

     !  Using periods, sites obtained from data file
     !     set up transmitter and receiver dictionaries
     call TXdictSetUp(nPer,periods,modes) 
     call RXdictSetUp(nSites,sites)
     call TypeDictSetup()

     select case (job)
     
     case (FORWARD_PRED)
        write(*,*) 'Calculating predicted data...'
        call fwdPred(sigma0,allData)
        ! write out impedances
        call write_Z(fidWrite,outputDataFile,nPer,periods,modes,   &
			nSites,sites,allData)

     case (FORWARD_SOLN)
        write(*,*) 'Calculating predicted data and saving the EM solution...'
        call create_EMsolnMTX(allData,eAll)
        call fwdPred(sigma0,allData,eAll)
        ! write out EM solutions
        call write_EMsolnMTX(fidWrite,EMsolnFile,eAll)
        ! write out all impedances
        call write_Z(fidWrite,outputDataFile,nPer,periods,modes,   &
			nSites,sites,allData)

     case (SENSITIVITIES)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
        header = 'Sensitivity Matrix' 
        call writeAll_Cond2D(fidWrite,SensFile,header,   &
                        allData%nData,sigma)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call read_Cond2D(fidRead,inputModelFile,dsigma,TEgrid)
        call Jmult(dsigma,sigma0,allData) 
        call write_Z(fidWrite,outputDataFile,nPer,periods,modes,   &
			nSites,sites,allData)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma) 
        call write_Cond2D(fidWrite,outputModelFile,dsigma,TEgrid)

     case (MULT_BY_J_T_MTX)
        write(*,*) 'Multiplying by J^T (all transmitters)...'
        call JmultT_MTX(sigma0,allData,sigma) 
        header = 'Sensitivity Matrix' 
        call writeAll_Cond2D(fidWrite,SensFile,header,   &
                        allData%nTx,sigma)

     case (MULT_BY_J_MTX)
        write(*,*) 'Multiplying by J (all transmitters)...'
        nTx = allData%nTx
        paramtype = ''
        allocate(sigma(nTx))
        do i = 1,nTx 
           call create_ModelParam(TEgrid,paramtype,sigma(i))
        enddo
        call readAll_Cond2D(fidRead,SensFile,nTx,header,sigma)
        call Jmult_MTX(sigma,sigma0,allData) 
        call write_Z(fidWrite,outputDataFile,nPer,periods,modes,   &
			nSites,sites,allData)
			
      case (FORWARD_PRED_DELTA)
        write(*,*) 'Calculating the predicted data for the perturbed model...'
        !  add sigma0 to epsilon * dsigma
        eps = 0.01
        call linComb_modelParam(ONE,sigma0,eps,dsigma,dsigma)
        call fwdPred(dsigma,allData)
        ! write out impedances computed with perturbed sigma
        call write_Z(fidWrite,outputDataFile,nPer,periods,modes,nSites,sites,allData)

     case (CREATE_DATA_ERRORS)
        write(*,*) 'Writing a data file with 5% errors in it...'
        error = 0.05
        call setError_dvecMTX(error,allData)
        ! write out impedances
        call write_Z(fidWrite,outputDataFile,nPer,periods,modes,   &
			nSites,sites,allData)
			     		
     case (NLCG_INVERSION)
        write(*,*) 'Starting the NLCG search...'
        lambda = ONE
        call NLCGsolver(allData,lambda,sigma0,sigma1)
        call write_Cond2D(fidWrite,outputModelFile,sigma1,TEgrid)
     end select

end program