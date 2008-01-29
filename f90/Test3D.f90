!*****************************************************************************
program Test3D
!  program for testing 3D forward modeling,
!   and sensitivity codes ....    Modified from TETMtest

     use modelparameter
     use dataspace
     use solnrhs
     use iomod
     use dataFunc
     use sensmatrix
     use mtinvsetup
     use nlcg
        
     implicit none
	
     !   Sticking with Egbert I/O for initial tests
     ! I/O units ... reuse generic read/write units if 
     !   possible
     integer (kind=4) :: fidRead = 1
     integer (kind=4) :: fidWrite = 2
     integer (kind=4) :: fidError = 99

     ! logical vbles to control which tests are done
     !   not all of these are implemented in SensMatrix now
     logical		:: IOtest = .false.
     logical		:: FWDtest = .false.
     logical		:: fwdPredTest = .false.
     logical		:: SensMatrixCalcTest = .false.
     logical		:: MultBySensTest = .false.
     logical		:: MultBySensTransTest = .false.
     logical		:: linPredTest = .false.
     logical		:: DirectSensTest = .false.
     logical		:: nlcgTest = .true.
     logical		:: logCond = .true.

     real (kind=selectedPrec), dimension(:), pointer :: periods
     real (kind=selectedPrec), dimension(:,:), pointer :: sites
     character*2, dimension(:), pointer	:: modes
     character*80	:: paramType = ''

     ! grid geometry data structure
     type(grid3d_t), target 	:: grid

     ! impedance data structure
     type(dvecMTX)		:: allData

     !  storage for the "background" conductivity parameter
     type(modelParam_t)		:: sigma0
     !  storage for the inverse solution
     type(modelParam_t)		:: sigma1
     !  storage for the full sensitivity matrix
     type(modelParam_t),dimension(:), pointer	:: sigma
     !  storage for a perturbation to conductivity
     type(modelParam_t)		:: dsigma
     type(emsolve_control)       :: solverParams    

     !  storage for EM solutions
     type(EMsolnMTX)		:: eAll

     integer(kind=4) :: Nzb, IER, nPer, i, iy,iz,iPer,nSigma,ierr
     integer(kind=4) :: nSites
     character*80 	cfile, gridType, header
     real(kind=selectedPrec)	:: epsilon = 0.01
     real(kind=selectedPrec)	:: lambda
     real(kind=selectedPrec)	:: error = 0.05

     if(logCond) then
        paramType = LOG_CELL
     else
        paramType = CELL
     endif
     
     ! Read input files ....
     ! (1)  at present hard-wired for file names
     ! In StartUp3D xml startup file is read to get some needed file names:
     !     Name of model grid/background conductivity file (Mackie
     !      3D format); solver control; (other files in startup not 
     !      presently used in this program; could be used to set
     !      outputs).   Then model file is read and grid geometry
     !      is set, as well as background conductivity.

     cfile = 'startup'
     call StartUp3D(cfile,ParamType,grid,solverParams,sigma0)
     !  complete preliminary grid calculations and setup 
     call gridCalcs(grid)

     !  set grid for higher level solver routines
     call set_SolnRHS_grid(grid)
     
     !  set solver control
     call setEMsolveControl(solverParams)
    
     ! (2) Read in data file (only needs to be a template for
     !      most tests -- periods/sites/mode -- except for
     !   multiplication by J^T; in this case a full data vector,
     !   not a template, is needed.
     if(MultBySensTransTest) then
        if(logCond) then
	   cfile = 'TestPredLog3D.imp'
        else
           cfile = 'TestPred3D.imp'
        endif
     else if (nlcgTest) then
        cfile = 'true3d.imp'		
     else
        cfile = 'Template_mini.imp'
     endif
		
     !call read_Z3D_ascii(fidRead,cfile,nPer,periods,nSites,sites,allData)
     call read_Z3D(fidRead,cfile,nPer,periods,nSites,sites,allData)

     if(IOtest) then
        call write_Z3D_ascii(fidWrite,'out.imp',nPer,periods, &
            nSites,sites,allData)
     end if
	 
     !  Using periods, sites obtained from data file
     !     set up transmitter and receiver dictionaries
     !   DICTIONARIES (and setup routines)
     !        ARE NOW IN DataFunc2d ... probably this is not
     !        a good place, since modeling modules that don't
     !        care about data are using TXdict now
     call TXdictSetUp(nPer,periods) 
     call RXdictSetUp(nSites,sites)
     call TypeDictSetup()

     ! (3) if necessary, read in dsigma (perturbation to background
     !       conductivity parameter
     if((MultBySensTest .or.DirectSensTest)) then
        if(logCond) then
            cfile = 'dTestLog3D.cpr'
        else
            cfile = 'dTest3D.cpr'
        endif
        call create_ModelParam(grid,paramType,dsigma)
        call read_Cond3D(fidRead,cfile,dsigma)
     endif

     if(FWDtest) then
        call create_EMsolnMTX(allData,eAll) 
        call fwdPred(sigma0,allData,eAll)
        ! write out EM solutions
        if(logCond) then
	    cfile = 'TestLog3D.sol'
	else
            cfile = 'Test3D.sol'
	endif
        !  output of multiple  solutions ... just for testing, done 
        !   with minimal metadata!


       call write_EMsolnMTX3D(fidWrite,cfile,eAll)

        ! write out all impedances
        if(logCond) then
	    cfile = 'TestPredLog3D.imp'
	else
            cfile = 'TestPred3D.imp'
	endif
!        call write_Z3D_ascii(fidWrite,cfile,nPer,periods,  &
        call write_Z3D(fidWrite,cfile,nPer,periods,  &
			nSites,sites,allData)

     endif

     if(fwdPredTest) then
        ! don't save solutions, just impedances
        call fwdPred(sigma0,allData)
        ! write out all impedances
	if(logCond) then
	    cfile = 'TestPredLog3D_A.imp'
	else
            cfile = 'TestPred3D_A.imp'
	endif
!	    call setError_dvecMTX(error,allData)
!        call write_Z3D_ascii(fidWrite,cfile,nPer,periods,  &
        call write_Z3D(fidWrite,cfile,nPer,periods,  &
			nSites,sites,allData)
     endif


     if(SensMatrixCalcTest) then
        call calcSensMatrix(allData,sigma0,sigma)
	if(logCond) then
	    cfile = 'TestLog3D.sns'
	else
            cfile = 'Test3D.sns'
	endif
        header = 'Sensitivity matrix test' 
        call writeAll_Cond3D(fidWrite,cfile,header,   &
			allData%nData,sigma)
     endif

     if(MultBySensTest) then
       call Jmult(dsigma,sigma0,allData) 
       if(logCond) then
	    cfile = 'TestSensMultLog3D.imp'
	else
            cfile = 'TestSensMult3D.imp'
	endif
        call write_Z3D(fidWrite,cfile,nPer,periods,  &
		nSites,sites,allData)
     endif

     if(MultBySensTransTest) then
        call create_EMsolnMTX(allData,eAll) 
        call fwdPred(sigma0,allData,eAll)
       call JmultT(sigma0,allData,dsigma,eAll)
!       call JmultT(sigma0,allData,dsigma)
       if(logCond) then
            cfile = 'TestSensMultTransLog3D.cpr'
        else
            cfile = 'TestSensMultTrans3D.cpr'
        endif
        call write_Cond3D(fidWrite,cfile,dSigma)
     endif

!     if(linPredTest) then
!        call sensMatMult(dsigma,sigma0,allData,epsilon) 
!        cfile = 'LinPredTest3D.imp'
!        call write_Z3D(fidWrite,cfile,nPer,periods,   &
!			nSites,sites,allData)
!     endif
     
     if(DirectSensTest) then
        !  add sigma0 to epsilon * dsigma
        call linComb_modelParam(ONE,sigma0,epsilon,dsigma,dsigma)
        call fwdPred(dsigma,allData)
        ! write out impedances computed with perturbed sigma
	if(logCond) then  
	   cfile = 'TestDeltaPredLog3D.imp'
	else   
	   cfile = 'TestDeltaPred3D.imp'
	endif
        call write_Z3D(fidWrite,cfile,nPer,periods,&  
                  nSites,sites,allData)
     endif
       

     if(nlcgTest) then
        write(*,*) 'Starting the NLCG search...'
        lambda = .1
        call NLCGsolver(allData,lambda,sigma0,sigma1)
        if(logCond) then
            cfile = 'solution3d_log.cpr'
        else
            cfile = 'solution3d.cpr'
        endif
        call write_Cond3D(fidWrite,cfile,sigma1)
     endif

     ! deallocate
     call deall_modelParam(sigma0)
     call deall_modelParam(dsigma)

end program
