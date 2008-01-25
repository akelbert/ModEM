!*****************************************************************************
program TETMtest
!  program for testing 2D TE/TM forward modeling,
!   and sensitivity codes, up through dcg tests (not yet all implemented)
!    Modified April 27 from FWDtestTE

!   All new module names ....
     use modelparameter
     use dataspace
     !use solnrhs
     use wsfwd2d
     use iomod
     use dataFunc
     use sensmatrix
!     use dcg
	 !use emsolver !temporarily added by AK for deall_EMsolnMTX
	 use grid2d !temporarily added by AK for deall_Grid2D
	 use NLCG
        
     implicit none
	
     !   Sticking with Egbert I/O for initial tests
     ! I/O units ... reuse generic read/write units if 
     !   possible
     integer (kind=4) :: fidRead = 1
     integer (kind=4) :: fidWrite = 2
     integer (kind=4) :: fidError = 99
	 integer (kind=4) :: j

     ! logical vbles to control which tests are done
     !   not all of these are implemented in SensMatrix2D now
     logical		:: IOtest = .false.
     logical		:: FWDtest = .false.
     logical		:: fwdPredTest = .true.
     logical		:: SensMatrixCalcTest = .false.
     logical		:: MultBySensTest = .false.
     logical		:: MultBySensTransTest = .false.
     logical		:: linPredTest = .false.
     logical		:: DirectSensTest = .false.
     logical		:: nlcgTest = .false.
     logical		:: logCond = .true.

     real (kind=selectedPrec), dimension(:), pointer :: periods
     real (kind=selectedPrec), dimension(:,:), pointer :: sites
     character*2, dimension(:), pointer	:: modes
     character*2	:: cMode = 'TE'
     character*80	:: paramType = ''

     ! grid geometry data structure
     type(grid2d_t), target 	:: TEgrid

     ! impedance data structure
     type(dvecMTX)		:: allData, allResid

     !  storage for the "background" conductivity parameter
     type(modelParam_t)		:: sigma0
     !  storage for the full sensitivity matrix
     type(modelParam_t),dimension(:), pointer	:: sigma
     !  storage for a perturbation to conductivity
     type(modelParam_t)		:: dsigma
     !  storage for the penalty functional derivative
     type(modelParam_t)		:: grad
     !  storage for the inverse solution
     type(modelParam_t)		:: sigma1
     !  the damping parameter
     real (kind=selectedPrec) :: lambda

     !  storage for EM solutions
     type(EMsolnMTX)		:: eAll

     integer (kind=4) :: Nzb, IER, nPer, i, iy,iz,iPer,nSigma
     integer (kind=4) :: nMode=1,nComp=2,nSites
     character*80 	cfile, gridType, header
     real(kind=selectedPrec)	:: epsilon = 0.01
     real(kind=selectedPrec)	:: error = 0.05

     ! Read input files ....
     !  at present hard-wired for file names
     ! (1) Read in numerical grid geometry
     cfile = '2D_ascii.grd'
     call read_Grid2D_ascii(fidRead,cfile,TEgrid)
     !  complete grid definition ... 
     call gridCalcs(TEgrid)

     !  set array size parameters in WS forward code module 
     !   these stay fixed for all forward modeling with this grid
     call setWSparams(TEgrid%Ny,TEgrid%Nz,TEgrid%Nza)

     !  set grid for higher level solver routines
     call set_SolnRHS_grid(TEgrid)
    
     ! (2) Read background conductivity parameter (allocate first,
     !    using size info obtained from grid ... this might change
     !    for different parameters!)
     if(logCond) then
         cfile = 'background_log.cpr'
     else
         cfile = 'background.cpr'
     endif
		 cfile = 'background_ascii.cpr'
     call create_modelParam(TEgrid,paramType,sigma0)
     call read_Cond2D_ascii(fidRead,cfile,sigma0)
	
     ! (3) Read in data file (only needs to be a template for
     !      most tests -- periods/sites/mode -- except for
     !   multiplication by J^T; in this case a full data vector,
     !   not a template, is needed.
     if(MultBySensTransTest .or. nlcgTest) then
	       cfile = 'true_'//cMode//'.imp'
     else
        cfile = 'template_'//cMode//'.imp'
     endif
		 cfile = 'data_te_real_error_conj'
     call read_Z_ascii(fidRead,cfile,nPer,periods,modes,nSites,sites,allData)

	   if(IOtest) then
	 	   call write_Z_ascii(fidWrite,'out.imp',nPer,periods,modes,nSites,sites,allData)
	   end if

     !  Using periods, sites obtained from data file
     !     set up transmitter and receiver dictionaries
     !   DICTIONARIES (and setup routines)
     !        ARE NOW IN DataFunc2d ... probably this is not
     !        a good place, since modeling modules that don't
     !        care about data are using TXdict now
     call TXdictSetUp(nPer,periods,modes) 
     call RXdictSetUp(nSites,sites)
     call TypeDictSetup()

     ! (4) if necessary, read in dsigma (perturbation to background
     !       conductivity parameter
     if((MultBySensTest .or.DirectSensTest).or. &
         (MultBySensTransTest)) then
        if(logCond) then
            cfile = 'dTestLog.cpr'
        else
            cfile = 'dTest.cpr'
        endif
        call create_ModelParam(TEgrid,paramtype,dsigma)
        call read_Cond2D(fidRead,cfile,dsigma)
     endif

     if(FWDtest) then
        write(*,*) 'Calculating predicted data and saving the EM solution...'
        call create_EMsolnMTX(allData,eAll) 
        call fwdPred(sigma0,allData,eAll)
        ! write out EM solutions
        if(logCond) then
	    cfile = 'TestLog'//cMode//'.sol'
	else
            cfile = 'Test'//cMode//'.sol'
	endif
        !  output of multiple  solutions ... just for testing, done 
        !   with minimal metadata!

        call write_EMsolnMTX(fidWrite,cfile,eAll)

        ! write out all impedances
        if(logCond) then
	    cfile = 'fwdPredLog'//cMode//'.imp'
	else
            cfile = 'fwdPred'//cMode//'.imp'
	endif
        call write_Z(fidWrite,cfile,nPer,periods,modes,   &
			nSites,sites,allData)

     endif

     if(fwdPredTest) then
        write(*,*) 'Calculating predicted data...'
        ! don't save solutions, just impedances        
        call fwdPred(sigma0,allData)
        ! write out all impedances
	if(logCond) then
	    cfile = 'fwdPredLog'//cMode//'_A.imp'
	else
            cfile = 'fwdPred'//cMode//'_A.imp'
	endif
		    call setError_dvecMTX(error,allData)
        call write_Z_ascii(fidWrite,cfile,nPer,periods,modes, &
			nSites,sites,allData)
     endif

     if(SensMatrixCalcTest) then
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
	if(logCond) then
	    cfile = 'TestLog'//cMode//'.sns'
	else
            cfile = 'Test'//cMode//'.sns'
	endif
        header = 'Sensitivity matrix test' 
        call writeAll_Cond2D(fidWrite,cfile,header,   &
			allData%nData,sigma)
     endif

     if(MultBySensTest) then
       write(*,*) 'Multiplying by J...'
       call Jmult(dsigma,sigma0,allData) 
       if(logCond) then
	    cfile = 'TestSensMultLog'//cMode//'.imp'
	else
            cfile = 'TestSensMult'//cMode//'.imp'
	endif
        call write_Z(fidWrite,cfile,nPer,periods,modes,nSites,sites,allData)
     endif

     if(MultBySensTransTest) then
       write(*,*) 'Multiplying by J^T...'
       if(fwdTest) then
         call JmultT(sigma0,allData,dsigma,eAll)
       else
         call JmultT(sigma0,allData,dsigma)
       end if
       if(logCond) then
            cfile = 'TestSensMultTransLog'//cMode//'.cpr'
        else
            cfile = 'TestSensMultTrans'//cMode//'.cpr'
        endif
        call write_Cond2D(fidWrite,cfile,dSigma)
     endif

!     if(linPredTest) then
!        call sensMatMult(dsigma,sigma0,allData,epsilon) 
!        cfile = 'LinPredTest'//cMode//'.imp'
!        call write_Z(fidWrite,cfile,nPer,periods,modes,nSites,sites,allData)
!     endif
     
     if(DirectSensTest) then
        write(*,*) 'Calculating the predicted data for the perturbed model...'
        !  add sigma0 to epsilon * dsigma
        call linComb_modelParam(ONE,sigma0,epsilon,dsigma,dsigma)
        call fwdPred(dsigma,allData)
        ! write out impedances computed with perturbed sigma
	if(logCond) then  
	   cfile = 'TestDeltaPredLog'//cMode//'.imp'
	else   
	   cfile = 'TestDeltaPred'//cMode//'.imp'
	endif
        call write_Z(fidWrite,cfile,nPer,periods,modes,nSites,sites,allData)
     endif
 
     if(nlcgTest) then
       write(*,*) 'Starting the NLCG search...'
        lambda = ONE
        call NLCGsolver(allData,lambda,sigma0,sigma1)
        if(logCond) then
            cfile = 'solution_'//cMode//'_log.cpr'
        else
            cfile = 'solution_'//cMode//'.cpr'
        endif
        call write_Cond2D_ascii(fidWrite,cfile,sigma1)
     endif
       
     ! deallocate
     if (associated(periods)) deallocate(periods)
     if (associated(modes)) deallocate(modes)
     if (associated(sites)) deallocate(sites)
     if (associated(sigma)) then
     	do j=1,size(sigma)
        	call deall_modelParam(sigma(j))
     	end do
     	deallocate(sigma)
     end if
     call deall_modelParam(sigma0)
     call deall_modelParam(dsigma)
	 call deall_DvecMTX(allData)
	 call deall_EMsolnMTX(eAll)
	 call deall_Grid2D(TEgrid)
	 call deall_RXdict
	 call deall_Dict

end program
