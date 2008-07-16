! *****************************************************************************
program Mod2DMT
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
     use ioascii
     use dataFunc
     use sensmatrix
     use main
     use nlcg
        
     implicit none

     ! file names
     character*80      :: rFile_Config = ''
     ! character*80      :: rFile_Grid = 'Test1.grd'
     character*80      :: rFile_Model = 'scratch/Input1.cpr'
     character*80      :: rFile_dModel = 'scratch/Input2.cpr'
     character*80      :: rFile_dModelMTX = 'scratch/Input2.sns'
     character*80      :: rFile_Data = 'scratch/Input.imp'
     character*80      :: wFile_Model = 'scratch/Out.cpr'
     character*80      :: wFile_dModel = 'scratch/Out.cpr'
     character*80      :: wFile_dModelMTX = 'scratch/Out.sns'
     character*80      :: wFile_Data = 'scratch/Out.imp'
     character*80      :: wFile_EMsoln = 'scratch/Out.sol'
     character*80      :: wFile_Sens = 'scratch/Out.sns'
     !type (ioFiles_t)  :: ioFiles

	 real					:: rtime  ! run time
	 real					:: ftime  ! run time per frequency
	 real					:: stime, etime ! start and end times
     integer, dimension(8)	:: tarray ! utility variable
        
     !  the damping parameter
     real(kind=selectedPrec)    :: lambda = 1.0
     real(kind=selectedPrec)    :: alpha = 0.1

     real(kind=selectedPrec)	:: eps
     real(kind=selectedPrec)	:: error
     
     real(kind=selectedPrec)	:: bg_cond_value
     real(kind=selectedPrec),dimension(:,:), pointer  :: bg_cond

     integer (kind=4) :: Nzb, IER, i, iy,iz,iPer,nSigma,nTx
     integer (kind=4) :: iargc,narg,k
     integer (kind=4) :: nMode=1,nComp=2
     character*80  gridType, header,arg, paramtype
     character*80, dimension(:), pointer :: temp
     
     call parseArgs('Mod2DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)
     
     call initGlobalData()
     
	 ! Start the (portable) clock
	 call date_and_time(values=tarray)
     stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
             
     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_Cond2D(fidWrite,cUserDef%wFile_Model,sigma0,TEgrid)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
				nSites,sites,allData)
		else if (write_model) then
        	write(*,*) 'Writing model and exiting...'
        	call write_Cond2D(fidWrite,cUserDef%wFile_Model,sigma0,TEgrid)		
		end if
        
     case (FORWARD)
        write(*,*) 'Calculating predicted data...'
        call fwdPred(sigma0,allData,eAll)
        if (write_EMsoln) then
        	! write out EM solutions
        	write(*,*) 'Saving the EM solution...'
        	call write_EMsolnMTX(fidWrite,cUserDef%wFile_EMsoln,eAll)
        end if
        ! write out all impedances
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
        header = 'Sensitivity Matrix' 
        call writeAll_Cond2D(fidWrite,cUserDef%wFile_Sens,header,   &
                        allData%nData,sigma)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call read_Cond2D(fidRead,cUserDef%rFile_dModel,dsigma,TEgrid)
        call Jmult(dsigma,sigma0,allData) 
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma) 
        call write_Cond2D(fidWrite,cUserDef%wFile_dModel,dsigma,TEgrid)

     case (MULT_BY_J_MTX)
        write(*,*) 'Multiplying by J (all transmitters)...'
        nTx = allData%nTx
        paramtype = ''
        allocate(sigma(nTx))
        do i = 1,nTx 
           call create_ModelParam(TEgrid,paramtype,sigma(i))
        enddo
        call readAll_Cond2D(fidRead,cUserDef%rFile_dModelMTX,nTx,header,sigma)
        call Jmult_MTX(sigma,sigma0,allData) 
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)
			
     case (MULT_BY_J_T_MTX)
        write(*,*) 'Multiplying by J^T (all transmitters)...'
        call JmultT_MTX(sigma0,allData,sigma) 
        header = 'Sensitivity Matrix' 
        call writeAll_Cond2D(fidWrite,cUserDef%wFile_dModelMTX,header,   &
                        allData%nTx,sigma)

     case (INVERSE_NLCG)
        write(*,*) 'Starting the NLCG search...'
        call NLCGsolver(allData,lambda,sigma0,sigma1,alpha)
        call write_Cond2D(fidWrite,cUserDef%wFile_Model,sigma1,TEgrid)
        if (write_data) then
        	call fwdPred(sigma1,allData)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,nSites,sites,allData)
        end if
        
     case default
     
        write(0,*) 'No job ',trim(cUserDef%job),' defined.'
     
     end select

	 ! cleaning up
	 call deallGlobalData()

	 call date_and_time(values=tarray)
	 etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	 rtime = etime - stime
	 write(0,*) ' elapsed time (mins) ',rtime/60.0

end program