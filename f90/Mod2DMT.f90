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
!     use modelspace
!     use dataspace
!     use solnrhs
!     use wsfwd2d
!     use ioascii
!     use dataFunc
     use sensmatrix
     use main
     use nlcg

     implicit none

     ! Character-based information specified by the user
     type (userdef_control)	:: cUserDef

     ! Variables required for storing the date and time
	 real					:: rtime  ! run time
	 real					:: ftime  ! run time per frequency
	 real					:: stime, etime ! start and end times
     integer, dimension(8)	:: tarray ! utility variable

     integer                :: nData
!     integer (kind=4) :: Nzb, IER, i, iy,iz,iPer,nSigma,nTx
!     integer (kind=4) :: iargc,narg,k
!     integer (kind=4) :: nMode=1,nComp=2
!     character*80  gridType, header,arg, paramtype
!     character*80, dimension(:), pointer :: temp

     call parseArgs('Mod2DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)

     call initGlobalData(cUserDef)

     ! set the grid for the numerical computations
     call setGrid(grid)

	 ! Start the (portable) clock
	 call date_and_time(values=tarray)
     stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
				nSites,sites,data_units,allData)
		else if (write_model) then
        	write(*,*) 'Writing model and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
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
			nSites,sites,data_units,allData)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
        nData = countData(allData)
        call writeVec_modelParam(nData,sigma,'Sensitivity matrix',cUserDef%wFile_Sens)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call Jmult(dsigma,sigma0,allData)
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,data_units,allData)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma)
        call write_modelParam(dsigma,cUserDef%wFile_dModel)

     case (INVERSE)
     	if (trim(cUserDef%search) == 'NLCG') then
        	write(*,*) 'Starting the NLCG search...'
        	call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,cUserDef%delta)
        else
        	write(*,*) 'Inverse search ',trim(cUserDef%search),' not yet implemented. Exiting...'
        	stop
        end if
        call write_modelParam(sigma1,cUserDef%wFile_Model)
        if (write_data) then
        	call fwdPred(sigma1,allData)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,nSites,sites,data_units,allData)
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
