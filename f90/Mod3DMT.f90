! *****************************************************************************
program Mod3DMT
!  program for running 3D MT forward, sensitivity and inverse modelling

     use sensmatrix
     use main
     use nlcg
     !use mtinvsetup

     implicit none

     ! Character-based information specified by the user
     type (userdef_control)	:: cUserDef

     ! Variables required for storing the date and time
	 real					:: rtime  ! run time
	 real					:: ftime  ! run time per frequency
	 real					:: stime, etime ! start and end times
     integer, dimension(8)	:: tarray ! utility variable

     integer                :: nData

     call parseArgs('Mod3DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)

     call initGlobalData(cUserDef)

	 ! Start the (portable) clock
	 call date_and_time(values=tarray)
     stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_modelParam(fidWrite,cUserDef%wFile_Model,sigma0)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,nSites,sites,siteids,data_units,compids,allData)
		else if (write_model) then
        	write(*,*) 'Writing model and exiting...'
        	call write_modelParam(fidWrite,cUserDef%wFile_Model,sigma0)
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
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,nSites,sites,siteids,data_units,compids,allData)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
        nData = countData(allData)
        call writeVec_modelParam(fidWrite,cUserDef%wFile_Sens,   &
                        nData,sigma,'Sensitivity matrix')

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call Jmult(dsigma,sigma0,allData)
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,nSites,sites,siteids,data_units,compids,allData)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma)
        call write_modelParam(fidWrite,cUserDef%wFile_dModel,dsigma)

     case (INVERSE)
     	if (trim(cUserDef%search) == 'NLCG') then
        	write(*,*) 'Starting the NLCG search...'
        	call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,cUserDef%delta)
        else
        	write(*,*) 'Inverse search ',trim(cUserDef%search),' not yet implemented. Exiting...'
        	stop
        end if
        call write_modelParam(fidWrite,cUserDef%wFile_Model,sigma1)
        if (write_data) then
        	call fwdPred(sigma1,allData)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,nSites,sites,siteids,data_units,compids,allData)
        end if

     case (TEST_COV)
        write(*,*) 'Multiplying input model parameter by covariance ...'
        call multBy_CmSqrt(sigma0,sigma1)
        call write_modelParam(fidWrite,cUserDef%wFile_Model,sigma1)

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
