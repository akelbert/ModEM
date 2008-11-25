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
!     use modelparameter
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

!     integer (kind=4) :: Nzb, IER, i, iy,iz,iPer,nSigma,nTx
!     integer (kind=4) :: iargc,narg,k
!     integer (kind=4) :: nMode=1,nComp=2
!     character*80  gridType, header,arg, paramtype
!     character*80, dimension(:), pointer :: temp

     call parseArgs('Mod2DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)

     call initGlobalData(cUserDef)

	 ! Start the (portable) clock
	 call date_and_time(values=tarray)
     stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_modelParam(fidWrite,cUserDef%wFile_Model,sigma0)
        	call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
				nSites,sites,allData)
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
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
        call writeVec_modelParam(fidWrite,cUserDef%wFile_Sens,   &
                        allData%nData,sigma,'Sensitivity matrix')

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call Jmult(dsigma,sigma0,allData)
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma)
        call write_modelParam(fidWrite,cUserDef%wFile_dModel,dsigma)

     case (MULT_BY_J_MTX)
        write(*,*) 'Multiplying by J (all transmitters)...'
        call Jmult_MTX(sigma,sigma0,allData)
        call write_Z(fidWrite,cUserDef%wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (MULT_BY_J_T_MTX)
        write(*,*) 'Multiplying by J^T (all transmitters)...'
        call JmultT_MTX(sigma0,allData,sigma)
        call writeVec_modelParam(fidWrite,cUserDef%wFile_dModelMTX,   &
                        allData%nTx,sigma,'J^T x d (all transmitters)')

     case (INVERSE_NLCG)
        write(*,*) 'Starting the NLCG search...'
        call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,cUserDef%delta)
        call write_modelParam(fidWrite,cUserDef%wFile_Model,sigma1)
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
