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

     use SensComp
     use SymmetryTest
     use Main
     use NLCG

     implicit none

     ! Character-based information specified by the user
     type (userdef_control)	:: cUserDef

     ! Variable required for storing the date and time
     type (timer_t)         :: timer

     integer                :: nData
!     integer (kind=4) :: Nzb, IER, i, iy,iz,iPer,nSigma,nTx
!     integer (kind=4) :: iargc,narg,k
!     integer (kind=4) :: nMode=1,nComp=2
!     character*80  gridType, header,arg, paramtype
!     character*80, dimension(:), pointer :: temp

     ! needed for symmetry testing
     type (solnVectorMTX_t) :: ePred
     type (rhsVectorMTX_t)  :: bAll

     call parseArgs('Mod2DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)

     call initGlobalData(cUserDef)

     ! set the grid for the numerical computations
     call setGrid(grid)

	 ! Start the (portable) clock
	 call reset_time(timer)

     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
        	call write_dataVectorMTX(allData,cUserDef%wFile_Data)
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
        	call write_solnVectorMTX(fidWrite,cUserDef%wFile_EMsoln,eAll)
        end if
        ! write out all impedances
        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcJ(allData,sigma0,sens)
        call write_sensMatrixMTX(sens,cUserDef%wFile_Sens)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call Jmult(dsigma,sigma0,allData)
        call write_dataVectorMTX(allData,cUserDef%wFile_Data)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma)
        call write_modelParam(dsigma,cUserDef%wFile_dModel)

     case (INVERSE)
     	if (trim(cUserDef%search) == 'NLCG') then
            ! sigma1 contains mHat on input (zero = starting from the prior)
        	write(*,*) 'Starting the NLCG search...'
        	sigma1 = dsigma
        	call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,cUserDef%rFile_invCtrl)
        else
        	write(*,*) 'Inverse search ',trim(cUserDef%search),' not yet implemented. Exiting...'
        	stop
        end if
        call write_modelParam(sigma1,cUserDef%wFile_Model)
        if (write_data) then
        	call write_dataVectorMTX(allData,cUserDef%wFile_Data)
        end if

     case (TEST_ADJ)
       select case (cUserDef%test)
           case('J')
               call Jtest(sigma0,dsigma,allData)
           case('L')
               !call fwdPred(sigma0,allData,eAll)
			   !ePred = eAll
			   call random_solnVectorMTX(eAll,cUserDef%delta)
			   call random_dataVectorMTX(allData,cUserDef%delta)
			   !call Ltest(sigma0,eAll,allData,ePred)
			   !call deall_solnVectorMTX(ePred)
               call Ltest(sigma0,eAll,allData)
           case('S')
               call random_solnVectorMTX(eAll,cUserDef%delta)
               bAll = eAll
               call Stest(sigma0,bAll,eAll)
           case('Q')
               call Qtest(sigma0,dsigma,allData)

           case('m')
               call random_modelParam(sigma0,cUserDef%delta)
               call write_modelParam(sigma0,cUserDef%wFile_Model)
               write_model = .false.
               write_data = .false.
           case('d')
               call random_dataVectorMTX(allData,cUserDef%delta)
               call write_dataVectorMTX(allData,cUserDef%wFile_Data)
               write_model = .false.
               write_data = .false.
           case('e')
               call random_solnVectorMTX(eAll,cUserDef%delta)
               call write_solnVectorMTX(ioWrite,cUserDef%wFile_EMsoln,eAll)
               write_model = .false.
               write_data = .false.

           case default
               write(0,*) 'Symmetry test for operator ',trim(cUserDef%test),' not yet implemented.'
       end select
	     if (write_model .and. write_data) then
	        write(*,*) 'Writing model and data files and exiting...'
	        call write_modelParam(dsigma,cUserDef%wFile_Model)
	        call write_dataVectorMTX(allData,cUserDef%wFile_Data)
	     else if (write_model) then
	        write(*,*) 'Writing model and exiting...'
	        call write_modelParam(dsigma,cUserDef%wFile_Model)
	     end if

     case default

        write(0,*) 'No job ',trim(cUserDef%job),' defined.'

     end select

	 ! cleaning up
	 call deallGlobalData()

	 call cleanUp()

	 write(0,*) ' elapsed time (mins) ',elapsed_time(timer)/60.0

end program
