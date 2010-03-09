! *****************************************************************************
program Mod3DMT
!  program for running 3D MT forward, sensitivity and inverse modelling

     use senscomp
     use main
     use nlcg
     use DCG
     !use mtinvsetup

#ifdef MPI     
     Use MPI_main 
#endif
     
     implicit none

     ! Character-based information specified by the user
     type (userdef_control)	:: cUserDef

     ! Variables required for storing the date and time
	 real					:: rtime  ! run time
	 real					:: ftime  ! run time per frequency
	 real					:: stime, etime ! start and end times
     integer, dimension(8)	:: tarray ! utility variable

#ifdef MPI 
            call  MPI_constructor 
	        write(6,*)'I am a PARALLEL version'        	        
#else
			write(6,*)'I am a SERIAL version'
#endif


        

     call parseArgs('Mod3DMT',cUserDef) ! OR readStartup(rFile_Startup,cUserDef)

          
	          

    
      call initGlobalData(cUserDef)
     ! set the grid for the numerical computations
      call setGrid(grid)
     
     

#ifdef MPI 
	if (taskid==0) then
	        file_id=2000
	        open(file_id,file='Nodes_Status.info')
	        write(file_id,*) 'Total Number of nodes= ', numworkers
	end if
	        write(6,*)'I am a PARALLEL version'

        if (taskid .ne. 0 ) then
          call Worker_job(sigma0,allData)
            if (trim(worker_job_task%what_to_do) .eq. 'STOPED')  then
               	 call deallGlobalData()
	             call cleanUp()
                 call MPI_destructor
              stop
            end if   
        end if
        	        
#else
			write(6,*)'I am a SERIAL version'
#endif      
     
     
     
#ifdef MPI
       !call Master_job_Distribute_userdef_control(cUserDef) 
       call Master_job_Distribute_Data_Size(allData,sigma0)
       call Master_job_Distribute_Data(allData)
       call Master_job_Distribute_Model(sigma0)       
#endif
     
	 ! Start the (portable) clock
	 call date_and_time(values=tarray)
     stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

     select case (cUserDef%job)

     case (READ_WRITE)
        if (write_model .and. write_data) then
        	write(*,*) 'Writing model and data files and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
        	call write_dataVecMTX(allData,cUserDef%wFile_Data)
		else if (write_model) then
        	write(*,*) 'Writing model and exiting...'
        	call write_modelParam(sigma0,cUserDef%wFile_Model)
		end if

     case (FORWARD)

     write(*,*) 'Calculating predicted data...'

#ifdef MPI   
   if (write_EMsoln) then 
        call Master_job_fwdPred(sigma0,allData,eAll)
        call Master_job_Collect_eAll(allData,eAll)
   else
        call Master_job_fwdPred(sigma0,allData)
   end if
        call Master_job_STOP_MESSAGE
#else          
        call fwdPred(sigma0,allData,eAll)
#endif
        
        if (write_EMsoln) then
        	! write out EM solutions
        	write(*,*) 'Saving the EM solution...'
        	call write_EMsolnMTX(fidWrite,cUserDef%wFile_EMsoln,eAll)
        end if
        ! write out all impedances
        call write_dataVecMTX(allData,cUserDef%wFile_Data)

     case (COMPUTE_J)
        write(*,*) 'Calculating the full sensitivity matrix...'
#ifdef MPI
        !call Master_job_COMPUTE_J(allData,sigma0,sens)
        call Master_job_STOP_MESSAGE 
#else                
        call calcSensMatrix(allData,sigma0,sens)
#endif
        call write_sensMatrixMTX(sens,cUserDef%wFile_Sens)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        
#ifdef MPI        
            !call Master_job_Jmult(dsigma,sigma0,allData)
            call Master_job_STOP_MESSAGE  
#else     
            call Jmult(dsigma,sigma0,allData)
#endif
        
        
        call write_dataVecMTX(allData,cUserDef%wFile_Data)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
#ifdef MPI
         call Master_job_JmultT(sigma0,allData,dsigma)
         call Master_job_STOP_MESSAGE       
#else
         call JmultT(sigma0,allData,dsigma)
#endif       
        
         call write_modelParam(dsigma,cUserDef%wFile_dModel)

     case (INVERSE)
     	if (trim(cUserDef%search) == 'NLCG') then
        	write(*,*) 'Starting the NLCG search...'
        	call NLCGsolver(allData,cUserDef%lambda,sigma0,sigma1,cUserDef%delta)
        
#ifdef MPI        	
        	call Master_job_STOP_MESSAGE
#endif         	
     	elseif (trim(cUserDef%search) == 'DCG') then
        	write(*,*) 'Starting the DCG search...'
        	call DCGsolver(allData,sigma0,sigma1)
            call write_modelParam(sigma1,cUserDef%wFile_Model)
        if (write_data) then
        	call write_dataVecMTX(allData,cUserDef%wFile_Data)
        end if
#ifdef MPI        	
        	call Master_job_STOP_MESSAGE
#endif        	
        	
       else
        	write(*,*) 'Inverse search ',trim(cUserDef%search),' not yet implemented. Exiting...'
        	stop
        end if
        call write_modelParam(sigma1,cUserDef%wFile_Model)
        if (write_data) then
        	call write_dataVecMTX(allData,cUserDef%wFile_Data)
        end if

     case (TEST_COV)
        write(*,*) 'Multiplying input model parameter by covariance ...'
        sigma1 = multBy_CmSqrt(sigma0)
        call write_modelParam(sigma1,cUserDef%wFile_Model)

     case default

        write(0,*) 'No job ',trim(cUserDef%job),' defined.'

     end select

#ifdef MPI       
		close(2000)
#endif

         
	 ! cleaning up
	 call deallGlobalData()
	 call cleanUp()
	 
#ifdef MPI
     	 call date_and_time(values=tarray)
	     etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	     rtime = etime - stime  
	 	 write(0,*) ' elapsed time (mins) ',rtime/60.0
	 call MPI_destructor
#else
	 call date_and_time(values=tarray)
	 etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	 rtime = etime - stime  
	 	 write(0,*) ' elapsed time (mins) ',rtime/60.0	 
#endif
	 


end program
 