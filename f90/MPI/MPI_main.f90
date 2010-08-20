
module MPI_main
#ifdef MPI
  
  use math_constants
  use utilities
  use datasens	 !!!!  inherits : dataspace, dataFunc, SolnSpace
  use SolverSens  !!!  inherits : modelspace, soln2d
  use ForwardSolver
  use sensComp
  use MPI_declaration
  use MPI_sub
  use main

Contains








!###########################################  MPI_initialization   ############################################################

Subroutine MPI_constructor(ctrl)

    implicit none
    type(userdef_control), intent(inout)   :: ctrl
    include 'mpif.h'
     
          call MPI_INIT( ierr )
          call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, ierr )
          numworkers = numtasks-1
          
call create_worker_job_task_mpi
call initUserCtrl(ctrl)
call create_userdef_control_MPI(ctrl)

End Subroutine MPI_constructor




!##############################################   Master_Job_FORWARD #########################################################



Subroutine Master_Job_fwdPred(sigma,d_pred,eAll)
    
    implicit none
    include 'mpif.h'
   type(modelParam_t), intent(in)	    :: sigma   
   type(dataVectorMTX_t), intent(inout)	:: d_pred
   type(solnVectorMTX_t), intent(inout), optional	:: eAll
   integer nTx,ndata,ndt
   logical keep_soln
   type(dataVector_t)                 :: d_temp_TX
   
   keep_soln = present(eAll)

   ! nTX is number of transmitters;
   nTx = d_pred%nTx
   if(associated(eAll_location)) then
             deallocate(eAll_location)
   end if
   allocate(eAll_location(nTx))
   do iper=1,nTx
   eAll_location(iper)=0
   end do
     starttime = MPI_Wtime() 
     
     ! First, distribute the current model to all workers
       call Master_job_Distribute_Model(sigma)
       call Master_job_Distribute_Data(d_pred)
     
               if(.not.d_pred%allocated) then
                  call errStop('data vector not allocated on input to fwdPred')
               end if

               if(present(eAll)) then
                  if(.not. eAll%allocated) then
                     call create_solnVectorMTX(nTx,eAll)
                  else if(d_pred%nTx .ne. eAll%nTx) then
                     call errStop('dimensions of eAll and d do not agree in fwdPred')
                  endif
               endif
               

               
               
                dest=0
                per_index=0
                worker_job_task%what_to_do='FORWARD'
                worker_job_task%several_Tx=.false.
                    do iper=1,d_pred%nTx
                          per_index=per_index+1
                          worker_job_task%per_index= per_index 
                          worker_job_task%keep_E_soln=keep_soln           
                          dest=dest+1                          
                          call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
                          write(file_id,*) 'FWD: Send Per. # ',per_index , 'to node # ',dest    
                          if (dest .ge. numworkers) then
                            goto 10
                          end if 
                          
                    end do

 
10                 continue              
                                   
        answers_to_receive = d_pred%nTx
        received_answers = 0
        

        
        do while (received_answers .lt. answers_to_receive)
                    ! Recv. node's Info.
                     call MPI_RECV(worker_job_task,1,worker_job_task_mpi, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                       
                       who=worker_job_task%taskid
                       which_per=worker_job_task%per_index
                       which_pol=worker_job_task%pol_index
                        
        ! Receive ONLY the predicted data from (who) node

               do ndt=1,d_pred%d(which_per)%ndt
	                  ndata=size(d_pred%d(which_per)%data(ndt)%value)                 
	                  call MPI_RECV(d_pred%d(which_per)%data(ndt)%value(1,1),ndata,MPI_DOUBLE_PRECISION, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
	                  call MPI_RECV(d_pred%d(which_per)%data(ndt)%error(1,1),ndata,MPI_DOUBLE_PRECISION, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
	                  call MPI_RECV(d_pred%d(which_per)%data(ndt)%errorBar,1,MPI_LOGICAL,  who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
               end do
                   eAll_location(which_per)=who

                  write(file_id,*)'FWD: Recv. Resp. for Per # ',which_per ,' from ', who
           
        ! Check if we send all transmitters, if not then send the next transmitter to the node who is free now....
        ! This part is very important if we have less nodes than transmitters.
                       if (Per_index .lt. d_pred%nTx) then
                           Per_index=Per_index+1
                           worker_job_task%per_index= per_index 
                           worker_job_task%what_to_do='FORWARD'
                           worker_job_task%keep_E_soln=keep_soln
                           worker_job_task%several_Tx=.true.
                           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,who, FROM_MASTER, MPI_COMM_WORLD, ierr) 
                           write(file_id,*) 'FWD_CONT: Send Per. # ',per_index , 'to node # ',who 
                       end if 
                       

                       
                       received_answers=received_answers+1        
        end do             
                
        write(file_id,*)'FWD: Finished calculating for (', d_pred%nTx , ') Transmitters '
                
                endtime=MPI_Wtime()
                time_used = endtime-starttime
        write(file_id,*)'FWD: TIME REQUIERED: ',time_used ,'s'

end subroutine Master_Job_fwdPred


!##############################################   Master_Job_Compute_J #########################################################

Subroutine Master_Job_Compute_J(d,sigma,dsigma)
    implicit none
    include 'mpif.h'
    
   type(dataVectorMTX_t), intent(inout)	:: d
   type(modelParam_t), intent(in)	:: sigma
   type(modelParam_t), pointer  :: dsigma(:)
   
   
   integer nTx,nTot,ii,j,row_index,total_stn,jj
   logical savedSolns
     
    
! nTX is number of transmitters;
   nTx = d%nTx
   ! nTot total number of data points
   nTot = countData(d) !d%Ndata 
   
    if(.not.associated(dsigma)) then
      ! allocate for sensitivity matrix
      allocate(dsigma(nTot))
      do j = 1,nTot
         ! this makes a copy of model param, of the same type
         !   as sigma0, then zeros it.
         call copy_ModelParam(dsigma(j),sigma0)
         call zero_ModelParam(dsigma(j))
      enddo
   endif  
   
   
     starttime = MPI_Wtime() 

               

               
               
                dest=0
                per_index=0
                worker_job_task%what_to_do='COMPUTE_J'
                    do iper=1,nTx
                          per_index=per_index+1
                          worker_job_task%per_index= per_index
                          stn_index=0
                        do istn=1,d%d(iper)%data(1)%nSite
                             stn_index=stn_index+1
                             worker_job_task%Stn_index=stn_index         
                             dest=dest+1
                             call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
                             write(file_id,1000) per_index , stn_index, dest    
                          
                          if (dest .ge. numworkers) then
                            goto 11
                          end if 
                       end do
                          
                    end do

 
11                 continue              
    
                                
        answers_to_receive =  nTot
        received_answers = 0
        do while (received_answers .lt. answers_to_receive)
        
                    ! Recv. node's status and id 
                     !call MPI_RECV(nodestate_vec(1),4,MPI_INTEGER, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                     call MPI_RECV(worker_job_task,1,worker_job_task_mpi, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                       
                       who=worker_job_task%taskid
                       which_per=worker_job_task%per_index
                       which_pol=worker_job_task%pol_index
                       which_stn=worker_job_task%stn_index
                        
        ! Receive dsigma(row_index) data from who
      
        
             do ii=1,2
                   if (which_per .gt. 1) then
                    total_stn=0
                    do jj= which_per-1,1,-1
                       total_stn=total_stn+d%d(jj)%data(1)%nSite
                    end do
                    
                    row_index=(total_stn*2) + ((which_stn-1)*2)+ii
                   else
                    row_index=((which_stn-1)*2)+ii
                   end if
                    
                  call create_modelParam_t_mpi(dsigma(row_index))
                  call MPI_RECV(dsigma(row_index),1 ,modelParam_t_mpi, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
                  call MPI_TYPE_FREE (modelParam_t_mpi, IERR) 
                  
 
             end do
                    



                  write(file_id,1001)received_answers,which_per , which_stn, who
         
                  
                  

        ! Check if we send all transmitters, if not then send the next transmitter to the node who is free now....
        !This part is very important if we have less nodes than transmitters.
  
  if (Per_index ==  nTx .and. stn_index == d%d(Per_index)%data(1)%nSite) goto 13 
   
 stn_index=stn_index+1 
 if ( stn_index .gt. d%d(Per_index)%data(1)%nSite ) then
          Per_index=Per_index+1
          stn_index=1 
 elseif ( stn_index .le. d%d(Per_index)%data(1)%nSite) then          
          Per_index=Per_index
 end if
 
 if (Per_index .gt. nTx ) goto 13
  

                                     
                           worker_job_task%per_index= per_index 
                           worker_job_task%Stn_index=stn_index 
                           worker_job_task%what_to_do='COMPUTE_J'
                           
                           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,who, FROM_MASTER, MPI_COMM_WORLD, ierr) 
                           write(file_id,1002) per_index , stn_index, who 
                      
                       
13     continue  
                       
                       received_answers=received_answers+2        
        end do             
 
 
!deallocate(dsigma_matrix_temp1)
                
                
        !DONE: Received soln for all transmitter from all nodes
        write(file_id,*)'COMPUTE_J: Finished calculating for (', d%nTx , ') Transmitters '
                
                endtime=MPI_Wtime()
                time_used = endtime-starttime
        write(file_id,*)'COMPUTE_J: TIME REQUIERED: ',time_used ,'s'
        
  1000 format ('SENS: Send Per # ',i6, ' and  Sta # ', i6 , ' to node # ',i6) 
  1001 format (i6,'SENS: Recv soln for Per # ',i6,' and Sta # ',i6,' from ', i6)
  1002 format( 'SENS_CONT: Send Per # ',i6, ' and  Sta # ', i6 , ' to node # ',i6)

end subroutine Master_Job_Compute_J



!##############################################    Master_job_JmultT #########################################################
Subroutine Master_job_JmultT(sigma,d,dsigma,eAll)

    implicit none
    include 'mpif.h'
    
   type(modelParam_t), intent(in)	:: sigma
   type(dataVectorMTX_t), intent(in)		:: d
   type(modelParam_t), intent(Out)  	:: dsigma
   type(solnVectorMTX_t), intent(in), optional	:: eAll
   type(modelParam_t)                   	:: dsigma_temp
   !type(modelParam_t), dimension(:), pointer 	:: Qcomb_matrix_temp
   type(modelParam_t)                       	:: Qcomb
   integer nTx,nTot,m_dimension,iDT,iTx
   logical savedSolns
   logical		:: calcSomeQ, firstQ
   type(solnVector_t)		:: e,e0
   type(rhsVector_t) 		:: comb
   character(3)         :: iterChar
   character*2          			:: mode
   savedSolns = present(eAll)
  ! nTX is number of transmitters;
   nTx = d%nTx
   ! nTot total number of data points
   nTot = countData(d) !d%Ndata 
   
   starttime = MPI_Wtime() 
 
   !First ditribute both model parameters and data
        call Master_job_Distribute_Model(sigma)
        call Master_job_Distribute_Data(d)
        
   dsigma_temp = sigma
   dsigma 	   = sigma
   call zero(dsigma_temp)
   call zero(dsigma)
                                
                dest=0
                per_index=0
                worker_job_task%what_to_do='JmultT'
                    do iper=1,d%nTx
                          per_index=per_index+1
                          worker_job_task%per_index= per_index 
                         if   (savedSolns) then 
                          dest=eAll_location(per_index)
                         else               
                          dest=dest+1
                         end if
                          
                          call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                   
                          write(file_id,*) 'JmultT: Send Per. # ',per_index , 'to node # ',dest    
                          if (dest .ge. numworkers) then
                            goto 20
                          end if 
                          
                    end do

 
20                 continue              
                                   
        answers_to_receive = d%nTx
        received_answers = 0
        

        
        do while (received_answers .lt. answers_to_receive)
                    ! Recv. node's status and id 
                     !call MPI_RECV(nodestate_vec(1),4,MPI_INTEGER, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                     call MPI_RECV(worker_job_task,1,worker_job_task_mpi, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                       
                       who=worker_job_task%taskid
                       which_per=worker_job_task%per_index
                       which_pol=worker_job_task%pol_index
                        
        ! Receive dsigma_temp from who
                  !call create_modelParam_t_mpi(dsigma_temp)
                  !call MPI_RECV(dsigma_temp,1 ,modelParam_t_mpi, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)  
                  !call linComb_modelParam(ONE,dsigma,ONE,dsigma_temp,dsigma)       
                  !call MPI_TYPE_FREE (modelParam_t_mpi, IERR)  
                  call create_model_param_place_holder(dsigma)
                  m_dimension=size(model_para_vec)
                  call MPI_RECV(model_para_vec(1),m_dimension ,MPI_REAL8, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)  
                  call set_model_para_values(dsigma_temp) 
                  call linComb_modelParam(ONE,dsigma,ONE,dsigma_temp,dsigma)
                     
                              
                  write(file_id,*)'JmultT: Recv. dsigma. for Per # ',which_per ,' from ', who

                 

        ! Check if we send all transmitters, if not then send the next transmitter to the node who is free now....
        !This part is very important if we have less nodes than transmitters.
                       if ( Per_index .lt. d%nTx) then
                           Per_index=Per_index+1
                           worker_job_task%per_index= per_index 
                           worker_job_task%what_to_do='JmultT'
                         if   (savedSolns) then 
                            dest=eAll_location(per_index)
                         else               
                          dest=who
                         end if
                           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)  
                           write(file_id,*) 'JmultT_CONT: Send Per. # ',per_index , 'to node # ',dest 
                       end if
                       

                       
                       received_answers=received_answers+1        
        end do             
    
 
                endtime=MPI_Wtime()
                time_used = endtime-starttime
        !DONE: Received soln for all transmitter from all nodes
        write(file_id,*)'JmultT: Finished calculating for (', d%nTx , ') Transmitters '
        endtime=MPI_Wtime()
        time_used = endtime-starttime
        write(file_id,*)'JmultT: TIME REQUIERED: ',time_used ,'s'



end Subroutine Master_job_JmultT

!##############################################    Master_job_Jmult #########################################################
Subroutine Master_job_Jmult(mHat,m,d,eAll)

    implicit none
    include 'mpif.h'
    
   type(dataVectorMTX_t), intent(inout)		:: d
   type(modelParam_t), intent(in)			:: mHat,m
   type(solnVectorMTX_t), intent(in), optional	:: eAll
   
   integer nTx,nTot,m_dimension,iDT,iTx,ndata,ndt
   logical savedSolns

   savedSolns = present(eAll)
  ! nTX is number of transmitters;
   nTx = d%nTx
   ! nTot total number of data points
   nTot = countData(d) !d%Ndata 
   starttime = MPI_Wtime() 
   
   ! First distribute m, mHat and d
    	    call Master_job_Distribute_Model(m,mHat) 
	        call Master_job_Distribute_Data(d)
	        
	              
                dest=0
                per_index=0
                worker_job_task%what_to_do='Jmult'
                    do iper=1,d%nTx
                          per_index=per_index+1
                          worker_job_task%per_index= per_index 
                         if (savedSolns) then 
                         dest=eAll_location(per_index)
                         else               
                         dest=dest+1
                         end if
                          
                          call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                   
                          write(file_id,*) 'Jmult: Send Per. # ',per_index , 'to node # ',dest    
                          if (dest .ge. numworkers) then
                            goto 20
                          end if 
                          
                    end do

 
20                 continue              
                                   
        answers_to_receive = d%nTx
        received_answers = 0
        

        
        do while (received_answers .lt. answers_to_receive)
                    ! Recv. node's status and id 
                     !call MPI_RECV(nodestate_vec(1),4,MPI_INTEGER, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                     call MPI_RECV(worker_job_task,1,worker_job_task_mpi, MPI_ANY_SOURCE, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)
                       
                       who=worker_job_task%taskid
                       which_per=worker_job_task%per_index
                       which_pol=worker_job_task%pol_index
                  
               do ndt=1,d%d(which_per)%ndt        
                  ndata=size(d%d(which_per)%data(ndt)%value) 
                  call MPI_RECV(d%d(which_per)%data(ndt)%value(1,1),ndata,MPI_DOUBLE_PRECISION, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
                  call MPI_RECV(d%d(which_per)%data(ndt)%error(1,1),ndata,MPI_DOUBLE_PRECISION, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
                  call MPI_RECV(d%d(which_per)%data(ndt)%errorbar,1,MPI_LOGICAL, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)      
              end do

                  write(file_id,*)'Jmult: Recv. data for Per # ',which_per ,' from ', who
     
        
        
        
        ! Check if we send all transmitters, if not then send the next transmitter to the node who is free now....
        !This part is very important if we have less nodes than transmitters.
                       if ( Per_index .lt. d%nTx) then
                           Per_index=Per_index+1
                           worker_job_task%per_index= per_index 
                           worker_job_task%what_to_do='Jmult'
                         if   (savedSolns) then 
                         dest=eAll_location(per_index)
                         else               
                          dest=who
                         end if
                           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)  
                           write(file_id,*) 'Jmult_CONT: Send Per. # ',per_index , 'to node # ',dest 
                       end if
                       

                       
                       received_answers=received_answers+1        
        end do  
                  
                
        !DONE: Received soln for all transmitter from all nodes
        write(file_id,*)'Jmult: Finished calculating for (', d%nTx , ') Transmitters '
                
                endtime=MPI_Wtime()
                time_used = endtime-starttime
        write(file_id,*)'Jmult: TIME REQUIERED: ',time_used ,'s'        
end Subroutine Master_job_Jmult

!############################################## Master_job_Distribute_Data #########################################################
Subroutine Master_job_Distribute_Data(d)
    implicit none
    include 'mpif.h'
    type(dataVectorMTX_t), intent(in)		:: d
    type(dataVectorMTX_t)           		:: d_temp
    integer nTx,nTot,ndata
    DOUBLE PRECISION ,pointer:: buffer(:)
     integer buffer_size,ndt
   ! nTX is number of transmitters;
   nTx = d%nTx
   ! nTot total number of data points
   nTot = countData(d)

        do dest=1,numworkers
            worker_job_task%what_to_do='Distribute Data'
            call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
        end do

call copy_dataVectorMTX(d_temp,d)


  do iper=1,d%nTx
       which_per=iper  
       do ndt=1,d_temp%d(which_per)%ndt
           ndata=size(d_temp%d(which_per)%data(ndt)%value)      
           call MPI_BCAST(d%d(which_per)%data(ndt)%value(1,1),ndata,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)
           call MPI_BCAST(d%d(which_per)%data(ndt)%error(1,1),ndata,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)
           call MPI_BCAST(d%d(which_per)%data(ndt)%errorBar,1,MPI_LOGICAL,0, MPI_COMM_WORLD,ierr)
       end do
           
 end do
 
         

 
  call deall_dataVectorMTX(d_temp)
  
end Subroutine Master_job_Distribute_Data


!############################################## Master_job_Distribute_Model #########################################################
Subroutine Master_job_Distribute_Model(sigma,delSigma)
    implicit none
    include 'mpif.h'
    type(modelParam_t), intent(in) 	:: sigma
    type(modelParam_t), intent(in), optional :: delSigma
    !local
    type(modelParam_t)          	:: sigma_temp
    Integer,pointer:: buffer(:)
    integer buffer_size,m_dimension
    logical send_delSigma
    
 send_delSigma = present(delSigma) 
 
 if (send_delSigma) then
         do dest=1,numworkers
            worker_job_task%what_to_do='Distribute delSigma'
            call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
        end do
         call create_model_param_place_holder(delSigma)
         call get_model_para_values(delSigma)      
         buffer_size=size(model_para_vec)
         call MPI_BCAST(model_para_vec(1),buffer_size,MPI_REAL8,0, MPI_COMM_WORLD,ierr)
        
 
        do dest=1,numworkers
            worker_job_task%what_to_do='Distribute Model'
            call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
        end do
 
         call create_model_param_place_holder(sigma)
         call get_model_para_values(sigma)      
         buffer_size=size(model_para_vec)
         call MPI_BCAST(model_para_vec(1),buffer_size,MPI_REAL8,0, MPI_COMM_WORLD,ierr)

else
        do dest=1,numworkers
            worker_job_task%what_to_do='Distribute Model'
            call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
        end do
 
        call create_model_param_place_holder(sigma)
         call get_model_para_values(sigma)      
        buffer_size=size(model_para_vec)
        call MPI_BCAST(model_para_vec(1),buffer_size,MPI_REAL8,0, MPI_COMM_WORLD,ierr)
end if
  
end Subroutine Master_job_Distribute_Model


!############################################## Master_job_Distribute_eAll #########################################################
Subroutine Master_job_Distribute_eAll(d,eAll)
    implicit none
    include 'mpif.h'
   type(dataVectorMTX_t), intent(in)		:: d
   type(solnVectorMTX_t), intent(in)  	:: eAll
       integer nTx,nTot
       
       
    nTx = d%nTx
    nTot = countData(d) !d%Ndata 
        do dest=1,numworkers
            worker_job_task%what_to_do='Distribute eAll'
            call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
        end do

  do iper=1,d%nTx
       which_per=iper
       do dest=1,numworkers
         call create_eAll_mpi(eAll)
         call MPI_SEND(eAll,1,eAll_mpi,dest, FROM_MASTER,MPI_COMM_WORLD, ierr)
         call MPI_TYPE_FREE (eAll_mpi, IERR)   
        end do
  end do  
end Subroutine Master_job_Distribute_eAll

!############################################## Master_job_Collect_eAll #########################################################
Subroutine Master_job_Collect_eAll(d,eAll)
    implicit none
    include 'mpif.h'
   type(dataVectorMTX_t), intent(in)		:: d
   type(dataVectorMTX_t)            		:: d_local

   type(solnVectorMTX_t)              	:: eAll_local
   type(solnVectorMTX_t)              	:: eAll_temp
   type(solnVectorMTX_t), intent(inout)	:: eAll
   type(solnVector_t)           		:: e0 
   integer nTx,nTot,iTx
       
    nTx = d%nTx

    
          

    
      if(.not. eAll%allocated) then
         call create_solnVectorMTX(d%nTx,eAll)
      else if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in Master_job_Collect_eAll')
      endif
     
      do iTx=1,nTx
         call create_solnVector(grid,iTx,e0)
         call copy_solnVector(eAll%solns(iTx),e0)
      end do
      
            

                      call create_solnVectorMTX(1,eAll_local)
                      call create_solnVector(grid,1,e0)
                      call copy_solnVector(eAll_local%solns(1),e0)

      do iper=1,d%nTx
            worker_job_task%what_to_do='Send eAll to Master'
            worker_job_task%per_index=iper
            who= eAll_location(iper)
            which_per=1
            write(6,*)'Will recv', iper, ' from ', who
            write(6,*) taskid,eAll%solns(1)%grid%nx
            
                   call MPI_SEND(worker_job_task,1,worker_job_task_mpi,who, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
                   call create_eAll_mpi(eAll_local)
                   call MPI_RECV(eAll_local%solns(1),1,eAll_mpi ,who, FROM_WORKER,MPI_COMM_WORLD,STATUS, ierr)
                   call copy_solnVector(eAll%solns(iper),eAll_local%solns(1))
                   call MPI_TYPE_FREE (eAll_mpi, IERR)
        
      end do
      
   
end Subroutine Master_job_Collect_eAll




!############################################## Master_job_Distribute_Data_Size #########################################################
Subroutine Master_job_Distribute_Data_Size(d,sigma0)
    implicit none
    include 'mpif.h'
    
    type(dataVectorMTX_t), intent(in)		:: d
    type(modelParam_t), intent(in) 	:: sigma0
   integer nTx,nTot
   logical savedSolns
     integer subversion
  integer version



     !call create_modelParam_t_mpi(sigma0)
     !modelParam_t_mpi_sing=modelParam_t_mpi

        do dest=1,numworkers
            worker_job_task%what_to_do='Distribute Data_Size'
            call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                     
        end do
        call create_model_param_place_holder(sigma0)


end Subroutine Master_job_Distribute_Data_Size
!############################################## Master_job_Distribute_userdef_control#########################################################
Subroutine Master_job_Distribute_userdef_control(ctrl)
    implicit none
    include 'mpif.h'
    
    type(userdef_control), intent(in)		:: ctrl
       write(6,*)'MASTER: ctrl%wFile_Sens',trim(ctrl%wFile_Sens) 
       write(6,*)'MASTER: ctrl%lambda',(ctrl%lambda)      
       write(6,*)'MASTER: ctrl%eps',(ctrl%eps)    
       write(6,*)'MASTER: ctrl%rFile_Cov',trim(ctrl%rFile_Cov)
       write(6,*)'MASTER: ctrl%search',trim(ctrl%search)
       write(6,*)'MASTER: ctrl%output_level',ctrl%output_level
      
        do dest=1,numworkers
           worker_job_task%what_to_do='Distribute userdef control'
           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)                                   
        end do
        
		!call create_userdef_control_MPI(ctrl)
		call MPI_BCAST(ctrl,1,userdef_control_MPI,0, MPI_COMM_WORLD,ierr)
		 




end Subroutine Master_job_Distribute_userdef_control
	
!##############################################    Master_job_Clean Memory ########################################################
Subroutine Master_job_Clean_Memory

    implicit none
    include 'mpif.h'
    
          write(file_id,*)'Sending Clean memory message to all nodes' 
      
       do dest=1,numworkers
           worker_job_task%what_to_do='Clean memory'
           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)  
           call MPI_RECV(worker_job_task,1,worker_job_task_mpi, dest, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)                                 
           write(file_id,*)'Node       :',worker_job_task%taskid, ' status=  ', trim(worker_job_task%what_to_do)          
                                 
        end do 

end Subroutine Master_job_Clean_Memory

!##############################################    Master_job_Stop_MESSAGE ########################################################
Subroutine Master_job_Stop_MESSAGE 

    implicit none
    include 'mpif.h'
    
          write(file_id,*)'FWD: Sending stop message to all nodes' 
      
       do dest=1,numworkers
           worker_job_task%what_to_do='STOP'
           !call MPI_SEND(per_index,1,MPI_INTEGER,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)  
           call MPI_SEND(worker_job_task,1,worker_job_task_mpi,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)  
           call MPI_RECV(worker_job_task,1,worker_job_task_mpi, dest, FROM_WORKER, MPI_COMM_WORLD, STATUS, ierr)                                 
           !call MPI_RECV(nodestate_vec(1),4,MPI_INTEGER, dest, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)
           write(file_id,*)'Node       :',worker_job_task%taskid, ' status=  ', trim(worker_job_task%what_to_do)          
                                 
        end do 

end Subroutine Master_job_Stop_MESSAGE













!############################################################   Worker_job :High Level Subroutine   #####################################################################
Subroutine Worker_job (sigma0,d)

    implicit none
    include 'mpif.h'
   
   
   
   type(modelParam_t)           	            :: sigma0
   type(modelParam_t)           	            :: sigma_temp,delSigma
   type(modelParam_t)                           :: Qcomb
   type(modelParam_t), pointer, dimension(:), save	:: sigma1
   type(modelParam_t), dimension(:), pointer  	:: dsigma11
    type(modelParam_t)                      	:: dsigma

   type(dataVectorMTX_t) ,intent(inout)    	   :: d
   type(dataVectorMTX_t)                	   :: measu_data
   type(dataVectorMTX_t)                :: d_local
   type(dataVector_t)                 :: d_temp_TX 
      
   type(solnVectorMTX_t)              :: eAll_local,eAll1  
   type(solnVectorMTX_t)              :: eAll,eAll_temp
   type(solnVector_t)           		:: e0  
   type(userdef_control)          :: ctrl
   Integer nTx,m_dimension,ndata,itx,ndt              
     character(80) 		  :: paramType

     real(kind=prec) :: vAir

      !call copy_dataVectorMTX(measu_data ,d)

nTx=d%nTx
recv_loop=0
per_index_vector=0

 do 
          recv_loop=recv_loop+1
          
          !Receive message including what to do and another info. requiered (i.e per_index_stn_index, ...,etc)
          call MPI_RECV(worker_job_task,1,worker_job_task_mpi,0, FROM_MASTER,MPI_COMM_WORLD, STATUS, ierr)

			if (taskid==1) then
			write(6,*) taskid,' TODO: ',trim(worker_job_task%what_to_do),(worker_job_task%keep_E_soln),worker_job_task%several_Tx
			end if
        
if (trim(worker_job_task%what_to_do) .eq. 'FORWARD') then

    if (.NOT. worker_job_task%several_Tx ) then  
         if (eAll1%allocated ) then
            call deall_solnVectorMTX(eAll1)
         end if
         call create_solnVectorMTX(nTx,eAll1)
      do iTx=1,nTx
         call create_solnVector(grid,iTx,e0)
         call copy_solnVector(eAll1%solns(iTx),e0)
      end do
      
    end if
                       per_index=worker_job_task%per_index
                       pol_index=worker_job_task%pol_index
                       worker_job_task%taskid=taskid
                       
                       ! This vector is used later to identify the postions of eAll%soln in eAll
                       !per_index_vector(per_index)=per_index_counter
  

                               

 
                      call  copy_dataVector(d_temp_TX ,measu_data%d(per_index))       
                      call  create_solnVector(grid,1,e0)
         
            ! Do the actual computation
                      call fwdPred_TX(sigma0,d_temp_TX,e0)
            ! Send Info. about the current slave.    
                      call MPI_SEND(worker_job_task,1,worker_job_task_mpi, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 

            ! Send predicted data                    
                   do nDt=1, d_temp_TX%ndt
                       ndata=size(d_temp_TX%data(nDt)%value)       
                       call MPI_SEND(d_temp_TX%data(nDt)%value(1,1),ndata,MPI_DOUBLE_PRECISION,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                       call MPI_SEND(d_temp_TX%data(nDt)%error(1,1),ndata,MPI_DOUBLE_PRECISION,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                       call MPI_SEND(d_temp_TX%data(nDt)%errorBar,1,MPI_LOGICAL,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                   end do
                   
            
            if (worker_job_task%keep_E_soln) then
                      !Keep soln for this period here, later the Master will:
                      ! - Collects eAll for all transmitters, if the used requier and output of eAll (Notice: this step can be avioded if we paralellize the IO stuff)
                      ! - Sends the coressponding  per_index to the worker who has the solution for that transmitter when computing JmutT or Jmult.
                      call copy_solnVector(eAll1%solns(per_index),e0)
                      !eAll%solns(per_index)%tx=per_index
                      eAll_exist=.true.
            end if
            
             
             
elseif (trim(worker_job_task%what_to_do) .eq. 'COMPUTE_J') then
      

                                   
                       per_index=worker_job_task%per_index
                       stn_index=worker_job_task%stn_index
                       worker_job_task%taskid=taskid
      
                      call copy_dataVectorMTX(d_local ,d)        
                      d_local%nTx           = 1
                      !d_local%Ndata         = 2
                      d_local%d(1)%data(1)%nSite    = 1
                      d_local%d(1)%tx       = per_index
                      d_local%d(1)%data(1)%rx(1)    = d%d(per_index)%data(1)%rx(stn_index)
                      d_local%d(1)%data(1)          = d%d(per_index)%data(1)
                      d_local%d(1)%data(1)%dataType = d%d(per_index)%data(1)%dataType
                      d_local%d(1)%data(1)%ncomp    = d%d(per_index)%data(1)%ncomp
                      d_local%d(1)%nDt              =1
                      
            ! Do the actual computation
                     ! call calcSensMatrix(d_local,sigma0,sigma1)
            ! Send Info. about the current slave.    
                      call MPI_SEND(worker_job_task,1,worker_job_task_mpi, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 

            
! Keep a copy of sigma0            
                       call copy_ModelParam(sigma_temp,sigma0)
 
                       call copy_ModelParam(sigma0,sigma1(1)) 
! Send sigma0 to the Master which is now sigma1(1)                        
                       call MPI_SEND(sigma0,1,modelParam_t_mpi_sing,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                       
                       call copy_ModelParam(sigma0,sigma1(2)) 
 ! Send sigma0 to the Master which is now sigma1(2)                        
                       call MPI_SEND(sigma0,1,modelParam_t_mpi_sing,0, FROM_WORKER,MPI_COMM_WORLD, ierr)

 ! Get back sigma0                        
                       call copy_ModelParam(sigma0,sigma_temp)
                  
elseif (trim(worker_job_task%what_to_do) .eq. 'JmultT') then
                                   
   dsigma = sigma0
   call zero(dsigma)
         
                       per_index=worker_job_task%per_index
                       worker_job_task%taskid=taskid
      
              ! Do the actual computation  
              
              if (eAll_exist) then         
                      call JmultT_TX(sigma0,d%d(per_index),dsigma,eAll1%solns(per_index))       
              else
                      call JmultT_TX(sigma0,d%d(per_index),dsigma)      
              end if               

            ! Send Info. about the current worker.    
                      call MPI_SEND(worker_job_task,1,worker_job_task_mpi, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 

                     
                     call create_model_param_place_holder(dsigma)
                     m_dimension=size(model_para_vec)  
                     call get_model_para_values(dsigma)
                     call MPI_SEND(model_para_vec(1),m_dimension,MPI_REAL8,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                     !sigma0 = sigma_temp
!Keep a copy of sigma0
                        !sigma_temp = sigma0
                        !sigma0     = dsigma                  
! Send sigma0, which is now dsigma11, to the Master
                       !call MPI_SEND(sigma0,1,modelParam_t_mpi_sing,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
! Get back sigma0                        
                       !sigma0 = sigma_temp
   per_index_vector=0     
                        
elseif (trim(worker_job_task%what_to_do) .eq. 'Jmult') then
                                   


         
                       per_index=worker_job_task%per_index
                       worker_job_task%taskid=taskid
          
              ! Do the actual computation                      

              if (eAll_exist) then        
                      call Jmult_TX(delSigma,sigma0,d%d(per_index),eAll%solns(per_index))       
              else
                      call Jmult_TX(delSigma,sigma0,d%d(per_index))    
              end if               

                      
            ! Send Info. about the current slave.    
                      call MPI_SEND(worker_job_task,1,worker_job_task_mpi, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 
             ! Send data
             do ndt=1,d%d(per_index)%ndt
                      ndata=size(d%d(per_index)%data(ndt)%value)         
                      call MPI_SEND(d%d(per_index)%data(ndt)%value(1,1),ndata,MPI_DOUBLE_PRECISION,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                      call MPI_SEND(d%d(per_index)%data(ndt)%error(1,1),ndata,MPI_DOUBLE_PRECISION,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                      call MPI_SEND(d%d(per_index)%data(ndt)%errorbar,1,MPI_LOGICAL,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
             end do
                      
   
                                             
elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Data') then
call copy_dataVectorMTX(d_local ,d) 
call copy_dataVectorMTX(d ,d_local) 
              do iper=1,d%nTx
                      which_per=iper
                    do ndt=1,d%d(which_per)%ndt                   
                      ndata=size(d%d(which_per)%data(ndt)%value)
                      call MPI_BCAST(d%d(which_per)%data(ndt)%value(1,1),ndata,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)
                      call MPI_BCAST(d%d(which_per)%data(ndt)%error(1,1),ndata,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)
                      call MPI_BCAST(d%d(which_per)%data(ndt)%errorBar,1,MPI_LOGICAL,0, MPI_COMM_WORLD,ierr)
                    end do   
                 
              end do   
!write(6,*)taskid,'reciev data'
   elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute eAll') then

   call create_solnVectorMTX(nTx,eAll)
              do iper=1,d%nTx
                  which_per=iper
                  call create_eAll_mpi(eAll)
                  call MPI_RECV(eAll,1,eAll_mpi ,0, FROM_MASTER,MPI_COMM_WORLD,STATUS, ierr)
                  call MPI_TYPE_FREE (eAll_mpi, IERR)                   
              end do  
         eAll_exist=.true.
              
elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Model') then
  
!Get the recent model from Master and save it in sigma0
m_dimension=size(model_para_vec)
   !call MPI_BCAST(sigma0,1,modelParam_t_mpi_sing,0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(model_para_vec(1),m_dimension,MPI_REAL8,0, MPI_COMM_WORLD,ierr) 
    call set_model_para_values(sigma0) 
     
elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute delSigma') then
  
!Get delsigma from master
m_dimension=size(model_para_vec)
   !call MPI_BCAST(sigma0,1,modelParam_t_mpi_sing,0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST(model_para_vec(1),m_dimension,MPI_REAL8,0, MPI_COMM_WORLD,ierr) 
    call copy_ModelParam(delSigma,sigma0)
    call set_model_para_values(delSigma)         


elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute userdef control') then
        
        !call initUserCtrl(ctrl)
		!call create_userdef_control_MPI(ctrl)
	  call MPI_BCAST(ctrl,1,userdef_control_MPI,0, MPI_COMM_WORLD,ierr)
	  if (taskid==1 ) then
       write(6,*)'WORKER: ctrl%wFile_Sens',trim(ctrl%wFile_Sens) 
       write(6,*)'WORKER: ctrl%lambda',(ctrl%lambda)      
       write(6,*)'WORKER: ctrl%eps',(ctrl%eps)    
       write(6,*)'WORKER: ctrl%rFile_Cov',trim(ctrl%rFile_Cov)
       write(6,*)'WORKER: ctrl%search',trim(ctrl%search)
       write(6,*)'WORKER: ctrl%output_level',ctrl%output_level
	  end if
	  
      call initGlobalData(ctrl)
      call setGrid(grid)
      call copy_dataVectorMTX(d ,allData)
      call copy_dataVectorMTX(measu_data ,d)

nTx=d%nTx

                        
elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Data_Size') then
   
!Create derived MPI data type called "modelParam_t_mpi", which has the same structures
!as the FORTRAN dirved data type "modelParam_t".
!This operation is called only ONE time at the begining of the prorgam.
!When sending /receiving the model parameter to/from the Master, use modelParam_t_mpi as an MPI derived
!data type. However, the model parameter that must be transfered must be first copied to sigma0, since
!"modelParam_t_mpi" is created using sigma0.  
        call create_model_param_place_holder(sigma0)
        
        !call create_modelParam_t_mpi(sigma0)
        !modelParam_t_mpi_sing=modelParam_t_mpi      

                                                
elseif (trim(worker_job_task%what_to_do) .eq. 'Send eAll to Master' ) then 



                       per_index=worker_job_task%per_index
                       worker_job_task%taskid=taskid
      

                      call create_solnVectorMTX(1,eAll_local)
                      call create_solnVector(grid,1,e0)
                      call copy_solnVector(eAll_local%solns(1),e0)
                      

                      which_per=1
                      call create_eAll_mpi(eAll_local)
                      call copy_solnVector(eAll_local%solns(1),eAll1%solns(per_index))                     
                      call MPI_SEND(eAll_local%solns(which_per),1,eAll_mpi,0, FROM_WORKER,MPI_COMM_WORLD, ierr)
                      call MPI_TYPE_FREE (eAll_mpi, IERR)
                       
elseif (trim(worker_job_task%what_to_do) .eq. 'Clean memory' ) then   
         
               
         call deall_modelParam(dsigma11(1))          
         deallocate (dsigma11)
         call deall_modelParam(sigma_temp)
         
         call deall_dataVectorMTX(d_local)
         call deall_dataVectorMTX(d_local)
         
         !call deall_solnVectorMTX(eAll_local1)    
         call deall_solnVectorMTX(eAll)
         call deall_solnVectorMTX(eAll_local)
         
         worker_job_task%what_to_do='Cleaned Memory and Waiting'
         worker_job_task%taskid=taskid
         call MPI_SEND(worker_job_task,1,worker_job_task_mpi, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 

                                   
elseif (trim(worker_job_task%what_to_do) .eq. 'STOP' ) then

                             worker_job_task%what_to_do='STOPED'
                             worker_job_task%taskid=taskid
                             call MPI_SEND(worker_job_task,1,worker_job_task_mpi, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 
                             !call MPI_TYPE_FREE (modelParam_t_mpi_sing, IERR)
                             exit

end if

if (taskid==1) then
write(6,*) taskid,' Finish: ',trim(worker_job_task%what_to_do)
end if

end do


End Subroutine Worker_job



subroutine create_dvecMTX_mpi(d)

    implicit none

    !include 'mpif.h'    
     type(dataVectorMTX_t), intent(in)	:: d
     integer ii,I1,extent_int,extent_real,extent_logic,extent,ndata
     
        


            ndata=size(d%d(which_per)%data(1)%value)
            nrx=size(d%d(which_per)%data(1)%rx)

    
    !goto 10

      typelist(0) = MPI_INTEGER
      !typelist(1) = MPI_INTEGER
      typelist(1) = MPI_INTEGER
      !typelist(3) = MPI_INTEGER
      typelist(2) = MPI_DOUBLE_PRECISION
      typelist(3) = MPI_DOUBLE_PRECISION
      typelist(4) = MPI_INTEGER
      typelist(5) = MPI_INTEGER
      typelist(6) = MPI_LOGICAL
      !typelist(9) = MPI_LOGICAL
      
      block_lengths(0) =2  
      !block_lengths(1) =1
      block_lengths(1) =2  
      !block_lengths(3) =1
      block_lengths(2) = ndata
      block_lengths(3) = ndata
      block_lengths(4) = nrx
      block_lengths(5) = 2
      block_lengths(6) = 2
      !block_lengths(9) = 1

    
      call MPI_Address(d%nTX,                                          address(0), ierr)
      !call MPI_Address(d%Ndata,                                        address(1), ierr)
      call MPI_Address(d%d(which_per),                                 address(1), ierr)
      call MPI_Address(d%d(which_per)%data(1)%nComp,                           address(2), ierr)
      !call MPI_Address(d%d(which_per)%nSite,                           address(4), ierr)
      call MPI_Address(d%d(which_per)%data(1)%value(1,1),                       address(3), ierr)
      call MPI_Address(d%d(which_per)%data(1)%error(1,1),                        address(4), ierr)
      call MPI_Address(d%d(which_per)%data(1)%rx(1),                           address(5), ierr)
      call MPI_Address(d%d(which_per)%tx,                              address(6), ierr)
      call MPI_Address(d%d(which_per)%allocated,                       address(7), ierr)
      !call MPI_Address(d%d(which_per)%errorBar,                        address(10), ierr)


    !write(6,*) taskid,ndata,nrx,which_per,address(0),address(1),address(2),address(3)

     do ii=0,6
       displacements(ii) = address(ii+1) - address(0)
     end do   
     
      call MPI_TYPE_STRUCT(7, block_lengths, displacements,typelist, dvecMTX_mpi, ierr)   
      call MPI_TYPE_COMMIT(dvecMTX_mpi, ierr)
  
 10 continue     
 
  

end subroutine create_dvecMTX_mpi

subroutine create_dvec_mpi(d)
    implicit none

    !include 'mpif.h'    
     type(dataVectorMTX_t), intent(in)	:: d
     integer ii,I1,extent_int,extent_real,extent_logic,ndata
     
        


            ndata=size(d%d(which_per)%data(1)%value)
            nrx=size(d%d(which_per)%data(1)%rx)
            
            
   offsets(0) = 0
   oldtypes(0) = MPI_INTEGER
   blockcounts(0) = 2

   call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
   offsets(1) = 2 * extent
   oldtypes(1) = MPI_DOUBLE_PRECISION
   blockcounts(1) = ndata

   call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
   offsets(2) = ndata* extent
   oldtypes(2) = MPI_DOUBLE_PRECISION
   blockcounts(2) = ndata  
   
      
   call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
   offsets(3) = ndata* extent
   oldtypes(3) = MPI_INTEGER
   blockcounts(3) = nrx  
   
   call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
   offsets(4) = nrx * extent
   oldtypes(4) = MPI_INTEGER
   blockcounts(4) = 2
   
   call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
   offsets(5) = 2 *extent
   oldtypes(5) = MPI_LOGICAL
   blockcounts(5) = 2     
   
   call MPI_TYPE_STRUCT(6, blockcounts, offsets, oldtypes,dvec_mpi, ierr)
   call MPI_TYPE_COMMIT(dvec_mpi, ierr)


end subroutine create_dvec_mpi


subroutine create_worker_job_task_mpi


    implicit none
    include 'mpif.h'
integer :: ii,intex,CHARACTERex
integer(MPI_ADDRESS_KIND) address1(0:20)


   offsets(0) = 0
   oldtypes(0) = MPI_CHARACTER
   blockcounts(0) =  80

   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(1) = (80* extent)+offsets(0)
   oldtypes(1) = MPI_INTEGER
   blockcounts(1) = 5
   		
   call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
   offsets(2) = (5* extent)+offsets(1)
   oldtypes(2) = MPI_LOGICAL
   blockcounts(2) = 2
   
       call MPI_TYPE_STRUCT(3, blockcounts, offsets, oldtypes,worker_job_task_mpi, ierr)
       call MPI_TYPE_COMMIT(worker_job_task_mpi, ierr)
      
   
     

      typelist(0) = MPI_CHARACTER
      typelist(1) = MPI_INTEGER
      typelist(2) = MPI_LOGICAL


      
      block_lengths(0) = 80  
      block_lengths(1) = 5
      block_lengths(2) = 2


      call MPI_Address(worker_job_task,                   address1(0), ierr)
      call MPI_Address(worker_job_task%what_to_do,        address1(1), ierr)
      call MPI_Address(worker_job_task%per_index,         address1(2), ierr)
      call MPI_Address(worker_job_task%keep_E_soln,        address1(3), ierr)
      

     do ii=0,2
       displacements(ii) = address1(ii+1) - address1(0)
     end do   
     
      !call MPI_TYPE_STRUCT(3, block_lengths, displacements,typelist, worker_job_task_mpi, ierr)   
      !call MPI_TYPE_COMMIT(worker_job_task_mpi, ierr) 
     


end subroutine create_worker_job_task_mpi

!##########################################################################
!Creating an MPI version of userdef_control_MPI 
subroutine create_userdef_control_MPI(ctrl) 
    implicit none
    include 'mpif.h'
    
	type(userdef_control), intent(in)   :: ctrl
	integer :: ii,intex,CHARACTERex,sum1
	integer(MPI_ADDRESS_KIND) address1(0:20)

	
   offsets(0) = 0
   oldtypes(0) = MPI_CHARACTER
   blockcounts(0) =  80*13

   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(1) = (80*13 * extent)+offsets(0)
   oldtypes(1) = MPI_REAL8
   blockcounts(1) = 2
   		
   call MPI_TYPE_EXTENT(MPI_REAL8, extent, ierr)
   offsets(2) = (2* extent)+offsets(1)
   oldtypes(2) = MPI_CHARACTER
   blockcounts(2) = 80*2
   
   call MPI_TYPE_EXTENT(MPI_CHARACTER, extent, ierr)
   offsets(3) = (80*2* extent)+offsets(2)
   oldtypes(3) = MPI_INTEGER
   blockcounts(3) = 1
   
     call MPI_TYPE_STRUCT(4, blockcounts, offsets, oldtypes,userdef_control_MPI, ierr)
     call MPI_TYPE_COMMIT(userdef_control_MPI, ierr)
    

      
   
      typelist(0) = MPI_CHARACTER                                                 
      typelist(1) = MPI_REAL8
      typelist(2) = MPI_CHARACTER
      typelist(3) = MPI_INTEGER

      
      block_lengths(0) = 80*13 
      block_lengths(1) = 2  
      block_lengths(2) = 80*2
      block_lengths(3) = 1
      
      call MPI_Address(ctrl,                	address1(0), ierr)
      call MPI_Address(ctrl%job,             	address1(1), ierr)
      call MPI_Address(ctrl%lambda,          	address1(2), ierr)
      call MPI_Address(ctrl%rFile_Cov,        	address1(3), ierr)
      call MPI_Address(ctrl%output_level,       address1(4), ierr)


     do ii=0,3
       displacements(ii) = address1(ii+1) - address1(0)
     end do   
     
     !call MPI_TYPE_STRUCT(4, block_lengths, displacements,typelist, userdef_control_MPI, ierr)   
     ! call MPI_TYPE_COMMIT(userdef_control_MPI, ierr) 






end  subroutine create_userdef_control_MPI





subroutine MPI_destructor
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)       
      call MPI_FINALIZE(ierr)

end subroutine MPI_destructor
#endif    
end module MPI_main
 
