
Module MPI_declaration
#ifdef MPI
     implicit none
include 'mpif.h'
! Decleration of general MPI stuff
Integer        :: taskid,numtasks,numworkers,ierr
character*80   :: todo,from_which_operation
Integer        :: file_id
Integer        :: iper, ipol, istn,per_index,pol_index,stn_index,answers_to_receive,received_answers
Integer        :: who, which_stn,which_per,which_pol,nodestate_vec(4),recv_loop 
Integer        :: vector_size
integer        :: size_x,size_y,size_z,nx,ip,nComp,nSite
integer      :: oldtypes(0:20), blockcounts(0:20),offsets(0:20) ,extent 
      integer      block_lengths(0:20)
     integer       displacements(0:20)
    !  integer (kind=MPI_Aint) displacements(0:20)
      integer      address(0:21)
      integer      typelist(0:21)
      
      
Integer        :: cvector_mpi_3D,gridDef3D_mpi,eAll_mpi,dvecMTX_mpi,dvec_mpi,grid_t_mpi,worker_job_task_mpi,modelParam_t_mpi,modelParam_t_mpi_sing,userdef_control_MPI
Integer        :: dvecMTX_mpi_vec(500),dvecMTX_mpi_vec1(500),dvec_mpi_vec(500)
 
      DOUBLE PRECISION      starttime,endtime,time_used,starttime_total,endtime_total
      
      Integer        MASTER, FROM_MASTER, FROM_WORKER,TAG
      INTEGER        STATUS(5)
      parameter      (MASTER = 0)
      parameter      (FROM_MASTER = 1)
      parameter      (FROM_WORKER = 2)
      parameter      (Tag=1)
      Integer        :: dest
       integer nrx
       
 type :: define_worker_job
SEQUENCE
     character*80  :: what_to_do='NOTHING'
     Integer       :: per_index,Stn_index,pol_index,data_type_index
     Integer       :: taskid
     logical       :: keep_E_soln=.false.
     logical       :: several_Tx=.false.
 end type define_worker_job     
      
 type(define_worker_job), save :: worker_job_task 
Integer, pointer, dimension(:) ::ndata_size,nrx_size
logical                        :: eAll_exist=.false.
Integer , pointer, dimension(:) :: eAll_location
Integer                        :: per_index_counter=0
Integer                        :: per_index_vector(100)
integer                        :: prec_MPI
real*8, pointer, dimension(:)        :: model_para_vec

#endif
 end  Module MPI_declaration
  
