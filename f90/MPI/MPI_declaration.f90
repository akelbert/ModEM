
Module MPI_declaration
#ifdef MPI
     implicit none
include 'mpif.h'

! Decleration of general MPI stuff
!******************************************************************** 
Integer        :: taskid,toatl_number_of_Proc,number_of_workers
Integer        :: file_id
Integer        :: MASTER, FROM_MASTER, FROM_WORKER,TAG,ierr,dest
INTEGER        :: STATUS(5)
parameter         (MASTER=0,FROM_MASTER=1,FROM_WORKER=2,Tag=1)
!******************************************************************** 




      
! Parameters required to create an MPI derived data types.
!********************************************************************      
Integer        :: cvector_mpi_3D,gridDef3D_mpi,eAll_mpi,dvecMTX_mpi 
Integer        :: dvec_mpi,grid_t_mpi,worker_job_task_mpi           
Integer        :: modelParam_t_mpi_sing,userdef_control_MPI
Integer        :: modelParam_t_mpi,extent
integer        :: oldtypes(0:20), blockcounts(0:20),offsets(0:20)  
integer        :: block_lengths(0:20)
integer        :: displacements(0:20)
integer        :: address(0:21)
integer        :: typelist(0:21)
!******************************************************************** 




        
! Parameters used in communication
!********************************************************************   
Integer        :: answers_to_receive,received_answers,recv_loop 
Integer        :: who, which_stn,which_per,which_dt,which_pol 
Integer , pointer, dimension(:)  :: eAll_location
logical                          :: eAll_exist=.false.
real*8,   pointer, dimension(:)  :: model_para_vec
character, pointer, dimension(:) :: eAll_para_vec   !! needed for MPI_pack/MPI_unpack; counted in bytes
Integer                          :: Nbytes
!********************************************************************     




! Time measuring   
!********************************************************************      
DOUBLE PRECISION    :: starttime,endtime,time_used
DOUBLE PRECISION    :: starttime_total,endtime_total 
!********************************************************************





! A drived data type to distribute Info. between Processors.
!********************************************************************        
type :: define_worker_job
     SEQUENCE
     character*80  :: what_to_do='NOTHING'
     Integer       :: per_index,Stn_index,pol_index,data_type_index
     Integer       :: taskid
     logical       :: keep_E_soln=.false.
     logical       :: several_Tx=.false.
 end type define_worker_job           
type(define_worker_job), save :: worker_job_task
!********************************************************************  


#endif
 end  Module MPI_declaration
  
