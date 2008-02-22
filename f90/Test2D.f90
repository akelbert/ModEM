! *****************************************************************************
program Test2D
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
     use iobinary
     use dataFunc
     use sensmatrix
     use nlcg
        
     implicit none
	
     ! I/O units ... reuse generic read/write units if 
     !   possible; for those kept open during program run, 
     !   reserve a specific unit here 
     integer (kind=4) :: fidRead = 1
     integer (kind=4) :: fidWrite = 2
     integer (kind=4) :: fidError = 99

     ! characters used to comunicate (from matlab)
     !   about which jobs is executed
     character*1		:: job = ' '
     character*1, parameter	:: FORWARD_PRED = 'F'
     character*1, parameter	:: FORWARD_SOLN = 'E'
     character*1, parameter	:: SENSITIVITIES = 'S'
     character*1, parameter	:: MULT_BY_J = 'G'
     character*1, parameter	:: MULT_BY_J_T = 'T'
     character*1, parameter	:: MULT_BY_J_MTX = 'N'
     character*1, parameter	:: MULT_BY_J_T_MTX = 'M'
     character*1, parameter	:: FORWARD_PRED_DELTA = 'D'
     character*1, parameter	:: CREATE_DATA_ERRORS = 'C'
     character*1, parameter	:: CREATE_BG_MODEL = 'B'
     character*1, parameter	:: NLCG_INVERSION = 'I'

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
        
     real (kind=selectedPrec), dimension(:), pointer	:: periods
     real (kind=selectedPrec), dimension(:,:), pointer	:: sites
     character*2, dimension(:), pointer     		:: modes

     ! grid geometry data structure
     type(grid2d_t), target 	:: TEgrid

     ! impedance data structure
     type(dvecMTX)		:: allData

     !  storage for the "background" conductivity parameter
     type(modelParam_t)		:: sigma0
     !  storage for the full sensitivity matrix
     type(modelParam_t),dimension(:),pointer	:: sigma
     !  storage for a perturbation to conductivity
     type(modelParam_t)		:: dsigma
     !  storage for the inverse solution
     type(modelParam_t)		:: sigma1

     !  storage for EM solutions
     type(EMsolnMTX)            :: eAll

     !  the damping parameter
     real(kind=selectedPrec)    :: lambda = 1.0
     real(kind=selectedPrec)    :: alpha = 0.1

     real(kind=selectedPrec)	:: eps
     real(kind=selectedPrec)	:: error
     
     real(kind=selectedPrec)	:: bg_cond_value
     real(kind=selectedPrec),dimension(:,:), pointer  :: bg_cond

     integer (kind=4) :: Nzb, IER, nPer, i, iy,iz,iPer,nSigma,nTx
     integer (kind=4) :: iargc,narg,k
     integer (kind=4) :: nMode=1,nComp=2,nSites
     character*80  gridType, header,arg, paramtype
     character*80, dimension(:), pointer :: temp

     !  parse command line ...  for now just uses first argument
     !   to set job
     narg = iargc()
     if(narg .gt. 0) then
        call getarg(1,arg)
        if(arg(1:1).eq.'-') job = arg(2:2)
     else
        write(*,*) 'Usage: Test2D -[option] [args]'
        write(*,*)
        write(*,*) ' Available options:'
        write(*,*) '-F  calculates the predicted data'
        write(*,*) '-E  calculates the predicted data and saves the EM solution'
        write(*,*) '-S  calculates and saves the full sensitivity matrix'
        write(*,*) '-G  multiplies a model by J to create a data vector'
        write(*,*) '-T  multiplies a data vector by J^T to create a model'
        write(*,*) '-N  evaluates sum( J m_i ) over transmitters to yield a data'
        write(*,*) '    vector; reads a sensitivity matrix from file to do this'
        write(*,*) '-M  multiples d_i by J_i^T separately for each transmitter, '
        write(*,*) '    to yield a bunch of models, one for each transmitter'
        write(*,*) '-D  calculated the predicted data for a perturbed model'
        !write(*,*) '-C  creates data errors in a crude way (inversion testing)'
        !write(*,*) '-B  creates background model for inversion from the grid'
        write(*,*) '-I  runs an NLCG inversion to yield an inverse model'
        write(*,*)
        write(*,*) ' Additional arguments:'
        write(*,*) '-F  rFile_Model rFile_Data wFile_Data'
        write(*,*) '-E  rFile_Model rFile_Data wFile_Data wFile_EMsoln'
        write(*,*) '-S  rFile_Model rFile_Data wFile_Sens'
        write(*,*) '-G  rFile_Model rFile_dModel rFile_Data wFile_Data'
        write(*,*) '-T  rFile_Model rFile_Data wFile_dModel'
        write(*,*) '-N  rFile_Model rFile_Data wFile_dModelMTX'
        write(*,*) '-M  rFile_Model rFile_dModelMTX rFile_Data wFile_Data'
        write(*,*) '-D  rFile_Model rFile_dModel rFile_Data wFile_Data'
        !write(*,*) '-C  rFile_Data wFile_Data'
        !write(*,*) '-B  rFile_Model wFile_Model bg_cond_value'
        write(*,*) '-I  rFile_Model rFile_Data wFile_Model wFile_Data lambda alpha'
        stop
     endif
     
     ! extract all following command line arguments
     allocate(temp(1:narg-1))
     do k = 1,narg-1
       call getarg(k+1,temp(k))
     end do

     ! save model file name in all cases
     
     select case (job)

      case (FORWARD_PRED, FORWARD_SOLN)
        ! F,E
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_Data = temp(2)
	    end if
	    if (narg > 3) then
	       wFile_Data = temp(3)
	    end if
	    if (narg > 4) then
	       wFile_EMsoln = temp(4)
	    end if

      case (FORWARD_PRED_DELTA)
        ! D
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
        if (narg > 2) then
	       rFile_dModel = temp(2)
	    end if
	    if (narg > 3) then
	       rFile_Data = temp(3)
	    end if
	    if (narg > 4) then
	       wFile_Data = temp(4)
	    end if
	    if (narg > 5) then
	       wFile_EMsoln = temp(5)
	    end if

      case (SENSITIVITIES) 
        ! S
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_Data = temp(2)
	    end if
	    if (narg > 3) then
	       wFile_Sens = temp(3)
	    end if

      case (MULT_BY_J) 
        ! G
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_dModel = temp(2)
	    end if
	    if (narg > 3) then
	       rFile_Data = temp(3)
	    end if
	    if (narg > 4) then
	       wFile_Data = temp(4)
	    end if

      case (MULT_BY_J_T)
        ! T
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_Data = temp(2)
	    end if
	    if (narg > 3) then
	       wFile_dModel = temp(3)
	    end if

      case (MULT_BY_J_T_MTX)
        ! N
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_Data = temp(2)
	    end if
	    if (narg > 3) then
	       wFile_dModelMTX = temp(3)
	    end if

      case (MULT_BY_J_MTX)
        ! M
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_dModelMTX = temp(2)
	    end if
	    if (narg > 3) then
	       rFile_Data = temp(3)
	    end if
	    if (narg > 4) then
	       wFile_Data = temp(4)
	    end if
			
      case (CREATE_DATA_ERRORS)
        ! C
	    if (narg > 1) then
	       rFile_Data = temp(1)
	    end if
	    if (narg > 2) then
	       wFile_Data = temp(2)
	    end if

      case (CREATE_BG_MODEL)
        ! B
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       wFile_Model = temp(2)
	    end if
        if (narg > 3) then
           read(temp(3),*) bg_cond_value
        else
           bg_cond_value = 1e-2
        end if
 	          
      case (NLCG_INVERSION)
        ! I
        if (narg > 1) then
	       rFile_Model = temp(1)
	    end if
	    if (narg > 2) then
	       rFile_Data = temp(2)
	    end if
	    if (narg > 3) then
	       wFile_Model = temp(3)
	    end if
	    if (narg > 4) then
	       wFile_Data = temp(4)
	    end if
	    if (narg > 5) then
          read(temp(5),*) lambda
        end if
        if (narg > 6) then
          read(temp(6),*) alpha
        end if
        
      case default
         call errStop('Unknown job. Please check your command line options')
      
     end select
     
     ! Read background conductivity parameter and grid
     if (len_trim(rFile_Model)>0) then
       call read_Cond2D(fidRead,rFile_Model,sigma0,TEgrid)
     
       !  set array size parameters in WS forward code module 
       !   these stay fixed for all forward modeling with this grid
       call setWSparams(TEgrid%Ny,TEgrid%Nz,TEgrid%Nza)

       !  set grid for higher level solver routines
       call set_SolnRHS_grid(TEgrid)
     else
       call warning('No input model parametrization')
     end if
	
     ! (3) Read in data file (only a template on input--periods/sites)
     if (len_trim(rFile_Data)>0) then
       call read_Z(fidRead,rFile_Data,nPer,periods,modes,nSites,sites,allData)
       !  Using periods, sites obtained from data file
       !     set up transmitter and receiver dictionaries
       call TXdictSetUp(nPer,periods,modes) 
       call RXdictSetUp(nSites,sites)
       call TypeDictSetup()
     else
       call warning('No input data file - unable to set up dictionaries')
     end if

	 ! Start the (portable) clock
	 call date_and_time(values=tarray)
     stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

     select case (job)

!     case (CREATE_BG_MODEL)
!        write(*,*) 'Creating background model for the inversion...'
!        allocate(bg_cond(TEgrid%Ny,TEgrid%Nz-TEgrid%Nza))
!        bg_cond = bg_cond_value
!        call create_modelParam(TEgrid,paramtype,sigma1,bg_cond)
!        call write_Cond2D(fidWrite,wFile_Model,sigma1,TEgrid)
!        deallocate(bg_cond)
     
     case (FORWARD_PRED)
        write(*,*) 'Calculating predicted data...'
        call fwdPred(sigma0,allData)
        ! write out impedances
        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (FORWARD_SOLN)
        write(*,*) 'Calculating predicted data and saving the EM solution...'
        call fwdPred(sigma0,allData,eAll)
        ! write out EM solutions
        call write_EMsolnMTX(fidWrite,wFile_EMsoln,eAll)
        ! write out all impedances
        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (SENSITIVITIES)
        write(*,*) 'Calculating the full sensitivity matrix...'
        call calcSensMatrix(allData,sigma0,sigma)
        header = 'Sensitivity Matrix' 
        call writeAll_Cond2D(fidWrite,wFile_Sens,header,   &
                        allData%nData,sigma)

     case (MULT_BY_J)
        write(*,*) 'Multiplying by J...'
        call read_Cond2D(fidRead,rFile_dModel,dsigma,TEgrid)
        call Jmult(dsigma,sigma0,allData) 
        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)

     case (MULT_BY_J_T)
        write(*,*) 'Multiplying by J^T...'
        call JmultT(sigma0,allData,dsigma) 
        call write_Cond2D(fidWrite,wFile_dModel,dsigma,TEgrid)

     case (MULT_BY_J_T_MTX)
        write(*,*) 'Multiplying by J^T (all transmitters)...'
        call JmultT_MTX(sigma0,allData,sigma) 
        header = 'Sensitivity Matrix' 
        call writeAll_Cond2D(fidWrite,wFile_dModelMTX,header,   &
                        allData%nTx,sigma)

     case (MULT_BY_J_MTX)
        write(*,*) 'Multiplying by J (all transmitters)...'
        nTx = allData%nTx
        paramtype = ''
        allocate(sigma(nTx))
        do i = 1,nTx 
           call create_ModelParam(TEgrid,paramtype,sigma(i))
        enddo
        call readAll_Cond2D(fidRead,rFile_dModelMTX,nTx,header,sigma)
        call Jmult_MTX(sigma,sigma0,allData) 
        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,   &
			nSites,sites,allData)
			
      case (FORWARD_PRED_DELTA)
        write(*,*) 'Calculating the predicted data for the perturbed model...'
        call read_Cond2D(fidRead,rFile_dModel,dsigma,TEgrid)
        !  add sigma0 to epsilon * dsigma
        eps = 0.01
        call linComb_modelParam(ONE,sigma0,eps,dsigma,dsigma)
        call fwdPred(dsigma,allData)
        ! write out impedances computed with perturbed sigma
        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,nSites,sites,allData)

!     case (CREATE_DATA_ERRORS)
!        write(*,*) 'Writing a data file with 5% errors in it...'
!        error = 0.05
!        call setError_dvecMTX(error,allData)
!        ! write out impedances
!        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,   &
!			nSites,sites,allData)
			     		
     case (NLCG_INVERSION)
        write(*,*) 'Starting the NLCG search...'
        call NLCGsolver(allData,lambda,sigma0,sigma1,alpha)
        call write_Cond2D(fidWrite,wFile_Model,sigma1,TEgrid)
        call fwdPred(sigma1,allData)
        call write_Z(fidWrite,wFile_Data,nPer,periods,modes,nSites,sites,allData)
     end select

	 call date_and_time(values=tarray)
	 etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	 rtime = etime - stime
	 write(0,*) ' elapsed time (mins) ',rtime/60.0

end program