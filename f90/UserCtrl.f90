! *****************************************************************************
module UserCtrl
  ! This module defines the derived data type structure with all filenames

  implicit none

  character*1, parameter    :: READ_WRITE = 'R'
  character*1, parameter	:: FORWARD = 'F'
  character*1, parameter	:: COMPUTE_J = 'J'
  character*1, parameter	:: MULT_BY_J = 'M'
  character*1, parameter	:: MULT_BY_J_T = 'T'
  character*1, parameter	:: INVERSE = 'I'
  character*1, parameter	:: TEST_COV = 'C'

  ! ***************************************************************************
  ! * input_info contains the list of all essential input information currently
  ! * read in from fn_startup.
  type :: userdef_control

	! Options: FORWARD, COMPUTE_J, MULT_BY_J, MULT_BY_J_T,
	!          NLCG, TEST_COV, READ_WRITE
	character(80)       :: job

	! File to set up inversion controls
	character(80)       :: rFile_invCtrl

	! File to set up forward solver controls
	character(80)       :: rFile_fwdCtrl

	! Input files
	character(80)       :: rFile_Grid, rFile_Model, rFile_Data
	character(80)       :: rFile_dModel

	! Output files
	character(80)       :: wFile_Grid, wFile_Model, wFile_Data
	character(80)       :: wFile_dModel
	character(80)       :: wFile_EMsoln, wFile_Sens

	! Specify damping parameter for the inversion
	real(8)             :: lambda

	! Initial step size parameter for inversions
	real(8)             :: delta

	! Specify covariance configuration
	character(80)       :: rFile_Cov

	! Choose the inverse search algorithm
	character(80)       :: search

	! Indicate how much output you want
	character(80)       :: verbose

  end type userdef_control

Contains

  subroutine initUserCtrl(ctrl)

  	type(userdef_control), intent(out)   :: ctrl

  	ctrl%job = ''
  	ctrl%rFile_invCtrl = ''
  	ctrl%rFile_fwdCtrl = ''
  	ctrl%rFile_Grid = ''
  	ctrl%wFile_Grid = ''
  	ctrl%rFile_Model = ''
  	ctrl%wFile_Model = ''
  	ctrl%rFile_Data = ''
  	ctrl%wFile_Data = ''
  	ctrl%rFile_dModel = ''
  	ctrl%wFile_dModel = ''
  	ctrl%wFile_EMsoln = ''
  	ctrl%wFile_Sens = ''
  	ctrl%lambda = 1
  	ctrl%delta = 1
  	ctrl%rFile_Cov = ''
  	ctrl%search = 'NLCG'
  	ctrl%verbose = 'regular'

  end subroutine initUserCtrl


  subroutine parseArgs(program,ctrl)

     character(1)     :: job
     integer (kind=4) :: iargc,narg,k
     character*80  gridType, header,arg, paramtype
     character*80, dimension(:), pointer :: temp
     character(*), intent(in)            :: program
     type(userdef_control), intent(out)  :: ctrl

     call initUserCtrl(ctrl)

     !  parse command line ...  for now just uses first argument
     !   to set job
     narg = iargc()
    if(narg .gt. 0) then
       k=1
       search_arg: &
		   do
		       k=k+1
		       if (k .eq. narg) exit  search_arg         
               call getarg ( k, arg )
                if(arg(1:1).eq.'-') then
                  narg=k-1
                  exit  search_arg
                end if
          end do search_arg     
     end if 

     if(narg .gt. 0) then
        call getarg(1,arg)
        if(arg(1:1).eq.'-') job = arg(2:2)
     else
        write(*,*) 'Usage: ',trim(program),' -[job] [args]'
        write(*,*)
        write(*,*) '[READ_WRITE]'
        write(*,*) ' -R  rFile_Model rFile_Data [wFile_Model wFile_Data]'
        write(*,*) '  Reads your input files and checks them for validity;'
        write(*,*) '  optionally also writes them out'
        write(*,*) '[FORWARD]'
        write(*,*) ' -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln]'
        write(*,*) '  Calculates the predicted data and saves the EM solution'
        write(*,*) '[COMPUTE_J]'
        write(*,*) ' -J  rFile_Model rFile_Data wFile_Sens'
        write(*,*) '  Calculates and saves the full J(acobian)'
        write(*,*) '[MULT_BY_J]'
        write(*,*) ' -M  rFile_Model rFile_dModel rFile_Data wFile_Data'
        write(*,*) '  Multiplies a model by J to create a data vector'
        write(*,*) '[MULT_BY_J_T]'
        write(*,*) ' -T  rFile_Model rFile_Data wFile_dModel'
        write(*,*) '  Multiplies a data vector by J^T to create a model'
        write(*,*) '[INVERSE]'
        write(*,*) ' -I NLCG rFile_Model rFile_Data wFile_Model [wFile_Data rFile_Cov lambda delta]'
        write(*,*) '  Runs an inverse search to yield an inverse model'
        write(*,*) '  Here, lambda =  the initial damping parameter'
        write(*,*) '        delta  =  the initial line search step size in the model units'
        write(*,*) '[TEST_COV]'
        write(*,*) ' -C  rFile_Model wFile_Model [rFile_Cov]'
        write(*,*)
        stop
     endif

     ! extract all following command line arguments
     allocate(temp(1:narg-1))
     do k = 1,narg-1
       call getarg(k+1,temp(k))
     end do

     ! set narg to the number of actual arguments (after the job is specified)
     narg = narg-1

     select case (job)

      case (READ_WRITE) !R
        if (narg < 2) then
           write(0,*) 'Usage: -R  rFile_Model rFile_Data [wFile_Model wFile_Data]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	    end if
	    if (narg > 2) then
	       ctrl%wFile_Model = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%wFile_Data = temp(4)
	    end if

      case (FORWARD) !F
        if (narg < 3) then
           write(0,*) 'Usage: -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_Data = temp(3)
	    end if
	    if (narg > 3) then
	       ctrl%wFile_EMsoln = temp(4)
	    end if

      case (COMPUTE_J) ! J
        if (narg < 3) then
           write(0,*) 'Usage: -J  rFile_Model rFile_Data wFile_Sens'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_Sens = temp(3)
	    end if

      case (MULT_BY_J) ! M
        if (narg < 4) then
           write(0,*) 'Usage: -M  rFile_Model rFile_dModel rFile_Data wFile_Data'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_dModel = temp(2)
	       ctrl%rFile_Data = temp(3)
	       ctrl%wFile_Data = temp(4)
	    end if

      case (MULT_BY_J_T) ! T
        if (narg < 3) then
           write(0,*) 'Usage: -T  rFile_Model rFile_Data wFile_dModel'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_dModel = temp(3)
	    end if

      case (INVERSE) ! I
        if (narg < 4) then
           write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data wFile_Model [wFile_Data rFile_Cov lambda delta]'
           write(0,*) 'Here, lambda =  the initial damping parameter'
           write(0,*) '      delta  =  the initial line search step size in the model units'
           stop
        else
           ctrl%search = temp(1)
           select case (ctrl%search)
           case ('NLCG','DCG','Hybrid')
              	! write(0,*) 'Inverse search ',trim(ctrl%search),' selected.'
           case default
				write(0,*) 'Unknown inverse search. Usage: -I [NLCG | DCG | Hybrid]'
				stop
           end select
	       ctrl%rFile_Model = temp(2)
	       ctrl%rFile_Data = temp(3)
	       ctrl%wFile_Model = temp(4)
	    end if
	    if (narg > 4) then
	       ctrl%wFile_Data = temp(5)
	    end if
	    if (narg > 5) then
	       ctrl%rFile_Cov = temp(6)
	    end if
	    if (narg > 6) then
          read(temp(7),*) ctrl%lambda
        end if
        if (narg > 7) then
          read(temp(8),*) ctrl%delta
        end if

      case (TEST_COV) ! C
        if (narg < 2) then
           write(0,*) 'Usage: -C  rFile_Model wFile_Model [rFile_Cov]'
           stop
        else
	    ctrl%rFile_Model = temp(1)
	    ctrl%wFile_Model = temp(2)
        end if
        if (narg > 2) then
            ctrl%rFile_Cov = temp(3)
        end if

      case default
         write(0,*) 'Unknown job. Please check your command line options'
         stop

     end select

     deallocate(temp)

     ! save this info for the main program
     ctrl%job = job

  end subroutine parseArgs

end module UserCtrl
