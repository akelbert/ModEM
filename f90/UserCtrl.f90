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
	!          INVERSE, TEST_COV, READ_WRITE
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

	! Misfit tolerance for the forward solver
	real(8)             :: eps

	! Specify covariance configuration
	character(80)       :: rFile_Cov

	! Choose the inverse search algorithm
	character(80)       :: search

	! Indicate how much output you want
	integer             :: output_level

  end type userdef_control

Contains

  ! Initialize userCtrl. The defaults are not empty strings.
  ! This is because the inquire statements with empty file names
  ! are not understood on IBM machines.

  subroutine initUserCtrl(ctrl)

  	type(userdef_control), intent(out)   :: ctrl

  	ctrl%job = ''
  	ctrl%rFile_invCtrl = 'n'
  	ctrl%rFile_fwdCtrl = 'n'
  	ctrl%rFile_Grid = 'n'
  	ctrl%wFile_Grid = 'n'
  	ctrl%rFile_Model = 'n'
  	ctrl%wFile_Model = 'n'
  	ctrl%rFile_Data = 'n'
  	ctrl%wFile_Data = 'n'
  	ctrl%rFile_dModel = 'n'
  	ctrl%wFile_dModel = 'n'
  	ctrl%wFile_EMsoln = 'n'
  	ctrl%wFile_Sens = 'n'
  	ctrl%lambda = 1
  	ctrl%eps = 1.0e-7
  	ctrl%rFile_Cov = 'n'
  	ctrl%search = 'NLCG'
  	ctrl%output_level = 3

  end subroutine initUserCtrl


  subroutine parseArgs(program,ctrl)

     character(1)     :: job
     integer (kind=4) :: iargc,narg,k
     integer          :: istat
     logical          :: exists
     character*80  gridType, header,arg, verbose, paramtype
     character*80, dimension(:), pointer :: temp
     character(*), intent(in)            :: program
     type(userdef_control), intent(out)  :: ctrl

     call initUserCtrl(ctrl)

     !  parse command line ...  for now just uses first argument
     !   to set job
     narg = iargc()

     !  quick fix against compilers which add additional parameters to mpirun;
     !   only read the first argument with a dash, or the argument -v followed
     !   by the verbose level; ignore the others
     verbose = 'regular'
     if(narg .gt. 1) then
       k=1
       search_arg: &
		   do
		       k=k+1
		       if (k .eq. narg) exit  search_arg
               call getarg ( k, arg )
                if(arg(1:1).eq.'-') then
                  if (arg(2:2).eq.'v') then
                    call getarg ( k+1, verbose )
                  end if
                  narg=k-1
                  exit  search_arg
                end if
          end do search_arg
     end if

	 !  set the level of output based on the user input
	 select case (verbose)
	 case ('debug')
	   print *,'Output all information including debugging lines.'
	   ctrl%output_level = 5
	 case ('full')
	   print *,'Output full information to screen and to files.'
	   ctrl%output_level = 4
	 case ('regular')
	   print *,'Output information to files, and progress report to screen (default).'
	   ctrl%output_level = 3
	 case ('compact')
	   print *,'Output information to files, and compact summary to screen.'
	   ctrl%output_level = 2
	 case ('result')
	   print *,'Output minimal information to files, and result to screen.'
	   ctrl%output_level = 1
	 case ('none')
	   print *,'Output nothing at all except result to screen and to files.'
	   ctrl%output_level = 0
	 case default
	   ctrl%output_level = 3
	 end select

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
        write(*,*) ' -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl]'
        write(*,*) '  Calculates the predicted data and saves the EM solution'
        write(*,*) '[COMPUTE_J]'
        write(*,*) ' -J  rFile_Model rFile_Data wFile_Sens [rFile_fwdCtrl]'
        write(*,*) '  Calculates and saves the full J(acobian)'
        write(*,*) '[MULT_BY_J]'
        write(*,*) ' -M  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
        write(*,*) '  Multiplies a model by J to create a data vector'
        write(*,*) '[MULT_BY_J_T]'
        write(*,*) ' -T  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
        write(*,*) '  Multiplies a data vector by J^T to create a model'
        write(*,*) '[INVERSE]'
        write(*,*) ' -I NLCG rFile_Model rFile_Data [lambda eps]'
        write(*,*) '  Here, lambda = the initial damping parameter for inversion'
        write(*,*) '           eps = misfit tolerance for the forward solver'
        write(*,*) 'OR'
        write(*,*) ' -I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]'
        write(0,*) '  Optionally, may also supply'
        write(0,*) '      the model covariance configuration file   [rFile_Cov]'
        write(0,*) '      the starting model parameter perturbation [rFile_dModel]'
        write(*,*) '  Runs an inverse search to yield an inverse model at every iteration'
        write(*,*) '[TEST_COV]'
        write(*,*) ' -C  rFile_Model wFile_Model [rFile_Cov]'
        write(*,*) '  Applies the model covariance to produce a smooth model output'
        write(*,*)
        write(*,*) 'Optional final argument -v [debug|full|regular|compact|result|none]'
        write(*,*) 'indicates the desired level of output to screen and to files.'
        write(*,*)
        stop
     endif

     ! extract all following command line arguments
     allocate(temp(1:narg-1))
     do k = 1,narg-1
       call getarg(k+1,temp(k))
     end do

     ! write feedback to standard output for records
     if (narg > 1) then
        write(6,*) 'Running job -',job,' with command line options:'
        do k = 1,narg-1
          write(6,'(x1a)',advance='no') trim(temp(k))
        end do
        write(6,*)
        write(6,*)
     end if

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
           write(0,*) 'Usage: -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl]'
           write(0,*)
           write(0,*) 'Here, rFile_fwdCtrl is the forward solver control file in the format'
           write(0,*)
           write(0,*) 'Number of QMR iters per divergence correction : 40'
           write(0,*) 'Maximum number of divergence correction calls : 20'
           write(0,*) 'Maximum number of divergence correction iters : 100'
           write(0,*) 'Misfit tolerance for EM forward solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for EM adjoint solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for divergence correction    : 1.0e-5'
           write(0,*)
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
           write(0,*) 'Usage: -J  rFile_Model rFile_Data wFile_Sens [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_Sens = temp(3)
	    end if

      case (MULT_BY_J) ! M
        if (narg < 4) then
           write(0,*) 'Usage: -M  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_dModel = temp(2)
	       ctrl%rFile_Data = temp(3)
	       ctrl%wFile_Data = temp(4)
	    end if

      case (MULT_BY_J_T) ! T
        if (narg < 3) then
           write(0,*) 'Usage: -T  rFile_Model rFile_Data wFile_dModel [rFile_fwdCtrl]'
           stop
        else
	       ctrl%rFile_Model = temp(1)
	       ctrl%rFile_Data = temp(2)
	       ctrl%wFile_dModel = temp(3)
	    end if

      case (INVERSE) ! I
        if (narg < 3) then
           write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [lambda eps]'
           write(0,*)
           write(0,*) 'Here, lambda = the initial damping parameter for inversion'
           write(0,*) '         eps = misfit tolerance for the forward solver'
           write(0,*) 'OR'
           write(0,*) 'Usage: -I NLCG rFile_Model rFile_Data [rFile_invCtrl rFile_fwdCtrl]'
           write(0,*)
           write(0,*) 'Here, rFile_invCtrl = the inversion control file in the format'
           write(0,*)
           write(0,*) 'Model and data output file name    : Example'
           write(0,*) 'Initial damping factor lambda      : 1.0'
           write(0,*) 'Initial search step in model units : 1.0'
           write(0,*) 'Restart when rms diff is less than : 2.0e-3'
           write(0,*) 'Exit search when rms is less than  : 1.05'
           write(0,*) 'Exit when lambda is less than      : 1.0e-4'
           write(0,*) 'Maximum number of iterations       : 120'
           write(0,*)
           write(0,*) '      rFile_fwdCtrl = the forward solver control file in the format'
           write(0,*)
           write(0,*) 'Number of QMR iters per divergence correction : 40'
           write(0,*) 'Maximum number of divergence correction calls : 20'
           write(0,*) 'Maximum number of divergence correction iters : 100'
           write(0,*) 'Misfit tolerance for EM forward solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for EM adjoint solver        : 1.0e-7'
           write(0,*) 'Misfit tolerance for divergence correction    : 1.0e-5'
           write(0,*)
           write(0,*) 'Optionally, may also supply'
           write(0,*)
           write(0,*) '      rFile_Cov     = the model covariance configuration file'
           write(0,*) '      rFile_dModel  = the starting model parameter perturbation'
           write(0,*)
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
	    end if
	    if (narg > 3) then
          read(temp(4),*,iostat=istat) ctrl%lambda
          if (istat .ne. 0) then
            ! check for the inverse solver configuration file
            ctrl%rFile_invCtrl = temp(4)
            inquire(FILE=ctrl%rFile_invCtrl,EXIST=exists)
            if (.not. exists) then
            	! problem - invalid argument
				write(0,*) 'Please specify a valid inverse control file or damping parameter'
				stop
            end if
          end if
        end if
        if (narg > 4) then
          read(temp(5),*,iostat=istat) ctrl%eps
          if (istat .ne. 0) then
            ! check for the forward solver configuration file
            ctrl%rFile_fwdCtrl = temp(5)
            inquire(FILE=ctrl%rFile_fwdCtrl,EXIST=exists)
            if (.not. exists) then
            	! problem - invalid argument
				write(0,*) 'Please specify a valid forward solver control file or misfit tolerance'
				stop
            end if
          end if
        end if
	    if (narg > 5) then
	       ctrl%rFile_Cov = temp(6)
	    end if
        if (narg > 6) then
            ! check for the optional starting model file
            ctrl%rFile_dModel = temp(7)
            inquire(FILE=ctrl%rFile_dModel,EXIST=exists)
            if (.not. exists) then
            	! problem - invalid argument
				write(0,*) 'Please specify a valid starting model file'
				stop
            end if
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
