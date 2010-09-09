! *****************************************************************************
module ForwardSolver
  !  High level interface/control module used by top level routines
  !   for initializing and using the solver.  The key public routines
  !   in this module have only generic (abstract data type) arguments
  !   and can thus be called from the top-level inversion routines.
  !
  !  All public routines must have the same generic names, parameter lists
  !   and abstract functionality; this implementation is for EARTH
  !
  !  Idea is to always call initSolver before calling fwdSolver,
  !   in particular before the first solution for each transmitter
  !   (frequency).  If called for the first time (in a program run,
  !   or after a call to exitSolver), full initialization
  !   (after deallocation/cleanup if required) is performed.

  use math_constants
  use SolnSpace
  use jacobian
  use transmitters
  use input
  use output
  use UserData
  use initFields
  use modelmap

  implicit none

  public :: initSolver, fwdSolve, sensSolve, exitSolver

  save
  private
  type(timer_t)     :: timer
  logical           :: solverInitialized = .false.
  logical           :: secondaryField = .false.
  logical           :: modelUpdated = .false.
  type(rscalar)     :: resist ! resistivity on the grid
  type(modelParam_t)    :: mPrev  ! store the previous model for efficiency
  type(solnVector_t)    :: hPrev  ! store the previous forward solver solution
  type(rhsVector_t)     :: b   ! needed to store the RHS for forward solver
  type(solnVectorMTX_t) :: h1d ! needed for secondary field formulation

Contains

  !**********************************************************************
  subroutine initSolver(iTx,m0,grid,h0,h,comb)
   !  Initializes both forward and sens solvers for transmitter iTx.

   integer, intent(in)                      :: iTx
   type(modelParam_t),intent(in), target	:: m0
   type(grid_t), intent(in), target         :: grid
   type(solnVector_t), intent(inout), optional :: h0
   type(solnVector_t), intent(inout), optional :: h
   type(rhsVector_t), intent(inout), optional :: comb
   !  local variables
   type(transmitter_t), pointer             :: freq
   logical                                  :: initFwd, initForSens

   initFwd = present(h0)
   initForSens = present(comb)

   freq => freqList%info(iTx)
   write(*,'(a46,es9.3,a5)') &
        'Initializing 3D SGFD global solver for period ',freq%period,' days'

   if(initFwd) then
     if(h0%allocated) then
        ! use the previous forward solution as starting solution
        ! for this frequency ... b should already exist
     else
        ! no forward solution computed yet... initialize
        call create_solnVector(grid,iTx,h0)
        b%nonzero_source = .true.
        call create_rhsVector(grid,iTx,b)
        call initialize_fields(h0%vec,b%source,grid)
     endif
   endif

   if(initForSens) then
        ! also initialize the optional outputs
      call create_solnVector(grid,iTx,h)
      comb%nonzero_source = .true.
      call create_rhsVector(grid,iTx,comb)
   endif

   secondaryField = freq%secondaryField

   if(.not. solverInitialized) then
      ! only do this when called for the first time
      mPrev = m0
      if(secondaryField) then
        call initField(cUserDef,grid,h1d)
      endif
      solverInitialized = .true.
   endif

   ! this will only be true when new model is supplied (m0 /= mPrev)
   modelUpdated = .true.

   if(modelUpdated) then
      ! compute the resistivity on the grid
      call create_rscalar(grid,resist,CENTER)
      call initModel(grid,m0,resist%v)
   endif

   ! reset timer
   call reset_time(timer)

  end subroutine initSolver

  !**********************************************************************
  subroutine fwdSolve(iTx,h)
!  Forward solver; uses b as the RHS. Similar to sensSolve, less general.
!  Temporary routine (will be merged with sensSolve).

   integer, intent(in)                        :: iTx
   type(solnVector_t), intent(inout)            :: h
   ! local variables
   type(transmitter_t), pointer                 :: freq
   real(kind=prec)                              :: omega
   logical                                      :: adjoint,sens

   call sensSolve(iTx,FWD,h)

  end subroutine fwdSolve

  !**********************************************************************
  subroutine sensSolve(iTx,FWDorADJ,h,comb)
!  Generic forward solver; if comb is specified, it is used as forcing;
!  otherwise, assume that the RHS is already stored in b.
!  Use the input e for starting solution
!  (normally, e will be zero on input for sensitivities).

   integer, intent(in)                        :: iTx
   character (len=3), intent(in)                :: FWDorADJ
   type(solnVector_t), intent(inout)            :: h
   type(rhsVector_t), intent(inout), optional    :: comb
   ! local variables
   type(transmitter_t), pointer                 :: freq
   real(kind=prec)                              :: omega
   type(cvector)                                :: source
   logical                                      :: adjoint,sens

   if(.not. present(comb)) then
    ! use b as RHS
    adjoint = .false.
    sens = .false.
    source = b%source
   else
    adjoint = (FWDorADJ .ne. FWD)
    sens = .true.
    source = comb%source
   endif

   freq => freqList%info(iTx)

   ! run FWD/ADJ solver
   write(*,'(a12,a3,a20,i4,a2,es12.6,a5)') &
    'Solving the ',FWDorADJ,' problem for period ',iTx,': ',1/freq%value,' secs'

   omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

   ! solve S_m <h> = <b> for vector <h>
   if (.not. adjoint) then
    !call operatorD_l_mult(source,h%grid)
   end if
   call operatorM(h%vec,source,omega,resist%v,h%grid,fwdCtrls,h%errflag,adjoint,sens)
   if (adjoint) then
    !call operatorD_l_divide(h%vec,h%grid)
   end if

   ! compute and output fields & C and D responses at cells
   call outputSolution(freq,h%vec,slices,h%grid,cUserDef,resist%v,'h')

   ! output full H-field cvector
   if (output_level > 3) then
      call outputField(freq,h%vec,cUserDef,'field')
   end if

   ! update pointer to the transmitter in solnVector
   h%tx = iTx

   ! clean up
   call deall(source)

   if (output_level > 1) then
      write (*,*) ' time taken (mins) ', elapsed_time(timer)/60.0
   end if


  end subroutine sensSolve


  !**********************************************************************
  subroutine exitSolver(e0,e,comb)
   ! Cleans up after fwdSolver and deallocates all solver data. Call
   ! initSolver to update solver data; call this to exit completely.
   ! Optionally, deallocates e0,e,comb
   type(solnVector_t), intent(inout), optional  :: e0
   type(solnVector_t), intent(inout), optional  ::e
   type(rhsVector_t), intent(inout), optional   ::comb

   ! local variables
   logical          :: initForSens

   initForSens = present(comb)

   if(present(e0)) then
      call deall_solnVector(e0)
   endif

   if(initForSens) then
      call deall_rhsVector(comb)
      call deall_solnVector(e)
   endif

   if(solverInitialized) then
      ! cleanup/deallocation routines for model operators
      if(secondaryField) then
        call deall_solnVectorMTX(h1d)
        secondaryField = .false.
      endif
      call deall_rhsVector(b)
      call deall_solnVector(hPrev)
      call deall_modelParam(mPrev)
      call deall_rscalar(resist)
      solverInitialized = .false.
      modelUpdated = .false.
   endif

   ! restart the clock
   call clear_time(timer)

  end subroutine exitSolver


end module ForwardSolver
