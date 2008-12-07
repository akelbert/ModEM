module emsolver

!  High level interface/control module used by top level routines
!   for initializing and using the solver.  The key public routines
!   in this module have only generic (abstract data type) arguments
!   and can thus be called from the top-level inversion routines.
!  A similar solver interface will be required for any specific use
!   of the top-level inversion modules
!  All public routines must have the same generic names, parameter lists
!   and abstract functionality as in this example for 2D MT.

use math_constants
use utilities
use datafunc ! inherit solnrhs soln2d modelparameter
use dataspace
use fwdtemod
use fwdtmmod
use solnrhs

implicit none
  
 type :: MTtx
  SEQUENCE
     !  An MT source is defined by frequency and boundary conditions
     !   at present there does not seem to be much need for BC info ... add
     !    if needed.  Other sorts of EM data may have more
     !    complex tx descriptions
     character(2)                :: mode = ''! = 'TE' or 'TM'
     ! angular frequency (radians/sec), and for convenience period (s)
     real(kind=selectedPrec)            :: omega = R_ZERO
     real(kind=selectedPrec)            :: period = R_ZERO
     ! index number to frequency/ period in solution file
     integer                    :: iPer
  end type MTtx

  ! transmitter dictionary for MT data will be an array of
  ! type mt_forcing (one element  for each frequency)
  ! Perhaps this should be moved to EMsolver module (and be private
  !    to this module)
  type (MTtx), pointer, save, public, dimension (:)   :: txDict


!    keep data structures used only by
!    routines in this module private
!   rhs data structures for solving forward, sensitivity probs
type(EMrhs), save, private		:: b0
!    keep track of which mode was solved for most recently
!    (to minimize reinitialization ... not clear this is needed!)
character*2, save, public		:: currentMode = '  '

! new type to store multiple EM solutions, one for each transmitter
!  NOTE: creation routine for this data type cannot be at the level
!   of the solnrhs module, since this depends on txDict (not available
!    for the lower level routine).  Might as well also leave the
!    type definition here, since you can't use the type until instances
!    can be created!
type :: EMsolnMTX
 SEQUENCE
  !  derived data type for storing solutions from multiple transmitters
  integer			:: nTx = 0
  type(EMsoln), pointer		:: solns(:)
  logical			:: allocated = .false.
end type EMsolnMTX

!  initialization routines (call Fwd version if no sensitivities are
!     are calculated).  Note that these routines are set up to
!    automatically manage memory and to figure out which initialization
!    (or reinitialization) steps are required (e.g., when the frequency
!    changes from the previous solver call, appropriate solver
!    coefficients are updated, matrices factored, etc.).  This
!    functionality needs to be maintained in implementations for new
!    problems! 
public initSolver

!  cleanup/deallocation routines
public exitSolver

! solver routines
public fwdSolve, sensSolve

! create/deallocate EMsoln
public create_EMsolnMTX, deall_EMsolnMTX

  !  SolnRHS_grid is used to define grid parameters for DataFunc and
  !    EMsolver modules.  Make a copy of the numerical
  !   grid geometry parameters in this module at the start of
  !   the inversion.

type(grid2d_t), target, save, private     :: SolnRHS_grid
  
Contains

!**********************************************************************

! Initializes and sets up transmitter dictionary for MT,
!  This is just a simple example of a routine for setting up the TX
!   dictionary; In this example we assume that there are nPer periods
!   for either TE or TM modes, or for both.
!  NOTE:   If TXdict if public, there is no reason for this to
!    be part of this module (but I leave it here for now!)

  subroutine TXdictSetUp(nTx,Periods,modes)

     integer, intent(in)         :: nTx
     real*8, intent(in)          :: periods(nTx)
     character*2, intent(in)     :: modes(nTx)
 
     ! local variables
     integer                     :: iTx

     allocate(txDict(nTx))
     do iTx = 1, nTx
        txDict(iTx)%period = Periods(iTx)
        txDict(iTx)%omega = (2*PI)/ txDict(iTx)%period
        txDict(iTx)%mode = modes(iTx)
     enddo

  end subroutine TXdictSetUp
  
! **************************************************************************
! Cleans up and deletes transmitter dictionary at end of program execution
  subroutine deall_txDict()

	integer     :: istat

    if (associated(txDict)) then
       deallocate(txDict,STAT=istat)
    end if

  end subroutine deall_txDict   
  
   !**********************************************************************
   subroutine initSolver(iTx,sigma,e0,e,comb)
   !   Initializes forward solver for
   !    transmitter iTx: in this instance TE or TM mode solvers 
   !    for the appropriate frequency depending on mode
   !   Idea is to call this before calling fwdSolve or sensSolve,
   !     in particular before the first solution for each transmitter
   !     (frequency).  If called for the first time (in a program run,
   !     or after a call to exitSolver), or if the previous call to
   !     this routine initialized for a different data type (TE vs. TM
   !     mode) full initialization (after deallocation/cleanup if required)
   !     is performed.
   !     
   !   iTx defines transmitter: for 2D MT, this provides info about 
   !       frequency and TE/TM mode; for 3D MT just frequency
   !   
   !   This now does all setup (including matrix factorization) for
   !     the appropriate mode/frequency
   !   NOTE: e and comb are optional calling arguments;
   !     both should be present if one is

   integer, intent(in)				:: iTx
   type(modelParam_t),intent(in), target		:: sigma
   !  following structures are initialized/created in this routine
   !	solution vector for forward problem
   type(EMsoln), intent(inout)			:: e0
   !	solution vector for sensitivity
   type(EMsoln), intent(inout), optional	:: e
   !	forcing for sensitivity
   type(EMrhs), intent(inout), optional		:: comb

   !  local variables
   integer					:: IER
   real(kind=selectedPrec)			:: period
   character*2          			:: mode
   logical					:: initForSens

   initForSens = present(comb)

   mode = txDict(iTx)%mode
   period = txDict(iTx)%period

   if(currentMode .ne. mode) then
      if(currentMode .ne. '  ') then
         !  not inital solution, but mode has changed from last
         !  sensitivity calculated: will need to reinitialize 
         !  solver + rhs/soln arrays ... first deallocate from
         !  previous mode 
         call deall_EMrhs(b0)
         call deall_EMsoln(e0)
         if(initForSens) then
            call deall_EMrhs(comb)
            call deall_EMsoln(e)
         endif
         select case(currentMode)
            case('TE')
               call Fwd2DdeallTE()
            case('TM')
               call Fwd2DdeallTM()
            case default
         end select
      endif
      currentMode = mode

      ! initialize coefficient matrix (frequency indpendent part)
      select case(mode)
         case('TE')
            call FWD2DSetupTE(SolnRHS_grid,sigma,IER)
            if(IER.lt.0) then
              call errStop('initializing for TE mode in initSolver')
            endif
         case('TM')
            call FWD2DSetupTM(SolnRHS_grid,sigma,IER)
            if(IER.lt.0) then
              call errStop('initializing for TM mode in initSolver')
            endif
         case default
            call errStop('mode must be TE or TM in initSolver')
      end select
   
      !  allocate for rhs for background, scratch sensitivity solutions
      b0%nonzero_source = .false.
      b0%nonzero_bc = .true.
      b0%adj = 'FWD'
      call create_EMrhs(SolnRHS_grid,mode,b0)

      !  allocate for background solution
      call create_EMsoln(SolnRHS_grid,mode,e0)

      if(initForSens) then
         !  allocate for sensitivity solution, RHS
         call create_EMsoln(solnRHS_grid,mode,e)
         comb%nonzero_source = .true.
         comb%nonzero_bc = .false.
         comb%adj = ''
         call create_EMrhs(SolnRHS_grid,mode,comb)
      endif
   endif

   !   complete initialization of coefficient matrix, then factor
   !   (factored matrix is saved in TE/TM mode modeling module
   !   This needs to be called before solving for a different frequency
   select case (mode)
      case ('TE')
         call UpdateFreqTE(period)
      case ('TM')
         call UpdateFreqTM(period)
      case default
   end select

   end subroutine initSolver

   !**********************************************************************
   subroutine exitSolver(e0,e,comb)
   !   deallocates rhs0, rhs and solver arrays
   type(EMsoln), intent(inout)			:: e0
   type(EMsoln), intent(inout), optional	::e
   type(EMrhs), intent(inout), optional		::comb

   ! local variables
   logical			:: initForSens

   initForSens = present(comb)

   call deall_EMrhs(b0)
   call deall_EMsoln(e0)

   if(initForSens) then
      call deall_EMrhs(comb)
      call deall_EMsoln(e)
   endif

   select case(currentMode)
      case('TE')
         call Fwd2DdeallTE()
      case('TM')
         call Fwd2DdeallTM()
      case default
   end select

   currentMode = '  '

   end subroutine exitSolver

   !**********************************************************************
   subroutine fwdSolve(iTx,e0)

   !  driver for 2d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b) is private to this module
   !  NOTE that this routine now does no initialization: this must be
   !   done with a call to initSolver before calling for the first
   !   time FOR EACH DATA TYPE/TRANSMITTER

   integer, intent(in)		:: iTx
   type(EMsoln), intent(inout)	:: e0

   ! local variables
   real(kind=selectedPrec)	:: period,omega
   integer			:: IER
   complex(kind=selectedPrec)	:: i_omega_mu
   character*2          	:: mode

   mode = txDict(iTx)%mode
   period = txDict(iTx)%period
   omega = txDict(iTx)%omega
   !  set period, complete setup of TE mode equation system
   i_omega_mu = cmplx(0.,isign*mu*omega,kind=selectedPrec)

   ! solve forward problem
   select case (mode)
      case ('TE')
         ! compute boundary conditions for 2D TE forward problem
         call SetBoundTE(period,b0%bc)
         !  solve 2D TE equations, return solution in e0
         call Fwd2DsolveTE(b0,e0%vec%v,IER)
         if(IER .LT. 0) then
            call errStop('solving TE mode equation in fwdSolve')
         endif
      case ('TM')
         ! compute boundary conditions for 2D TM forward problem
         call SetBoundTM(period,b0%bc)
         !  solve 2D TM equations, return solution in e0
         call Fwd2DsolveTM(b0,e0%vec%v,IER)
         if(IER .LT. 0) then
            call errStop('solving TM mode equation in fwdSolve')
         endif
      case default
   end select
   !  set transmitter index for this solution 
   !          (sigma is set at initialization)
   e0%tx = iTx
   e0%mode = mode
   e0%period = period
   e0%omega = omega

   end subroutine fwdSolve

   !**********************************************************************
   subroutine sensSolve(iTx,FWDorADJ,comb,e)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine now does no initialization: this must be
   !   done with a call to initSolver before calling for the first
   !   time FOR EACH DATA TYPE/TRANSMITTER

   integer, intent(in)          	:: iTx
   character*3, intent(in)		    :: FWDorADJ
   type(EMrhs), intent(inout)		:: comb
   type(EMsoln), intent(inout)		:: e

   ! local variables
   integer      :: IER

   comb%adj = FWDorADJ
   if(txDict(iTx)%mode.eq.'TE') then
       call Fwd2DsolveTE(comb,e%vec%v,IER)

       if(IER .LT. 0) then
          call errStop('solving TE mode equation in sensSolve')
       endif
    else
       call Fwd2DsolveTM(comb,e%vec%v,IER)
       if(IER .LT. 0) then
          call errStop('solving TM mode equation in sensSolve')
       endif
   endif
   e%tx = iTx
   end subroutine sensSolve   

!**********************************************************************
   subroutine create_EMsolnMTX(d,eAll)

      type(dvecMTX),intent(in)          :: d
      type(EMsolnMTX), intent(inout)     :: eAll

      !  local variables
      integer                           :: j
      character(2)                      :: mode

      call deall_EMsolnMTX(eAll)

      eAll%nTx = d%nTx
      allocate(eAll%solns(d%nTx))
      eAll%allocated = .true.
      do j = 1,d%nTx
         mode = txDict(d%d(j)%tx)%mode
         call create_EMsoln(SolnRHS_grid,mode,eAll%solns(j))
      enddo

   end subroutine create_EMsolnMTX
   
   !**********************************************************************
   subroutine deall_EMsolnMTX(eAll)

      type(EMsolnMTX), intent(inout)     :: eAll

      !  local variables
      integer                           :: j

	  do j = 1,eAll%nTx
	  	call deall_EMsoln(eAll%solns(j))
	  end do

      if (associated(eAll%solns)) deallocate(eAll%solns)
      eAll%allocated = .false.

   end subroutine deall_EMsolnMTX

!**********************************************************************
    subroutine set_SolnRHS_grid(grid)
!    Call this routine to set basic grid geometry parameters
!       before using any other routines in this module (most depend
!       on saved SolnRHS_grid to define grid geometry)

       type (grid2d_t), intent(in)     :: grid

       SolnRHS_grid = grid

    end subroutine set_SolnRHS_grid
    
!**********************************************************************
    subroutine delete_SolnRHS_grid
!    Call this routine when SolnRHS_grid is no longer needed

       call deall_grid2d(SolnRHS_grid)

    end subroutine delete_SolnRHS_grid
    
end module emsolver
