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
use datafunc ! inherit solnrhs soln2d modelparameter datagridinfo
use dataspace
use fwdtemod
use fwdtmmod
use solnrhs

implicit none
  
!    keep data structures used only by
!    routines in this module private
!   rhs data structures for solving forward, sensitivity probs
type(EMrhs), save, private		:: b0
!    keep track of which mode was solved for most recently
!    (to minimize reinitialization ... not clear this is needed!)
character*2, save, public		:: currentMode = '  '

! new type to store multiple EM solutions, one for each transmitter
!  NOTE: creation routine for this data type cannot be at the level
!   of the solnrhs module, since this depends on typeDict (not available
!    for the lower level routine).  Might as well also leave the
!    type definition here, since you can't use the type until instances
!    can be created!
type :: EMsolnMTX
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

Contains

   !**********************************************************************
   subroutine initSolver(iDT,iTx,sigma,e0,e,comb)
   !   Initializes forward solver for data type iDt,
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
   !   iDt defines "data type" ... for 2D MT this provide info
   !       about TE/TM mode;
   !   iTx defines transmitter: for MT, just frequency 
   !   
   !   This now does all setup (including matrix factorization) for
   !     the appropriate mode/frequency
   !   NOTE: e and comb are optional calling arguments;
   !     both should be present if one is

   integer, intent(in)				:: iDT
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
   character*80         			:: gridType
   character*2          			:: mode
   logical					:: initForSens

   initForSens = present(comb)

   mode = typeDict(iDT)%mode
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
            gridType = NODE
            if(IER.lt.0) then
              call errStop('initializing for TE mode in initSolver')
            endif
         case('TM')
            call FWD2DSetupTM(SolnRHS_grid,sigma,IER)
            gridType = NODE_EARTH
            if(IER.lt.0) then
              call errStop('initializing for TM mode in initSolver')
            endif
         case default
            call errStop('mode must be TE or TM in initSolver')
      end select
   
      !  allocate for rhs for background, scratch sensitivity solutions
      b0%mode = mode
      b0%nonzero_source = .false.
      b0%nonzero_bc = .true.
      b0%adj = 'FWD'
      call create_EMrhs(SolnRHS_grid,gridType,b0)

      !  allocate for background solution
      call create_EMsoln(SolnRHS_grid,gridType,e0)

      if(initForSens) then
         !  allocate for sensitivity solution, RHS
         call create_EMsoln(solnRHS_grid,gridType,e)
         comb%mode = mode
         comb%nonzero_source = .true.
         comb%nonzero_bc = .false.
         comb%adj = ''
         call create_EMrhs(SolnRHS_grid,gridType,comb)
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
   subroutine fwdSolve(iTx,iDT,e0)

   !  driver for 2d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b) is private to this module
   !  NOTE that this routine now does no initialization: this must be
   !   done with a call to initSolver before calling for the first
   !   time FOR EACH DATA TYPE/TRANSMITTER

   integer, intent(in)		:: iTx,iDT
   type(EMsoln), intent(inout)	:: e0

   ! local variables
   real(kind=selectedPrec)	:: period,omega
   integer			:: IER
   complex(kind=selectedPrec)	:: i_omega_mu
   character*2          	:: mode

   mode = typeDict(iDT)%mode
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
   subroutine sensSolve(iTx,iDT,FWDorADJ,comb,e)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine now does no initialization: this must be
   !   done with a call to initSolver before calling for the first
   !   time FOR EACH DATA TYPE/TRANSMITTER

   integer, intent(in)          	:: iTx,iDT
   character*3, intent(in)		    :: FWDorADJ
   type(EMrhs), intent(inout)		:: comb
   type(EMsoln), intent(inout)		:: e

   ! local variables
   integer      :: IER

   comb%adj = FWDorADJ
   if(typeDict(iDT)%mode.eq.'TE') then
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
      character*80                      :: gridType

      call deall_EMsolnMTX(eAll)

      eAll%nTx = d%nTx
      allocate(eAll%solns(d%nTx))
      eAll%allocated = .true.
      do j = 1,d%nTx
         if(typeDict(d%d(j)%dataType)%mode .eq. 'TE') then
            gridType = NODE
         else
            gridType = NODE_EARTH
         endif
         call create_EMsoln(SolnRHS_grid,gridType,eAll%solns(j))
      enddo

   end subroutine create_EMsolnMTX
   
   !**********************************************************************
   subroutine deall_EMsolnMTX(eAll)

      type(EMsolnMTX), intent(inout)     :: eAll

      !  local variables
      integer                           :: j
      character*80                      :: gridType

	  do j = 1,eAll%nTx
	  	call deall_EMsoln(eAll%solns(j))
	  end do

      if (associated(eAll%solns)) deallocate(eAll%solns)
      eAll%allocated = .false.

   end subroutine deall_EMsolnMTX

end module emsolver
