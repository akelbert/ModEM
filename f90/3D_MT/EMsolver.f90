module emsolver

!  High level interface/control module used by top level routines
!   for initializing and using the solver.  The key public routines
!   in this module have only generic (abstract data type) arguments
!   and can thus be called from the top-level inversion routines.
!  A similar solver interface will be required for any specific use
!   of the top-level inversion modules
!  All public routines must have the same generic names, parameter lists
!   and abstract functionality; this implementation is for 3D-MT
!
!  Main difference from 2D: no need to keep track of TE/TM modes

use math_constants
use datafunc
use dataspace
use emsolve3d
use boundary_ws

implicit none
  
!    keep data structures used only by
!    routines in this module private
!   rhs data structures for solving forward, sensitivity probs
type(RHS), save, private		:: b0

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
public fwdSolve, sensSolve, create_EMsolnMTX 

logical, save, private		:: modelDataInitialized = .false.
!  logical, save, private		:: sigmaNotCurrent = .true.

Contains
   
   !**********************************************************************
   subroutine initSolver(iDT,iTx,sigma,e0,e,comb)
   !   Initializes forward solver for data type iDt, transmitter iTx.
   !     Idea is to call this before calling fwdSolve or sensSolve,
   !     in particular before the first solution for each transmitter
   !     (frequency).  If called for the first time (in a program run,
   !     or after a call to exitSolver), full initialization 
   !     (after deallocation/cleanup if required) is performed.
   !     
   !   iDt defines "data type" ... for 2D MT this provide info
   !       about TE/TM mode;
   !   iTx defines transmitter: for MT, just frequency 
   !   
   !   This now does all setup (including matrix factorization) for
   !     the appropriate mode/frequency
   !   NOTE: e and comb are optional calling arguments;
   !     both should be present if one is
   !   
   !   For 3DMT initialization
   !     of the solver does not depend on iDT which is not
   !     used here.  (Used to distinguish TE/TM in 2D)

   integer, intent(in)				:: iDT
   integer, intent(in)				:: iTx
   type(modelParam_t),intent(in), target		:: sigma
   !  following structures are initialized
   !	solution vector for forward problem
   type(EMsoln), intent(inout)			:: e0
   !	solution vector for sensitivity
   type(EMsoln), intent(inout), optional	:: e
   !	forcing for sensitivity
   type(EMrhs), intent(inout), optional		:: comb

   !  local variables
   integer		:: IER,k
   character*80         :: gridType
   logical		:: initForSens,sigmaNotCurrent

   initForSens = present(comb)
   
	 !  allocate for scratch EMrhs structure for background, sensitivity
   b0%nonzero_Source = .false.
   b0%nonzero_bc = .true.
   b0%adj = 'FWD'
   b0%sparse_Source = .false.
   call create_RHS(SolnRHS_grid,b0)

   !  allocate for background solution
   call create_EMsoln(SolnRHS_grid,e0)
   e0%sigma => sigma

   if(initForSens) then
      !  allocate for sensitivity solution, RHS
      call create_EMsoln(solnRHS_grid,e)
      e%sigma => sigma
      do k = 1,comb%nPol

        comb%b(k)%nonzero_source = .true.
        comb%b(k)%nonzero_bc = .false.
      !   assuming here that we don't use sparse storage ... we could!
        comb%b(k)%sparse_Source = .false.
        comb%b(k)%adj = ''
      enddo
      call create_EMrhs(SolnRHS_grid,comb)
   endif

   if(.NOT.modelDataInitialized) then
   !   Initialize modelData, setup model operators
      call ModelDataInit(SolnRHS_grid)
      call ModelOperatorSetup()
      modelDataInitialized = .true.
   endif

!    the following needs work ... want to avoid reinitializing
!     operator coefficients when conductivity does not change;
!     need to have a way to reset sigmaNotCurrent to false when
!     conductivity changes (one idea: assign a random number whenever
!     a conductivity parameter is modified (by any of the routines in
!     module ModelParameter); store this in the modelOperator module (which
!     is where updateCond sits) and have updateCond compare the random key
!     with what is stored)
!  if(sigmaNotCurrent) then
       call updateCond(sigma)
!      sigmaNotCurrent = .false.
!   endif

   ! This needs to be called before solving for a different frequency
   !!!!!!!  BUT AFTER UPDATECOND !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call UpdateFreq(txDict(iTx)%omega)

   end subroutine initSolver

   !**********************************************************************
   subroutine exitSolver(e0,e,comb)
   !   deallocates b0, comb, e0, e and solver arrays
   type(EMsoln), intent(inout)			:: e0
   type(EMsoln), intent(inout), optional	::e
   type(EMrhs), intent(inout), optional		::comb

   ! local variables
   logical			:: initForSens

   initForSens = present(comb)

   call deall_RHS(b0)
   call deall_EMsoln(e0)

   if(initForSens) then
      call deall_EMrhs(comb)
      call deall_EMsoln(e)
   endif

   !  Need cleanup/deallocation routines for model operators
   end subroutine exitSolver

   !**********************************************************************
   subroutine fwdSolve(iTx,iDT,e0)

   !  driver for 3d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b0) is generated locally--i.e.
   !   boundary conditions are set internally (NOTE: could use transmitter
   !   dictionary to indireclty provide information about boundary
   !    conditions.  Presently we set BC using WS approach.   
   !  NOTE that this routine calls UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.
   !  NOTE: iDT is actually not used in the 3DMT version, as solver does not
   !       depend on data type

   integer, intent(in)		:: iTx,iDT
   type(EMsoln), intent(inout)	:: e0

   ! local variables
   real(kind=selectedPrec)	:: period, omega
   integer			:: IER,iMode
   complex(kind=selectedPrec)	:: i_omega_mu

   omega = txDict(iTx)%omega
   period = txDict(iTx)%period
   !  set period, complete setup of 3D EM equation system
   i_omega_mu = cmplx(0.,isign*mu*omega,kind=selectedPrec)

   !  complete operator intialization, for this frequency
   !  call UpdateFreq(txDict(iTx)%omega)
   !  loop over polarizations
   do iMode = 1,2
      ! compute boundary conditions for polarization iMode
      !   Need to figure out mapping of conductivity needed here
      !    (have a pointer to modelParam in e0)!
      call BC_x0_WS(imode,period,e0%grid,e0%sigma, & 
		e0%pol(imode),b0%bc)
			write(*,'(a12,a3,a18,es12.6,a10,i2)') 'Solving the ','FWD', &
				' problem for freq ',omega/(2*PI),' & mode # ',imode
      call FWDsolve3D(b0,omega,e0%pol(imode))
   enddo
   !  set transmitter index for this solution 
   !          (sigma is set at initialization)
   e0%tx = iTx
   e0%period = period
   e0%omega = omega

   end subroutine fwdSolve

   !**********************************************************************
   subroutine sensSolve(iTx,iDT,FWDorADJ,comb,e)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine DOES NOT call UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.
   !  This final initialization step must (at present) be done by 
   !    calling fwdSolve before calling this routine.

   integer, intent(in)          	:: iTx,iDT
   character*3, intent(in)		:: FWDorADJ
   type(EMrhs), intent(inout)		:: comb
   type(EMsoln), intent(inout)		:: e

   ! local variables
   integer      			:: IER,iMode
   real(kind=selectedPrec) 		:: omega, period

   omega = txDict(iTx)%omega
   period = txDict(iTx)%period
   !  zero starting solution, solve for both modes
   call zero_EMsoln(e)
   do iMode = 1,2
      comb%b(imode)%adj = FWDorADJ
			write(*,'(a12,a3,a18,es12.6,a10,i2)') 'Solving the ',FWDorADJ, &
				' problem for freq ',omega/(2*PI),' & mode # ',imode
      call FWDsolve3d(comb%b(imode),omega,e%pol(imode))
   enddo

   e%tx = iTx
	 e%period = period
	 e%omega = omega

   end subroutine sensSolve   

!**********************************************************************
   subroutine create_EMsolnMTX(d,eAll)

      type(dvecMTX),intent(in)          :: d
      type(EMsolnMTX), intent(inout)     :: eAll

      !  local variables
      integer                           :: j
      character*80                      :: gridType

      eAll%nTx = d%nTx
      allocate(eAll%solns(d%nTx))
      eAll%allocated = .true.
      do j = 1,d%nTx
         call create_EMsoln(SolnRHS_grid,eAll%solns(j))
      enddo

   end subroutine create_EMsolnMTX

end module emsolver
