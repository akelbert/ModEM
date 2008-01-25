module datagridinfo

use math_constants
use grid3d


 type :: MTtx
     !  An MT source is defined by frequency and boundary conditions
     !   at present there does not seem to be much need for BC info ... add
     !    if needed.  Other sorts of EM data may have more
     !    complex tx descriptions
     ! angular frequency (radians/sec), and for convenience period (s)
     real(kind=selectedPrec)            :: omega = R_ZERO
     real(kind=selectedPrec)            :: period = R_ZERO
     ! index number to frequency/ period in solution file
     integer                    :: iPer
  end type MTtx

  ! transmitter dictionary txDict for 3D-MT data will be an array of
  ! type MTtx (one element  for each frequency)
  ! Perhaps this should be moved to EMsolver module (and be private
  !    to that module?)   
  ! NOTE: could have multiple transmitter dictionaries, each of which
  !    could constist of elements of different types; nothing about
  !    the dictionary or the elements that it consists of is used
  !    in higher level routines
  type (MTtx), pointer, save, public, dimension (:)   :: txDict

  type :: dataType
!
     !  stores information about the "data type" -- which could include
     !   information that is relevant to transmitter, receiver, or both
     !   E.g., for 2D MT, "mode" (TE/TM) is relevant to both receiver and
     !    transmitter.  For 3D MT, this can be used to distinguish between
     !    full impedance tensors, impedance tensors+vertical field TFs,
     !    off-diagonal impedance only, interstation TFs, etc.
     !    Different data types may correspond to different EM solutions
     !     (TE vs. TM; joint inversion of data from multiple geophysical
     !        techniques), or not (different receiver configurations in 3D MT).
     !    In some cases, multiple transmitter dictionaries may be needed (and
     !     which to use would be determined in, e.g., EMsolve).
     !    Similarly, in some cases multiple receiver dictionaries may or may not
     !     be required when there are multiple dataTypes
     !       (e.g., interstation TFs require more information, to
     !        specify a pair of site locations; adding more local TF components
     !         would not require multiple receiver dictionaries)
     !     dictionaries may be required 
     !
     !   The following two attributes (essentially receiver attributes)
     !      must be defined for all data types.  These are accessed 
     !       and used by the top-level inversion routines:
     logical                    :: isComplex = .false.
     logical                    :: calcQ = .false.

     character*80               :: name = ''
     !  for 3D MT data types will initially be used only to distinguish
     !    between local transfer function types 
     !         (later maybe add interstation TFs)
     integer			:: tfType, nComp
!     
  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), pointer, save, public, dimension(:) :: typeDict

  ! obviously much more needs to be done to streamline addition of data types ... 
  integer, parameter    :: Full_Impedance = 1
  integer, parameter    :: Impedance_Plus_Hz = 2
  integer, parameter    :: Off_Diagonal_Impedance = 3
  

  !  SolnRHS_grid is used to define grid parameters for DataFunc and
  !    EMsolver modules.  Make a copy of the numerical
  !   grid geometry parameters in this module at the start of
  !   the inversion.

  type(grid3d_t), target, save         :: SolnRHS_grid

contains

!**********************************************************************
    subroutine set_SolnRHS_grid(grid)
!    Call this routine to set basic grid geometry parameters
!       before using any other routines in this module (most depend
!       on saved SolnRHS_grid to define grid geometry)

       type (grid3d_t), intent(in)     :: grid

       SolnRHS_grid = grid

    end subroutine set_SolnRHS_grid
!**********************************************************************

! Initializes and sets up transmitter dictionary for MT,
!  This is just a simple example of a routine for setting up the TX
!   dictionary; In this example we assume that there are nPer periods
!   for either TE or TM modes, or for both.
!  NOTE:   If TXdict if public, there is no reason for this to
!    be part of this module (but I leave it here for now!)

  subroutine TXdictSetUp(nTx,Periods)

     integer, intent(in)         :: nTx
     real(kind=selectedPrec), intent(in)          :: periods(nTx)
 
     ! local variables
     integer                     :: iTx

     allocate(txDict(nTx))
     do iTx = 1, nTx
        txDict(iTx)%period = Periods(iTx)
        txDict(iTx)%omega = (2*PI)/ txDict(iTx)%period
     enddo

  end subroutine TXdictSetUp

!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine TypeDictSetup()

     allocate(typeDict(3))

     typeDict(Full_Impedance)%name = 'Full Impedance'
     typeDict(Full_Impedance)%isComplex = .true.
     typeDict(Full_Impedance)%calcQ     = .false.
     typeDict(Full_Impedance)%tfType     = Full_Impedance
     typeDict(Full_Impedance)%nComp     = 8

     typeDict(Impedance_Plus_Hz)%name = 'Full Impedance Plus Hz'
     typeDict(Impedance_Plus_Hz)%isComplex = .true.
     typeDict(Impedance_Plus_Hz)%calcQ     = .false.
     typeDict(Impedance_Plus_Hz)%tfType     = Impedance_Plus_Hz
     typeDict(Impedance_Plus_Hz)%nComp     = 12

     typeDict(Off_Diagonal_Impedance)%name = 'Off Diagonal Impedance'
     typeDict(Off_Diagonal_Impedance)%isComplex = .true.
     typeDict(Off_Diagonal_Impedance)%calcQ     = .false.
     typeDict(Off_Diagonal_Impedance)%tfType     = Off_Diagonal_Impedance
     typeDict(Off_Diagonal_Impedance)%nComp     = 4

  end subroutine TypeDictSetUp

! **************************************************************************
  subroutine deall_Dict

     deallocate(txDict)
     deallocate(typeDict)

  end subroutine deall_Dict

end module datagridinfo
