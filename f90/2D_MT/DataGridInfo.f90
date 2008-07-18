module datagridinfo

use grid2d

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

  ! transmitter dictionary for MT data will be an array of
  ! type mt_forcing (one element  for each frequency)
  ! Perhaps this should be moved to EMsolver module (and be private
  !    to this module)
  type (MTtx), allocatable, save, public, dimension (:)   :: txDict

  type :: dataType

     !  stores information about the "data type"
     !   The following two attributes must be defined for all
     !    data types; these are accessed and used by the top-level
     !    inversion routines.
     logical                    :: isComplex = .false.
     logical                    :: calcQ = .false.
     !    Other attributes might be different (different number,
     !        different names, types, etc.) for  different applications.
     character*80               :: name = ''
     character*2                :: mode = ''! = 'TE' or 'TM'
     !  could add rxDictNumber to keep track of reciever dictionary
     !  number used for this dataType (only 1 receiver dictionary now,
     !   so this is omitted)

  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), allocatable, save, public, dimension(:) :: typeDict

  ! add data types here ... this all needs work!
  integer, parameter    :: TE_Impedance = 1
  integer, parameter    :: TM_Impedance = 2

  !  SolnRHS_grid is used to define grid parameters for DataFunc and
  !    EMsolver modules.  Make a copy of the numerical
  !   grid geometry parameters in this module at the start of
  !   the inversion.

  type(grid2d_t), target, save         :: SolnRHS_grid

contains

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
     enddo

  end subroutine TXdictSetUp

!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine TypeDictSetup()

     allocate(typeDict(2))
     typeDict(TE_Impedance)%name = 'TE Impedance'
     typeDict(TE_Impedance)%isComplex = .true.
     typeDict(TE_Impedance)%calcQ     = .false.
     typeDict(TE_Impedance)%mode     = TE
     typeDict(TM_Impedance)%name = 'TM Impedance'
     typeDict(TM_Impedance)%isComplex = .true.
     typeDict(TM_Impedance)%calcQ     = .true.
     typeDict(TM_Impedance)%mode     = TM

  end subroutine TypeDictSetUp

! **************************************************************************
  subroutine deall_Dict

     deallocate(txDict)
     deallocate(typeDict)

  end subroutine deall_Dict

end module datagridinfo
