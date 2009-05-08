! *****************************************************************************
module transmitters
  ! This module contains the transmitter dictionary (txDict) for 3D MT

  use math_constants

  implicit none

  public			:: TXdictSetUp, deall_TXdict

  type :: MTtx
     !  An MT source is defined by frequency and boundary conditions
     !   at present there does not seem to be much need for BC info ... add
     !    if needed.  Other sorts of EM data may have more
     !    complex tx descriptions
     character(2)                :: mode = ''! = 'TE' or 'TM'
     ! angular frequency (radians/sec), and for convenience period (s)
     real(kind=prec)            :: omega = R_ZERO
     real(kind=prec)            :: period = R_ZERO
     ! index number to frequency/ period in solution file
     integer                    :: iPer
   end type MTtx

   ! transmitter dictionary for MT data will be an array of
   ! type mt_forcing (one element  for each frequency)
   ! Perhaps this should be moved to EMsolver module (and be private
   !    to this module)
   type (MTtx), pointer, save, public, dimension (:)   :: txDict

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

end module transmitters
