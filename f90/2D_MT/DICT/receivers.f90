! *****************************************************************************
module receivers
  ! This module contains the receiver dictionary (rxDict) for 2D MT

  use math_constants

  implicit none

  public			:: RXdictSetUp, deall_RXdict

  type :: MTrx
     ! x gives location of EM measurements
     !  multiple receiver dictionaries can be defined, and
     !   different dictionaries can be used for different data types
     !  Additonal elements of MTrx data type can be added to
     !   accomodate additional data types
     real(kind=prec)			::  x(2)
  end type MTrx

  ! receiver dictionary for 2D MT data will be an array of
  !  type MTrx (one element of the array for each site)
  !  Two components of MTrx%x are position along the profile,
  !     and vertical position (generally on the surface)
  !  Note that the receiver dictionary is only used inside the
  !    data functional module
  type (MTrx), pointer, save, public, dimension(:) :: rxDict


Contains

  ! **************************************************************************
  ! Initializes and sets up receiver dictionary
  ! Now the reciever dictionary only contains the location of the point obs
  subroutine RXdictSetUp(nSites,siteLocations)
    !  siteLocatins(2,nSites) is array of measurement locations (x,z)
    !   corresponding to grid  (normally z = 0 for flat Earth surface)

    integer, intent(in)	 		:: nSites
    real(kind=prec), intent(in)	:: siteLocations(2,nSites)

    !  local variables
    integer                             :: i

    allocate(rxDict(nSites))
    do i = 1,nSites
       rxDict(i)%x = siteLocations(:,i)
    enddo

  end subroutine RXdictSetUp

  ! **************************************************************************
  ! Cleans up and deletes receiver dictionary at end of program execution
  subroutine deall_RXdict()

	integer     :: istat

    if (associated(rxDict)) then
       deallocate(rxDict,STAT=istat)
    end if

  end subroutine deall_RXdict


end module receivers
