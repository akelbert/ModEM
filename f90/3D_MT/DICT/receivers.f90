! *****************************************************************************
module receivers
  ! This module contains the receiver dictionary (rxDict) for 3D MT

  use math_constants

  implicit none

  public			:: RXdictSetUp, deall_RXdict

  !  multiple receiver dictionaries can be defined, as
  !   different sorts of meta-data may be required for different data
  !   types; e.g., intersite TFs would require two locations.
  !   Note that the location must be specified relative to the same
  !     coordinate system used for the grid, with consistent units:
  !     in 3D MT, the origin of the grid is given relative to
  !     a refernce origin, with everything in meters.  The data site
  !     locations should also be given in meters, relative to this
  !     same physical origin.  All 3 components are required, to support
  !     observations anywhere within the solution domain.
  !  Additonal data types, to be used as elements of additional
  !    dictionaries can be added to accomodate additional data types
  type :: MTrx
     ! x gives location of EM measurements;
  	 ! x(1) points North, x(2) points East, x(3) points down
     real (kind=prec)                   ::  x(3)
     ! optional site ID; needed for input and output only
     character(80)                      ::  id=''
  end type MTrx

  ! receiver dictionary for 3D MT data will be an array of
  !  type MTrx (one element of the array for each site)
  !  Note that receiver dictionaries are only used inside the
  !    data functional module, and can thus be private
  type (MTrx), pointer, save, public, dimension(:) :: rxDict


Contains

  ! **************************************************************************
  ! Initializes and sets up receiver dictionary
  ! The reciever dictionary contains sparse vectors required
  ! for magnetic and electric field vector evaluation
  subroutine RXdictSetUp(nSites,siteLocations,siteIDs)

    integer, intent(in)	 		:: nSites
    real(kind=prec), intent(in)	:: siteLocations(3,nSites)
    character(*), intent(in), optional  :: siteIDs(nSites)
	character(3) :: id

    !  local variables
    integer		:: i

    allocate(rxDict(nSites))

    do i = 1,nSites
    	rxDict(i)%x = siteLocations(:,i)
		if (present(siteIDs)) then
			rxDict(i)%id = siteIDs(i)
		else
			write(id,'(i3.3)') i
			rxDict(i)%id = id
		end if
    enddo

  end subroutine RXdictSetUp

  ! **************************************************************************
  ! Cleans up and deletes receiver dictionary at end of program execution
  subroutine deall_rxDict()

	integer     :: istat

    if (associated(rxDict)) then
       deallocate(rxDict,STAT=istat)
    end if

  end subroutine deall_rxDict

end module receivers
