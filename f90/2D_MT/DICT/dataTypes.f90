! *****************************************************************************
module dataTypes
  ! This module contains the data type dictionary (typeDict) for 2D MT

  use math_constants
  use utilities

  implicit none

  type :: dataType

     !  stores information about the "data type"
     !   The following two attributes must be defined for all
     !    data types; these are accessed and used by the top-level
     !    inversion routines.
     logical                    :: isComplex = .false.
     logical                    :: calcQ = .false.
     !    Other attributes might be different (different number,
     !        different names, types, etc.) for  different applications.
     ! character(2)                :: mode = ''! = 'TE' or 'TM'
     character(80)               :: name = ''
     ! the number of components in the data type
     integer           			:: nComp
     !  could add rxDictNumber to keep track of reciever dictionary
     !  number used for this dataType (only 1 receiver dictionary now,
     !   so this is omitted)
     ! id(nComp) are the text comments that describe the components; these are
     ! initialized, but can be overwritten by the info from the data file.
     character(15), pointer, dimension(:) :: id

  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), pointer, save, public, dimension(:) :: typeDict

  ! add data types here ... this all needs work!
  integer, parameter    :: TE_Impedance = 1
  integer, parameter    :: TM_Impedance = 2

Contains

!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine setup_typeDict()

  	 integer     :: istat

     allocate(typeDict(2),STAT=istat)
     typeDict(TE_Impedance)%name = 'TE Impedance'
     typeDict(TE_Impedance)%isComplex = .true.
     typeDict(TE_Impedance)%calcQ     = .false.
     typeDict(TE_Impedance)%nComp     = 2
     allocate(typeDict(TE_Impedance)%id(2),STAT=istat)
     typeDict(TE_Impedance)%id(1)     = 'Re'
     typeDict(TE_Impedance)%id(2)     = 'Im'
     typeDict(TM_Impedance)%name = 'TM Impedance'
     typeDict(TM_Impedance)%isComplex = .true.
     typeDict(TM_Impedance)%calcQ     = .true.
     typeDict(TM_Impedance)%nComp     = 2
     allocate(typeDict(TM_Impedance)%id(2),STAT=istat)
     typeDict(TM_Impedance)%id(1)     = 'Re'
     typeDict(TM_Impedance)%id(2)     = 'Im'

  end subroutine setup_typeDict

! **************************************************************************
! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_typeDict()

	integer     :: j,istat

	if (associated(typeDict)) then

	   do j = 1,size(typeDict)
	      if (associated(typeDict(j)%id)) then
	         deallocate(typeDict(j)%id,STAT=istat)
	      end if
	   end do

       deallocate(typeDict,STAT=istat)

    end if

  end subroutine deall_typeDict

!**********************************************************************
! Computes the value by which the data must be multiplied to convert
! from the old units to the new units.
! The units may be any of the following.
! 1) SI units for E/B: [V/m]/[T] (used in ModEM code)
! 2) practical units for E/B: [mV/km]/[nT]
! 3) SI units for E/H: [V/m]/[A/m] = Ohm

  function ImpUnits(oldUnits,newUnits) result (SI_factor)

	character(*), intent(in)    :: oldUnits, newUnits
	real(kind=prec)             :: SI_factor
	! local
	real(kind=prec)             :: factor1, factor2

	! first convert the old units to [V/m]/[T]
	if (index(oldUnits,'[V/m]/[T]')>0) then
	   ! SI units for E/B
	   factor1 = ONE
	else if (index(oldUnits,'[mV/km]/[nT]')>0) then
	   ! practical units for E/B
	   factor1 = ONE * 1000.0
	else if ((index(oldUnits,'[V/m]/[A/m]')>0) .or. (index(oldUnits,'Ohm')>0)) then
	   ! SI units for E/H
	   factor1 = ONE * 1000.0 * 10000.0/(4*PI) ! approx. 796000.0
	else
	   call errStop('Unknown input units in ImpUnits: '//trim(oldUnits))
	end if

	! now convert [V/m]/[T] to the new units
	if (index(newUnits,'[V/m]/[T]')>0) then
	   ! SI units for E/B
	   factor2 = ONE
	else if (index(newUnits,'[mV/km]/[nT]')>0) then
	   ! practical units for E/B
	   factor2 = ONE / (1000.0)
	else if ((index(newUnits,'[V/m]/[A/m]')>0) .or. (index(newUnits,'Ohm')>0)) then
	   ! SI units for E/H
	   factor2 = ONE / (1000.0 * 10000.0/(4*PI))
	else
	   call errStop('Unknown output units in ImpUnits: '//trim(newUnits))
	end if

	SI_factor = factor1 * factor2

  end function ImpUnits

!**********************************************************************
! Figures out the data type from the component ids

  function ImpType(mode) result (dataType)

    character(*), intent(in)    :: mode
	integer	             	 	:: dataType

    if (index(mode,'TE')>0) then
        dataType = TE_Impedance
    else if (index(mode,'TM')>0) then
        dataType = TM_Impedance
    else
	   call errStop('Unknown mode in ImpType: '//trim(mode))
    end if

  end function ImpType

end module dataTypes
