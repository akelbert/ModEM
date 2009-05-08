! *****************************************************************************
module dataTypes
  ! This module contains the data type dictionary (typeDict) for 2D MT

  use math_constants

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
     !  could add rxDictNumber to keep track of reciever dictionary
     !  number used for this dataType (only 1 receiver dictionary now,
     !   so this is omitted)

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
  subroutine TypeDictSetup()

     allocate(typeDict(2))
     typeDict(TE_Impedance)%name = 'TE Impedance'
     typeDict(TE_Impedance)%isComplex = .true.
     typeDict(TE_Impedance)%calcQ     = .false.
     typeDict(TM_Impedance)%name = 'TM Impedance'
     typeDict(TM_Impedance)%isComplex = .true.
     typeDict(TM_Impedance)%calcQ     = .true.

  end subroutine TypeDictSetUp

! **************************************************************************
! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_typeDict()

	integer     :: istat

    if (associated(typeDict)) then
       deallocate(typeDict,STAT=istat)
    end if

  end subroutine deall_typeDict

end module dataTypes
