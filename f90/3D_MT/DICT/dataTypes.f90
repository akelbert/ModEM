! *****************************************************************************
module dataTypes
  ! This module contains the data type dictionary (typeDict) for 3D MT

  use math_constants

  implicit none

  public			:: typeDictSetUp, deall_typeDict

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
  type :: dataType

     ! The following two attributes (essentially receiver attributes)
     ! must be defined for all data types.  These are accessed
     ! and used by the top-level inversion routines.
     !
     ! AK: It is debatable whether this is a good idea... I suggest rethinking:
     ! isComplex can easily be stored in dataVec, and calcQ ... well, needs to
     ! be done differently! I think that the sole purpose of the dataType
     ! should be to identify the mapping from the EM solution to the data vector.
     logical           :: isComplex = .false.
     logical           :: calcQ = .false.

     ! a user-friendly description of this data type
     character(80)     :: name = ''

     ! for 3D MT data types will initially be used only to distinguish
     ! between local transfer function types (later maybe add interstation TFs)
     integer		   :: tfType

     ! the number of components in the data type
     integer           :: nComp

     ! id(nComp) are the text comments that describe the components; these are
     ! initialized, but can be overwritten by the info from the data file.
     !
     ! AK: I suspect the way this is currently set up implies that the default order
     ! of components must be preserved! This needs more attention; until then,
     ! the user should specify the components in the default order.
     character(15), pointer, dimension(:) :: id

     ! the units of the data type
     ! (obviously, the way this is currently formulated, Impedance_Plus_Hz will
     ! be partly non-dimensional; but this works for practical purposes)
     character(80)     :: units = ''

  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), pointer, save, public, dimension(:) :: typeDict

  integer, parameter   :: Full_Impedance = 1
  integer, parameter   :: Impedance_Plus_Hz = 2
  integer, parameter   :: Off_Diagonal_Impedance = 3


Contains


!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine TypeDictSetUp()

  	 integer     :: istat

     allocate(typeDict(3),STAT=istat)

     typeDict(Full_Impedance)%name = 'Full Impedance'
     typeDict(Full_Impedance)%isComplex = .true.
     typeDict(Full_Impedance)%calcQ     = .false.
     typeDict(Full_Impedance)%tfType     = Full_Impedance
     typeDict(Full_Impedance)%units   = '[V/m]/[T]'
     typeDict(Full_Impedance)%nComp     = 8
     allocate(typeDict(Full_Impedance)%id(8),STAT=istat)
     typeDict(Full_Impedance)%id(1) = 'Re(Zxx)'
     typeDict(Full_Impedance)%id(2) = 'Im(Zxx)'
     typeDict(Full_Impedance)%id(3) = 'Re(Zxy)'
     typeDict(Full_Impedance)%id(4) = 'Im(Zxy)'
     typeDict(Full_Impedance)%id(5) = 'Re(Zyx)'
     typeDict(Full_Impedance)%id(6) = 'Im(Zyx)'
     typeDict(Full_Impedance)%id(7) = 'Re(Zyy)'
     typeDict(Full_Impedance)%id(8) = 'Im(Zyy)'

     typeDict(Impedance_Plus_Hz)%name = 'Full Impedance Plus Hz'
     typeDict(Impedance_Plus_Hz)%isComplex = .true.
     typeDict(Impedance_Plus_Hz)%calcQ     = .false.
     typeDict(Impedance_Plus_Hz)%tfType     = Impedance_Plus_Hz
     typeDict(Impedance_Plus_Hz)%units   = '[V/m]/[T]'
     typeDict(Impedance_Plus_Hz)%nComp     = 12
     allocate(typeDict(Full_Impedance)%id(12),STAT=istat)
     typeDict(Full_Impedance)%id(1)  = 'Re(Zxx)'
     typeDict(Full_Impedance)%id(2)  = 'Im(Zxx)'
     typeDict(Full_Impedance)%id(3)  = 'Re(Zxy)'
     typeDict(Full_Impedance)%id(4)  = 'Im(Zxy)'
     typeDict(Full_Impedance)%id(5)  = 'Re(Zyx)'
     typeDict(Full_Impedance)%id(6)  = 'Im(Zyx)'
     typeDict(Full_Impedance)%id(7)  = 'Re(Zyy)'
     typeDict(Full_Impedance)%id(8)  = 'Im(Zyy)'
     typeDict(Full_Impedance)%id(9)  = 'Re(Tx)'
     typeDict(Full_Impedance)%id(10) = 'Im(Tx)'
     typeDict(Full_Impedance)%id(11) = 'Re(Ty)'
     typeDict(Full_Impedance)%id(12) = 'Im(Ty)'

     typeDict(Off_Diagonal_Impedance)%name = 'Off Diagonal Impedance'
     typeDict(Off_Diagonal_Impedance)%isComplex = .true.
     typeDict(Off_Diagonal_Impedance)%calcQ     = .false.
     typeDict(Off_Diagonal_Impedance)%tfType     = Off_Diagonal_Impedance
     typeDict(Off_Diagonal_Impedance)%units  = '[V/m]/[T]'
     typeDict(Off_Diagonal_Impedance)%nComp     = 4
     allocate(typeDict(Full_Impedance)%id(4),STAT=istat)
     typeDict(Full_Impedance)%id(1) = 'Re(Zxy)'
     typeDict(Full_Impedance)%id(2) = 'Im(Zxy)'
     typeDict(Full_Impedance)%id(3) = 'Re(Zyx)'
     typeDict(Full_Impedance)%id(4) = 'Im(Zyx)'

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
