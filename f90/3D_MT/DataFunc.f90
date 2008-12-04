! *****************************************************************************
module dataFunc
  ! 3D MT data functionals:  This module contains
  !  (1) receiver dictionary (rxDict), containing sparse representation of
  !    data functionals required to define MT impedance or other interpretation
  !    parameters, including linear and linearized functionals.
  !  (2) routines for evaluation of impedances, and ultimately other
  !       interpretation parameters
  !  (3) routines to compute data functionals for linearized
  !       impedances,  and ultimately other interpretation paramters
  !   The idea:
  !     -> first initialize rxDict by calling
  !          appropriate initialization/setup routines
  !           (included here, assuming simple cases only)
  !     -> data are stored in structures (defined in module data_vector)
  !        which contain indices into transmitter and receiver dictionaries
  !        in addition to actual data values.  These indices are used by
  !        the data function computation routines to compute predicted data.
  !
  !
  !  This module is specific to 3D MT; similar modules will need to be written
  !     to implement data functionals for other problems

  use interpeb
  use solnrhs
  use math_constants

  implicit none

  public			:: RXdictSetup, deall_RXdict
  !   Names of these routines must be as here, as these are called by
  !    top-level inversion routines
  public                        :: nonLinDataFunc, linDataFunc
  public                        :: EMSparseQtoModelParam

  type :: MTrx
     ! x gives location of EM measurements
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
     real(kind=selectedPrec)                    ::  x(3)
     ! optional site ID; needed for input and output only
     character(80)                              ::  id=''
  end type MTrx

  ! receiver dictionary for 3D MT data will be an array of
  !  type MTrx (one element of the array for each site)
  !  Note that receiver dictionaries are only used inside the
  !    data functional module, and can thus be private
  type (MTrx), pointer, save, private, dimension(:) :: rxDict

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

     character(80)              :: name = ''
     !  for 3D MT data types will initially be used only to distinguish
     !    between local transfer function types
     !         (later maybe add interstation TFs)
     integer			:: tfType, nComp

     ! AK: id(nComp) are the text comments that describe the components
     !character(80), pointer, dimension(:)   :: id
     ! however, in order to implement this everything has to be done differently.
     ! I don't like the whole dictionary business, it's just not logical... one
     ! should be able to set up the dvec uniquely from the dictionaries and the
     ! data values; the way this is done is backwards, and requires supplementary
     ! temporary information in the upper level routines to set everything up
     ! or to gather the information together, which really defeats the whole purpose
     ! of this headache. It just hasn't been thought through. For now, we stick
     ! with the typeDict, and assume our data is in the correct order etc.
     character(200)             :: comment = ''
     ! AK: let us also store the units for a dataType here (obviously, the way
     ! this is currently formulated, Impedance_Plus_Hz will be partly non-dimensional)
     character(80)              :: units = ''
!
  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), pointer, save, public, dimension(:) :: typeDict

  ! obviously much more needs to be done to streamline addition of data types ...
  integer, parameter    :: Full_Impedance = 1
  integer, parameter    :: Impedance_Plus_Hz = 2
  integer, parameter    :: Off_Diagonal_Impedance = 3


Contains

  ! **************************************************************************
  ! Initializes and sets up receiver dictionary
  ! The reciever dictionary contains sparse vectors required
  ! for magnetic and electric field vector evaluation
  subroutine RXdictSetUp(nSites,siteLocations,siteIDs)

    integer, intent(in)	 		:: nSites
    real(kind=selectedPrec), intent(in)	:: siteLocations(3,nSites)
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

!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine TypeDictSetup()

  	 integer     :: istat

     allocate(typeDict(3),STAT=istat)

     typeDict(Full_Impedance)%name = 'Full Impedance'
     typeDict(Full_Impedance)%isComplex = .true.
     typeDict(Full_Impedance)%calcQ     = .false.
     typeDict(Full_Impedance)%tfType     = Full_Impedance
     typeDict(Full_Impedance)%nComp     = 8
     typeDict(Full_Impedance)%units   = '[V/m]/[T]'
     write(typeDict(Full_Impedance)%comment,'(8a15)') &
     	'Re(Zxx)','Im(Zxx)','Re(Zxy)','Im(Zxy)','Re(Zyx)','Im(Zyx)','Re(Zyy)','Im(Zyy)'

     typeDict(Impedance_Plus_Hz)%name = 'Full Impedance Plus Hz'
     typeDict(Impedance_Plus_Hz)%isComplex = .true.
     typeDict(Impedance_Plus_Hz)%calcQ     = .false.
     typeDict(Impedance_Plus_Hz)%tfType     = Impedance_Plus_Hz
     typeDict(Impedance_Plus_Hz)%nComp     = 12
     typeDict(Impedance_Plus_Hz)%units   = '[V/m]/[T]'
     write(typeDict(Impedance_Plus_Hz)%comment,'(12a15)') &
     	'Re(Zxx)','Im(Zxx)','Re(Zxy)','Im(Zxy)','Re(Zyx)','Im(Zyx)','Re(Zyy)','Im(Zyy)','Re(Tx)','Im(Tx)','Re(Ty)','Im(Ty)'

     typeDict(Off_Diagonal_Impedance)%name = 'Off Diagonal Impedance'
     typeDict(Off_Diagonal_Impedance)%isComplex = .true.
     typeDict(Off_Diagonal_Impedance)%calcQ     = .false.
     typeDict(Off_Diagonal_Impedance)%tfType     = Off_Diagonal_Impedance
     typeDict(Off_Diagonal_Impedance)%nComp     = 4
     typeDict(Off_Diagonal_Impedance)%units  = '[V/m]/[T]'
     write(typeDict(Off_Diagonal_Impedance)%comment,'(4a15)') &
     	'Re(Zxy)','Im(Zxy)','Re(Zyx)','Im(Zyx)'

  end subroutine TypeDictSetUp

! **************************************************************************
! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_typeDict()

	integer     :: istat

    if (associated(typeDict)) then
       deallocate(typeDict,STAT=istat)
    end if

  end subroutine deall_typeDict

!******************************************************************************
  subroutine nonLinDataFunc(ef,Sigma,iDT,iRX,Z,Binv)
  ! given electric field solutions (both modes--and note
  !    that the solution knows about the transmitter used),
  ! and indices into data types and receiver dictionaries for one
  ! data vector compute the complex impedance tensor.
  ! Binv is optional output argument, used needed for linearized
  ! impedance calculation in this module (not used by higher levels)

  implicit none
  type (EMsoln), intent(in)		:: ef
  type (modelParam_t), intent(in) :: Sigma ! used to compute ef
  integer, intent(in)			:: iDT
  integer, intent(in) 			:: iRX
  complex(kind=selectedPrec), intent(inout)	:: Z(*)

  ! Definition of the impedance elements:
  !   iDT=Full_Impedance
  ! Z(1) = Zxx; Z(2) = Zxy; Z(3) = Zyx; Z(4) = Zyy
  !   iDT=Impedance_Plus_Hz
  ! Z(5) = Zzx; Z(6) = Zzy in addition
  !   iDT=Off_Diagonal_Impedance
  !  Z(1) = Zxy, Z(2) = Zyx

  !  optional argument, useful for linearized impedance
  complex(kind=selectedPrec), intent(out), optional	:: Binv(2,2)

  !  local variables
  integer			:: iMode, i,j,xyz,ij
  real(kind=selectedPrec)	:: omega,x(3)
  complex(kind=selectedPrec)	:: BB(3,2),EE(2,2)
  complex(kind=selectedPrec)	:: det,i_omega,ctemp
  type(sparsevecC)		:: Lex,Ley,Lbx,Lby,Lbz
  logical			:: ComputeHz

  !  probably should dependence on omega into BinterpSetup, as in 2D!
  omega = ef%omega
  x = rxDict(iRX)%x
  ComputeHz = (typeDict(iDT)%tfType.EQ.Impedance_Plus_Hz)

  ! First set up interpolation functionals for Ex, Ey, Bx, By, Bz
  xyz = 1
  call EinterpSetUp(ef%grid,x,xyz,Lex)
  xyz = 2
  call EinterpSetUp(ef%grid,x,xyz,Ley)
  xyz = 1
  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
  xyz = 2
  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
  if(ComputeHz) then
     xyz = 3
     call BfromESetUp(ef%grid,omega,x,xyz,Lbz)
  endif

  ! loop over modes
  do iMode = 1,2
      ! electric fields
      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
      ! magnetic fields
      ! NOTE: sparse vectors containing data functionals for magnetic
      ! fields ARE NOW divided by i_omega
      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
      if(ComputeHz) then
        BB(3,iMode) = dotProd_noConj_scvector_f(Lbz,ef%pol(iMode))
      endif
  enddo

  ! compute impedances  (invert horizontal H matrix
  ! using Kramer's rule)
  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
  ctemp = BB(1,1)
  BB(1,1) =  BB(2,2)/det
  BB(2,2) =  ctemp/det
  BB(1,2) = -BB(1,2)/det
  BB(2,1) = -BB(2,1)/det

  selectcase (iDT)
     case(Full_Impedance)
        do j = 1,2
           do i = 1,2
              ij = 2*(i-1)+j
              Z(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
           enddo
        enddo
     case(Impedance_Plus_Hz)
        do j = 1,2
           do i = 1,2
              ij = 2*(i-1)+j
              Z(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
           enddo
	   !  vertical field TF
           ij = 4+j
           Z(ij) = BB(3,1)*BB(1,j)+BB(3,2)*BB(2,j)
        enddo
     case(Off_Diagonal_Impedance)
        Z(1) = EE(1,1)*BB(1,2)+EE(1,2)*BB(2,2)
        Z(2) = EE(2,1)*BB(1,1)+EE(2,2)*BB(2,1)
  end select

  if(present(Binv)) then
     Binv = BB(1:2,:)
  endif

  ! clean up
  call deall_sparsevecc(Lex)
  call deall_sparsevecc(Ley)
  call deall_sparsevecc(Lbx)
  call deall_sparsevecc(Lby)
  call deall_sparsevecc(Lbz)

  end subroutine nonLinDataFunc
!
!****************************************************************************
  subroutine linDataFunc(e0,Sigma0,iDT,iRX,L,Q)
  !  given input background electric field solution (both modes; e0),
  !  indices into data type/receiver dictionaries
  !  compute array of sparse complex vectors giving coefficients
  !  of linearized impedance functional (complex representation)
  !  For data types with attribute computeQ = .true. also returns
  !     derivative of data functional with respect to model paramters
  !         (NOT YET IMPLEMENTED FOR 3D MT!!!!)

  type (EMsoln), intent(in)		   :: e0
  type (modelParam_t), intent(in)  :: Sigma0
  integer, intent(in)			   :: iDT, iRX
  !   NOTE: Lz and Qz have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE:  in principal the comparable input arguments in
  !        the 2D program should also be of type EMsparse!
  type(EMsparse), intent(inout)		:: L(*),Q(*)

  !  local variables
  complex(kind=selectedPrec)	:: Binv(2,2)
  complex (kind=selectedPrec)	:: Z(6), i_omega,c1,c2
  real(kind=selectedPrec)	:: x(3),omega
  type(sparsevecc)		:: L1
  integer			:: i,j,k,nComp,IJ(2,6),xyz,n
  type(sparsevecC)		:: Lex,Ley,Lbx,Lby,Lbz
  logical			:: ComputeHz

  omega = e0%omega
  x = rxDict(iRX)%x

  !  set up which components are needed,  ... and ! evaluate
  !   impedance, Binv for background solution
  !          ... appear to need full impedance for offdiagonal
  select case(iDT)
     case(Full_Impedance)
        nComp = 4;
        ComputeHz = .false.
        do j = 1,2
           do i = 1,2
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
           enddo
        enddo
        Call nonlinDataFunc(e0,Sigma0,Full_Impedance,iRX,Z,Binv)
     case(Impedance_Plus_Hz)
        nComp = 6;
        ComputeHz = .true.
        do j = 1,2
           do i = 1,3
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
           enddo
        enddo
        Call nonlinDataFunc(e0,Sigma0,Impedance_Plus_Hz,iRX,Z,Binv)
     case(Off_Diagonal_Impedance)
        nComp = 2
        ComputeHz = .false.
        IJ(1,1) = 1
        IJ(2,1) = 2
        IJ(1,2) = 2
        IJ(2,2) = 1
        Call nonlinDataFunc(e0,Sigma0,Full_Impedance,iRX,Z,Binv)
  endselect

  ! Then set up interpolation functionals for Ex, Ey, Bx, By, Bz
  xyz = 1
  call EinterpSetUp(e0%grid,x,xyz,Lex)
  xyz = 2
  call EinterpSetUp(e0%grid,x,xyz,Ley)
  xyz = 1
  call BfromESetUp(e0%grid,omega,x,xyz,Lbx)
  xyz = 2
  call BfromESetUp(e0%grid,omega,x,xyz,Lby)
  if(ComputeHz) then
     xyz = 3
     call BfromESetUp(e0%grid,omega,x,xyz,Lbz)
  endif

  !  compute sparse vector representations of linearized functionals
  do n = 1,nComp
     !  i runs over rows of TF matrix, j runs over columns of TF
     i = IJ(1,n)
     j = IJ(2,n)
     c1 = Z(2*(i-1)+1)
     c2 = Z(2*(i-1)+2)
     Call linComb_sparsevecc(Lbx,c1,Lby,c2,L1)
     do k = 1,2
        !  k defines which mode the linearized functional is
        !   to be applied to
        c1 = Binv(k,j)
        c2 = -c1
        if(i.eq.1) then
           !  component in x row of impedance tensor
           Call linComb_sparsevecc(Lex,c1,L1,c2,L(n)%L(k))
        elseif(i.eq.2) then
           !  component in y row of impedance tensor
           Call linComb_sparsevecc(Ley,c1,L1,c2,L(n)%L(k))
        else
           !  component in Bz row (vertical field TF)
           Call linComb_sparsevecc(Lbz,c1,L1,c2,L(n)%L(k))
        endif
     enddo
  enddo

  ! clean up
  call deall_sparsevecc(L1)
  call deall_sparsevecc(Lex)
  call deall_sparsevecc(Ley)
  call deall_sparsevecc(Lbx)
  call deall_sparsevecc(Lby)
  call deall_sparsevecc(Lbz)

  end subroutine linDataFunc
!****************************************************************************
 subroutine EMSparseQtoModelParam(cs,Q,sigma0,dsigmaReal,dSigmaImag)

   !   adds cs*Q to (dsigmaReal,dSigmaImag)
   !   cs is a complex constant, Q is a sparse scalar field defined on
   !     grid cells (but represented with EMsparse data object ... the
   !     xyz component indicators are ignored).  dsigmaReal/dsigmaImag
   !     are used to collect sums of real and imaginary parts; dsigmaImag
   !     is optional.
   !
   !  Mostly this is just a wrapper for QtoModelParam.
   !   QtoModelParam has to be part of the modelParameter module,
   !   in order to keep modelParam attributes private.
   !  To avoid reference to objects below the level of generic SolnRHS
   !   module objects in higher level inversion routines it is necessary
   !   to wrap QtoModelParam with this simple routine.
   !
   !  Now, this also is used to multiply Q by cs before adding.

   complex(kind=selectedPrec),intent(in)	:: cs
   type (EMsparse), intent(in)                  :: Q
   type (modelParam_t), intent(in)                :: sigma0
   type (modelParam_t), intent(inout)             :: dsigmaReal
   type (modelParam_t), intent(inout),optional    :: dsigmaImag

   !  local variables
   integer		:: k
   type(sparsevecc)	:: Ltemp

   if(present(dSigmaImag)) then
      do k = 1,Q%nPol
         call scMult_sparsevecc(cs,Q%L(k),Ltemp)
         call QtoModelParam(Ltemp,sigma0,dsigmaReal,dSigmaImag)
      enddo
   else
      do k = 1,Q%nPol
         call scMult_sparsevecc(cs,Q%L(k),Ltemp)
         call QtoModelParam(Ltemp,sigma0,dsigmaReal)
      enddo
   endif

   call deall_sparsevecc(Ltemp)

   end subroutine EMSparseQtoModelParam

end module dataFunc
