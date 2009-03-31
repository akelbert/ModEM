! *****************************************************************************
module datafunc
  ! 2D MT data functionals ... now for TE + TM impedance ... could add
  !  This module contains 
  !  (1) receiver dictionary (rxDict); representation of
  !    data functionals for MT impedance are eliminated
  !  (2) transmitter dictionary (txDict), information about sources,
  !      now including TE/TM mode as well as frequency
  !  (3) data type dictionary typeDict
  !   DICTIONARIES ARE AN ARRAY OF reciever/transmitter descriptors
  !       THESE AR NOW PUBLIC ... this is all still under development
  !  (4) routines for evaluation of impedances, and ultimately other
  !	  interpretation parameters
  !  (5) routines to compute data functionals for linearized
  !       impedances,  and ultimately other interpretation paramters
  !   The idea: 
  !     -> rxDict and txDict are still private to this module
  !     -> first initialize txDict and rxDict by calling
  !          appropriate initialization/setup routines
  !           (included here, assuming simple cases only)
  !     -> data are stored in structures (defined in module data_vector)
  !        which contain indices into transmitter and receiver dictionaries
  !        in addition to actual data values.  These indices are used by
  !        the data function computation routines to compute predicted data.
  !
  ! This module is specific to 2D MT; similar modules would need to be written
  !     to implement data functionals for other problems

  use math_constants
  use utilities
  use interpeb        !  basic interpolation routines for 2D TE and TM
                        !    solution grids; allow computation of both E and B
                        !    at an arbitrary point in either grid
  use modelparameter, only:	rhoC	!  model parameterization
					!  dependent function needed
					!  for TM mode data
					!  functionals
  use solnrhs
 
  implicit none

  public			:: RXdictSetup, deall_RXdict
  !   Names of these routines must be as here, as these are called by
  !    top-level inversion routines
  public			:: nonLinDataFunc, linDataFunc
  public			:: EMSparseQtoModelParam

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
  type (MTrx), pointer, save, private, dimension(:) :: rxDict

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

!******************************************************************************
  subroutine nonLinDataFunc(ef,Sigma,iDT,iRX,Z)

  !   2D MT impedance; gets mode from solution
  ! given electric field solution  (stored as type cvector)
  !   + indices into receiver dictionary 
  ! compute complex scalar imepdance
  !   This now creates the required sparse vectors 
  !      (either for TE or TM, as appropriate) for impedance evaluation

  type (EMsoln_t), intent(in)			:: ef
  ! model parameter used to computed ef
  type (modelParam_t), intent(in)   :: Sigma
  ! indicies into data type and receiver dictionaries
  integer, intent(in)				:: iDt, iRX
  !  complex data returned; note that (a) this will always be complex
  !   (even if data is real ... then only real part is used by calling
  !    routine)  and (b) Z will always be an array in the calling
  !    routine (as in top-level inversion routines) and so needs to
  !    be treated as an array here, even if there is only one element.
  !   As an example: to add tippers to TE mode, dimension on Z 
  !     will have to be changed to 2!
  complex(kind=prec), intent(inout)	:: Z(1)

  !  local variables
  type(sparsevecc)		:: Lb,Le
  complex(kind=prec)	:: B,E
  real(kind=prec)	:: omega, x(2)
  logical			:: Conj_Case = .false.
  character*2			:: mode
  character*80			:: msg

  !  get mode, frequency for transmitter used to compute solution ef
  mode = ef%mode
  omega =  ef%omega
  ! get location from reciever dictionary
  x = rxDict(iRX)%x

  if(mode.eq.'TE') then
     ! electric field
     call NodeInterpSetup2D(ef%grid,x,mode,Le)
     ! magnetic field
     call BfromESetUp_TE(ef%grid,x,omega,Lb)


  elseif(mode.eq.'TM') then
     ! magnetic field
     call NodeInterpSetup2D(ef%grid,x,mode,Lb)
     ! electric field
     call EfromBSetUp_TM(ef%grid,x,omega,Sigma,Le)
  else
     call errStop('option not available in nonLinDataFunc')
  endif

  !  Using sparse vector representations of data functionals, 
  !          compute impedance
  E = dotProd_scvector(Le,ef%vec,Conj_case)  
  ! magnetic field
  B = dotProd_scvector(Lb,ef%vec,Conj_case)  
  
  ! impedance is trivial for 2D!
  Z = E/B

  ! clean up
  call deall_sparsevecc(Le)
  call deall_sparsevecc(Lb)

  end subroutine nonLinDataFunc

!****************************************************************************
  subroutine linDataFunc(e0,Sigma0,iDT,iRX,Lz,Qz)
  !  given input background electric field solution,
  !  index into receiver dictionary for a single site (iRX)
  !  compute sparse complex vector giving coefficients
  !  of linearized impedance functional (complex representation)
  !  For TM mode solution also returns sparse vector Q (model
  !     paramter space) for derivative of data functional with
  !     respect to model paramters; Q is not referenced for TE data
  !   NOTE: Lz and Qz have to be declared as arrays for 
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  
  !  electric field solutions are stored as type EMsoln
  type (EMsoln_t), intent(in)		:: e0
  ! model parameter used to computed e0
  type (modelParam_t), intent(in)   :: Sigma0
  ! indicies into data type and receiver dictionaries
  integer, intent(in)			:: iRX,iDT
  !  Lz, Qz will always be arrays in the calling
  !    routine (as in top-level inversion routines) and so need to
  !    be treated as arrays here, even if there is only one element.
  !   As an example: to add tippers to TE mode, dimension on LZ will
  !    have to be changed to 2. 
  type(EMsparse_t), intent(inout)		:: Lz(1), Qz(1)

  !  local variables
  complex (kind=prec)		:: B,E,c_E,c_B
  type(sparsevecc)			:: Le,Lb
  real (kind=prec)		:: x(2), omega
  character*2				:: mode
  logical				:: Conj_case = .false.

  !  get mode, frequency for transmitter used to compute solution ef
  mode = e0%mode
  omega =  e0%omega
  ! get location from reciever dictionary
  x = rxDict(iRX)%x

  !   evaluate E, B, for background solution
  if(mode.eq.'TE') then
     ! electric field
     call NodeInterpSetup2D(e0%grid,x,mode,Le)
     ! magnetic field
     call BfromESetUp_TE(e0%grid,x,omega,Lb)
  elseif(mode.eq.'TM') then
     ! magnetic field
     call NodeInterpSetup2D(e0%grid,x,mode,Lb)
     ! electric field
     call EfromBSetUp_TM(e0%grid,x,omega,Sigma0,Le,e0%vec,Qz(1)%L)
  else
     call errStop('option not available in linDataFunc')
  endif
  
  !  Compute electric, magnetic field for background soln.
  E = dotProd_scvector(Le,e0%vec,Conj_case)  
  B = dotProd_scvector(Lb,e0%vec,Conj_case)  

  !  compute sparse vector representations of linearized 
  !    impedance functional; coefficients depend on E, B
  c_E = C_ONE/B
  c_B = -E/(B*B)
  !  Nominally, Lz = c_E*Le + c_B*Lb
  ! However, note that Lz is of type EMsparse (here just a wrapped
  !    version of a sparsevecc object)
  call linComb_sparsevecc(Le,c_E,Lb,c_B,Lz(1)%L) 

  if(mode .eq. 'TM') then
     !  also need to multiply parameter space sparse vector
     !    (derivative of Le coefficients wrt cell resistivities)
     !   by 1/B  (actually 1/B == 1 !)
     Qz(1)%L%C = c_E*Qz(1)%L%C
  endif

  ! clean up
  call deall_sparsevecc(Le)
  call deall_sparsevecc(Lb)

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

   complex(kind=prec),intent(in)	:: cs
   type (EMsparse_t), intent(in)                  :: Q
   type (modelParam_t), intent(in)                :: sigma0
   type (modelParam_t), intent(inout)             :: dsigmaReal
   type (modelParam_t), intent(inout),optional    :: dsigmaImag
   
   !  local variables
   type(sparsevecc)	:: Ltemp

   if(present(dSigmaImag)) then
         call scMult_sparsevecc(cs,Q%L,Ltemp)
         call QtoModelParam(Ltemp,sigma0,dsigmaReal,dSigmaImag)
   else
         call scMult_sparsevecc(cs,Q%L,Ltemp)
         call QtoModelParam(Ltemp,sigma0,dsigmaReal)
   endif

   call deall_sparsevecc(Ltemp)

   end subroutine EMSparseQtoModelParam
   

end module dataFunc
