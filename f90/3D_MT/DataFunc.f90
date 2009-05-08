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

  use emfieldinterp
  use solnrhs
  use receivers
  use transmitters
  use datatypes

  implicit none

  !   Names of these routines must be as here, as these are called by
  !    top-level inversion routines
  public                        :: nonLinDataFunc, linDataFunc


Contains

!******************************************************************************
  subroutine nonLinDataFunc(ef,Sigma,iDT,iRX,Z,Binv)
  ! given electric field solutions (both modes--and note
  !    that the solution knows about the transmitter used),
  ! and indices into data types and receiver dictionaries for one
  ! data vector compute the complex impedance tensor.
  ! Binv is optional output argument, used needed for linearized
  ! impedance calculation in this module (not used by higher levels)

  implicit none
  type (EMsoln_t), intent(in)		:: ef
  type (modelParam_t), intent(in) :: Sigma ! used to compute ef
  integer, intent(in)			:: iDT
  integer, intent(in) 			:: iRX
  complex(kind=prec), intent(inout)	:: Z(*)

  ! Definition of the impedance elements:
  !   iDT=Full_Impedance
  ! Z(1) = Zxx; Z(2) = Zxy; Z(3) = Zyx; Z(4) = Zyy
  !   iDT=Impedance_Plus_Hz
  ! Z(5) = Zzx; Z(6) = Zzy in addition
  !   iDT=Off_Diagonal_Impedance
  !  Z(1) = Zxy, Z(2) = Zyx

  !  optional argument, useful for linearized impedance
  complex(kind=prec), intent(out), optional	:: Binv(2,2)

  !  local variables
  integer			:: iMode, i,j,xyz,ij
  real(kind=prec)	:: omega,x(3)
  complex(kind=prec)	:: BB(3,2),EE(2,2)
  complex(kind=prec)	:: det,i_omega,ctemp
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

  type (EMsoln_t), intent(in)		   :: e0
  type (modelParam_t), intent(in)  :: Sigma0
  integer, intent(in)			   :: iDT, iRX
  !   NOTE: Lz and Qz have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE:  in principal the comparable input arguments in
  !        the 2D program should also be of type EMsparse!
  type(EMsparse_t), intent(inout)		:: L(*),Q(*)

  !  local variables
  complex(kind=prec)	:: Binv(2,2)
  complex (kind=prec)	:: Z(6), i_omega,c1,c2
  real(kind=prec)	:: x(3),omega
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

end module dataFunc
