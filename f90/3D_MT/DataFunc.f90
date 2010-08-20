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
  use SolnSpace
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
  type (solnVector_t), intent(in)		:: ef
  type (modelParam_t), intent(in) :: Sigma ! used to compute ef
  integer, intent(in)			:: iDT
  integer, intent(in) 			:: iRX
  complex(kind=prec), intent(inout)	:: Z(*)

  ! Definition of the impedance elements:
  !   iDT=Full_Impedance
  ! 			Z(1) = Zxx; Z(2) = Zxy; Z(3) = Zyx; Z(4) = Zyy
  !   iDT=Impedance_Plus_Hz
  ! 			Z(5) = Zzx; Z(6) = Zzy in addition
  !   iDT=Off_Diagonal_Impedance
  !  			Z(1) = Zxy, Z(2) = Zyx
  !   iDT=Full_Vertical_Components
  !  			Z(1) = Tx, Z(2) = Ty
  !   iDT=Full_Interstation_TF
  ! 		 	Z(1) = Mxx; Z(2) = Mxy; Z(3) = Myx; Z(4) = Myy
  !   iDT=Off_Diagonal_Rho_Phi (Added on behalf of Kristina Tietze, GFZ-Potsdam)
  !  			Z(1) = Rhoxy , Z(2) = Phixy, Z(3) = Rhoyx, Z(4) = Phiyx
  
  !  optional argument, useful for linearized impedance
  complex(kind=prec), intent(out), optional	:: Binv(2,2)

  !  local variables
  integer			:: iMode, i,j,xyz,ij
  real(kind=prec)	:: omega,x(3),x_ref(3)
  complex(kind=prec)    :: tempZ(2)
  complex(kind=prec)	:: BB(3,2),EE(2,2),RR(2,2)
  complex(kind=prec)	:: det,i_omega,ctemp
  type(sparsevecC)		:: Lex,Ley,Lbx,Lby,Lbz,Lrx,Lry
  logical			:: ComputeHz,ComputeE

  !  probably should dependence on omega into BinterpSetup, as in 2D!
  omega = txDict(ef%tx)%omega


  
  selectcase (iDT)
     case(Full_Impedance)
               x     = rxDict(iRX)%x         !Local site position (x,y,z)
		     ! First set up interpolation functionals for Ex, Ey
			  xyz = 1
			  call EinterpSetUp(ef%grid,x,xyz,Lex)
			  xyz = 2
			  call EinterpSetUp(ef%grid,x,xyz,Ley)
			 ! Then set up interpolation functionals for Bx, By      
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  ! loop over modes
			  do iMode = 1,2
			      ! electric fields
			      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
			      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			 end do
			 !invert horizontal B matrix using Kramer's rule.
			  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
			  ctemp = BB(1,1)
			  BB(1,1) =  BB(2,2)/det
			  BB(2,2) =  ctemp/det
			  BB(1,2) = -BB(1,2)/det
			  BB(2,1) = -BB(2,1)/det
			  	      
		        do j = 1,2
		           do i = 1,2
		              ij = 2*(i-1)+j
		              Z(ij) = EE(i,1)*BB(1,j)+EE(i,2)*BB(2,j)
		           enddo
		        enddo
		        
     case(Impedance_Plus_Hz)
              x     = rxDict(iRX)%x          !Local site position (x,y,z)
		     ! First set up interpolation functionals for Ex, Ey
			  xyz = 1
			  call EinterpSetUp(ef%grid,x,xyz,Lex)
			  xyz = 2
			  call EinterpSetUp(ef%grid,x,xyz,Ley)
			 ! Then set up interpolation functionals for Bx, By, Bz      
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  xyz = 3
     		  call BfromESetUp(ef%grid,omega,x,xyz,Lbz) 
			  ! loop over modes
			  do iMode = 1,2
			      ! electric fields
			      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
			      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			      BB(3,iMode) = dotProd_noConj_scvector_f(Lbz,ef%pol(iMode))
			 end do
			 !invert horizontal B matrix using Kramer's rule.
			  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
			  ctemp = BB(1,1)
			  BB(1,1) =  BB(2,2)/det
			  BB(2,2) =  ctemp/det
			  BB(1,2) = -BB(1,2)/det
			  BB(2,1) = -BB(2,1)/det
			    
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
              x     = rxDict(iRX)%x          !Local site position (x,y,z)
		     ! First set up interpolation functionals for Ex, Ey
			  xyz = 1
			  call EinterpSetUp(ef%grid,x,xyz,Lex)
			  xyz = 2
			  call EinterpSetUp(ef%grid,x,xyz,Ley)
			 ! Then set up interpolation functionals for Bx, By      
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  ! loop over modes
			  do iMode = 1,2
			      ! electric fields
			      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
			      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			 end do
			 !invert horizontal B matrix using Kramer's rule.
			  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
			  ctemp = BB(1,1)
			  BB(1,1) =  BB(2,2)/det
			  BB(2,2) =  ctemp/det
			  BB(1,2) = -BB(1,2)/det
			  BB(2,1) = -BB(2,1)/det
			  
			  
	        Z(1) = EE(1,1)*BB(1,2)+EE(1,2)*BB(2,2)
	        Z(2) = EE(2,1)*BB(1,1)+EE(2,2)*BB(2,1)
        
     case(Full_Vertical_Components)
               x     = rxDict(iRX)%x          !Local site position (x,y,z)
              !  Vertical field TF
			 ! First set up interpolation functionals for Bx, By, Bz      
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  xyz = 3
     		  call BfromESetUp(ef%grid,omega,x,xyz,Lbz) 
			  ! loop over modes
			  do iMode = 1,2
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			      BB(3,iMode) = dotProd_noConj_scvector_f(Lbz,ef%pol(iMode))
			 end do
			 !invert horizontal B matrix using Kramer's rule.
			  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
			  ctemp = BB(1,1)
			  BB(1,1) =  BB(2,2)/det
			  BB(2,2) =  ctemp/det
			  BB(1,2) = -BB(1,2)/det
			  BB(2,1) = -BB(2,1)/det
			         		  
       
              Z(1) = BB(3,1)*BB(1,1)+BB(3,2)*BB(2,1)
              Z(2) = BB(3,1)*BB(1,2)+BB(3,2)*BB(2,2)
     case(Full_Interstation_TF)
              x     = rxDict(iRX)%x          !Local site position (x,y,z)
              x_ref = rxDict(iRX)%r          !Reference site position (x,y,z)
  			 ! First set up interpolation functionals for Bx, By at local site    
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)   
		     !Then set up interpolation functionals for Bx, By at the referance site 
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x_ref,xyz,Lrx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x_ref,xyz,Lry)
			    do iMode = 1,2
			      ! magnetic fields at local station
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			      ! magnetic fields, at the REFERANCE station
			      RR(1,iMode) = dotProd_noConj_scvector_f(Lrx,ef%pol(iMode))
			      RR(2,iMode) = dotProd_noConj_scvector_f(Lry,ef%pol(iMode))
			    end do
			  ! Compute the inverse of RR using Kramer's rule
			  det = RR(1,1)*RR(2,2)-RR(1,2)*RR(2,1)
			  ctemp = RR(1,1)
			  RR(1,1) =  RR(2,2)/det
			  RR(2,2) =  ctemp/det
			  RR(1,2) = -RR(1,2)/det
			  RR(2,1) = -RR(2,1)/det
			  ! Z = BB * RR^-1
			         do j = 1,2
			           do i = 1,2
			              ij = 2*(i-1)+j
			              Z(ij) = BB(i,1)*RR(1,j)+BB(i,2)*RR(2,j)
			           enddo
			        enddo
       case(Off_Diagonal_Rho_Phi)
               x     = rxDict(iRX)%x          !Local site position (x,y,z)
		     ! First set up interpolation functionals for Ex, Ey
			  xyz = 1
			  call EinterpSetUp(ef%grid,x,xyz,Lex)
			  xyz = 2
			  call EinterpSetUp(ef%grid,x,xyz,Ley)
			 ! Then set up interpolation functionals for Bx, By      
			  xyz = 1
			  call BfromESetUp(ef%grid,omega,x,xyz,Lbx)
			  xyz = 2
			  call BfromESetUp(ef%grid,omega,x,xyz,Lby)
			  ! loop over modes
			  do iMode = 1,2
			      ! electric fields
			      EE(1,iMode) =  dotProd_noConj_scvector_f(Lex,ef%pol(iMode))
			      EE(2,iMode) =  dotProd_noConj_scvector_f(Ley,ef%pol(iMode))
			      ! magnetic fields
			      BB(1,iMode) = dotProd_noConj_scvector_f(Lbx,ef%pol(iMode))
			      BB(2,iMode) = dotProd_noConj_scvector_f(Lby,ef%pol(iMode))
			 end do
			 !invert horizontal B matrix using Kramer's rule.
			  det = BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1)
			  ctemp = BB(1,1)
			  BB(1,1) =  BB(2,2)/det
			  BB(2,2) =  ctemp/det
			  BB(1,2) = -BB(1,2)/det
			  BB(2,1) = -BB(2,1)/det
			        
		       tempZ(1) = EE(1,1)*BB(1,2)+EE(1,2)*BB(2,2)
		       tempZ(2) = EE(2,1)*BB(1,1)+EE(2,2)*BB(2,1)
		       Z(1)   = abs(tempZ(1))**2*MU_0/omega
		       z(2)   = atan2(ISIGN*imag(tempZ(1)),real(tempZ(1)))*R2D
		       Z(3)   = abs(tempZ(2))**2*MU_0/omega
		       Z(4)   = atan2(ISIGN*imag(tempZ(2)),real(tempZ(2)))*R2D                
    end select

  if(present(Binv)) then
      if(typeDict(iDT)%tfType .eq. Full_Interstation_TF) then
         Binv = RR(1:2,:)
      else
         Binv = BB(1:2,:)
      end if   
  endif

  ! clean up
  call deall_sparsevecc(Lex)
  call deall_sparsevecc(Ley)
  call deall_sparsevecc(Lbx)
  call deall_sparsevecc(Lby)
  call deall_sparsevecc(Lbz)
  call deall_sparsevecc(Lrx)
  call deall_sparsevecc(Lry)
  
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

  type (solnVector_t), intent(in)		   :: e0
  type (modelParam_t), intent(in)  :: Sigma0
  integer, intent(in)			   :: iDT, iRX
  !   NOTE: Lz and Qz have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE:  in principal the comparable input arguments in
  !        the 2D program should also be of type sparseVector!
  type(sparseVector_t), intent(inout)		:: L(*),Q(*)

  !  local variables
  complex(kind=prec)	:: Binv(2,2)
  complex (kind=prec)	:: Z(6), i_omega,c1,c2
  real(kind=prec)	:: x(3),x_ref(3),omega
  type(sparsevecc)		:: L1
  integer			:: i,j,k,nComp,IJ(3,6),xyz,n, predictedComp
  type(sparsevecC)		:: Lex,Ley,Lbx,Lby,Lbz,Lrx,Lry
  logical			:: ComputeHz

  omega = txDict(e0%tx)%omega
  	 x     = rxDict(iRX)%x
     x_ref = rxDict(iRX)%r          !Reference site position (x,y,z)

  !  set up which components are needed,  ... and ! evaluate
  !   impedance, Binv for background solution
  !          ... appear to need full impedance for offdiagonal

  !   Some modifications to allow for other TFs ... e.g., the case
  !     of Hz TFs only (changes also slightly simplify generalization
  !      to interstation TFs) :  increase first dimension of array IJ
  !      to 3: (IJ(1,:) = row index in TF matrix Z;
  !             IJ(2,:) = column index in TF martrix X;
  !             IJ(3,:) = predicted field component ...
  !                     Ex = 1; Ey =2; Bz = 3; (Bx = 4; By = 5,  at referance site)  
  !						(can add more cases)
  !
  select case(iDT)
     case(Full_Impedance)
        nComp = 4
        ComputeHz = .false.
        do j = 1,2
           do i = 1,2
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
              IJ(3,2*(i-1)+j) = i
           enddo
        enddo
        Call nonlinDataFunc(e0,Sigma0,Full_Impedance,iRX,Z,Binv)
     case(Impedance_Plus_Hz)
        nComp = 6
        ComputeHz = .true.
        do j = 1,2
           do i = 1,3
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
              IJ(3,2*(i-1)+j) = i
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
        IJ(3,1) = 1
        IJ(3,2) = 2
        Call nonlinDataFunc(e0,Sigma0,Full_Impedance,iRX,Z,Binv)
     case(Full_Vertical_Components)
        nComp = 2
        ComputeHz = .true.
        IJ(1,1) = 1
        IJ(1,2) = 1
        IJ(2,1) = 1
        IJ(2,2) = 2
        IJ(3,1) = 3
        IJ(3,2) = 3
        Call nonlinDataFunc(e0,Sigma0,Full_Vertical_Components,iRX,Z,Binv)
     case(Full_Interstation_TF)
        nComp = 4
        ComputeHz = .false.
        do j = 1,2
           do i = 1,2
              IJ(1,2*(i-1)+j) = i
              IJ(2,2*(i-1)+j) = j
              IJ(3,2*(i-1)+j) = i+3
           enddo
        enddo
        Call nonlinDataFunc(e0,Sigma0,Full_Interstation_TF,iRX,Z,Binv) 
     endselect

  ! Then set up interpolation functionals for Ex, Ey, Bx, By, Bz, (Bx,By at referance station)
  xyz = 1
  call EinterpSetUp(e0%grid,x,xyz,Lex)
  xyz = 2
  call EinterpSetUp(e0%grid,x,xyz,Ley)
  xyz = 1
  call BfromESetUp(e0%grid,omega,x,xyz,Lbx)
  xyz = 2
  call BfromESetUp(e0%grid,omega,x,xyz,Lby)
  xyz = 1
  call BfromESetUp(e0%grid,omega,x_ref,xyz,Lrx)
  xyz = 2
  call BfromESetUp(e0%grid,omega,x_ref,xyz,Lry)
  if(ComputeHz) then
     xyz = 3
     call BfromESetUp(e0%grid,omega,x,xyz,Lbz)
  endif

  !  compute sparse vector representations of linearized functionals
  do n = 1,nComp
     !  i runs over rows of TF matrix, j runs over columns of TF
     i = IJ(1,n)
     j = IJ(2,n)
     predictedComp = IJ(3,n)
     c1 = Z(2*(i-1)+1)
     c2 = Z(2*(i-1)+2)
    if(typeDict(iDT)%tfType .eq. Full_Interstation_TF) then 
      Call linComb_sparsevecc(Lrx,c1,Lry,c2,L1)
    else 
      Call linComb_sparsevecc(Lbx,c1,Lby,c2,L1)
    end if  
     do k = 1,2
        !  k defines which mode the linearized functional is
        !   to be applied to
        c1 = Binv(k,j)  !In case of interstaion TF, Binv = RRinv.
        c2 = -c1
        if(predictedComp.eq.1) then
           !  component in x row of impedance tensor
           Call linComb_sparsevecc(Lex,c1,L1,c2,L(n)%L(k))
        elseif(predictedComp.eq.2) then
           !  component in y row of impedance tensor
           Call linComb_sparsevecc(Ley,c1,L1,c2,L(n)%L(k))
        elseif(predictedComp.eq.3) then
           !  component in Bz row (vertical field TF)
           Call linComb_sparsevecc(Lbz,c1,L1,c2,L(n)%L(k))
        elseif(predictedComp.eq.4) then
           !  component in x row (interstation TF)
           Call linComb_sparsevecc(Lrx,c1,L1,c2,L(n)%L(k))
        elseif(predictedComp.eq.5) then
           !  component in y row (interstation TF)
           Call linComb_sparsevecc(Lry,c1,L1,c2,L(n)%L(k))        
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
  call deall_sparsevecc(Lrx)
  call deall_sparsevecc(Lry)
  
  end subroutine linDataFunc

end module dataFunc
