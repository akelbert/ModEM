module MeasComb
!  higher level module that calls routines in dataFunc2D to
!   --> evaluate non-linear and linear data functionals for all sites
!       in a dvec object
!   --> construct the comb for forcing the adjoint model for the linear
!       combination of data functionals corresponding to a dvec object
!   This module is generic, and should work for a broad range of different
!     problems (in particular, 2D and 3D MT).  Much of the code in these
!     routines is for handling two general efficiency issues that arise
!     with frequency domain EM data: Even at a single site there can be multiple 
!     components that can be evaluated from a single EMsoln object: e.g.,
!     in 2D MT, TE mode impedance and Tipper; or for 3D 4 components of
!     impedance tensor; and data can be complex or real.  These modules
!     handle these different cases efficiently, so that these details 
!     are hidden from manipulations in the sensitivity calculation module
!   This module also handles the sparse vectors Q, which define the
!     direct sensitivity of the data functionals to model parameters.
 
!*****************************************************************************
use math_constants
use utilities
use dataspace
use modelparameter
use solnrhs
use dataFunc

implicit none

 public 		:: linDataMeas,  linDataComb, dataMeas

Contains

  subroutine linDataMeas(e0,ef,d,dSigma)
  ! given measured and  background electric field solutions (ef,e0)
  ! evaluate linearized data functionals for all sites represented
  !  in a dvec object.

  !  electric field solutions are stored as type EMsoln
  type (EMsoln), intent(in)			:: ef,e0
  ! d provides indices into receiver dictionary on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dvec), intent(inout)            	:: d
  !  optional input of conductivity perturbation dSigma used to
  !   solve for ef, needed if data functional depends on input parameter
  !   and this routine is used to multiply J (sensitivity matrix) times
  !    a model parameter vector
  type (modelParam_t), intent(in), optional	:: dSigma	

  !  local variables
  complex(kind=selectedPrec)	:: Z
  type(EMsparse), pointer	:: Lz(:), Qz(:)
  logical               	:: Conj_Case = .false.
  logical               	:: calcQ
  integer               	:: iSite, ncomp, itx, &
					nFunc, iDT, iComp, iFunc
  character*80          	:: msg
  type (modelParam_t)		:: sigmaQreal, sigmaQimag

  !  first check consistency of e0, ef, d
  !             ... all should point to same transmitter!
  if((d%tx.ne.e0%tx) .or. (d%tx.ne.ef%tx) .or. &
          (e0%tx.ne.ef%tx)) then
     msg = 'Error: transmitter incompatability in linDataMeas'
     call errStop(msg)
  endif

  itx = d%tx
  iDT = d%dataType
  calcQ = typeDict(iDT)%calcQ
  !  calcQ is true when data functional coefficients depend on
  !   model parameters 
  if(calcQ .and. .not.present(dSigma)) then
     msg = 'Error: dSigma required as input to linDataMeas for this data type'
     call errStop(msg)
  endif
  ncomp = d%ncomp
  if(typeDict(iDT)%isComplex) then
     !  data are complex; one sensitivity calculation can be
     !   used for both real and imaginary parts
     if(mod(ncomp,2).ne.0) then
        msg = 'Error: for complex data # of components must be even'
        call errStop(msg)
     endif
     nFunc = ncomp/2
  else
     !  data are treated as real: full sensitivity computation is required
     !   for each component
     nFunc = ncomp
  endif
  allocate(Lz(nFunc))
  allocate(Qz(nFunc))

  !  loop over sites
  do iSite = 1,d%nSite
     ! compute sparse vector representations of linearized
     ! data functionals for transfer function (TE impedance)
     ! elements at one site
     !   Note that elements of QZ are created inside linDataFunc only
     !    if required
     !   linDataFunc returns one Lz (and Qz if appropriate)
     !    for each of nFunc functionals
     call linDataFunc(e0,iDT,d%rx(iSite),Lz,Qz)
     if((iSite .eq. 1).and.calcQ) then
        call copy_ModelParam(sigmaQreal,e0%sigma)
        if(typeDict(iDT)%isComplex) then
           call copy_ModelParam(sigmaQimag,e0%sigma)
        endif
     endif 
     iComp = 1
     do iFunc  = 1, nFunc
        Z = dotProd_EMsparseEMsoln(Lz(iFunc),ef,Conj_Case)
        if(typeDict(iDT)%isComplex) then
           d%data(iComp,iSite) = real(Z)
           iComp = iComp + 1
           d%data(iComp,iSite) = imag(Z)
           iComp = iComp + 1
           if(calcQ) then
              call zero_ModelParam(sigmaQreal)
              call zero_ModelParam(sigmaQimag)
              call EMSparseQtoModelParam(C_ONE,Qz(iFunc), &
			e0%sigma,sigmaQreal,sigmaQimag)
              d%data(iComp-2,iSite) = d%data(iComp-2,iSite) &
	      		+ dotProd_modelParam(sigmaQreal,dSigma)	
              d%data(iComp-1,iSite) = d%data(iComp-1,iSite) &
       			+ dotProd_modelParam(sigmaQimag,dSigma)				
           endif
        else
           d%data(iComp,iSite) = real(Z)
           iComp = iComp + 1
           if(calcQ) then
              call zero_ModelParam(sigmaQreal)
              call EMSparseQtoModelParam(C_ONE,Qz(iFunc),e0%sigma,sigmaQreal)
              d%data(iComp-1,iSite) = d%data(iComp-1,iSite) +  &
			dotProd_modelParam(sigmaQreal,dSigma)
           endif
        endif
     enddo
  enddo
  !  deallocate local arrays
  do iFunc = 1,nFunc
     call deall_EMSparse(Lz(iFunc))
     if(calcQ) then
        call deall_EMSparse(Qz(iFunc))
        call deall_modelParam(sigmaQreal)
        if(typeDict(iDT)%isComplex) then
           call deall_modelParam(sigmaQimag)
        endif
     endif 
  enddo
  deallocate(Lz)
  deallocate(Qz)

  end subroutine linDataMeas

!*****************************************************************************
  subroutine linDataComb(e0,d,comb,Qcomb)

  ! given background electric field solution
  ! and a dvec object (element of data space containing
  ! MT data for one frequency, one or more sites) compute adjoint
  ! of measurement operator: i.e., the comb constructed from the scaled
  ! superposition of data kernals (scaled by conjugate of data values ...
  !  e.g., if d contains residuals, this can be used to set up for
  !  gradient calculation).
  ! Also returns conductivity parameter data structure corresponding
  !   to part of sensitivity due to derivative of data functionals
  !   with respect to model parameters
  !
  !  NOTE: this will not zero comb (or Qcomb): repeated calls with different
  !  instances of dvec will add new comb elements to input comb
  !  NOTE: we are only supporting full storage sources in comb;
  !    the elements of comb should be allocated and zeroed before calling
  ! As with linDataPred, all of the receiver, transmitter, and data type
  !   information is obtained from the dictionaries using indices 
  !    stored in d%rx, d%tx, d%dataType

  !  electric field solutions are stored as type EMsoln
  type (EMsoln), intent(in)     :: e0
  ! d provides indices into receiver dictionary on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dvec), intent(in)               	:: d
  ! Output 
  type (EMrhs), intent(inout)          		:: comb
  ! Optional output
  type (modelParam_t), intent(inout),optional	:: Qcomb

  !  local variables
  integer               		:: iSite, iTx, iDT, nComp, &
					 nFunc, iComp, iFunc
  complex(kind=selectedPrec)		:: Z
  logical               		:: calcQ
  type(EMSparse),pointer    	:: Lz(:),Qz(:)
  character*80				:: msg

  itx = d%tx
  iDT = d%dataType
  calcQ = typeDict(iDT)%calcQ
  !  calcQ is true when data functional coefficients depend on
  !   model parameters   
  nComp = d%nComp
  if(typeDict(iDT)%isComplex) then
     !  data are complex; one sensitivity calculation can be
     !   used for both real and imaginary parts
     if(mod(ncomp,2).ne.0) then
        msg = 'Error: for complex data # of components must be even'
        call errStop(msg)
     endif
     nFunc = ncomp/2
  else
     !  data are treated as real: full sensitivity computation is required
     !   for each component
     nFunc = ncomp
  endif
  allocate(Lz(nFunc))
  allocate(Qz(nFunc))
  
  !  loop over sites
  do iSite = 1,d%nSite
     ! compute sparse vector representations of linearized
     ! data functionals for transfer function
     ! elements at one site
     call linDataFunc(e0,iDT,d%rx(iSite),Lz,Qz)
     iComp = 1
     do iFunc  = 1, nFunc
        if(typeDict(iDT)%isComplex) then
           !  move real data in dvec into complex conjugate of TF (impedance)
           !  multiply this by data kernel for complex impedance ...
           !      (take real part in parameter space)
           Z = cmplx(d%data(iComp,iSite),-d%data(iComp+1,iSite),8)
           iComp = iComp+2
        else
           Z = cmplx(d%data(iComp,iSite),0.0,8)
           iComp = iComp+1
        endif
        call add_EMsparseEMrhs(Z,Lz(iFunc),comb)
        if(calcQ) then
           ! adds to input Qcomb (type modelParam)
           call EMSparseQtoModelParam(Z,Qz(iFunc),e0%sigma,Qcomb)
        endif
     enddo
   enddo

  !  deallocate local arrays
  do iFunc = 1,nFunc
     call deall_EMSparse(Lz(iFunc))
     if(calcQ) then
        call deall_EMSparse(Qz(iFunc))
     endif 
  enddo
  deallocate(Lz)
  deallocate(Qz)

  end subroutine linDataComb

!****************************************************************************
  subroutine dataMeas(ef,d)
  ! given solution for a single TX ef compute predicted data 
  !  at all sites, returning result in d
  !   Data type information is obtained from typeDict, 
  !   using dictionary indices stored in d%dataType
  !   Calls nonLinDataFunc to do the actual impedance calculation
  
  type (EMsoln), intent(in)     :: ef
  ! d provides indices into receiver dictionary on
  ! input. Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dvec), intent(inout)    :: d

  !  local variables
  integer               ::  iSite, ncomp,nFunc,iDT,iComp, iFunc
  complex (kind=selectedPrec), pointer, dimension(:)      ::  Z
  character(len=80)		:: msg

  
  iDT = d%dataType
  ncomp = d%ncomp
  if(typeDict(iDT)%isComplex) then
     !  data are complex; one sensitivity calculation can be
     !   used for both real and imaginary parts
     if(mod(ncomp,2).ne.0) then
        msg = 'Error: for complex data # of components must be even'
        call errStop(msg)
     endif
     nFunc = ncomp/2
  else
     !  data are treated as real
     nFunc = ncomp
  endif
  allocate(Z(nFunc))

  ! loop over sites
  do iSite = 1,d%nSite
     call nonLinDataFunc(ef,iDT,d%rx(iSite),Z)
     !  copy responses in Z (possibly complex, possibly multiple
     !         components) into dvec object d
     !  Loop over components
     iComp = 0
     do iFunc  = 1, nFunc
        if(typeDict(iDT)%isComplex) then
           iComp = iComp + 1
           d%data(iComp,iSite) = real(Z(iFunc))
           iComp = iComp + 1
           d%data(iComp,iSite) = imag(Z(iFunc))
        else
           iComp = iComp + 1
           d%data(iComp,iSite) = real(Z(iFunc))
        endif  
     enddo  
  enddo
  
  ! predicted data have no error bars defined
  d%errorBar = .false.

  end subroutine dataMeas

!****************************************************************************

end module MeasComb
