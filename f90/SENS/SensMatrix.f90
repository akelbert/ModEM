module sensMatrix

use math_constants
use utilities
use meascomb	 !!!!  inherits : dataspace, dataFunc, solnrhs
use modelparamsens  !!!  inherits : modelspace, soln2d
use emsolver

implicit none

public 	:: calcSensMatrix, Jmult, JmultT, fwdPred, JmultT_MTX, &
		Jmult_MTX

Contains

   !**********************************************************************
   subroutine calcSensMatrix(d,sigma0,dsigma)
   !  Calculate sensitivity matrix for data in d
   !  Approaching a generic form that will work for any problem;
   !   Documentation not edited.  Code not debugged!
   !
   !  Mostly self-contained, but BEFORE CALLING:
   !   1)  read in and initialize grid (sigma0 points
   !        to this grid)
   !   2)  call setWSparams to set grid dimensions inside
   !         WS forward modeling module
   !   3)  set up transmitter and receiver dictionaries;
   !       "pointers" to entries in these dictionaries are
   !        attached to multi-transmitter data vector d .
   !        Note that the actual frequencies used for the
   !        solver are extracted from the transmitter dictionary
   !   IT IS NOT necessary to initialize for modeling
   !     before calling
   !
   !   d is the input data vector, here just used to identify
   !     receiver transmitter pairs to compute sensitivities for
   type(dvecMTX), intent(inout)	:: d
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the output array of data sensitivities,
   !   one for each element in the data array.  Each sensitivity
   !    is an element of type modelParam, an abstract
   !    data type that defines the unknow conductivity
   type(modelParam_t), pointer   :: dsigma(:)

   !  local variables
   type(EMsoln)		:: e,e0
   type(EMrhs)		:: comb
   integer 		:: j,nTx,k,nSite,nTot,ii,iTx, &
				iDT,nfunc,ncomp,iRx,iFunc
   type(EMsparse), pointer	:: L(:),Q(:)
   logical 		:: calcQ

   ! nTX is number of transmitters; nTot is number of real data
   nTx = d%nTx
   nTot = d%Ndata

   if(.not.associated(dsigma)) then
      ! allocate for sensitivity matrix
      allocate(dsigma(nTot))
      do j = 1,nTot
         ! this makes a copy of model param, of the same type
         !   as sigma0, then zeros it.
         call copy_ModelParam(dsigma(j),sigma0)
         call zero_ModelParam(dsigma(j))
      enddo
   endif

   ! loop over frequencies, computing all sensitivities for
   !   one frequency
   ii = 0
   do j = 1,nTx
      !   get indices into dictionaries for this dvec
      iTx = d%d(j)%tx
      iDT = d%d(j)%dataType
      calcQ = typeDict(iDT)%calcQ

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,e0,e,comb)

      !  solve forward problem; result is stored in e0
      call fwdSolve(iTx,e0)

      ! get data type info for this dvec: here we allow for either
      !  real or complex data, though we always store as real
      !  Complex data come in pairs, stored as two succesive
      !  "components" in the dvec structure.
      !  "isComplex" is an attribute of a dvec structure; can
      !  be different for different dvecs within the overall data
      !  vector
      ncomp = d%d(j)%ncomp
      nSite = d%d(j)%nSite
      if(typeDict(iDT)%isComplex) then
         !  data are complex; one sensitivity calculation can be
         !   used for both real and imaginary parts
         if(mod(ncomp,2).ne.0) then
            call errStop('for complex data # of components must be even in CalcSensMatrix')
         endif
         nFunc = ncomp/2
      else
         !  data are treated as real: full sensitivity computation is required
         !   for each component
         nFunc = ncomp
      endif
      allocate(L(nFunc))
      allocate(Q(nFunc))

      ! loop over sites, computing sensitivity for all components
      !   for each site
      do k = 1,nSite

         iRx = d%d(j)%rx(k)

         !  compute linearized data functional(s) : L
         call linDataFunc(e0,sigma0,iDT,iRx,L,Q)

         ! loop over functionals  (for TE/TM impedances nFunc = 1)
         do iFunc = 1,nFunc

            ! solve forward problem for each of nFunc functionals
            call zero_EMrhs(comb)
            call add_EMsparseEMrhs(C_ONE,L(iFunc),comb)

            call sensSolve(iTx,TRN,comb,e)

            ! multiply by P^T (and add Q^T if appropriate)
            if(typeDict(iDT)%isComplex) then
               ii = ii + 2
               call PmultT(e0,sigma0,e,dsigma(ii-1),dsigma(ii))
               if(calcQ) then
                  call EMSparseQtoModelParam(C_ONE,Q(iFunc),sigma0, &
			dsigma(ii-1),dsigma(ii))
               endif
            else
               ii = ii+1
               call PmultT(e0,sigma0,e,dsigma(ii))
               if(calcQ) then
                  call EMSparseQtoModelParam(C_ONE,Q(iFunc),sigma0,dsigma(ii))
               endif
            endif
         enddo  ! iFunc
      enddo
      deallocate(L)
      deallocate(Q)
   enddo

   call exitSolver(e0,e,comb)

   end subroutine calcSensMatrix

   !**********************************************************************
   subroutine Jmult(delSig,sigma0,d,eAll)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for all transmitters in a datavector (i.e., multiple dvec objects)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  indicies into dictionaries are attached to data vector d .
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions)
   !
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   delSig is the input conductivity parameter perturbation
   type(modelParam_t), intent(in)	:: delSig
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dvecMTX), intent(inout)		:: d
   type(EMsolnMTX), intent(inout), optional	:: eAll

   !  local variables
   type(EMsoln)		:: e,e0
   type(EMrhs)		:: comb
   integer 		:: j,iTx,iDT
   logical		:: savedSolns

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in Jmult')
      endif
   endif

   ! loop over frequencies : solve forward system twice,
   !    to compute background, perturbation solutions
   !    apply data functionals (linearized about background soln)
   !    to perturbation solution
   do j = 1,d%nTx

      !   get indices into dictionaries for this dvec
      iTx = d%d(j)%tx
      iDT = d%d(j)%dataType

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,e0,e,comb)

      if(savedSolns) then
         !e0 = eAll%solns(j)
         call copy_EMsoln(e0,eAll%solns(j))
      else
         !  solve forward problem; result is stored in e0
         call fwdSolve(iTx,e0)
      endif

      !  compute rhs (stored in comb) for forward sensitivity
      !  calculation, using conductivity perturbations and
      !  background soln:
      call Pmult(e0,sigma0,delSig,comb)

      ! solve forward problem with source in comb
      call sensSolve(iTx,FWD,comb,e)

      ! finally apply linearized data functionals
      if(TypeDict(iDT)%calcQ) then
         call linDataMeas(e0,sigma0,e,d%d(j),delSig)
      else
	 call linDataMeas(e0,sigma0,e,d%d(j))
      endif

   enddo

   !  clean up
   call exitSolver(e0,e,comb)

   end subroutine Jmult

   !**********************************************************************
   subroutine JmultT(sigma0,d,dsigma,eAll)

   !  Transpose of Jmult mujltiplied by data vector d; output is a
   !      single conductivity parameter in dsigma
   !
   !   sigma0 is background conductivity parameter
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   IN THE PROPER ORDER (at present) !!!!
   !   (This can be used to avoid recomputing forward solutions,
   !    e.g., in a CG solution scheme)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  "pointers" to dictionary entries are attached to multi-transmitter
   !   data vector d

   type(modelParam_t), intent(in)	:: sigma0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dvecMTX), intent(in)		:: d
   !   delSig is the output conductivity parameter
   type(modelParam_t), intent(Out)  	:: dsigma
   type(EMsolnMTX), intent(in), optional	:: eAll

   !  local variables
   type(EMsoln)		:: e,e0
   type(EMrhs) 		:: comb
   type(modelParam_t)	:: sigmaTemp, Qcomb
   integer 		:: j,iTx,iDT
   logical		:: calcSomeQ, firstQ
   logical		:: savedSolns

   calcSomeQ = .false.
   firstQ = .true.
   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in JmultT')
      endif
   endif

   call copy_ModelParam(dsigma,sigma0)
   call zero_ModelParam(dsigma)

   ! loop over transmitters:
   do j = 1,d%nTx

      !   get indices into dictionaries for this dvec
      iTx = d%d(j)%tx
      iDT = d%d(j)%dataType

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,e0,e,comb)
      if(j.eq.1) then
        call copy_ModelParam(sigmaTemp,sigma0)
	call zero_ModelParam(sigmaTemp)
      endif

      if(savedSolns) then
         !e0 = eAll%solns(j)
	     call copy_EMsoln(e0,eAll%solns(j))
      else
         !  solve forward problem; result is stored in e0
         call fwdSolve(iTx,e0)
      endif

      ! set up comb using linearized data functionals
      !  ... for dvecs with data functionals depending on conductivity
      !   parameter also compute analagous comb in parameter space
      if(typeDict(iDT)%calcQ) then
         if(firstQ) then
            !  first transmitter for which Q must be calculated:
            !   ==> allocate and zero Qcomb (use copy so that paramtype
            !         is set correctly)
            call copy_ModelParam(Qcomb,sigmaTemp)
            !  NOTE: linDataComb ADDS to Qcomb, not overwrites
            !   ==> only zero Qcomb for first transmitter requiring Q
            call zero_ModelParam(Qcomb)
            !  set flags indicating that Q is now non-zero
            calcSomeQ = .true.
            firstQ = .false.
          endif
            !  BUT: linDataComb overwrites comb ... so zero this
            !       for every transmitter
          call zero_EMrhs(comb)
          call linDataComb(e0,sigma0,d%d(j),comb,Qcomb)
      else
          call zero_EMrhs(comb)
          call linDataComb(e0,sigma0,d%d(j),comb)
      endif
      ! solve forward problem with source in comb
      call sensSolve(iTx,TRN,comb,e)

      ! map from nodes to conductivity space ... result in sigmaTemp
      !  HERE PmultT overwrites sigmaTemp
      !  NOTE: here we throw away imaginary part, even for complex
      !     data (conceivably might want to save this in some cases!)
      call PmultT(e0,sigma0,e,sigmaTemp)
      call linComb_modelParam(ONE,dsigma,ONE,sigmaTemp,dsigma)
   enddo

   if(calcSomeQ) then
      !  add Qcomb
      call linComb_ModelParam(ONE,dsigma,ONE,Qcomb,dsigma)
      call deall_modelParam(Qcomb)
   endif

   !  clean up
   call exitSolver(e0,e,comb)
   call deall_modelParam(sigmaTemp)

   end subroutine JmultT

   !**********************************************************************
   subroutine Jmult_MTX(delSig,sigma0,d,eAll)

   !  Calculate product of sensitivity matrix and a separate
   !    model parameter for each transmitter in a datavector
   !    (i.e., multiple dvec objects)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  indicies into dictionaries are attached to data vector d .
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions)
   !
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   delSig is the array of input conductivity parameter perturbations
   type(modelParam_t), dimension(:), pointer	:: delSig
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dvecMTX), intent(inout)		:: d
   type(EMsolnMTX), intent(inout), optional	:: eAll

   !  local variables
   type(EMsoln)		:: e,e0
   type(EMrhs)		:: comb
   integer 		:: j,iTx,iDT,nTx
   logical		:: savedSolns

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in Jmult_MTX')
      endif
   endif

   nTx = d%nTx
   if(.not.associated(delSig)) then
      ! allocate for sensitivity matrix
      allocate(delSig(nTx))
      do j = 1,nTx
         ! this makes a copy of model param, of the same type
         !   as sigma0, then zeros it.
         call copy_ModelParam(delSig(j),sigma0)
         call zero_ModelParam(delSig(j))
      enddo
   endif

   ! loop over frequencies : solve forward system twice,
   !    to compute background, perturbation solutions
   !    apply data functionals (linearized about background soln)
   !    to perturbation solution
   do j = 1,nTx

      !   get indices into dictionaries for this dvec
      iTx = d%d(j)%tx
      iDT = d%d(j)%dataType

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,e0,e,comb)

      if(savedSolns) then
         ! e0 = eAll%solns(j)
         call copy_EMsoln(e0,eAll%solns(j))
      else
         !  solve forward problem; result is stored in e0
         call fwdSolve(iTx,e0)
      endif

      !  compute rhs (stored in comb) for forward sensitivity
      !  calculation, using conductivity perturbations and
      !  background soln:
      call Pmult(e0,sigma0,delSig(j),comb)

      ! solve forward problem with source in comb
      call sensSolve(iTx,FWD,comb,e)

      ! finally apply linearized data functionals
      if(TypeDict(iDT)%calcQ) then
         call linDataMeas(e0,sigma0,e,d%d(j),delSig(j))
      else
	 call linDataMeas(e0,sigma0,e,d%d(j))
      endif

   enddo

   !  clean up
   call exitSolver(e0,e,comb)

   end subroutine Jmult_MTX

   !**********************************************************************
   subroutine JmultT_MTX(sigma0,d,dsigma,eAll)

   !  Transpose of Jmult mujltiplied by data vector d; output is an
   !       array of conductivity parameters in allocatble array dsigma,
   !       one for each transmitter.  The sum of dsigma over transmitters
   !       is the same as the output of JmultT
   !
   !   sigma0 is background conductivity parameter
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   IN THE PROPER ORDER (at present) !!!!
   !   (This can be used to avoid recomputing forward solutions,
   !    e.g., in a CG solution scheme)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  "pointers" to dictionary entries are attached to multi-transmitter
   !   data vector d

   type(modelParam_t), intent(in)	:: sigma0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dvecMTX), intent(in)		:: d
   !   delSig is the output conductivity parameter
   type(modelParam_t), dimension(:), pointer  	:: dsigma
   type(EMsolnMTX), intent(in), optional	:: eAll

   !  local variables
   type(EMsoln)		:: e,e0
   type(EMrhs) 		:: comb
   type(modelParam_t)	:: Qcomb
   logical		:: calcSomeQ = .false.
   integer 		:: j,iTx,iDT,nTx
   logical		:: savedSolns

   ! nTX is number of transmitters;
   nTx = d%nTx

   if(.not.associated(dsigma)) then
      ! allocate for sensitivity matrix
      allocate(dsigma(nTx))
      do j = 1,nTx
         ! this makes a copy of model param, of the same type
         !   as sigma0, then zeros it.
         call copy_ModelParam(dsigma(j),sigma0)
         call zero_ModelParam(dsigma(j))
         iDT = d%d(j)%dataType
         if(typeDict(iDT)%calcQ) then
            calcSomeQ = .true.
         endif
      enddo
   endif

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in JmultT_MTX')
      endif
   endif

   if(calcSomeQ) then
      !  will need model parameter storage in Qcomb; use copy to create
      call copy_ModelParam(Qcomb,sigma0)
   endif

   ! loop over transmitters:
   do j = 1,d%nTx

      !   get indices into dictionaries for this dvec
      iTx = d%d(j)%tx
      iDT = d%d(j)%dataType

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,e0,e,comb)

      if(savedSolns) then
         ! e0 = eAll%solns(j)
	     call copy_EMsoln(e0,eAll%solns(j))
      else
         !  solve forward problem; result is stored in e0
         call fwdSolve(iTx,e0)
      endif

      ! set up comb using linearized data functionals
      !  ... for dvecs with data functionals depending on conductivity
      !   parameter also compute analagous comb in parameter space
      call zero_EMrhs(comb)
      if(typeDict(iDT)%calcQ) then
         call linDataComb(e0,sigma0,d%d(j),comb,Qcomb)
      else
         call linDataComb(e0,sigma0,d%d(j),comb)
      endif

      ! solve forward problem with source in comb
      call sensSolve(iTx,TRN,comb,e)

      ! map from nodes to conductivity space ... result in sigmaTemp
      !  HERE PmultT overwrites sigmaTemp
      !  NOTE: here we throw away imaginary part, even for complex
      !     data (conceivably might want to save this in some cases!)
      call PmultT(e0,sigma0,e,dsigma(j))
      if(typeDict(iDT)%calcQ) then
         !  add Qcomb
         call linComb_ModelParam(ONE,dsigma(j),ONE,Qcomb,dsigma(j))
      endif
   enddo

   !  clean up
   call exitSolver(e0,e,comb)
   if(calcSomeQ) then
      call deall_modelParam(Qcomb)
   endif

   end subroutine JmultT_MTX

   !**********************************************************************
   subroutine fwdPred(sigma,d,eAll)

   !  Calculate predicted data for dvecMTX object d
   !    and for conductivity parameter sigma
   !  Optionally returns array of EM solutions eAll, one for
   !    each transmitter (if eAll is present)
   !
   !  First need to set up transmitter, receiver, dataType dictionaries;
   !  "pointers" to dictionaries are attached to multi-transmitter
   !   data vector d .
   !
   !   sigma is input conductivity parameter
   type(modelParam_t), intent(in)	:: sigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver/transmitter
   type(dvecMTX), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(EMsolnMTX), intent(inout), optional	:: eAll

   ! local variables
   type(EMsoln)				:: e0
   integer				:: iTx,j

   if(.not.d%allocated) then
      call errStop('data vector not allocated on input to fwdPred')
   end if

   if(present(eAll)) then
      if(.not. eAll%allocated) then
         !call errStop('eAll has to be allocated on input to fwdPred')
         call create_EMsolnMTX(d,eAll)
      else if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in fwdPred')
      endif
   endif

   ! loop over transmitters: solve forward system for each,
   !    apply (non-linear) data functionals
   ! write(6,*) 'd%nTx = ',d%nTx

   ! set errorBar=.false. since predicted data do not have
   !    well-defined error bars
   d%errorBar = .false.

   do j = 1,d%nTx
      !   get indices into dictionaries for this dvec
      iTx = d%d(j)%tx

      !  do any necessary initialization for transmitter iTx
      call initSolver(iTx,sigma,e0)

      ! compute forward solution
      call fwdSolve(iTx,e0)

      ! apply data functionals
      call dataMeas(e0,sigma,d%d(j))

      if(present(eAll)) then
         call copy_EMsoln(eAll%solns(j),e0)
      endif

   enddo

   ! deallocate and clean up
   call exitSolver(e0)

   end subroutine fwdPred
   !**********************************************************************

end module sensMatrix
