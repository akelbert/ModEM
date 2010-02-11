module sensComp

  use math_constants
  use utilities
  use datasens	 !!!!  inherits : dataspace, dataFunc, solnrhs
  use modelsens  !!!  inherits : modelspace, soln2d
  use emsolver
  use sensmatrix

  implicit none

  public 	:: calcSensMatrix, setGrid, cleanUp
  public	:: Jmult,   Jmult_TX
  public	:: JmultT,  JmultT_TX
  public	:: fwdPred, fwdPred_TX

  ! numerical discretization used to compute the EM solution
  !  (may be different from the grid stored in model parameter)
  type (grid_t), target, save, private     :: grid

  ! temporary EM fields, that are saved for efficiency - to avoid
  !  memory allocation & deallocation for each transmitter
  type(EMsoln_t), save, private		:: e,e0
  type(EMrhs_t) , save, private		:: comb


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
   type(dataVecMTX_t), intent(in)	:: d
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the output array of data sensitivities,
   !   one for each element in the data array.  Each sensitivity
   !    is an element of type modelParam, an abstract
   !    data type that defines the unknow conductivity
   type(modelParam_t), pointer   :: dsigma(:)

   !  local variables
   type(dataVec_t)  :: dVec
   integer 		:: i,j,nTx,k,nSite,nTot,ii,iTx, &
				iDT,nfunc,ncomp,iRx,iFunc
   type(EMsparse_t), pointer	:: L(:),Q(:)
   logical 		:: calcQ

   ! nTot is number of real data
   nTot = countData(d)

   if(.not.associated(dsigma)) then
      ! allocate for sensitivity matrix
      allocate(dsigma(nTot))
      do j = 1,nTot
         ! this makes a copy of model param, of the same type
         !   as sigma0, then zeros it.
         dsigma(j) = sigma0
         call zero(dsigma(j))
      enddo
   endif

   ! nTX is number of transmitters
   nTx = d%nTx

   ! loop over frequencies, computing all sensitivities for
   !   one frequency
   ii = 0
   do j = 1,nTx

      ! identify the transmitter
      iTx = d%d(j)%tx

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,grid,e0,e,comb)

      !  solve forward problem; result is stored in e0
      call fwdSolve(iTx,e0)

      ! now loop over data types
      do i = 1,d%d(j)%nDt

        ! make a local copy of the dataVec for this transmitter and dataType
        dVec = d%d(j)%data(i)

        ! get index into dataType dictionary for this dataVec
        iDT = dVec%dataType
        calcQ = typeDict(iDT)%calcQ

        ! get data type info for this dataVec: here we allow for either
        !  real or complex data, though we always store as real
        !  Complex data come in pairs, stored as two succesive
        !  "components" in the dataVec structure.
        !  "isComplex" is an attribute of a dataVec structure; can
        !  be different for different dataVecs within the overall data
        !  vector
        nComp = dVec%nComp
        nSite = dVec%nSite
        if(dVec%isComplex) then
           !  data are complex; one sensitivity calculation can be
           !   used for both real and imaginary parts
           if(mod(nComp,2).ne.0) then
              call errStop('for complex data # of components must be even in CalcSensMatrix')
           endif
           nFunc = nComp/2
        else
           !  data are treated as real: full sensitivity computation is required
           !   for each component
           nFunc = nComp
        endif
        allocate(L(nFunc))
        allocate(Q(nFunc))

        ! loop over sites, computing sensitivity for all components for each site
        do k = 1,nSite

           iRx = dVec%rx(k)

           !  compute linearized data functional(s) : L
           call linDataFunc(e0,sigma0,iDT,iRx,L,Q)

           ! loop over functionals  (for TE/TM impedances nFunc = 1)
           do iFunc = 1,nFunc

              ! solve forward problem for each of nFunc functionals
              call zero_EMrhs(comb)
              call add_EMsparseEMrhs(C_ONE,L(iFunc),comb)

              call sensSolve(iTx,TRN,comb,e)

              ! multiply by P^T (and add Q^T if appropriate)
              if(dVec%isComplex) then
                 ii = ii + 2
                 call PmultT(e0,sigma0,e,dsigma(ii-1),dsigma(ii))
                 if(calcQ) then
                    call QaddT(C_ONE,Q(iFunc),sigma0,dsigma(ii-1),dsigma(ii))
                 endif
              else
                 ii = ii+1
                 call PmultT(e0,sigma0,e,dsigma(ii))
                 if(calcQ) then
                    call QaddT(C_ONE,Q(iFunc),sigma0,dsigma(ii))
                 endif
              endif

              ! deallocate temporary sparse vectors
              call deall_EMsparse(L(iFunc))
              call deall_EMsparse(Q(iFunc))

           enddo  ! iFunc
        enddo  ! sites

        deallocate(L)
        deallocate(Q)

      enddo  ! dataType's
   enddo  ! tx

   call deall_dataVec(dvec)

  end subroutine calcSensMatrix


  !**********************************************************************
  subroutine Jmult_TX(dsigma,sigma0,d,emsoln)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for one transmitter (but possibly for multiple data types)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  indicies into dictionaries are attached to data vector d .
   !
   !  If optional input parameter e1 is present, it must contain
   !   solutions for this transmitter for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions)
   !
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the input conductivity parameter perturbation
   type(modelParam_t), intent(in)	:: dsigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVecTX_t), intent(inout)		:: d
   type(EMsoln_t), intent(inout), optional	:: emsoln

   !  local variables
   integer 		:: i,j,iTx,iDT
   logical		:: savedSolns

      savedSolns = present(emsoln)
      if(savedSolns) then
         if(d%tx .ne. emsoln%tx) then
            call errStop('EM soln and d do not agree in Jmult')
         endif
      endif

      !   get index into transmitter dictionary for this dataVecTX
      iTx = d%tx

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,grid,e0,e,comb)

      if(savedSolns) then
         ! e0 = emsoln
         call copy_EMsoln(e0,emsoln)
      else
         ! solve forward problem; result is stored in e0
         call fwdSolve(iTx,e0)
      endif

      !  compute rhs (stored in comb) for forward sensitivity
      !  calculation, using conductivity perturbations and
      !  background soln:
      call Pmult(e0,sigma0,dsigma,comb)

      ! solve forward problem with source in comb
      call sensSolve(iTx,FWD,comb,e)

      ! now, loop over data types
      do i = 1,d%nDt
         ! get index into the dataType dictionary
         iDT = d%data(i)%dataType
         ! finally apply linearized data functionals
         if(typeDict(iDT)%calcQ) then
            call linDataMeas(e0,sigma0,e,d%data(i),dsigma)
         else
	        call linDataMeas(e0,sigma0,e,d%data(i))
         endif
      enddo

  end subroutine Jmult_TX


  !**********************************************************************
  subroutine Jmult(dsigma,sigma0,d,eAll)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for all transmitters in a datavector (i.e., multiple dataVec objects)
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
   !   dsigma is the input conductivity parameter perturbation
   type(modelParam_t), intent(in)	:: dsigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVecMTX_t), intent(inout)		:: d
   type(EMsolnMTX_t), intent(inout), optional	:: eAll

   !  local variables
   integer 		:: j
   logical		:: savedSolns

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in Jmult_MTX')
      endif
   endif

   ! loop over frequencies : solve forward system twice,
   !    to compute background, perturbation solutions
   !    apply data functionals (linearized about background soln)
   !    to perturbation solution
   do j = 1,d%nTx

    if(savedSolns) then
   	   ! compute d = J x m for a single transmitter
   	   call Jmult_TX(dsigma,sigma0,d%d(j),eAll%solns(j))
    else
       ! do not pass the EM soln to Jmult - it will be computed
   	   call Jmult_TX(dsigma,sigma0,d%d(j))
    endif

   enddo  ! tx


  end subroutine Jmult


  !**********************************************************************
  subroutine JmultT_TX(sigma0,d,dsigma,emsoln)

   !  Transpose of Jmult mujltiplied by data vector d for one transmitter;
   !   output is a single conductivity parameter in dsigma
   !
   !   sigma0 is background conductivity parameter
   !
   !  If optional input parameter e0 is present, it must contain
   !   solutions for this transmitter for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions,
   !    e.g., in a CG solution scheme)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  "pointers" to dictionary entries are attached to data vector d

   type(modelParam_t), intent(in)	:: sigma0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVecTX_t), intent(in)		:: d
   !   dsigma is the output conductivity parameter
   type(modelParam_t), intent(inout)  	:: dsigma
   type(EMsoln_t), intent(in), optional	:: emsoln

   !  local variables
   type(modelParam_t)	:: sigmaTemp, Qcomb
   integer 		:: i,j,iTx,iDT
   logical		:: calcSomeQ, firstQ
   logical		:: savedSolns

      calcSomeQ = .false.
      firstQ = .true.
      savedSolns = present(emsoln)
      if(savedSolns) then
         if(d%tx .ne. emsoln%tx) then
            call errStop('EM soln and d do not agree in JmultT')
         endif
      endif

      dsigma = sigma0
      call zero(dsigma)

      sigmaTemp = sigma0
      call zero(sigmaTemp)

      !   get index into transmitter dictionary for this dataVecTX
      iTx = d%tx

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,grid,e0,e,comb)

      if(savedSolns) then
         ! e0 = e1
         call copy_EMsoln(e0,emsoln)
      else
         ! solve forward problem; result is stored in e0
         call fwdSolve(iTx,e0)
      endif

      ! now, loop over data types
      do i = 1,d%nDt

         ! get index into dataType dictionary
         iDT = d%data(i)%dataType

         ! set up comb using linearized data functionals
         !  ... for dataVecs with data functionals depending on conductivity
         !   parameter also compute analagous comb in parameter space
         if(typeDict(iDT)%calcQ) then
            if(firstQ) then
               !  first transmitter for which Q must be calculated:
               !   ==> allocate and zero Qcomb (use copy so that paramtype
               !         is set correctly)
               !  NOTE: linDataComb ADDS to Qcomb, not overwrites
               !   ==> only zero Qcomb for first transmitter requiring Q
               Qcomb = sigmaTemp
               call zero(Qcomb)
               !  set flags indicating that Q is now non-zero
               calcSomeQ = .true.
               firstQ = .false.
            endif
            !  BUT: linDataComb overwrites comb ... so zero this
            !       for every transmitter
            call zero_EMrhs(comb)
            call linDataComb(e0,sigma0,d%data(i),comb,Qcomb)
         else
            call zero_EMrhs(comb)
            call linDataComb(e0,sigma0,d%data(i),comb)
         endif
         ! solve forward problem with source in comb
         call sensSolve(iTx,TRN,comb,e)

         ! map from nodes to conductivity space ... result in sigmaTemp
         !  HERE PmultT overwrites sigmaTemp
         !  NOTE: here we throw away imaginary part, even for complex
         !     data (conceivably might want to save this in some cases!)
         call PmultT(e0,sigma0,e,sigmaTemp)
         call linComb_modelParam(ONE,dsigma,ONE,sigmaTemp,dsigma)

      enddo  ! dataType's

	  ! In theory, linDataComb should compute Qcomb for each transmitter
	  ! and data type, and add it to dsigma; however, this is currently
	  ! implemented such that linDataComb *adds* to Qcomb, each time it
	  ! is called. So for now, we add the total Qcomb once at the end.
      if(calcSomeQ) then
         !  add Qcomb
         call linComb_ModelParam(ONE,dsigma,ONE,Qcomb,dsigma)
      endif

      !  clean up
      call deall_modelParam(Qcomb)
      call deall_modelParam(sigmaTemp)

  end subroutine JmultT_TX


  !**********************************************************************
  subroutine JmultT(sigma0,d,dsigma,eAll)

   !  Transpose of Jmult multiplied by data vector d for all transmitters;
   !   output is a single conductivity parameter in dsigma
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
   type(dataVecMTX_t), intent(in)		:: d
   !   dsigma is the output conductivity parameter
   type(modelParam_t), intent(out)  	:: dsigma
   type(EMsolnMTX_t), intent(in), optional	:: eAll

   !  local variables
   type(modelParam_t)	:: sigmaTemp
   integer 		:: j
   logical		:: savedSolns

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in JmultT_MTX')
      endif
   endif

   dsigma = sigma0
   call zero(dsigma)

   ! loop over transmitters
   do j = 1,d%nTx

    if(savedSolns) then
   	   ! compute sigmaTemp = JT x d for a single transmitter
   	   call JmultT_TX(sigma0,d%d(j),sigmaTemp,eAll%solns(j))
    else
       ! do not pass the EM soln to JmultT - it will be computed
   	   call JmultT_TX(sigma0,d%d(j),sigmaTemp)
    endif
   	! ... add to dsigma
   	call linComb_modelParam(ONE,dsigma,ONE,sigmaTemp,dsigma)

   enddo  ! tx

   !  clean up
   call deall_modelParam(sigmaTemp)

  end subroutine JmultT


  !**********************************************************************
  subroutine fwdPred_TX(sigma,d,emsoln)

   !  Calculate predicted data for single-transmitter dataVecTX object d
   !    and for conductivity parameter sigma
   !  Also returns the EM solution e for this transmitter
   !
   !  First need to set up transmitter, receiver, dataType dictionaries;
   !  "pointers" to dictionaries are attached to data vector d .
   !
   !   sigma is input conductivity parameter
   type(modelParam_t), intent(in)	:: sigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver/transmitter
   type(dataVecTX_t), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(EMsoln_t), intent(inout)	:: emsoln

   ! local variables
   integer				:: iTx,i,j

      if(.not.d%allocated) then
         call errStop('data vector not allocated on input to fwdPred')
      end if

      ! get index into transmitter dictionary for this dataVecTX
      iTx = d%tx

      !  do any necessary initialization for transmitter iTx
      call initSolver(iTx,sigma,grid,emsoln)

      ! compute forward solution
      call fwdSolve(iTx,emsoln)

      ! cycle over data types
      do i = 1,d%nDt

         ! set errorBar=.false. since predicted data do not have
         ! well-defined error bars (important for the inversion)
         d%data(i)%errorBar = .false.

         ! apply data functionals
         call dataMeas(emsoln,sigma,d%data(i))

      enddo

  end subroutine fwdPred_TX


  !**********************************************************************
  subroutine fwdPred(sigma,d,eAll)

   !  Calculate predicted data for dataVecMTX object d
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
   type(dataVecMTX_t), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(EMsolnMTX_t), intent(inout), optional	:: eAll

   ! local variables
   integer				:: j

   if(.not.d%allocated) then
      call errStop('data vector not allocated on input to fwdPred')
   end if

   if(present(eAll)) then
      if(.not. eAll%allocated) then
         call create_EMsolnMTX(d%nTx,eAll)
      else if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in fwdPred')
      endif
   endif

   ! loop over transmitters: solve forward system for each,
   !    apply (non-linear) data functionals
   do j = 1,d%nTx

      call fwdPred_TX(sigma,d%d(j),e0)

      if(present(eAll)) then
         call copy_EMsoln(eAll%solns(j),e0)
      endif

   enddo

  end subroutine fwdPred


  !**********************************************************************
  subroutine setGrid(newgrid)

   !  Use to set and/or update the numerical grid, that is then used
   !   all computations in this module;
   !   This is not a pointer target.
   !  Might also have to run exitSolver at this point, if we are updating
   !   the grid during an inversion; that restarts the EMsolver module.

   type(grid_t), intent(in)     :: newgrid

   grid = newgrid

  end subroutine setGrid


  !**********************************************************************
  subroutine cleanUp()

   ! Subroutine to deallocate all memory stored in this module

   call exitSolver(e0,e,comb)
   call deall_grid(grid)

  end subroutine cleanUp

end module sensComp
