module sensComp

  use griddef
  use data_vectors
  use dataFunc
  use modelmap
  use jacobian
  use output
  use initFields
  use dataMisfit
  use boundaries
  use dataspace
  use SolnSpace

  implicit none

  public 	:: calcSensMatrix, Jmult, JmultT, fwdPred, setGrid, cleanUp, deall_sensMatrix

  ! numerical discretization used to compute the EM solution
  !  (may be different from the grid stored in model parameter)
  ! currently this is inherited from global
  ! type(grid_t), target, save, private     :: grid

  ! utility variables necessary to time the computations;
  !  including a private variable rtime that stores total run time
  real, save, private						:: rtime  ! run time
  real, save, private						:: ftime  ! run time per frequency
  real, save, private						:: stime, etime ! start and end times
  integer, dimension(8), private			:: tarray ! utility variable

  ! ***************************************************************************
  ! * type sensMatrix_t stores the full sensitivity matrix or the
  ! * partial derivatives with respect to the original model parameters
  ! * for all frequencies and all data functionals

  type :: sensMatrix_t

	type(modelParam_t), pointer, dimension(:,:)	 :: dm	  !(nfreq,nfunc)
	logical                                      :: allocated=.false.

  end type sensMatrix_t

Contains

   !**********************************************************************
   subroutine calcSensMatrix(d,sigma0,dsigma)
   !  Calculate sensitivity matrix for data in d
   !
   !   d is the input data vector, here just used to identify
   !     receiver transmitter pairs to compute sensitivities for
   type(dataVectorMTX_t), intent(in)	:: d
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the output array of data sensitivities,
   !   one for each element in the data array.  Each sensitivity
   !    is an element of type modelParam, an abstract
   !    data type that defines the unknow conductivity
   type(modelParam_t), pointer   :: dsigma(:)


   end subroutine calcSensMatrix

   !**********************************************************************
   subroutine Jmult(dsigma,sigma0,d,eAll)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for all transmitters in a datavector (i.e., multiple dataVec objects)
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
   type(dataVectorMTX_t), intent(inout)		:: d
   type(solnVectorMTX_t), intent(inout), optional	:: eAll



   end subroutine Jmult

   !**********************************************************************
   subroutine JmultT(m0,d,dm,H,dR)

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

   type(modelParam_t), intent(in)	:: m0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVectorMTX_t), intent(in)		:: d
   !   dsigma is the output conductivity parameter
   type(modelParam_t), intent(out)  	:: dm
   type(solnVectorMTX_t), intent(in), optional	:: H
   type(sensMatrix_t), intent(inout), optional :: dR

    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
    type (cvector)							:: Hj,B,F,Hconj,dH,dE,Econj
	type (rvector)							:: dE_real
	type (rscalar)							:: rho,drho
	type (sparsevecc)						:: Hb
	type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	type (modelParam_t)						:: dmisfit
	integer									:: ifreq,ifunc
	logical									:: adjoint,delta
	character(1)							:: cfunc
	character(80)							:: fn_err
	logical		:: savedSolns

	savedSolns = present(H)
	if(savedSolns) then
		if(d%nTx .ne. H%nTx) then
			call errStop('dimensions of H and d do not agree in JmultT')
		endif
	endif

	if (present(dR)) then
		if (.not. dR%allocated) then
			allocate(dR%dm(nfreq,nfunc),STAT=istat)
			dR%allocated = .true.
		endif
	endif

	dm = m0
	call zero(dm)

	! Temporary model parameter
	dmisfit = m0
	call zero(dmisfit)

    ! Allocate the resistivity vectors
	call create_rscalar(grid,rho,CENTER)
	call create_rscalar(grid,drho,CENTER)

	! Compute model information everywhere else in the domain
	call initModel(grid,m0,rho%v)

	! Start the (portable) clock
	call date_and_time(values=tarray)

	do ifreq=1,nfreq

	  stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  if(savedSolns) then
		Hj = H%solns(ifreq)
	  else
		call initialize_fields(Hj,B)

	  	write(6,*) 'Solving 3D forward problem for freq ',ifreq,freq%value

	  	! solve A <h> = <b> for vector <h>
	  	adjoint=.FALSE.
	  	call operatorM(Hj,B,omega,rho%v,grid,fwdCtrls,errflag,adjoint)

	  	! compute and output fields & C and D responses at cells
	  	call outputSolution(freq,Hj,slices,grid,cUserDef,rho%v,'h')
	  end if

	  !print *, 'Model initialized',H%solns(ifreq)%allocated,Hj%allocated

	  do ifunc=1,nfunc

		write(6,*)
		write(6,*) 'Multiplying by J^T for ',&
				  trim(TFList%info(ifunc)%name), ' responses...',ifreq,freq%value

		! $G_\omega r_\omega$
  		call operatorG(d%v(ifreq,ifunc,:),Hj,F)

		!print *, 'operator G successful',Hj%x


		! $M*^{-1}_{\rho,-\omega} ( G_\omega r_\omega )$
		adjoint = .TRUE.
		delta = .TRUE.
		! Forcing term F should not contain any non-zero boundary values
		call operatorM(dH,F,omega,rho%v,grid,fwdCtrls,errflag,adjoint,delta)
		!call create_cvector(grid,dH,EDGE)
		!dH%x = C_ONE
		!dH%y = C_ONE
		!dH%z = C_ONE
		! call outputSolution(freq,dH,slices,grid,cUserDef,rho,'dh')

		! Pre-divide the interior components of dH by elementary areas
		call operatorD_Si_divide(dH,grid)

		! $C D_{S_i}^{-1} M*^{-1}_{\rho,-\omega} ( G_\omega r_\omega )$
		call createBC(Hb,grid)
		call insertBC(Hb,dH)
		call operatorC(dH,dE,grid)

		!	cfunc = trim(TFList%info(ifunc)%name)
		!	fn_err = trim(outFiles%fn_err)//trim(cfunc)
		!	call initFileWrite(fn_err,ioERR)
		!	do i= 1,grid%nx
		!	  do j =1,grid%ny
		!		do k =1,grid%nz
		!		  if (dreal(dE%y(i,j,k)) > 1.0) then
		!			write(ioERR,*) i,j,k, dreal(dH%y(i,j,k)), dreal(dE%y(i,j,k))
		!		  end if
		!		end do
		!	  end do
		!	end do
		!	close(ioERR)

		! $\bar{\e} = C \bar{\h}$
		Hconj = conjg(Hj)
		call operatorD_l_mult(Hconj,grid)
		call operatorC(Hconj,Econj,grid)

		! $D_{\bar{\e}} C D_{S_i}^{-1} M*^{-1}_{\rho,-\omega} ( G_\omega r_\omega )$
		call diagMult(Econj,dE,dE) !dE = Econj * dE

		! $\Re( D_{\bar{\e}} C D_{S_i}^{-1} M*^{-1}_{\rho,-\omega} ( G_\omega r_\omega ) )$
		dE_real = real(dE)

		! $L^T \delta{R}$
		call operatorLt(drho%v,dE_real,grid)

		! $P^T L^T \delta{R}$
		call operatorPt(drho,dmisfit)

		! this line is needed to counter some mistake in the above... should be -dmisfit
		call scMult(MinusONE,dmisfit,dmisfit)

		! add to the total model parametrization
		call linComb(ONE,dm,ONE,dmisfit,dm)

		! save the misfit derivatives if required
		if (present(dR)) then
			dR%dm(ifreq,ifunc) = dmisfit
		end if

	  end do

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  print *
	  rtime = rtime + ftime

	end do

	call deall_modelParam(dmisfit)
	call deall_rscalar(rho)
	call deall_rscalar(drho)
	call deall_sparsevecc(Hb)
	call deall_cvector(Hj)
	call deall_cvector(B)
	call deall_cvector(F)
	call deall_cvector(dH)
	call deall_cvector(Hconj)
	call deall_cvector(dE)
	call deall_cvector(Econj)
	call deall_rvector(dE_real)

   end subroutine JmultT


   !**********************************************************************
   subroutine fwdPred(m,d,H)

   !  Calculate predicted data for dataVectorMTX object d
   !    and for conductivity parameter sigma
   !  Optionally returns array of EM solutions eAll, one for
   !    each transmitter (if eAll is present)
   !
   !  First need to set up transmitter, receiver, dataType dictionaries;
   !  "pointers" to dictionaries are attached to multi-transmitter
   !   data vector d .
   !
   !   sigma is input conductivity parameter
   type(modelParam_t), intent(in)	:: m
   !   d is the computed (output) data vector, also used to identify
   !     receiver/transmitter
   type(dataVectorMTX_t), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(solnVectorMTX_t), intent(inout), optional	:: H
   ! local variables
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
    type (cvector)							:: Hj,B,F,Hconj,B_tilde,dH,dE,Econj,Bzero,dR
	type (rvector)							:: dE_real
	type (rscalar)							:: rho
	type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	integer									:: ifreq,ifunc
	logical									:: adjoint,delta
	character(1)							:: cfunc
	character(80)							:: fn_err
	type(dataVectorMTX_t)                      :: dat,psi
	integer                                 :: nfreq,nfunc,nobs

    if(.not.d%allocated) then
       call errStop('data vector not allocated on input to fwdPred')
    end if

    if(present(H)) then
       if(.not. H%allocated) then
          call create_solnVectorMTX(d%nTx,H,grid)
       else if(d%nTx .ne. H%nTx) then
          call errStop('dimensions of H and d do not agree in fwdPred')
       end if
    end if

    ! Allocate temporary data and response storage
    nfreq = freqList%n
    nfunc = TFList%n
    nobs  = obsList%n

    dat = d
    psi = d
    call zero(psi)

    ! Allocate the resistivity vector
	call create_rscalar(grid,rho,CENTER)

	! Compute model information everywhere else in the domain
	call initModel(grid,m,rho%v)

	! Start the (portable) clock
	call date_and_time(values=tarray)

	call initialize_fields(Hj,B)

	do ifreq=1,freqList%n

	  stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  write(6,*) 'Solving 3D forward problem for freq ',ifreq,freq%value

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  call operatorM(Hj,B,omega,rho%v,grid,fwdCtrls,errflag,adjoint)
	  !call create_cvector(grid,Hj,EDGE)
	  !Hj%x = C_ONE
	  !Hj%y = C_ONE
	  !Hj%z = C_ONE

	  ! compute and output fields & C and D responses at cells
	  call outputSolution(freq,Hj,slices,grid,cUserDef,rho%v,'h')

	  ! output full H-field cvector
	  if (output_level > 3) then
	  	call outputField(freq,Hj,cUserDef,'field')
	  end if

	  ! compute and output C and D responses at observatories
	  call calcResponses(freq,Hj,dat,psi)
	  call outputResponses(freq,psi,freqList,TFList,obsList,outFiles,dat)

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  print *
	  rtime = rtime + ftime

	  if(present(H)) then
	     H%solns(ifreq) = Hj
	     H%tx(ifreq) = ifreq
	     H%errflag(ifreq) = errflag
	  end if

	end do

	d = psi

	call deall_dataVectorMTX(dat)
	call deall_dataVectorMTX(psi)
	call deall_rscalar(rho)
	call deall_cvector(Hj)
	call deall_cvector(B)

   end subroutine fwdPred

   !**********************************************************************
   subroutine setGrid(newgrid)

   !  Use to set and/or update the numerical grid, that is then used
   !   all computations in this module;
   !   This is not a pointer target.

   type(grid_t), intent(in)     :: newgrid
   ! local
   integer                      :: istat

    grid = newgrid

    ! Update the interpolation weights of observatories to the new list
    call initObsList(grid,obsList)

   	! Too complicated to rewrite input to all subroutines that use x,y,z,nx,ny,nz
	! in terms of the grid variable, but that would be the way to do it in the
	! future. Currently, use this patch instead. x,y,z stored in module griddef
	nx = grid%nx; ny = grid%ny; nz = grid%nz
	nzEarth = grid%nzEarth; nzAir = grid%nzAir
	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x; y = grid%y; z = grid%z

   end subroutine setGrid

   !**********************************************************************
   subroutine cleanUp()

   !  delete the private grid variable

   ! local
   integer                      :: istat

    call deall_grid(grid)

	deallocate(x,y,z,STAT=istat)

   end subroutine cleanUp

   !**********************************************************************
   subroutine deall_sensMatrix(dsigma)

   !  delete the private grid variable
   type (sensMatrix_t), intent(inout)   :: dsigma
   ! local
   integer                      :: i,j,istat

	if (dsigma%allocated) then
		do i = 1,size(dsigma%dm,1)
			do j = 1,size(dsigma%dm,2)
				call deall_modelParam(dsigma%dm(i,j))
			enddo
		enddo
	endif

	if (associated(dsigma%dm)) deallocate(dsigma%dm,STAT=istat)
	dsigma%allocated = .false.

   end subroutine deall_sensMatrix

end module sensComp
