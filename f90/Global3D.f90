
program earth

  use main
  use output
  use senscomp
  use nlcg

  use global ! for symmetry testing only

  ! testing
  use ringcurrent
  use dimensions
  use field_vectors
  use data_vectors
  implicit none

  real										:: runtime
  character(80)								:: fn_startup=''
  character(80)								:: fn_model='rho.out'
  integer									:: i,j,k,ios,ii
  integer									:: istat
  real(8),dimension(:),allocatable			:: da, R1, R2
  real(8),dimension(:,:),allocatable			:: dR
  complex(8),dimension(:,:,:),allocatable	:: psi1,psi2,dpsi1,dpsi2
  integer									:: ifreq,ifunc,iobs
  real(8), dimension(2)                         :: value1, value2, delta
  real(8)                                       :: eps

  runtime = 0.0d0

  fn_startup='fwd_startup'
  !fn_startup=''
  eps = 0.01

  call InitGlobalData(fn_startup)

  call outputModel(outFiles%fn_model,grid,rho)

  call InitGlobalArrays()

  if (cUserDef%calculate == 'nothing') then
    call write_modelParam(param,'testOutput.prm')
	call DeleteGlobalData()
	call DeleteGlobalArrays()
	print *,'Successfully initialized the model. Exiting...'
	stop
  end if

  ! Model initialisation complete, and boundary conditions have been set.
  !	Now solve forward problem for each frequency, in turn, and for each
  ! compute the derivative of the penalty functional or the Jacobian,
  !	if required

  select case ( cUserDef%calculate )

  case ('inverse')

	print *,'Running NLCG inversion.'

	call calc_inverse(runtime)

  case ('original')

	print *,'Calculating responses using the original forward solver.'

	call calc_original(runtime)

  case ('secondary')

	print *,'Calculating solution using the secondary field formulation.'

	call calc_secondary_fields(runtime)

  case ('responses')

	print *,'Calculating responses only.'

	call calc_responses(runtime)

  case ('derivative')

	print *,'Calculating responses and the penalty functional derivative.'

	call calc_derivative(runtime)

  case ('jacobian')

	print *,'Calculating responses and the full Jacobian.'

	call calc_jacobian(runtime)

  case ('symmetry')

	print *,'Perform the symmetry test and exit.'

	!call calc_symmetry(runtime)
	call calc_symmetric_operators(runtime)

  case ('test_derivative')

	print *,'Run a test to validate derivative computations.'

	! this includes 'frozen' variables, but we're just testing
	allocate(da(param%nc))
	call random_number(da)
	da = da * 0.005

    open(unit=50,file='da.dat',status='unknown')
	write(50,*) 'da = ',da
	close(50)

	allocate(dR(TFList%n,param%nc),R1(TFList%n),R2(TFList%n))

  ! Run Jacobian computations for the unperturbed parametrization (a)
	call calc_derivative(runtime)
        R1 = misfitValue
        dR = dmisfitValue

	call DeleteGlobalData()
	call DeleteGlobalArrays()

	call InitGlobalData(fn_startup,da)
	call InitGlobalArrays()

  ! Run Jacobian computations after the perturbation (a+da)
	call calc_derivative(runtime)
        R2 = misfitValue

        do ifunc=1,TFList%n
           write(0,'(a40)') trim(TFList%info(ifunc)%name)//'-response: R1-R2, dR'
           write(0,'(2g15.7)') R2(ifunc)-R1(ifunc),dot_product(da(:),dR(ifunc,:))
        end do


        deallocate(da,dR,R1,R2)

  case ('validate')

	print *,'Run full test to validate Jacobian computations.'

	allocate(psi1(nfreq,nfunc,nobs),dpsi1(nfreq,nfunc,nobs))
	allocate(psi2(nfreq,nfunc,nobs),dpsi2(nfreq,nfunc,nobs))

	allocate(da(param%nc))
	call random_number(da)
	da = da * eps
   write(0,*) 'Model perturbation used is uniform random times ',eps
   open(unit=50,file='da.dat',status='unknown')
	write(50,*) 'da = ',da
	close(50)

  ! Run Jacobian computations for the unperturbed parametrization (a)
	call calc_jacobian(runtime)

	write(0,*) 'Data functionals 1: written into file psi1.dat'
    open(unit=51,file='psi1.dat',status='unknown')
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  write(51,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,value: ',ifreq,ifunc,iobs,&
				psi%v(ifreq,ifunc,iobs)%resp%value
		end do
	  end do
	end do
	close(51)
	write(0,*) 'Derivative 1 = Jacobian x perturbation: written into file dpsi1.dat'
    open(unit=52,file='dpsi1.dat',status='unknown')
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  dpsi1(ifreq,ifunc,iobs) = dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		  write(52,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,delta: ',ifreq,ifunc,iobs,&
				dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		end do
	  end do
	end do
	close(52)

	write(0,*) 'Data functionals 1: '
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  write(0,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,value: ',ifreq,ifunc,iobs,&
				psi%v(ifreq,ifunc,iobs)%resp%value
		end do
	  end do
	end do
	write(0,*) 'Derivative 1 = Jacobian x perturbation: '
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  dpsi1(ifreq,ifunc,iobs) = dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		  write(0,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,delta: ',ifreq,ifunc,iobs,&
				dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		end do
	  end do
	end do

	psi1 = psi%v(:,:,:)%resp%value

	call DeleteGlobalData()
	call DeleteGlobalArrays()

	call InitGlobalData(fn_startup,da)
	call InitGlobalArrays()

  ! Run Jacobian computations after the perturbation (a+da)
	call calc_jacobian(runtime)

	write(0,*) 'Data functionals 2: written into file psi2.dat'
    open(unit=51,file='psi2.dat',status='unknown')
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  write(51,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,value: ',ifreq,ifunc,iobs,&
				psi%v(ifreq,ifunc,iobs)%resp%value
		end do
	  end do
	end do
	close(51)
	write(0,*) 'Derivative 2 = Jacobian x perturbation: written into file dpsi2.dat'
    open(unit=52,file='dpsi2.dat',status='unknown')
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  dpsi1(ifreq,ifunc,iobs) = dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		  write(52,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,delta: ',ifreq,ifunc,iobs,&
				dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		end do
	  end do
	end do
	close(52)

	write(0,*) 'Data functionals 2: '
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  write(0,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,value: ',ifreq,ifunc,iobs,&
				psi%v(ifreq,ifunc,iobs)%resp%value
		end do
	  end do
	end do
	write(0,*) 'Derivative 2 = Jacobian x perturbation: '
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
		  dpsi2(ifreq,ifunc,iobs) = dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		  write(0,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,delta: ',ifreq,ifunc,iobs,&
				dot_product(da(:),sens%da(ifreq,ifunc,iobs,:))
		end do
	  end do
	end do

	psi2 = psi%v(:,:,:)%resp%value


	write(0,*) 'Psi2 - Psi1 / derivative (should be 1): '
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
                  value1(1) = dreal(psi1(ifreq,ifunc,iobs))
                  value1(2) = dimag(psi1(ifreq,ifunc,iobs))
                  value2(1) = dreal(psi2(ifreq,ifunc,iobs))
                  value2(2) = dimag(psi2(ifreq,ifunc,iobs))
                  delta(1) = dreal(dpsi1(ifreq,ifunc,iobs))
                  delta(2) = dimag(dpsi1(ifreq,ifunc,iobs))
		  write(0,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,value: ',ifreq,ifunc,iobs,&
                              (value2-value1)/delta
		end do
	  end do
	end do

	write(0,*) 'Psi2 - Psi1 / avg_derivative (should be 1): '
	dpsi1 = (dpsi1+dpsi2)/2
	do ifreq=1,nfreq
	  do ifunc=1,nfunc
		do iobs=1,nobs
                  value1(1) = dreal(psi1(ifreq,ifunc,iobs))
                  value1(2) = dimag(psi1(ifreq,ifunc,iobs))
                  value2(1) = dreal(psi2(ifreq,ifunc,iobs))
                  value2(2) = dimag(psi2(ifreq,ifunc,iobs))
                  delta(1) = dreal(dpsi1(ifreq,ifunc,iobs))
                  delta(2) = dimag(dpsi1(ifreq,ifunc,iobs))
		  write(0,'(a30,3i8,2g15.7)') 'ifreq,ifunc,iobs,value: ',ifreq,ifunc,iobs,&
                              (value2-value1)/delta
		end do
	  end do
	end do

	deallocate(psi1,dpsi1,psi2,dpsi2)

  case default

	print *,'Please specify exactly which quantities you would like to calculate.'
    print *,'Exiting...'

  end select


  continue

  write(0,*) ' elapsed time (mins) ',runtime/60.0

  call CloseOutFiles()

  call DeleteGlobalArrays()
  call DeleteGlobalData()


end program earth


  ! ***************************************************************************
  ! * calc_inverse runs the NLCG inversion to produce an inverse model

  subroutine calc_inverse(rtime)

  	use global
  	use nlcg

	real, intent(inout)						:: rtime  ! run time
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    type (modelParam_t)                     :: invparam ! inverse model

	! Start the (portable) clock
	call date_and_time(values=tarray)
	stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

  	call NLCGsolver(dat,cUserDef%damping,param,invparam,cUserDef%fn_invctrl)

	call date_and_time(values=tarray)
	etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	rtime = etime - stime

	call deall_modelParam(invparam)

  end subroutine calc_inverse


  ! ***************************************************************************
  ! * calc_original is a subroutine to calculate specified responses only for
  ! * the specified model, list of frequencies and possibly geometric locations,
  ! * output them and exit; we use the original forward solver by Uyeshima,
  ! * Schultz, Toh & Fujii

  subroutine calc_original(rtime)

	use global
	use field_vectors
	use ringcurrent
	use maxwells
	use output_orig
	use transmitters
  	implicit none

	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency

    complex(8)                              :: cdat,ddat  ! responses
	complex(8),allocatable,dimension(:)		:: vectorh,vectorb

	! These are the boundary conditions, that are never changed from
	! their original values. However, if Hx,Hy,Hz are *not* global,
	! the boundary conditions need to be explicitly inserted in them
	! before running SolveMaxwells - this is necessary to keep the
	! divergence correction valid. In theory, divergence correction
	! should not depend on the boundary conditions at all. However,
	! due to the bug in the realization, it does depend on z-component.
	! Currently, no work with the boundary values is required since
	! hx,hy,hz are indeed global.
	! type (sparsevecc)						:: bvH

    ! Average C response defined here
	complex(8),allocatable,dimension(:)		:: avg_cresp
	integer									:: istat,ifreq

	! Start the (portable) clock
	call date_and_time(values=tarray)

  ! Dynamic allocation of all global vectors
	allocate(vectorh(np2),vectorb(np2),STAT=istat)

  ! Initialize all global vectors here
	vectorh = C_ZERO
	vectorb = C_ZERO

    !call initialize_fields_vec(vectorh,vectorb,bvH)
	allocate(hx(np1),hy(np1),hz(np1),STAT=istat)
	! Initialise the H field at the outer boundary, then make it
	! depend linearly on radius throughout the computational domain
	! such that all field components are zero at the CMB. The P10
	! source at the boundary has been modified to the more accurate
	! initialization by Ikuko Fujii allowing for the ringcurrent.
	call set_initial_p10(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	!call set_initial_ringcurrent(nx,ny,nz,x,y,z,hx,hy,hz,np1)
	!call initBoundaryValues(bvH,hx,hy,hz,grid)

	! initialize vector <h> of A <h> = <b>
	call copyd3_d1_b(nx,ny,nz,hx,hy,hz,vectorh,x,y,z)

	! compute <b> of A <h> = <b>
	call calcb(nx,ny,nz,hx,hy,hz,vectorb,rho,x,y,z)


	! ---- Initialise I/O

	! call initFileWrite(outFiles%fn_jxjyjz,ioJ)
	! call initFileWrite(outFiles%fn_err,ioERR)
	! call initFileWrite(outFiles%fn_hxhyhz,ioH)
	call initFileWrite(outFiles%fn_cresp,ioC)
	call initFileWrite(outFiles%fn_dresp,ioD)
	call initFileWrite(outFiles%fn_hx,ioHx)
	call initFileWrite(outFiles%fn_hy,ioHy)
	call initFileWrite(outFiles%fn_hz,ioHz)

	allocate(avg_cresp(nfreq),STAT=istat)
	avg_cresp = C_ZERO

	do ifreq=1,freqList%n

	  stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	  omega  = 2.0d0*pi*freqList%info(ifreq)%value     ! angular frequency (radians/sec)

	  write(6,*) 'Solving 3D forward problem for freq ',ifreq,freqList%info(ifreq)%value

	  ! Initialize the fields required for divergence computations
	  !hx = (0.0d0,0.0d0)
	  !hy = (0.0d0,0.0d0)
	  !hz = (0.0d0,0.0d0)
	  !call insertBoundaryValues(bvH,hx,hy,hz) ! only z-component used - see comments

	  ! solve A <h> = <b> for vector <h>
	  call mult_vec_by_l(nx,ny,nz,vectorh,x,y,z)
	  call SolveMaxwells(vectorh,vectorb,omega,rho,fwdCtrls,errflag)

	  ! <h> = l*(Hx,Hy,Hz) -> <Hx>,<Hy>,<Hz>
	  call divide_vec_by_l(nx,ny,nz,vectorh,x,y,z)
	  call copyd1_d3_b(nx,ny,nz,hx,hy,hz,vectorh,x,y,z)
	  !call insertBoundaryValues(bvH,hx,hy,hz)

	  ! Controlled output
	  if (errflag==0) then

		 !call allh_out(nx,ny,nz,Hx,Hy,Hz,x,y,z,nzAir)
		 call calc_response(omega,nx,ny,nz,hx,hy,hz,x,y,z,nzAir,cdat,ddat)
		 !call allj_out(ioJ,nx,ny,nz,nzAir,omega,x,y,z,rho,hx,hy,hz)

	  else if (errflag==1) then

		 call calc_response(omega,nx,ny,nz,hx,hy,hz,x,y,z,nzAir,cdat,ddat)
		 print*,'error not reduced to specified level by ',fwdCtrls%nrelmax, &
				 ' iterations'
		 print*,'fields are being output, and program will continue'

	  else

   		 print *,'Forward solver has terminated abnormally at freq=',freqList%info(ifreq)%value
		 print *,'Exiting...'
		 exit

	  end if


	  avg_cresp(ifreq) = cdat

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  rtime = rtime + ftime

	end do

	continue

    call avg_cresp_out(freqList%info%value,nfreq,outFiles%fn_avg_cresp,avg_cresp)
	deallocate(avg_cresp,STAT=istat)
	deallocate(vectorh,vectorb,STAT=istat)
	deallocate(hx,hy,hz)

  end subroutine calc_original	! calc_original


  ! ***************************************************************************
  ! * calc_responses is a subroutine to calculate specified responses only for
  ! * the specified model, list of frequencies and possibly geometric locations,
  ! * output them and exit.

  subroutine calc_responses(rtime)

	use griddef
	use data_vectors
	use global
	use dataFunc
	use jacobian
	use output
	use initFields
	use dataMisfit
	use boundaries	!for testing
	use senscomp
    use transmitters
	implicit none

	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
	type (solnVectorMTX_t)                      :: H
    !type (cvector)							:: H,B,F
	!type (sparsevecc)						:: Hb
	!type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	integer									:: ifreq
	!logical									:: adjoint

	! Start the (portable) clock
	call date_and_time(values=tarray)
	stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

    ! Call the forward solver for all frequencies
    psi = dat
    call fwdPred(param,psi,H)

	do ifreq=1,freqList%n

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  call calcResiduals(freq,dat,psi,res)

	  ! If computing different kinds of misfits, be consistent;
	  ! write different kinds into different data structures
	  call calcMisfit(freq,res,misfit,misfit%name)

	!call OutputResiduals(freq,res,TFList,obsList)

	end do

	if (output_level>0) then
	write(0,*)
  	do i=1,freqList%n
	  do j=1,TFList%n
    write(6,'(a12,i3,a2,a35,i2,a2,g15.7)') 'Misfit for ',misfit%ndat(i,j),&
			trim(TFList%info(j)%name),' responses, computed for frequency ',i,' :',&
			misfit%value(i,j)/(2*misfit%ndat(i,j))
	  end do
	end do
	end if

	call misfitSumUp(res,misfit,misfitValue)

	!call outputMisfit(param,misfit,misfitValue,cUserDef)

	continue
	call deall_solnVectorMTX(H)

	call date_and_time(values=tarray)
	etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	rtime = etime - stime


  end subroutine calc_responses	! calc_responses


  ! ***************************************************************************
  ! * calc_secondary_fields is a subroutine to calculate fields and responses
  ! * using secondary field formulation; reading the primary field (1D solution)
  ! * from a file, or computing it internally (not implemented yet).

  subroutine calc_secondary_fields(rtime)

	use griddef
	use data_vectors
	use global
	use dataFunc
	use jacobian
	use input
	use output
	use initFields
	use dataMisfit
	use boundaries	!for testing
	use senscomp
    use transmitters
    use dataTypes
	implicit none

	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
	type (solnVectorMTX_t)						:: H1D
	type (cvector)                     	    :: Hj, Bj, dH, F
	type (rvector)							:: drhoF
	type (rscalar)							:: rho1D, drho
	type (modelParam_t)                     :: param1D
    !type (cvector)							:: H,B,F
	!type (sparsevecc)						:: Hb
	!type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	integer									:: ifreq
	logical									:: adjoint, delta

	! Start the (portable) clock
	call date_and_time(values=tarray)
	stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	! Don't need to initialize fields: initial dH and boundary conditions are zero
	! call initialize_fields(dH,Bj)

	! Make 1D parametrization out of full, and map to grid (primary cell centers)
	param1D = getRadial(param)
	call create_rscalar(grid,rho1D,CENTER)
	call initModel(grid,param1D,rho1D%v)

	! Take the difference on the grid, to avoid the problem with zero resistivity
	! call linComb_rscalar(ONE,rho,MinusONE,rho1D,drho)
	call create_rscalar(grid,drho,CENTER)
	drho%v = rho - rho1D%v

	! Map the resistivity vector to primary cell faces (dual edges)
	call operatorL(drho%v,drhoF,grid)

	! Read the primary field solution, computed for a 1D model
	! (alternatively, compute it locally)
	call initField(cUserDef,grid,H1D)

	do ifreq=1,freqList%n

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  ! Compute the RHS = - del x drho (del x H)
	  call operatorlC(H1D%solns(ifreq),F,grid)
	  call diagMult(drhoF,F,F)
	  call operatorCt(F,Bj,grid)
	  call operatorD_Si_divide(Bj,grid)
	  call linComb(C_MinusONE,Bj,C_ZERO,Bj,Bj)

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  delta=.TRUE.
	  call operatorM(dH,Bj,omega,rho,grid,fwdCtrls,errflag,adjoint,delta)

	  ! Full solution for one frequency is the sum H1D + dH
	  Hj = dH
	  call linComb_cvector(C_ONE,H1D%solns(ifreq),C_ONE,dH,Hj)

	  ! compute and output fields; C and D responses at cells
	  call outputSolution(freq,Hj,slices,grid,cUserDef,rho,'h')

	  ! output full H-field cvector
	  if (output_level > 3) then
	  	call outputField(freq,Hj,cUserDef,'field')
	  end if

	  ! compute and output C and D responses at observatories
	  call calcResponses(freq,Hj,dat,psi)
	  call outputResponses(freq,psi,freqList,TFList,obsList,outFiles,dat)

	  call calcResiduals(freq,dat,psi,res)

	  ! If computing different kinds of misfits, be consistent;
	  ! write different kinds into different data structures
	  call calcMisfit(freq,res,misfit,misfit%name)

	!call OutputResiduals(freq,res,TFList,obsList)

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  rtime = rtime + ftime

	end do

	if (output_level>0) then
	write(0,*)
  	do i=1,freqList%n
	  do j=1,TFList%n
    write(6,'(a12,i3,a2,a35,i2,a2,g15.7)') 'Misfit for ',misfit%ndat(i,j),&
			trim(TFList%info(j)%name),' responses, computed for frequency ',i,' :',&
			misfit%value(i,j)/(2*misfit%ndat(i,j))
	  end do
	end do
	end if

	call misfitSumUp(res,misfit,misfitValue)

	!call outputMisfit(param,misfit,misfitValue,cUserDef)

	continue
	call deall_cvector(Hj)
	call deall_cvector(Bj)
	call deall_cvector(dH)
	call deall_cvector(F)
	call deall_rvector(drhoF)
	call deall_rscalar(drho)
	call deall_rscalar(rho1D)
	call deall_modelParam(param1D)
	call deall_solnVectorMTX(H1D)

	call date_and_time(values=tarray)
	etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	rtime = etime - stime


  end subroutine calc_secondary_fields	! calc_secondary_fields


  ! ***************************************************************************
  ! * calc_derivative is a subroutine to calculate specified responses and the
  ! * full derivative of the penalty functional with respect to the original
  ! * model parameters for the list of frequencies, and exit. Designed for use
  ! * by gradient-based inverse solvers.

  subroutine calc_derivative(rtime)

	use data_vectors
	use global
	use jacobian
	use boundaries
	use output
	use initFields
	use dataFunc
	use dataMisfit
	use senscomp
    use transmitters
    use dataTypes
	implicit none

	type (sensMatrix_t)					    :: dR
	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
	type (solnVectorMTX_t)                      :: H
    !type (cvector)							:: H,B,F,Hconj,B_tilde,dH,dE,Econj,Bzero,dR
	!type (rvector)							:: dE_real
	!type (rscalar)							:: drho
	!type (sparsevecc)						:: Hb
	!type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	type (modelParam_t)						:: dmisfit,dmisfitSmooth
	integer									:: ifreq,ifunc
	logical									:: adjoint,delta
	character(1)							:: cfunc
	character(80)							:: fn_err
	real(8)									:: SS, RMS
	integer									:: Ndata

	! Start the (portable) clock
	call date_and_time(values=tarray)
	stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

    ! Call the forward solver for all frequencies
    psi = dat
    call fwdPred(param,psi,H)

	do ifreq=1,freqList%n

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  ! res = psi-dat (note: multiply by +2 to obtain derivative)
	  call calcResiduals(freq,dat,psi,res)
	  !call OutputResiduals(freq,res,TFList,obsList)

	  ! If computing different kinds of misfits, be consistent;
	  ! write different kinds into different data structures
	  !call calcMisfit(freq,res,misfit,misfit%name)

	  ! compute the weighted residuals and start the derivative computations
	  call calcResiduals(freq,dat,psi,wres,weighted=.TRUE.)

	end do

!	! initialize res
!	res = psi
!	! compute residual: res = dat-psi
!	call linComb(ONE,dat,MinusONE,psi,res)
!	! normalize residuals, compute sum of squares
!	wres = res
!	call normalizeData(wres,2)
!	SS = dotProd(res,wres)
!	Ndata = countData(res)
!	RMS = sqrt(SS/Ndata)
!	print *, 'Total RMS squared misfit for ',Ndata,' data values = ',RMS**2

	! Call multiplication by J^T
	call JmultT(param,wres,dmisfit,H,dR)

	!dmisfit = -2. * dmisfit
	call scMult(TWO,dmisfit,dmisfit)

	dmisfitSmooth = multBy_CmSqrt(dmisfit)

	! We can no longer store misfit%dRda - misfit derivative for each frequency,
	! functional and model parameter. They are now summed up internally in JmultT.
	! We might want to add this functionality back!!!
	do ifreq=1,freqList%n
	  do ifunc=1,nfunc

		call scMult(TWO,dR%dm(ifreq,ifunc),dR%dm(ifreq,ifunc))
	    dmisfitSmooth = multBy_CmSqrt(dR%dm(ifreq,ifunc))

		call getParamValues_modelParam(dmisfitSmooth,misfit%dRda(ifreq,ifunc,:))

	  end do
	end do

	do ifreq=1,freqList%n

	  freq = freqList%info(ifreq)

	  ! If computing different kinds of misfits, be consistent;
	  ! write different kinds into different data structures
	  call calcMisfit(freq,res,misfit,misfit%name)

	end do

	! Output the summary and exit
	call misfitSumUp(res,misfit,misfitValue,dmisfitValue)
	!call misfitSumUp(res,misfit,misfitValue)

	if (output_level>0) then
	write(0,*)
  	do i=1,freqList%n
	  do j=1,TFList%n
    write(6,'(a12,i3,a3,a35,i2,a2,g15.7)') 'Misfit for ',misfit%ndat(i,j),&
			trim(TFList%info(j)%name),' responses, computed for frequency ',i,' :',&
			misfit%value(i,j)/(2*misfit%ndat(i,j))
	  end do
	end do
	end if

	!do i=1,freqList%n
	!  do j=1,TFList%n
	!	write(0,'(a40,2i3,i6,g15.7)') 'ifreq,ifunc,ndat,misfit = ',i,j,&
	 ! 				misfit%ndat(i,j),misfit%value(i,j)/(2*misfit%ndat(i,j))
	!  end do
	!end do

	!call outputMisfit(param,misfit,misfitValue,cUserDef)
	!call outputDerivative(param,misfit,dmisfitValue,cUserDef)

	continue
	call deall_sensMatrix(dR)
	call deall_solnVectorMTX(H)
	call deall_modelParam(dmisfit)
	call deall_modelParam(dmisfitSmooth)

	call date_and_time(values=tarray)
	etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	rtime = etime - stime

  end subroutine calc_derivative  ! calc_derivative


  ! ***************************************************************************
  ! * calc_jacobian is a subroutine to calculate specified responses and the
  ! * full Jacobian of these responses with respect to the original model
  ! * parameters, for each frequency in turn. Designed to use for the detailed
  ! * sensitivity analysis at a chosen point in the model parametrization domain,
  ! * such as at a local minimum.

  subroutine calc_jacobian(rtime)

	use data_vectors
	use global
	use jacobian
	use boundaries
	use output
	use initFields
	use dataFunc
	use dataMisfit
    use transmitters
    use dataTypes
    use receivers
	implicit none

	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j
    type (cvector)							:: H,B,Hconj,B_tilde,dH,dE,Econj,F
	type (sparsevecc)						:: bv_H,bv_dH,Hb
	type (rvector)							:: dE_real,dE_imag
	type (sparsevecc)						:: g_sparse
	type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	type (modelParam_t)						:: rsens,isens,rsensTemp,isensTemp
	integer									:: ifreq,ifunc,iobs,ilayer
	logical									:: adjoint,delta
	real(8),dimension(ncoeff)				:: da

	! Start the (portable) clock
	call date_and_time(values=tarray)
    stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	!write(0,*) param%p(:)%value

	call initialize_fields(H,B)

	do ifreq=1,freqList%n

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  write(6,*) 'Solving 3D forward problem for freq ',ifreq,freq%value

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  call operatorM(H,B,omega,rho,grid,fwdCtrls,errflag,adjoint)

	  ! compute and output fields; C and D responses at cells
	  call outputSolution(freq,H,slices,grid,cUserDef,rho,'h')
	  !call outputH('jacobian',H,grid)
	  !call outputCellResp('jacobian',H)

	  ! compute and output C and D responses at observatories
	  call calcResponses(freq,H,dat,psi)
	  call outputResponses(freq,psi,freqList,TFList,obsList,outFiles,dat)
	  ! compute C and D residuals
	  call calcResiduals(freq,dat,psi,res)

	  ! if computing different kinds of misfits, be consistent;
	  ! write different kinds into different data structures
	  call calcMisfit(freq,res,misfit,misfit%name)


	  do ifunc=1,nfunc

		do iobs=1,nobs

		  if (.not.obsList%info(iobs)%defined) then
			sens%da_real(ifreq,ifunc,iobs,:) = 0.0d0
			sens%da_imag(ifreq,ifunc,iobs,:) = 0.0d0
			cycle
		  end if

		  ! compute $\g_j$ and set $\tilde{\b} = \g_j$
		  call Lrows(TFList%info(ifunc),obsList%info(iobs),H,g_sparse)
		  call create_cvector(grid,F,EDGE)
		  call add_scvector(C_ONE,g_sparse,F)

		  ! $M*^{-1}_{\rho,-\omega} ( g_j )$
		  adjoint = .TRUE.
		  delta = .TRUE.
		  call create_cvector(grid,dH,EDGE)
		  ! Forcing term F should not contain any non-zero boundary values
		  call operatorM(dH,F,omega,rho,grid,fwdCtrls,errflag,adjoint,delta)
		  !call outputH('jacobian_dh',dH,grid)
		  !call outputCellResp('jacobian_dh',dH)
		  ! call outputSolution(freq,dH,slices,grid,cUserDef,rho,'dh')

		  ! Pre-divide the interior components of dH by elementary areas
		  call operatorD_Si_divide(dH,grid)

		  ! $C D_{S_i}^{-1} M*^{-1}_{\rho,\omega} ( g_j)$
		  call createBC(Hb,grid)
		  call insertBC(Hb,dH)
		  call operatorC(dH,dE,grid)

		  ! $\bar{\e} = C \bar{\h}$
		  Hconj = conjg(H)
		  call operatorD_l_mult(Hconj,grid)
		  call operatorC(Hconj,Econj,grid)

		  ! $D_{\bar{\e}} C D_{S_i}^{-1} M*^{-1}_{\rho,\omega} ( g_j )$
		  call diagMult(Econj,dE,dE) !dE = Econj * dE

		  ! $L^T D_{\bar{\e}} ... M*^{-1}_{\rho,\omega} g_j$
		  dE_real = real(dE)
		  call operatorLt(sens%drho_real%v,dE_real,grid)

		  dE_imag = imag(dE)
		  call operatorLt(sens%drho_imag%v,dE_imag,grid)

		  ! $P^T L^T D_{\bar{\e}} C A^{-1}_{\rho,-\omega} g_j$
		  call operatorPt(sens%drho_real,rsensTemp)
		  call operatorPt(sens%drho_imag,isensTemp)

		  rsens = multBy_CmSqrt(rsensTemp)
		  isens = multBy_CmSqrt(isensTemp)

		  call getParamValues_modelParam(rsens,sens%da_real(ifreq,ifunc,iobs,:))
		  call getParamValues_modelParam(isens,sens%da_imag(ifreq,ifunc,iobs,:))


		  ! account for the fact that this is $\pd{\bar{\psi}^j_\omega}{\a}$
		  sens%da_imag = - sens%da_imag

		  sens%da(ifreq,ifunc,iobs,:) = - dcmplx(sens%da_real(ifreq,ifunc,iobs,:),sens%da_imag(ifreq,ifunc,iobs,:))

		end do

                call outputJacobian(freq,psi,sens,freqList,TFList,obsList,param,cUserDef)
		  !print *, 'ifreq,ifunc,iobs,dpsi/da: ',ifreq,ifunc,iobs,&
		!		  sens%da(ifreq,ifunc,iobs,:)

	  end do

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  rtime = rtime + ftime

	end do

	call misfitSumUp(res,misfit,misfitValue)

	!call outputMisfit(param,misfit,misfitValue,cUserDef)

	continue

	call deall_modelParam(rsensTemp)
	call deall_modelParam(rsens)
	call deall_modelParam(isensTemp)
	call deall_modelParam(isens)
	call deall_sparsevecc(bv_H)
	call deall_sparsevecc(bv_dH)
	call deall_sparsevecc(Hb)
	call deall_sparsevecc(g_sparse)

  end subroutine calc_jacobian	! calc_jacobian


  ! ***************************************************************************
  ! * calc_symmetry is a subroutine to run a symmetry test for reciprocity of
  ! * the solution by computing the derivative the direct way, rather than adjoint;
  ! * this method may take a long time to complete if many parameters are used.

  subroutine calc_symmetry(rtime)

	use data_vectors
	use global
	use jacobian
	use boundaries
	use output
	use initFields
	use dataFunc
	use dataMisfit
    use transmitters
    use dataTypes
	implicit none

	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
    type (cvector)							:: H,B,F,Hconj,B_tilde,dH,dE,Econj,Bzero,dR
	type (rvector)							:: dE_real
	type (rscalar)							:: drho
	type (sparsevecc)						:: Hb
	type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	type (modelCoeff_t)						:: coeff
	integer									:: ifreq,ifunc,index
	logical									:: adjoint,delta
	character(1)							:: cfunc
	character(80)							:: fn_err

	! Start the (portable) clock
	call date_and_time(values=tarray)
    stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	call initialize_fields(H,B)

	call create_rscalar(grid,drho,CENTER)

	do ifreq=1,freqList%n

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  write(6,*) 'Solving 3D forward problem for freq ',ifreq,freq%value

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  call operatorM(H,B,omega,rho,grid,fwdCtrls,errflag,adjoint)
	  !call create_cvector(grid,H,EDGE)
	  !H%x = C_ONE
	  !H%y = C_ONE
	  !H%z = C_ONE

	  ! compute and output fields & C and D responses at cells
	  call outputSolution(freq,H,slices,grid,cUserDef,rho,'h')

	  ! compute and output C and D responses at observatories
	  call calcResponses(freq,H,dat,psi)
	  call outputResponses(freq,psi,freqList,TFList,obsList,outFiles,dat)

	  ! compute and output C and D residuals
	  call calcResiduals(freq,dat,psi,res)
	  !call outputResiduals(freq,res,TFList,obsList,outFiles)

	  ! if computing different kinds of misfits, be consistent;
	  ! write different kinds into different data structures
	  call calcMisfit(freq,res,misfit,misfit%name)

	  ! compute the weighted residuals and start the derivative computations
	  call calcResiduals(freq,dat,psi,wres,weighted=.TRUE.)

	  do ifunc=1,nfunc

		print *
		print *, 'Starting the derivative computations for ',&
				  trim(TFList%info(ifunc)%name), ' responses...'


		!---------------------------------------------------------
		! Direct computations

		do index = 1,param%nc

		  coeff = getCoeff_modelParam(param,index)

		  if (coeff%frozen) then
			cycle
		  end if

		  drho = vectorPj(param,index)

		  call operatorL(drho%v,dE_real,grid)
		  !dE = dE_real

		  ! ${\e} = C {\h}$
		  call create_cvector(grid,Hconj,EDGE)
		  Hconj = H
		  call operatorD_l_mult(Hconj,grid)
		  call operatorC(Hconj,Econj,grid)

		  call diagMult(Econj,dE_real,dE) !dE = Econj * dE_real

		  call operatorCt(dE,F,grid)
		  call createBC(Hb,grid)
		  call insertBC(Hb,F)

		  call operatorD_Si_divide(F,grid)

		  ! $M*^{-1}_{\rho,-\omega} ( G_\omega r_\omega )$
		  adjoint = .FALSE.
		  delta = .TRUE.
		  ! Forcing term F should not contain any non-zero boundary values
		  call operatorM(dH,F,omega,rho,grid,fwdCtrls,errflag,adjoint,delta)

		  call Lmult(dH,H,psi%v(ifreq,ifunc,:))

		  misfit%dRda(ifreq,ifunc,index) &
			= -2. * dreal(sum(conjg(wres%v(ifreq,ifunc,:)%resp%value) * psi%v(ifreq,ifunc,:)%resp%value))

		end do
		! End of direct computations
		!---------------------------------------------------------

	  end do

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  rtime = rtime + ftime

	end do

	! Output the summary and exit
	call misfitSumUp(res,misfit,misfitValue,dmisfitValue)

	!call outputMisfit(param,misfit,misfitValue,cUserDef)
	!call outputDerivative(param,misfit,dmisfitValue,cUserDef)

	continue

	call deall_rscalar(drho)
	call deall_sparsevecc(Hb)

  end subroutine calc_symmetry  ! calc_symmetry


  ! ***************************************************************************
  ! * calc_symmetric_operators is a subroutine to run a symmetry test for each
  ! * of the operators involved in the calculation

  subroutine calc_symmetric_operators(rtime)

	use data_vectors
	use global
	use jacobian
	use output
	use modelmap
	use initFields
	use dataFunc
	use dataMisfit
	use sg_spherical
    use transmitters
    use dataTypes
	implicit none

	real, intent(inout)						:: rtime  ! run time
	real									:: ftime  ! run time per frequency
	real									:: stime, etime ! start and end times
    integer, dimension(8)					:: tarray ! utility variable
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,kk
    type (cvector)							:: H,H1,H2,H3,B,F,F1,F2,F3,F1conj,F2conj
	type (rvector)							:: r1,r2
	type (cvector)							:: e1,e2
	type (rscalar)							:: drho,drho1,drho2
	type (sparsevecc)						:: g_sparse
	type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	type (modelParam_t)						:: dparam
  	logical									:: adjoint,delta
	integer									:: ifreq,ifunc
	integer, dimension(3)					:: k
	!real(8), dimension(grid%nx+1,grid%ny,grid%nz)		:: random_x
	!real(8), dimension(grid%nx,grid%ny+1,grid%nz)		:: random_y
	!real(8), dimension(grid%nx,grid%ny,grid%nz+1)		:: random_z
	real(8), dimension(grid%nx,grid%ny+1,grid%nz+1)		:: random_x
	real(8), dimension(grid%nx+1,grid%ny,grid%nz+1)		:: random_y
	real(8), dimension(grid%nx+1,grid%ny+1,grid%nz)		:: random_z
	complex(8), dimension(3,3)				:: value
	real(8)									:: value0,value1,value2
	logical, dimension(:,:,:), allocatable	:: mask
	!logical									:: verbose

	! Start the (portable) clock
	call date_and_time(values=tarray)
    stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	allocate(mask(grid%nx,grid%ny,grid%nz))
	mask = .TRUE.

	! Symmetry test for operators P (initModel) and Pt
	print *, 'Symmetry test for operators P (initModel) and Pt'
	call create_rscalar(grid,drho1,CENTER)
	call create_rscalar(grid,drho2,CENTER)
	drho1%v = ONE
	call operatorP(param,drho2)
	value0 = dot_product(pack(drho1%v,mask),pack(drho2%v,mask))
	value1 = sum(drho1%v(:,:,grid%nzAir+1:nz-1) * drho2%v(:,:,grid%nzAir+1:nz-1))
	print *, "Value 1 = ",value1
	dparam = param
	call operatorPt(drho1,dparam)
	!value2 = dot_product(param%p%value,dparam%p%value)
	value2 = dotProd_modelParam_f(param,dparam)
	print *, "Value 2 = ",value2
	!stop
	do i= grid%nx,grid%nx
	  do j =grid%ny,grid%ny
		do kk =1, grid%nz
		  if (drho1%v(i,j,kk)*drho2%v(i,j,kk) /= 0.0) then
			!print *, 'rho:',i,j,kk,rho(i,j,kk),drho2%v(i,j,kk)
		  end if
		end do
	  end do
	end do
	!print *,param%p%value
	!print *,dparam%p%value


	! Symmetry test for operators C and Ct
	print *, 'Symmetry test for operators C and Ct'
	call create_cvector(grid,f1,FACE)
	call random_number(random_x)
	call random_number(random_y)
	call random_number(random_z)
	f1%x = C_ONE
	f1%y = C_ONE
	f1%z = C_ONE
	!f1%z(:,:,1) = random_z(:,:,1)
	!f1%z(:,:,4) = random_z(:,:,4)
	!f1%z(:,:,nz+1) = random_z(:,:,nz+1)
	!call validate_cvector(f1,verbose)
	call create_cvector(grid,e1,EDGE)
	e1%x = random_x
	e1%y = random_y
	e1%z = random_z
	!e1%z(1,1,9) = random_z(1,1,9)
	!call validate_cvector(e1)
	call operatorC(e1,f2,grid)
	!call validate_cvector(f2,verbose)
	!f2%x(grid%nx+1,:,:) = C_ZERO
	!value1 = dotProd_noConj(f1,f2)
	value1 = dotProd(f1,f2)
	print *, "Value 1 = ",value1
	call operatorCt(f1,e2,grid)
	!call validate_cvector(e2,verbose)
	!e2%y(grid%nx+1,:,:) = C_ZERO
	!e2%z(grid%nx+1,:,:) = C_ZERO
	!e2%z(2:grid%nx,1,:) = C_ZERO
	!e2%z(2:grid%nx,ny+1,:) = C_ZERO
	!value2 = dotProd_noConj(e1,e2)
	value2 = dotProd(e1,e2)
	print *, "Value 2 = ",value2


	! Symmetry test for operators G and Gt
	print *, 'Symmetry test for operators G and Gt'
	call create_cvector(grid,H,EDGE)
	H%x = (2.0,1.0)
	H%y = (0.5,3.0)
	H%z = C_ONE
	freq = freqList%info(1)
	call calcResponses(freq,H,dat,psi)
	!dat%v(1,1,1)%resp%value = (1000000.,-200000)
	!dat%v(1,1,2)%resp%value = (1500000.,-250000)
  	call LmultT(dat%v(1,1,:),H,E2)
	!print *,dat%v(1,1,:)%resp%value
	value1 = dotProd(e1,e2)
	print *, "Value 1 = ",value1
  	call Lmult(E1,H,psi%v(1,1,:))
	!print *,psi%v(1,1,:)%resp%value
	value2 = sum(conjg(dat%v(1,1,:)%resp%value) * psi%v(1,1,:)%resp%value)
	print *, "Value 2 = ",value2


	!stop

	! Symmetry test for operators L and Lt
	print *, 'Symmetry test for operators L and Lt'
	call create_rvector(grid,r1,FACE)
	r1%x = 2.0
	r1%y = 3.5
	r1%z = 5.0
	!e1%z(2,4,6) = 5.0
	!e1%x(1,:,:) = R_ZERO
	!e1%x(grid%nx,:,:) = R_ZERO
	drho1%v = ONE*2.0 !rho
	call operatorL(drho1%v,r2,grid)
	!call operatorL(rho,e2,grid)
	r2%x(grid%nx+1,:,:) = R_ZERO
	!e1%x(grid%nx+1,:,:) = R_ZERO
	value1 = dotProd(r1,r2)
	print *, "Value 1 = ",value1
	call operatorLt(drho2%v,r1,grid)
	value0 = 0.0d0
	do i= 1,grid%nx
	  do j =1,grid%ny
		do kk =1, grid%nz
		  !value0 = value0 + drho1%v(i,j,kk)*drho2%v(i,j,kk)
		  if (drho1%v(i,j,kk)*drho2%v(i,j,kk) > 0.0) then
			!print *, 'rho:',i,j,kk, drho2%v(i,j,kk)
		  end if
		  if (r1%z(i,j,kk)*r2%z(i,j,kk) > 0.0) then
			!print *, 'e: ',i,j,kk, r2%z(i,j,kk)
		  end if
		end do
	  end do
	end do
	value0 = dot_product(pack(drho1%v,mask),pack(drho2%v,mask))
	!value2 = sum(drho%v(:,:,grid%nzAir+1:nz) * drho2%v(:,:,grid%nzAir+1:nz))
	print *, "Value 2 = ",value0


	!stop

	! Symmetry test for operators M  and M*

	print *, 'Symmetry test for operators M  and M*'

	call initialize_fields(H,B)


	k(1) = 8
	k(2) = 3
	k(3) = 11

	do i=1,3
	  if (.not.obsList%info(k(i))%defined) then
		write(0,*) 'Error: observatory ',k(i),' not defined'
		stop
	  end if
	end do

	ifreq=1
	ifunc=1

	  omega  = 2.0d0*pi*freqList%info(ifreq)%value     ! angular frequency (radians/sec)

	  write(6,*) 'Solving 3D forward problem for freq ',ifreq,freqList%info(ifreq)%value

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  call operatorM(H,B,omega,rho,grid,fwdCtrls,errflag,adjoint)

	  print *, 'Starting the symmetry test for ',k(1), ' and ',k(2)


	  ! compute $\g_j$ and set $\tilde{\b} = \g_j$
	  call Lrows(TFList%info(ifunc),obsList%info(k(1)),H,g_sparse)
	  call create_cvector(grid,F1,EDGE)
	  call add_scvector(C_ONE,g_sparse,F1)

	  ! compute $\g_j$ and set $\tilde{\b} = \g_j$
	  call Lrows(TFList%info(ifunc),obsList%info(k(2)),H,g_sparse)
	  call create_cvector(grid,F2,EDGE)
	  call add_scvector(C_ONE,g_sparse,F2)

	  write(6,*) 'Solving 3D forward problem for obs.',1,' index ',k(1)

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  call operatorM(H1,F1,omega,rho,grid,fwdCtrls,errflag,adjoint)

	  value(1,2) = dotProd(F2,H1)
	  print *, 'Solution: 1, ',value(1,2)

	  write(6,*) 'Solving 3D adjoint problem for obs.',2,' index ',k(2)

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.TRUE.
	  call operatorM(H2,F2,omega,rho,grid,fwdCtrls,errflag,adjoint)

	  value(2,1) = dotProd(F1,H2)
	  print *, 'Solution: 2, ',value(2,1)


	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  rtime = rtime + ftime

	  print *, 'The fwd and adj results should be complex conjugates of each other...'

	  print *, 'Solution g_2^* Mfwd^{-1} g_1: ',value(1,2)
	  print *, 'Solution g_1^* Madj^{-1} g_2: ',value(2,1)

	continue

	call deall_modelParam(dparam)
	call deall_rscalar(drho)
	call deall_sparsevecc(g_sparse)

  end subroutine calc_symmetric_operators  ! calc_symmetric_operators
