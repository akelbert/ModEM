! *****************************************************************************
module main
	! These subroutines are called from the main program earth3d only

  use input
  use iospec
  use iotypes
  use interp
  use math_constants
  use model_operators
  use modeldef
  use modelmap
  use global
  use grid_orig
  use field_vectors
  use data_vectors
  use sg_scalar
  implicit none


Contains

  ! ***************************************************************************
  ! * InitGlobalData is the routine to call to initialize all derived data types
  ! * and other variables defined in modules basics, modeldef, datadef and
  ! * in this module. These include:
  ! * 1) constants
  ! * 2) grid information: nx,ny,nz,x,y,z, cell centre coordinates
  ! * 3) model information: layers, coefficients, rho on the grid
  ! * 4) periods/freq info: nfreq, freq
  ! * 5) forward solver controls
  ! * 6) output file names
  ! * The information about the grid, frequency and resistivity map is saved
  ! * in the module global, and the grid is also saved in module griddef.
  ! * Frequencies are listed in ascending order (by value)
  ! * Observatories are listed in ascending order (alphabetic)

  subroutine InitGlobalData(fn_startup,da)

	implicit none
    character(80), intent(inout)				:: fn_startup
	integer										:: i,ios=0,istat=0
	character(100)								:: label
	real(8),dimension(:),intent(in),optional	:: da
	type (param_info)							:: p_delta
	type (coeff_info)							:: coeff

	!--------------------------------------------------------------------------
	! Initialize the math/physics constants defined in module basics
	!call init_const(pi,d2r,r2d,mu0,rair,Rearth)
	!--------------------------------------------------------------------------
	! Initialize the user-defined switches and file names
    call readStartFile(fn_startup,cUserDef)
	!--------------------------------------------------------------------------
	! Initialize preconditioning information
	call initFunctional(cUserDef,misfitType)
	!--------------------------------------------------------------------------
	! Read and compute grid information
	call initGrid(cUserDef,grid)
	!--------------------------------------------------------------------------
	! Read and save thin shell conductance information, if present
	call initCrust(cUserDef,grid,crust)
	!--------------------------------------------------------------------------
	! Read main parametrization
	call initParamY(cUserDef,p_input)
	!--------------------------------------------------------------------------
	! Check whether base parametrization exists
	inquire(FILE=cUserDef%fn_param0,EXIST=exists)
	!--------------------------------------------------------------------------
	! If exists, read it; else create a skeleton from main parametrization
	if (exists) then
	  call initParamY(cUserDef,p0_input,p0=.TRUE.)
		! NOTE: regularization information is taken from the prior model!!!
		if (.not.VerifyParamY(p_input,p0_input)) then
			write(0,*) 'Warning: Using the layer structure from the prior model'
			p_input%L = p0_input%L
			if (.not.VerifyParamY(p_input,p0_input)) then
				write(0,*) 'Warning: Base parametrization incompatible with main model'
				stop
			end if
		end if
	else
	  write(0,*) 'Warning: No base parametrization found; zero model will be used'
	  p0_input = CreateEmptyParamY(p_input)
	end if
	!--------------------------------------------------------------------------
	! Compute the correction (only needed if run for a test perturbation)
	if (present(da)) then
	  p_delta = CreateEmptyParamY(p_input)
	  call FillParamValuesY(p_delta,da)
	  p_input = p_input + p_delta
	  call DeleteParamY(p_delta)
	  !param%p(:)%value = param%p(:)%value + da(:)
	end if
	!--------------------------------------------------------------------------
	! 'Smooth' the parametrization by applying inverted regularization operator
	param = SmoothSqrtParamY(p_input)
	!--------------------------------------------------------------------------
	! Compute parametrization to use (for model norm we will still use p_input)
	param = param + p0_input
	!--------------------------------------------------------------------------
	! Allocate the resistivity vector
	allocate(rho(grid%nx,grid%ny,grid%nz))
	!--------------------------------------------------------------------------
	! Initialize model information in air and possibly crust, if given
	call initResist(grid,crust,rho)
	!--------------------------------------------------------------------------
	! Compute model information everywhere else in the domain
	call initModel(grid,param,rho)
	!--------------------------------------------------------------------------
	! Read the information about the frequencies
	call initFreq(cUserDef,freqList)
	!--------------------------------------------------------------------------
	! Initialize observatory locations required for output
	call initCoords(cUserDef,obsList)
	!--------------------------------------------------------------------------
	! Compute the observatory locations on the grid and interpolation weights
	call initObsList(grid,obsList)
	!--------------------------------------------------------------------------
	! Initialize transfer function information
	call initTF(cUserDef,TFList)
	!--------------------------------------------------------------------------
	! Read radii at which we will output the full solution
	call initSlices(cUserDef,slices)
	!--------------------------------------------------------------------------
	! Read forward solver controls to store in fwdCtrls
	call initControls(cUserDef,fwdCtrls)
	!--------------------------------------------------------------------------
	! Initialize the file names to store in outFiles
    call initOutput(cUserDef,outFiles)
	!--------------------------------------------------------------------------
	nfreq=freqList%n
	nobs =obsList%n
	nfunc=TFList%n
	ncoeff=param%nc

	select case (cUserDef%verbose)
	case ('debug')
	  print *,'Output all information including debugging lines.'
	  output_level = 5
	case ('full')
	  print *,'Output full information to screen and to files.'
	  output_level = 4
	case ('regular')
	  print *,'Output information to files, and progress report to screen (default).'
	  output_level = 3
	case ('compact')
	  print *,'Output information to files, and compact summary to screen.'
	  output_level = 2
	case ('result')
	  print *,'Output information to files, and result to screen.'
	  output_level = 1
	case ('none')
	  print *,'Output nothing at all except result to screen and to files.'
	  output_level = 0
	case default
	  output_level = 3
	end select

	!--------------------------------------------------------------------------
	! Helpful output
	call PrintParamInfoY(p0_input,output_level-1,"Prior model m_0 = ")
	call PrintParamInfoY(p_input,output_level-1,"Input model \hat{m} = ")
	call PrintParamInfoY(param,output_level,"Final model m = ")
	!--------------------------------------------------------------------------

	! If this information is required, initialize data functionals
	if (cUserDef%calculate == 'original') then
	else
	end if

	! If this information is required, initialize data functionals
	if (cUserDef%calculate == 'original') then
	else
	  allocate(dat%v(nfreq,nfunc,nobs),dat%n(nfreq,nfunc),STAT=istat)
	  allocate(psi%v(nfreq,nfunc,nobs),psi%n(nfreq,nfunc),STAT=istat)
	  allocate(res%v(nfreq,nfunc,nobs),res%n(nfreq,nfunc),STAT=istat)
	  allocate(ndat(nfreq,nfunc),STAT=istat)
	  allocate(misfitValue(nfunc))
	  call initData(cUserDef,dat,obsList,freqList,TFList)
	  call initMisfit(misfitType,TFList,freqList,dat,misfit)
	  if (cUserDef%calculate /= 'responses') then
		allocate(dmisfitValue(nfunc,ncoeff))
		allocate(wres%v(nfreq,nfunc,nobs),wres%n(nfreq,nfunc),STAT=istat)  ! weighted residuals
		allocate(misfit%dRda(nfreq,nfunc,ncoeff),STAT=istat)
		call create_rscalar(grid,sens%drho_real,CENTER)
		call create_rscalar(grid,sens%drho_imag,CENTER)
		allocate(sens%da_real(nfreq,nfunc,nobs,ncoeff),STAT=istat)
		allocate(sens%da_imag(nfreq,nfunc,nobs,ncoeff),STAT=istat)
		allocate(sens%da(nfreq,nfunc,nobs,ncoeff),STAT=istat)
	  end if
	end if


    ! Output the frequencies
    write(0,*) 'nfreq =',nfreq
    do i=1,nfreq
	  write(0,*)' freq(',trim(freqList%info(i)%code),')=',freqList%info(i)%value
	end do

	! Too complicated to rewrite input to all subroutines that use x,y,z,nx,ny,nz
	! in terms of the grid variable, but that would be the way to do it in the
	! future. Currently, use this patch instead. x,y,z stored in module griddef
	nx = grid%nx; ny = grid%ny; nz = grid%nz
	nzEarth = grid%nzEarth; nzAir = grid%nzAir
	allocate(x(nx),y(ny+1),z(nz+1))
	x = grid%x; y = grid%y; z = grid%z

	return

  end subroutine InitGlobalData	! InitGlobalData


  ! ***************************************************************************
  ! * DeleteGlobalData deallocates all allocatable data defined globally.
  subroutine DeleteGlobalData()

	integer	:: istat

	! Deallocate global variables that have been allocated by InitGlobalData()
	call DeleteParamY(param)
	call DeleteParamY(p_input)
	call DeleteParamY(p0_input)
	call DeleteParamY(p_smooth)
	call DeleteParamY(p_diff)
	deallocate(x,y,z,STAT=istat)
	deallocate(grid%x,grid%y,grid%z,STAT=istat)
	deallocate(grid%ph,grid%th,grid%r,STAT=istat)
	deallocate(freqList%info,STAT=istat)
	deallocate(rho,STAT=istat)
	deallocate(misfit%value,STAT=istat)
	deallocate(misfit%ndat,STAT=istat)
	deallocate(misfit%weight,STAT=istat)
	deallocate(misfit%dRda,STAT=istat)
	deallocate(sens%da_real,STAT=istat)
	deallocate(sens%da_imag,STAT=istat)
	deallocate(sens%da,STAT=istat)
	call deall_rscalar(sens%drho_real)
	call deall_rscalar(sens%drho_imag)
	deallocate(dat%v,dat%n,STAT=istat)
	deallocate(psi%v,psi%n,STAT=istat)
	deallocate(res%v,res%n,STAT=istat)
	deallocate(wres%v,wres%n,STAT=istat)

	if (allocated(misfitValue)) then
	  deallocate(misfitValue)
	end if
	if (allocated(dmisfitValue)) then
	  deallocate(dmisfitValue)
	end if

  end subroutine DeleteGlobalData	! DeleteGlobalData


  ! ***************************************************************************
  ! * InitDivergenceCorrection allocates and initializes with zeros all global
  ! * field vectors that are required to keep the divergence correction valid
  subroutine initDivergenceCorrection()

	!use field_vectors
	use dimensions

	integer					  :: istat

	allocate(hx(np1),hy(np1),hz(np1),STAT=istat)
	hx = (0.0d0,0.0d0)
	hy = (0.0d0,0.0d0)
	hz = (0.0d0,0.0d0)
	allocate(sx(np1),sy(np1),sz(np1),STAT=istat)
	sx = (0.0d0,0.0d0)
	sy = (0.0d0,0.0d0)
	sz = (0.0d0,0.0d0)
	allocate(divr(np1),divi(np1),STAT=istat)
	divr = R_ZERO
	divi = R_ZERO
	allocate(divsr(np1),divsi(np1),STAT=istat)
	divsr = R_ZERO
	divsi = R_ZERO

  end subroutine initDivergenceCorrection


  ! ***************************************************************************
  ! * InitGlobalArrays allocates and initializes with zeros all global arrays
  ! * with an allocatable attribute defined in various modules. No allocatable
  ! * arrays should be defined in the main program. Use modules instead.
  subroutine InitGlobalArrays()

	use dimensions
	use boundaries
	use ringcurrent
	!use initFields
	integer	:: istat

	np1=   nx*(ny+1)*(nz+1)	! number of scalars on grid (with some liberal extras)
	np2=3 *nx*(ny+1)*(nz+1)	! number of vectors on grid (with some liberal extras)
	np3=   nx*(ny-1)+2		! number of surface nodes (exact)
  	np4=9 *nx*(ny+1)*(nz+1)
	np5=21*nx*(ny+1)*(nz+1)
	np6=4 *nx*(ny+1)*(nz+1)

	bv1= 2*(ny-1)*(nx-1) + (ny-1)
	bv2= 2*nx*(ny-1) + nx
	bv3= nx*(ny-1)+2

  ! Dynamic allocation of all global vectors
!	allocate(vectorh(np2),vectorb(np2),STAT=istat)

  ! Initialize all global vectors here
!	vectorh = C_ZERO
!	vectorb = C_ZERO

	call initDivergenceCorrection()

  end subroutine InitGlobalArrays ! InitGlobalArrays


  ! ***************************************************************************
  ! * DeleteGlobalArrays does exactly what it is meant to do: deallocates
  ! * all allocatable arrays. No checks are currently performed, since currently
  ! * all memory is deallocated in this subroutine. If memory is ever deallocated
  ! * otherwise, if(allocated(var)) checks are required.
  subroutine DeleteGlobalArrays()

	integer	:: istat

	! Deallocate fields and temporary variables
	deallocate(hx,hy,hz,STAT=istat)
	deallocate(sx,sy,sz,STAT=istat)
	deallocate(divr,divi,STAT=istat)
	deallocate(divsr,divsi,STAT=istat)
!	deallocate(vectorh,vectorb,STAT=istat)

  end subroutine DeleteGlobalArrays	! DeleteGlobalArrays


end module main
