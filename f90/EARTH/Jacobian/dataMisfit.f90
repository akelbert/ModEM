! *****************************************************************************
module dataMisfit
  ! Module containing the definitions and subroutines related to the penalty
  ! functional only. Only use it if we want to compute the data misfit.

  use global
  use responses
  use functionals
  use model_operators
  use dataspace
  implicit none


Contains

  ! ***************************************************************************
  ! * calcResponses is the subroutine to compute the values of data functionals

  subroutine calcResponses(freq,H,dat,psi)

	!uses: TFList,obsList

	type (transmitter_t), intent(in)		:: freq
	type (dataVecMTX_t), intent(in)			:: dat
	type (dataVecMTX_t), intent(inout)		:: psi
	type (cvector), intent(in)				:: H
	integer									:: i,j,k

    psi%nTx = dat%nTx
	psi%n = dat%n

	i = freq%i

	do j=1,TFList%n
	  select case ( TFList%info(j)%name )
	  case ('C')
		do k=1,nobs
		  psi%v(i,j,k)%resp%value = dataResp(C_ratio,obsList%info(k),H)
		  psi%v(i,j,k)%resp%err = 0.0d0
		  psi%v(i,j,k)%resp%exists = dat%v(i,j,k)%resp%exists
		  psi%v(i,j,k)%func => TFList%info(j)
		  psi%v(i,j,k)%obs => obsList%info(k)
		  psi%v(i,j,k)%freq => freqList%info(i)
		end do
	  case ('D')
		do k=1,nobs
		  psi%v(i,j,k)%resp%value = dataResp(D_ratio,obsList%info(k),H)
		  psi%v(i,j,k)%resp%err = 0.0d0
		  psi%v(i,j,k)%resp%exists = dat%v(i,j,k)%resp%exists
		  psi%v(i,j,k)%func => TFList%info(j)
		  psi%v(i,j,k)%obs => obsList%info(k)
		  psi%v(i,j,k)%freq => freqList%info(i)
		end do
	  case default
		write(0,*) 'Warning: no transfer functions specified to compute'
		stop
	  end select
	  do k=1,obsList%n
		if (.not.obsList%info(k)%defined) then
		  psi%v(i,j,k)%resp%exists = .FALSE.
		  psi%n(i,j) = psi%n(i,j) - 1
		end if
	  end do
	end do

  end subroutine calcResponses  ! calcResponses


  ! ***************************************************************************
  ! * calcResiduals is the subroutine to compute the data residuals from the
  ! * data and the corresponding responses obtained by first computing the fields.

  subroutine calcResiduals(freq,dat,psi,res,weighted)

	type (transmitter_t), intent(in)			:: freq
	type (dataVecMTX_t), intent(in)				:: dat,psi
	type (dataVecMTX_t), intent(inout)			:: res
	logical, intent(in), optional				:: weighted
	integer										:: i,j,k
	complex(8)									:: value
	real(8)										:: error

	i = freq%i

	! Cycle over functional types and the observatories for this frequency
	do j=1,TFList%n
	  do k=1,obsList%n

  		res%v(i,j,k)%func => dat%v(i,j,k)%func
  		res%v(i,j,k)%obs => dat%v(i,j,k)%obs
  		res%v(i,j,k)%freq => dat%v(i,j,k)%freq
		! The residual exists if both the data entry and the response make sense
		if (dat%v(i,j,k)%resp%exists.and.psi%v(i,j,k)%resp%exists) then
		  res%v(i,j,k)%resp%exists = .TRUE.
		  res%v(i,j,k)%resp%value = psi%v(i,j,k)%resp%value - dat%v(i,j,k)%resp%value
		  res%v(i,j,k)%resp%err = dat%v(i,j,k)%resp%err
		else
		  res%v(i,j,k)%resp%exists = .FALSE.
		  res%v(i,j,k)%resp%value = C_ZERO
		  res%v(i,j,k)%resp%err = LARGE_REAL
		end if

	  end do
	  ! Count the number of existing residuals
	  res%n(i,j) = count(res%v(i,j,:)%resp%exists)
	end do

	! If instead weighted responses are required, divide by error squared
	if (present(weighted)) then
	  if (weighted) then
		do j=1,TFList%n
		  do k=1,obsList%n
			value = res%v(i,j,k)%resp%value
			error = res%v(i,j,k)%resp%err
			res%v(i,j,k)%resp%value = value/(error**2)
			res%v(i,j,k)%resp%err = ONE
		  end do
		end do
	  end if
	end if

  end subroutine calcResiduals  ! calcResiduals


  ! ***************************************************************************
  ! * calcMisfit performs the misfit computations per frequency

  subroutine calcMisfit(freq,res,misfit,name)

	type (transmitter_t), intent(in)				:: freq
	type (misfit_t), intent(inout)			:: misfit
	type (dataVecMTX_t), intent(in)				:: res
	character(80), intent(in)					:: name
	real(8)										:: rval,ival,error
	integer										:: i,j,k

	i = freq%i

	!call OutputResiduals(freq,res,TFList,obsList)
	do j=1,TFList%n
	  do k=1,obsList%n
		if (.not.res%v(i,j,k)%resp%exists) then
		  cycle
		end if

		rval = dreal(res%v(i,j,k)%resp%value)
		ival = dimag(res%v(i,j,k)%resp%value)
		error = res%v(i,j,k)%resp%err

		misfit%ndat(i,j) = res%n(i,j)

		! Calculate the misfit
		select case (trim(name))
		case ('Mean Squared')
		  misfit%value(i,j) = &
			 misfit%value(i,j) + (rval/error)**2 + (ival/error)**2
		case default
		  write(0,*) 'Warning: (compute_misfit) unknown penalty functional ',trim(name)
		  return
		end select
	  end do
	  !write(0,'(a40,2i3,i6,g15.7)') 'ifreq,ifunc,ndat,misfit = ',i,j,&
	  !				misfit%ndat(i,j),misfit%value(i,j)/(2*misfit%ndat(i,j))
	end do

  end subroutine calcMisfit	! calcMisfit


  ! ***************************************************************************
  ! * misfitSumUp

  subroutine misfitSumUp(res,misfit,misfitValue,dmisfitValue)

	! this is temporary: output to files will be moved to output.f90
	use iotypes
	use iospec
	! uses: TFList,p_input,param
	type (misfit_t), intent(inout)			  :: misfit
	type (dataVecMTX_t), intent(in)				  :: res
	type (modelParam_t)							  :: dparam,weighted_norm
  	real(8), dimension(:), intent(out)			  :: misfitValue  !nfunc
	real(8), dimension(:,:), intent(out),optional :: dmisfitValue !nfunc,ncoeff
	type (modelCoeff_t), dimension(param%nc)	  :: norm,grad,point
	character(80)                    :: comment
	real(8)										  :: v1,v2
	integer										  :: j,l

	misfitValue = 0.0d0
	write(0,*) 'Penalty functional used: ',trim(misfit%name)
	do j=1,TFList%n
	  misfitValue(j) = sum(misfit%value(:,j))/(2*sum(misfit%ndat(:,j)))
	  write(0,*) 'Total data misfit for ',&
				  trim(TFList%info(j)%name),' response = ',misfitValue(j)
	end do

	v1=dot_product(misfit%weight,misfitValue)
	v2=dotProd_modelParam_f(p_input,p_input)
	v2=v2*misfit%damping

	write(0,*) 'Weighted model norm = ',v2

	write(0,*) 'Total data misfit   = ',v1+v2

	do j=1,TFList%n
	  misfitValue(j) = misfit%weight(j)*misfitValue(j)+v2
	end do

	!---------------------------------------------------------------------------
	! Important output to file; this will be a separate subroutine in output.f90
	call initFileWrite(cUserDef%fn_misfit, ioOut)

	write(ioOut,'(a10)',advance='no') 'MISFIT = '
	write (ioOut,'(g17.9)') v1+v2

	close(ioOut)
	!---------------------------------------------------------------------------

	! Output derivative summary

	dparam = param

	if(.not.present(dmisfitValue)) then
	  return
	end if

	do j=1,TFList%n
		write(comment,*) 'Total derivative for ',trim(TFList%info(j)%name),' response = '
!	  if (output_level>2) then
!	  write(0,*) 'Total derivative for ',&
!				  trim(TFList%info(j)%name),' response = '
!	  end if
	  do l=1,param%nc
		dmisfitValue(j,l) = sum(misfit%dRda(:,j,l))/(2*sum(misfit%ndat(:,j)))
	  end do
	  call fillParamValues_modelParam(dparam,dmisfitValue(j,:))
		call print_modelParam(dparam,output_level,trim(comment))
		call getCoeffArray_modelParam(dparam,grad)
!	  do l=1,param%nc
!		grad(l)=getCoeff_modelParam(dparam,l)
!		if (output_level>1) then
!		  write (0,'(2i8,g17.9)') grad(l)%L%num,grad(l)%F%num,grad(l)%value
!	        end if
!	  end do

	end do

	!if (output_level>1) then
	!  write(0,*) 'Weighted model norm = '
	!end if

	!do l=1,param%nc
	!  norm(l)=getCoeff_modelParam(p_input,l)
	!  if (output_level>1) then
	!    write (0,'(2i8,g17.9)') norm(l)%L%num,norm(l)%F%num,2.0d0*misfit%damping*norm(l)%value
	!  end if
	!end do
	weighted_norm = p_input
	call getCoeffArray_modelParam(p_input,norm)
	call scMult_modelParam(TWO*misfit%damping,p_input,weighted_norm)
	call print_modelParam(weighted_norm,output_level,"Weighted model norm derivative = ")

	! Compute total gradient including model norm term
	do l=1,param%nc
	  grad(l)%value=dot_product(misfit%weight,dmisfitValue(:,l))
	end do
	!dparam = dparam + weighted_norm

	norm%value = 2.0d0 * misfit%damping * norm%value
	grad%value = grad%value + norm%value
	call fillParam_modelParam(dparam,grad)


	call print_modelParam(dparam,output_level,"Total derivative = ")
!	do l=1,param%nc
!	  write (0,'(2i8,g17.9)') grad(l)%L%num,grad(l)%F%num,grad(l)%value
!	end do

	do j=1,TFList%n
	  dmisfitValue(j,:) = misfit%weight(j)*dmisfitValue(j,:)+norm%value
	end do

	call getCoeffArray_modelParam(p_input,point)

	!---------------------------------------------------------------------------
	! Important output to file; this will be a separate subroutine in output.f90
	call initFileWrite(cUserDef%fn_gradient, ioOut)

	write(ioOut,'(a15)',advance='no') 'DERIVATIVE = '
	do l=1,param%nc
	  if(grad(l)%frozen) then
		cycle
	  end if
	  write (ioOut,'(g17.9)',advance='no') grad(l)%value
	end do

	close(ioOut)
	!---------------------------------------------------------------------------
	! Important output to file; this will be a separate subroutine in output.f90
	call initFileWrite(cUserDef%fn_point, ioOut)

	write(ioOut,'(i1)') 1
	do l=1,param%nc
	  if(point(l)%frozen) then
		cycle
	  end if
	  write (ioOut,'(g17.9)',advance='no') point(l)%value
	end do

	close(ioOut)
	!---------------------------------------------------------------------------
	! Free memory
	call deall_modelParam(dparam)
	call deall_modelParam(weighted_norm)

  end subroutine misfitSumUp  ! misfitSumUp


end module dataMisfit
