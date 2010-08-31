! *****************************************************************************
module output
  ! Module containing the subroutines for output only

  ! NEED TO CLEAN UP THIS MODULE

  use iotypes
  use iospec
  use datadef
  use griddef
  use sg_vector
  use sg_sparse_vector
  use modeldef
  !use global
  use responses
  use functionals
  use dataspace
  use dataTypes
  use transmitters
  use receivers
  implicit none


Contains

  ! ***************************************************************************
  ! * Output full solution into separate files called "slices": each file
  ! * contains values for a single frequency and a single radius

  subroutine outputField(freq,H,cUserDef,extension)

	type (input_info), intent(in)			  :: cUserDef
	type (transmitter_t), intent(in)			  :: freq
	type (cvector), intent(in)				  :: H
	integer									  :: i,j,k,n,ios,istat
	character(3)							  :: ichar
	character(10)							  :: echar
	character(*), intent(in), optional		  :: extension
	character(80)							  :: fn_output

	write (ichar,'(i3.3)') freq%i
	if (present(extension)) then
	  echar = trim(extension)
	else
	  echar = 'H'
	end if
	fn_output = trim(cUserDef%modelname)//'_'//trim(freq%code)//'.'//trim(echar)
	!print *, "Output solution to ",fn_output

	open(ioWRITE,file=fn_output,status='unknown',form='formatted',iostat=ios)
	write(ioWRITE,'(a33,f0.3,a6)') "# Full H-field output for period ",freq%period," days."
	call write_cvector(ioWRITE,H)
	close(ioWRITE)

  end subroutine outputField	! outputField


  ! ***************************************************************************
  ! * Output full solution into separate files called "slices": each file
  ! * contains values for a single frequency and a single radius

  subroutine outputSolution(freq,H,slices,grid,cUserDef,rho,extension)

	type (input_info), intent(in)			  :: cUserDef
	type (transmitter_t), intent(in)			  :: freq
	type (Rad_List), intent(in)				  :: slices
	type (grid_t), intent(in)			  :: grid
	type (cvector), intent(in)				  :: H
	type (solution_t)							  :: Hij
	type (receiver_t)							  :: obs
	real(8), dimension(:,:,:), intent(in)	  :: rho	!(nx,ny,nz)
	complex(8)								  :: Hx,Hy,Hz
	integer									  :: i,j,k,n,ios,istat
	character(3)							  :: ichar
	character(10)							  :: echar
	character(*), intent(in), optional		  :: extension
	character(80)							  :: fn_output

	write (ichar,'(i3.3)') freq%i
	if (present(extension)) then
	  echar = trim(extension)
	else
	  echar = 'h'
	end if
	fn_output = trim(cUserDef%modelname)//'_'//trim(freq%code)//'.'//trim(echar)
	!print *, "Output solution to ",fn_output

	open(ioWRITE,file=fn_output,status='unknown',form='formatted',iostat=ios)
	write(ioWRITE,'(a85)') "#Radius ij_code GM_Lon GM_Lat rho(i,j,k) rho(i,j,k-1) Hx Hy Hz C_ratio D_ratio (km)";
	write(ioWRITE,'(a8,f0.3,a6,i6,a4,i6,a4,i6)') &
	  ' period=',freq%period,' nrad=',slices%n,' nx=',grid%nx,' ny=',grid%ny

	do n = 1,slices%n
	  ! Check that radius is valid
	  if ((slices%r(n) < grid%r(grid%nz)).or.(slices%r(n) > grid%r(2))) then
		write(*,*) 'Warning: output radius ',slices%r(n),' too close to a boundary; will not output data'
		cycle
	  end if
	  ! If it is, find the closest grid radius and compute the solution
	  call fieldValue_ij(slices%r(n),H,Hij)
	  do j = 2,grid%ny-1
		do i = 1,grid%nx
		  k = Hij%o(i,j)%k
		  obs = Hij%o(i,j)
		  Hx = Hij%x(i,j)
		  Hy = Hij%y(i,j)
		  Hz = Hij%z(i,j)
		  write(ioWRITE,'(f12.3,a10,14g15.7)') &
			  obs%rad,trim(obs%code),&
			  obs%lon,obs%lat,&
			  rho(i,j,k),rho(i,j,k-1),&
			  Hx,Hy,Hz,&
			  m2km * C_ratio(Hx,Hy,Hz,obs%colat*d2r),&
			  m2km * D_ratio(Hx,Hy,Hz,obs%colat*d2r)
		end do
	  end do
	end do

	close(ioWRITE)

	! Deallocate solution
	deallocate(Hij%x,Hij%y,Hij%z,STAT=istat)
	do i = 1,grid%nx
	  do j = 1,grid%ny
	    call deall_sparsevecc(Hij%o(i,j)%Lx)
	    call deall_sparsevecc(Hij%o(i,j)%Ly)
	    call deall_sparsevecc(Hij%o(i,j)%Lz)
	  end do
	end do
	deallocate(Hij%o,STAT=istat)

  end subroutine outputSolution	! outputSolution


  ! ***************************************************************************
  ! * OutputResponses writes the chosen kind of responses calculated at every
  ! * observatory location to an output file

  subroutine outputResponses(freq,psi,freqList,TFList,obsList,outFiles,dat)

	type (transmitter_t), intent(in)					:: freq
	type (dataVectorMTX_t), intent(in)					:: psi
	type (dataVectorMTX_t), intent(in), optional		    :: dat
	type (TF_List), intent(in)						:: TFList
	type (Obs_List), intent(in)						:: obsList
	type (Freq_List), intent(in)					:: freqList
	type (output_info), intent(in)					:: outFiles
	character(80)									:: fn_response
	complex(8), dimension(:), allocatable			:: Resp,RespRatio
	complex(8), dimension(:), allocatable			:: FieldData
	real(8), dimension(:), allocatable			:: FieldError
	real(8)											:: rval,ival,err,rms
	integer											:: i,j,k
	integer											:: ios,istat

	i = freq%i

	allocate(Resp(obsList%n),RespRatio(obsList%n),FieldData(obsList%n),FieldError(obsList%n), STAT=istat)

	do j=1,TFList%n

	  select case ( trim(TFList%info(j)%name) )
	  case ('C')
		fn_response = outFiles%fn_cdat
		RespRatio(:) = psi%v(i,j,:)%resp%value
		FieldData(:) = dat%v(i,j,:)%resp%value
		FieldError(:) = dat%v(i,j,:)%resp%err
		do k=1,obsList%n
			Resp(k) = RespRatio(k) * dtan(obsList%info(k)%colat*d2r) * m2km
		end do

	  case ('D')
		fn_response = outFiles%fn_ddat
		RespRatio(:) = psi%v(i,j,:)%resp%value
		FieldData(:) = dat%v(i,j,:)%resp%value
		FieldError(:) = dat%v(i,j,:)%resp%err
		do k=1,obsList%n
			Resp(k) = RespRatio(k) * dsin(obsList%info(k)%colat*d2r) * m2km
		end do

	  case default
		write(0,*) 'Warning: unknown transfer function: ',&
		  trim(TFList%info(j)%name)
		cycle
	  end select

!          inquire(FILE=fn_response,EXIST=exists)
!	  if ((.not.exists).and.(i==1)) then
!		open(ioResp,file=fn_response,status='unknown',form='formatted',iostat=ios)
!	    write(ioResp,*) "#Output of earth3d (for ",freqList%n," frequency values)";
!		write(ioResp,*) "#Period Code GM_Lon GM_Lat Real(km) Imag(km) Error(km)";
!	  else if(exists .and.(i==1)) then
!                write(0,*) 'Warning: Response file already exists. Appending...'
!		open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
!          else
!		open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
!	  end if

		! NB: No appending for now to make it easier to debug and cleaner. Jan 19, 2007
	  if (i==1) then
		  open(ioResp,file=fn_response,status='unknown',form='formatted',iostat=ios)
	    write(ioResp,*) "#Output of earth3d (for ",freqList%n," frequency values)";
		  write(ioResp,*) "#Period Code GM_Lon GM_Lat Real(km) Imag(km) Error(km)";
    else
		  open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
	  end if
		!write(ioResp,*) "freq = ",freq%value
		do k=1,size(RespRatio)
		  if (.not.obsList%info(k)%defined) then
			cycle
		  end if
		  write(ioResp,'(f8.3,a12,4g15.7)',advance='no') &
			  freqList%info(i)%period,&
			  trim(obsList%info(k)%code),&
			  obsList%info(k)%lon,obsList%info(k)%lat,Resp(k)
		  if(present(dat)) then
			if (dat%v(i,j,k)%resp%exists) then
			  rval = dreal(FieldData(k)-RespRatio(k))
			  ival = dimag(FieldData(k)-RespRatio(k))
			  err = FieldError(k)
			  rms = ((rval/err)**2 + (ival/err)**2)/2
			  write(ioResp,'(g15.7)') rms
			else
			  write(ioResp,'(g15.7)') 999999.9
			end if
		  else
			write(ioResp,*)
		  end if
		end do
	  close(ioResp)

!	  if (fn_response == '') then
!		do k=1,size(Resp)
!		  write(*,*) trim(obsList%info(k)%code),Resp(k) * m2km
!		end do
!	  end if

	end do

	deallocate(Resp)
	return

  end subroutine outputResponses  ! outputResponses


  ! ***************************************************************************
  ! * OutputMisfit outputs misfit and data sensitivities information

  subroutine outputMisfit(param,misfit,misfitValue,cUserDef)

	type (input_info), intent(in)			  :: cUserDef
	type (misfit_t), intent(in)			  :: misfit
	type (modelParam_t), intent(in)			  :: param
	real(8),dimension(:),intent(in)			  :: misfitValue
    character(80)                             :: fn_misfit

        fn_misfit = cUserDef%fn_misfit

	call initFileWrite(fn_misfit, ioWRITE)

	write(ioWRITE,'(a10)',advance='no') 'MISFIT = '
	write (ioWRITE,'(g15.7)') dot_product(misfit%weight,misfitValue)

	close(ioWRITE)

  end subroutine outputMisfit ! outputMisfit


  ! ***************************************************************************
  ! * OutputDerivative outputs the data sensitivities information

  subroutine outputDerivative(param,misfit,dmisfitValue,cUserDef)

	type (input_info), intent(in)			  :: cUserDef
	type (misfit_t), intent(in)			  :: misfit
	type (modelParam_t), intent(in)			  :: param
	real(8),dimension(:,:),intent(in)		  :: dmisfitValue
	integer									  :: i,j,l
    character(80)                             :: fn_misfit

        fn_misfit = cUserDef%fn_gradient

	call initFileWrite(fn_misfit, ioWRITE)

	write(ioWRITE,'(a15)',advance='no') 'DERIVATIVE = '
	do l=1,param%nc
	  !write (ioWRITE,'(g15.7)',advance='no') misfitInfo%d%c(j,i)%value
	  write (ioWRITE,'(g15.7)',advance='no') dot_product(misfit%weight,dmisfitValue(:,l))
	end do

	close(ioWRITE)

  end subroutine outputDerivative ! outputDerivative


  ! *****************************************************************************
  ! * OutputResiduals writes out the non-averaged residuals for all observatories
  ! * for each frequency and data functional

  subroutine outputResiduals(freq,res,TFList,obsList,outFiles)

	type (transmitter_t), intent(in)					:: freq
	type (dataVectorMTX_t), intent(in)					:: res
	type (TF_List), intent(in)						:: TFList
	type (Obs_List), intent(in)						:: obsList
	type (output_info), intent(in)					:: outFiles
	character(80)									:: fname
	integer											:: i,j,k,ios
	real(8)											:: rval,ival,error

	i = freq%i

	fname = outFiles%fn_residuals

	if (i == 1) then
	  open(ioWRITE,file=fname,status='unknown',form='formatted',iostat=ios)
	  write(ioWRITE,'(a40)') 'ifreq,ifunc,code,res,err,res/err = ';
	else
	  open(ioWRITE,file=fname,position='append', form='formatted',iostat=ios)
	end if

	do j = 1,TFList%n
	  do k = 1,obsList%n

		if (.not.res%v(i,j,k)%resp%exists) then
		  cycle
		end if

		rval = dreal(res%v(i,j,k)%resp%value)
		ival = dimag(res%v(i,j,k)%resp%value)
		error = res%v(i,j,k)%resp%err

  		write(ioWRITE,'(2i5,a10,5g15.7)') i,j,&
		  adjustr(trim(res%v(i,j,k)%obs%code)),&
		  rval,ival,error,rval/error,ival/error

	  end do
	end do

  end subroutine outputResiduals	! outputResiduals


  ! ***************************************************************************
  ! * OutputSens outputs the full data sensitivities information

  subroutine outputSens(misfitInfo,param,freqList,TFList)

	type (misfitInfo_t),dimension(:,:),intent(in)	  :: misfitInfo
	type (modelParam_t), intent(in)			  :: param
	type (Freq_List), intent(in)			  :: freqList
	type (TF_List), intent(in)				  :: TFList
	integer									  :: ifreq,ifunc,ilayer,icoeff,ivar
	real(8)									  :: misfitValue,dmisfitValue
        character(80)                                     :: fn_misfit

        fn_misfit = 'info.out'

	call initFileWrite(fn_misfit, ioWRITE)

	write(ioWRITE,*) '# ifreq,ifunc,ilayer,icoeff,code,coeff,misfit,derivative:'
	do ifreq=1,freqList%n
	  do ifunc=1,TFList%n
		do ilayer=1,param%nL
		  do icoeff=1,param%nF

			if(param%c(ilayer,icoeff)%frozen) then
			  cycle
			end if

			misfitValue = misfitInfo(ifreq,ifunc)%v/(2*misfitInfo(ifreq,ifunc)%n)
			dmisfitValue = misfitInfo(ifreq,ifunc)%d%c(ilayer,icoeff)%value

			write(ioWRITE,*) ifreq,ifunc,ilayer,icoeff,param%c(ilayer,icoeff)%code,&
					   param%c(ilayer,icoeff)%value,&
					   misfitValue,dmisfitValue

		  end do
		end do
	  end do
	end do

  end subroutine outputSens ! outputSens


  ! ***************************************************************************
  ! * OutputJacobian outputs the full Jacobian matrix

  subroutine outputJacobian(freq,psi,sens,freqList,TFList,obsList,param,cUserDef)

    use model_operators

	type (transmitter_t), intent(in)					:: freq
        type (sensitivity_t), intent(in)                                  :: sens
	type (dataVectorMTX_t), intent(in)					:: psi
	type (modelParam_t), intent(in), optional		    :: param
	type (TF_List), intent(in)						:: TFList
	type (Obs_List), intent(in)						:: obsList
	type (Freq_List), intent(in)					:: freqList
	!type (output_info), intent(in)					:: outFiles
	type (modelCoeff_t)								:: coeff
	type (input_info), intent(in)			  :: cUserDef
	character(3)							  :: ichar
	character(10)							  :: echar
	character(80)									:: fn_jacobian
	complex(8), dimension(:), allocatable			:: Resp
	real(8)											:: rval,ival,error,rms
	integer											:: i,j,k,l
	integer											:: ios,istat

	i = freq%i

	!write (ichar,'(i3.3)') freq%i

	allocate(Resp(obsList%n), STAT=istat)

	do j=1,TFList%n

	  select case ( trim(TFList%info(j)%name) )
	  case ('C')
                echar = 'cj'
	  case ('D')
                echar = 'dj'
	  case default
		write(0,*) 'Warning: unknown transfer function: ',&
		  trim(TFList%info(j)%name)
		cycle
	  end select

          fn_jacobian = trim(cUserDef%modelname)//'_'//trim(freq%code)//'.'//trim(echar)
	  print *, "Outputting Jacobian to ",fn_jacobian

	  Resp(:) = psi%v(i,j,:)%resp%value

          open(ioSens,file=fn_jacobian,status='unknown',form='formatted',iostat=ios)
          !write(ioSens,*) "#Full Jacobian, output of earth3d (for ",freqList%n," frequency values)";
          !write(ioSens,*) "#Period Code GM_Lon GM_Lat Layer Coeff.code Coeff.value Re(psi) Im(psi) Re(dpsi) Im(dpsi)";

          do k=1,size(Resp)
            if (.not.obsList%info(k)%defined) then
              cycle
            end if
            do l=1,param%nc
               	coeff = getCoeff_modelParam(param,l)

                if (coeff%frozen) then
                   cycle
                end if

                write(ioSens,'(f8.3,a12,2g15.7,i4,a12,g15.7,2g15.7,2g15.7)',advance='yes') &
                     freqList%info(i)%period,&
                     trim(obsList%info(k)%code),&
                     obsList%info(k)%lon,obsList%info(k)%lat,&
                     coeff%L%num,coeff%code,coeff%value,&
                     Resp(k),&
                     sens%da(i,j,k,l)

            end do
          end do
	  close(ioSens)

	end do

	deallocate(Resp)
	return

  end subroutine outputJacobian ! outputJacobian



  ! ***************************************************************************
  ! * OutputModel writes the resistivities on the grid to an output file

  subroutine outputModel(fn_model,mygrid,rho)

	character(80), intent(in)		  :: fn_model
	type (grid_t), intent(in)	  :: mygrid
	real(8), dimension(:,:,:), intent(in)	:: rho	!(nx,ny,nz)
	real(8)													:: lon,lat,depth
	integer							  :: ios,i,j,k

	if (fn_model /= '') then
	  open(ioMdl,file=fn_model,status='unknown',form='formatted',iostat=ios)
                write(ioMdl,*) '# lon(i), lat(j), depth(k), rho(ijk) of cell node'
                write(ioMdl,'(3i3)') mygrid%nx,mygrid%ny,mygrid%nzEarth
		do k=mygrid%nzAir+1,mygrid%nz
		  do i=1,mygrid%nx
			do j=1,mygrid%ny
				lon = mygrid%ph(i)*r2d
				lat = 90.0d0-mygrid%th(j)*r2d
				depth = EARTH_R - mygrid%r(k)
			  write(ioMdl,'(3g15.7)',advance='no') lon,lat,depth
			  write(ioMdl,'(g15.7)') rho(i,j,k)
			end do
		  end do
		end do
	  close(ioMdl)
	else
	  do k=mygrid%nzAir,mygrid%nz
		i=mygrid%nx
		j=1
		write(*,'(a10,i3,i3,i3)',advance='no') 'i,j,k = ',i,j,k
		write(*,*) ' rho(ijk) = ',rho(i,j,k)
	  end do
	  do k=mygrid%nzAir,mygrid%nz
		write(*,*) 'k,radius(k) = ',k,mygrid%r(k)
	  end do
	end if

  end subroutine outputModel  ! outputModel



end module output
