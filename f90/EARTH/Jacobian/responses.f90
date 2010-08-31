! *****************************************************************************
module responses
	! Module responses contains the variables and the basic subroutines to do
	! with the data responses and the penalty functional and their derivatives

  use SolnSpace
  use receivers
  use global, only: grid
  implicit none


Contains

  ! ***************************************************************************
  ! * dataResp performs the calculation that returns a response of a
  ! * chosen kind at a chosen observatory location, given H defined on grid
  ! * We are likely to need the responses at observatories (use this subroutine)
  ! * and at grid nodes on the surface. The latter is very easy to compute.
  ! * Use e.g. C(H%x(i,j,nzAir+1),H%y(i,j,nzAir+1),H%z(i,j,nzAir+1),grid%th(j)).
  ! * For responses at locations other than the observatory, construct an
  ! * observatory type (receiver_t) at the required location and call this routine.
  ! * Interpolation is automatic through LocateReceiver->ComputeInterpWeights.

  function dataResp(func,obs,H) result(Resp)

	complex(8), external							:: func
	type (receiver_t), intent(inout)					:: obs
	type (cvector), intent(in)						:: H
	complex(8)										:: Resp

	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz


	if (.not.obs%located) then
	  call LocateReceiver(grid,obs)
	end if

	if (.not.obs%defined) then
	  !write(0,*) 'Observatory ',trim(obs%code),' is not defined at dataResp; hence ignore'
	  return
	end if

	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H)
	Hz = dotProd_noConj(Lz,H)

	! Calculate the response at this location
	Resp = func(Hx,Hy,Hz,obs%colat*d2r)

	!write(*,'(a40,2i6,5g15.7)') 'i,j,lon,colat,Hx,Hy,Hz = ',&
	!	obs%i,obs%j,obs%lon,obs%colat,Hx,Hy,Hz

	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)

  end function dataResp ! response


  ! ***************************************************************************
  ! * dataResp_ijk performs a calculation such that the output is
  ! * the response at a cell (defined at the midpoint of the interval from the
  ! * (i,j,k)'th to the (i,j+1,k)'th node), given H defined on grid
  ! * With the current interpolation routine, k is assumed to be the index
  ! * of the Earth's surface

  function dataResp_ijk(func,i,j,k,H) result(Resp)

	complex(8), external							:: func
	integer, intent(in)								:: i,j,k
	type (cvector), intent(in)						:: H
	complex(8)										:: Resp
	type (receiver_t)									:: obs
	real(8)											:: rad,lon,colat
	character(80)									:: name

	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz

	rad = grid%r(k)
	lon = grid%ph(i)*r2d
	colat = (grid%th(j)+grid%th(j+1))*r2d/2

	write(name,'(2i5)') i,j
	name = trim(name)

	! Build an observatory at this cell with code = name
	call CreateReceiver(obs,rad,lon,colat,name)

	! Compute interpolation parameters
	call LocateReceiver(grid,obs)

	if (.not.obs%defined) then
	  return
	end if

	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H)
	Hz = dotProd_noConj(Lz,H)

	! Calculate the response at this location
	Resp = func(Hx,Hy,Hz,obs%colat*d2r)

	!write(*,'(a40,2i6,5g15.7)') 'i,j,lon,colat,Hx,Hy,Hz = ',&
	!	obs%i,obs%j,obs%lon,obs%colat,Hx,Hy,Hz
	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)

  end function dataResp_ijk ! dataResp_ijk


  ! ***************************************************************************
  ! * fieldValue_ij computes a set of magnetic field solutions
  ! * at the cells (defined at the midpoint of the interval from the
  ! * (i,j,k)'th to the (i,j+1,k)'th node), given H defined on grid, at
  ! * the specified grid radius

  subroutine fieldValue_ij(rad,H,Hij)

	real(8), intent(in)								:: rad
	type (cvector), intent(in)						:: H
	type (solution_t), intent(out)					:: Hij
	real(8)											:: lon,colat
	character(80)									:: name
	integer											:: i,j,istat
	integer											:: nx,ny

	nx = grid%nx
	ny = grid%ny

	allocate(Hij%x(nx,ny),Hij%y(nx,ny),Hij%z(nx,ny),Hij%o(nx,ny),STAT=istat)

	Hij%x = C_ZERO
	Hij%y = C_ZERO
	Hij%z = C_ZERO

	do j = 2,ny-1
	  do i = 1,nx

		!Use if obs is at midpoint between (i,j,k)'th to the (i,j+1,k)'th nodes
		!lon = grid%ph(i)*r2d
		!Use if obs is at the center of the (i,j,k)'th cell
		if (i == nx) then
		  lon = (grid%ph(i)+2*PI)*r2d/2
		else
		  lon = (grid%ph(i)+grid%ph(i+1))*r2d/2
		end if
		colat = (grid%th(j)+grid%th(j+1))*r2d/2

		write(name,'(2i4)') i,j
		name = trim(name)

		! Build an observatory at this cell with code = name
		call CreateReceiver(Hij%o(i,j),rad,lon,colat,name)

		! NB: The receiver is shifted down to the nearest grid radius!!!
		call LocateReceiver(grid,Hij%o(i,j))

		if (Hij%o(i,j)%located) then
		  Hij%x(i,j) = dotProd_noConj(Hij%o(i,j)%Lx,H) ! dotProd_noConj_scvector_f
		  Hij%y(i,j) = dotProd_noConj(Hij%o(i,j)%Ly,H)
		  Hij%z(i,j) = dotProd_noConj(Hij%o(i,j)%Lz,H)
		end if

	  end do
	end do

  end subroutine fieldValue_ij ! fieldValue_ij


end module responses
