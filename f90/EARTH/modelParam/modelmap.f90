! *****************************************************************************
module modelmap
  ! Module to initialize the model on the grid once the grid and parametrization
  ! information is available

  !use basics
  use griddef
  use modeldef
  use paramfunc
  implicit none


Contains


  ! ***************************************************************************
  ! * initResist initializes resistivity structure with the known values above
  ! * mantle; if crust distribution is known, this is also used; if not, this
  ! * case is handled by setting grid%nzCrust=grid%nzAir, so that no resistivity
  ! * values below air layers are initialized here.

  subroutine initResist(grid,crust,resist)

    type (grid_t), intent(in)					:: grid
	type (modelShell_t), intent(in)					:: crust
	real(8), dimension(:,:,:), intent(out)			:: resist !(nx,ny,nz)
	real(8)											:: crust_depth
	integer											:: i,j,k

	resist = R_ZERO

	forall (i=1:grid%nx, j=1:grid%ny, k=1:grid%nzAir)
	  resist(i,j,k) = 1/SIGMA_AIR
	end forall

	crust_depth = KM2M*(EARTH_R-CRUST_R)

	do i=1,grid%nx
	  do j=1,grid%ny
		do k=grid%nzAir+1,grid%nzCrust ! if no crust, nzCrust == nzAir
		  if(crust%cond(i,j) <= R_ZERO) then
			write(0, *) 'Error: (initResist) negative or infinite resistivity at',i,j,k
			stop
		  end if
		  resist(i,j,k) = crust_depth/crust%cond(i,j)
		end do
	  end do
	end do


  end subroutine initResist	! initResist


  ! ***************************************************************************
  ! * initModel is the routine that generates the 3-D resistivity map on the
  ! * grid using the information stored in the grid and the parametrization

  subroutine initModel(grid,param,resist)

    type (grid_t), intent(in)					:: grid
    type (modelParam_t), intent(inout)				:: param
	!type (rscalar), intent(out)					:: resist
	real(8), dimension(:,:,:), intent(inout)		:: resist !(nx,ny,nz)
	integer											:: i,j,k,l,istat
	integer											:: iL,ip
	real(8)											:: value,coeff
	type (modelLayer_t), pointer						:: this_layer
	type (modelFunc_t)								:: func
	type (modelPoint_t)								:: point

	! First initialize resistivity in air and possibly crust, if given
	call initResist(grid,param%crust,resist)

	! Test to make sure no grid is defined outside the layered region
	if (grid%r(grid%nz+1) < param%L(param%nL)%lbound) then
	  param%L(param%nL)%lbound = grid%r(grid%nz+1)
	end if

	do k=grid%nzCrust+1,grid%nz

	  ! Find current layer by locating the upper boundary of a cell
	  do l=1,param%nL
		if (in_layer(grid%r(k),param%L(l))) then
		  this_layer => param%L(l)
		  exit
		end if
	  end do

	  iL = this_layer%num

	  do i=1,grid%nx
		do j=1,grid%ny

		  point%phi   = (grid%ph(i) + grid%ph(i+1))/2
		  point%theta = (grid%th(j) + grid%th(j+1))/2
		  point%r     = (grid%r(k) + grid%r(k+1))/2

		  ! Sum up the coeffs * F_at_point in the given layer
		  value = 0.0d0
		  do ip=1,param%nF

			coeff = param%c(iL,ip)%value
			func = param%F(ip)
			if(coeff /= 0.0d0) then
			  value = value + coeff * F_at_point(func,point)
			end if

		  end do

		  if (this_layer%if_log) then
			resist(i,j,k) = exp(log(10.0d0)*value)
		  else
			resist(i,j,k) = value
		  end if

		  if (resist(i,j,k) <= 0.0d0) then
			write(0, *) 'Error: (initModel) negative or zero resistivity at',i,j,k
			stop
		  end if

		end do
	  end do
	end do


  end subroutine initModel	! initModel


end module modelmap
