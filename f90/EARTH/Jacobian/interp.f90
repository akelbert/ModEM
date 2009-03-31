! *****************************************************************************
module interp
  ! Module containing the Earth's surface interpolation subroutines

  use modeldef
  use datadef
  use sg_vector
  use sg_sparse_vector
  use utilities
  implicit none


Contains

  ! ***************************************************************************
  ! * initObsList initializes the interpolation parameters for a pre-defined
  ! * list of receivers; call in the initialization block, so that we know which
  ! * receivers are defined while initializing the data

  subroutine initObsList(grid,obsList)

	type (grid_t), intent(in)				:: grid
	type (Obs_List), intent(inout)				:: obsList
	integer										:: nobs,i

	nobs = obsList%n
	do i=1,nobs
	  call LocateReceiver(grid,obsList%info(i))
	  !write(*,'(a40,a6,2i6,2g15.7)') 'code,i,j,lon,colat = ',obsList%info(i)%code,&
	  !	obsList%info(i)%i,obsList%info(i)%j,obsList%info(i)%lon,obsList%info(i)%colat
	end do

  end subroutine initObsList


  ! ***************************************************************************
  ! * ComputeInterpWeights computes the sparse vector L at a single observatory,
  ! * such that H%x @ obs = L%x * H etc, using Bilinear Interpolation.
  ! * Currently the weights are such that they H @ obs would compute as
  ! x:
  !   If observatory is below the mid-longitude of the cell face,
  !	  surfH%x(iobs) = (1-q2)*((1-p1)*H%x(i,j,k) + p1*H%x(i,j+1,k))
  !					+ q2 * ((1-p1)*H%x(i-1,j,k) + p1*H%x(i-1,j+1,k))
  !   If observatory is above the mid-longitude of the cell face,
  !	  surfH%x(iobs) = (1-q2)*((1-p1)*H%x(i,j,k) + p1*H%x(i,j+1,k))
  !					+ q2 * ((1-p1)*H%x(i+1,j,k) + p1*H%x(i+1,j+1,k))
  !	  Else
  !	  surfH%x(iobs) = (1-p1)*H%x(i,j,k) + p1*H%x(i,j+1,k)
  !	y:
  !   If observatory is below the mid-latitude of the cell face,
  !	  surfH%y(iobs) = (1-p2)*((1-q1)*H%y(i,j,k) + q1*H%y(i+1,j,k))
  !					+ p2 * ((1-q1)*H%y(i,j-1,k) + q1*H%y(i+1,j-1,k))
  !   If observatory is above the mid-latitude of the cell face,
  !	  surfH%y(iobs) = (1-p2)*((1-q1)*H%y(i,j,k) + q1*H%y(i+1,j,k))
  !					+ p2 * ((1-q1)*H%y(i,j+1,k) + q1*H%y(i+1,j+1,k))
  !	  Else
  !	  surfH%y(iobs) = (1-q1)*H%y(i,j,k) + q1*H%y(i+1,j,k)
  ! z:
  !	  surfH%z(iobs) =  p1*q1 * (H%z(i,j,k-1)+H%z(i,j,k))/2 &
  !					+ p1*(1-q1) * (H%z(i+1,j,k-1)+H%z(i+1,j,k))/2 &
  !					+ q1*(1-p1) * (H%z(i,j+1,k-1)+H%z(i,j+1,k))/2 &
  !					+ (1-p1)*(1-q1) * (H%z(i+1,j+1,k-1)+H%z(i+1,j+1,k))/2
  !
  ! e.g. for z-component the weights are
  !	  w%z(i,j,k-1) = (1-p)*(1-q)/2
  !	  w%z(i,j,k) = (1-p)*(1-q)/2
  !	  w%z(i+1,j,k-1) = (1-p)*q/2
  !	  w%z(i+1,j,k) = (1-p)*q/2
  !	  w%z(i,j+1,k-1) = (1-q)*p/2
  !	  w%z(i,j+1,k) = (1-q)*p/2
  !	  w%z(i+1,j+1,k-1) = p*q/2
  !	  w%z(i+1,j+1,k) = p*q/2

  subroutine ComputeInterpWeights(grid,obs)

	type (grid_t) , intent(in)					:: grid
	type (receiver)	 , intent(inout)				:: obs
	type (sparsevecc)								:: Lx,Ly,Lz
	integer											:: i,j,k,index
	real(8)											:: p1,p2,q1,q2
	integer, parameter								:: x=1,y=2,z=3

	  i = obs%i
	  j = obs%j
	  k = obs%k
	  p1 = obs%p_crn  ! N-S distance from the cell corner
	  p2 = obs%p_ctr  ! N-S distance from the cell center
	  q1 = obs%q_crn  ! W-E distance from the cell corner
	  q2 = obs%q_ctr  ! W-E distance from the cell center

	  ! Create a sparse vector Lx and insert weights
	  if (clean(q1) > 0.5d0) then  ! above the mid-longitude
		call create_sparsevecc(4,Lx,EDGE)
		call newValueR_sparsevecc(Lx,1,(1-q2)*(1-p1),i,j,k,x)
		call newValueR_sparsevecc(Lx,2,(1-q2)*p1,i,j+1,k,x)
		call newValueR_sparsevecc(Lx,3,q2*(1-p1),i+1,j,k,x)
		call newValueR_sparsevecc(Lx,4,q2*p1,i+1,j+1,k,x)
	  else if (clean(q1) < 0.5d0) then ! below the mid-longitude
		call create_sparsevecc(4,Lx,EDGE)
		call newValueR_sparsevecc(Lx,1,(1-q2)*(1-p1),i,j,k,x)
		call newValueR_sparsevecc(Lx,2,(1-q2)*p1,i,j+1,k,x)
		call newValueR_sparsevecc(Lx,3,q2*(1-p1),i-1,j,k,x)
		call newValueR_sparsevecc(Lx,4,q2*p1,i-1,j+1,k,x)
	  else	! at the mid-longitude
		call create_sparsevecc(2,Lx,EDGE)
		call newValueR_sparsevecc(Lx,1,1-p1,i,j,k,x)
		call newValueR_sparsevecc(Lx,2,p1,i,j+1,k,x)
	  end if

	  ! Create a sparse vector Ly and insert weights
	  if (clean(p1) > 0.5d0) then  ! above the mid-latitude
		call create_sparsevecc(4,Ly,EDGE)
		call newValueR_sparsevecc(Ly,1,(1-p2)*(1-q1),i,j,k,y)
		call newValueR_sparsevecc(Ly,2,(1-p2)*q1,i+1,j,k,y)
		call newValueR_sparsevecc(Ly,3,p2*(1-q1),i,j+1,k,y)
		call newValueR_sparsevecc(Ly,4,p2*q1,i+1,j+1,k,y)
	  else if (clean(p1) < 0.5d0) then ! below the mid-latitude
		call create_sparsevecc(4,Ly,EDGE)
		call newValueR_sparsevecc(Ly,1,(1-p2)*(1-q1),i,j,k,y)
		call newValueR_sparsevecc(Ly,2,(1-p2)*q1,i+1,j,k,y)
		call newValueR_sparsevecc(Ly,3,p2*(1-q1),i,j-1,k,y)
		call newValueR_sparsevecc(Ly,4,p2*q1,i+1,j-1,k,y)
	  else	! at the mid-latitude
		call create_sparsevecc(2,Ly,EDGE)
		call newValueR_sparsevecc(Ly,1,1-q1,i,j,k,y)
		call newValueR_sparsevecc(Ly,2,q1,i+1,j,k,y)
	  end if

	  ! Create a sparse vector Lz and insert weights
	  call create_sparsevecc(8,Lz,EDGE)
	  call newValueR_sparsevecc(Lz,1,(1-p1)*(1-q1)/2,i,j,k-1,z)
	  call newValueR_sparsevecc(Lz,2,(1-p1)*(1-q1)/2,i,j,k,z)
	  call newValueR_sparsevecc(Lz,3,(1-p1)*q1/2,i+1,j,k-1,z)
	  call newValueR_sparsevecc(Lz,4,(1-p1)*q1/2,i+1,j,k,z)
	  call newValueR_sparsevecc(Lz,5,p1*(1-q1)/2,i,j+1,k-1,z)
	  call newValueR_sparsevecc(Lz,6,p1*(1-q1)/2,i,j+1,k,z)
	  call newValueR_sparsevecc(Lz,7,p1*q1/2,i+1,j+1,k-1,z)
	  call newValueR_sparsevecc(Lz,8,p1*q1/2,i+1,j+1,k,z)

	  ! To avoid problems on the zero longitude meridian
	  do index=1,Lx%nCoeff
		if (Lx%i(index) <= 0) then
		  Lx%i(index) = grid%nx - Lx%i(index)
		else if (Lx%i(index) > grid%nx) then
		  Lx%i(index) = Lx%i(index) - grid%nx
		end if
	  end do
	  do index=1,Ly%nCoeff
		if (Ly%i(index) <= 0) then
		  Ly%i(index) = grid%nx - Ly%i(index)
		else if (Ly%i(index) > grid%nx) then
		  Ly%i(index) = Ly%i(index) - grid%nx
		end if
	  end do
	  do index=1,Lz%nCoeff
		if (Lz%i(index) <= 0) then
		  Lz%i(index) = grid%nx - Lz%i(index)
		else if (Lz%i(index) > grid%nx) then
		  Lz%i(index) = Lz%i(index) - grid%nx
		end if
	  end do

	  ! Save the interpolation weights in the receiver variable
	  obs%Lx = Lx
	  obs%Ly = Ly
	  obs%Lz = Lz
	  return

  end subroutine ComputeInterpWeights	! ComputeInterpWeights


  ! ***************************************************************************
  ! * CreateReceiver creates an observatory with a specified name at a specified
  ! * location: longitude and colatitude (in degrees)

  subroutine CreateReceiver(o,rad,lon,colat,name)

	type (receiver), intent(out)	:: o
	real(8), intent(in)				:: rad,lon,colat
	character(80), intent(in)		:: name

	o%code = trim(name)
	o%rad = rad
	o%lon = lon
	o%colat = colat
	o%lat = 90.0 - colat
	o%defined = .TRUE.

  end subroutine CreateReceiver	! CreateReceiver


  ! ***************************************************************************
  ! * LocateReceiver is a subroutine that locates an observatory on the surface
  ! * grid. It computes the linear weights that will be used for
  ! * Bilinear Interpolation.
  ! *
  ! * Comment:
  ! * An alternative would be a subroutine that locates an observatory on the
  ! * surface grid, using the method of Spherical Linear Interpolation, devised
  ! * in Ken Shoemake, �Animating Rotation with Quaternion Curves�,
  ! * Computer Graphics, Volume 19, Number 3, 1985 (proof by Glassner, 1999).
  ! * In contrast to the Linear Interpolation, we would then compute T0,T1
  ! * such that h = T0*h0 + T1*h1, so that
  ! * T0 = \dfrac{\sin{(1-t)\Omega}}{\sin{Omega}} and
  ! * T1 = \dfrac{\sin{t\Omega}}{\sin{Omega}}.
  ! * Equivalently,
  ! * T0 = \dfrac{\sin{\Omega-\alpha}}{\sin{Omega}} and
  ! * T1 = \dfrac{\sin{\alpha}}{\sin{Omega}},
  ! *
  ! * In general, let $o_0$ and $o_1$ be the first and last points of the arc on
  ! * the geodesic (a great circle) which passes through our observatory, and
  ! * let $o$ denote the observatory itself. Consider them to be radial vectors
  ! * from the center of the Earth to the surface. Let $\Omega$ be the angle
  ! * subtended by the arc, while $\alpha$ be the angle corresponding to the
  ! * shortest distance from $o_0$ to $o$.
  ! *
  ! * Define the parameter $t\in [0,1]$, such that $t = \dfrac{\alpha}{\Omega}$.
  ! *
  ! * Then t and 1-t are effectively the weights required for interpolation,
  ! * except that they should be inverted: h = (1-t)*h0 + t*h1,
  ! * where t is the parameter in [0,1] (either p or q),
  ! * and h,h0,h1 are any values (either vectors or scalars), defined in the
  ! * points o,o0,o1 respectively.
  ! *
  ! * This subroutine outputs the grid indices of the cell on the upper face of
  ! * which the observatory is located. It also outputs (0,0)<=(p,q)<=(1,1),
  ! * where p and q are such that if
  ! * lxox = the shortest length between (i,j,k)'th and (i,j+1,k)'th x-edges
  ! * lxo  = the shortest length from the (i,j,k)'th x-edge to the observatory
  ! * lyoy = the shortest length between (i,j,k)'th and (i+1,j,k)'th y-edges
  ! * lyo  = the shortest length from the (i,j,k)'th y-edge to the observatory
  ! * then
  ! * lxox = grid%r(k) * (grid%th(j+1) - grid%th(j))
  ! * lxo  = grid%r(k) * (th - grid%th(j))
  !	* lyoy = grid%r(k) * dsin(th) * (grid%ph(i+1) - grid%ph(i))
  !	* lyo  = grid%r(k) * dsin(th) * (ph - grid%ph(i))
  ! * and p_crn,q_crn are the proportions:
  !	* p_crn = lxo/lxox
  !	* q_crn = lyo/lyoy
  ! *
  ! * Also, for bilinear interpolation we will require the proportions of the
  ! * perpendicular distances along the lines through the cell face centers.
  ! * These are p_ctr and q_ctr. Best described by a drawing that I will include
  ! * in my thesis.
  ! *
  ! * There might be a better way of doing the interpolation...

  subroutine LocateReceiver(grid,o)

	type (grid_t), intent(in)	:: grid
	type (receiver), intent(inout)	:: o
	integer							:: i,j,k
	real(8)							:: ph,th,r
	real(8)							:: p_crn,q_crn ! corner as origin
	real(8)							:: p_ctr,q_ctr ! center as origin
	real(8)							:: ph_crn,ph_crn_p,ph_crn_m
	real(8)							:: th_crn,th_crn_p,th_crn_m
	real(8)							:: ph_ctr,ph_ctr_p,ph_ctr_m
	real(8)							:: th_ctr,th_ctr_p,th_ctr_m
	real(8),dimension(:),allocatable:: ph_grid,th_grid
	integer							:: nx,ny,nz
	integer							:: istat

	if (.not.o%defined) then
	  !write(0,*) 'Error: (LocateReceiver) receiver ',trim(o%code),' not defined'
	  return
	else if (o%located) then
	  return
	end if

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz

	allocate(ph_grid(0:nx+2),th_grid(1:ny+1),STAT=istat)
	ph_grid(1:nx) = grid%ph(1:nx)
	ph_grid(0) = grid%ph(nx) - 2*pi	! required at i == 1
	ph_grid(nx+1) = grid%ph(1) + 2*pi ! required at i == nx-1 or nx
	ph_grid(nx+2) = grid%ph(2) + 2*pi ! required at i == nx
	th_grid(1:ny+1) = grid%th(1:ny+1)

	!k = grid%nzAir+1	! on the air/Earth boundary
	r = o%rad
	k = maxNode(r,grid%r(1:nz+1))

	! first element in the array has to correspond to zero longitude
	ph = o%lon * d2r
	i = minNode(ph,ph_grid(1:nx))

	! first element in the array has to correspond to zero co-latitude
	th = o%colat * d2r
	j = minNode(th,th_grid(1:ny))

	! Special case: receiver too close to one of the poles
	! At poles Hx is not defined, Hy ambiguous, tan(th) = 0
	if ((j == 1).or.(j == ny)) then
	  write(0,*) 'Warning: (LocateReceiver) receiver ',trim(o%code),' too close to poles'
	  o%i = i
	  o%j = j
	  o%k = k
	  o%defined = .FALSE. ! fields not defined, so ignore the receiver
	  return
	end if

	! Special case: receiver is located at equator (or vero close)
	! At equator tan(th)->inf,cos(th)=0, hence C- and D-responses aren't useful
!	if ((th >= pi/2-EPS_GRID).and.(th <= pi/2+EPS_GRID)) then
!	  write(0,*) 'Warning: (LocateReceiver) receiver ',trim(o%code),' located at equator'
!	  o%i = i
!	  o%j = j
!	  o%k = k
!	  o%defined = .FALSE.
!	  return
!	end if

	! Set angular distances from origin to cell corners
	ph_crn = ph_grid(i)
	ph_crn_p = ph_grid(i+1)
	ph_crn_m = ph_grid(i-1)
	th_crn = th_grid(j)
	th_crn_p = th_grid(j+1)
	th_crn_m = th_grid(j-1)

	! Set angular distances from origin to cell face centers
	ph_ctr = (ph_grid(i)+ph_grid(i+1))/2
	ph_ctr_p = (ph_grid(i+1)+ph_grid(i+2))/2
	ph_ctr_m = (ph_grid(i-1)+ph_grid(i))/2
	th_ctr = (th_grid(j)+th_grid(j+1))/2
	th_ctr_p = (th_grid(j+1)+th_grid(j+2))/2
	th_ctr_m = (th_grid(j-1)+th_grid(j))/2

	p_crn = clean((th - th_crn)/(th_crn_p - th_crn))
	q_crn = clean((ph - ph_crn)/(ph_crn_p - ph_crn))

	! Expressions against machine error problems
	! p_crn = dnint(p_crn*LARGE_REAL)/LARGE_REAL
	! q_crn = dnint(q_crn*LARGE_REAL)/LARGE_REAL

	if ((p_crn < 0.0d0).OR.(p_crn > 1.0d0)) then
	  write(0, *) 'Warning: (LocateReceiver) the value of p_crn is ',p_crn,' at ',trim(o%code)
	  if ((p_crn < 0.0d0).AND.(p_crn > - EPS_GRID)) then
		p_crn = 0.0d0
	  else
		stop
	  end if
	else if ((q_crn < 0.0d0).OR.(q_crn > 1.0d0)) then
	  write(0, *) 'Warning: (LocateReceiver) the value of q_crn is ',q_crn,' at ',trim(o%code)
	  if ((q_crn < 0.0d0).AND.(q_crn > - EPS_GRID)) then
		q_crn = 0.0d0
	  else
		stop
	  end if
	end if

	if (p_crn .lt. 0.5d0) then
	  p_ctr = clean((th_ctr - th)/(th_ctr - th_ctr_m))
	else if (p_crn .gt. 0.5d0) then
	  p_ctr = clean((th - th_ctr)/(th_ctr_p - th_ctr))
	else
	  p_ctr = 0.0d0
	end if

	if (q_crn .lt. 0.5d0) then
	  q_ctr = clean((ph_ctr - ph)/(ph_ctr - ph_ctr_m))
	else if (q_crn .gt. 0.5d0) then
	  q_ctr = clean((ph - ph_ctr)/(ph_ctr_p - ph_ctr))
	else
	  q_ctr = 0.0d0
	end if

	! Expressions against machine error problems
	! p_ctr = dnint(p_ctr*LARGE_REAL)/LARGE_REAL
	! q_ctr = dnint(q_ctr*LARGE_REAL)/LARGE_REAL

	! print *, 'i,j,p_crn,q_crn,p_ctr,q_ctr = ',i,j,p_crn,q_crn,p_ctr,q_ctr


	if ((p_ctr < 0.0d0).OR.(p_ctr > 1.0d0)) then
	  write(0, *) 'Warning: (LocateReceiver) the value of p_ctr is ',p_ctr,' at ',trim(o%code)
	  if ((p_ctr < 0.0d0).AND.(p_ctr > - EPS_GRID)) then
		p_ctr = 0.0d0
	  else
		stop
	  end if
	else if ((q_ctr < 0.0d0).OR.(q_ctr > 1.0d0)) then
	  write(0, *) 'Warning: (LocateReceiver) the value of q_ctr is ',q_ctr,' at ',trim(o%code)
	  if ((q_ctr < 0.0d0).AND.(q_ctr > - EPS_GRID)) then
		q_ctr = 0.0d0
	  else
		stop
	  end if
	end if

	! Proportions computed. Now save them in the receiver and exit.
	o%i = i
	o%j = j
	o%k = k
	o%p_crn = p_crn
	o%p_ctr = p_ctr
	o%q_crn = q_crn
	o%q_ctr = q_ctr

	! Initialize Lx,Ly,Lz
	call ComputeInterpWeights(grid,o)

	o%located = .TRUE.

	deallocate(ph_grid,th_grid)

	return

  end subroutine LocateReceiver	! LocateReceiver



end module interp
