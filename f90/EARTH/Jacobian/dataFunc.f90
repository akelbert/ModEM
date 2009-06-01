! *****************************************************************************
module dataFunc
  ! Module containing the subroutines for a-posteriori analysis of the output
  ! vector h=lH from the routine SolveMaxwells

  use global
  use modeldef
  use datadef
  use griddef
  use sg_vector
  use sg_sparse_vector
  use interp
  use dataspace
  implicit none


Contains

  ! ***************************************************************************
  ! * operatorG acts on an input vector dvec length m = # of complex data values,
  ! * to produce a full complex vector defined on edges of the grid.
  ! * Matrix G consists of m vectors g(1),..,g(m), each defined on grid edges.
  ! * The output vector is G * dvec, computed by matrix rules of multiplication.
  ! * The reason we need to pass dvec with intent(inout) is that the observatory
  ! * info in dvec could be updated through the LocateReceiver subroutine.
  ! * An observatory is 'not defined' if it is too close to one of the poles.
  ! * The vector dvec should not contain data for these observatories in the first
  ! * place. These are just the checks to avoid segmentation fault.
  ! * Output vector should not contain any non-zero boundary values.

  subroutine operatorG(dvec,H,outE)

	! uses: funcList, obsList, freqList, grid
	type (dataValue_t), dimension(:), intent(in)		:: dvec	! residuals
	type (cvector), intent(in)							:: H
	type (cvector), intent(out)							:: outE
	type (sparsevecc)									:: g_sparse
	integer												:: j,istat
	integer												:: m !# of complex data pts

	! Can choose whether we want the input vector on which G acts to contain
	! a single type of responses, or a mixed set of different types; no limitations
	! or specific requirements on the number of data points in the input vector.
	m = size(dvec)

	! Create the full vector on edges that will be the output
	call create_cvector(grid,outE,EDGE)

	! Perform matrix multiplication column by column
	do j = 1,m
  	  if (.not.dvec(j)%obs%defined) then
		!write(0,*) 'Observatory ',trim(dvec(j)%obs%code),' is not defined at operatorG'
		cycle
	  end if
	  ! Compute the column vector of matrix G, which is g_j
	  call compute_g(dvec(j)%func,dvec(j)%obs,H,g_sparse)
	  ! Output vector x = G*r = \sum_{j=1}^m {g_j r_j}
	  call add_scvector(dvec(j)%resp%value,g_sparse,outE)
	end do

	call deall_sparsevecc(g_sparse)
	return	! return the outE

  end subroutine operatorG	! operatorG


  ! ***************************************************************************
  ! * operatorGt acts on an input vector defined on edges, to produce an output
  ! * data vector of length m; $G^* h = (g_1^* h,..,g_m^* h)^T$ is a vector of
  ! * linearised data functionals.

  subroutine operatorGt(inE,H,dvec)

	! uses: funcList, obsList, freqList, grid
	type (dataValue_t), dimension(:), intent(inout)		:: dvec	! residuals
	type (cvector), intent(in)							:: H
	type (cvector), intent(in)							:: inE
	type (sparsevecc)									:: g_sparse
	integer												:: j,istat
	integer												:: m !# of complex data pts

	! Can choose whether we want the input vector on which G acts to contain
	! a single type of responses, or a mixed set of different types; no limitations
	! or specific requirements on the number of data points in the input vector.
	m = size(dvec)

	! Perform matrix multiplication column by column
	do j = 1,m
  	  if (.not.dvec(j)%obs%defined) then
		!write(0,*) 'Observatory ',trim(dvec(j)%obs%code),' is not defined at operatorG'
		cycle
	  end if
	  ! Compute the column vector of matrix G, which is g_j
	  call compute_g(dvec(j)%func,dvec(j)%obs,H,g_sparse)
	  ! Compute the j'th entry of the data vector, g_j^* inE
	  dvec(j)%resp%value = dotProd_scvector_f(g_sparse,inE)
	end do

	call deall_sparsevecc(g_sparse)
	return	! return the outE

  end subroutine operatorGt	! operatorGt



  ! ***************************************************************************
  ! * compute_g is a subroutine to output a full vector g_j defined on edges,
  ! * such that for a single frequency and a single observatory,
  ! * $\pd{psi_{\omega}^j}{veca} = g_j^* \pd{vecH}{veca}$

  subroutine compute_g(dataType,obs,H,g_sparse)

	! uses: grid
	type (functional_t), intent(in)					:: dataType
	type (receiver_t), intent(in)					:: obs
	type (cvector), intent(in)						:: H
	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz
	complex(8)										:: pd_Hx,pd_Hy,pd_Hz
	type (sparsevecc)								:: gc_sparse  ! g*
	type (sparsevecc), intent(out)					:: g_sparse	! g
	real(8)											:: EARTH_R

	if (.not.obs%located) then
	  write(0,*) 'Error: Observatory ',trim(obs%code),' not yet located in compute_g'
	  stop
	  !call LocateReceiver(grid,obs)
	end if

	if (.not.obs%defined) then
	  write(0,*) 'Error: Observatory ',trim(obs%code),' is not defined at compute_g'
	  stop
	end if


	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H)
	Hz = dotProd_noConj(Lz,H)

	EARTH_R = grid%z(grid%nzAir+1) * m2km

	if (dataType%name == 'C') then

!	  pd_Hx = C_ZERO
!	  pd_Hy = - km2m * (EARTH_R/2) * dtan(obs%colat*d2r) * Hz/(Hy*Hy)
!	  pd_Hz =   km2m * (EARTH_R/2) * dtan(obs%colat*d2r) * 1/Hy
	  pd_Hx = C_ZERO
	  pd_Hy = - km2m * (EARTH_R/2) * Hz/(Hy*Hy)
	  pd_Hz =   km2m * (EARTH_R/2) * 1/Hy

	  call linComb_sparsevecc(Ly,pd_Hy,Lz,pd_Hz,gc_sparse)

	else if (dataType%name == 'D') then

!	  pd_Hx =   km2m * (EARTH_R/2) * dsin(obs%colat*d2r) * 1/Hy
!	  pd_Hy = - km2m * (EARTH_R/2) * dsin(obs%colat*d2r) * Hx/(Hy*Hy)
!	  pd_Hz = C_ZERO
	  pd_Hx =   km2m * (EARTH_R/2) * 1/Hy
	  pd_Hy = - km2m * (EARTH_R/2) * Hx/(Hy*Hy)
	  pd_Hz = C_ZERO

	  call linComb_sparsevecc(Ly,pd_Hy,Lx,pd_Hx,gc_sparse)

	else

	  write(0, *) 'Error: (compute_g) unknown data functional', dataType%name
	  stop

	end if

	! This was g* (conjugated), now we want g
	g_sparse = conjg(gc_sparse)

	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)
	call deall_sparsevecc(gc_sparse)

  end subroutine compute_g	! compute_g



end module dataFunc
