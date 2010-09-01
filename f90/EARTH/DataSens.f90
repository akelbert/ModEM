! *****************************************************************************
! Temporary module. Will be deleted when we switch to the modular DataSpace.
module DataSens
  ! Lmult & LmultT

  use DataSpace
  use DataFunc
  use receivers
  implicit none


Contains

  ! ***************************************************************************
  ! * LmultT acts on an input vector dvec length m = # of complex data values,
  ! * to produce a full complex vector defined on edges of the grid.
  ! * Matrix G consists of m vectors g(1),..,g(m), each defined on grid edges.
  ! * The output vector is G * dvec, computed by matrix rules of multiplication.
  ! * The reason we need to pass dvec with intent(inout) is that the observatory
  ! * info in dvec could be updated through the LocateReceiver subroutine.
  ! * An observatory is 'not defined' if it is too close to one of the poles.
  ! * The vector dvec should not contain data for these observatories in the first
  ! * place. These are just the checks to avoid segmentation fault.
  ! * Output vector should not contain any non-zero boundary values.

  subroutine LmultT(H,m0,dvec,outE) !LmultT(dvec,H,outE)

	! uses: funcList, obsList, freqList, grid
  ! background model parameter
  type (modelParam_t), intent(in)  :: m0
  !  electric field solutions are stored as type solnVector
  type (solnVector_t), intent(in)     :: H
  ! d provides indices into receiver dictionary on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dataValue_t), dimension(:), intent(in)		:: dvec	! residuals
  !type (dataVector_t), intent(in)               :: dvec
  ! Output
  type (rhsVector_t), intent(inout)                 :: outE

	!type (cvector), intent(in)							:: H
	!type (cvector), intent(inout)							:: outE
	!type (sparsevecc)									:: g_sparse
    type (sparseVector_t)                                   :: g_sparse(1)
	integer												:: j,istat
	integer												:: m !# of complex data pts
    type(receiver_t), pointer                           :: obs

	! Can choose whether we want the input vector on which G acts to contain
	! a single type of responses, or a mixed set of different types; no limitations
	! or specific requirements on the number of data points in the input vector.
	m = size(dvec)

	! Create the full vector on edges that will be the output
	call create_cvector(grid,outE%source,EDGE)

	! Perform matrix multiplication column by column
	do j = 1,m
      obs => obsList%info(dvec(j)%rx)
	   print *, 'LmultT for observatory ',trim(obs%code)
  	  if (.not.obs%defined) then
		!write(0,*) 'Observatory ',trim(obs%code),' is not defined at LmultT'
		cycle
	  end if
	  ! Compute the column vector of matrix G, which is g_j
      call Lrows(H,m0,dvec(j)%dataType,dvec(j)%rx,g_sparse)
	  ! Output vector x = G*r = \sum_{j=1}^m {g_j r_j}
	  call add_scvector(dvec(j)%resp%value,g_sparse(1)%L,outE%source)
	end do

	call deall_sparsevecc(g_sparse(1)%L)
	return	! return the outE

  end subroutine LmultT	! LmultT


  ! ***************************************************************************
  ! * Lmult acts on an input vector defined on edges, to produce an output
  ! * data vector of length m; $G^* h = (g_1^* h,..,g_m^* h)^T$ is a vector of
  ! * linearised data functionals.

  subroutine Lmult(H,m0,inE,dvec) !Lmult(inE,H,dvec)

	! uses: funcList, obsList, freqList, grid
  !  electric field solutions are stored as type solnVector
  type (solnVector_t), intent(in)           :: H,inE
  ! d provides indices into receiver and data type dictionaries on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dataValue_t), dimension(:), intent(inout)		:: dvec	! residuals
  !type (dataVector_t), intent(inout)      :: dvec
  !  background model parameter
  type (modelParam_t), intent(in)           :: m0


	!type (cvector), intent(in)							:: H
	!type (cvector), intent(in)							:: inE
    !type (sparsevecc)                                   :: g_sparse
	type (sparseVector_t)								:: g_sparse(1)
	integer												:: j,istat
	integer												:: m !# of complex data pts
	type(receiver_t), pointer                           :: obs

	! Can choose whether we want the input vector on which G acts to contain
	! a single type of responses, or a mixed set of different types; no limitations
	! or specific requirements on the number of data points in the input vector.
	m = size(dvec)

	! Perform matrix multiplication column by column
	do j = 1,m
	  obs => obsList%info(dvec(j)%rx)
  	  if (.not.obs%defined) then
		!write(0,*) 'Observatory ',trim(obs%code),' is not defined at LmultT'
		cycle
	  end if
	  ! Compute the column vector of matrix G, which is g_j
	  call Lrows(H,m0,dvec(j)%dataType,dvec(j)%rx,g_sparse)
	  ! Compute the j'th entry of the data vector, g_j^* inE
	  dvec(j)%resp%value = dotProd_scvector_f(g_sparse(1)%L,inE%vec)
	end do

	call deall_sparsevecc(g_sparse(1)%L)
	return	! return the outE

  end subroutine Lmult	! Lmult

end module DataSens
