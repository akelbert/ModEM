! *****************************************************************************
module initFields
  ! Module containing the subroutines to initialize the field vectors.
  ! No memory allocation is done here. Vectors are assumed to already exist.

  use grid_orig
  use field_vectors
  use dimensions
  use ringcurrent
  use coreFwd
  use global
  use wrapper
  use sg_sparse_vector
  implicit none


Contains


  ! ***************************************************************************
  ! * initialize_fields is the routine to compute the initial values of the
  ! * magnetic fields based on the boundary conditions. Currently there is an
  ! * ambiguity that has to be dealt with, in that the _vertical_ outer boundary
  ! * information is not saved in the transformation Hx/Hy/Hz -> vectorh.
  ! * Hence the converse transformation vectorh -> Hx/Hy/Hz does not reconstruct
  ! * the original boundary values, but rather sets them to zero.
  ! * If Hx/Hy/Hz have not had the 'save' attribute, the boundary information in
  ! * them is soon lost. To avoid this problem we apply a patch by saving the
  ! * upper and lower boundary fields in bvH. 
  ! * Then Hx/Hy/Hz vectors should be identical to their initialized values after
  ! * the following two-line transformation:
  !	* call copyd1_d3_b(nx,ny,nz,Hx,Hy,Hz,hvec,x,y,z)
  ! * call insertBoundaryValues(bvH,Hx,Hy,Hz)

  subroutine initialize_fields_vec(hvec,bvec,bvH)


	complex(8), dimension(:), intent(inout)	:: hvec,bvec
	type (sparsevecc), intent(out)			:: bvH
	integer									:: istat
	
	!allocate(Hx(np1),Hy(np1),Hz(np1),STAT=istat)

	! Initialise the H field at the outer boundary, then make it
	! depend linearly on radius throughout the computational domain
	! such that all field components are zero at the CMB. The P10
	! source at the boundary has been modified to the more accurate
	! initialization by Ikuko Fujii allowing for the ringcurrent. 
	call set_initial_p10(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	! call set_initial_ringcurrent(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	call initBoundaryValues(bvH,Hx,Hy,Hz,grid)

	! save the boundary values of Hz
	! call insertBoundaryValues(bvH,Hx,Hy,Hz)

	! initialize vector <h> of A <h> = <b>
	call copyd3_d1_b(nx,ny,nz,Hx,Hy,Hz,hvec,x,y,z)

	! compute <b> of A <h> = <b>
	call calcb(nx,ny,nz,Hx,Hy,Hz,bvec,rho,x,y,z)

	!deallocate(Hx,Hy,Hz)

	!do i=1,bvH%ncoeff
	  !print *,'bvH: ',bvH%i(i),bvH%j(i),bvH%k(i),bvH%xyz(i),bvH%c(i)
	!end do
	!call create_cvector(grid, Bzero, EDGE)
	!call compare_fields(B, Bzero)


  end subroutine initialize_fields_vec



  ! ***************************************************************************
  ! * initialize_fields is the routine to compute the initial values of the
  ! * magnetic fields based on the boundary conditions
  subroutine initialize_fields(H,F)

	!use field_vectors

	type (cvector), intent(inout)			:: H,F
	type (sparsevecc)						:: Hb
	complex(8), dimension(:), allocatable	:: hvec
	!complex(8), dimension(:), allocatable   :: Hx,Hy,Hz
	integer									:: istat,i

	! Local allocations; perhaps not required if the variables become global
	allocate(hvec(np2),STAT=istat)
	!allocate(hx(np1),hy(np1),hz(np1),STAT=istat)

	!---- Initialise the H field at the outer boundary, then make it
	!     depend linearly with radius throughout the computational domain
	!     such that all field components are zero at the CMB. The P10
	!		source at the boundary has been modified to the more accurate
	!		initialization by Ikuko Fujii allowing for the ringcurrent. 
	call set_initial_p10(nx,ny,nz,x,y,z,Hx,Hy,Hz,np1)
	! NB --- This is now done in InitGlobalArrays
	! call set_initial_ringcurrent(nx,ny,nz,x,y,z,hx,hy,hz,np1)
	call copyd3_d1_b(nx,ny,nz,hx,hy,hz,hvec,x,y,z)	  ! extract interior components
	call initBoundaryValues(Hb,hx,hy,hz,grid)		  ! extract boundary components

	! Initialize the magnetic field
	call copyd1_d3_d(hvec,H,grid,bc=Hb)

	! Initialize the forcing with the boundary magnetic fields
	call create_cvector(grid,F,EDGE)
	call add_scvector(C_ONE,Hb,F)

	!deallocate(hx,hy,hz)
	deallocate(hvec)
  
  end subroutine initialize_fields




end module initFields
