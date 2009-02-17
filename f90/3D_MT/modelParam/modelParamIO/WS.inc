! I/O routines for 3D_MT modelParam - extended Weerachai Siripunvaraporn's format

	!******************************************************************
	subroutine read_modelParam_WS(fid,cfile,m,paramType,grid)

      ! opens cfile on unit fid, writes out object of
      ! type modelParam in Weerachai Siripunvaraporn's format, closes file

      integer, intent(in)		           :: fid
      character(*), intent(in)             :: cfile
      type(modelParam_t), intent(out)	   :: m
      character(80), intent(out)           :: paramType
      type(grid3d_t), intent(inout)        :: grid
      ! local variables
      type(rscalar)                        :: rho, ccond
      character(80)                        :: someChar=''
      integer                              :: Nx, Ny, NzEarth, NzAir, someIndex=0, i, j, k, istat
      real(8)                              :: alpha=3.

	  paramType = ''

      ! Open file and read grid
      open(unit=fid, file=cfile, status='old')

      ! First read the comment line
      read(fid,'(a80)') someChar

	  ! Now read the second line with the grid dimensions
      read(fid,'(a80)') someChar
      read(someChar,*) Nx, Ny, NzEarth, someIndex

      if (someIndex /= 0) then
          write(0, *) 'Mapping not supported yet in read_modelParam_WS'
          stop
      end if

	  ! By default assume 'LINEAR RHO' - Weerachai's linear resistivity format
      if (index(someChar,'LOGE')>0) then
         paramType = 'LOGE'
      else
         paramType = 'LINEAR'
      end if

	  ! No information about the air layers in file. Hardcoded here.
	  nzAir = 6

	  ! Initialize and read in the grid
      call create_Grid3D(Nx,Ny,NzAir,NzEarth,grid)

      read(fid,*) (grid%dx(j),j=1,grid%nx)
      read(fid,*) (grid%dy(j),j=1,grid%ny)
      read(fid,*) (grid%dz(j),j=grid%nzAir+1,grid%nzAir+grid%nzEarth)

      !   Following is Kush's approach to setting air layers:
      ! mirror imaging the dz values in the air layer with respect to
      ! earth layer as far as we can using the following formulation
      ! air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom)
      i = grid%nzAir+1
      j = 0
      do k = grid%nzAir, 1, -1
          j = j + 1
          grid%dz(k) = ((alpha)**(j-1))*grid%dz(i)
          i = i + 1
      end do

      ! the topmost air layer has to be at least 30 km
      if (grid%dz(1).lt.30000) then
          grid%dz(1) = 30000
      end if

	  ! Create and read in the resistivity values
      call create_rscalar(grid,rho,CELL_EARTH)

 	  ! Read the resistivity
 	  do k=1,nzEarth
 	  	do j=1,Ny
 	  		do i=Nx,1,-1
 	  			read(fid,*,iostat=istat) rho%v(i,j,k)
 	  		end do
 	  	end do
 	  end do

 	  ! Note that our standard subroutine doesn't work with Weerachai's
 	  ! real value format. It is still better than either Mackie's or WS's...
      ! call read_rscalar(fid,rho)

      ! Convert from conductivity to resistivity
      call copy_rscalar(ccond,rho)
      if (index(paramType,trim(LOGE))>0) then
      	ccond%v = - rho%v
      else if (index(paramType,trim(LINEAR))>0) then
       	ccond%v = ONE/rho%v
 	  end if

      ! In case the grid origin is stored next (in metres!)...
      read(fid,*,iostat=istat) grid%ox, grid%oy, grid%oz

	  ! Defaults to the grid centre at the Earth's surface
      if (istat /= 0) then
          grid%ox = -sum(grid%dx)/2.0
          grid%oy = -sum(grid%dy)/2.0
          grid%oz = 0.0
      endif

      read(fid,*,iostat=istat) grid%rotdeg
      if (istat /= 0) then
          grid%rotdeg = 0.0
      end if

      close(fid)

	  ! Finally create the model parameter
	  call create_modelParam(grid,paramType,m,ccond)

	  ! Convert modelParam to logarithmic for computations
	  paramType = LOGE
	  call setType_modelParam(m,paramType)

	  ! now done with the rscalars, so deallocate
	  call deall_rscalar(rho)
	  call deall_rscalar(ccond)

	end subroutine read_modelParam_WS

	!******************************************************************
	subroutine write_modelParam_WS(fid,cfile,m,comment)

      ! opens cfile on unit fid, writes out object of
      ! type modelParam in Weerachai Siripunvaraporn's format, closes file

      integer, intent(in)		           :: fid
      character(*), intent(in)             :: cfile
      type(modelParam_t), intent(in)	   :: m
      character(*), intent(in), optional   :: comment
      ! local variables
      type(grid3d_t)                       :: grid
      type(rscalar)                        :: rho,ccond
      character(80)                        :: paramType
      integer                              :: Nx, Ny, NzEarth, i, j, k, istat

	  ! Convert modelParam to linear for output
	  paramType = LINEAR
	  call getValue_modelParam(m,paramType,ccond)

	  grid = ccond%grid

      open(unit=fid, file=cfile, status='unknown')

	  if (present(comment)) then
	  	write(fid,*) '# ',trim(comment)
	  else
      	write(fid,*) '# Model written by ModEM in WS format'
      end if

      ! Write grid geometry definitions
      Nx=grid%nx
      Ny=grid%ny
      NzEarth=grid%nz - grid%nzAir

      write(fid,'(4i5)',advance='no') Nx,Ny,NzEarth,0
      write(fid,*) trim(paramType)

      ! Write grid spacings
      do j=1,grid%nx
      	write(fid,'(f12.3)',advance='no') grid%dx(j)
	  end do
	  write(fid,*)
	  do j=1,grid%ny
      	write(fid,'(f12.3)',advance='no') grid%dy(j)
	  end do
	  write(fid,*)
	  do j=grid%nzAir+1,grid%nz
      	write(fid,'(f12.3)',advance='no') grid%dz(j)
	  end do
	  write(fid,*)

      ! Convert conductivity to resistivity
      call copy_rscalar(rho,ccond)
      if (index(paramType,trim(LOGE))>0) then
      	rho%v = - ccond%v
      else if (index(paramType,trim(LINEAR))>0) then
       	rho%v = ONE/ccond%v
 	  end if

 	  ! Write the resistivity
 	  write(fid,*)
 	  do k=1,nzEarth
 	  	do j=1,Ny
 	  		do i=Nx,1,-1
 	  			write(fid,'(es13.5)',iostat=istat,advance='no') rho%v(i,j,k)
 	  		end do
 	  		write(fid,*)
 	  	end do
 	  	write(fid,*)
 	  end do

 	  ! Note that our standard subroutine doesn't work with Weerachai's
 	  ! real value format. It is still better than either Mackie's or WS's...
	  ! call write_rscalar(fid,rho)

      ! Also write the grid origin (in metres!) and rotation (in degrees)...
      write(fid,'(3f9.3)',iostat=istat) grid%ox, grid%oy, grid%oz
      write(fid,'(f9.3)',iostat=istat)  grid%rotdeg
      close(fid)

	  ! Now done with ccond, rho and grid, so deallocate
	  call deall_rscalar(rho)
	  call deall_rscalar(ccond)
	  call deall_grid3d(grid)

      end subroutine write_modelParam_WS
