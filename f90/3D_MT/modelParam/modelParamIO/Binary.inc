! Binary IO routines for 3D_MT modelParam (Gary Egbert's format)

  !******************************************************************
   subroutine write_modelParam_binary(fid,cfile,m)

   !  open cfile on unit fid, writes out object of
   !   type modelParam in standard *binary* format (comparable
   !   format written by write_Cond3D), then close file
   ! NOTE: uses grid for the model size; needs to be thought
   !       about if we decide to decouple grid and the model

      integer, intent(in)               :: fid
      character(*), intent(in)  :: cfile
      type(modelParam_t), intent(in)              :: m

      type(rscalar)           :: ccond
      type(grid3d_t)          :: grid
      real(kind=prec) :: AirCond
      character(80)           :: paramType

      paramType = LOGE

      call getValue_modelParam(m,paramType,ccond,AirCond)

      grid = ccond%grid

      open(unit=fid, file=cfile, form='unformatted')
      write(fid) grid%nx,grid%ny,grid%nzEarth,grid%nzAir
      write(fid) grid%dx
      write(fid) grid%dy
      write(fid) grid%dz(grid%NzAir+1:grid%Nz)
      write(fid) paramType
      write(fid) AirCond
      write(fid) ccond%v
      close(fid)

      call deall_rscalar(ccond)
      call deall_grid3d(grid)

      end subroutine write_modelParam_binary

  !******************************************************************
   subroutine read_modelParam_binary(fid,cfile,m,paramType,grid)

   !  open cfile on unit fid, writes out object of
   !   type modelParam in standard *binary* format (comparable
   !   format written by write_Cond3D), then close file

      integer, intent(in)                  :: fid
      character(*), intent(in) 		       :: cfile
      type(modelParam_t), intent(out)      :: m
      character(80), intent(in)            :: paramType
      type(grid3d_t), intent(inout), optional :: grid

      integer		:: NzAir,Nz,Nx,Ny,NzEarth
      real(kind=prec) 	:: AirCond
      type(rscalar)             :: ccond
      character(80)             :: gridType,rParamType

      open(unit=fid, file=cfile, form='unformatted',status='old')

      read(fid) Nx,Ny,NzEarth,NzAir

      call create_grid3d(Nx,Ny,NzAir,NzEarth,grid)
      read(fid)  grid%dx
      read(fid)  grid%dy
      read(fid)  grid%dz(NzAir+1:NzAir+NzEarth)
      call gridCalcs(grid)

      read(fid) rParamType

      gridType = CELL_EARTH
      call create_rscalar(grid,ccond,gridType)

      read(fid) AirCond
      read(fid) ccond%v

      close(fid)

      ! move cell conductivities read into rscalar object into a modelParam
	  ! object ... this dance needed to keep modelParam attributes private
	  !   First need to create model parameter
	  call create_modelParam(grid,rParamType,m,ccond,AirCond)

	  ! convert modelParam to the required paramType
	  call setType_modelParam(m,paramType)

	  ! now done with ccond, so deallocate
	  call deall_rscalar(ccond)

      end subroutine read_modelParam_binary

     !******************************************************************
      subroutine writeVec_modelParam_binary(fid,cfile,nSigma,m,header)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file

      integer, intent(in)               :: fid,nSigma
      character(*), intent(in)          :: cfile, header
      type(modelParam_t), intent(in)      :: m(nSigma)

      integer 		:: i,Nz,NzAir

      open(unit=fid, file=cfile, form='unformatted')
      write(fid) header
      write(fid) nSigma

      do i = 1,nSigma
         Nz = m(i)%grid%Nz
         NzAir = m(i)%grid%Nz-m(i)%NzEarth
         write(fid) m(i)%paramType
         write(fid) m(i)%Nx,m(i)%Ny,m(i)%NzEarth
         write(fid) m(i)%grid%dx
         write(fid) m(i)%grid%dy
         write(fid) m(i)%grid%dz(NzAir+1:Nz)
         write(fid) m(i)%AirCond
         write(fid) m(i)%CellCond%v
      enddo
      close(fid)
      end subroutine writeVec_modelParam_binary

     !******************************************************************
      subroutine readVec_modelParam_binary(fid,cfile,nSigma,m,header)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file

      integer, intent(in)               :: fid
      integer, intent(in)              :: nSigma
      character(*), intent(in)          :: cfile
      character(80), intent(out)        :: header
      type(modelParam_t), intent(inout)   :: m(nSigma)

      integer 		:: i,nS,Nx,Ny,Nz,NzAir,NzEarth,istat

      open(unit=fid, file=cfile, form='unformatted')
      read(fid) header
      read(fid) nS

      if(nS .NE. nSigma) then
          call errStop('size of sigma does not agree with contents of file in readAll_Cond2D')
      endif

      do i = 1,nSigma
         read(fid) m(i)%paramType
         read(fid) Nx,Ny,NzEarth
         if((m(i)%Nx .NE. Nx) .OR. (m(i)%Ny .NE. Ny) .OR. (m(i)%NzEarth .NE. NzEarth)) then
            close(fid)
            call errStop('Size of cond does not agree with contents of file in readAll_Cond2D')
         else
            Nz = m(i)%grid%Nz
            NzAir = m(i)%grid%Nz-m(i)%NzEarth
            read(fid) m(i)%grid%dx
            read(fid) m(i)%grid%dy
            read(fid) m(i)%grid%dz(NzAir+1:Nz)
            read(fid) m(i)%AirCond
            read(fid) m(i)%CellCond%v
         endif
      enddo
      close(fid)
      end subroutine readVec_modelParam_binary
