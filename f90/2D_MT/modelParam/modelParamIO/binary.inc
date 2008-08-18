! Binary IO routines for modelParam (Gary Egbert's format)

      !****************************************************************** 
      subroutine write_modelParam_binary(fid,cfile,cond,grid)

      !  open cfile on unit fid, writes out object of
      !   type modelParam in standard format (readable by matlab
      !   routine readCond_2D.m), closes file

      integer, intent(in)		:: fid
      character*80, intent(in)	:: cfile
      type(modelParam_t), intent(in)		:: cond 
      type(grid2d_t), intent(in)            :: grid
      real (kind=selectedPrec), allocatable :: v(:,:)
      Integer                     ::k,j,Nz,Nza,NzEarth,Ny,stat
      character*80                :: paramType
      real (kind=selectedPrec)    :: AirCond

      Ny = grid%Ny
      NzEarth = grid%Nz-grid%Nza
      allocate(v(Ny,NzEarth))
 
      ! extract conductivity values from the modelParam structure
      call getValue_modelParam(cond,v,paramType,AirCond)
            
      open(unit=fid, file=cfile, form='unformatted')
      ! write grid dimensions ...
      write(fid) grid%Ny,grid%Nz,grid%Nza

      ! write grid spacings
      write(fid) grid%Dy
      write(fid) grid%Dz
      
      ! write model parameter
      write(fid) paramType
      write(fid) Ny,NzEarth
      write(fid) v
      write(fid) AirCond
      close(fid)

      deallocate(v)
       
      end subroutine write_modelParam_binary

     !******************************************************************
      subroutine read_modelParam_binary(fid,cfile,cond,grid)

      !  open cfile on unit fid, read in object of
      !   type modelParam in standard format 
      !   must either know or read in the grid
 
      integer, intent(in)		:: fid
      character(*), intent(in)	:: cfile
      type(modelParam_t), intent(inout)		:: cond 
      type(grid2d_t), intent(inout)         :: grid
      real (kind=selectedPrec), allocatable :: v(:,:)

      !  local variables
      integer 		:: Ny, Nz, Nza, NzEarth
      character*80	:: paramType
      real (kind=selectedPrec) :: AirCond

      open(unit=fid, file=cfile, form='unformatted',status='OLD')
      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      read(fid) Ny,Nz,Nza

      ! then allocate for grid
      call create_Grid2D(Ny,Nz,Nza,grid)
      
      !    read in grid spacings
      read(fid) grid%Dy
      read(fid) grid%Dz

      ! complete grid definition
      call gridCalcs(grid)

      read(fid) paramType
      read(fid) Ny,NzEarth

      if((grid%Ny .NE. Ny).OR.(grid%Nz - grid%Nza .NE. NzEarth)) then
         close(fid)
         call errStop('Grid size does not agree with model parameter size in readCond_2D')
      end if
      
      ! deallocate modelParam if it is allocated
      call deall_modelParam(cond)
 
      ! read info from file
      allocate(v(Ny,NzEarth)) 
      read(fid) v
      read(fid) AirCond
      close(fid)

      ! write into the modelParam structure
      call create_modelParam(grid,paramType,cond,v,AirCond)
      
      ! clean up
      deallocate(v)
      
      end subroutine read_modelParam_binary

      !******************************************************************
      subroutine writeVec_modelParam_binary(fid,cfile,nSigma,sigma,header)

      !  open cfile on unit fid, writes out nSigma objects of
      !   type modelParam , closes file

      integer, intent(in)		:: fid,nSigma
      character(*), intent(in)		:: cfile, header
      type(modelParam_t), intent(in)	:: sigma(nSigma) 

      integer i
      character(80) temp

      temp = header

      open(unit=fid, file=cfile, form='unformatted')
      write(fid) temp
      write(fid) nSigma
      do i = 1,nSigma
         write(fid) sigma(i)%paramType
         write(fid) sigma(i)%Ny,sigma(i)%NzEarth
         write(fid) sigma(i)%v
         write(fid) sigma(i)%AirCond
      enddo
      close(fid)
      end subroutine writeVec_modelParam_binary

     !******************************************************************
      subroutine readVec_modelParam_binary(fid,cfile,nSigma,sigma,header)

      !  open cfile on unit fid, read nSigma objects of
      !   type modelParam , closes file
      !  sigma(nsigma) must be allocated before calling

      integer, intent(in)		:: fid
      character(*), intent(in)		:: cfile
      integer, intent(in)		:: nSigma
      character(80), intent(out)		:: header
      type(modelParam_t), intent(inout) 	:: sigma(nsigma)

      ! local variables
      integer i, nS, Ny,NzEarth

      open(unit=fid, file=cfile, form='unformatted')
      read(fid) header
      read(fid) nS
      if(nS .NE. nSigma) then
          call errStop('size of sigma does not agree with contents of file in readAll_Cond2D')
      endif

      do i = 1,nSigma
         read(fid) sigma(i)%paramType
         read(fid) Ny,NzEarth
         if((sigma(i)%Ny .NE. Ny).OR. (sigma(i)%NzEarth .NE. NzEarth)) then
            close(fid)
            call errStop('Size of cond does not agree with contents of file in readAll_Cond2D')
         else
            read(fid) sigma(i)%v
            read(fid) sigma(i)%AirCond
         endif
      enddo
      close(fid)

      end subroutine readVec_modelParam_binary