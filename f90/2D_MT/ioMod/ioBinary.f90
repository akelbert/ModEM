! *****************************************************************************
!  Module for basic input and output of standard data structures
!  for 2D MT modeling and inversion code 
module ioBinary
   use math_constants
   use emfield
   use dataspace
   use emsolver
   use utilities
   implicit none

   !  routines that are public
   public	::  read_cvector,write_cvector,write_EMsolnMTX, &
		read_Grid2D, write_Z, read_Z, &
		read_cond2d, write_cond2d

   Contains
     !******************************************************************
      subroutine read_cvector(fid,cfile,vec)

      !  open cfile on unit fid, read in object of
      !   type cvector in standard format 
      !   vec must be allocated before calling
 

      integer, intent(in)		:: fid
      character*80, intent(in)	:: cfile
      character*80		:: msg
      type(cvector), intent(inout)		:: vec 

      !  local variables
      integer 		:: N1, N2
      character*80	:: gridType


      if(vec%allocated) then
         open(unit=fid, file=cfile, form='unformatted')
         read(fid) gridType
         read(fid) N1,N2
         if((vec%N1 .NE. N1).OR.(vec%N2 .NE. N2)) then
            close(fid)
            msg = 'Size of vec does not agree with contents of file'
            call errStop(msg)
         else
            vec%gridType = gridType
            read(fid) vec%v
            close(fid)
         endif
      else
         msg = 'vec must be allocated before call to read_cvector'
         call errStop(msg)
      endif
      end subroutine read_cvector

     !******************************************************************
      subroutine openR_cvector(fid,cfile,header,nVec)

      !  open cfile on unit fid for reading multiple cvector objects
      !   reads header, number of vectors in file

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      character*80, intent(out)		:: header
      integer, intent(out)		:: nVec

      open(unit=fid, file=cfile, form='unformatted',status = 'old')
      read(fid) header
      read(fid) nVec
      end subroutine openR_cvector

     !******************************************************************
      subroutine cvectorRead(fid,vec)

      !  reads one object of type cvector from unit fid 
      !  by default, just reads next object ... might modify to
      !   allow skipping to read an arbitrary record number

      integer, intent(in)		:: fid
      type(cvector), intent(inout)	:: vec 
   
      ! local variables:
      character*80	:: gridType, msg
      integer		:: N1,N2
 
      if(vec%allocated) then
         read(fid) gridType
         read(fid) N1,N2
         if((vec%N1 .NE. N1).OR.(vec%N2 .NE. N2)) then
            msg = 'Size of vec does not agree with contents of file'
            call errStop(msg)
         else
            vec%gridType = gridType
            read(fid) vec%v
         endif
      else
         msg = 'vec must be allocated before call to cvectorRead'
         call errStop(msg)
      endif

      end subroutine cvectorRead

     !******************************************************************
      subroutine openW_cvector(fid,cfile,header,nVec)

      !  open cfile on unit fid for writing multiple cvector objects
      !   writes header, number of vectors to be put in file

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      character*80, intent(in)		:: header
      integer, intent(in)		:: nVec

      open(unit=fid, file=cfile, form='unformatted',status='unknown')
      write(fid) header
      write(fid) nVec
      end subroutine openW_cvector

     !******************************************************************
      subroutine cvectorWrite(fid,vec)

      !  writes one object of type cvector to unit fid, already opened

      integer, intent(in)		:: fid
      type(cvector), intent(in) 		:: vec 
 
      write(fid) vec%gridType
      write(fid) vec%N1,vec%N2
      write(fid) vec%v

      end subroutine cvectorWrite

     !******************************************************************
      subroutine write_cvector(fid,cfile,vec)

      !  open cfile on unit fid, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      type(cvector), intent(in)		:: vec 

      open(unit=fid, file=cfile, form='unformatted',status='unknown')
      write(fid) vec%gridType
      write(fid) vec%N1,vec%N2
      write(fid) vec%v
      close(fid)
      end subroutine write_cvector

     !******************************************************************
      subroutine write_EMsolnMTX(fid,cfile,eAll)

      !  open cfile on unit fid, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file
      !  NOT coded at present to specifically write out TE/TM
      !    solutions, periods, etc. (can get this infor from
      !    eAll%solns(j)%tx, but only with access to TXdict.

      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      type(EMsolnMTX), intent(in)		:: eAll 
  
      integer		:: j
    
      open(unit=fid, file=cfile, form='unformatted',status='unknown')
      
      write(fid) eAll%nTx
      do j = 1,eAll%nTx
         write(fid) eAll%solns(j)%vec%gridType
         write(fid) eAll%solns(j)%vec%N1,eAll%solns(j)%vec%N2
         write(fid) eAll%solns(j)%vec%v
      enddo
      close(fid)
      return
      end subroutine write_EMsolnMTX

     !**********************************************************************
      subroutine read_Grid2D(fid,cfile,grid)
     !  reads in basic grid, allocating for Dy, Dz 
      integer, intent(in)		:: fid
      character*80, intent(in)		:: cfile
      type(grid2d_t), intent(inout)	:: grid

      ! local variables:
      integer				:: Ny, Nz, Nza
      
      open(unit=fid,file=cfile,form='unformatted',status='OLD')
      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      read(fid) Ny,Nz,Nza

      ! then allocate for grid
      call create_Grid2D(Ny,Nz,Nza,grid)

      !    read in grid spacings
      read(fid) grid%Dy
      read(fid) grid%Dz
      close(fid)
      end subroutine read_Grid2D


     !**********************************************************************
      subroutine write_Z(fid,cfile,nTx,periods,modes,nSites,sites,allData)
     ! writes impedance file, including list of periods, siteLocations
     !   NOTE: this assumes that the arrays "sites" and "periods" are
     !    essentially identical to the receiver and transmitter dictionaries
     !    (in which case, why have both these arrays and the dicts?)

      integer, intent(in)       :: fid
      character(*), intent(in)  :: cfile
      integer, intent(in)       :: nTx,nSites
      real(kind =8),intent(in)   :: periods(nTx)
      character*2, intent(in)	:: modes(nTx)
      real(kind = 8),intent(in)   :: sites(2,nSites)
      type(dvecMTX), intent(in)      :: allData
      real(kind = 8), dimension(:,:), allocatable :: siteTemp

     ! local variables
      integer   :: ns,iTx,k,j

      open(unit=fid,file=cfile,form='unformatted',status='unknown')
      write(fid) allData%nTx
      ! loop over periods
      do iTx = 1,allData%nTx
         ns = allData%d(iTx)%nSite
         ! write period, number of sites for this period
         write(fid) periods(iTx),modes(iTx),ns
         ! write site locations for this period
         allocate(siteTemp(2,ns))
         do k = 1,ns
            do j = 1,2
               siteTemp(j,k) = sites(j,allData%d(iTx)%rx(k))
            enddo
         enddo
         write(fid) siteTemp
         write(fid) allData%d(iTx)%data
         write(fid) allData%d(iTx)%err
         deallocate(siteTemp)
      enddo
      close(fid)
      end subroutine write_Z

     !**********************************************************************
      subroutine read_Z(fid,cfile,nTx,periods,modes,nSites,sites,allData)
     ! reads in data file, returns list of periods, modes, siteLocations, and
     !   sets up data vector structure, including data and error bars
      integer, intent(in)       :: fid
      character(*), intent(in)  :: cfile
      integer, intent(out)      :: nTx,nSites
      real(kind = 8),dimension(:), pointer     :: periods
      real(kind = 8),dimension(:,:), pointer   :: sites
      character*2, dimension(:), pointer	:: modes
      type(dvecMTX), intent(inout)   :: allData

     ! local variables
      integer   :: nComp = 2
      integer   :: ns,iTx,k,l,j,Ndata
      real(kind = 8), allocatable, dimension(:,:) :: siteTemp,siteTempAll
      logical   :: newSite

      open(unit=fid,file=cfile,form='unformatted',status='old')
      read(fid) nTx
      allocate(periods(nTx))
      allocate(modes(nTx))
      allocate(allData%d(nTx))
      allData%allocated = .true.
      allData%errorBar = .true.
      allData%nTx = nTx
     ! loop over transmitters (periods/modes)
      Ndata = 0
      do iTx = 1,nTx
         ! read in number of sites for this period
         read(fid) periods(iTx),modes(iTx),ns
         ! read in site locations
         allocate(siteTemp(2,ns))
         read(fid) siteTemp
         ! create dvec object, read in data
         allData%d(iTx)%errorBar = .true.
         call create_Dvec(nComp,ns,allData%d(iTx))
         Ndata  = Ndata + nComp*ns
         allData%d(iTx)%tx = iTx
	 if(modes(iTx) .eq. 'TM') then
	    allData%d(iTx)%datatype = 2
	 else
	    allData%d(iTx)%datatype = 1
	 endif
         read(fid) allData%d(iTx)%data
         read(fid) allData%d(iTx)%err
         if(iTx .eq. 1) then
           ! allocate temporary storage for full sites list
           ! (this might not always work ... I am assuming that the
           !        max number of sites is ns from period one * nTx)
            allocate(siteTempAll(2,ns*nTx))
            nSites = ns
            do k = 1,ns
               siteTempAll(1,k) = siteTemp(1,k)
               siteTempAll(2,k) = siteTemp(2,k)
               allData%d(iTx)%rx(k) = k
            enddo
         else
            ! check to see if site locations are already in list
            !  if not add to list; in any event set reciever "pointer" rx
            do k = 1,ns
               newSite = .true.
               do l = 1,nSites
                  if((siteTemp(1,k).eq.siteTempAll(1,l)).and.  &
		  	(siteTemp(2,k).eq.siteTempAll(2,l))) then
                     newSite = .false.
                     allData%d(iTx)%rx(k) = l
                     exit
                  endif
               enddo
               if(newSite) then
                  nSites = nSites+1
                  siteTempAll(1,nSites) = siteTemp(1,k)
                  siteTempAll(2,nSites) = siteTemp(2,k)
                  allData%d(iTx)%rx(k) = nSites
               endif
            enddo
         endif
         deallocate(siteTemp)
      enddo
      allData%Ndata = Ndata
      ! copy list of unique sites into "sites" array
      allocate(sites(2,nSites))
      do k = 1,nSites
         do j = 1,2
            sites(j,k) = siteTempAll(j,k)
         enddo
      enddo
      close(fid)
      end subroutine read_Z

    !******************************************************************
      subroutine write_Cond2D(fid,cfile,cond,grid)

      !  open cfile on unit fid, writes out object of
      !   type modelParam in standard format (readable by matlab
      !   routine readCond2D.m), closes file

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
      write(fid) paramType
      write(fid) Ny,NzEarth
      write(fid) v
      write(fid) AirCond
      close(fid)

      deallocate(v)
       
      end subroutine write_Cond2D

     !******************************************************************
      subroutine read_Cond2D(fid,cfile,cond,grid)

      !  open cfile on unit fid, read in object of
      !   type modelParam in standard format 
      !   must either know or read in the grid
 
      integer, intent(in)		:: fid
      character(*), intent(in)	:: cfile
      type(modelParam_t), intent(inout)		:: cond 
      type(grid2d_t), intent(inout)         :: grid
      real (kind=selectedPrec), allocatable :: v(:,:)

      !  local variables
      integer 		:: Ny, NzEarth
      character*80	:: paramType
      real (kind=selectedPrec) :: AirCond

      open(unit=fid, file=cfile, form='unformatted',status='OLD')
      read(fid) paramType
      read(fid) Ny,NzEarth

      if((grid%Ny .NE. Ny).OR.(grid%Nz - grid%Nza .NE. NzEarth)) then
         close(fid)
         call errStop('Grid size does not agree with contents of file in read_cond2d')
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
      
      end subroutine read_Cond2D

end module ioBinary
