! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: 2D MT

  use math_constants
  use file_units
  use dataspace
  use transmitters
  use receivers
  use datatypes

  implicit none

  private

  interface read_dataVectorMTX
	MODULE PROCEDURE read_Z
  end interface

  interface write_dataVectorMTX
	MODULE PROCEDURE write_Z
  end interface


  public     :: read_dataVectorMTX, write_dataVectorMTX


  character(100), private, save :: info_in_file
  integer,        private, save :: sign_in_file
  character(20),  private, save :: units_in_file


Contains

!**********************************************************************
! writes impedance file in ASCII format

   subroutine write_Z(allData,cfile)

      character(*), intent(in)					:: cfile
      type(dataVectorMTX_t), intent(in)			:: allData
      ! local variables
      type(dataBlock_t)                           :: data
      character(400)                            :: info
      integer									:: i,j,k,istat
      integer                                   :: nBlocks,iTx,iDt
      real(kind=prec)    						:: SI_factor
      logical                                   :: conjugate
      character(10)                             :: tab='           '


      info = 'Impedance responses from ModEM'
      SI_factor = ImpUnits('[V/m]/[T]',units_in_file)
      if (sign_in_file == ISIGN) then
        conjugate = .false.
      else if (abs(sign_in_file) == 1) then
        conjugate = .true.
      end if

      open(unit=ioDat,file=cfile,form='formatted',status='unknown')
      write(ioDat,'(a12)',advance='no') 'Description:'
      write(ioDat,*) trim(info)
      write(ioDat,'(a6)',advance='no') 'Units:'
      write(ioDat,*) trim(units_in_file)
      write(ioDat,'(a17,i3)') 'Sign convention: ',sign_in_file
      write(ioDat,*)

      nBlocks = countDataBlock(allData)
      write(ioDat,'(i5)') nBlocks

      ! loop over periods
      do iTx = 1,allData%nTx
         do iDt = 1,allData%d(iTx)%nDt

             data = allData%d(iTx)%data(iDt)

	         ! write period, number of sites and components for this period and data type
	         write(ioDat,'(es12.6,a5,i5)') txDict(iTx)%period,txDict(iTx)%mode,data%nSite
			 ! write latitude, longitude and elevation
			 do i = 1,2
			   do j = 1,data%nSite
	         	  write(ioDat,'(f14.3)',advance='no') rxDict(data%rx(j))%x(i)
			   enddo
	           write(ioDat,*)
			 enddo
			 write(ioDat,'(a8)',advance='no') ' '
			 ! write the descriptor
			 do k = 1,data%nComp
			 	write(ioDat,'(a15)',advance='no') trim(typeDict(data%dataType)%id(k))
			 enddo
			 write(ioDat,*)

	         ! conjugate data as necessary
	         if (conjugate) then
	         	do k=2,data%nComp,2
	         		data%value(k,:) = - data%value(k,:)
	         	end do
	         end if

	         ! write data
	         do j = 1,data%nSite
	         	! Note: temporarily, we write site id's according to their number;
	         	! in the future, they will be stored in the receiver dictionary
	         	write(ioDat,'(a10)',advance='no') trim(rxDict(data%rx(j))%id)
      	        !  note that each data field in a dataVec (i.e., allData%d(iTx)%data)
	         	!   is a real array of size (nComp,nSite)
	         	do k = 1,data%nComp
	         		write(ioDat,'(es15.6)',advance='no') data%value(k,j)*SI_factor
	         	enddo
	         	write(ioDat,*)
	         	write(ioDat,'(a10)',advance='no') tab
	         	do k = 1,data%nComp
	         	    if (data%errorBar) then
	         		   write(ioDat,'(es15.6)',advance='no') data%error(k,j)*SI_factor
	         		else
	         		   write(ioDat,'(es15.6)',advance='no') R_ZERO
	         		endif
	         	enddo
	         	write(ioDat,*)
	         enddo
         enddo
      enddo
      close(ioDat)
      call deall_dataBlock(data)

   end subroutine write_Z

!******************************************************************
! reads in the ASCII data file, sets up all dictionaries
! and the allData structure, including data and error bars

   subroutine read_Z(allData,cfile)

     character(*), intent(in)  				:: cfile
     type(dataVectorMTX_t), intent(inout)   	:: allData
     ! local variables
     type(dataBlock_t), dimension(:), pointer :: Data
     integer      							:: nBlocks,nTx,nDt,nSite,nComp,nData
     integer                                :: iTx,iRx,iDt,i,j,k,istat
     character(10), pointer,dimension(:)    :: siteIDs
     real(kind=prec), pointer, dimension(:,:) :: siteLoc
     real(kind=prec)                        :: SI_factor, Period
     logical      							:: conjugate, errorBar, isComplex
     character(400) 						:: temp,header
     character(2)                           :: Mode

      ! First, set up the data type dictionary, if it's not in existence yet
      call setup_typeDict()

      ! Now, read the data file
      open(unit=ioDat,file=cfile,form='formatted',status='old')
      read(ioDat,'(a13,a100)') temp,info_in_file
      read(ioDat,'(a7,a20)') temp,units_in_file

      SI_factor = ImpUnits(units_in_file,'[V/m]/[T]')

      read(ioDat,'(a17,i3)') temp,sign_in_file
      read(ioDat,*)

      if (sign_in_file == ISIGN) then
        conjugate = .false.
      else if (abs(sign_in_file) == 1) then
        conjugate = .true.
      else
        call errStop('Unknown sign convention in the data file '//cfile)
      end if

      ! Usually, nBlocks = nTx. However, in the rare case when there are more than
      ! one data type per frequency, this number can be bigger. We first read all
      ! blocks in the data file, then create the allData structure and dictionaries.
      read(ioDat,*) nBlocks

      allocate(Data(nBlocks),STAT=istat)
      nData = 0

      do i = 1,nBlocks
         ! read in period, number of components and sites for this dataVec
         read(ioDat,*) Period,Mode,nSite

         ! read in site locations
         allocate(siteLoc(nSite,2),siteIDs(nSite),STAT=istat)

         read(ioDat,*) (siteLoc(j,1),j=1,nSite)
         read(ioDat,*) (siteLoc(j,2),j=1,nSite)

         nComp = 2

         ! create dataVec object, read in data
         isComplex = .true.
         errorBar = .true.
         call create_dataBlock(nComp,nSite,Data(i),isComplex,errorBar)
         nData  = nData + nComp*nSite

         read(ioDat,'(a400)') header ! header line: need to parse this

         Data(i)%dataType = ImpType(Mode)

         do k=1,nSite
         	read(ioDat,*) siteIDs(k), (Data(i)%value(j,k),j=1,nComp)
         	read(ioDat,*)             (Data(i)%error(j,k),j=1,nComp)
         end do

         ! convert data to SI units
         Data(i)%value = SI_factor * Data(i)%value
         Data(i)%error = SI_factor * Data(i)%error

         ! conjugate data as necessary
         if (conjugate) then
         	do j=2,nComp,2
         		Data(i)%value(j,:) = - Data(i)%value(j,:)
         	end do
         end if

         ! Update the transmitter dictionary and the index (sets up if necessary)
         iTx = update_txDict(Period,Mode)
         Data(i)%tx = iTx

         ! Now, update the receiver dictionary and indices
         if (i .eq. 1) then
         	call setup_rxDict(nSite,siteLoc,siteIDs)
         	do k=1,nSite
         		Data(i)%rx(k) = k
         	end do
         else
         	do k=1,nSite
         		iRx = update_rxDict(siteLoc(k,:),siteIDs(k))
         		Data(i)%rx(k) = iRx
         	end do
         end if
         deallocate(siteLoc,siteIDs,STAT=istat)

	  end do

	  close(ioDat)

	  ! Finally, set up allData
	  nTx = size(txDict)
      call create_dataVectorMTX(nTx,allData)
      do iTx = 1,nTx
         ! count number of data types for this transmitter
         nDt = 0
         do i = 1,nBlocks
         	if (Data(i)%tx == iTx) then
         		nDt = nDt+1
         	end if
         end do
         ! initialize the data vectors for each
         call create_dataVector(nDt,allData%d(iTx))
         allData%d(iTx)%tx = iTx
         iDt = 1
         do i = 1,nBlocks
         	if (Data(i)%tx == iTx) then
         		allData%d(iTx)%data(iDt) = Data(i)
         		iDt = iDt+1
         	end if
         end do
         allData%d(iTx)%allocated = .true.
      end do
      allData%allocated = .true.

      ! Deallocate the temporary Data array
	  do i = 1,nBlocks
	     call deall_dataBlock(Data(i))
	  end do
	  deallocate(Data,STAT=istat)


  end subroutine read_Z


!**********************************************************************
  subroutine setError_dataBlock(err,d)
   !  whether error bars exist or not, sets new error bars
	 !  according to the input parameter err, which specifies the
	 !  fractional error in the output; any errors stored in d
	 !  are overwritten
	 ! TEMPORARY - need to be replaced with something data type specific

   real (kind=prec), intent(in) :: err
   type(dataBlock_t),intent(inout)             :: d

   !  local variables
   integer      			:: nComp, nSite

   if(.not.d%allocated) then
      call errStop('data vector not allocated in setError_dataVec')
   endif

   nComp = d%nComp
   nSite = d%nSite
   if (.not. d%errorBar) then
     allocate(d%error(nComp,nSite))
	 d%errorBar = .true.
   end if

   ! we should really be using dataType%isComplex information
	 ! create data errors for complex data, but we currently do not
	 ! have access to the dataType dictionary from here (I think);
	 ! so this is a temporary solution

	 d%error = err * abs(d%value)

  end subroutine setError_dataBlock

!**********************************************************************
  subroutine setError_dataVectorMTX(err,d)
   !  whether error bars exist or not, sets new error bars
	 !  according to the input parameter err, which specifies the
	 !  fractional error in the output; any errors stored in d
	 !  are overwritten

	 real (kind=prec), intent(in) :: err
   type(dataVectorMTX_t),intent(inout)          :: d

   !  local variables
   integer      			:: nTx, nData, j, k

   if(.not.d%allocated) then
      call errStop('data vector in setError_dataVectorMTX not allocated')
   endif

   do j = 1, d%nTx
     do k = 1, d%d(j)%nDt
        call setError_dataBlock(err,d%d(j)%data(k))
     enddo
   enddo

 end subroutine setError_dataVectorMTX


end module DataIO
