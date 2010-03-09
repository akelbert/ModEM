! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: 3D MT

  use math_constants
  use file_units
  use dataspace
  use transmitters
  use receivers
  use datatypes
  use utilities

  implicit none

  private

  interface read_dataVecMTX
	MODULE PROCEDURE read_Z
  end interface
  
  interface read_dataVecMTX_list_format
  	MODULE PROCEDURE read_Z_list_format
  end interface
  		
  		
  interface write_dataVecMTX
	MODULE PROCEDURE write_Z
  end interface


  public     :: read_dataVecMTX, write_dataVecMTX,read_dataVecMTX_list_format


  character(100), private, save :: info_in_file
  integer,        private, save :: sign_in_file
  character(20),  private, save :: units_in_file


Contains
!**********************************************************************   
! reads in the ASCII data file (List format), sets up all dictionaries
! and the allData structure, including data and error bars
   subroutine read_Z_list_format(allData,cfile)
 
     character(*), intent(in)  				:: cfile
     type(dataVecMTX_t), intent(inout)   	:: allData
     ! local variables
     integer      							:: nTx,nRx,nDt,nComp,iComp,dataType,nSite,max_nDt,max_nComp
     integer                                :: iTx,iRx,iDt,i,j,k,istat,iSite
     real(kind=prec), dimension(3) 			:: siteLoc_local,siteLoc_ref
     real(kind=prec)						:: SI_factor
     real(kind=prec)                        :: Period
     logical      							:: conjugate, errorBar, isComplex
     character(400)                         :: Data_Type_Key_word
     character(400)                         :: temp,line,stn_id_local,stn_id_ref
     logical								:: new_site,new_period
     integer 								:: stat,string_position
     real(kind=prec), pointer, dimension(:)	:: temp_data,temp_error
     character(100)                         :: args(50)
	 Integer, pointer, dimension(:,:)		:: stn_vec
     Integer, pointer, dimension(:)			::nSite_vec,data_type_vec	 
     integer								:: nargs,ios,per_counter,data_type
   
     real(kind=prec), pointer, dimension(:,:,:,:) :: global_data,global_data_temp    
	 real(kind=prec), pointer, dimension(:,:,:,:) :: global_error,global_error_temp 
	 type(dataVecMTX_t)                    :: allData_temp   

     
          ! First, set up the data type dictionary, if it's not in existence yet
      call setup_typeDict()
	  max_nDt=size(typeDict)
	  max_nComp=Maxval(typeDict(:)%nComp)
      allocate(temp_data(max_nComp),temp_error(max_nComp))
      ! Now, read the data file
      open(unit=ioDat,file=cfile,status='old')
         
      read(ioDat,'(a13,a100)') temp,info_in_file
      read(ioDat,'(a17,i3)') temp,sign_in_file
      write(6,*) 'sign_in_file ',sign_in_file
       if (sign_in_file == ISIGN) then
        conjugate = .false.
      else if (abs(sign_in_file) == 1) then
        conjugate = .true.
      else
        call errStop('Unknown sign convention in the data file '//cfile)
      end if

Reading_all_file_loop: &
   do
             read(ioDat,'(A)', iostat=stat) line
              if(stat /= 0) exit Reading_all_file_loop     !End of the file
              
              line=trim(line) 
		      call COMPACT(line)                           !Removes spaces between arrguments of the line
		      call parse(line,' ',args,nargs)              !Counts and finds out the arrguments and number of arrguemnts in the line
		      
		      if (nargs==0) then 						   ! If its an empty line---> ignor it
		      !write(6,*) '###################'
		        per_counter=0 
		        cycle
		      end if
              !Check if it is a header line
              string_position=findstr(trim(args(1)),'#')
              if (string_position .eq. 1 )then           !If the first character in the first argument is '#', ---> Header line (or comments line),ignor it
              !write(6,*) '###################'
                per_counter=0
                cycle 
              end if
              ! It is a line that contains data
                      per_counter=per_counter+1
 		              line=trim(line) 
				      call COMPACT(line)                           !Removes spaces between arrguments of the line
				      call parse(line,' ',args,nargs)              !Counts and finds out the arrguments and number of arrguemnts in the line             
                            Data_Type_Key_word=args(1)
					        read(args(2),*)Period
					        stn_id_local=args(3)
					        read(args(4:6),*)siteLoc_local(:)
					        
  					        if (Data_Type_Key_word =='Full_Interstation_TF') then
  					            data_type = Full_Interstation_TF
					            nComp=typeDict(data_type)%nComp
					            stn_id_ref=args(7)
					            read(args(8:10),*)siteLoc_ref(:)
					            units_in_file=args(11)
						        read(args(12:19),*)temp_data(1:8)
						        read(args(20:23),*)(temp_error(j),j=1,8,2)
						        do j=2,8,2
						         temp_error(j)= temp_error(j-1)	
						        end do
						        call  update_rxDict(siteLoc_local,stn_id_local,iRx,new_site,siteLoc_ref,stn_id_ref)  					        
					        elseif (Data_Type_Key_word =='Full_Impedance') then
							    data_type = Full_Impedance	
					            nComp=typeDict(data_type)%nComp
					            units_in_file=args(7)
						        read(args(8:15),*)temp_data(1:8)
						        read(args(16:19),*)(temp_error(j),j=1,8,2)
						        do j=2,8,2
						         temp_error(j)= temp_error(j-1)	
						        end do
						         call  update_rxDict(siteLoc_local,stn_id_local,iRx,new_site)       
					        elseif (Data_Type_Key_word =='Full_Vertical_Components') then
						        data_type = Full_Vertical_Components
						        nComp=typeDict(data_type)%nComp
					            units_in_file=args(7)
						        read(args(8:11),*)temp_data(1:4)
						        read(args(12:13),*)(temp_error(j),j=1,4,2)
						        do j=2,4,2
						         temp_error(j)= temp_error(j-1)	
						        end do  
						         call  update_rxDict(siteLoc_local,stn_id_local,iRx,new_site) 
					        end if
                            call  update_txDict(Period,2,iTx,new_period)
                            							
                             if(new_period .or. new_site ) then
							     nTx=size(global_data,1)
								 nRx=size(global_data,2)
						         allocate(global_data_temp(nTx,nRx,max_nDt,max_nComp))
								 allocate(global_error_temp(nTx,nRx,max_nDt,max_nComp))
								 global_data_temp=0.0
								 global_data_temp=global_data
								 global_error_temp=global_error

							     if (associated (global_data))  deallocate(global_data)
								 if (associated (global_error))  deallocate(global_error)
							     nTx=size(txDict)
								 nRx=size(rxDict)								 
                                 allocate(global_data(nTx,nRx,max_nDt,max_nComp))
								 allocate(global_error(nTx,nRx,max_nDt,max_nComp))
								 global_data=0.0
								 do i=1,size(global_data_temp,1)
								   do j=1,size(global_data_temp,2)
								      do k=1,max_nDt
									    global_data(i,j,k,:)=global_data_temp(i,j,k,:)
										global_error(i,j,k,:)=global_error_temp(i,j,k,:)
									  end do
                                   end do
								 end do
								 deallocate(global_data_temp,global_error_temp)
                             end if
							         if (conjugate) then
										do j=2,nComp,2
											temp_data(j)= - temp_data(j)
										end do
									 end if							 
                            call get_nComp_DT(trim(Data_Type_Key_word),dataType,nComp)
							SI_factor = ImpUnits(units_in_file,typeDict(dataType)%units(1)) !Assuming there is NO mixing in the data types: EACH DATA TYPE HAS ITS OWN UNITS
                            global_data(iTx,iRx,dataType,1:8)  = SI_factor*temp_data(1:8)
							global_error(iTx,iRx,dataType,1:8) = SI_factor*temp_error(1:8)
                            if(stat /= 0) exit Reading_all_file_loop     !End of the file

  end do Reading_all_file_loop
    close(ioDat)
	
	nTx =size(txDict)
	nRx= size(rxDict)

    
	 allocate(nSite_vec(max_nDt))
     allocate(data_type_vec(max_nDt))
	 allocate(stn_vec(nRx,max_nDt))
	 
	call create_dataVecMTX(nTx,allData) 
  do iTx=1,nTx
        nDt=0 
     do iDt=1,max_nDt
        nSite_vec(iDt)=0
  		do iRx=1,nRx
 			  if( global_data(iTx,iRx,iDt,1) .ne. 0  ) then
 			   nSite_vec(iDt)=nSite_vec(iDt)+1
			   stn_vec(nSite_vec(iDt),iDt)=iRx
 			  end if 
 		end do
			  if (nSite_vec(iDt) .ne. 0) then
			  nDt=nDt+1
              data_type_vec(nDt)=iDt
			  end if
 		end do
          call create_dataVecTX(nDt,allData%d(iTx))
		  do iDt=1,nDt
		  	isComplex = .true.
			errorBar = .true.
			nComp=typeDict(data_type_vec(iDt))%nComp
			nSite=nSite_vec(data_type_vec(iDt))
		    call create_dataVec(nComp,nSite,allData%d(iTx)%data(iDt),isComplex,errorBar)			
			     do iSite=1,nSite
				    iRx=stn_vec(isite,data_type_vec(iDt))
				    do icomp=1,nComp
				      allData%d(iTx)%data(iDT)%value(icomp,iSite) =  global_data(iTx,iRx,data_type_vec(iDt),icomp)
				      allData%d(iTx)%data(iDt)%error(icomp,iSite) =  global_error(iTx,iRx,data_type_vec(iDt),icomp)
                    end do				   
				   allData%d(iTx)%data(iDt)%rx(iSite) = iRx   
				 end do  
				 allData%d(iTx)%data(iDt)%dataType=data_type_vec(iDt)
				 allData%d(iTx)%data(iDt)%tx=iTx
				
             end do
			 allData%d(iTx)%allocated = .true.	
			 allData%d(iTx)%tx=iTx	 
			 allData%d(iTx)%nDt=nDt
  end do
  allData%allocated = .true.
  
  ! this here just to test the output with the old format
  units_in_file= '[V/m]/[A/m]'
                         
     
   end subroutine read_Z_list_format     
!**********************************************************************
! writes impedance file in ASCII format (old format)

   subroutine write_Z(allData,cfile)

      character(*), intent(in)					:: cfile
      type(dataVecMTX_t), intent(in)			:: allData
      ! local variables
      type(dataVec_t)                           :: data
      character(400)                            :: info
      integer									:: i,j,k,istat
      integer                                   :: nBlocks,iTx,iDt
      real(kind=prec), pointer, dimension(:)    :: SI_factor
      logical                                   :: conjugate
      character(10)                             :: tab='           '


      info = 'Impedance responses from ModEM'
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

      nBlocks = countDataVec(allData)
      write(ioDat,'(i5)') nBlocks

      ! loop over periods
      do iTx = 1,allData%nTx
         do iDt = 1,allData%d(iTx)%nDt

             data = allData%d(iTx)%data(iDt)

	         ! write period, number of sites and components for this period and data type
	         write(ioDat,'(es12.6,2i5)') txDict(iTx)%period,data%nComp,data%nSite
			 ! write latitude, longitude and elevation
			 do i = 1,3
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

             ! compute the units conversion factors
             allocate(SI_factor(data%nComp),STAT=istat)
             do k=1,data%nComp
                SI_factor(k) = ImpUnits(typeDict(data%dataType)%units(k),units_in_file)
             end do

	         ! write data
	         do j = 1,data%nSite
	         	! site id's are stored in the receiver dictionary
	         	write(ioDat,'(a10)',advance='no') trim(rxDict(data%rx(j))%id)
      	        !  note that each data field in a dataVec (i.e., allData%d(iTx)%data)
	         	!   is a real array of size (nComp,nSite)
	         	do k = 1,data%nComp
	         		write(ioDat,'(es15.6)',advance='no') data%value(k,j)*SI_factor(k)
	         	enddo
	         	write(ioDat,*)
	         	write(ioDat,'(a10)',advance='no') tab
	         	do k = 1,data%nComp
	         	    if (data%errorBar) then
	         		   write(ioDat,'(es15.6)',advance='no') data%error(k,j)*SI_factor(k)
	         		else
	         		   write(ioDat,'(es15.6)',advance='no') R_ZERO
	         		endif
	         	enddo
	         	write(ioDat,*)
	         enddo
	         deallocate(SI_factor,STAT=istat)
         enddo
      enddo
      close(ioDat)
      call deall_dataVec(data)

   end subroutine write_Z

!******************************************************************
! reads in the ASCII data file, sets up all dictionaries
! and the allData structure, including data and error bars

   subroutine read_Z(allData,cfile)

     character(*), intent(in)  				:: cfile
     type(dataVecMTX_t), intent(inout)   	:: allData
     ! local variables
     type(dataVec_t), dimension(:), pointer :: Data
     integer      							:: nBlocks,nTx,nDt,nSite,nComp,nData
     integer                                :: iTx,iRx,iDt,i,j,k,istat
     character(10), pointer,dimension(:)    :: siteIDs
     real(kind=prec), pointer, dimension(:,:) :: siteLoc
     real(kind=prec), pointer, dimension(:) :: SI_factor
     real(kind=prec)                        :: Period
     logical      							:: conjugate, errorBar, isComplex
     character(400) 						:: temp,header

      ! First, set up the data type dictionary, if it's not in existence yet
      call setup_typeDict()

      ! Now, read the data file
      open(unit=ioDat,file=cfile,form='formatted',status='old')
      read(ioDat,'(a13,a100)') temp,info_in_file
      read(ioDat,'(a7,a20)') temp,units_in_file
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
         read(ioDat,*) Period,nComp,nSite

         ! read in site locations
         allocate(siteLoc(nSite,3),siteIDs(nSite),STAT=istat)

         read(ioDat,*) (siteLoc(j,1),j=1,nSite)
         read(ioDat,*) (siteLoc(j,2),j=1,nSite)
         read(ioDat,*) (siteLoc(j,3),j=1,nSite)

         ! create dataVec object, read in data
         isComplex = .true.
         errorBar = .true.
         call create_dataVec(nComp,nSite,Data(i),isComplex,errorBar)
         nData  = nData + nComp*nSite

         read(ioDat,'(a400)') header ! header line: need to parse this

         ! set data type and its SI units
         Data(i)%dataType = ImpType(nComp,header)
         allocate(SI_factor(nComp),STAT=istat)
         do j=1,nComp
            SI_factor(j) = ImpUnits(units_in_file,typeDict(Data(i)%dataType)%units(j))
         end do

         do k=1,nSite
         	read(ioDat,*) siteIDs(k), (Data(i)%value(j,k),j=1,nComp)
         	read(ioDat,*)             (Data(i)%error(j,k),j=1,nComp)
         end do

         ! convert data to SI units
         do j=1,nComp
            Data(i)%value(j,:) = SI_factor(j) * Data(i)%value(j,:)
            Data(i)%error(j,:) = SI_factor(j) * Data(i)%error(j,:)
         end do

         ! conjugate data as necessary
         if (conjugate) then
         	do j=2,nComp,2
         		Data(i)%value(j,:) = - Data(i)%value(j,:)
         	end do
         end if

         ! Update the transmitter dictionary and the index (sets up if necessary)
         call update_txDict(Period,2,iTx)
         Data(i)%tx = iTx

         ! Now, update the receiver dictionary and indices
         if (i .eq. 1) then
         	call setup_rxDict(nSite,siteLoc,siteIDs)
         	do k=1,nSite
         		Data(i)%rx(k) = k
         	end do
         else
         	do k=1,nSite
         		call update_rxDict(siteLoc(k,:),siteIDs(k),iRx)
         		Data(i)%rx(k) = iRx
         	end do
         end if
         deallocate(siteLoc,siteIDs,SI_factor,STAT=istat)

	  end do

	  close(ioDat)

	  ! Finally, set up allData
	  nTx = size(txDict)
      call create_dataVecMTX(nTx,allData)
      do iTx = 1,nTx
         ! count number of data types for this transmitter
         nDt = 0
         do i = 1,nBlocks
         	if (Data(i)%tx == iTx) then
         		nDt = nDt+1
         	end if
         end do
         ! initialize the data vectors for each
         call create_dataVecTX(nDt,allData%d(iTx))
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
	     call deall_dataVec(Data(i))
	  end do
	  deallocate(Data,STAT=istat)


  end subroutine read_Z


!**********************************************************************
  subroutine setError_dataVec(err,d)
   !  whether error bars exist or not, sets new error bars
	 !  according to the input parameter err, which specifies the
	 !  fractional error in the output; any errors stored in d
	 !  are overwritten
	 ! TEMPORARY - need to be replaced with something data type specific

   real (kind=prec), intent(in) :: err
   type(dataVec_t),intent(inout)             :: d

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

  end subroutine setError_dataVec

!**********************************************************************
  subroutine setError_dataVecMTX(err,d)
   !  whether error bars exist or not, sets new error bars
	 !  according to the input parameter err, which specifies the
	 !  fractional error in the output; any errors stored in d
	 !  are overwritten

	 real (kind=prec), intent(in) :: err
   type(dataVecMTX_t),intent(inout)          :: d

   !  local variables
   integer      			:: nTx, nData, j, k

   if(.not.d%allocated) then
      call errStop('data vector in setError_dataVecMTX not allocated')
   endif

   do j = 1, d%nTx
     do k = 1, d%d(j)%nDt
        call setError_dataVec(err,d%d(j)%data(k))
     enddo
   enddo

 end subroutine setError_dataVecMTX

end module DataIO
