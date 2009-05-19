! *****************************************************************************
module input
  ! This module defines input file information and input routines

  use math_constants
  use griddef
  use utilities
  use modeldef
  use datadef
  use iotypes
  implicit none

  integer, parameter							:: ioStartup=101
  integer, parameter							:: ioMdl=1
  integer, parameter							:: ioGrd=2
  integer, parameter							:: ioShell=3
  integer, parameter							:: ioPrm=4
  integer, parameter							:: ioPt=32
  integer, parameter							:: ioPer=31
  integer, parameter							:: ioCtrl=16
  integer, parameter							:: ioCond=23
  integer, parameter							:: ioDat=17
  integer, parameter							:: ioObs=18
  integer, parameter							:: ioFunc=19
  integer, parameter							:: ioRad=15

  logical										:: exists ! for I/O inquiries
  character(100)								:: label  ! first line in files

Contains

  ! ***************************************************************************
  ! * readStartFile reads the filename from the screen, if it is not specified.
  ! * This file contains essential information for the program to run.

  subroutine readStartFile(fn_startup,cUserDef)

    character(80), intent(inout)  		:: fn_startup
	type (input_info), intent(out)		:: cUserDef
    integer								:: ios
	character(20)						:: string

    ! passed an empty string
    if (fn_startup == '') then
       ! Prompt the user to give the name for the startup file
       write(0, *) '**********************************************************'
       write(0, *) 'Please, ENTER THE START FILE:'
       read(*, '(a80)') fn_startup
    end if

    open (unit=ioStartup,file=fn_startup,status='old',iostat=ios)

    if(ios/=0) then
       write(0,*) 'Error opening file:', fn_startup
    endif

    ! This is the list of options specified in the startup file

    read (ioStartup,'(a17,a80)') string,cUserDef%paramname;
    read (ioStartup,'(a17,a80)') string,cUserDef%modelname;
    read (ioStartup,'(a17,a80)') string,cUserDef%verbose;
    read (ioStartup,'(a17,a80)') string,cUserDef%calculate;
    read (ioStartup,'(a17,g15.7)') string,cUserDef%damping;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_grid;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_shell;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_period;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_coords;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_func;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_ctrl;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_slices;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_param0;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_param;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_cdata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_ddata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_misfit;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_gradient;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_point;

    close(ioStartup)

	cUserDef%step_size = 1.0d0

  end subroutine readStartFile  ! readStartFile


  ! ***************************************************************************
  ! * initGrid reads the modelfile fn_grid to store the grid inGrid only.
  ! * Traditionally, for global spherical grid, assume the following directions:
  ! * x -> phi (longitude, varies from 0 to 360)
  ! * y -> theta (co-latitude = 90 - latitude; varies from 0 to 180)
  ! * z -> -r (radius from Earth's centre, r=2871.0 at CMB,
  ! *									     6371.0 at Earth/air interface)
  ! *
  ! * For consistency with the original Randy Mackie, 06-27-85, 3-D model, and
  ! * also for consistency with the current forward solver subroutines, we are
  ! * keeping the following grid structure in this forward solver:
  ! *      line 1: dimensions of model (x,y,z)
  ! *      line 2: x(*) in degrees (interval)
  ! *      line 3: y(*) in degrees (position from n-pole)
  ! *      line 4: z(*) in km (distance from center of the earth, decreasing)

  subroutine initGrid(cUserDef,mygrid)

    type (input_info), intent(in)					:: cUserDef
    type (grid_t) , intent(out)					:: mygrid
    integer				                            :: ios,istat,i
    integer				                            :: nx,ny,nz,nzAir,nzCrust
	real(8), dimension(:), allocatable				:: x,y,z


	open(ioGrd,file=cUserDef%fn_grid,status='old',iostat=ios)

    write(6,*) 'Reading from the grid file ',cUserDef%fn_grid
	read(ioGrd,*) nx,ny,nz
	! model grid and resistivity memory allocation
	allocate(x(nx+1),y(ny+1),z(nz+1))

	! read the x-intervals and y in degrees, z in km from the top of the air layer down
    read(ioGrd,*) x(1:nx)
    read(ioGrd,*) y(1:ny+1)
    read(ioGrd,*) z(1:nz+1)

	close(ioGrd)

	! round vertical grid values to nearest meter to prevent precision errors
	do i=1,nz+1
	  z(i)=nearest_meter(z(i))
	end do

	nzAir = 0
	do i=1,nz
	  if (clean(z(i)) > EARTH_R + EPS_GRID) then
		nzAir = nzAir + 1
	  end if
	end do

	nzCrust = 0
	inquire(FILE=cUserDef%fn_shell,EXIST=exists)
	! If no thin shell information present, assume no crust for model computations
	if (.not.exists) then
	  write(0,*) 'Warning: No thin shell conductance distribution specified; assume no crust'
	  nzCrust = nzAir
	else
	  do i=1,nz
		if (clean(z(i)) > CRUST_R + EPS_GRID) then
		  nzCrust = nzCrust + 1
		end if
	  end do
	end if

	! fill in the grid structure
	allocate(mygrid%x(nx),mygrid%y(ny+1),mygrid%z(nz+1),STAT=istat)
	mygrid%nx = nx
	mygrid%ny = ny
	mygrid%nz = nz
	mygrid%nzAir = nzAir
	mygrid%nzEarth = nz - nzAir
	mygrid%nzCrust = nzCrust
	mygrid%nzMantle = nz - nzCrust
    mygrid%x(1:nx)   = x(1:nx)*d2r
    mygrid%y(1:ny+1) = y(1:ny+1)*d2r
    mygrid%z(1:nz+1) = z(1:nz+1)*1000.0D0

	! fill in the coordinates of cell centres; first redefine x
	allocate(mygrid%ph(1:nx+1),mygrid%th(ny+1),mygrid%r(nz+1),STAT=istat)
	x(1) = 0.0d0
	do i=1,nx
	  x(i+1) = x(i)+mygrid%x(i)
	end do
	!x(nx+1) = x(1) ! don't do that - problems with interpolation!!!!
	y(1:ny+1) = mygrid%y(1:ny+1)
	z(1:nz+1) = mygrid%z(1:nz+1)/1000.0d0

	! now save the cell node coordinates in radians and km, respectively
	phi: do i=1,nx+1
	  mygrid%ph(i) = x(i)
	end do phi

	theta: do i=1,ny+1
	  mygrid%th(i) = y(i)
	end do theta

	radius: do i=1,nz+1
	  mygrid%r(i) = z(i)
	end do radius

	deallocate(x,y,z)

	return

  end subroutine initGrid	! initGrid


  ! ***************************************************************************
  ! * initCrust reads the modelfile fn_shell to store S-conductance of Earth's
  ! * crust in GM coordinates. It also stores the depth of the crust in km.
  ! * This should be called after initializing the grid, since it uses the grid
  ! * dimensions to determine the number of entries in file.

  subroutine initCrust(cUserDef,mygrid,mycrust)

    type (input_info), intent(in)					:: cUserDef
    type (grid_t) , intent(in)					:: mygrid
    type (modelShell_t) , intent(out)					:: mycrust
    integer				                            :: ios,istat,i


	inquire(FILE=cUserDef%fn_shell,EXIST=exists)
	if(.not.exists) then
	  mycrust%allocated = .false.
	  return
	end if

	allocate(mycrust%cond(mygrid%nx,mygrid%ny),STAT=istat)

	open(ioShell,file=cUserDef%fn_shell,status='old',iostat=ios)

    write(6,*) 'Reading from the thin sheet conductance file ',cUserDef%fn_shell

	do i=1,mygrid%nx
	  read(ioShell,*,iostat=ios) mycrust%cond(i,1:mygrid%ny)
	end do

	close(ioShell)

	mycrust%allocated = .true.
	return

  end subroutine initCrust	! initCrust

  ! ***************************************************************************
  ! * initFreq reads the file fn_period (currently periods in days, could be
  ! * modified in the future) to store the freq or angular frequencies omega
  subroutine initFreq(cUserDef,myfreq)

	! Not sure yet, whether freq, period or omega will be required
	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (Freq_List), intent(out)							:: myfreq
    integer				                                    :: num
    integer				                                    :: i,j,ios=0
    real(8)								                    :: tmp
    real(8), dimension(:), allocatable				        :: value,days
    character(80)                                                       :: basename

    open(ioPer,file=cUserDef%fn_period,status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the periods file ',cUserDef%fn_period
    read(ioPer,'(a)') label
    ! write(6,*) label

	read(ioPer,*) num
	allocate(value(num),days(num))
	do i=1,num
      read(ioPer,*) days(i) ! reading period in *days*
	  value(i)=1/(days(i)*24*60*60)	! turn into freq.
	  value(i)=clean(value(i))
	end do

!        read(ioPer,*) value(i) ! reading frequency value
!	  days(i)=(1/value(i))/(24*60*60)	! turn into period
!	end do

	close(ioPer)

	! sort the values in ascending order
    do i=1,num
      do j=i+1,num
		if(value(i)>value(j)) then
          tmp=value(i)
          value(i)=value(j)
          value(j)=tmp
		  tmp=days(i)
          days(i)=days(j)
          days(j)=tmp
        end if
      end do
    end do

	!omega(1:num)=2.0d0*pi*value(1:num)
	allocate(myfreq%info(num))
	myfreq%n=num
	myfreq%info(1:num)%value=value(1:num)
	myfreq%info(1:num)%period=days(1:num)
	do i=1,num
	  myfreq%info(i)%i = i
	end do

	deallocate(value)

    basename=cUserDef%fn_period(1:index(cUserDef%fn_period,'.')-1)

    inquire(FILE=trim(basename)//'.codes',EXIST=exists)
    ! If file is not present, set the codes to the frequency numbers
    if (.not.exists) then
       do i=1,num
          write (myfreq%info(i)%code,'(i3.3)') i
       end do
    else
       open(ioPer,file=trim(basename)//'.codes',status='old',form='formatted',iostat=ios)
       read(ioPer,'(a)') label
       read(ioPer,*) num !should be the same as number of frequencies
       do i=1,num
          read(ioPer,*) myfreq%info(i)%code
       end do
       close(ioPer)
    end if

    return

  end subroutine initFreq  !	initFreq


  ! ***************************************************************************
  ! * initCoords reads the file fn_coords that contains the observatory codes
  ! * and locations. We use this information to output responses at these points,
  ! * also to store the data, if computing the Jacobian

  subroutine initCoords(cUserDef,myobs)

	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (Obs_List), intent(out)							:: myobs
    integer				                                    :: num
    integer				                                    :: i,j,ios=0

    open(ioObs,file=cUserDef%fn_coords,status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the coordinates file ',cUserDef%fn_coords
    read(ioObs,'(a)') label
    ! write(6,*) label

	read(ioObs,*) num
	allocate(myobs%info(num))
	do i=1,num

      read(ioObs,*) myobs%info(i)%code,myobs%info(i)%colat,myobs%info(i)%lon
	  myobs%info(i)%lat = 90.0d0 - myobs%info(i)%colat
	  myobs%info(i)%defined = .TRUE.

	end do

	close(ioObs)

	myobs%n = num


    return

  end subroutine initCoords !	initCoords


  ! ***************************************************************************
  ! * initTF reads the file fn_func that contains the information about the
  ! * number and the types of data functionals to use
  subroutine initTF(cUserDef,TFList)

	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (TF_List), intent(out)								:: TFList
    integer				                                    :: num,i,ios

    open(ioFunc,file=cUserDef%fn_func,status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the transfer functions file ',cUserDef%fn_func
    read(ioFunc,'(a)') label

	read(ioFunc,*) num
	allocate(TFList%info(num))
	do i=1,num

      read(ioFunc,*) TFList%info(i)%name,TFList%info(i)%nComp,TFList%info(i)%w

	end do

	close(ioFunc)

	TFList%n = num

  end subroutine initTF	! initTF


  ! *****************************************************************************
  ! * initSlices reads the file fn_slices that contains the information about the
  ! * number and the values of grid radii at which we output full solution
  subroutine initSlices(cUserDef,slices)

	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (Rad_List), intent(out)							:: slices
    integer				                                    :: num,i,ios

	inquire(FILE=cUserDef%fn_slices,EXIST=exists)
	! If file is present, initialize the output radii ("slices")
	if (.not.exists) then
	  write(0,*) 'Warning: No radii specified; we will only output surface solution'
	  allocate(slices%r(1))
	  slices%n = 1
	  slices%r(1) = EARTH_R
	end if

    open(ioRad,file=cUserDef%fn_slices,status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the output radii file ',cUserDef%fn_slices
    read(ioRad,'(a)') label

	read(ioRad,*) num
	allocate(slices%r(num))
	do i=1,num

      read(ioRad,*) slices%r(i)

	end do

	close(ioRad)

	slices%n = num

  end subroutine initSlices	! initSlices


  ! ***************************************************************************
  ! * initMisfit initializes the values required to compute the full penalty
  ! * functional
  subroutine initMisfit(misfitType,TFList,freqList,dat,misfit)

	type (misfitDef_t), intent(in)						:: misfitType
	type (dataVecMTX_t), intent(in)						:: dat
	type (Freq_List), intent(in)						:: freqList
	type (TF_List), intent(in)							:: TFList
	type (misfit_t), intent(out)						:: misfit
	!integer, dimension(:,:), intent(out)				:: ndat
	integer												:: ifunc,ifreq,iobs
	integer												:: istat
	real(8)												:: total_weight

	allocate(misfit%value(freqList%n,TFList%n),STAT=istat)
	allocate(misfit%ndat(freqList%n,TFList%n),STAT=istat)
	allocate(misfit%weight(TFList%n),STAT=istat)

	misfit%name = misfitType%name
	misfit%damping = misfitType%mu

	!ndat = 0
	!do ifunc = 1,TFList%n
  	  ! compute the total number of observations of each type for all frequencies
	!  do ifreq = 1,freqList%n
	!	ndat(ifreq,ifunc) = ndat(ifreq,ifunc) + count(dat%v(ifreq,ifunc,:)%resp%exists)
	!  end do
	!end do

	misfit%ndat = dat%n
	misfit%value = 0.0d0

	total_weight = dot_product(TFList%info%w,sum(dat%n,1))
	do ifunc = 1,TFList%n
	  misfit%weight(ifunc) = TFList%info(ifunc)%w * sum(dat%n(:,ifunc))/total_weight
	end do

  end subroutine initMisfit

  ! ***************************************************************************
  ! * initFunctional specifies basic information about the data functional

  subroutine initFunctional(cUserDef,misfitType)

	implicit none
	type (input_info), intent(in)					:: cUserDef
	type (misfitDef_t), intent(out)					:: misfitType

	misfitType%name = 'Mean Squared'  ! Default data functional

	misfitType%mu   = cUserDef%damping

  end subroutine initFunctional	! initFunctional


  ! ***************************************************************************
  ! * initData initializes the data according to the information stored in
  ! * the list of transfer functions TFList
  subroutine initData(cUserDef,mydat,myobs,myfreq,myfunc)

	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (dataVecMTX_t), intent(inout)						:: mydat
	type (Obs_List), target, intent(inout)					:: myobs
	type (Freq_List), target, intent(inout)					:: myfreq
	type (TF_List), target, intent(in)						:: myfunc
	real(8)                         :: theta
	integer													:: i,j,k

    mydat%nTx = myfreq%n

	do j=1,myfunc%n

	  do i=1,myfreq%n
		do k=1,myobs%n
		  mydat%v(i,j,k)%freq => myfreq%info(i)
		  mydat%v(i,j,k)%obs => myobs%info(k)
		  mydat%v(i,j,k)%func => myfunc%info(j)
		end do
	  end do

	  select case ( myfunc%info(j)%name )

	  case ('C')
		call initDataList(cUserDef%fn_cdata,mydat%v(:,j,:),myobs,myfreq)
		if (count(mydat%v(:,j,:)%resp%exists)==0) then
  		  write(0,*) 'Warning: Data misfit will not be calculated for C responses: data file not available'
		end if
	  do i=1,myfreq%n
		  do k=1,myobs%n
				theta = mydat%v(i,j,k)%obs%colat*d2r;
				if ((theta >= pi/2-EPS_GRID).and.(theta <= pi/2+EPS_GRID)) then
				  write(0,*) 'Error: (initData) receiver ',trim(mydat%v(i,j,k)%obs%code),' located too close to equator;'
				  write(0,*) 'Error: (initData) please delete this receiver from the list and try again.'
					exit
				end if
		    mydat%v(i,j,k)%resp%value = mydat%v(i,j,k)%resp%value/dtan(theta)
		    mydat%v(i,j,k)%resp%err = mydat%v(i,j,k)%resp%err/dtan(theta)
		  end do
	  end do

	  case ('D')
		call initDataList(cUserDef%fn_ddata,mydat%v(:,j,:),myobs,myfreq)
		if (count(mydat%v(:,j,:)%resp%exists)==0) then
  		  write(0,*) 'Warning: Data misfit will not be calculated for D responses: data file not available'
		end if
	  do i=1,myfreq%n
		  do k=1,myobs%n
				theta = mydat%v(i,j,k)%obs%colat*d2r;
		    mydat%v(i,j,k)%resp%value = mydat%v(i,j,k)%resp%value/dsin(theta)
		    mydat%v(i,j,k)%resp%err = mydat%v(i,j,k)%resp%err/dsin(theta)
		  end do
	  end do

	  case default

		write(0,*) 'Please specify correct transfer functions in the file ',trim(cUserDef%fn_func)
		stop

	  end select

	  do i=1,myfreq%n
		mydat%n(i,j) = count(mydat%v(i,j,:)%resp%exists)
	  end do

	end do

  end subroutine initData ! initData

  ! ***************************************************************************
  ! * initData reads the file fn_data (in future could be several files) that
  ! * contain all the information about the available data: the values of the
  ! * responses, data errors corresponding to an observatory and a frequency.
  ! * The logical values stored in global variables freq and obs are also
  ! * initialized here. If in future they are not needed, they can easily be
  ! * deleted.
  ! * Currently only using C-responses! If we ever want to read and access
  ! * D- and J-responses, rename this subroutine to initDataC, file to fn_dataC,
  ! * and write analogous subroutines to initialize these values
  ! * (can be done by copying and pasting, change C to D or J).
  ! * Until this is done, do not even try to access dat%D and dat%J !!!!!!
  subroutine initDataList(fname,mydat,myobs,myfreq)

	implicit none
	character(*), intent(in)					:: fname
	type (dataValue_t), dimension(:,:), intent(inout)		:: mydat
	type (Obs_List), target, intent(inout)				:: myobs
	type (Freq_List), target, intent(inout)				:: myfreq
    integer				                                :: ifreq,iobs,count
    character(100)				                        :: label
    integer				                                :: i,j,ios=0,istat=0
	real(8)								:: const,large
	real(8)								:: freq,days,prev
	character(3)							:: code
	real(8)								:: lon,lat
	real(8)								:: creal,cimag,cerr
	integer								:: nfreq,nobs
	logical								:: new

	large = 2.0e15
	const = 1.0d-2
	prev = 0.0d0

	!nfreq = size(mydat,1)
	!allocate(f(nfreq),STAT=istat)
	!f(:) = mydat(:,1)%freq%info%value
	!nobs  = size(mydat,2)
	!allocate(o(nobs),STAT=istat)
	!o(:) = mydat(1,:)%obs%info%code

	! Initialize the data functionals
	mydat(:,:)%resp%exists = .FALSE.
	mydat(:,:)%resp%value = dcmplx(0.0d0,0.0d0)
	mydat(:,:)%resp%err = large

	inquire(FILE=trim(fname),EXIST=exists)
	! If data file is present, initialize data, functionals and residuals
	if (.not.exists) then
	  return
	end if

    open(ioDat,file=trim(fname),status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the data file ',trim(fname)
    read(ioDat,'(a)') label
    write(6,*) label
    read(ioDat,'(a)') label

    count=0


	READ_DATA: do

	  !read(ioDat,'(a12,1x)',advance='no',iostat=ios)  label !'# freq, obs: '
	  !read(ioDat,*,iostat=ios) freq, num

	  if (ios /= 0) exit

	  read(ioDat,*,iostat=ios) days,code,lon,lat,creal,cimag,cerr

	  freq  = 1/(days*24*60*60)
	  new = .FALSE.
	  if (abs((prev-freq)/freq) > const) then
		new = .TRUE.
		prev = freq
	  end if

	  ! Find the relevant frequency in the list
	  if (new) then
		ifreq = 0
		do i = 1,myfreq%n
  		  if (abs((myfreq%info(i)%value-freq)/myfreq%info(i)%value) < const) then
			print *, 'Reading data for the frequency ',i, freq !myfreq%info(i)%value
			ifreq = i
			exit
		  end if
		end do
	  end if
	  ! If the frequency is not in our list, skip the data
	  if (ifreq == 0) then
		!print *, 'Frequency ',freq,' is not in the frequency list; ignore'
		!do j=1,num
		  !read(ioDat,*,iostat=ios)
		!end do
		cycle READ_DATA
	  end if
	  ! Find the relevant observatory in the list
	  iobs = 0
	  do j = 1,myobs%n
  		if (trim(code) == trim(myobs%info(j)%code)) then
		  iobs = j
		  count = count + 1
		  exit
		end if
	  end do
	  ! If observatory is not in our list, ignore it
	  if (iobs == 0) then
		!print *, 'Observatory ',trim(code),' in file ',trim(fname),' is not in our list'
		!print *, 'This is an error; exiting...'
		!stop
		cycle READ_DATA
	  end if

	  ! If everything is correct, save the data into the variable mydat(ifreq,:)
	  !mydat(ifreq,iobs)%freq => myfreq%info(ifreq)
	  !mydat(ifreq,iobs)%obs	 => myobs%info(iobs)
	  !mydat(ifreq,iobs)%dataFunc	 => myfunc%info(ifunc)
	  mydat(ifreq,iobs)%resp%value = km2m * dcmplx(creal,cimag)
	  mydat(ifreq,iobs)%resp%err = km2m * cerr
	  mydat(ifreq,iobs)%resp%exists = .TRUE.

	end do READ_DATA

	close(ioDat)

	if (count==0) then
  	  write(0,*) 'Warning: No data matches given frequencies and observatories'
	else
  	  write(0,*) 'Number of complex data values in total: ',count
	end if

	! We haven't assigned the targets to some pointers, which is not good, so do it
	!do i = 1, myfreq%n
	 ! do j = 1, myobs%n
	!	if (.not.mydatC(i,j)%resp%exists) then
	!	  mydatC(i,j)%freq => myfreq%info(i)
	!	  mydatC(i,j)%obs  => myobs%info(j)
	!	end if
	 ! end do
	!end do

    return

  end subroutine initDataList !	initData


  ! ***************************************************************************
  ! * initData reads the file fn_data (in future could be several files) that
  ! * contain all the information about the available data: the values of the
  ! * responses, data errors corresponding to an observatory and a frequency.
  ! * The logical values stored in global variables freq and obs are also
  ! * initialized here. If in future they are not needed, they can easily be
  ! * deleted.
  ! * Currently only using C-responses! If we ever want to read and access
  ! * D- and J-responses, rename this subroutine to initDataC, file to fn_dataC,
  ! * and write analogous subroutines to initialize these values
  ! * (can be done by copying and pasting, change C to D or J).
  ! * Until this is done, do not even try to access dat%D and dat%J !!!!!!
  subroutine initDataC(cUserDef,mydatC,myobs,myfreq)

	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (dataValue_t), dimension(:,:), intent(out)			:: mydatC
	type (Obs_List), target, intent(inout)					:: myobs
	type (Freq_List), target, intent(inout)					:: myfreq
    integer				                                    :: ifreq,iobs,num
    character(100)				                            :: label
    integer				                                    :: i,j,ios=0,istat=0
	real(8)													:: const,large
	real(8)													:: freq
	character(3)											:: code
	real(8)													:: lon,lat
	real(8)													:: creal,cimag,cerr

	! Initialize logicals if_C (if_D, if_J) in data types (transmitter_t) and (receiver_t)
	!do i = 1, myfreq%n
	!  do j = 1, myobs%n
	!	allocate(myobs%info(j)%if_C(i),STAT=istat)
	!	allocate(myfreq%info(i)%if_C(j),STAT=istat)
	!	myobs%info(j)%if_C(i) = .FALSE.
	!	myfreq%info(i)%if_C(j) = .FALSE.
	!  end do
	!end do

	large = 2.0e15

	! Initialize the data functionals
	mydatC(:,:)%resp%exists = .FALSE.
	mydatC(:,:)%resp%value = dcmplx(0.0d0,0.0d0)
	mydatC(:,:)%resp%err = large

    open(ioDat,file=cUserDef%fn_cdata,status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the data file ',cUserDef%fn_cdata
    read(ioDat,'(a)') label
    write(6,*) label

	READ_DATA: do

	  read(ioDat,'(a12,1x)',advance='no',iostat=ios)  label !'# freq, obs: '
	  read(ioDat,*,iostat=ios) freq, num

	  if (ios /= 0) exit

	  const = 1.0d-2

	  ! Before reading the values for this frequency, check whether we need them
	  ifreq = 0

	  do i = 1,myfreq%n
  		if (abs((myfreq%info(i)%value-freq)/myfreq%info(i)%value) < const) then
		  print *, 'Reading data for the frequency ',i, freq !myfreq%info(i)%value
		  ifreq = i
		  exit
		end if
	  end do
	  ! If the frequency is not in our list, skip the data
	  if (ifreq == 0) then
		print *, 'Frequency ',freq,' is not in the frequency list; ignore'
		do j=1,num
		  read(ioDat,*,iostat=ios)
		end do
		cycle READ_DATA
	  end if

	  ! If the answer is yes, read them into the variable mydat(ifreq,:)
	  do i=1,num
		read(ioDat,*,iostat=ios) code,lon,lat,creal,cimag,cerr
		! - note - dividing by tan(colat) to obtain ratios
		creal = creal/dtan((90-lat)*d2r)
		cimag = cimag/dtan((90-lat)*d2r)
		cerr = cerr/dtan((90-lat)*d2r)
		! - note - end
		iobs = 0
		do j = 1,myobs%n
  		  if (trim(code) == trim(myobs%info(j)%code)) then
			iobs = j
			exit
		  end if
		end do
		! If observatory is not in our list, output an error - should not happen!..
		if (iobs == 0) then
		  print *, 'Observatory ',trim(code),' in file ',trim(cUserDef%fn_cdata),' is not in our list'
		  print *, 'This is an error; exiting...'
		  stop
		else
		  mydatC(ifreq,iobs)%freq => myfreq%info(ifreq)
		  mydatC(ifreq,iobs)%obs	 => myobs%info(iobs)
		  !mydatC(ifreq,iobs)%dataFunc	 => myfunc%info(ifunc)
		  mydatC(ifreq,iobs)%resp%value = km2m * dcmplx(creal,cimag)
		  mydatC(ifreq,iobs)%resp%err = km2m * cerr
		  mydatC(ifreq,iobs)%resp%exists = .TRUE.
		end if

	  end do

	end do READ_DATA

	close(ioDat)

	! Update logicals if_C (if_D, if_J) in data types (transmitter_t) and (receiver_t)
	do i = 1, myfreq%n
	  do j = 1, myobs%n
		if (mydatC(i,j)%resp%exists) then
		  !myobs%info(j)%if_C(i) = .TRUE.
		  !myfreq%info(i)%if_C(j) = .TRUE.
		else ! We haven't assigned the targets to our pointers, which is not good, so do it
		  mydatC(i,j)%freq => myfreq%info(i)
		  mydatC(i,j)%obs  => myobs%info(j)
		end if
	  end do
	end do

    return

  end subroutine initDataC !	initData

  ! ***************************************************************************
  ! * initData reads the file fn_data (in future could be several files) that
  ! * contain all the information about the available data: the values of the
  ! * responses, data errors corresponding to an observatory and a frequency.
  ! * The logical values stored in global variables freq and obs are also
  ! * initialized here. If in future they are not needed, they can easily be
  ! * deleted.
  subroutine initDataD(cUserDef,mydatD,myobs,myfreq)

	implicit none
    type (input_info), intent(in)							:: cUserDef
	type (dataValue_t), dimension(:,:), intent(out)			:: mydatD
	type (Obs_List), target, intent(inout)					:: myobs
	type (Freq_List), target, intent(inout)					:: myfreq
    integer				                                    :: ifreq,iobs,num
    character(100)				                            :: label
    integer				                                    :: i,j,ios=0,istat=0
	real(8)													:: const,large
	real(8)													:: freq
	character(3)											:: code
	real(8)													:: lon,lat
	real(8)													:: dreal,dimag,derr

	! Initialize logicals if_C (if_D, if_J) in data types (transmitter_t) and (receiver_t)
	!do i = 1, myfreq%n
	!  do j = 1, myobs%n
	!	allocate(myobs%info(j)%if_D(i),STAT=istat)
	!	allocate(myfreq%info(i)%if_D(j),STAT=istat)
	!	myobs%info(j)%if_D(i) = .FALSE.
	!	myfreq%info(i)%if_D(j) = .FALSE.
	!  end do
	!end do

	large = 2.0e15

	! Initialize the data functionals
	mydatD(:,:)%resp%exists = .FALSE.
	mydatD(:,:)%resp%value = dcmplx(0.0d0,0.0d0)
	mydatD(:,:)%resp%err = large

    open(ioDat,file=cUserDef%fn_ddata,status='old',form='formatted',iostat=ios)

    write(6,*) 'Reading from the data file ',cUserDef%fn_ddata
    read(ioDat,'(a)') label
    write(6,*) label

	READ_DATA: do

	  read(ioDat,'(a12,1x)',advance='no',iostat=ios)  label !'# freq, obs: '
	  read(ioDat,*,iostat=ios) freq, num

	  if (ios /= 0) exit

	  const = 1.0d-2

	  ! Before reading the values for this frequency, check whether we need them
	  ifreq = 0

	  do i = 1,myfreq%n
  		if (abs((myfreq%info(i)%value-freq)/myfreq%info(i)%value) < const) then
		  print *, 'Reading data for the frequency ',i, freq !myfreq%info(i)%value
		  ifreq = i
		  exit
		end if
	  end do
	  ! If the frequency is not in our list, skip the data
	  if (ifreq == 0) then
		print *, 'Frequency ',freq,' is not in the frequency list; ignore'
		do j=1,num
		  read(ioDat,*,iostat=ios)
		end do
		cycle READ_DATA
	  end if

	  ! If the answer is yes, read them into the variable mydat(ifreq,:)
	  do i=1,num
		read(ioDat,*,iostat=ios) code,lon,lat,dreal,dimag,derr
		! - note - dividing by sin(colat) to obtain ratios
		dreal = dreal/dsin((90-lat)*d2r)
		dimag = dimag/dsin((90-lat)*d2r)
		derr = derr/dsin((90-lat)*d2r)
		! - note - end
		iobs = 0
		do j = 1,myobs%n
  		  if (trim(code) == trim(myobs%info(j)%code)) then
			iobs = j
			exit
		  end if
		end do
		! If observatory is not in our list, output an error - should not happen!..
		if (iobs == 0) then
		  print *, 'Observatory ',trim(code),' in file ',trim(cUserDef%fn_ddata),' is not in our list'
		  print *, 'This is an error; exiting...'
		  stop
		else
		  mydatD(ifreq,iobs)%freq => myfreq%info(ifreq)
		  mydatD(ifreq,iobs)%obs	 => myobs%info(iobs)
		  !mydatD(ifreq,iobs)%dataFunc	 => myfunc%info(ifunc)
		  mydatD(ifreq,iobs)%resp%value = km2m * dcmplx(dreal,dimag)
		  mydatD(ifreq,iobs)%resp%err = km2m * derr
		  mydatD(ifreq,iobs)%resp%exists = .TRUE.
		end if

	  end do

	end do READ_DATA

	close(ioDat)

	! Update logicals if_C (if_D, if_J) in data types (transmitter_t) and (receiver_t)
	do i = 1, myfreq%n
	  do j = 1, myobs%n
		if (mydatD(i,j)%resp%exists) then
		  !myobs%info(j)%if_D(i) = .TRUE.
		  !myfreq%info(i)%if_D(j) = .TRUE.
		else ! We haven't assigned the targets to our pointers, which is not good, so do it
		  mydatD(i,j)%freq => myfreq%info(i)
		  mydatD(i,j)%obs  => myobs%info(j)
		end if
	  end do
	end do

    return

  end subroutine initDataD !	initData



  ! ***************************************************************************
  ! * initModelParam reads the parametrization info to store in the derived data
  ! * type variables (modelCoeff_t, modelLayer_t) in the special case when the
  ! * parametrization is in terms of spherical harmonics of various degree/order
  ! * per layer, each of the parameters being either variable (keyword 'range') or
  ! * constant (keyword 'const'). We keep this information internally in such
  ! * a format, that generalizations of this case would be easy to implement.
  ! * In the parametrization type variables, each layer has the same (maximum)
  ! * degree and order of spherical harmonics, and hence identical numbers of
  ! * coefficients to store. Out of these coefficients, those that have been
  ! * defined in the script bear logical keyword 'exists'. Out of those, variable
  ! * coefficients have logical 'frozen==.FALSE.'
  ! * Obviously, the information on the range is not required for the forward
  ! * solver to operate. This is provided for the inversion, which will share
  ! * the same input format for now.

  subroutine initModelParam(cUserDef,myparam,p0)

	use model_operators

    type (input_info), intent(in)					:: cUserDef
    type (modelParam_t), intent(inout)					:: myparam
	logical, intent(in), optional		:: p0
    integer								:: ilayer,i,j,k,n,l,m
	integer								:: nF,nL
    integer								:: sum,sum0,degree
    integer								:: ios,istat
    real(8)								:: upperb,lowerb,width,depth,alpha,beta
    character(5)							:: if_log_char
    character(5)							:: if_var_char
    logical								:: if_log, if_fixed
    character(80)							:: prmname, string
    real(8)                                :: v,min,max

    lowerb = EARTH_R
    depth = 0.0d0
    width = 0.0d0
    sum = 0
    sum0 = 0

	if(present(p0)) then
	  if(p0) then
		open(ioPrm,file=cUserDef%fn_param0,status='old',form='formatted',iostat=ios)
		write(6,*) 'Reading from the parametrization file ',cUserDef%fn_param0
	  end if
	else
	  open(ioPrm,file=cUserDef%fn_param,status='old',form='formatted',iostat=ios)
	  write(6,*) 'Reading from the parametrization file ',cUserDef%fn_param
	end if

    read(ioPrm,'(a8,a80)') string,prmname

	if (index(prmname,'harmonic')==0) then
       write(0, *) 'Error: (initModelParam) not a spherical harmonic parametrization'
       stop
	else
	  i = index(prmname,'layers')
	  read(prmname(i+7:len(prmname)),'(i2)') nL
	  i = index(prmname,'degree')
	  read(prmname(i+7:len(prmname)),'(i2)') degree
	  write(6,*) prmname
	end if


	call create_modelParam(myparam,nL,degree)

    !write(6,'(a50,i3)') 'Number of layers in script: ',myparam%nL


    ! Unwind spherical harmonic parametrization

    do n=1,nL

	   read(ioPrm, *)

       upperb = lowerb

 	   ! If you wish to read the line in sections, use advance='no' specifier
     !  read(ioPrm,'(a3,1x,a6,1x,i2,a6,g15.7)',iostat=ios) &
			!if_log_char,string,degree,string,depth

			!read(ioPrm,'(a3)',iostat=ios,advance='no') if_log_char

			read(ioPrm,'(a80)',iostat=ios) string
			i = index(string,'degree')
			j = index(string,'layer')
			k = index(string,'reg')


			if (k==0) then
					! no regularisation specified for this layer
					alpha = 0.0d0
					beta  = 1.0d0
					k = len(string)
			else
					read(string(k+4:len(string)),*) alpha,beta
			end if

			read(string(1:i-1),*) if_log_char
			read(string(i+7:j),*) degree
			read(string(j+6:k),*) depth

       read(ioPrm,'(a80)',iostat=ios) string

       if (if_log_char == 'log') then  ! log means log_{10}
          if_log = .TRUE.
       else
          if_log = .FALSE.
       end if
	   lowerb = EARTH_R - depth

	   call setLayer_modelParam(myparam,n,upperb,lowerb,alpha,beta,if_log)

       sum0=sum0+sum
       sum=0
       do l=0,degree
          sum = sum + (2*l+1)
       end do
       !write(6,'(a46,i2,a2,i3)') 'Number of coefficients in layer ',n,': ',sum

       do i=1,sum
          read(ioPrm,*,iostat=ios) l,m,v,min,max,if_var_char
		  if (if_var_char == 'range') then
			if_fixed = .FALSE.
		  else if (if_var_char == 'const') then
			if_fixed = .TRUE.
		  else
			write(0, *) 'Error: (initModelParam) wrong character constant for ',n,l,m
			stop
		  end if
		  call setCoeffValue_modelParam(myparam,n,l,m,v,min,max,if_fixed)
		  !print *,'Values: ',l,m,v,min,max,if_fixed
       end do

    end do

    ! this check is not necessary, but helps to debug parametrization scripts
    !write(6,'(a50,i3)') 'Number of variable parameters in script: ',count(.not.myparam%c%frozen)
	write(6,*)

    close(ioPrm)

    return

  end subroutine initModelParam	! initModelParam


  ! ***************************************************************************
  ! * initControls reads the file fn_ctrl and sets the values of main control
  ! * parameters for the forward solver, stored in fwdCtrls

  subroutine initControls(cUserDef,fwdCtrls)

	implicit none
	type (input_info), intent(in)						:: cUserDef
	type (fwdCtrl_t), intent(out)						:: fwdCtrls

	  open(ioCtrl,file=cUserDef%fn_ctrl,form='formatted',status='old')

      write(6,*) 'Reading from the forward solver controls file ',cUserDef%fn_ctrl
      read(ioCtrl,*) fwdCtrls%ipotloopmax
      read(ioCtrl,*) fwdCtrls%errend
      read(ioCtrl,*) fwdCtrls%nrelmax
      read(ioCtrl,*) fwdCtrls%n_reldivh
      read(ioCtrl,*) fwdCtrls%ipot0,fwdCtrls%ipotint,fwdCtrls%ipot_max

      close(ioCtrl)

  end subroutine initControls	! initControls


  ! ***************************************************************************
  ! * initOutput initialises all the output file names stored in outFiles
  ! * Assuming type (output_info) outFiles and type (input_info) cUserDef are
  ! * not available.
  subroutine initOutput(cUserDef,outFiles)

	implicit none
    type (input_info) ,intent(in)					:: cUserDef
    type (output_info),intent(out)					:: outFiles
	integer											:: i

    if (cUserDef%modelname == '') then
       write(0, *) 'modelname not specified yet for FileInfoInit'
       stop
    end if

	i=index(cUserDef%modelname,' ')
	i=i-1

	outFiles%fn_hx	  =cUserDef%modelname(1:i)//'.hx';
	outFiles%fn_hy	  =cUserDef%modelname(1:i)//'.hy';
	outFiles%fn_hz	  =cUserDef%modelname(1:i)//'.hz';
	outFiles%fn_jxjyjz =cUserDef%modelname(1:i)//'.jxjyjz';
	outFiles%fn_hxhyhz =cUserDef%modelname(1:i)//'.hxhyhz';
	outFiles%fn_err	  =cUserDef%modelname(1:i)//'.err';
	outFiles%fn_cresp  =cUserDef%modelname(1:i)//'.cresp';
	outFiles%fn_dresp  =cUserDef%modelname(1:i)//'.dresp';
	outFiles%fn_cdat  =cUserDef%modelname(1:i)//'.cout';
	outFiles%fn_ddat  =cUserDef%modelname(1:i)//'.dout';
	outFiles%fn_cjac  =cUserDef%modelname(1:i)//'.cj';
	outFiles%fn_djac  =cUserDef%modelname(1:i)//'.dj';
	outFiles%fn_model  =cUserDef%modelname(1:i)//'.rho';
	outFiles%fn_residuals  =cUserDef%modelname(1:i)//'.res';
	outFiles%fn_bv	  =cUserDef%modelname(1:i)//'.bv';

	outFiles%fn_avg_cresp='earth_response.cdat';
	outFiles%fn_avg_dresp='earth_response.ddat';


  end subroutine initOutput	! initOutput

end module input
