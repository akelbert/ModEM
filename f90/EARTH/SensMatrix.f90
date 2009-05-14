module sensMatrix

use griddef
use data_vectors
use global, only: cUserDef, fwdCtrls, outFiles, obsList, TFList, freqList, slices
use dataFunc
use jacobian
use output
use initFields
use dataMisfit

implicit none

public 	:: calcSensMatrix, Jmult, JmultT, fwdPred, setGrid

! numerical discretization used to compute the EM solution
!  (may be different from the grid stored in model parameter)
type(grid_t), target, save, private     :: grid

! utility variables necessary to time the computations;
!  including a public variable rtime that stores total run time
real, save, public						:: rtime  ! run time
real, save, private						:: ftime  ! run time per frequency
real, save, private						:: stime, etime ! start and end times
integer, dimension(8), private			:: tarray ! utility variable

Contains

   !**********************************************************************
   subroutine calcSensMatrix(d,sigma0,dsigma)
   !  Calculate sensitivity matrix for data in d
   !
   !   d is the input data vector, here just used to identify
   !     receiver transmitter pairs to compute sensitivities for
   type(dataVecMTX_t), intent(in)	:: d
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the output array of data sensitivities,
   !   one for each element in the data array.  Each sensitivity
   !    is an element of type modelParam, an abstract
   !    data type that defines the unknow conductivity
   type(modelParam_t), pointer   :: dsigma(:)


   end subroutine calcSensMatrix

   !**********************************************************************
   subroutine Jmult(dsigma,sigma0,d,eAll)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for all transmitters in a datavector (i.e., multiple dataVec objects)
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions)
   !
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the input conductivity parameter perturbation
   type(modelParam_t), intent(in)	:: dsigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVecMTX_t), intent(inout)		:: d
   type(EMsolnMTX_t), intent(inout), optional	:: eAll



   end subroutine Jmult

   !**********************************************************************
   subroutine JmultT(sigma0,d,dsigma,eAll)

   !  Transpose of Jmult mujltiplied by data vector d; output is a
   !      single conductivity parameter in dsigma
   !
   !   sigma0 is background conductivity parameter
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   IN THE PROPER ORDER (at present) !!!!
   !   (This can be used to avoid recomputing forward solutions,
   !    e.g., in a CG solution scheme)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  "pointers" to dictionary entries are attached to multi-transmitter
   !   data vector d

   type(modelParam_t), intent(in)	:: sigma0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVecMTX_t), intent(in)		:: d
   !   dsigma is the output conductivity parameter
   type(modelParam_t), intent(Out)  	:: dsigma
   type(EMsolnMTX_t), intent(in), optional	:: eAll



   end subroutine JmultT


   !**********************************************************************
   subroutine fwdPred(m,d,H)

   !  Calculate predicted data for dataVecMTX object d
   !    and for conductivity parameter sigma
   !  Optionally returns array of EM solutions eAll, one for
   !    each transmitter (if eAll is present)
   !
   !  First need to set up transmitter, receiver, dataType dictionaries;
   !  "pointers" to dictionaries are attached to multi-transmitter
   !   data vector d .
   !
   !   sigma is input conductivity parameter
   type(modelParam_t), intent(in)	:: m
   !   d is the computed (output) data vector, also used to identify
   !     receiver/transmitter
   type(dataVecMTX_t), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(EMsolnMTX_t), intent(inout), optional	:: H
   ! local variables
    integer	                                :: errflag	! internal error flag
	real(8)									:: omega  ! variable angular frequency
	integer									:: istat,i,j,k
    type (cvector)							:: Hj,B,F,Hconj,B_tilde,dH,dE,Econj,Bzero,dR
	type (rvector)							:: dE_real
	type (rscalar)							:: drho
	type (sparsevecc)						:: Hb
	type (functional_t)						:: dataType
	type (transmitter_t)					:: freq
	type (modelParam_t)						:: dmisfit,dmisfit2
	integer									:: ifreq,ifunc
	logical									:: adjoint,delta
	character(1)							:: cfunc
	character(80)							:: fn_err
	type(dataVecMTX_t)                      :: dat,psi
	integer                                 :: nfreq,nfunc,nobs

    if(.not.d%allocated) then
       call errStop('data vector not allocated on input to fwdPred')
    end if

    if(present(H)) then
       if(.not. H%allocated) then
          call create_EMsolnMTX(d%nTx,H,grid)
       else if(d%nTx .ne. H%nTx) then
          call errStop('dimensions of H and d do not agree in fwdPred')
       end if
    end if

    ! Allocate temporary data and response storage
    nfreq = freqList%n
    nfunc = TFList%n
    nobs  = obsList%n
	allocate(dat%v(nfreq,nfunc,nobs),dat%n(nfreq,nfunc),STAT=istat)
	allocate(psi%v(nfreq,nfunc,nobs),psi%n(nfreq,nfunc),STAT=istat)
	dat%v = d%v
	dat%n = d%n

    ! Allocate the resistivity vector
	call create_rscalar(grid,rho,CENTER)

	! Compute model information everywhere else in the domain
	call initModel(grid,m,rho%v)

	! Start the (portable) clock
	call date_and_time(values=tarray)

	call initialize_fields(Hj,B)

	do ifreq=1,freqList%n

	  stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

	  freq = freqList%info(ifreq)

	  omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

	  write(6,*) 'Solving 3D forward problem for freq ',ifreq,freq%value

	  ! solve A <h> = <b> for vector <h>
	  adjoint=.FALSE.
	  call operatorM(Hj,B,omega,rho%v,grid,fwdCtrls,errflag,adjoint)
	  !call create_cvector(grid,Hj,EDGE)
	  !Hj%x = C_ONE
	  !Hj%y = C_ONE
	  !Hj%z = C_ONE

	  ! compute and output fields & C and D responses at cells
	  call outputSolution(freq,Hj,slices,grid,cUserDef,rho%v,'h')

	  ! compute and output C and D responses at observatories
	  call calcResponses(freq,Hj,dat,psi)
	  call outputResponses(freq,psi,freqList,TFList,obsList,outFiles,dat)

	  call date_and_time(values=tarray)
	  etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
	  ftime = etime - stime
	  print *,'Time taken (secs) ',ftime
	  print *
	  rtime = rtime + ftime

	  if(present(H)) then
	     H%solns(ifreq) = Hj
	     H%tx(ifreq) = ifreq
	     H%errflag(ifreq) = errflag
	  end if

	end do

    d%v = psi%v
    d%n = psi%n

	deall_rscalar(rho)

   end subroutine fwdPred

   !**********************************************************************
   subroutine setGrid(newgrid)

   !  Use to set and/or update the numerical grid, that is then used
   !   all computations in this module;
   !   This is not a pointer target.

   type(grid_t), intent(in)     :: newgrid
   ! local
   integer                      :: istat

    grid = newgrid

   	! Too complicated to rewrite input to all subroutines that use x,y,z,nx,ny,nz
	! in terms of the grid variable, but that would be the way to do it in the
	! future. Currently, use this patch instead. x,y,z stored in module griddef
	nx = grid%nx; ny = grid%ny; nz = grid%nz
	nzEarth = grid%nzEarth; nzAir = grid%nzAir
	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = grid%x; y = grid%y; z = grid%z

   end subroutine setGrid

end module sensMatrix
