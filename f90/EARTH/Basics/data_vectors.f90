! *****************************************************************************
module data_vectors
	! Module data_vectors contains the variables and the basic subroutines to do
	! with the data vectors and the penalty functional and their derivatives

  use datadef
  use modeldef
  implicit none

  ! ***************************************************************************
  ! * storing data and transfer functions, dim(nfreq,nfunc,nobs)
  ! * Keep in mind that in general nobs(ifreq) <= nobs; use logical indicators. 
  !type (dataValue_t), allocatable, dimension(:,:,:)	:: dat,psi,res,wres
  type (dataVecMTX_t), save				:: dat,psi,res,wres

  ! ***************************************************************************
  ! * storing data misfit for each freq and func, dim(nfreq,nfunc)
  ! * total number of observations per data type & misfit weights, dim(nfunc)
  !real(8), allocatable, dimension(:,:), save		:: misfit
  integer, allocatable, dimension(:,:), save		:: ndat  
  !real(8), allocatable, dimension(:), save		:: weight
  real(8), allocatable, dimension(:), save		:: misfitValue
  real(8), allocatable, dimension(:,:), save		:: dmisfitValue
  type (misfit_t),save				:: misfit
  
  type (sensitivity_t),save				:: sens

  !type (misfitInfo_t),dimension(:,:),allocatable,save   :: misfitInfo ! (nfreq,nfunc)

end module data_vectors
