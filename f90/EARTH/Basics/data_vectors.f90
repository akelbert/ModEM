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
  !type (data_value), allocatable, dimension(:,:,:)	:: dat,psi,res,wres
  type (data_array), save				:: dat,psi,res,wres

  ! ***************************************************************************
  ! * storing data misfit for each freq and func, dim(nfreq,nfunc)
  ! * total number of observations per data type & misfit weights, dim(nfunc)
  !real(8), allocatable, dimension(:,:), save		:: misfit
  integer, allocatable, dimension(:,:), save		:: ndat  
  !real(8), allocatable, dimension(:), save		:: weight
  real(8), allocatable, dimension(:), save		:: misfitValue
  real(8), allocatable, dimension(:,:), save		:: dmisfitValue
  type (data_misfit),save				:: misfit
  
  type (sensitivity),save				:: sens

  !type (misfit_info),dimension(:,:),allocatable,save   :: misfitInfo ! (nfreq,nfunc)

end module data_vectors
