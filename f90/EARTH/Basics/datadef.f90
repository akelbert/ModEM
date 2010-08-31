! *****************************************************************************
module datadef
	! This module contains the type definitions for the structures that contain
	! the full information about the data and data functionals, frequencies and
	! other data-related information

  use griddef
  use sg_scalar
  use sg_sparse_vector
  implicit none

  ! --------------------------------------
  ! Define in a globally accessible module
  ! integer										:: nobs
  ! integer										:: nfreq
  ! integer										:: nfunc
  ! type (dataValue_t), dimension(nfreq,nobs)	:: dat
  ! type (Obs_List)								:: obsList
  ! type (Freq_List)							:: freqList
  ! type (TF_List)								:: funcList
  ! Instances of data_func: datC, psiC, resC etc
  ! dat is data with errors
  ! psi are computed functionals
  ! res are residuals
  ! These are vectors dimension(nfreq,nobs)
  ! Use datC%resp%exists == .FALSE. to indicate
  ! that data is not present at a given
  ! frequency + observatory pair.
  ! --------------------------------------


  ! ***************************************************************************
  ! * type sensitivity_t contains the full information about the data sensitivities
  ! * with respect to the original model parameters and to
  ! * each cell resistivity for all frequencies
  ! * used in the Jacobian calculations only
  type :: sensitivity_t

	complex(8), pointer, dimension(:,:,:,:)	:: da  !(nfreq,nfunc,nobs,nvar)
	real(8), pointer, dimension(:,:,:,:)	:: da_real  !(nfreq,nfunc,nobs,nvar)
	real(8), pointer, dimension(:,:,:,:)	:: da_imag  !(nfreq,nfunc,nobs,nvar)
	type (rscalar)							:: drho_real
	type (rscalar)							:: drho_imag
	!real(8), pointer, dimension(:,:,:)		:: drho_real  !(nx,ny,nz) - single freq, all func
	!real(8), pointer, dimension(:,:,:)		:: drho_imag  !(nx,ny,nz) - single freq, all func

  end type sensitivity_t


  ! ***************************************************************************
  ! * type misfit_t contains the full information about the misfit and its'
  ! * partial derivatives with respect to the original model parameters and to
  ! * each cell resistivity for all frequencies
  ! * used in the Jacobian and derivative calculations
  type :: misfit_t

	character(80)							:: name
	real(8)                   :: damping
	real(8), pointer, dimension(:,:)		:: value  !(nfreq,nfunc)
	integer, pointer, dimension(:,:)		:: ndat	  !(nfreq,nfunc)
	real(8), pointer, dimension(:)			:: weight !nfunc
	real(8), pointer, dimension(:,:,:)	    :: dRda	  !(nfreq,nfunc,ncoeff)
	!real(8), pointer, dimension(:,:,:,:)	:: dRda	  !(nfreq,nfunc,nlayer,nparam)

  end type misfit_t


end module datadef
