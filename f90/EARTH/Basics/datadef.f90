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


  ! ***************************************************************************
  ! * contains the list of transfer functions
  type :: TF_List

	integer										:: n
	type (functional_t), pointer, dimension(:)	:: info	!nfunc

  end type TF_List


  ! ***************************************************************************
  ! * type functional_t contains the definitions of a data functional used in
  ! * the expression for the penalty functional; we will further specify the rules
  ! * by which they are computed, perhaps by pointing at a function
  type :: functional_t

	! Here specify rules used to compute this data functional
	character(80)						    :: name
    ! nComp is number of complex EM components observed (e.g., 3 for
    ! MT (with vertical field TFs))
	integer									:: nComp
	! The weight of this data type in the calculation of misfit
	real(8)									:: w

  end type functional_t

  ! ***************************************************************************
  ! * type response_t contains the full information about a single response at a
  ! * single frequency and observatory: value, error, anything else
  type :: response_t

	complex(8)								:: value
	real(8)									:: err
	logical									:: exists

  end type response_t


  ! ***************************************************************************
  ! * type dataValue_t contains the definition of a single data functionals used
  ! * in the expression for the penalty functional; they are specified for a
  ! * particular type of response and a frequency + observatory pair
  type :: dataValue_t

	type (transmitter_t), pointer			    :: freq
	type (receiver_t),	pointer			        :: obs
	type (functional_t),	pointer				:: func
	type (response_t)							:: resp
	! For the more general case, define dimension(freq%nMode,dataType%nComp)
	! type (response_t), pointer, dimension(:,:) :: resp

  end type dataValue_t


  ! ***************************************************************************
  ! * type dataVecMTX_t contains the 3-D array of all data functionals we define;
  ! * it should be considered an ordered three-dimensional list of data
  ! * For some observatories data "do not exist" - either because they were
  ! * not physically present in the input data file, or since the observatory
  ! * is decided to be too close to one of the poles or the equator, so that
  ! * the responses at these locations would not be reliable. In these cases
  ! * we set v(i,j,k)%resp%exists to .FALSE.
  ! * The variable n is the number of "existing" data values for a given freq.
  ! * and functional type. It can also be computed by count(v(i,j,:)%resp%exists);
  ! * however, it is handy to have these values at hand.
  type :: dataVecMTX_t

	type (dataValue_t), pointer, dimension(:,:,:) :: v	!nfreq,nfunc,nobs
	integer, pointer, dimension(:,:)			  :: n	!nfreq,nfunc
	! Total number of frequencies stored in this data vector
	integer                                       :: ntx
	logical                                       :: allocated=.FALSE.

  end type dataVecMTX_t


  ! ***************************************************************************
  ! * list of radii ("slices") at which we output the full solution
  type :: Rad_List

	integer										:: n
	real(8), pointer, dimension(:)				:: r	!nrad

  end type Rad_List


  ! ***************************************************************************
  ! * full solution information for a single 'slice' at a fixed radius (in km)
  type :: solution_t

	type (receiver_t), pointer, dimension(:,:)	:: o
	complex(8), pointer, dimension(:,:)			:: x,y,z	!nx,ny

  end type solution_t


  ! ***************************************************************************
  ! * contains the full information about the observatories only
  type :: Obs_List

	integer										:: n
	type (receiver_t), pointer, dimension(:)		:: info	!nobs

  end type Obs_List


  ! ***************************************************************************
  ! * type receiver_t contains the information about a single observatory; we
  ! * define as many of them as there are observatories (nobs)
  ! * This definition of the receiver has the limitation that it is assumed to
  ! * be located on an ij-plane. It can be extended to being located at any
  ! * point in the domain, if required, by modifying the LocateReceiver code,
  ! * and making the ComputeInterpWeights a 3-D interpolation routine.
  type :: receiver_t

	! observatory code that usually has three letters, but may have more
	character(80)							:: code
	! observatory location: co-latitude and longitude in degrees
	real(8)									:: colat,lat,lon
	! observatory location: radius in km (default is EARTH_R)
	real(8)									:: rad = EARTH_R
	! once you define a receiver, need to set defined to TRUE
	logical									:: defined=.FALSE.

	! these values specify the location of the receiver relative to the grid,
	! required for interpolation. They only need to be computed once.
	integer									:: i,j,k
	! the proportions for the distance from the cell corner, for linear interp.
	real(8)									:: p_crn,q_crn
	! the proportions for the distance from the cell center, for bilinear interp.
	real(8)									:: p_ctr,q_ctr
	! the location of the observatory relative to the center of the face (which side
	! from the mid-face) can be derived from comparing the values p_crn and q_crn to 1/2.

	! vectors that store the weights for interpolation
	type (sparsevecc)						:: Lx,Ly,Lz

	! indicator located is set to TRUE once the above values are computed
	logical									:: located=.FALSE.

  end type receiver_t


  ! ***************************************************************************
  ! * storing the frequency dimensions and values in ascending order
  type :: Freq_List

	integer										:: n
	type (transmitter_t), pointer, dimension(:)	:: info	!nfreq

  end type Freq_List


  ! ***************************************************************************
  ! * type transmitter_t contains the information about the frequencies; and also
  ! * the info on data functionals required to calculate for each frequency
  type :: transmitter_t

	! the code may be needed to indicate frequency number in the full data set
    character(80)							:: code
	! frequency value
	real(8)		  							:: value
	! period in days used for input/output only
	real(8)		  							:: period ! (in days)
	! order index of this frequency that is used for the output
	integer									:: i
    ! nMode is number of "modes" for transmitter (e.g., 2 for MT)
	integer									:: nMode

  end type transmitter_t


end module datadef
