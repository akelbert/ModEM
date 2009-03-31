module ModelSpace

  use griddef
  use math_constants
  use utilities
  use sg_scalar
  use sg_vector
  use sg_sparse_vector
  use modelparameter

  implicit none

  type :: modelVec_t

     private
     type (modelParam_t)         :: m
     type (grid_t), pointer      :: grid
     logical                     :: loge = .true.
     logical			         :: allocated = .false.

  end type modelVec_t

  interface assignment (=)
  	MODULE PROCEDURE copy_modelVec
  end interface

  interface operator (*)
  	MODULE PROCEDURE scMult_modelVec
  end interface

  interface dotProd
  	MODULE PROCEDURE dotProd_modelVec
  end interface

  interface linComb
  	MODULE PROCEDURE linComb_modelVec
  end interface

  interface zero
  	MODULE PROCEDURE zero_modelVec
  end interface

  interface deall
  	MODULE PROCEDURE deall_modelVec
  end interface

  interface read
  	MODULE PROCEDURE zero_modelVec
  end interface

  interface write
  	MODULE PROCEDURE deall_modelVec
  end interface

  public  :: create_modelVec, deall_modelVec, zero_modelVec, copy_modelVec
  public  :: linComb_modelVec, dotProd_modelVec, scMult_modelVec
  public  :: read_modelVec, write_modelVec, readVec_modelVec, writeVec_modelVec

Contains

  !**********************************************************************
  ! creates mVec from a model parameter;
  ! for internal use within the modelParam group of modules only

  subroutine create_modelVec(grid,m,mVec,linear)

     implicit none
     type (grid_t), intent(in), target	:: grid
     type (modelParam_t), intent(in)        :: m
     type (modelVec_t), intent(out)         :: mVec
     logical, intent(in), optional          :: linear

     if(mVec%allocated) then
        call deall_modelVec(mVec)
     endif

     mVec%allocated = .true.
     mVec%m = m
     mVec%grid => grid

     if (.not. present(linear)) then
     	mVec%loge = .true.
     else
     	mVec%loge = .not. linear
     endif

  end subroutine create_modelVec

  !************************************************************

  subroutine deall_modelVec(mVec)

     type (modelVec_t), intent(inout)   :: mVec

     if(mVec%allocated) then
        call deall_modelParam(mVec%m)
        nullify(mVec%grid)
     endif

     mVec%allocated = .false.

  end subroutine deall_modelVec

  !**********************************************************************

  subroutine zero_modelVec(mVec)

     type(modelVec_t), intent(inout)    :: mVec

     call zero_modelParam(mVec%m)

  end subroutine zero_modelVec

  !**********************************************************************

  subroutine copy_modelVec(mOut,mIn)

     type (modelVec_t), intent(in)       :: mIn
     type (modelVec_t), intent(inout)    :: mOut

     if (mOut%allocated) then
     	call deall_modelVec(mOut)
     endif

     mOut%m = mIn%m
     mOut%grid => mIn%grid
     mOut%loge = mIn%loge
     mOut%allocated = .true.

  end subroutine copy_modelVec

  !**********************************************************************
  !  forms the linear combination of model vectors
  !    m = a1*m1 + a2*m2
  !  where a1 and a2 are real constants and m1 and m2
  !   are model vectors
  !   output m may overwrite m1 or m2

  subroutine linComb_modelVec(a1,mVec1,a2,mVec2,mVec)


     real (kind=prec), intent(in)       :: a1,a2
     type (modelVec_t), intent(in)		:: mVec1,mVec2
     type (modelVec_t), intent(inout)	:: mVec

     if (mVec%allocated) then
     	call deall_modelVec(mVec)
     endif

     if ((.not. mVec1%allocated) .or. (.not. mVec2%allocated)) then
     	call errStop('input model vectors not allocated in linComb_modelVec')
     endif

     if (mVec1%loge /= mVec2%loge) then
     	call errStop('input model vectors incompatible in linComb_modelVec')
     endif

     call linComb_modelParam(a1, mVec1%m, a2, mVec2%m, mVec%m)

     mVec%grid => mVec1%grid
     mVec%loge = mVec1%loge
     mVec%allocated = .true.

  end subroutine linComb_modelVec

  !**********************************************************************

  function dotProd_modelVec(mVec1,mVec2) result(r)

     real (kind=prec)                   :: r
     type (modelVec_t), intent(in)      :: mVec1,mVec2

     r = dotProd_modelParam_f(mVec1%m, mVec2%m)

  end function dotProd_modelVec

  !**********************************************************************
  !  computes mOut = a * mIn for modelVec object m and real scalar a

  function scMult_modelVec(a,mIn) result (mOut)

    real (kind=prec), intent(in)		:: a
    type (modelVec_t), intent(in)	    :: mIn
    type (modelVec_t)                   :: mOut

    if (.not. mIn%allocated) then
       call errStop('input not allocated on call to scMult_modelVec')
    endif

	call linComb_modelVec(R_ZERO, mIn, a, mIn, mOut)

  end function scMult_modelVec


  !**********************************************************************
  !  we have to assume that modelParam and grid are in different files,
  !  since they are conceptually different entities;
  !  if they are stored in the same file (e.g. formats of Randie Mackie
  !  or Weerachai Siripunvaraporn), write two separate routines to extract
  !  the model parameter and the grid

  subroutine read_modelVec(mVec, cfile_model, cfile_grid)

    type (modelVec_t), intent(out)	        :: mVec
    character (:), intent(in)		        :: cfile_model
    character (:), intent(in), optional     :: cfile_grid
    ! local
    integer                                 :: status

	if (present(cfile_grid)) then
		call read_grid(mVec%grid, cfile_grid, status)
	else
		call read_grid(mVec%grid, cfile_model, status)
	endif

	call read_modelParam(mVec%m, cfile_model, status)

	mVec%loge = .true.

	if (status == 0) then
		mVec%allocated = .true.
	endif

  end function read_modelVec


  !**********************************************************************
  !  we have to assume that modelParam and grid are in different files,
  !  since they are conceptually different entities;
  !  if they are stored in the same file (e.g. formats of Randie Mackie
  !  or Weerachai Siripunvaraporn), we just write out the model param.
  !  of course, in this case that will have to contain some grid info!

  subroutine write_modelVec(mVec, cfile_model, cfile_grid)

    type (modelVec_t), intent(out)	        :: mVec
    character (:), intent(in)		        :: cfile_model
    character (:), intent(in), optional     :: cfile_grid
    ! local
    integer                                 :: status

	if (present(cfile_grid)) then
		call write_grid(mVec%grid, cfile_grid, status)
	endif

	call write_modelParam(mVec%m, cfile_model, status)

  end function write_modelVec


end module ModelSpace
