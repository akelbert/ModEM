module ModelSpace

  use griddef
  use math_constants
  use utilities
  use modelparameter ! inherits EMfield

  implicit none

  type :: modelVec_t

     !private
     type (modelParam_t)         :: m
     type (grid_t), pointer      :: grid
     logical			         :: allocated = .false.

  end type modelVec_t

  interface assignment (=)
  	MODULE PROCEDURE copy_modelVec
  end interface

  interface operator (*)
  	MODULE PROCEDURE scMult_modelVec_f
  end interface

  interface dotProd
  	MODULE PROCEDURE dotProd_modelVec_f
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

  ! generic
  private :: create_modelVec
  public  :: deall_modelVec, zero_modelVec, copy_modelVec
  public  :: linComb_modelVec, dotProd_modelVec_f, scMult_modelVec_f
  public  :: read_modelVec, write_modelVec, readVec_modelVec, writeVec_modelVec
  public  :: getGrid_modelVec
  ! not generic, but mappings from model to grid have to exist
  public  :: ModelToOneCell, ModelToCell, dCellToModel
  public  :: dModelToEdge, dEdgeToModel, dModelToNode, dNodeToModel

 ! definitions for CmSqrt: must be consistent with the include file below
  include "modelParam/modelCov/Diffusion.hd"

Contains

  !  The included file must contain subroutines create_CmSqrt, deall_CmSqrt, multBy...
  include "modelParam/modelCov/Diffusion.inc"

  !**********************************************************************
  ! creates mVec from a model parameter;
  ! for internal use within the modelParam group of modules only

  subroutine create_modelVec(grid,m,mVec)

     implicit none
     type (grid_t), intent(in), target	 :: grid
     type (modelParam_t), intent(in)     :: m
     type (modelVec_t), intent(out)      :: mVec

     if(mVec%allocated) then
        call deall_modelVec(mVec)
     endif

     mVec%allocated = .true.
     mVec%m = m
     mVec%grid => grid

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

     call linComb_modelParam(a1, mVec1%m, a2, mVec2%m, mVec%m)

     mVec%grid => mVec1%grid
     mVec%allocated = .true.

  end subroutine linComb_modelVec

  !**********************************************************************

  function dotProd_modelVec_f(mVec1,mVec2) result(r)

     real (kind=prec)                   :: r
     type (modelVec_t), intent(in)      :: mVec1,mVec2

     r = dotProd_modelParam_f(mVec1%m, mVec2%m)

  end function dotProd_modelVec_f

  !**********************************************************************
  !  computes mOut = a * mIn for modelVec object m and real scalar a

  function scMult_modelVec_f(a,mIn) result (mOut)

    real (kind=prec), intent(in)		:: a
    type (modelVec_t), intent(in)	    :: mIn
    type (modelVec_t)                   :: mOut

    if (.not. mIn%allocated) then
       call errStop('input not allocated on call to scMult_modelVec_f')
    endif

	call linComb_modelVec(R_ZERO, mIn, a, mIn, mOut)

  end function scMult_modelVec_f

  !**********************************************************************
  !  extracts the grid from a modelVec object; this is only needed since
  !  the attributes are private (in F2003, can declare modelParam private
  !  and grid and allocated attributes could be public)

  subroutine getGrid_modelVec(mIn,grid)

    type (modelVec_t), intent(in)     :: mIn
    type (grid_t), intent(out)        :: grid

    if (.not. mIn%allocated) then
       call warning('model vector not allocated on call to getGrid_modelVec')
    endif

    grid = mIn%grid

  end subroutine getGrid_modelVec

  ! I/O interfaces

  !**********************************************************************
  !  we have to assume that modelParam and grid are in different files,
  !  since they are conceptually different entities;
  !  if they are stored in the same file (e.g. formats of Randie Mackie
  !  or Weerachai Siripunvaraporn), write two separate routines to extract
  !  the model parameter and the grid

  subroutine read_modelVec(mVec, cfile_model, cfile_grid)

    type (modelVec_t), intent(out)	        :: mVec
    character(*), intent(in)		        :: cfile_model
    character(*), intent(in), optional      :: cfile_grid
    ! local
    integer                                 :: status

	if (present(cfile_grid)) then
		call read_grid(mVec%grid, cfile_grid, status)
	else
		call read_grid(mVec%grid, cfile_model, status)
	endif

	call read_modelParam(mVec%m, cfile_model, status)

	if (status == 0) then
		mVec%allocated = .true.
	endif

  end subroutine read_modelVec


  !**********************************************************************
  !  we have to assume that modelParam and grid are in different files,
  !  since they are conceptually different entities;
  !  if they are stored in the same file (e.g. formats of Randie Mackie
  !  or Weerachai Siripunvaraporn), we just write out the model param.
  !  of course, in this case that will have to contain some grid info!

  subroutine write_modelVec(mVec, cfile_model, cfile_grid)

    type (modelVec_t), intent(out)	        :: mVec
    character(*), intent(in)		        :: cfile_model
    character(*), intent(in), optional      :: cfile_grid
    ! local
    integer                                 :: status

	if (present(cfile_grid)) then
		call write_grid(mVec%grid, cfile_grid, status)
	endif

	call write_modelParam(mVec%m, cfile_model, status)

  end subroutine write_modelVec

  !**********************************************************************
  !  temporary read routine for e.g. the full Jacobian; need to come up
  !  with a data type to store a matrix of model parameters

  subroutine readVec_modelVec(mlen, mVec, header, cfile_model, cfile_grid)

    integer, intent(in)                     :: mlen
    type (modelVec_t), pointer, intent(out)	:: mVec(:)
    character(*), intent(out)               :: header
    character(*), intent(in)		        :: cfile_model
    character(*), intent(in), optional      :: cfile_grid
    ! local
    type (modelParam_t), pointer            :: m(:)
    type (grid_t)                           :: grid
    integer                                 :: i, status

    if (.not. associated(mVec)) then
        allocate(mVec(mlen), STAT=status)
    endif

	if (present(cfile_grid)) then
		call read_grid(grid, cfile_grid, status)
	else
		call read_grid(grid, cfile_model, status)
	endif

    allocate(m(mlen), STAT=status)

	call readVec_modelParam(mlen, m, header, cfile_model, status)

    do i = 1,mlen
        mVec(i)%grid = grid
        mVec(i)%m = m(i)
	    if (status == 0) then
		    mVec(i)%allocated = .true.
	    endif
    enddo

    deallocate(m, STAT=status)

  end subroutine readVec_modelVec


  !**********************************************************************
  !  temporary read routine for e.g. the full Jacobian; need to come up
  !  with a data type to store a matrix of model parameters

  subroutine writeVec_modelVec(mlen, mVec, header, cfile_model, cfile_grid)

    integer, intent(in)                     :: mlen
    type (modelVec_t), pointer, intent(out)	:: mVec(:)
    character(*), intent(in)                :: header
    character(*), intent(in)		        :: cfile_model
    character(*), intent(in), optional      :: cfile_grid
    ! local
    type (modelParam_t), pointer            :: m(:)
    integer                                 :: i, status

    if (.not. associated(mVec)) then
        call errStop('the sensitivity matrix is not allocated in writeVec_modelVec')
    elseif (mlen /= size(mVec)) then
        call errStop('the size of sensitivity matrix is not as expected in writeVec_modelVec')
    endif

	if (present(cfile_grid)) then
		call write_grid(mVec(1)%grid, cfile_grid, status)
	endif

    allocate(m(mlen), STAT=status)
    do i = 1,mlen
        m(i) = mVec(i)%m
    enddo

	call writeVec_modelParam(mlen, m, header, cfile_model, status)

    deallocate(m, STAT=status)

  end subroutine writeVec_modelVec


  ! ModelMap interfaces
  !
  ! These can be generalized to work for any problem. For now, they are
  ! not generic and depend on the implementation of ModelMap module.
  ! These mapping will use both the grid and the model parameter in a
  ! modelVec. However, for the 3D MT problem, grid is already part of
  ! the model parameter. For simplicity, leaving it the way it is.

  !**********************************************************************
  ! Non-linear mapping from model parameter to a single grid cell

  function ModelToOneCell(mVec, j, k) result (r)

     type (modelVec_t), intent(in)           :: mVec
     integer, intent(in)                     :: j,k
     real (kind=prec)                        :: r

     r = rhoC(mVec%m,j,k)

  end function ModelToOneCell

  !**********************************************************************
  ! Non-linear mapping from model parameter to grid cell centers

  subroutine ModelToCell(mVec, Ny, Nz, cCond)

     type (modelVec_t), intent(in)           :: mVec
     integer, intent(in)	                 :: Ny,Nz
     real (kind=prec), intent(inout)         :: cCond(Ny,Nz)

     call CondParamToArray(mVec%m,Ny,Nz,cCond)

  end subroutine ModelToCell

  !**********************************************************************
  ! Linear mapping from model parameter to grid edges d\pi/dm

  subroutine dModelToEdge(mVec, eCondY, eCondZ, mVec0)

     type (modelVec_t), intent(in)           :: mVec
     type (cvector), intent(out)             :: eCondY, eCondZ
     type (modelVec_t), intent(in)           :: mVec0

     call CellToEdge(mVec%m, mVec0%m, eCondY, eCondZ)

  end subroutine dModelToEdge

  !**********************************************************************
  ! Linear mapping from grid edges to model parameter (d\pi/dm)^T
  ! (adjoint of dModelToEdge)

  subroutine dEdgeToModel(eCondY, eCondZ, mVec, mVec0)

     type (cvector), intent(in)              :: eCondY, eCondZ
     type (modelVec_t), intent(out)          :: mVec
     type (modelVec_t), intent(in)           :: mVec0

     call EdgeToCell(eCondY, eCondZ, mVec0%m, mVec%m)
     mVec%grid = mVec0%grid
     mVec%allocated = .true.

  end subroutine dEdgeToModel

  !**********************************************************************
  ! Linear mapping from model parameter to grid nodes - a different
  ! version of d\pi/dm

  subroutine dModelToNode(mVec, cCond, mVec0)

     type (modelVec_t), intent(in)           :: mVec
     type (cvector), intent(out)             :: cCond
     type (modelVec_t), intent(in)           :: mVec0

     call CellToNode(mVec%m, cCond, mVec0%m)

  end subroutine dModelToNode

  !**********************************************************************
  ! Linear mapping from grid nodes to model parameter - version of
  ! (d\pi/dm)^T (adjoint of dNodeToModel)

  subroutine dNodeToModel(cCond, mVec, mVec0)

     type (cvector), intent(in)              :: cCond
     type (modelVec_t), intent(out)          :: mVec
     type (modelVec_t), intent(in)           :: mVec0

     call NodeToCell(cCond, mVec%m, mVec0%m)
     mVec%grid = mVec0%grid
     mVec%allocated = .true.

  end subroutine dNodeToModel

  !****************************************************************************
  subroutine dCellToModel(Qj,sigma0,dsigmaReal,dSigmaImag)

   !  given input sparse vector defined on cells,
   !   compute vector defined on parameter space
   !   output result as dsigmaReal, dsigmaImag
   !   i.e. multiplies a sparse vector by (d\pi/dm)^T
   !   used to compute and multiply by Q^T
   !  a wrapper for SparseCelltoModelParam

   type (sparsevecc), intent(in)                  :: Qj
   type (modelVec_t), intent(in)                :: sigma0
   type (modelVec_t), intent(out)             :: dsigmaReal
   type (modelVec_t), intent(out),optional    :: dsigmaImag

   if(present(dSigmaImag)) then
         call SparseCelltoModelParam(Qj,sigma0%m,dsigmaReal%m,dSigmaImag%m)
         dSigmaReal%grid = sigma0%grid
         dSigmaReal%allocated = .true.
         dSigmaImag%grid = sigma0%grid
         dSigmaImag%allocated = .true.
   else
         call SparseCelltoModelParam(Qj,sigma0%m,dsigmaReal%m)
         dSigmaReal%grid = sigma0%grid
         dSigmaReal%allocated = .true.
   endif

  end subroutine dCellToModel

end module ModelSpace
