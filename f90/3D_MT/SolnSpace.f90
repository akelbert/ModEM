module SolnSpace
!   higher level module to define EM solution and RHS objects
!   plus basic methods, linear algebra, dot products
!   3D MT version
!
! Defines: EMsoln, EMsparse, EMrhs
! Uses: EMfield

use math_constants
use utilities
use sg_vector
use sg_boundary
use sg_sparse_vector
use transmitters

implicit none

interface assignment (=)
   MODULE PROCEDURE copy_EMsoln
   MODULE PROCEDURE copy_EMsolnMTX
   MODULE PROCEDURE copy_EMsparse
   !MODULE PROCEDURE copy_EMrhs - doesn't exist yet
end interface

  type :: EMsoln_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!
    !!   However, the specific implementation in this module
    !!    will depend on the specific problem.  This version is for 3D MT,
    !!  and collects solutions for both "polarizations" in a single
    !!   structure.  Thus the impedance operator acts on an object of this type.

    !! e is array of electric field solution vectors for nPol source
    !! polarizations; 2 for 3D MT - one electrical field solution for each
    !! allows a variable number of polarizations to accommodate other uses
    !! e.g. active source applications
    integer					:: nPol = 2
    type(cvector), pointer  :: pol(:)

    !! tx points to information in the transmitter dictionary about the source
    !!   used to compute the solution, e.g. omega/period;
    !!   do not duplicate it here to avoid potential problems
    integer 			:: tx = 0

    !! grid is a pointer to numerical discretization stored in SensMatrix
    type(grid_t), pointer	:: grid

    !! allocated when the EMsoln was created but not yet deallocated
    logical			:: allocated = .false.

    !! avoid memory leaks: set this to true for function outputs only
    logical			:: temporary = .false.

  end type EMsoln_t

  type :: EMsolnMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer						:: nTx = 0
    type(EMsoln_t), pointer		:: solns(:)
    logical						:: allocated = .false.
    logical						:: temporary = .false.
  end type EMsolnMTX_t

  type :: EMsparse_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!   Here we need nPol sparse vectors, one for each polarization
    !!    to represent linear data functionals on objects of type EMsoln
    !!  this doesn't really need the grid... not used anywhere at the moment
    !!  but then the generic create interface shouldn't use it either
    !!  unless sparse vectors require the grid pointer for something.
    integer						:: nPol = 2
    integer						:: nCoeff = 0
    type(sparsevecc), pointer	:: L(:)
    logical						:: allocated = .false.
    logical						:: temporary = .false.
    type(grid_t), pointer       :: grid
    integer						:: tx = 0
  end type EMsparse_t

  type :: RHS_t
     ! merges internal sources and boundary conditions into a single
     ! data structure; this is a compact representation of the right
     !  hand side for the induction equations for ONE mode

     character*3		:: adj = ''
     character*2		:: polName ! Ex or Ey
     logical                    :: nonzero_BC     = .false.
     logical                    :: nonzero_Source = .false.
     logical                    :: sparse_Source  = .false.
     logical			:: allocated      = .false.
    logical					:: temporary = .false.
     type (cvector) 		:: s
     type (sparsevecc) 		:: sSparse
     type (cboundary) 		:: bc
     type(grid_t), pointer	:: grid
     integer                :: tx = 0
  end type RHS_t

  type :: EMrhs_t
     ! rhs data structure for multiple polarizations, the abstract
     !  full rhs used for the abstract full EMsoln
     ! pointer to the grid needed for cleaner routines in this module.

     integer				:: nPol = 2
     type (RHS_t), pointer	:: b(:)
     logical				:: allocated = .false.
     logical				:: temporary = .false.
     type(grid_t), pointer	:: grid
     integer				:: tx = 0
  end type EMrhs_t

contains

!**********************************************************************
!           Basic EMsoln methods
!**********************************************************************

     subroutine create_EMsoln(grid,iTx,e)

     !  generic routine for creating the EMsoln type for 3D problems:
     !  number of polarization obtained from the transmitter dictionary

       implicit none
       type(grid_t), intent(in), target	    :: grid
       integer, intent(in)                  :: iTx
       type (EMsoln_t), intent(inout)		:: e

       ! local variables
       integer				:: k,istat

       if (e%allocated) then
          if (associated(e%grid, target=grid) .and. (e%tx == iTx)) then
             ! do nothing
             return
          else
             call deall_EMsoln(e)
          end if
       end if

       e%nPol = txDict(iTx)%nPol
       allocate(e%pol(e%nPol), STAT=istat)
       do k = 1,e%nPol
          call create_cvector(grid,e%pol(k),EDGE)
       enddo
       e%tx = iTx
       e%grid => grid

	   e%allocated = .true.

     end subroutine create_EMsoln

     !************************************************************
     subroutine deall_EMsoln(e)

       !  3D  version
       implicit none
       type (EMsoln_t), intent(inout)   :: e

       ! local variables
       integer				:: k, istat

       if (associated(e%pol)) then
          do k = 1,e%nPol
             call deall_cvector(e%pol(k))
          enddo
          deallocate(e%pol, STAT=istat)
       endif

       if(associated(e%grid)) then
           nullify(e%grid)
       endif

       e%allocated = .false.

     end subroutine deall_EMsoln

     !************************************************************
     subroutine copy_EMsoln(eOut,eIn)

       !  3D  version
       implicit none
       type (EMsoln_t), intent(in)	:: eIn
       type (EMsoln_t), intent(inout)	:: eOut

       ! local variables
       integer				:: k

       if (.not. eIn%allocated) then
         call errStop('input EM soln not allocated yet in copy_EMsoln')
       endif

       call create_EMsoln(eIn%grid,eIn%tx,eOut)

       do k = 1,eIn%nPol
          call copy_cvector(eOut%pol(k),eIn%pol(k))
       enddo

       eOut%allocated = eIn%allocated

       !if (eIn%temporary) then
       !   call deall_EMsoln(eIn)
       !endif

     end subroutine copy_EMsoln

     !**********************************************************************
     subroutine zero_EMsoln(e)
     !  zeros a solution space object

       type(EMsoln_t), intent(inout)	:: e

       ! local variables
       integer				:: k

       do k = 1,e%nPol
          call zero_cvector(e%pol(k))
       enddo

     end subroutine zero_EMsoln


!**********************************************************************
!           Basic EMsolnMTX methods
!**********************************************************************

   subroutine create_EMsolnMTX(nTx,eAll)

      integer, intent(in)               :: nTx
      type(EMsolnMTX_t), intent(inout)  :: eAll

      !  local variables
      integer                           :: istat

      if (eAll%allocated) then
         if (eAll%nTx == nTx) then
            ! do nothing
            return
         else
            call deall_EMsolnMTX(eAll)
         end if
      end if

      eAll%nTx = nTx
      allocate(eAll%solns(nTx), STAT=istat)
      eAll%allocated = .true.

   end subroutine create_EMsolnMTX

   !**********************************************************************
   subroutine deall_EMsolnMTX(eAll)

      type(EMsolnMTX_t), intent(inout)     :: eAll

      !  local variables
      integer                           :: j, istat

	  do j = 1,eAll%nTx
	  	call deall_EMsoln(eAll%solns(j))
	  end do

      if (associated(eAll%solns)) deallocate(eAll%solns, STAT=istat)
      eAll%allocated = .false.

   end subroutine deall_EMsolnMTX

   !************************************************************
   subroutine copy_EMsolnMTX(eOut,eIn)

       !  3D  version
       implicit none
       type (EMsolnMTX_t), intent(in)	:: eIn
       type (EMsolnMTX_t), intent(inout)	:: eOut

       ! local variables
       integer				:: j

       if (.not. eIn%allocated) then
         call errStop('input multi-transmitter EM soln not allocated yet in copy_EMsolnMTX')
       endif

       call create_EMsolnMTX(eIn%nTx,eOut)

       do j = 1,eIn%nTx
          call copy_EMsoln(eOut%solns(j),eIn%solns(j))
       enddo

       eOut%allocated = eIn%allocated

       !if (eIn%temporary) then
       !   call deall_EMsolnMTX(eIn)
       !endif

   end subroutine copy_EMsolnMTX

!**********************************************************************
!           Basic EMsparse methods
!**********************************************************************
!    don't really use linear combinations ... so not implemented

    subroutine create_EMsparse(grid,iTx,LC,nCoeff)

      !  generic routine for creating the EMsparse type for 3D problems:
      !  number of polarization obtained from the transmitter dictionary

       implicit none
       type(grid_t), intent(in), target	    :: grid
       integer, intent(in)                  :: iTx
       type (EMsparse_t), intent(inout)		:: LC
       integer, intent(in), optional		:: nCoeff

       ! local variables
       integer				:: nc,k,istat

       if (present(nCoeff)) then
          nc = nCoeff
       else
          nc = 0 ! will reallocate memory later in the program
       end if

       if (LC%allocated) then
          ! it is safest and quite efficient to reallocate sparse vectors
          call deall_EMsparse(LC)
       end if

       LC%nPol = txDict(iTx)%nPol
       allocate(LC%L(LC%nPol), STAT=istat)
       do k = 1,LC%nPol
          call create_sparsevecc(nc,LC%L(k),EDGE)
       enddo

       LC%nCoeff = nc
       LC%grid => grid
       LC%tx = iTx
	   LC%allocated = .true.

    end subroutine create_EMsparse

    !************************************************************
    subroutine deall_EMsparse(LC)

      ! 3D version
      type (EMsparse_t), intent(inout)   	:: LC

      ! local variables
      integer				:: k,istat

      if (associated(LC%L)) then
         do k = 1,LC%nPol
            call deall_sparsevecc(LC%L(k))
         enddo
         deallocate(LC%L, STAT=istat)
         nullify(LC%grid)
      endif

   end subroutine deall_EMsparse

   !************************************************************
   subroutine copy_EMsparse(eOut,eIn)

       !  3D  version
       implicit none
       type (EMsparse_t), intent(in)	:: eIn
       type (EMsparse_t), intent(inout)	:: eOut

       ! local variables
       integer				:: k

       if (.not. eIn%allocated) then
         call errStop('input EM sparse vector not allocated yet in copy_EMsparse')
       endif

       call create_EMsparse(eIn%grid,eIn%tx,eOut,eIn%nCoeff)

       do k = 1,eIn%nPol
          call copy_sparsevecc(eOut%L(k),eIn%L(k))
       enddo

       eOut%allocated = eIn%allocated

       !if (eIn%temporary) then
       !   call deall_EMsparse(eIn)
       !endif

   end subroutine copy_EMsparse

!**********************************************************************
!           combined EMsoln/EMsparse methods
!**********************************************************************
!  subroutine add_EMsparseEMsoln(cs,SV,FV)  does not appear to be needed

   function dotProd_EMsparseEMsoln(SV,FV,Conj_Case) result(c)

       type (EMsparse_t), intent(in)             :: SV  ! sparse vector
       type (EMsoln_t), intent(in)               :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)		:: c
       integer					:: k

       c = C_ZERO
       if(conj_case) then
          do k = 1,SV%nPol
             c = c + dotProd_scvector_f(SV%L(k),FV%pol(k))
          enddo
       else
          do k = 1,SV%nPol
             c = c + dotProd_noConj_scvector_f(SV%L(k),FV%pol(k))
          enddo
       endif

   end function dotProd_EMsparseEMsoln


!**********************************************************************
!           RHS methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid before calling
     !   NO gridType needed for 3DMT

     subroutine create_RHS(grid,iTx,b)
     !   3D version  ...
     !     does not create sparse vectors if sparsesource = .true.
     !       (this would in any event require knowing number of
     !         non-zero coefficients to allow for)
		 ! NOTE: Does not do anything if b%allocated is .true.

       type(grid_t), intent(in),target	:: grid
       integer, intent(in)              :: iTx
       type (RHS_t), intent(inout)   	:: b

			 if (b%allocated) then
					! do nothing - exit the create subroutine
					return
			 endif

       if (b%nonzero_BC) then
          ! create boundary condition data structures for each polarization
          !  NOTE: get rid of "used for" in BC type def
          call create_cboundary(grid,b%bc)
          b%allocated = .true.
       endif

       if (b%nonzero_source) then
          if(b%sparse_source) then
            !   can't create sparse vector without knowing how
            !    many components it has; nothing to do until we
            !    actually use (e.g., add another sparse vector)
          else
              ! Initialize full array for storage of source
              call create_cvector(grid, b%s, EDGE)
              b%allocated = .true.
          endif
       endif

       b%tx = iTx
       b%grid => grid

     end subroutine create_RHS

    !************************************************************
     subroutine deall_RHS(b)

       type (RHS_t), intent(inout)   :: b

       if (b%nonzero_BC) then
          call deall_cboundary(b%bc)

       endif

       if (b%nonzero_source) then
          if(b%sparse_source) then
             call deall_sparsevecc(b%sSparse)
          else
             call deall_cvector(b%s)
          endif
       endif

       nullify(b%grid)
       b%allocated = .false.

     end subroutine deall_RHS

     !**********************************************************************


     !**********************************************************************
     subroutine zero_RHS(b)
     !  zeros a RHS object

       type(RHS_t), intent(inout)	:: b

       if(b%nonzero_source .and. b%allocated) then
          if(b%sparse_source) then
             !  if sparse vector is zeroed, all components are
             !    deleted ...
             call deall_sparsevecc(b%sSparse)
          else
             call zero_cvector(b%s)
          endif
       else if(b%nonzero_bc .and. b%allocated) then
          call zero_cboundary(b%bc)
       else
          if(.not.b%allocated) then
             call errStop('Input not yet allocated in zero_RHS')
          endif
       endif
     end subroutine zero_RHS


!**********************************************************************
!           EMrhs methods
!**********************************************************************

     subroutine create_EMrhs(grid,iTx,b)
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid before calling
     !   3D version  ...
     !     does not create sparse vectors if sparsesource = .true.
     !       (this would in any event require knowing number of
     !         non-zero coefficients to allow for)

       type(grid_t), intent(in),target	:: grid
       integer, intent(in)              :: iTx
       type (EMrhs_t), intent(inout)   	:: b

       integer				:: k,istat

       if (b%allocated) then
          if (associated(b%grid, target=grid) .and. (b%tx == iTx)) then
             ! do nothing
             return
          else
             call deall_EMrhs(b)
          end if
       end if

       b%nPol = txDict(iTx)%nPol
       allocate(b%b(b%nPol), STAT=istat)
       do k = 1,b%nPol
          call create_RHS(grid,iTx,b%b(k))
       enddo

       b%allocated = .true.

     end subroutine create_EMrhs

    !************************************************************
     subroutine deall_EMrhs(b)

       type (EMrhs_t), intent(inout)   :: b

       integer			:: k,istat

       if (associated(b%b)) then
          do k = 1,b%nPol
             call deall_RHS(b%b(k))
          enddo
          deallocate(b%b, STAT=istat)
       endif

       b%allocated = .false.

     end subroutine deall_EMrhs

!     !************************************************************
!     ! need copy_RHS for this... too complicated, will write this
!     ! when it is needed
!     subroutine copy_EMrhs(bOut,bIn)
!
!       !  3D  version
!       implicit none
!       type (EMrhs_t), intent(in)	:: bIn
!       type (EMrhs_t), intent(inout)	:: bOut
!
!       ! local variables
!       integer				:: k
!
!       if (.not. bIn%allocated) then
!         call errStop('input EM RHS not allocated yet in copy_EMrhs')
!       endif
!
!       call create_EMrhs(bIn%grid,bIn%tx,bOut)
!
!       do k = 1,bIn%nPol
!          call copy_RHS(bOut%b(k),bIn%b(k))
!       enddo
!
!       bOut%allocated = bIn%allocated
!
!       if (bIn%temporary) then
!          call deall_EMrhs(bIn)
!       endif
!
!     end subroutine copy_EMrhs

     !**********************************************************************
     subroutine zero_EMrhs(b)
     !  zeros a EMrhs object

       type(EMrhs_t), intent(inout)	:: b

       integer			:: k

       do k = 1,b%nPol
          call zero_RHS(b%b(k))
       enddo
     end subroutine zero_EMrhs

     !**********************************************************************

     subroutine add_EMsparseEMrhs(cs,SV,comb)

     !   Forms the linear combination cs*SV+comb where:
     !     cs is complex
     !     SV is an EMsparse object
     !     comb is an EMrhs object
     !   Result is returned in comb
     !
     !   In this implementation, an EMrhs object

       complex(kind=prec), intent(in)  :: cs
       type (EMsparse_t), intent(in)             :: SV  ! sparse vector
       type (EMrhs_t), intent(inout)             :: comb  ! full vector

     !  local variables
       type(sparsevecc)			:: temp
       integer				:: k

       do k = 1,comb%nPol
          comb%b(k)%nonzero_source = .true.
          if(comb%b(k)%sparse_Source) then
             ! use sparse vector storage for output
             if(comb%b(k)%sSparse%allocated) then
                !  add to contentes of comb sparse vector and cs*SV
                call linComb_sparsevecc(comb%b(k)%sSparse,C_ONE, &
			SV%L(k),cs,temp)
             else
                call scMult_sparsevecc(cs,SV%L(k),temp)
             endif
             call copy_sparsevecc(temp,comb%b(k)%sSparse)
          else
             !  might want to check for allocation of comb%b(k)%s first
             call add_scvector(cs,SV%L(k),comb%b(k)%s)
          endif
       enddo

       call deall_sparsevecc(temp)

     end subroutine add_EMsparseEMrhs

end module SolnSpace
