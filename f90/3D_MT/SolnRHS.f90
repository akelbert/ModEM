module solnrhs
!   higher level module to define EM solution and RHS objects
!   plus basic methods, linear algebra, dot products
!   3D MT version
!
! Defines: EMsoln, EMsparse, EMrhs
! Uses: EMfield, ModelParameter

use math_constants
use utilities
use modelparameter
use sg_vector
use sg_boundary
use sg_sparse_vector

implicit none

  type :: EMsoln
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!
    !!   However, the specific implementation in this module
    !!    will depend on the specific problem.  This version is for 3D MT,
    !!  and collects solutions for both "polarizations" in a single
    !!   structure.  Thus the impedance operator acts on an object of this type.

    !! e is array of electric field solution vectors for two source
    !! polarizations; one electrical field solution for each
    !!  might generalize to allow a variable number of polarizations
    integer			:: nPol = 2
    type(cvector), dimension(2)  	:: pol

    !! omega, period, tx are information about the source used to compute
    !!   the solution
    real(kind=selectedPrec)	:: omega = R_ZERO
    real(kind=selectedPrec)	:: period = R_ZERO
    integer 			:: tx = 0

    !! sigma, grid are pointers to the model parameter and grid used 
    type(modelParam_t), pointer	:: sigma
    type(grid3d_t), pointer	:: grid

		!! allocated when the EMsoln was created but not yet deallocated
    logical			:: allocated = .false.

  end type EMsoln

  type :: EMsparse
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!   Here we need two sparse vectors, one for each polarization
    !!    to represent linear data functionals on objects of type EMsoln
    integer			:: nPol = 2
    type(sparseVecC)		:: L(2)
  end type EMsparse

  type :: RHS
     ! merges internal sources and boundary conditions into a single
     ! data structure; this is a compact representation of the right
     !  hand side for the induction equations for ONE mode

     character*3		:: adj = ''
     character*2		:: polName ! Ex or Ey
     logical                    :: nonzero_BC
     logical                    :: nonzero_Source = .false.
     logical                    :: sparse_Source = .false.
     logical			:: allocated = .false.
     type (cvector) 		:: s
     type (sparsevecc) 		:: sSparse
     type (cboundary) 		:: bc
     type(grid3d_t), pointer	:: grid
  end type RHS

  type :: EMrhs
     ! rhs data structure for multiple polarizations, the abstract
     !  full rhs used for the abstract full EMsoln

     integer			:: nPol = 2
     type (RHS) 		:: b(2)
     logical			:: allocated = .false.
  end type EMrhs

contains

!**********************************************************************
!           Basic EMsoln methods
!**********************************************************************

     subroutine create_EMsoln(grid,e)

     !  3DMT  version:  NO gridType needed
     !  does not set transmitter or pointer to conductivity
       implicit none
       type(grid3d_t), intent(in), target	:: grid
       type (EMsoln), intent(inout)		:: e
      
       ! local variables
       integer				:: k

			 if (e%allocated) then
				! do nothing
				return
			 end if

       e%nPol = 2
       do k = 1,e%nPol 
          call create_cvector(grid,e%pol(k),EDGE)
       enddo
       e%grid => grid

			 e%allocated = .true.

     end subroutine create_EMsoln

     !************************************************************
     subroutine deall_EMsoln(e)
     !  3D-MT  version
       implicit none
       type (EMsoln), intent(inout)   :: e

       ! local variables
       integer				:: k

       do k = 1,e%nPol 
          call deall_cvector(e%pol(k))
       enddo

       if(associated(e%grid)) then
           nullify(e%grid)
       endif
       if(associated(e%sigma)) then
           nullify(e%sigma)
       endif

			 e%allocated = .false.

     end subroutine deall_EMsoln

     !************************************************************
     subroutine copy_EMsoln(eOut,eIn)
     !  3D-MT  version
       implicit none
       type (EMsoln), intent(in)	:: eIn
       type (EMsoln), intent(inout)	:: eOut
       
       ! local variables
       integer				:: k
       !  should have some error checking for eIn ...
       do k = 1,eIn%nPol 
          call copy_cvector(eOut%pol(k),eIn%pol(k))
       enddo
       
       eOut%omega = eIn%omega
       eOut%period = eIn%period
       eOut%tx = eIn%tx
       eOut%sigma => eIn%sigma
       eOut%grid => eIn%grid

     end subroutine copy_EMsoln

     !**********************************************************************
     subroutine zero_EMsoln(e)
     !  zeros a solution space object
   
       type(EMsoln), intent(inout)	:: e

       ! local variables
       integer				:: k

       do k = 1,e%nPol 
          call zero_cvector(e%pol(k))
       enddo

     end subroutine zero_EMsoln

!**********************************************************************
!           Basic EMsparse methods
!**********************************************************************
!  NOTE:  create_EMsparse appears to be unnecessary!
!     also don't really use linear combinations ...

     !************************************************************
     subroutine deall_EMsparse(LC)

       ! 3D version
       type (EMsparse), intent(inout)   	:: LC

       ! local variables
       integer				:: k
       
       LC%nPol = 2
       do k = 1,LC%nPol 
          call deall_sparsevecc(LC%L(k))
       enddo

     end subroutine deall_EMsparse
!**********************************************************************
!           combined EMsoln/EMsparse methods
!**********************************************************************
     function dotProd_EMsparseEMsoln(SV,FV,Conj_Case) result(c)

       type (EMsparse), intent(in)             :: SV  ! sparse vector
       type (EMsoln), intent(in)               :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=selectedPrec)		:: c
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

     !  subroutine add_EMsparseEMsoln(cs,SV,FV)  does not appear to be needed

!**********************************************************************
!           RHS methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid before calling
     !   NO gridType needed for 3DMT

     subroutine create_RHS(grid,b)
     !   3D version  ... 
     !     does not create sparse vectors if sparsesource = .true.
     !       (this would in any event require knowing number of
     !         non-zero coefficients to allow for)  
		 ! NOTE: Does not do anything if b%allocated is .true.    
     
       type(grid3d_t), intent(in),target	:: grid
       type (rhs), intent(inout)   	:: b

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

       b%grid => grid

     end subroutine create_RHS

    !************************************************************
     subroutine deall_RHS(b)

       type (rhs), intent(inout)   :: b

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
   
       type(RHS), intent(inout)	:: b

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
             msg = 'Input not yet allocated in zero_RHS'
             call errStop(msg)
          endif
       endif
     end subroutine zero_RHS
!**********************************************************************
!           EMrhs methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid before calling
     !   NO gridType needed for 3DMT

     subroutine create_EMrhs(grid,b)
     !   3D version  ... 
     !     does not create sparse vectors if sparsesource = .true.
     !       (this would in any event require knowing number of
     !         non-zero coefficients to allow for)      
     
       type(grid3d_t), intent(in),target	:: grid
       type (EMrhs), intent(inout)   	:: b

       integer				:: k

       if (b%allocated) then
          ! do nothing - exit the create subroutine
          return
       endif

       do k = 1,b%nPol
          call create_RHS(grid,b%b(k)) 
       enddo

       b%allocated = .true.

     end subroutine create_EMrhs

    !************************************************************
     subroutine deall_EMrhs(b)

       type (EMrhs), intent(inout)   :: b

       integer			:: k

       do k = 1,b%nPol
             call deall_RHS(b%b(k))
       enddo

       b%allocated = .false.

     end subroutine deall_EMrhs

     !**********************************************************************

     subroutine add_EMsparseEMrhs(cs,SV,comb)

     !   Forms the linear combination cs*SV+comb where:
     !     cs is complex
     !     SV is an EMsparse object
     !     comb is an EMrhs object
     !   Result is returned in comb
     !
     !   In this implementation, an EMrhs object

       complex(kind=selectedPrec), intent(in)  :: cs
       type (EMsparse), intent(in)             :: SV  ! sparse vector
       type (EMrhs), intent(inout)             :: comb  ! full vector

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

     end subroutine add_EMsparseEMrhs

     !**********************************************************************
     subroutine zero_EMrhs(b)
     !  zeros a EMrhs object
   
       type(EMrhs), intent(inout)	:: b

       integer			:: k

       do k = 1,b%nPol
          call zero_RHS(b%b(k))
       enddo
     end subroutine zero_EMrhs

end module solnrhs
