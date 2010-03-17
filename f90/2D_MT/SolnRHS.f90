module solnrhs
!   higher level module to define EM solution and EMrhs objects
!   plus basic methods, linear algebra, dot products
!
! Defines: EMsoln, EMsparse, EMrhs
! Uses: EMfield

use emfield
use transmitters

implicit none


 type :: EMsoln_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    !!
    !!   However, the specific implementation in this module
    !!    will depend on the specific problem.  This version is for 2D MT.
    !!   the basic solution object
    type(cvector)			:: vec

    !!  Mode (TE or TM) and omega/period are stored in the transmitter
    !!  dictionary. We avoid multiple potential problems by not duplicating
    !!  this information in EM solution (even when we have it available).
    !!  Then, if we need access to mode or frequency, we refer to the dictionary.
    type(grid_t), pointer		:: grid
    integer 				:: tx = 0
    logical                 :: allocated = .false.
  end type EMsoln_t

  type :: EMsolnMTX_t
    !! Generic solution type for storing solutions from multiple transmitters
    integer			:: nTx = 0
    type(EMsoln_t), pointer		:: solns(:)
    logical			:: allocated = .false.
  end type EMsolnMTX_t

  type :: EMsparse_t
    !!   Generic solution type, same name must be used to allow
    !!   use of higher level inversion modules on different problems.
    type(sparsevecc)			:: L
  end type EMsparse_t


  type :: EMrhs_t
     !!   right hand side for solving both TE and TM mode equations,
     !!   forward or adjoint problems (mode & omega stored in txDict)
     character*3			:: adj = ''
     logical           			:: nonzero_source = .false.
     logical				:: nonzero_bc = .false.
     logical				:: allocated = .false.
     type(cvector)			:: source
     complex(kind=prec), pointer, dimension(:)	:: bc
     type(grid_t), pointer		:: grid
     integer                :: tx = 0
  end type EMrhs_t

contains

!**********************************************************************
!           Basic EMsoln methods
!**********************************************************************
     subroutine create_EMsoln(grid,iTx,e)
       ! no need to pass the mode to create_EMsoln:
       ! we now get it from the transmitter dict.
       ! now the interface is generic,
       ! and could be used as such

       implicit none
       type(grid_t), intent(in), target	:: grid
       integer, intent(in)              :: iTx
       type (EMsoln_t), intent(inout)	:: e

       ! local
       character(2)      :: mode
       character(80)     :: gridType

       mode = txDict(iTx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
       else
          call errStop('Unknown mode in create_EMsoln')
       endif

       call create_cvector(grid,gridType,e%vec)
       e%tx = iTx
       e%grid => grid
       e%allocated = .true.

     end subroutine create_EMsoln

     !************************************************************
     subroutine deall_EMsoln(e)
       implicit none
       type (EMsoln_t), intent(inout)   :: e

       call deall_cvector(e%vec)
       if(associated(e%grid)) then
           nullify(e%grid)
       endif
       e%allocated = .false.

     end subroutine deall_EMsoln

     !************************************************************
     subroutine copy_EMsoln(eOut,eIn)

       implicit none
       type (EMsoln_t), intent(in)	:: eIn
       type (EMsoln_t), intent(inout)	:: eOut

       !  should have some error checking for eIn ...
       call copy_cvector(eOut%vec,eIn%vec)

       eOut%tx = eIn%tx
       eOut%grid => eIn%grid
       eOut%allocated = eIn%allocated

     end subroutine copy_EMsoln

     !**********************************************************************
     function dotProd_EMsoln(e1,e2,conj_Case) result(c)
     !   dot product of two solution space objects
       type(EMsoln_t), intent(in)		:: e1,e2
       logical, intent(in)		:: conj_case
       complex(kind=prec)	:: c

       c = dotProd_cvector(e1%vec,e2%vec,conj_case)

     end function dotProd_EMsoln

     !**********************************************************************
     subroutine zero_EMsoln(e)
     !  zeros a solution space object

       type(EMsoln_t), intent(inout)	:: e

       e%vec%v = C_ZERO

     end subroutine zero_EMsoln


!**********************************************************************
!           Basic EMsolnMTX methods
!**********************************************************************

   subroutine create_EMsolnMTX(nTx,eAll)

      integer, intent(in)               :: nTx
      type(EMsolnMTX_t), intent(inout)  :: eAll

      ! local variables
      integer      :: istat

      call deall_EMsolnMTX(eAll)

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

!**********************************************************************
!           Basic EMsparse methods
!**********************************************************************

     subroutine create_EMsparse(grid,iTx,LC,nCoeff)

       type(grid_t), intent(in)	:: grid
       integer, intent(in)		:: iTx
       type(EMsparse_t)			:: LC
       integer, optional, intent(in)	:: nCoeff

       ! local
       character(2)      :: mode
       character(80)     :: gridType
       integer			 :: nc

       if (present(nCoeff)) then
          nc = nCoeff
       else
          nc = 0 ! will reallocate memory later in the program
       end if

       mode = txDict(iTx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
       else
          call errStop('Unknown mode in create_EMsparse')
       endif

       call create_sparsevecc(grid,gridType,nc,LC%L)

     end subroutine create_EMsparse

     !************************************************************
     subroutine deall_EMsparse(LC)

       type (EMsparse_t), intent(inout)   	:: LC

       call deall_sparsevecc(LC%L)

     end subroutine deall_EMsparse

     !**********************************************************************
     subroutine linComb_EMsparse(Lin1,c1,Lin2,c2,Lout)
       ! linear combination of two EMsparse objects

       type (EMsparse_t), intent(in)		:: Lin1,Lin2
       complex (kind=prec), intent(in)	:: c1,c2
       type (EMsparse_t), intent(inout)		:: Lout

       call linComb_sparsevecc(Lin1%L,c1,Lin2%L,c2,Lout%L)

     end subroutine linComb_EMsparse

!**********************************************************************
!           combined EMsoln/EMsparse methods
!**********************************************************************
     function dotProd_EMsparseEMsoln(SV,FV,Conj_Case) result(c)

       type (EMsparse_t), intent(in)             :: SV  ! sparse vector
       type (EMsoln_t), intent(in)               :: FV  ! full vector
       logical, intent(in)                     :: conj_Case ! = .true.
       complex(kind=prec)		:: c

       c = dotProd_scvector(SV%L,FV%vec,Conj_Case)

     end function dotProd_EMsparseEMsoln

     !**********************************************************************
     subroutine add_EMsparseEMsoln(cs,SV,FV)

        complex(kind=prec), intent(in)	:: cs
        type (EMsparse_t), intent(in)	:: SV  ! sparse vector
        type (EMsoln_t), intent(inout)	:: FV  ! full vector

       call add_scvector(cs,SV%L,FV%vec)

     end subroutine add_EMsparseEMsoln

!**********************************************************************
!           EMrhs methods
!**********************************************************************
     !  allocates and initializes arrays for the "rhs" structure
     !   set pointer to grid + mode/source fields of rhs before calling
     subroutine create_EMrhs(grid,iTx,b)

       type(grid_t), intent(in),target	:: grid
       integer, intent(in)              :: iTx
       type (EMrhs_t), intent(inout)   	:: b

       !  local variables
       character(2)     :: mode
       character(80)    :: gridType
       integer ::       Nz,Ny,Nza,Nzi,nBC

       Nz = grid%Nz
       Ny = grid%Ny
       Nza = grid%Nza

       mode = txDict(iTx)%mode

       if(mode .eq. 'TE') then
          gridType = NODE
          Nzi = Nz
       else if(mode .eq. 'TM') then
          gridType = NODE_EARTH
          Nzi = Nz-Nza
       else
          call errStop('Unknown mode in create_EMrhs')
       endif

       nBC = 2*(Ny+1)+2*(Nzi+1)

       if (b%nonzero_BC) then
          ! Initialize array for boundary condition data
          allocate(b%BC(nBC))
       endif

       if (b%nonzero_source) then
           ! Initialize full array for storage of source
           call create_cvector(grid,gridType,b%source)
       endif

       b%tx = iTx
       b%grid => grid
       b%allocated = .true.

     end subroutine create_EMrhs

    !************************************************************
     subroutine deall_EMrhs(b)

       type (EMrhs_t), intent(inout)   :: b

       if (.not.(b%allocated)) then
         return
       endif

       if (b%nonzero_BC) then
          deallocate(b%BC)
          nullify(b%BC)
       endif

       if (b%nonzero_source) then
          call deall_cvector(b%source)
       endif
       b%allocated = .false.

     end subroutine deall_EMrhs

     !**********************************************************************
     subroutine add_EMsparseEMrhs(cs,SV,FV)

        complex(kind=prec), intent(in)	:: cs
        type (EMsparse_t), intent(in)		:: SV  ! sparse vector
        type (EMrhs_t), intent(inout)		:: FV  ! full vector

       call add_scvector(cs,SV%L,FV%source)
       FV%nonzero_source = .true.

     end subroutine add_EMsparseEMrhs

     !**********************************************************************
     subroutine zero_EMrhs(e)
     !  zeros a solution space object

       type(EMrhs_t), intent(inout)	:: e

       if(e%nonzero_source .and. e%allocated) then
          e%source%v = C_ZERO
       else
          if(e%nonzero_bc .and. e%allocated) then
             e%bc = C_ZERO
          endif
       endif
     end subroutine zero_EMrhs

end module solnrhs
