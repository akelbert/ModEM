! *****************************************************************************
module sg_scalar
  ! This module creates data types for scalar fields defined on
  ! scalar (center or corner) nodes of a staggered Cartesian grid, along with 
  ! basic algebraic (scalar) operations. Such basic operations are: allocation,
  ! deallocation, initialization, copying, and algebraic operations (linear 
  ! combinations, scalar products, dot products).
  ! Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes) modules.

  use math_constants		! math/ physics constants
  use grid3d
  implicit none

  ! Generic interfaces are done through subroutines
  ! creates scalar (center or corner) nodes
  INTERFACE create
     module procedure create_rscalar
     module procedure create_cscalar
  END INTERFACE

  ! deallocates the scalar (center or corner) nodes
  INTERFACE deall
     module procedure deall_rscalar
     module procedure deall_cscalar
  END INTERFACE

  ! a real scalar multiplies the scalar (center or corner) nodes
  INTERFACE scMult
     module procedure scMult_rscalar
     module procedure scMult_cscalar
  END INTERFACE

  INTERFACE scMultAdd
     module procedure scMult_cscalar
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_cscalar
  END INTERFACE

  ! adds the scalar (center or corner) nodes
  INTERFACE add
     module procedure add_rscalar
     module procedure add_cscalar
  END INTERFACE

  ! subtracts the scalar (center or corner) nodes
  INTERFACE subtract
     module procedure subtract_rscalar
     module procedure subtract_cscalar
  END INTERFACE

  ! pointwise scalar multiplication of scalar (center or corner) nodes and
  ! pointwise complex-real (mixed) multiplication of scalar (center or corner)
  ! nodes
  ! Both are scalar data types
  INTERFACE diagMult
     module procedure diagMult_rscalar
     module procedure diagMult_cscalar
     module procedure diagMult_rcscalar
     module procedure diagMult_crscalar
  END INTERFACE

  ! zeros the scalar (center or corner) nodes
  INTERFACE zero
     module procedure zero_rscalar
     module procedure zero_cscalar
  END INTERFACE

  ! computes the dot product
  INTERFACE dotProd
     MODULE PROCEDURE dotProd_rscalar_f
     MODULE PROCEDURE dotProd_cscalar_f
  END INTERFACE


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_cscalar
     MODULE PROCEDURE copy_rscalar
  END INTERFACE

  ! Interface operators are done through functions
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_rscalar_f
     MODULE PROCEDURE add_cscalar_f
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtract_rscalar_f
     MODULE PROCEDURE subtract_cscalar_f
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE scMult_rscalar_f
     MODULE PROCEDURE scMult_cscalar_f
     MODULE PROCEDURE diagMult_rscalar_f
     MODULE PROCEDURE diagMult_cscalar_f
     MODULE PROCEDURE diagMult_rcscalar_f
     MODULE PROCEDURE diagMult_crscalar_f 
  END INTERFACE


  public		::   create_rscalar,  create_cscalar, &
       deall_rscalar, deall_cscalar, &
       copy_rscalar, copy_cscalar, &
       zero_rscalar, zero_cscalar, &
       scMult_rscalar, scMult_cscalar, &
       scMult_rscalar_f, scMult_cscalar_f, &
       add_rscalar, add_cscalar, &
       add_rscalar_f, add_cscalar_f, &
       subtract_rscalar, subtract_cscalar, &
       subtract_rscalar_f, subtract_cscalar_f, &
       diagMult_rscalar, diagMult_cscalar, &
       diagMult_rscalar_f, diagMult_cscalar_f, &
       diagMult_rcscalar, diagMult_crscalar, &
       diagMult_rcscalar_f, diagMult_crscalar_f, &
       dotProd_rscalar_f, dotProd_cscalar_f, &
       ! the two routines below are not part of any interface
  linComb_cscalar, scMultAdd_cscalar, CornerVolume

  ! ***************************************************************************
  ! type scalar defines scalar for either edge or face in a staggered grid as
  ! a complex field
  type :: cscalar

     ! store the intention of the use as a character string: 'Center', 
     ! CORNER, CELL_EARTH
     character (len=80)	                               :: gridType

     ! Note that the arrays are defined through dynamic memory allocation  
     complex(kind=selectedPrec), pointer, dimension(:,:,:)    :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction: 
     integer                                           :: nx = 0, ny = 0, nz = 0

     ! allocated:  .true.  v array has been allocated
     logical		                               :: allocated = .false.

     ! pointer to parent grid
     type (grid3d_t), pointer                              :: grid
     
  end type cscalar


  ! ***************************************************************************
  ! type scalar defines scalar for either edge or face in a staggered grid as
  ! a real field
  type :: rscalar

     ! store the intention of the use as a character string: 'Center'
     ! CORNER, CELL_EARTH
     character (len=80)	                           	:: gridType

     ! Typical usage:  conductivity averaged on centers of
     ! staggered grid
     ! v: dimension Nx, Ny, Nz
     ! Note that the arrays are defined through dynamic memory allocation  
     real(kind=selectedPrec), pointer, dimension(:,:,:)        :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction: 
     integer                                            :: nx = 0, ny = 0, nz =0

     ! allocated:  .true.  v array has been allocated
     logical		                                :: allocated = .false.

     ! pointer to parent grid
     type (grid3d_t), pointer                               :: grid

  end type rscalar


Contains
  ! CREATE GRID_scalar
  ! * subroutine create_rscalar(igrid, E, gridType)
  ! * subroutine create_cscalar(igrid, E, gridType)

  ! DEALLOCATE GRID_scalar
  ! * subroutine deall_rscalar(E)
  ! * subroutine deall_cscalar(E)

  ! COPY GRID_scalar:  (=)
  ! * subroutine copy_rscalar(E2,E1)
  ! * subroutine copy_cscalar(E2,E1)

  ! ZERO GRID_scalar  
  ! * subroutine zero_rscalar(E)
  ! * subroutine zero_cscalar(E)

  ! SCALAR MULTIPLICATION:
  ! * subroutine scMult_rscalar(c, E1, E2)
  ! * subroutine scMult_cscalar(c, E1, E2)

  ! SCALAR MULTIPLICATION: FUNCTION VERSION (c*E)
  ! * function scMult_rscalar_f(c, E1) result(E2)
  ! * function scMult_cscalar_f(c, E1) result(E2)

  ! scalar SUM:
  ! * subroutine add_rscalar(E1, E2, E3)
  ! * subroutine add_cscalar(E1, E2, E3)

  ! scalar SUM: FUNCTION VERSION (E3 = E1+E2)
  ! * function add_rscalar_f(E1, E2) result(E3)
  ! * function add_cscalar_f(E1, E2) result(E3)

  ! scalar SUBTRACT:
  ! * subroutine subtract_rscalar(E1, E2, E3)
  ! * subroutine subtract_cscalar(E1, E2, E3)

  ! scalar SUBTRACT: FUNCTION VERSION (E3 = E1-E2)
  ! * function subtract_rscalar_f(E1, E2) result(E3)
  ! * function subtract_cscalar_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF scalars: 
  ! * subroutine diagMult_rscalar(E1, E2, E3)
  ! * subroutine diagMult_cscalar(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF scalars: FUNCTION VERSION (c*E)
  ! * function diagMult_rscalar_f(E1, E2) result(E3)
  ! * function diagMult_cscalar_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF scalars and scalars:
  ! one combination is real and complex scalars and another is complex
  ! and real scalars as the sequence for inputs
  ! * subroutine diagMult_crscalar(E1, E2, E3)
  ! * subroutine diagMult_rcscalar(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF scalarS and SCALARS: FUNCTION VERSION (c*E)
  ! one combination is real and complex scalars and another is complex
  ! and real scalars as the sequence for inputs
  ! * function diagMult_crscalar_f(E1, E2) result(E3)
  ! * function diagMult_rcscalar_f(E1, E2) result(E3)

  ! scalar DOT PRODUCT
  ! * function dotProd_rscalar(E1, E2) result(r)
  ! * function dotProd_cscalar(E1, E2) result(c)

  ! COMPLEX LINEAR COMBINATION ... formally redundent with above
  ! only doing ones that we are sure to want
  ! * subroutine linCom_cscalar(inc1, E1, inc2, E2, E3)
  ! * subroutine scMultAdd_S_node(c, E1, E2)  (E2 = E2+c*E1)

  ! The algebraic routines expect all input and output
  ! variables to be of the correct type, already allocated,
  ! and of the correct size.  


  !****************************************************************************
  ! create_rscalar creates variable of derived type rscalar,
  ! using grid definition in structure "grid" ;
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage
  subroutine create_rscalar(igrid, E, gridType)

    implicit none
    type(grid3d_t), target, intent(in)    :: igrid
    ! the grid for which an scalar (center or corner) node field is being 
    ! initialized
    type (rscalar), intent(out)        :: E

    integer                            :: status,nx,ny,nz,nzEarth

    character (len=80)                 :: gridType

    if(E%allocated) then
       ! first deallocate memory for v 
       deallocate(E%v, STAT=status)
    end if

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    nzEarth = igrid%nzEarth
    E%nx = nx
    E%ny = ny
    E%nz = nz
    E%gridType = gridType

    ! allocate memory for v
    ! E%allocated will be true if all allocations succeed
    if (E%gridType == CENTER) then
       allocate(E%v(nx,ny,nz), STAT=status)
    else if (E%gridType == CORNER) then
       allocate(E%v(nx+1,ny+1,nz+1), STAT=status)
    else if(E%gridType == CELL_EARTH) then
       E%nz = nzEarth
       allocate(E%v(nx,ny,nzEarth), STAT=status)
    else
       write (0, *) 'not a known tag: create_rscalar'    
    end if
    E%allocated = status .EQ. 0

    if (E%allocated) then
       E%v = 0.0
    end if

  end subroutine create_rscalar  ! create_rscalar

  !****************************************************************************
  ! create_cscalar creates variable of derived type cscalar,
  ! using grid definition in structure "grid" ;
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage
  subroutine create_cscalar(igrid, E, gridType)

    implicit none
    type(grid3d_t), target, intent(in)     :: igrid
    ! the grid for which an scalar (center or corner) node field is being 
    ! initialized
    type (cscalar), intent(out)         :: E

    integer                             :: status,nx,ny,nz

    character (len=80)                  :: gridType

    if(E%allocated) then
       ! first deallocate memory for v 
       deallocate(E%v,STAT=status)
    end if

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    E%nx = nx
    E%ny = ny
    E%nz = nz
    E%gridType = gridType

    ! allocate memory for v ; 
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    if (E%gridType == CENTER) then
       allocate(E%v(nx,ny,nz), STAT=status)
    else if (E%gridType == CORNER) then
       allocate(E%v(nx+1,ny+1,nz+1), STAT=status)
    else 
       write (0, *) 'not a known tag: create_cscalar'    
    end if
    E%allocated = E%allocated .and. (status .EQ. 0)

    if (E%allocated) then
       E%v = C_ZERO
    end if
    ! print *, 'E%allocated', E%allocated

  end subroutine create_cscalar  ! create_cscalar


  !****************************************************************************
  ! deall_rscalar destoys variable of derived type rscalar,
  ! deallocating memory
  subroutine deall_rscalar(E)

    implicit none
    type (rscalar)  :: E
    integer	    :: status

    if(E%allocated) then
       ! deallocate memory for v 
       deallocate(E%v,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_rscalar  ! deall_rscalar


  !****************************************************************************
  ! deall_cscalar destoys variable of derived type cscalar,
  ! deallocating memory
  subroutine deall_cscalar(E)

    implicit none
    type (cscalar)  :: E
    integer	    :: status

    ! deallocate memory for v 
    if(E%allocated) then
       ! deallocate memory for v 
       deallocate(E%v,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_cscalar  ! deall_cscalar


  !****************************************************************************
  ! copy_rscalar makes an exact copy of derived data type 
  ! rscalar;   NOTE: first argument is output
  subroutine copy_rscalar(E2, E1)

    implicit none
    type (rscalar), intent(in)       :: E1
    type (rscalar), intent(inout)    :: E2
    integer	                     :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_rscalar'
    else

       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%v = E1%v
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage in copy_rscalar'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for v 
             deallocate(E2%v,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_rscalar(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          E2%v = E1%v
          E2%gridType = E1%gridType

       end if

    end if

  end subroutine copy_rscalar  ! copy_rscalar


  !****************************************************************************
  ! copy_cscalar makes an exact copy of derived data type 
  ! cscalar; 
  subroutine copy_cscalar(E2, E1) 

    implicit none
    type (cscalar), intent(in)            :: E1
    type (cscalar), intent(inout)         :: E2
    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_cscalar'
    else

       if((E2%nx == E1%nx).and.(E2%ny == E1%ny).and.(E2%nz == E1%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%v = E1%v
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage in copy_cscalar'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for v 
             deallocate(E2%v, STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_cscalar(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          E2%v = E1%v
          E2%gridType = E1%gridType

       end if

    end if

  end subroutine copy_cscalar  ! copy_cscalar


  !****************************************************************************
  ! zero_rscalar zeros variable of derived data type 
  ! rscalar;
  subroutine zero_rscalar(E)

    implicit none
    type (rscalar), intent(inout)   :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_rscalar: E not allocated'
    else

       E%v = R_ZERO

    end if

  end subroutine zero_rscalar


  !****************************************************************************
  ! zero_cscalar zeros variable of derived data type 
  ! cscalar;

  subroutine zero_cscalar(E)

    implicit none
    type (cscalar), intent(inout) :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_cscalar: E not allocated'
    else

       E%v = C_ZERO

    end if

  end subroutine zero_cscalar ! zero_cscalar


  !****************************************************************************
  ! scMult_cscalar multiplies scalar stored as devired data type
  ! cscalar with a complex scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_cscalar(c, E1, E2)

    implicit none
    complex(kind=selectedPrec), intent(in)                      :: c          
    ! a complex scalar to be multiplied with
    type (cscalar), intent(in)                       :: E1            
    type (cscalar), intent(inout)                    :: E2 

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for scMult_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component
             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage in scMult_cscalar'
          end if

       else
          write(0, *) 'Error:scMult_cscalar: scalars not same size'

       end if
    end if

  end subroutine scMult_cscalar ! scMult_cscalar


  !****************************************************************************
  ! scMult_cscalar_f multiplies scalar stored as devired data type
  ! cscalar with a complex scalar; function version
  function scMult_cscalar_f(c, E1) result(E2)

    implicit none
    complex(kind=selectedPrec), intent(in)                      :: c          
    ! a complex scalar to be multiplied with
    type (cscalar), intent(in)                       :: E1            
    type (cscalar)                                   :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_cscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component
             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage in scMult_cscalar_f'
          end if

       else

          write(0, *) 'Error:scMult_cscalar_f: scalars not same size'

       end if
    end if

  end function scMult_cscalar_f ! scMult_cscalar_f


  ! ***************************************************************************
  ! scMult_rscalar multiplies scalar stored as devired data type
  ! rscalar with a real scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_rscalar(c, E1, E2)

    implicit none
    real (kind=selectedPrec), intent(in)                         :: c          
    ! a real scalar to be multiplied with
    type (rscalar), intent(in)                       :: E1            
    type (rscalar), intent(inout)                    :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for v-component
             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage for scMult_rscalar'
          end if

       else

          write(0, *) 'Error:scMult_rscalar: scalars not same size'

       end if
    end if

  end subroutine scMult_rscalar ! scMult_rscalar


  !****************************************************************************
  ! scMult_rscalar_f multiplies scalar stored as devired data type
  ! rscalar with a real scalar; function version
  function scMult_rscalar_f(c, E1) result(E2)

    implicit none
    real (kind=selectedPrec), intent(in)                         :: c          
    ! a complex scalar to be multiplied with
    type (rscalar), intent(in)                       :: E1            
    type (rscalar)                                   :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for v-component
             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage for scMult_rscalar_f'
          end if

       else

          write(0, *) 'Error:scMult_rscalar_f: scalars not same size'

       end if
    end if

  end function scMult_rscalar_f ! scMult_rscalar_f


  !****************************************************************************
  ! add_rscalar adds scalars stored as devired data type
  ! rscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3     
    
   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rscalar'
       stop
    endif        

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_rscalar'
          end if

       else

          write(0, *) 'Error:add_rscalar: scalars not same size'

       end if
    end if

  end subroutine add_rscalar ! add_rscalar


  !****************************************************************************
  ! add_rscalar_f adds scalars stored as devired data type
  ! rscalar with ; function version
  function add_rscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar)                           :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_rscalar_f'
          end if

       else

          write(0, *) 'Error:add_rscalar_f: scalars not same size'

       end if
    end if

  end function add_rscalar_f ! add_rscalar_f


  !****************************************************************************
  ! add_cscalar adds scalars stored as devired data type
  ! cscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_cscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cscalar'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for add_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_cscalar'
          end if

       else

          write(0, *) 'Error:add_cscalar: scalars not same size'

       end if
    end if

  end subroutine add_cscalar ! add_cscalar


  !****************************************************************************
  ! add_cscalar_f adds scalars stored as devired data type
  ! cscalar with ; function version
  function add_cscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar)                           :: E3             

     if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_cscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_cscalar_f'
          end if

       else

          write(0, *) 'Error:add_cscalar_f: scalars not same size'

       end if
    end if

  end function add_cscalar_f ! add_cscalar_f


  !****************************************************************************
  ! subtract_rscalar subtracts scalars stored as devired data type
  ! rscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_rscalar'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_rscalar'
          end if

       else

          write(0, *) 'Error:add_rscalar: scalars not same size'

       end if
    end if

  end subroutine subtract_rscalar ! subtract_rscalar


  !****************************************************************************
  ! subtract_rscalar_f subtracts scalars stored as devired data type
  ! rscalar with ; function version
  function subtract_rscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar)                           :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_rscalar_f'
          end if

       else

          write(0, *) 'Error:subtract_rscalar_f: scalars not same size'

       end if
    end if

  end function subtract_rscalar_f ! subtract_rscalar_f


  !****************************************************************************
  ! subtract_cscalar subtracts scalars stored as devired data type
  ! cscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_cscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cscalar'
       stop
    endif  

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for subtract_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_cscalar'
          end if

       else

          write(0, *) 'Error:subtract_cscalar: scalars not same size'

       end if
    end if

  end subroutine subtract_cscalar ! subtract_cscalar


  !****************************************************************************
  ! subtract_cscalar_f subtracts scalars stored as devired data type
  ! cscalar with ; function version
  function subtract_cscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar)                           :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_cscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_cscalar_f'
          end if

       else

          write(0, *) 'Error:subtract_cscalar_f: scalars not same size'

       end if
    end if

  end function subtract_cscalar_f ! subtract_cscalar_f


  !****************************************************************************
  ! diagMult_rscalar multiplies two scalars E1, E2 stored as devired data 
  ! type rscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rscalar'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rscalar'
          end if

       else

          write(0, *) 'Error:diagMult_rscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_rscalar ! diagMult_rscalar


  !****************************************************************************
  ! diagMult_rscalar_f multiplies two scalars E1, E2 stored as devired 
  ! data type rscalar pointwise; function version
  function diagMult_rscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_rscalar_f: scalars not same size'

       end if
    end if

  end function diagMult_rscalar_f ! diagMult_rscalar_f


  !****************************************************************************
  ! diagMult_cscalar multiplies two scalars E1, E2 stored as devired data 
  ! type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_cscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar), intent(inout)            :: E3
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cscalar'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diaMult_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_cscalar'
          end if

       else

          write(0, *) 'Error:diagMult_cscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_cscalar ! diagMult_cscalar


  !****************************************************************************
  ! diagMult_cscalar_f multiplies two scalars E1, E2 stored as devired 
  ! data  type cscalar pointwise; function version
  function diagMult_cscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar)                           :: E3

   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_cscalar_f'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_cscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_cscalar_f: scalars not same size'

       end if
    end if

  end function diagMult_cscalar_f ! diagMult_cscalar_f


  !****************************************************************************
  ! diagMult_crscalar multiplies scalar E1 with scalar E2 stored as 
  ! devired type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_crscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1
    type (rscalar), intent(in)               :: E2    
    type (cscalar), intent(inout)            :: E3
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crscalar'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_crscalar'
          end if

       else

          write(0, *) 'Error:diagMult_crscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_crscalar ! diagMult_crscalar


  !****************************************************************************
  ! diagMult_crscalar_f multiplies scalar E1 with scalar E2 stored as 
  ! derived data type cscalar pointwise; function version
  function diagMult_crscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1
    type (rscalar), intent(in)               :: E2
    type (cscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crscalar_f'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_crscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_Node_MixedCR_f: scalars not same size'

       end if
    end if

  end function diagMult_crscalar_f ! diagMult_crscalar_f


  !****************************************************************************
  ! diagMult_rcscalar multiplies scalar E1 with scalar E2 stored as 
  ! derived type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rcscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1
    type (cscalar), intent(in)               :: E2    
    type (cscalar), intent(inout)            :: E3

   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcscalar'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rcscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rcscalar'
          end if

       else

          write(0, *) 'Error:diagMult_rcscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_rcscalar ! diagMult_rcscalar


  !****************************************************************************
  ! diagMult_rcscalar_f multiplies scalar E1 with scalar E2 stored as 
  ! derived data type cscalar pointwise; function version
  function diagMult_rcscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1
    type (cscalar), intent(in)               :: E2
    type (cscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcscalar_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if RHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rcscalar_f'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_Node_MixedRC_f'
          end if

       else

          write(0, *) 'Error:diagMult_Node_MixedRC_f: scalars not same size'

       end if
    end if

  end function diagMult_rcscalar_f ! diagMult_rcscalar_f


  !****************************************************************************
  ! dotProd_rscalar computes dot product of two vecors stored
  ! as derived data type rscalar, returning a real number
  function dotProd_rscalar_f(E1, E2) result(r)

    implicit none
    type (rscalar), intent(in)   :: E1, E2
    real (kind=selectedPrec)		     :: r 

    r = 0.0

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_rscalar'
       stop
    endif  

    ! Check whether both input scalars are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if (E1%gridType == E2%gridType) then

          r = r + sum(E1%v * E2%v)

       else
          write (0, *) 'not compatible usage for dotProd_rscalar'
       end if

    else

       write(0, *) 'Error:dotProd_rscalar: scalars not same size'

    end if

  end function dotProd_rscalar_f  ! dotProd_rscalar


  !****************************************************************************
  ! dotProd_cscalar computes dot product of two vecors stored
  ! as derived data type cscalar, returning a complex number
  function dotProd_cscalar_f(E1, E2) result(c)

    implicit none
    type (cscalar), intent(in)       :: E1, E2
    complex(kind=selectedPrec)		     :: c

    c = C_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_cscalar'
       stop
    endif  

    ! Check whether both input scalars are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if(E1%gridType == E2%gridType) then

          c = c + sum(conjg(E1%v) * E2%v)

       else
          write (0, *) 'not compatible usage for dotProd_cscalar'
       end if

    else

       write(0, *) 'Error:dotProd_cscalar: scalars not same size'

    end if

  end function dotProd_cscalar_f ! dotProd_cscalar


  !****************************************************************************
  ! linComb_cscalar computes linear combination of two scalars
  ! stored as derived data type cscalar; subroutine, not a function
  ! both input scalars must have the same dimension
  subroutine linComb_cscalar(inc1, E1, inc2, E2, E3)

    implicit none
    !   input scalars
    type (cscalar), intent(in)             :: E1, E2     
    !  input complex scalars
    complex (kind=8), intent(in)           :: inc1, inc2
    type (cscalar), intent(inout)          :: E3
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for linComb_cscalar'
       stop
    endif                

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for linComb_cscalar'
    else

       ! Check whether all scalars are of the same size
       if ((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! form linear combinatoin
             E3%v = inc1*E1%v + inc2*E2%v 

          else
             write (0, *) 'not compatible usage for linComb_cscalar'
          end if

       else

          write(0, *) 'Error:linComb_cscalar:  scalars not same size'

       end if
    end if

  end subroutine linComb_cscalar ! linComb_cscalar


  !****************************************************************************
  ! scMultadd_cscalar multiplies scalar E1 stored as derived data type
  ! cscalar with a complex scalar c, adding result to output scalar E2
  subroutine scMultAdd_cscalar(c, E1, E2)

    implicit none
    complex(kind=selectedPrec), intent(in)                      :: c          
    ! a complex scalar to be multiplied with
    type (cscalar), intent(in)                       :: E1            
    type (cscalar), intent(inout)                    :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultAdd_cscalar'
       stop
    endif  

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultAdd_cscalar'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component
             E2%v = E2%v + E1%v * c

          else
             write (0, *) 'not compatible usage for scMultAdd_cscalar'
          end if

       else

          write(0, *) 'Error:scMultAdd_cscalar: scalars not same size'

       end if
    end if

  end subroutine scMultAdd_cscalar ! scMultAdd_cscalar

  ! *************************************************************************
  ! * CornerVolume creates volume elements centered around the corners of
  ! * the grid, and stores them as real scalars with gridType=CORNER.
  
  subroutine CornerVolume(inGr, cV)

    implicit none
    type (grid3d_t), intent(in)          :: inGr  ! input grid
    type (rscalar), intent(inout)      :: cV    ! center volume as output
    integer                            :: ix, iy, iz        
    ! dummy variables

    ! Checks whether the size is the same
    if ((inGr%nx == cV%nx).and.&
         (inGr%ny == cV%ny).and.&
         (inGr%nz == cV%nz)) then

       if (cV%gridType == CORNER) then

          ! center volume is only using the internal corner nodes
          do ix = 2, inGr%nx
             do iy = 2, inGr%ny
                do iz = 2, inGr%nz

                   ! note that we are multiplying
                   ! using the distances with corner of a cell as a center
                   cV%v(ix, iy, iz) = inGr%delX(ix)*inGr%delY(iy)*&
                        inGr%delZ(iz)

                enddo
             enddo
          enddo

       else
          write (0, *) 'CornerVolume: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-grid size and center volume are not the same size'
    endif

  end subroutine CornerVolume

end module sg_scalar ! sg_scalar
