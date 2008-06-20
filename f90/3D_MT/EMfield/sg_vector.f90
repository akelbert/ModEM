! *****************************************************************************
module sg_vector
  ! This module creates data types for vector fields defined on
  ! edge/ face nodes of a staggered Cartesian grid, along with basic
  ! algebraic (vector space) operations. Example of these basic operations
  ! are allocation, deallocation, intialization, copying, and algebriac
  ! operations (linear combinations, scalar products, dot products).
  ! Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes) modules.

  use math_constants		! math/ physics constants
  use grid3d
  implicit none

  ! Generic interfaces are done through subroutines
  ! creates edge/ face nodes
  INTERFACE create
     module procedure create_rvector
     module procedure create_cvector
  END INTERFACE

  ! deallocates the edge/ face nodes
  INTERFACE deall
     module procedure deall_rvector
     module procedure deall_cvector
  END INTERFACE
  
  ! scalar value multiplies the edge/ face nodes
  INTERFACE scMult
     module procedure scMult_rvector
     module procedure scMult_cvector
     module procedure scMultReal_cvector
  END INTERFACE

  INTERFACE scMultAdd
     module procedure scMultAdd_cvector
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_cvector
  END INTERFACE

  ! adds the edge/ face nodes
  INTERFACE add
     module procedure add_rvector
     module procedure add_cvector
  END INTERFACE

  ! subtracts the edge/ face nodes
  INTERFACE subtract
     module procedure subtract_rvector
     module procedure subtract_cvector
  END INTERFACE

  ! pointwise vector (two vector data types) multiplication of edge/ face 
  ! nodes
  ! and pointwise real-complex (mixed) multiplication of edge/ face nodes
  ! Both are vector data types
  INTERFACE diagMult
     module procedure diagMult_rvector
     module procedure diagMult_cvector
     module procedure diagMult_rcvector
     module procedure diagMult_crvector
  END INTERFACE

  ! pointwise real-complex (mixed) division of edge/ face nodes
  ! Both are vector data types
  INTERFACE diagDiv
     module procedure diagDiv_rcvector
     module procedure diagDiv_crvector
  END INTERFACE

  ! zeros the edge/ face nodes
  INTERFACE zero
     module procedure zero_rvector
     module procedure zero_cvector
  END INTERFACE

  INTERFACE dotProd
     MODULE PROCEDURE dotProd_rvector_f
     MODULE PROCEDURE dotProd_cvector_f
  END INTERFACE

  INTERFACE dotProd_noConj
     MODULE PROCEDURE dotProd_noConj_cvector_f
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_cvector
     MODULE PROCEDURE copy_rvector
  END INTERFACE

  ! Interface operators are done through functions
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_rvector_f
     MODULE PROCEDURE add_cvector_f
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtract_rvector_f
     MODULE PROCEDURE subtract_cvector_f
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE scMult_rvector_f
     MODULE PROCEDURE scMult_cvector_f
     MODULE PROCEDURE scMultReal_cvector_f
     MODULE PROCEDURE diagMult_rvector_f
     MODULE PROCEDURE diagMult_cvector_f
     MODULE PROCEDURE diagMult_rcvector_f
     MODULE PROCEDURE diagMult_crvector_f 
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE diagDiv_rcvector_f
     MODULE PROCEDURE diagDiv_crvector_f 
  END INTERFACE

  public		::   create_rvector,  create_cvector, &
       deall_rvector, deall_cvector, &
       copy_rvector, copy_cvector, &
       zero_rvector, zero_cvector, &
       scMult_rvector, scMult_cvector, scMultReal_cvector, &
       scMult_rvector_f, scMult_cvector_f, scMultReal_cvector_f, &
       add_rvector, add_cvector, &
       add_rvector_f, add_cvector_f, &
       subtract_rvector, subtract_cvector, &
       subtract_rvector_f, subtract_cvector_f, &
       diagMult_rvector, diagMult_cvector, &
       diagMult_rvector_f, diagMult_cvector_f, &
       diagMult_rcvector, diagMult_crvector, &
       diagMult_rcvector_f, diagMult_crvector_f, &
       diagDiv_rcvector, diagDiv_crvector, &
       diagDiv_rcvector_f, diagDiv_crvector_f, &
       dotProd_rvector_f, dotProd_cvector_f, &
       getReal_cvector, getImag_cvector, &
       ! the two routines below are not part of any interface
  	linComb_cvector, scMultAdd_cvector


  ! ***************************************************************************
  ! type vector defines vector for either edge or face in a staggered grid as
  ! a complex field
  type :: cvector

     ! store the intention of the use in a character string defined
     ! in Grid3D as a parameter: EDGE or FACE are two possibilities
     character (len=80)	                             :: gridType

     ! Typical usage:  electrical fields on cell edges of
     ! staggered grid
     ! For example, in an EDGE, the dimensions would be
     ! x: edge nodes in x-direction: dimension Nx, Ny+1, Nz+1
     ! y: edge nodes in y-direction: dimension Nx+1, Ny, Nz+1
     ! z: edge nodes in z-direction: dimension Nx+1, Ny+1, Nz
     ! Note that the arrays are defined through dynamic memory allocation  
     complex(kind=selectedPrec), pointer, dimension(:,:,:)  :: x, y, z

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction: 
     integer                                          :: nx = 0, ny = 0, nz = 0

     ! allocated:  .true.  x, y, z arrays have been allocated
     logical		                              :: allocated = .false.

     ! pointer to parent grid
     type (grid3d_t), pointer                             :: grid		

  end type cvector
  
  ! ***************************************************************************
  ! type vector defines vector for either edge or face in a staggered grid as
  ! a real field
  type :: rvector

     ! store the intention of the use in a character string defined
     ! in Grid3D as a parameter: EDGE or FACE are two possibilities
     character (len=80)	                              :: gridType

     ! Typical usage:  conductivity averaged on cell edges of
     ! staggered grid
     ! x: edge nodes in x-direction: dimension Nx, Ny+1, Nz+1
     ! y: edge nodes in y-direction: dimension Nx+1, Ny, Nz+1
     ! z: edge nodes in z-direction: dimension Nx+1, Ny+1, Nz
     ! Note that the arrays are defined through dynamic memory allocation  
     real(kind=selectedPrec), pointer, dimension(:,:,:) :: x,y,z

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction: 
     integer                                          :: nx = 0, ny = 0, nz = 0


     ! allocated:  .true.  x, y, z arrays have been allocated
     logical		                              :: allocated = .false.

     ! pointer to parent grid
     type (grid3d_t), pointer                             :: grid

  end type rvector
  
Contains
  ! CREATE GRID_edge/ face VECTORS
  ! * subroutine create_rvector(igrid, E, gridType)
  ! * subroutine create_cvector(igrid, E, gridType)

  ! DEALLOCATE GRID_edge/ face VECTORS
  ! * subroutine deall_rvector(E)
  ! * subroutine deall_cvector(E)

  ! COPY GRID_edge/ face VECTORS :  (=)
  ! * subroutine copy_rvector(E2,E1)
  ! * subroutine copy_cvector(E2,E1)

  ! ZERO GRID_edge/ face VECTORS 
  ! * subroutine zero_rvector(E)
  ! * subroutine zero_cvector(E)

  ! SCALAR MULTIPLICATION:
  ! * subroutine scMult_rvector(c, E1, E2)
  ! * subroutine scMult_cvector(c, E1, E2)
  ! * subroutine scMultReal_cvector(c, E1, E2)

  ! SCALAR MULTIPLICATION: FUNCTION VERSION (c*E)
  ! * function scMult_rvector_f(c, E1) result(E2)
  ! * function scMult_cvector_f(c, E1) result(E2)
  ! * function scMultReal_cvector_f(c, E1) result(E2)

  ! VECTOR SUM:
  ! * subroutine add_rvector(E1, E2, E3)
  ! * subroutine add_cvector(E1, E2, E3)

  ! VECTOR SUM: FUNCTION VERSION (E3 = E1+E2)
  ! * function add_rvector_f(E1, E2) result(E3)
  ! * function add_cvector_f(E1, E2) result(E3)

  ! VECTOR SUBTRACT:
  ! * subroutine subtract_rvector(E1, E2, E3)
  ! * subroutine subtract_cvector(E1, E2, E3)

  ! VECTOR SUBTRACT: FUNCTION VERSION (E3 = E1-E2)
  ! * function subtract_rvector_f(E1, E2) result(E3)
  ! * function subtract_cvector_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF VECTORS: 
  ! * subroutine diagMult_rvector(E1, E2, E3)
  ! * subroutine diagMult_cvector(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF VECTORS: FUNCTION VERSION (c*E)
  ! * function diagMult_rvector_f(E1, E2) result(E3)
  ! * function diagMult_cvector_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF VECTORS and SCALARS:
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * subroutine diagMult_crvector(E1, E2, E3)
  ! * subroutine diagMult_rcvector(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF VECTORS and SCALARS: FUNCTION VERSION (c*E)
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * function diagMult_crvector_f(E1, E2) result(E3)
  ! * function diagMult_rcvector_f(E1, E2) result(E3)

  ! POINTWISE DIVISION OF VECTORS and SCALARS:
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * subroutine diagDiv_crvector(E1, E2, E3)
  ! * subroutine diagDiv_rcvector(E1, E2, E3)

  ! POINTWISE DIVISION OF VECTORS and SCALARS: FUNCTION VERSION (c*E)
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * function diagDiv_crvector_f(E1, E2) result(E3)
  ! * function diagDiv_rcvector_f(E1, E2) result(E3)

  ! VECTOR DOT PRODUCT
  ! * function dotProd_rvector(E1, E2) result(r)
  ! * function dotProd_cvector(E1, E2) result(c)

  ! COMPLEX LINEAR COMBINATION ... formally redundent with above
  ! only doing ones that we are sure to want
  ! * subroutine linCom_cvector(inc1, E1, inc2, E2, E3)
  ! * subroutine scMultAdd_V_node(c, E1, E2)  (E2 = E2+c*E1)

  ! The algebraic routines expect all input and output
  ! variables to be of the correct type, already allocated,
  ! and of the correct size.  


  !****************************************************************************
  ! create_rvector creates variable of derived type rvector,
  ! using grid definition in structure "grid" ;
  ! allocates memory in x,y,z component arrays
  ! gridType is a character string to describe intended usage
  subroutine create_rvector(igrid, E, gridType)

    implicit none
    type(grid3d_t), target, intent(in)    :: igrid
    ! the grid for which an edge/ face node field is being initialized
    type (rvector), intent(out)        :: E

    integer                            :: status,nx,ny,nz

    character (len=80), intent(in)     :: gridType

    if(E%allocated) then
       ! first deallocate memory for x,y,z 
       deallocate(E%x, E%y, E%z, STAT=status)
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

    ! gridType
    E%gridType = gridType

    ! allocate memory for x,y,z 
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    if (E%gridType == EDGE) then
       allocate(E%x(nx,ny+1,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%y(nx+1,ny,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%z(nx+1,ny+1,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
    else if (E%gridType == FACE) then
       allocate(E%x(nx+1,ny,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%y(nx,ny+1,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%z(nx,ny,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
    else 
       write (0, *) 'not a known tag'
    end if

    if (E%allocated) then
       E%x = 0.0
       E%y = 0.0
       E%z = 0.0
    end if

  end subroutine create_rvector  ! create_rvector


  !****************************************************************************
  ! create_cvector creates variable of derived type cvector,
  ! using grid definition in structure "grid" ;
  ! allocates memory in x,y,z component arrays
  ! gridType is a character string to describe intended usage
  subroutine create_cvector(igrid, E, gridType)

    implicit none
    type(grid3d_t), target, intent(in)     :: igrid
    ! the grid for which an edge/ face node field is being initialized
    type (cvector), intent(out)         :: E

    integer                             :: status,nx,ny,nz

    character (len=80), intent(in)      :: gridType

    if(E%allocated) then
       ! first deallocate memory for x,y,z 
       deallocate(E%x, E%y, E%z,STAT=status)
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
    ! print *, 'nx, ny, nz', nx, ny, nz

    ! gridType
    E%gridType = gridType

    ! allocate memory for x,y,z ; 
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    if (E%gridType == EDGE) then
       allocate(E%x(nx,ny+1,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%y(nx+1,ny,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%z(nx+1,ny+1,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
    else if (E%gridType == FACE) then
       allocate(E%x(nx+1,ny,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%y(nx,ny+1,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%z(nx,ny,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
    else 
       write (0, *) 'not a known tag'
    end if

    if (E%allocated) then
       E%x = C_ZERO
       E%y = C_ZERO
       E%z = C_ZERO
    end if
    ! print *, 'E%allocated', E%allocated

  end subroutine create_cvector  ! create_cvector


  !****************************************************************************
  ! deall_rvector destoys variable of derived type rvector,
  ! deallocating memory
  subroutine deall_rvector(E)

    implicit none
    type (rvector)  :: E
    integer	    :: status

    if(E%allocated) then
       ! deallocate memory for x,y,z 
       deallocate(E%x, E%y, E%z,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_rvector  ! deall_rvector


  !****************************************************************************
  ! deall_cvector destoys variable of derived type cvector,
  ! deallocating memory
  subroutine deall_cvector(E)

    implicit none
    type (cvector)  :: E
    integer	    :: status

    ! deallocate memory for x,y,z 
    if(E%allocated) then
       ! deallocate memory for x,y,z 
       deallocate(E%x, E%y, E%z,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_cvector  ! deall_cvector


  !****************************************************************************
  ! copy_rvector makes an exact copy of derived data type 
  ! rvector;   NOTE: first argument is output
  subroutine copy_rvector(E2, E1)

    implicit none
    type (rvector), intent(in)       :: E1
    type (rvector), intent(inout)    :: E2
    integer	                     :: status
 
    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_rvector'
    else

       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%x = E1%x
             E2%y = E1%y
             E2%z = E1%z
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage for copy_rvector'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for x,y,z 
             deallocate(E2%x, E2%y, E2%z,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_rvector(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          E2%x = E1%x
          E2%y = E1%y
          E2%z = E1%z
          E2%gridType = E1%gridType

       end if

    end if

  end subroutine copy_rvector  ! copy_rvector


  !****************************************************************************
  ! copy_cvector makes an exact copy of derived data type 
  ! cvector; 
  subroutine copy_cvector(E2, E1) 

    implicit none
    type (cvector), intent(in)            :: E1
    type (cvector), intent(inout)         :: E2
    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_cvector'
    else

       if((E2%nx == E1%nx).and.(E2%ny == E1%ny).and.(E2%nz == E1%nz)) then

          if  (E1%gridType == E2%gridType) then

             ! just copy components
             E2%x = E1%x
             E2%y = E1%y
             E2%z = E1%z
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage for copy_cvector'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for x,y,z 
             deallocate(E2%x, E2%y, E2%z,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_cvector(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          E2%x = E1%x
          E2%y = E1%y
          E2%z = E1%z
          E2%gridType = E1%gridType

       end if

    end if

  end subroutine copy_cvector  ! copy_cvector


  !****************************************************************************
  ! zero_rvector zeros variable of derived data type 
  ! rvector;
  subroutine zero_rvector(E)

    implicit none
    type (rvector), intent(inout)   :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_rvector: E not allocated'
    else

       E%x = R_ZERO
       E%y = R_ZERO
       E%z = R_ZERO

    end if

  end subroutine zero_rvector


  !****************************************************************************
  ! zero_cvector zeros variable of derived data type 
  ! cvector;
  subroutine zero_cvector(E)

    implicit none
    type (cvector), intent(inout) :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_cvector: E not allocated'
    else

       E%x = C_ZERO
       E%y = C_ZERO
       E%z = C_ZERO

    end if

  end subroutine zero_cvector ! zero_cvector


  !****************************************************************************
  ! scMult_cvector multiplies vector stored as derived data type
  ! cvector with a complex scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_cvector(c, E1, E2)

    implicit none
    complex(kind=selectedPrec), intent(in)                      :: c          
    ! a complex scalar to be multiplied with
    type (cvector), intent(in)                       :: E1            
    type (cvector), intent(inout)                    :: E2 

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cvector'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for scMult_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_cvector'
          end if

       else
          write(0, *) 'Error:scMult_cvector: vectors not same size'

       end if
    end if

  end subroutine scMult_cvector ! scMult_cvector


  !****************************************************************************
  ! scMult_cvector_f multiplies vector stored as derived data type
  ! cvector with a complex scalar; function version
  function scMult_cvector_f(c, E1) result(E2)

    implicit none
    complex(kind=selectedPrec), intent(in)                      :: c          
    ! a complex scalar to be multiplied with
    type (cvector), intent(in)                       :: E1            
    type (cvector)                                   :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cvector_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_cvector_f'
          end if

       else

          write(0, *) 'Error:scMult_cvector_f: vectors not same size'

       end if
    end if

  end function scMult_cvector_f ! scMult_cvector_f


  !****************************************************************************
  ! scMultReal_cvector multiplies vector stored as derived data type
  ! cvector with a real scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMultReal_cvector(c, E1, E2)

    implicit none
    real (kind=selectedPrec), intent(in)                         :: c          
    ! a real scalar to be multiplied with
    type (cvector), intent(in)                       :: E1            
    type (cvector), intent(inout)                    :: E2 

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultReal_cvector'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for scMultReal_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMultReal_cvector'
          end if

       else
          write(0, *) 'Error:scMultReal_cvector: vectors not same size'

       end if
    end if

  end subroutine scMultReal_cvector ! scMultReal_cvector


  !****************************************************************************
  ! scMult_cvector_f multiplies vector stored as derived data type
  ! cvector with a real scalar; function version
  function scMultReal_cvector_f(c, E1) result(E2)

    implicit none
    real (kind=selectedPrec), intent(in)			     :: c          
    ! a real scalar to be multiplied with
    type (cvector), intent(in)                       :: E1            
    type (cvector)                                   :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultReal_cvector_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultReal_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMultReal_cvector_f'
          end if

       else

          write(0, *) 'Error:scMultReal_cvector_f: vectors not same size'

       end if
    end if

  end function scMultReal_cvector_f ! scMultReal_cvector_f


  ! ***************************************************************************
  ! scMult_rvector multiplies vector stored as derived data type
  ! rvector with a real scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_rvector(c, E1, E2)

    implicit none
    real (kind=selectedPrec), intent(in)                         :: c          
    ! a real scalar to be multiplied with
    type (rvector), intent(in)                       :: E1            
    type (rvector), intent(inout)                    :: E2             

   if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rvector'
       stop
    endif
    
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_rvector'
          end if

       else

          write(0, *) 'Error:scMult_rvector: vectors not same size'

       end if
    end if

  end subroutine scMult_rvector ! scMult_rvector


  !****************************************************************************
  ! scMult_rvector_f multiplies vector stored as derived data type
  ! rvector with a real scalar; function version
  function scMult_rvector_f(c, E1) result(E2)

    implicit none
    real (kind=selectedPrec), intent(in)                         :: c          
    ! a complex scalar to be multiplied with
    type (rvector), intent(in)                       :: E1            
    type (rvector)                                   :: E2 
    
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rvector_f'
       stop
    endif            

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_rvector_f'
          end if

       else

          write(0, *) 'Error:scMult_rvector_f: vectors not same size'

       end if
    end if

  end function scMult_rvector_f ! scMult_rvector_f


  !****************************************************************************
  ! add_rvector adds vectors stored as derived data type
  ! rvector with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rvector'
       stop
    endif        

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_rvector'
          end if

       else

          write(0, *) 'Error:add_rvector: vectors not same size'

       end if
    end if

  end subroutine add_rvector ! add_rvector


  !****************************************************************************
  ! add_rvector_f adds vectors stored as derived data type
  ! rvector with ; function version
  function add_rvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector)                           :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rvector_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_rvector_f'
          end if

       else

          write(0, *) 'Error:add_rvector_f: vectors not same size'

       end if
    end if

  end function add_rvector_f ! add_rvector_f


  !****************************************************************************
  ! add_cvector adds vectors stored as derived data type
  ! cvector with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_cvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for add_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_cvector'
          end if

       else

          write(0, *) 'Error:add_cvector: vectors not same size'

       end if
    end if

  end subroutine add_cvector ! add_cvector


  !****************************************************************************
  ! add_cvector_f adds vectors stored as derived data type
  ! cvector with ; function version
  function add_cvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector)                           :: E3             

     if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cvector_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_cvector_f'
          end if

       else

          write(0, *) 'Error:add_cvector: vectors not same size'

       end if
    end if

  end function add_cvector_f ! add_cvector_f


  !****************************************************************************
  ! subtract_rvector subtracts vectors stored as derived data type rvector with
  ! ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_rvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_rvector'
          end if

       else

          write(0, *) 'Error: subtract_rvector: vectors not same size'

       end if
    end if

  end subroutine subtract_rvector ! subtract_rvector


  !****************************************************************************
  ! subtract_rvector_f subtracts vectors stored as derived data type
  ! rvector with ; function version
  function subtract_rvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector)                           :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rvector_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_rvector_f'
          end if

       else

          write(0, *) 'Error:subtract_rvector_f: vectors not same size'

       end if
    end if

  end function subtract_rvector_f ! subtract_rvector_f


  !****************************************************************************
  ! subtract_cvector subtracts vectors stored as derived data type
  ! cvector with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_cvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector), intent(inout)            :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cvector'
       stop
    endif  

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for subtract_cvector'
    else
    
       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_cvector'
          end if

       else

          write(0, *) 'Error:subtract_cvector: vectors not same size'

       end if
    end if

  end subroutine subtract_cvector ! subtract_cvector


  !****************************************************************************
  ! subtract_cvector_f subtracts vectors stored as derived data type
  ! cvector with ; function version
  function subtract_cvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector)                           :: E3             

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cvector_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_cvector_f'
          end if

       else

          write(0, *) 'Error: subtract_cvector_f: vectors not same size'

       end if
    end if

  end function subtract_cvector_f ! subtract_cvector_f


  !****************************************************************************
  ! diagMult_rvector multiplies two vectors E1, E2 stored as derived data 
  ! type rvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector), intent(inout)            :: E3
    
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rvector'
          end if

       else

          write(0, *) 'Error:diagMult_rvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_rvector ! diagMult_rvector


  !****************************************************************************
  ! diagMult_rvector_f multiplies two vectors E1, E2 stored as derived 
  ! data type rvector pointwise; function version
  function diagMult_rvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rvector_f'
       stop
    endif 

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_rvector_f: vectors not same size'

       end if
    end if

  end function diagMult_rvector_f ! diagMult_rvector_f


  !****************************************************************************
  ! diagMult_cvector multiplies two vectors E1, E2 stored as derived data 
  ! type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_cvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diaMult_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_cvector'
          end if

       else

          write(0, *) 'Error:diagMult_cvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_cvector ! diagMult_cvector


  !****************************************************************************
  ! diagMult_cvector_f multiplies two vectors E1, E2 stored as derived 
  ! data  type cvector pointwise; function version
  function diagMult_cvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector)                           :: E3

   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cvector_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if RHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'RHS was not allocated for diagMult_cvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_cvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_cvector_f: vectors not same size'

       end if
    end if

  end function diagMult_cvector_f ! diagMult_cvector_f


  !****************************************************************************
  ! diagMult_crvector multiplies complex vector E1 with scalar vector E2 
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_crvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2    
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crvector'
    else
    
       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_crvector'
          end if

       else

          write(0, *) 'Error:diagMult_crvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_crvector ! diagMult_crvector


  !****************************************************************************
  ! diagMult_crvector_f multiplies complex vector E1 with real vector 
  ! E2 stored as derived data type cvector pointwise; function version
  function diagMult_crvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crvector_f'
       stop
    endif  

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_crvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_Node_MixedCR_f: vectors not same size'

       end if
    end if

  end function diagMult_crvector_f ! diagMult_crvector_f


  !****************************************************************************
  ! diagMult_rcvector multiplies real vector E1 with complex vector E2 
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rcvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2    
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rcvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rcvector'
          end if

       else

          write(0, *) 'Error:diagMult_rcvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_rcvector ! diagMult_rcvector


  !****************************************************************************
  ! diagMult_rcvector_f multiplies real vector E1 with complex vector E2 
  ! stored as derived data type cvector pointwise; function version
  function diagMult_rcvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcvector_f'
       stop
    endif 

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if RHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'RHS was not allocated for diagMult_rcvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rcvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_Node_MixedRC_f: vectors not same size'

       end if
    end if

  end function diagMult_rcvector_f ! diagMult_rcvector_f


  !****************************************************************************
  ! diagDiv_crvector divides complex vector E1 with scalar vector E2 
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_crvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2    
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_crvector'
       stop
    endif 
    
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_crvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for x,y,z-components
             E3%x = E1%x / E2%x
             E3%y = E1%y / E2%y
             E3%z = E1%z / E2%z

          else
             write (0, *) 'not compatible usage for diagDiv_crvector'
          end if

       else

          write(0, *) 'Error:diagDiv_crvector: vectors not same size'

       end if
    end if

  end subroutine diagDiv_crvector ! diagDiv_crvector


  !****************************************************************************
  ! diagDiv_crvector_f divides complex vector E1 with real vector 
  ! E2 stored as derived data type cvector pointwise; function version
  function diagDiv_crvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_crvector_f'
       stop
    endif 

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_crvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for x,y,z-components
             E3%x = E1%x / E2%x
             E3%y = E1%y / E2%y
             E3%z = E1%z / E2%z

          else
             write (0, *) 'not compatible usage for diagDicrvector_f'
          end if

       else

          write(0, *) 'Error:diagDicrvector_f: vectors not same size'

       end if
    end if

  end function diagDiv_crvector_f ! diagDiv_crvector_f


  !****************************************************************************
  ! diagDiv_rcvector divides real vector E1 with complex vector E2 
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_rcvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2    
    type (cvector), intent(inout)            :: E3

       if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rcvector'
       stop
    endif  

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rcvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise divides for x,y,z-components
             E3%x = E2%x / E1%x
             E3%y = E2%y / E1%y
             E3%z = E2%z / E1%z

          else
             write (0, *) 'not compatible usage for diagDiv_rcvector'
          end if

       else

          write(0, *) 'Error:diagDiv_rcvector: vectors not same size'

       end if
    end if

  end subroutine diagDiv_rcvector ! diagDiv_rcvector


  !****************************************************************************
  ! diagDiv_rcvector_f divides real vector E1 with complex vector E2 
  ! stored as derived data type cvector pointwise; function version
  function diagDiv_rcvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rcvector_f'
       stop
    endif 

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rcvector_f'
    else
    
       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise divides for x,y,z-components
             E3%x = E2%x / E1%x
             E3%y = E2%y / E1%y
             E3%z = E2%z / E1%z

          else
             write (0, *) 'not compatible usage for diagDiv_rcvector_f'
          end if

       else

          write(0, *) 'Error:diagDircvector_f: vectors not same size'

       end if
    end if

  end function diagDiv_rcvector_f ! diagDiv_rcvector_f


  !****************************************************************************
  ! dotProd_rvector computes dot product of two vecors stored
  ! as derived data type rvector, returning a real number
  function dotProd_rvector_f(E1, E2) result(r)

    implicit none
    type (rvector), intent(in)   :: E1, E2
    real (kind=selectedPrec)		     :: r 

    r = 0.0
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_rvector'
       stop
    endif  

    ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if ((E1%gridType == E2%gridType)) then

          r = r + sum(E1%x * E2%x)
          r = r + sum(E1%y * E2%y)
          r = r + sum(E1%z * E2%z)

       else
          write (0, *) 'dotProd_rvector: not compatible usage'
       end if

    else

       write(0, *) 'Error:dotProd_rvector: vectors not same size'

    end if

  end function dotProd_rvector_f  ! dotProd_rvector

  !****************************************************************************
  ! dotProd_cvector computes dot product of two vecors stored
  ! as derived data type cvector, returning a complex number
  function dotProd_cvector_f(E1, E2) result(c)

    implicit none
    type (cvector), intent(in)       :: E1, E2
    complex(kind=selectedPrec)		     :: c

    c = C_ZERO
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_cvector'
       stop
    endif  

    ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if ((E1%gridType == E2%gridType)) then

          c = c + sum(conjg(E1%x) * E2%x)
          c = c + sum(conjg(E1%y) * E2%y)
          c = c + sum(conjg(E1%z) * E2%z)

       else
          write (0, *) 'dotProd_cvector: not compatible usage'
       end if

    else

       write(0, *) 'Error:dotProd_cvector: vectors not same size'

    end if

  end function dotProd_cvector_f ! dotProd_cvector

  !****************************************************************************
  ! dotProd_noConj_cvector computes dot product of two vecors stored
  ! as derived data type cvector, returning a complex number
  !  IN THIS VERSION CONJUGATES ARE NOT USED

  function dotProd_noConj_cvector_f(E1, E2) result(c)

    implicit none
    type (cvector), intent(in)       :: E1, E2
    complex(kind=selectedPrec)		     :: c

    c = C_ZERO
    
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_noConj_cvector'
       stop
    endif  

    ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
       if ((E1%gridType == E2%gridType)) then
          c = c + sum(E1%x * E2%x)
          c = c + sum(E1%y * E2%y)
          c = c + sum(E1%z * E2%z)
       else
          write (0, *) 'dotProd_noConj_cvector: not compatible usage'
       end if
    else
       write(0, *) 'Error:dotProd_noConj_cvector: vectors not same size'
    end if

  end function dotProd_noConj_cvector_f ! dotProd_cvector


  !****************************************************************************
  ! linComb_cvector computes linear combination of two vectors
  ! stored as derived data type cvector; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_cvector(inc1, E1, inc2, E2, E3)

    implicit none
    !   input vectors
    type (cvector), intent(in)             :: E1, E2     
    !  input complex scalars
    complex (kind=8), intent(in)           :: inc1, inc2
    type (cvector), intent(inout)          :: E3              

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for linComb_cvector'
       stop
    endif                

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for linComb_cvector'
    else

       ! Check whether all vectors are of the same size
       if ((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! form linear combinatoin
             E3%x = inc1*E1%x + inc2*E2%x  
             E3%y = inc1*E1%y + inc2*E2%y
             E3%z = inc1*E1%z + inc2*E2%z

          else
             write (0, *) 'not compatible usage for linComb_cvector'
          end if

       else

          write(0, *) 'Error:linComb_cvector:  vectors not same size'

       end if
    end if

  end subroutine linComb_cvector ! linComb_cvector


  !****************************************************************************
  ! scMultadd_cvector multiplies vector E1 stored as derived data type
  ! cvector with a complex scalar c, adding result to output vector E2
  subroutine scMultAdd_cvector(c, E1, E2)

    implicit none
    complex(kind=selectedPrec), intent(in)                      :: c          
    ! a complex scalar to be multiplied with
    type (cvector), intent(in)                       :: E1            
    type (cvector), intent(inout)                    :: E2             

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultAdd_cvector'
       stop
    endif  

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultAdd_cvector'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if ((E1%gridType == E2%gridType)) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E2%x + E1%x * c
             E2%y = E2%y + E1%y * c
             E2%z = E2%z + E1%z * c

          else
             write (0, *) 'not compatible usage for scMultAdd_cvector'
          end if

       else

          write(0, *) 'Error:scMultAdd_cvector: vectors not same size'

       end if
    end if

  end subroutine scMultAdd_cvector ! scMultAdd_cvector

!***************************************************************

  subroutine getReal_cvector(eC,eR)

  !   takes real part of cvector object eC
  !   for now: let's assume that eR and eC are consistent,
  !   and that eR has been allocated ... clean up the present
  !   garbage implmentation of error checking later!

  type(cvector), intent(in)		:: eC
  type(rvector), intent(inout)		:: eR

  eR%x = real(eC%x)
  eR%y = real(eC%y)
  eR%z = real(eC%z)

  end subroutine getReal_cvector

!***************************************************************

  subroutine getImag_cvector(eC,eR)

  !   takes real part of cvector object eC
  !   for now: let's assume that eR and eC are consistent,
  !   and that eR has been allocated ... clean up the present
  !   garbage implmentation of error checking later!

  type(cvector), intent(in)		:: eC
  type(rvector), intent(inout)		:: eR

  eR%x = imag(eC%x)
  eR%y = imag(eC%y)
  eR%z = imag(eC%z)

  end subroutine getImag_cvector
  
!***************************************************************

  ! *************************************************************************
  ! * EdgeVolume creates volume elements centered around the edges of
  ! * the grid, and stores them as real vectors with gridType=EDGE.
  ! *
  ! * A diagonal matrix multiplication of the edge volume with the difference
  ! * equations enables us to make a symmetrical matrix. Remember,
  ! * the electrical fields are defined on the center of the edges, therefore,
  ! * the edge volume is centered about the electrical field measurement.

  subroutine EdgeVolume(inGr, eV)

    implicit none
    type (grid3d_t), intent(in)             :: inGr     ! input model
    type (rvector), intent(inout)         :: eV       ! edge volume
    integer                               :: ix, iy, iz        
    ! dummy variables

    ! Checks whether the size is the same
    if ((inGr%nx == eV%nx).and.&
         (inGr%ny == eV%ny).and.&
         (inGr%nz == eV%nz)) then

       if (eV%gridType == EDGE) then

          ! edge volume are made for all the edges
          ! for x-components
          do ix = 1,inGr%nx
             do iy = 1,inGr%ny+1
                do iz = 1,inGr%nz+1

                   ! eV%x values are centered within dx.
                   eV%x(ix, iy, iz) = inGr%dx(ix)*inGr%delY(iy)*&
                        inGr%delZ(iz)

                enddo
             enddo
          enddo

          ! edge volume are made for all the edges
          ! for y-components
          do ix = 1,inGr%nx+1
             do iy = 1,inGr%ny
                do iz = 1,inGr%nz+1

                   ! eV%y values are centered within dy.
                   eV%y(ix, iy, iz) = inGr%delX(ix)*inGr%dy(iy)*&
                        inGr%delZ(iz)

                enddo
             enddo
          enddo

          ! edge volume are made for all the edges
          ! for z-components
          do ix = 1,inGr%nx+1
             do iy = 1,inGr%ny+1
                do iz = 1,inGr%nz

                   ! eV%z values are centered within dz. 
                   eV%z(ix, iy, iz) = inGr%delX(ix)*inGr%delY(iy)*&
                        inGr%dz(iz)

                enddo
             enddo
          enddo

       else
          write (0, *) 'EdgeVolume: not compatible usage for existing data types'
       end if


    else
       write(0, *) 'Error-grid size and edge volume are not the same size'
    endif

  end subroutine EdgeVolume  ! EdgeVolume

end module sg_vector ! sg_vector
