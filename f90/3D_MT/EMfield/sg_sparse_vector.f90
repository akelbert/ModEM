! *****************************************************************************
module sg_sparse_vector
  ! Sparse complex vector operations. At present, sparse vectors are used for
  ! representation of data functionals. Only complex vectors are supported, and
  ! only those routines needed for modeling and initial work on inversion have
  ! been developed.
  ! Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes)
  !  modules.

  use math_constants
  use grid3d
  use sg_vector
  implicit none

  INTERFACE create
     module procedure create_sparsevecc
  END INTERFACE

  INTERFACE deall
     module procedure deall_sparsevecc
  END INTERFACE
  
  INTERFACE scMult
     module procedure scMult_sparsevecc
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_sparsevecc
  END INTERFACE

  INTERFACE add
     module procedure add_scvector
  END INTERFACE

  INTERFACE dotProd
     MODULE PROCEDURE dotProd_scvector_f
     !MODULE PROCEDURE dotProd_csvector_f
  END INTERFACE

  INTERFACE dotProd_noConj
     MODULE PROCEDURE dotProd_noConj_scvector_f
     !MODULE PROCEDURE dotProd_noConj_csvector_f
  END INTERFACE
	
  interface assignment (=)
     MODULE PROCEDURE copy_sparsevecc
     MODULE PROCEDURE copy_csvector
  end interface

  public	:: sparsevecc
  public	:: create_sparsevecc, deall_sparsevecc 
  public  :: copy_csvector
  public	:: copy_sparsevecc, linComb_sparsevecc, scMult_sparsevecc
  public	:: dotProd_scvector_f, dotProd_csvector_f
  public	:: add_scvector


!**************************************************************************
  type :: sparsevecc

     ! complex vector defined on edge/ face nodes: 
     ! intention of use as a character string: 'Edge', or 'Face'
     character (len=80)	                             	:: gridType=''
     ! nCoeff is number of non-zero nodes
     integer 						:: nCoeff  = 0
     ! xyz = 1,2,3 refers to x, y or z components, 
     ! i,j,k are arrays of indices that defines grid location
     integer , pointer, dimension(:) 		:: i,j,k,xyz
     ! c is complex array of coefficients
     complex (kind=selectedPrec), pointer, dimension(:) 	:: c
     ! has sparse vector been allocated?
     logical					:: allocated = .false.
     ! pointer to the parent grid
     type (grid3d_t), pointer                 	:: grid

  end type sparsevecc


Contains

  ! **********************************************************************
  ! delete/ deallocate the sparse vector
  subroutine deall_sparsevecc(oldLC)

    implicit none
    type (sparsevecc), intent(inout)                :: oldLC
    integer                                        :: status

       deallocate(oldLC%i,STAT=status)
       deallocate(oldLC%j, STAT=status)
       deallocate(oldLC%k, STAT=status)
       deallocate(oldLC%xyz, STAT=status)
       deallocate(oldLC%c, STAT=status)
       oldLC%gridType = ''
       oldLC%allocated = .false.

  end subroutine deall_sparsevecc


  ! **********************************************************************
  ! create an object of type sparsevecc of length nCoeff
  subroutine create_sparsevecc(nCoeff,newLC, gridType)

    implicit none	
    integer, intent(in) 			:: nCoeff
    type (sparsevecc), intent(inout) 		:: newLC
    character (len=80), intent(in)     		:: gridType
    integer					:: status
    
    ! the old baggage is out of the door
    if(newLC%allocated) then
       deallocate(newLC%i, STAT = status)
       deallocate(newLC%j, STAT = status)
       deallocate(newLC%k, STAT = status)
       deallocate(newLC%xyz, STAT = status)
       deallocate(newLC%c, STAT = status)
       newLC%gridType = ''
       newLC%allocated = .false.
    endif

    newLC%allocated = .true.
    allocate(newLC%i(nCoeff),STAT=status) 
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%j(nCoeff),STAT=status) 
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%k(nCoeff),STAT=status) 
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%xyz(nCoeff),STAT=status) 
    newLC%allocated = newLC%allocated .and. (status .eq. 0)
    allocate(newLC%c(nCoeff),STAT=status) 
    newLC%allocated = newLC%allocated .and. (status .eq. 0)

    newLC%nCoeff = nCoeff
    newLC%i = 0
    newLC%j = 0
    newLC%k = 0
    newLC%xyz = 0
    newLC%c = C_ZERO
    newLC%gridType = gridType

  end subroutine create_sparsevecc


  ! **********************************************************************
  ! this will copy a sparse complex vector from SV1 to SV2 ...
  ! interface to =
  ! basically like copy commands for vectors, scalars, BC
  ! check for szie consistency (nCoeff), reallocate output if needed
  ! note that before allocation nCoeff = 0
  ! Remember, SV2 = SV1
  subroutine copy_sparsevecc(SV2,SV1)

    implicit none
    type (sparsevecc), target, intent(in)	:: SV1
    type (sparsevecc), intent(inout)		:: SV2
    integer	                     		:: status

    ! check to see if RHS (SV1) is active (allocated)
    if(.not.SV1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_sparsevecc'
       stop
    else

       ! happen to have the same specs
       if(SV1%nCoeff == SV2%nCoeff) then

          if (SV1%gridType == SV2%gridType) then

             ! just copy the components
             SV2%i = SV1%i
             SV2%j = SV1%j
             SV2%k = SV1%k
             SV2%xyz = SV1%xyz
             SV2%c = SV1%c
	     SV2%grid => SV1%grid

	  else
             write (0, *) 'not compatible usage for copy_sparsevecc'
          end if

          ! SV1 and SV2 do not have the same number of coefficients
       else

	  ! internal memory allocation is strongly discouraged. But this
	  ! is an exception
          if(SV2%allocated) then
             ! first deallocate memory for all the components 
             deallocate(SV2%i, SV2%j, SV2%k, SV2%xyz, SV2%c ,STAT=status)
	     SV2%gridType = ''
       	     SV2%allocated = .false.
          end if

          !  then allocate SV2 as correct szie ...
          Call create_sparsevecc(SV1%nCoeff, SV2, SV1%gridType)
          !   .... and copy SV1
          SV2%nCoeff = SV1%nCoeff
	  SV2%i = SV1%i
	  SV2%j = SV1%j
	  SV2%k = SV1%k
	  SV2%xyz = SV1%xyz
	  SV2%c = SV1%c
	  SV2%grid => SV1%grid

       end if

    end if

  end subroutine copy_sparsevecc


  ! **********************************************************************
  ! linear combination of two sparse vectors, output as a sparse vector
  ! allocates (or reallocates) output sparse vector Loc3
  subroutine linComb_sparsevecc(Lic1,ic1,Lic2,ic2,Loc3)

    implicit none
    type (sparsevecc), target, intent(in)	:: Lic1,Lic2
    type (sparsevecc), intent(inout)		:: Loc3
    complex (kind=8), intent(in)		:: ic1,ic2

    integer					:: n,m,nm
    integer					:: nCoeffSum
    integer, pointer, dimension(:)  	:: Lic1oldIndex
    integer	                     		:: status

    allocate(Lic1oldIndex(Lic2%nCoeff), STAT = status)

    ! it all depends on how many nodes are common
    if(Loc3%allocated) then
       deallocate(Loc3%i, STAT = status)
       deallocate(Loc3%j, STAT = status)
       deallocate(Loc3%k, STAT = status)
       deallocate(Loc3%xyz, STAT = status)
       deallocate(Loc3%c, STAT = status)
       Loc3%gridType = ''
       Loc3%allocated = .false.
    endif

    if (Lic1%gridType == Lic2%gridType) then

       ! count common indices
       nCoeffSum = Lic1%nCoeff+Lic2%nCoeff 
       do m = 1,Lic2%nCoeff
          Lic1oldIndex(m) = 0
          do n = 1,Lic1%nCoeff

             if((Lic1%xyz(n).eq.Lic2%xyz(m)).and.(Lic1%i(n).eq.Lic2%i(m)).and.  &
                  (Lic1%j(n).eq.Lic2%j(m)).and.(Lic1%k(n).eq.Lic2%k(m))) then
                nCoeffSum = nCoeffSum-1
                Lic1oldIndex(m) = n
                exit
             endif

          enddo
       enddo

    else
       write (0, *) 'not compatible usage for LinCompSparseVecC'
    end if

    Call create_sparsevecc(nCoeffsum,Loc3, Lic1%gridType)
    Loc3%grid => Lic1%grid    
    nm = Lic1%nCoeff
    Loc3%i(1:nm) = Lic1%i
    Loc3%j(1:nm) = Lic1%j
    Loc3%k(1:nm) = Lic1%k
    Loc3%xyz(1:nm) = Lic1%xyz
    Loc3%c(1:nm) = ic1*Lic1%c

    do m = 1,Lic2%nCoeff
       ! if none of them are common, just concatenate
       if(Lic1oldIndex(m).eq.0) then
          nm = nm+1
          Loc3%i(nm) = Lic2%i(m)
          Loc3%j(nm) = Lic2%j(m)
          Loc3%k(nm) = Lic2%k(m)
          Loc3%xyz(nm) = Lic2%xyz(m)
          Loc3%c(nm) = ic2*Lic2%c(m)
       else
          Loc3%c(Lic1oldIndex(m)) =  Loc3%c(Lic1oldIndex(m))+ic2*Lic2%c(m)
       endif
    enddo

  end subroutine linComb_sparsevecc

  ! **********************************************************************

  subroutine scMult_sparsevecc(cs,Lin,Lout)
  ! multiply comples sparse vector by a complex scalar cs
  !  output (Lout) can overwrite input (Lin)
  !  allocates for output if necessary

    implicit none
    complex(kind=selectedPrec), intent(in)	:: cs
    type (sparsevecc), intent(in)		:: Lin
    type (sparsevecc), intent(inout)            :: Lout

    ! local variables
    integer                                     :: n,m,nm
    integer                                     :: nCoeffSum
    integer, pointer, dimension(:)          :: Lic1oldIndex
    integer                                     :: status

    !  make sure Lout is allocated and of the correct size
    if(Lout%allocated) then
       if(Lout%nCoeff .ne. Lin%nCoeff) then
          call deall_sparsevecc(Lout)
          call create_sparsevecc(Lin%nCoeff,Lout,Lin%gridType)
       endif
    else
       call create_sparsevecc(Lin%nCoeff,Lout,Lin%gridType)
    endif

    Lout%i = Lin%i
    Lout%j = Lin%j
    Lout%k = Lin%k
    Lout%xyz = Lin%xyz
    Lout%c = cs*Lin%c
    Lout%grid => Lin%grid

    end subroutine scMult_sparsevecc

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  function dotProd_scvector_f(SV,V) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=selectedPrec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_scvector'
       stop
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_scvector'
       stop
    endif

    ! sum over  non-zero terms in sparse vector (conjugate sparse)
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + conjg(SV%c(i)) * V%x(xi, yi, zi)

             ! dealing with y-component 
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + conjg(SV%c(i)) * V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + conjg(SV%c(i)) * V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_scvector'
          stop
       endif

    enddo

  end function dotProd_scvector_f

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  function dotProd_csvector_f(V,SV) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=selectedPrec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_csvector'
       stop
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_csvector'
       stop
    endif

    ! sum over  non-zero terms in sparse vector (conjugate full)
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * conjg(V%x(xi, yi, zi))

             ! dealing with y-component 
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * conjg(V%y(xi, yi, zi))

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * conjg(V%z(xi, yi, zi))
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_csvector'
          stop
       endif

    enddo

  end function dotProd_csvector_f

  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  !   FOR THIS VERSION FIRST VECTOR IS NOT CONJUGATED
  function dotProd_noConj_scvector_f(SV,V) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=selectedPrec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for SparseFullDotProdC'
       stop
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for SparseFullDotProdC'
       stop
    endif

    ! sum over  non-zero terms in sparse vector
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%x(xi, yi, zi)

             ! dealing with y-component 
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for SparseFullDotProdC'
          stop
       endif

    enddo

  end function dotProd_noConj_scvector_f


  ! **********************************************************************
  ! compute complex dot product between a sparse vector SV and a vector of
  ! type cvector ... result in c
  !   FOR THIS VERSION FIRST VECTOR IS NOT CONJUGATED
  function dotProd_noConj_csvector_f(V,SV) result(c)

    implicit none
    type (sparsevecc), intent(in)		:: SV
    type (cvector), intent(in)			:: V
    complex(kind=selectedPrec)				:: c
    integer					:: i
    integer					:: xi, yi, zi

    c = C_ZERO

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_csvector'
       stop
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for dotProd_csvector'
       stop
    endif

    ! sum over  non-zero terms in sparse vector
    ! (need to check xyz the component)
    ! Remember, xyz = 1,2,3 refers to x, y or z components
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%x(xi, yi, zi)

             ! dealing with y-component 
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             c = c + SV%c(i) * V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for dotProd_csvector'
          stop
       endif

    enddo

  end function dotProd_noConj_csvector_f


  ! **********************************************************************
  ! compute sum cs*SV + V where SV is a complex sparse vector. V is a full 
  ! complex vector of type cvector, and cs is a complex scalar; result 
  ! overwrites V. Can also be used to construct a complex full vector 
  ! from a sparse complex vector if cs = (1.0, 0.0) and V = 0 (initially)
  ! or copy from a sparse vector to a full vector 
  subroutine add_scvector(cs,SV,V)

    implicit none
    type (sparsevecc), intent(in)	:: SV
    type (cvector), intent(inout)	:: V
    complex(kind=selectedPrec), intent(in)		:: cs
    integer				:: i
    integer				:: xi, yi, zi

    if((.not.SV%allocated).or.(.not.V%allocated)) then
       write(0,*) 'RHS not allocated yet for add_scvector'
       stop
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for add_scvector'
       stop
    endif

    ! loop over non-zero terms in sparse vector, adding to
    ! corresponding terms in full vector 
    ! (need to check component xyz ...)
    do i = 1,SV%nCoeff

       ! generic test for both edge and face (all the components)
       if ((SV%i(i).le.V%grid%nx+1).or.(SV%j(i).le.V%grid%ny+1).or.&
            (SV%k(i).le.V%grid%nz+1)) then

          ! dealing with x-components
          if (SV%xyz(i) == 1) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             V%x(xi, yi, zi) = cs*SV%c(i) + V%x(xi, yi, zi)

             ! dealing with y-component 
          else if (SV%xyz(i) == 2) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             V%y(xi, yi, zi) = cs*SV%c(i) + V%y(xi, yi, zi)

             ! dealing with z-component
          else if (SV%xyz(i) == 3) then
             xi = SV%i(i)
             yi = SV%j(i)
             zi = SV%k(i)
             V%z(xi, yi, zi) = cs*SV%c(i) + V%z(xi, yi, zi)
          end if

       else
          write(0,*) 'IJK out of bounds for add_scvector'
          stop
       endif

    enddo

  end subroutine add_scvector


  ! **********************************************************************
  ! copy from a full vector to a sparse vector, this routine has quite a 
  ! limited functionality as it assumes one knows the total number of
  ! non-zero vectors in full vector description
  subroutine copy_csvector(SV,V)

    implicit none
    type (sparsevecc), intent(inout)	:: SV
    type (cvector), target, intent(in)	:: V
    integer				:: i
    integer				:: xi, yi, zi

    if(.not.V%allocated) then
       write(0,*) 'RHS not allocated yet for copy_csvector'
       stop
    endif

    if(.not.SV%allocated) then
       write(0,*) 'LHS not allocated yet for copy_csvector'
       stop
    endif

    if (SV%gridType /= V%gridType) then
       write(0,*) 'not compatible usage for copy_csvector'
       stop
    endif

    if (V%gridType == EDGE) then 
    i = 0
    ! for x - component
    do xi = 1, V%nx
       do yi = 1, V%ny+1
	  do zi = 1, V%nz+1
	     if (V%x(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   stop
                end if
                SV%xyz(i) = 1
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%x(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for y - component
    do xi = 1, V%nx+1
       do yi = 1, V%ny
	  do zi = 1, V%nz+1
	     if (V%y(xi, yi, zi) /= 0.0) then
                i = i + 1 
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   stop
                end if
                SV%xyz(i) = 2
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%y(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for z - component
    do xi = 1, V%nx+1
       do yi = 1, V%ny+1
	  do zi = 1, V%nz
	     if (V%z(xi, yi, zi) /= 0.0) then
                i = i + 1 
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   stop
                end if
                SV%xyz(i) = 3
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%z(xi, yi, zi)
             end if
          end do
       end do
    end do

 else if (V%gridType == FACE) then
    ! for x - component
    do xi = 1, V%nx+1
       do yi = 1, V%ny
	  do zi = 1, V%nz
	     if (V%x(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   stop
                end if
                SV%xyz(i) = 1
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%x(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for y - component
    do xi = 1, V%nx
       do yi = 1, V%ny+1
	  do zi = 1, V%nz
	     if (V%y(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   stop
                end if
                SV%xyz(i) = 2
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%y(xi, yi, zi)
             end if
          end do
       end do
    end do
    ! for z - component
    do xi = 1, V%nx
       do yi = 1, V%ny
	  do zi = 1, V%nz+1
	     if (V%z(xi, yi, zi) /= 0.0) then
                i = i + 1
                if (i > SV%nCoeff) then
                   write(0, *) 'outside sparse vector nCoeff: copy_csvector'
                   stop
                end if
                SV%xyz(i) = 3
                SV%i(i) = xi
                SV%j(i) = yi
                SV%k(i) = zi
                SV%c(i) = V%z(xi, yi, zi)
             end if
          end do
       end do
    end do

 else 

    write (0, *) 'Vector (full) use not proper in copy_csvector' 

 end if 
 
 SV%grid => V%grid

end subroutine copy_csvector

end module sg_sparse_vector
