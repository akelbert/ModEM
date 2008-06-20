module grid2d
   !  Defines grid data structure, basic grid methods
   
   use math_constants

   ! Possible grid types, on which cvector is defined. Viewing the grid
   ! as 2D, NODES and EDGES correspond to the corners and the edges of
   ! all square grid elements, and CELLS corresponds to the centers of
   ! these squares. EARTH means the air is excluded from the grid.
   character(len=80), parameter		:: CELL = 'CELL' 
   character(len=80), parameter		:: NODE = 'NODE'
   character(len=80), parameter		:: EDGE = 'EDGE'
   character(len=80), parameter		:: NODE_EARTH = 'NODE EARTH'
   character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'
   character(len=80), parameter		:: EDGE_EARTH = 'EDGE EARTH'

   type ::  grid2d_t
      !  grid2d_t is derived data type used to store grid geometry
      integer   :: Nz = 0
      integer   :: Ny = 0
      integer   :: Nza = 0
      real (kind=selectedPrec), pointer, dimension(:) :: Dy,Dz
      real (kind=selectedPrec), pointer, dimension(:) :: Dely,Delz
      real (kind=selectedPrec), pointer, dimension(:) :: yNode,zNode
      real (kind=selectedPrec), pointer, dimension(:) :: yCenter,zCenter
   end type grid2d_t

   Contains

     !************************************************************************
     subroutine create_Grid2D(Ny,Nz,Nza,grid)
       !  creates finite differences grid2d_t structure of 
       !  size Nz x Ny, allocates arrays
       !
       implicit none
       integer, intent(in)		:: Nz,Ny,Nza
       type (grid2d_t) , intent(inout)	:: grid

       grid%Nza = Nza
       grid%Nz = Nz
       grid%Ny = Ny
       allocate(grid%Dz(Nz))
       allocate(grid%Dy(Ny))
       allocate(grid%Delz(Nz+1))
       allocate(grid%Dely(Ny+1))
       allocate(grid%zNode(Nz+1))
       allocate(grid%yNode(Ny+1))
       allocate(grid%zCenter(Nz))
       allocate(grid%yCenter(Ny))

     end subroutine create_Grid2D 
     
     !************************************************************************
     subroutine deall_Grid2D(grid)
       !  deallocates finite differences grid2d_t structure
       !
       implicit none
       type (grid2d_t) , intent(inout)	:: grid

       if (associated(grid%Dz)) deallocate(grid%Dz)
       if (associated(grid%Dy)) deallocate(grid%Dy)
       if (associated(grid%Delz)) deallocate(grid%Delz)
       if (associated(grid%Dely)) deallocate(grid%Dely)
       if (associated(grid%zNode)) deallocate(grid%zNode)
       if (associated(grid%yNode)) deallocate(grid%yNode)
       if (associated(grid%zCenter)) deallocate(grid%zCenter)
       if (associated(grid%yCenter)) deallocate(grid%yCenter)

     end subroutine deall_Grid2D 
     
     !***********************************************************************
     ! gridCalcs grid: after allocation, read in Dy, Dz, then call:
     subroutine gridCalcs(grid)
       implicit none
       type (grid2d_t) , intent(inout)	:: grid
       !  local variables
       integer ::	Nz,Ny,iy,iz,Nza

       Nz = grid%Nz
       Nza = grid%Nza
       Ny = grid%Ny
       do iy = 2,Ny
          grid%Dely(iy) = (grid%Dy(iy-1)+grid%Dy(iy))/2.
       enddo
       grid%Dely(1) = grid%Dely(1)/2.
       grid%Dely(Ny+1) = grid%Dely(Ny)/2.
       do iz = 2,Nz
            grid%Delz(iz) = (grid%Dz(iz-1)+grid%Dz(iz))/2.
        enddo
        grid%Delz(1) = grid%Delz(1)/2.
        grid%Delz(Nz+1) = grid%Delz(Nz)/2.

        !  for now assume y0 = 0, z0 = 0 (left edge, Earth surface)
        grid%yNode(1) = 0
        grid%zNode(1) = 0
        do iy = 1,Ny
           grid%yNode(iy+1) = grid%yNode(iy) + grid%Dy(iy)
        enddo
        ! start summing from top of domain (includes air)
        do iz = 1,Nz
           grid%zNode(iz+1) = grid%zNode(iz) + grid%Dz(iz)
        enddo
        !  next construct positions of cell centers
        do iy = 1,Ny
           grid%yCenter(iy) = (grid%yNode(iy)+grid%yNode(iy+1))/2.;
        enddo
        do iz = 1,grid%Nz
           grid%zCenter(iz) = (grid%zNode(iz)+grid%zNode(iz+1))/2.;
        enddo
     end subroutine gridCalcs
end module grid2D
