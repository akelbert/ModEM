! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Generic interface required to call this from Level I inversion routines...
  ! Version: Global 3D

  use UserData
  use input
  use output

  implicit none

  private

  public     :: write_dataVectorMTX


Contains

!**********************************************************************
! writes global responses file in ASCII format

   subroutine write_dataVectorMTX(allResp,cfile)

      character(*), intent(in)                      :: cfile
      type(dataVectorMTX_t), intent(in)             :: allResp
      ! local variables
      integer									:: i
      character(80)								:: strtmp

      do i=1,allData%nTx

		write(strtmp,*) trim(cfile),'.cout'
		outFiles%fn_cdat = strtmp

		write(strtmp,*) trim(cfile),'.dout'
		outFiles%fn_ddat = strtmp

		call outputResponses(allResp%d(i),outFiles,allData%d(i))

      end do

   end subroutine write_dataVectorMTX


end module DataIO
