! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: Global 3D

  use data_vectors
  use global
  use input
  use output

  implicit none

  private

  public     :: write_dataVecMTX


Contains

!**********************************************************************
! writes global responses file in ASCII format

   subroutine write_dataVecMTX(allData,cfile)

      character(*), intent(in)					:: cfile
      type(dataVecMTX_t), intent(in)			:: allData
      ! local variables
      type (transmitter_t)						:: freq
      integer									:: ifreq
      character(80)								:: strtmp

      do ifreq=1,freqList%n

		freq = freqList%info(ifreq)

		write(strtmp,*) trim(cfile),'.cout'
		outFiles%fn_cdat = strtmp

		write(strtmp,*) trim(cfile),'.dout'
		outFiles%fn_ddat = strtmp

		call outputResponses(freq,allData,freqList,TFList,obsList,outFiles,dat)

      end do

   end subroutine write_dataVecMTX


end module DataIO
