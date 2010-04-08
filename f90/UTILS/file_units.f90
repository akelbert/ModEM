! *****************************************************************************
module file_units
  ! This module defines I/O file unit information, used throughout the code

  implicit none

  ! User-defined parameter that sets the level of output
  integer, save									:: output_level=0

  ! Startup and control files
  integer, parameter							:: ioStartup=101
  integer, parameter							:: ioFwdCtrl=105
  integer, parameter							:: ioInvCtrl=106

  ! Grid, model, data
  integer, parameter							:: ioPrm=1
  integer, parameter							:: ioGrd=2
  integer, parameter							:: ioDat=3
  integer, parameter							:: ioMdl=4

  ! Generic error, read and write units
  integer, parameter							:: ioERR=7
  integer, parameter							:: ioREAD=8
  integer, parameter							:: ioWRITE=9

  ! Dictionary files, if they exist
  integer, parameter							:: ioRX=17
  integer, parameter							:: ioTX=18
  integer, parameter							:: ioDT=19

  ! Fields, currents etc
  integer, parameter	  				        :: ioH=10
  integer, parameter	  				        :: ioE=11
  integer, parameter	  				        :: ioJ=12

  ! Others, as necessary
  integer, parameter							:: ioShell=28

  integer, parameter							:: ioRad=15
  integer, parameter					        :: ioSens=22
  integer, parameter                            :: ioC=23, ioD=24

end module file_units
