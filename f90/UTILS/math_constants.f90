! *****************************************************************************
! math_constants.f90 defines math constants and other parameters

module math_constants
  implicit none

  !  use this to select singal or double precision
  integer, parameter                    :: SP = selected_real_kind(6,37)
  integer, parameter                    :: DP = selected_real_kind(15,307)
  integer, parameter			:: selectedPrec = DP

  real (kind=selectedPrec), parameter	:: PI = 3.14159265_selectedPrec
  real (kind=selectedPrec), parameter	:: MU = PI*.0000004_selectedPrec

  ! Conductivity of the air for computational purposes
  real (kind=selectedPrec), parameter	:: SIGMA_AIR = 1.0e-10

  ! Useful conversion constants
  real (kind=selectedPrec), parameter	:: D2R = PI/180._selectedPrec 
  real (kind=selectedPrec), parameter	:: R2D = 180._selectedPrec/PI 
  real (kind=selectedPrec), parameter	:: KM2M = 1000.0_selectedPrec
  real (kind=selectedPrec), parameter	:: M2KM = 0.001_selectedPrec
  
  ! Sign convention
  integer, parameter                    :: isign = -1

  real(kind=selectedPrec),parameter	:: TWO = 2.0_selectedPrec
  real(kind=selectedPrec),parameter	:: ONE = 1.0_selectedPrec
  real(kind=selectedPrec),parameter	:: THREE = 3.0_selectedPrec
  real(kind=selectedPrec),parameter	:: EIGHT = 8.0_selectedPrec
  real(kind=selectedPrec),parameter	:: R_ZERO = 0.0_selectedPrec
  real(kind=selectedPrec),parameter	:: MinusONE = -1.0_selectedPrec
  real(kind=selectedPrec),parameter	:: MinusTWO = -2.0_selectedPrec
  
  ! Real precision constants / tolerance
  real(kind=selectedPrec),parameter	:: LARGE_REAL = 1.0e20
  real(kind=selectedPrec),parameter :: TOL4= 0.0001_dp
  real(kind=selectedPrec),parameter :: TOL6= 0.000001_dp
  real(kind=selectedPrec),parameter :: TOL8= 0.00000001_dp

  complex(kind=selectedPrec), parameter	:: C_ONE = (1.0_selectedPrec,0.0_selectedPrec)
  complex(kind=selectedPrec), parameter :: C_ZERO = (0.0_selectedPrec, 0.0_selectedPrec)
  complex(kind=selectedPrec), parameter	:: C_MinusOne = (-1.0_selectedPrec, 0.0_selectedPrec)

  character(len=80), parameter		:: CELL = 'CELL'
  character(len=80), parameter		:: NODE = 'NODE'
  character(len=80), parameter		:: FACE = 'FACE'
  character(len=80), parameter		:: EDGE = 'Edge'
  character(len=80), parameter		:: CENTER = 'Center'
  character(len=80), parameter		:: Corner = 'Corner'
  character(len=80), parameter		:: NODE_EARTH = 'NODE EARTH'
  character(len=80), parameter		:: FACE_EARTH = 'FACE EARTH'
  character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'
  character(len=80), parameter		:: LOG_CELL = 'LOG_CELL'

  character(len=2), parameter		:: TE = 'TE'
  character(len=2), parameter		:: TM = 'TM'

  character(len=3), parameter		:: FWD = 'FWD'
  character(len=3), parameter		:: TRN = 'TRN'
  character(len=3), parameter		:: ADJ = 'ADJ'

end module math_constants !math_constants
