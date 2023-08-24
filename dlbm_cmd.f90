module dlbm_cmd
use dlbm, only: wp
implicit none

type :: cmd_settings
  integer :: N 
    ! Number of lattice cells
  real(wp) :: length
    ! Initial length (half-thickness or radius)
  real(wp) :: X0
    ! Initial moisture concentration on a dry basis
  real(wp) :: Xeq
    ! Equilibrium moisture concentration on a dry basis

  real(wp) :: max_time

  integer :: nsteps
  integer :: nout
end type

contains

end module