module constants

use types, only: dp
implicit none

private
public pi, e_, i_, kB, R_, NA, mass, ff_parameters

real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e_   = 2.7182818284590452353602874713527_dp
real(dp), parameter :: kB   = 1.380649E-23_dp
real(dp), parameter :: NA   = 6.02214076E-23_dp
real(dp), parameter :: R_   = 8.314472_dp
complex(dp), parameter :: i_ = (0, 1)

! Mass
!                                                He      Ne      Ar      Kr      Xe         conversion to kg
real(dp), dimension(1:5), parameter :: mass = (/4.0026, 20.180, 39.948, 83.798, 131.29/) * 1.660538782E-27

! Potential/Force parameters
real(dp), dimension(0:6,1:5), parameter :: ff_parameters = reshape((/ &
! He
24238.01564, -5.934943188, -11.17928721, -1.821078761, -0.7762283939, 0.3703465767, -3.210773756, &
! Ne
81648.44026, -9.829947619, -6.482831445, -0.5073208921, -0.4906026951, 0.9921472732, 0., &
! Ar
65214.64725, -9.452343340, -19.42488828, -1.958381959, -1.854049555, 0.7454617542, 0., &
! Kr
60249.13228, -9.456080572, -24.40996013, -2.18227961, -1.959180470, 0.874092399, 0., &
! Xe
44977.3164, -9.121814449, -29.63636182, -2.278991444, -1.876430370, 0.8701531593, 0./), (/7,5/))


end module
