module forcefields

use types, only: dp
implicit none
!private
!public spring_f, lj_f, ff
contains
    
    
function spring_f(particle1, particle2, nparticles, r, force_const)
  integer, intent(in) :: nparticles, particle1, particle2
  real(dp), intent(in) :: force_const
  real(dp), dimension(3, nparticles, -1:2), intent(in) :: r

  !real(dp), dimension(3, nparticles, 1) :: calc_f
  real(dp) :: l, l0, f
  real(dp), dimension(3) :: d
  real(dp), dimension(3) :: spring_f

  l0 = 1. ! Equilibrium lengh

  ! Calculate 3d-distances
  d = r(:, particle2, 0) - r(:, particle1, 0)

  ! Calculate net distance between particles
  l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)

  ! Calculate force (spring)
  f = ( force_const * (l - l0))
  spring_f(:) = f * (d(:)/l)

end function spring_f
    
    
function lj_f(particle1, particle2, nparticles, r, sigma, epsilon)
  integer, intent(in) :: nparticles, particle1, particle2
  real(dp), intent(in) :: epsilon, sigma ! sigma = l0, epsilon = depth
  real(dp), dimension(3, nparticles, -1:2), intent(in) :: r

  !real(dp), dimension(3, nparticles, 1) :: calc_f
  real(dp) :: l, l0, f
  real(dp), dimension(3) :: d
  real(dp), dimension(3) :: lj_f

  l0 = 1. ! Equilibrium lengh

  ! Calculate 3d-distances
  d = r(:, particle2, 0) - r(:, particle1, 0)

  ! Calculate net distance between particles
  l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)

  ! Calculate force lj-potential
  f = 4 * epsilon * (-((12 * sigma ** 12)/(l ** 13) + ((6 * sigma ** 6)/(l ** 7))))
  lj_f(:) = f * (d(:)/l)
end function lj_f



function ff(particle1, particle2, nparticles, r, l, sigma, epsilon)
  integer, intent(in) :: nparticles, particle1, particle2

  real(dp), intent(in) :: epsilon, sigma, l ! sigma = (P = 0), epsilon = depth
  real(dp), dimension(3, nparticles, -1:2), intent(in) :: r

  !real(dp), dimension(3, nparticles, 1) :: calc_f
  real(dp) :: l0
  real(dp), dimension(3) :: d
  real(dp), dimension(3) :: ff

  l0 = 1. ! Equilibrium lengh

  ! Calculate xyz-distances
  d = r(:, particle2, 0) - r(:, particle1, 0)

  ! Calculate net distance between particles: is done in main iterationloop

  ! Decide if l < cutoff is done in main iterationloop

  ff(:) = 4 * epsilon * (-((12 * sigma ** 12)/(l ** 13) + ((6 * sigma ** 6)/(l ** 7)))) * (d(:)/l)


  
end function ff      

    
end module forcefields