module forcefields

use types, only: dp
implicit none

contains

    
function spring_f(d, l, l0, force_const)
  real(dp), intent(in) ::l, l0, force_const
  real(dp), dimension(3), intent(in) :: d
  real(dp), dimension(3) :: spring_f

  spring_f(:) = (force_const * (l - l0)) * (d(:)/l)

end function spring_f
    

function lj_f(d, l, sigma, epsilon)
  real(dp), intent(in) :: epsilon, sigma, l
  real(dp), dimension(3), intent(in) :: d
  real(dp), dimension(3) :: lj_f

  lj_f(:) = 4 * epsilon * (-((12 * sigma ** 12)/(l ** 13) + ((6 * sigma ** 6)/(l ** 7)))) * (d(:)/l)  
end function lj_f


function min_f(d, l, l0)
  real(dp), intent(in) :: l, l0 ! l0 gives desired mean distance
  real(dp), dimension(3), intent(in) :: d
  real(dp), dimension(3) :: min_f

  min_f(:) = 4 * (-((12 * l0 ** 12)/(l ** 13) + ((6 * l0 ** 6)/(l ** 7)))) * (d(:)/l)  
end function min_f


function fitted_lj_f(d, r_in, a) !a0, a1, a2, a3, a4, a5, a6)
  ! Potential from Two-body interatomic potentials for He, Ne, Ar, Kr, and Xe from ab 
  ! initio data, K. Deiters1 and R. J. Sadus, J. Chem. Phys. 150, 134504 (2019)
  real(dp), intent(in), dimension(0:6) :: a ! a0, a1, a2, a3, a4, a5, a6
  real(dp), intent(in) :: r_in ! input r in m
  real(dp) :: r, kB = 0.13806504e-22
  real(dp), dimension(3), intent(in) :: d
  real(dp), dimension(3) :: fitted_lj_f

  r = r_in * 1E9 ! Conversion to nm

  ! Force derivative of SAAPx potential
  ! Force converted into Newton
  fitted_lj_f(:) = -(9*(a(0)*exp(a(6)*r**2 + a(1)*r)/r + a(2)*exp(a(3)*r) + a(4))*kB*a(5)*r**5) / ((a(5)*r**6 + 1)**2) &
  + (3*kB*(-a(0)*exp(a(6)*r**2 + a(1)*r)/r**2 + a(0)*(2*a(6)*r + a(1))*exp(a(6)*r**2 + a(1)*r)/r &
  + a(2)*a(3)*exp(a(3)*r)))/(2*(a(5)*r**6 + 1)) * (d(:) / r)  
end function fitted_lj_f


end module forcefields