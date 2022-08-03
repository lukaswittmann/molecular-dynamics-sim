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


end module forcefields