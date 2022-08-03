module utils

use types, only: dp
implicit none

contains

subroutine savetxt(filename, d)
    character(len=*), intent(in) :: filename  ! File to save the array to
    real(dp), intent(in) :: d(:, :)           ! The 2D array to save

    integer :: s, i
    open(newunit=s, file=filename, status="replace")
    do i = 1, size(d, 1)
        write(s, *) d(i, :)
    end do
    close(s)
end subroutine savetxt

subroutine save_md_array(filename, d) ! Saves each particles coords in each line
    character(len=*), intent(in) :: filename  ! File to save the array to
    real(dp), intent(in) :: d(:, :)           ! The 2D array to save

    integer :: s, i
    open(newunit=s, file=filename, status="replace")
    do i = 1, size(d, 2)
        write(s, *) d(:, i) ! xyz, particle
    end do
    close(s)
end subroutine save_md_array

subroutine print_rv(r, v, nparticles)
    integer, intent(in) :: nparticles
    real(dp), dimension(3,nparticles), intent(in) :: r, v
    print *, '=============================== Positions:  =============================='
    print *, r
    print *, '=============================== Velocities: =============================='
    print *, v
    print *, '=========================================================================='

end subroutine print_rv

subroutine save_trajectory(filename, d) ! Saves each particles coords in each line
    character(len=*), intent(in) :: filename  ! File to save the array to
    real(dp), intent(in) :: d(:)           ! The 2D array to save

    integer :: s, i
    open(newunit=s, file=filename, access="sequential",  position="append")

    !do i = 1, size(d, 1)
        write(s, '(*(G0.6,:,"'//achar(9)//'"))') d(:) ! xyz, particle
    !end do

    close(s)
end subroutine save_trajectory

! -------------------------------------------------------------------------------------------- !

function rand_r(nparticles, boxsize)
  ! Generate random starting positions
    integer, intent(in) :: nparticles
    real(dp), dimension(3) :: boxsize
    real(dp), dimension(3,nparticles,-1:2) :: rand_r
    integer :: iter
    real(dp) :: randx, randy, randz
    do iter = 1, nparticles
        call random_number(randx)
        call random_number(randy)
        call random_number(randz)
        rand_r(1, iter, :) = randx * boxsize(1)
        rand_r(2, iter, :) = randy * boxsize(2)
        rand_r(3, iter, :) = randz * boxsize(3)
    end do   
end function rand_r

function rand_v(nparticles, factor)
    ! Generate random starting positions
      integer, intent(in) :: nparticles, factor
      real(dp), dimension(3,nparticles) :: rand_v
      integer :: iter
      real(dp) :: randvx, randvy, randvz
      do iter = 1, nparticles
          call random_number(randvx)
          call random_number(randvy)
          call random_number(randvz)
          rand_v(1, iter) = (randvx - 0.5) * (factor / 100)
          rand_v(2, iter) = (randvy - 0.5) * (factor / 100)
          rand_v(3, iter) = (randvz - 0.5) * (factor / 100)
      end do   
end function rand_v

end module utils
