program main
  
use types, only: dp
use forcefields
use utils
implicit none

integer, parameter :: nparticles = 50
real(dp), parameter :: dt = 0.1, t = 500.
integer, parameter :: maxiter = INT(t/dt + 1)
integer :: beginning, rate, end, iter, particle1, particle2, interaction, particle
character (len=10) :: file_name

real(dp), dimension(3, nparticles, 1) :: properties ! 1: mass  2: ..   3: .. [x,y,z] for anisotropy
real(dp), dimension(3, nparticles, -1:2) :: r       ! ([x,y,z], particle, [t-dt, t, t+dt, t+2dt])
real(dp), dimension(3, nparticles) :: v = 0         ! ([x,y,z], particle)
real(dp), dimension(3, nparticles,1) :: f = 0       ! ([f(x),f(y),f(z)], particle)


! Generate random starting positions and initial velocities
r(:,:,:) = rand_r(nparticles, 10) 
!v(:,:) = rand_v(nparticles, 100)
!r(:,:,-1) = r(:,:,0) - v(:,:) * dt
call print_final_rv(r(:,:,0:0), v, nparticles) 


properties(:,:,1) = 1. ! Set mass

! Basic md-run information
print *, "t =", t, "dt =", dt, "maxiter = ", maxiter
call system_clock(beginning, rate)


! Integration loop
integrator: do iter = 0, maxiter

f = 0
do particle1 = 1, nparticles          ! Iteration over each particle
  do particle2 = 1, nparticles        ! Interaction with each other particle
    if (particle1 /= particle2) then  ! No selfinteraction
      f(:,particle1,1) = f(:,particle1,1) + lj_f(particle1, particle2, nparticles, r, 1.6_dp, 1.0_dp)
    end if  
  end do
end do

  ! Do Verlet integration step using
  r(:,:,1) = 2 * r(:,:,0) - r(:,:,-1) + (f(:,:,1) / properties(:,:,1)) * (dt ** 2)

  ! Calculate velocities by numerical derivation
  v = (r(:,:,2) - r(:,:,0)) / (2 * dt)

  ! Move array to t = t+dt for next integration step
  r = cshift(r,  1, DIM=3)
    
! Save trajectories and velocities
!if (mod(iter, 100)==0) then
!  do particle = 1, nparticles
!    write (file_name,'(i0)') particle
!    call save_trajectory("export/r_"//trim(file_name)//".txt", r(:,particle,0))
    !call save_trajectory("export/v_1.txt", v(:,particle))
    !call save_trajectory("export/f_1.txt", f(:,particle,1))
!  end do
!end if

end do integrator
  
print *, " "
call print_final_rv(r, v, nparticles) 

call system_clock(end)
print *, "Elapsed time: ", real(end - beginning) / real(rate), "seconds."


end program main
