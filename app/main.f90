program main
  
use types, only: dp
use forcefields
use utils
implicit none

integer, parameter :: nparticles = 20
real(dp), parameter :: dt = 0.001, t = 100., cutoff = 10., sigma = 1.6, epsilon = 0.1
integer, parameter :: maxiter = INT(t/dt + 1)
integer :: beginning, rate, end, iter, particle1, particle2, interaction, particle, dim
character (len=10) :: file_name

real(dp), dimension(3) :: boxsize = (/10.,10.,10./) ! define xyz boxsize
real(dp), dimension(3, nparticles, 1) :: properties ! 1: mass  2: ..   3: .. [x,y,z] for anisotropy
real(dp), dimension(3, nparticles, -1:2) :: r       ! ([x,y,z], particle, [t-dt, t, t+dt, t+2dt])
real(dp), dimension(3, nparticles) :: v = 0         ! ([x,y,z], particle)
real(dp), dimension(3, nparticles,1) :: f = 0       ! ([f(x),f(y),f(z)], particle)
real(dp), dimension(3) :: d, d0
real(dp) :: l ! Needed for cutoff and calculation of force

! Generate random starting positions in boxsize and initial velocities
r(:,:,:) = rand_r(nparticles, boxsize) 
!v(:,:) = rand_v(nparticles, 100)
!r(:,:,-1) = r(:,:,0) - v(:,:) * dt
call print_final_rv(r(:,:,0:0), v, nparticles) 



properties(:,:,1) = 1. ! Set mass

! Basic md-run information
print *, "t =", t, "dt =", dt
print *, "maxiter = ", maxiter, "nparticles =", nparticles
call system_clock(beginning, rate)

! -------------------------------------------------------------------------------------------- !

! Integration loop
integrator: do iter = 0, maxiter


! Particle teleporter to fulfill periodic boundary condition
do particle = 1, nparticles ! Iterate over all particles
  do dim = 1, 3 ! Iterate over all 3d space

    if (r(dim, particle, 0) > boxsize(dim)) then
      r(dim, particle, :) = r(dim, particle, :) - boxsize(dim)
    else if (r(dim, particle, 0) < 0.0_dp) then
      r(dim, particle, :) = r(dim, particle, :) + boxsize(dim)
    end if

  end do
end do


f = 0
do particle1 = 1, nparticles   ! Interaction for each particle
  do particle2 = 1, nparticles ! ..with each other particle
    if (particle1 /= particle2) then  

      ! Realparticle interaction
      d0 = r(:, particle2, 0) - r(:, particle1, 0) ! Realspace
      d = d0

      ! Calculate distance, decide if l < cutoff
      l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)
      if (l < cutoff) then
        f(:,particle1,1) = f(:,particle1,1) &
        + ff(particle1, particle2, nparticles, d, l, sigma, epsilon)
      end if

      do dim = 1, 3
        d = d0 ! Realspace
        d(dim) = d(dim) + boxsize(dim) ! Mirrorparticle interactions +xyz

        l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)
        if (l < cutoff) then
          f(:,particle1,1) = f(:,particle1,1) &
          + ff(particle1, particle2, nparticles, d, l, sigma, epsilon)
        end if
      end do

      do dim = 1, 3
        d = d0 ! Realspace
        d(dim) = d(dim) - boxsize(dim) ! Mirrorparticle interactions -xyz
        
        l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)
        if (l < cutoff) then
          f(:,particle1,1) = f(:,particle1,1) &
          + ff(particle1, particle2, nparticles, d, l, sigma, epsilon)
        end if
      end do
       
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
if (mod(iter, 100)==0) then
  do particle = 1, nparticles
    write (file_name,'(i0)') particle
    call save_trajectory("export/r_"//trim(file_name)//".txt", r(:,particle,0))
    !call save_trajectory("export/v_1.txt", v(:,particle))
    !call save_trajectory("export/f_1.txt", f(:,particle,1))
  end do
end if

end do integrator

! -------------------------------------------------------------------------------------------- !

  
print *, " "
call print_final_rv(r, v, nparticles) 

call system_clock(end)
print *, "Elapsed time: ", real(end - beginning) / real(rate), "seconds."


end program main
