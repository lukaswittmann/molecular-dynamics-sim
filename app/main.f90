program main
  
use types, only: dp
use forcefields
use utils
use energy_minimization
implicit none

integer, parameter :: nparticles = 50
real(dp), parameter :: dt = 0.001, t = 100., cutoff = 8., sigma = 1., epsilon = 1.  ! sigma = (P = 0), epsilon = depth
integer, parameter :: maxiter = INT(t/dt + 1)
integer :: beginning, rate, end, iter, particle1, particle2, interaction, particle, dim, dim1, dim2, dim3
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
v(:,:) = rand_v(nparticles, 1)
r(:,:,-1) = r(:,:,0) - v(:,:) * dt
!call print_rv(r(:,:,0:0), v, nparticles) 

properties(:,:,1) = 1. ! Set mass

! Energy minimization
print '(/,A,/)', "** ENERGY MINIMIZATION **"
r = minimize_energy(r, nparticles, properties, boxsize, sigma, cutoff)

! Basic md-run information
print *, "t =", t, "dt =", dt
print *, "maxiter = ", maxiter, "nparticles =", nparticles
call system_clock(beginning, rate)
print '(/,A,/)', "** MD RUN **"

! -------------------------------------------------------------------------------------------- !

! Integration loop
integrator: do iter = 0, maxiter

  ! "Particle-teleporter" to fulfill periodic boundary condition
  do particle = 1, nparticles ! Iterate over all particles
    do dim = 1, 3 ! Iterate over 3d space

      if (r(dim, particle, 0) > boxsize(dim)) then
        r(dim, particle, :) = r(dim, particle, :) - boxsize(dim)
      else if (r(dim, particle, 0) < 0.0_dp) then
        r(dim, particle, :) = r(dim, particle, :) + boxsize(dim)
      end if

    end do
  end do

  ! Calculation of particle interactions
  f = 0
  do particle1 = 1, nparticles   ! Interaction for each particle
    do particle2 = 1, nparticles ! ..with each other particle
      if (particle1 /= particle2) then  

        d0 = r(:, particle2, 0) - boxsize(:) - r(:, particle1, 0)

        do dim1 = 1,3
          do dim2 = 1,3
            do dim3 = 1,3

              ! Define new position of image and for dim1=dim2=dim3=1 of real particle interaction
              d = d0
              d(1) = d(1) + (dim1 - 1) * boxsize(1)
              d(2) = d(2) + (dim2 - 1) * boxsize(2)
              d(3) = d(3) + (dim3 - 1) * boxsize(3)

              l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)

              if (l < cutoff) then
                f(:,particle1,1) = f(:,particle1,1) + lj_f(d, l, sigma, epsilon)
              end if

            end do
          end do
        end do

      end if
    end do
  end do

  ! Do Verlet integration step
  r(:,:,1) = 2 * r(:,:,0) - r(:,:,-1) + (f(:,:,1) / properties(:,:,1)) * (dt ** 2)

  ! Calculate velocities by numerical derivation
  v = (r(:,:,1) - r(:,:,-1)) / (2 * dt) ! v(t+dt)

  ! Move array to t = t+dt for next integration step
  r = cshift(r,  1, DIM=3)

  if (mod(iter, 1000)==0) then
    ! Calculate temperature (kinetic energy) in [a.u]
    print *, "Energy =",  (properties(1,1,1) / 2 * sum(v ** 2)), (real(iter)/real(maxiter)*100), "%"
  end if
    
  ! Save trajectories, velocities and forces every .. iterations
  !if (mod(iter, 20)==0) then
  !  do particle = 1, nparticles
  !    write (file_name,'(i0)') particle
  !    call save_trajectory("export/r_"//trim(file_name)//".txt", r(:,particle,0))
  !    call save_trajectory("export/v_"//trim(file_name)//".txt", v(:,particle))
  !    call save_trajectory("export/f_"//trim(file_name)//".txt", f(:,particle,1))
  !  end do
  !end if

  ! Checks if simulation failed every .. iterations
  if (mod(iter, 1000)==0) then
    if ((maxval(r(:,:,-1)) > maxval(boxsize(:))) .or. (minval(r(:,:,-1)) < 0)) then
      print '(/,A,/)', "ERROR: Calculation out of bounds, aborting..."
      exit integrator
    end if
  end if

end do integrator

! -------------------------------------------------------------------------------------------- !
  
!call print_rv(r, v, nparticles) 

call system_clock(end)
print *, "Elapsed time: ", real(end - beginning) / real(rate), "seconds."

end program main
