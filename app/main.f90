program main
  
use types
use forcefields
use utils
use energy_minimization
use constants
implicit none

! Md run settings
integer, parameter :: nparticles = 512
real(dp), parameter :: dt = (2 * 1E-15), t = (10. * 1E-12) ! dt in fs to s, t in ps to s
integer, parameter :: maxiter = INT(t/dt)
real(dp), dimension(3) :: boxsize = (/3.2,3.2,3.2/) * 1E-9 ! Boxsize in nm to meter
real(dp), dimension(3) :: boxsize_half = (/1.6,1.6,1.6/) * 1E-9 ! Cutoff in nm to m
real(dp) :: cutoff = 1.6

! Particle settings
integer, parameter :: atom_type = 3 ! 1: He, 2: Ne, 3: Ar, 4: Kr. 5: Xe
! real(dp), parameter :: sigma = 3.40100e-01 , epsilon = 9.78638e-01 ! Ar

! Thermostat settings
real(dp), parameter :: tau_T = (1 * 1E-12) ! relaxation time of thermostat in ps to s
real(dp) :: E, T_inst, T0 = 300. ! in K

! Barostat settings
real(dp), parameter :: tau_P =  (1. * 1E-8) ! relaxation time of thermostat in ns to s
real(dp) :: Xi, mu, V_inst, P_inst, P0 = (1. * 1E5)  ! in Bar to Pa

! Runtime variables
integer ::  beginning, rate, end, iter, particle1, particle2, interaction, particle, dim, dim1, dim2, dim3
character (len=10) :: file_name
real(dp), dimension(3, nparticles, 1) :: properties ! 1: mass  2: ..   3: .. [x,y,z] for anisotropy
real(dp), dimension(3, nparticles, -1:1) :: r       ! ([x,y,z], particle, [t-dt, t, t+dt])
real(dp), dimension(3, nparticles) :: v = 0         ! ([x,y,z], particle)
real(dp), dimension(3, nparticles, nparticles, 1) :: r_ij = 0 ! Distance triangle matrix
real(dp), dimension(3, nparticles, nparticles, 1) :: f_ij = 0 ! Force triangle matrix
real(dp), dimension(3) :: d, d0
real(dp) :: timing, l

! Generate random starting positions in boxsize and initial velocities
r(:,:,:) = rand_r(nparticles, boxsize) 
v(:,:) = rand_v(nparticles, int(1 * 1E4))
r(:,:,-1) = r(:,:,0) - v(:,:) * dt
!call print_rv(r(:,:,0:0), v, nparticles) 

! Energy minimization
!print '(/,A,/)', "** ENERGY MINIMIZATION **"
!r = minimize_energy(r, nparticles, properties, boxsize, 5E-11_dp, boxsize_half)

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

! Calculation of particle distances with pbc and interactions
r_ij = 0
f_ij = 0
do particle1 = 1, nparticles   ! Interaction for each particle
  do particle2 = (particle1 + 1), nparticles ! ..with each other particle

    ! Calculate xyz distances of particle1 and 2
    r_ij(:, particle1, particle2, 1) = r(:, particle2, 0) - r(:, particle1, 0)

    ! Dimension x
    if (r_ij(1, particle1, particle2, 1) < -1 * boxsize_half(1)) then
      r_ij(1, particle1, particle2, 1) = r_ij(1, particle1, particle2, 1) + boxsize_half(1)
    else if (r_ij(1, particle1, particle2, 1) > boxsize_half(1)) then
      r_ij(1, particle1, particle2, 1) = r_ij(1, particle1, particle2, 1) - boxsize_half(1)
    end if

    ! Dimension y
    if (r_ij(2, particle1, particle2, 1) < -1 * boxsize_half(2)) then
      r_ij(2, particle1, particle2, 1) = r_ij(2, particle1, particle2, 1) + boxsize_half(2)
    else if (r_ij(2, particle1, particle2, 1) > boxsize_half(2)) then
      r_ij(2, particle1, particle2, 1) = r_ij(2, particle1, particle2, 1) - boxsize_half(2)
    end if

    ! Dimension z
    if (r_ij(3, particle1, particle2, 1) < -1 * boxsize_half(3)) then
      r_ij(3, particle1, particle2, 1) = r_ij(3, particle1, particle2, 1) + boxsize_half(3)
    else if (r_ij(3, particle1, particle2, 1) > boxsize_half(3)) then
      r_ij(3, particle1, particle2, 1) = r_ij(3, particle1, particle2, 1) - boxsize_half(3)
    end if

    ! Mirror distances to other triangle matrix
    r_ij(:, particle2, particle1, 1) = - r_ij(:, particle1, particle2, 1)

    ! Calculate total distance for cutoff
    l = SQRT(r_ij(1, particle1, particle2, 1) ** 2 + r_ij(2, particle1, particle2, 1) ** 2 + r_ij(3, particle1, particle2, 1) ** 2)

    ! Calculate interactions if l < cutoff
    if (l < cutoff) then
      f_ij(:, particle1, particle2, 1) = (fitted_lj_f(r_ij(:, particle1, particle2, 1), l, ff_parameters(:,atom_type)) / 2)
      f_ij(:, particle2, particle1, 1) = - f_ij(:, particle1, particle2, 1)
    end if

    ! Calculate Virial Xi_ext doublesum, non ideal part of pressure
    Xi = Xi + 2 * sum(abs(r_ij(:, particle1, particle2, 1)) * abs(f_ij(:, particle1, particle2, 1)), DIM=1)

  end do
end do

! Do Verlet integration step
  r(:,:,1) = 2 * r(:,:,0) - r(:,:,-1) + (sum(f_ij(:, :, :, 1), DIM=3) / mass(atom_type)) * (dt ** 2)

! Calculate velocities by numerical derivation
  v = (r(:,:,1) - r(:,:,-1)) / (2 * dt) ! v(t+dt)

! Move array to t = t+dt for next integration step
  r = cshift(r,  1, DIM=3)

! Thermostat
  ! Calculate E and instantaneous temperature
  E = (mass(atom_type) / 2) * sum((v(1,:) ** 2 + v(2,:) ** 2 + v(3,:) ** 2))
  T_inst = ((2 * E) / (3 * kB * nparticles))
  ! Weak coupling thermostat
  v = v * SQRT(1 + (dt / tau_T) * (T0 / T_inst - 1))
  ! Apply velocity to r(t-dt) for algorithm
  r(:,:,-1) = r(:,:,0) - (v(:,:) * dt)

! Barostat
  ! Calculate instantaneous volume and pressure
  V_inst = boxsize(1) * boxsize(2) * boxsize(3)
  P_inst = ((2 * E) / (3 * V_inst)) + (1 / (3 * V_inst)) * Xi
  ! Scaling factor
  mu = (1 - (dt / tau_P) * (P0 - P_inst)) ** (1./3.)
  ! Scale coordinates and boxsize by mu
  r = r * mu
  boxsize = boxsize * mu
  boxsize_half = boxsize_half * mu
  cutoff = cutoff * mu
  Xi = 0

! Console notifier during run
  if (mod(iter, 100)==0) then
    print *, " "
    print *, "Energy (J)",  E, "Temperature (K)", T_inst
    print *, "Volume (nm^3)",  (V_inst * 1E9 ** 3), "Pressure (Pa)", P_inst
    print *, (real(iter)/real(maxiter)*100), "% ", (iter * dt * 1E12), "ps", iter, "/", maxiter
    !print *, "Percentage non-ideal:", (((1 / (3 * V_inst)) * Xi)/(((2 * E) / (3 * V_inst)) + (1 / (3 * V_inst)) * Xi))
    print *, " "
  end if

! Save trajectories, velocities, forces and thermodynamics data every .. iterations
  ! if (mod(iter, 10)==0) then
  !   ! do particle = 1, nparticles
  !   !   write (file_name,'(i0)') particle
  !   !   call save_trajectory("export/r_"//trim(file_name)//".txt", r(:,particle,0))
  !   !   call save_trajectory("export/v_"//trim(file_name)//".txt", v(:,particle))
  !   !   call save_trajectory("export/f_"//trim(file_name)//".txt", f(:,particle,1))
  !   ! end do
  !   !call save_td_data("export/time.txt", (iter*dt))
  !   !call save_td_data("export/T.txt", T_inst)
  !   !call save_td_data("export/V.txt", V_inst)
  !   !call save_td_data("export/P.txt", P_inst)
  !   call save_trajectory("export/Ekin.txt", (0.5 * mass(atom_type) * (SQRT(v(1,:) ** 2 + v(2,:) ** 2 + v(3,:) ** 2)) ** 2))
  ! end if

! Checks if simulation failed every .. iterations
  if (mod(iter, 1000)==0) then
    if ((maxval(r(:,:,-1)) > maxval(boxsize(:))) .or. (minval(r(:,:,-1)) < 0)) then
      print *, " "
      print *, "ERROR: Calculation out of bounds, aborting after", iter, "iterations.."
      print *, " "
      exit integrator
    end if
  end if

end do integrator

! -------------------------------------------------------------------------------------------- !

call system_clock(end)
timing = real(end - beginning) / real(rate)
print *, "Elapsed time: ", timing, "seconds"
print *, "with", (maxiter/timing), "iter/s and", (t/timing*1E15), "fs/s"

end program main
