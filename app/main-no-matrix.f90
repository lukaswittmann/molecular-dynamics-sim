program main
  
use types
use forcefields
use utils
use energy_minimization
use constants
implicit none

! Md run settings
integer, parameter :: nparticles = 250
real(dp), parameter :: dt = (0.002 * 1E-12), t = (200. * 1E-12) ! dt in ps to s, t in ps to s
real(dp), parameter :: cutoff = (5. * 1E-9) ! in nm to m
integer, parameter :: maxiter = INT(t/dt)
real(dp), dimension(3) :: boxsize = (/10.,10.,10./) * 1E-9 ! Boxsize in nm to meter

! Particle settings
integer, parameter :: atom_type = 3 ! 1: He, 2: Ne, 3: Ar, 4: Kr. 5: Xe
real(dp), parameter :: sigma = 3.40100e-01 , epsilon = 9.78638e-01 ! Ar

! Thermostat settings
real(dp), parameter :: tau_T = (1 * 1E-12) ! relaxation time of thermostat in ps to s
real(dp) :: E, T_inst, T0 = 300. ! in K

! Barostat settings
real(dp), parameter :: tau_P =  (1. * 1E-6) ! relaxation time of thermostat in ns to s
real(dp) :: mu, V_inst, P_inst, P0 = (1. * 1E5)  ! in Bar to Pa

! Runtime variables
integer :: beginning, rate, end, iter, particle1, particle2, interaction, particle, dim, dim1, dim2, dim3
character (len=10) :: file_name
real(dp), dimension(3, nparticles, 1) :: properties ! 1: mass  2: ..   3: .. [x,y,z] for anisotropy
real(dp), dimension(3, nparticles, -1:2) :: r       ! ([x,y,z], particle, [t-dt, t, t+dt, t+2dt])
real(dp), dimension(3, nparticles) :: v = 0         ! ([x,y,z], particle)
real(dp), dimension(3, nparticles,1) :: f = 0       ! ([f(x),f(y),f(z)], particle)
real(dp), dimension(3) :: d, d0
real(dp) :: l, test_R = NA * kB

! Generate random starting positions in boxsize and initial velocities
r(:,:,:) = rand_r(nparticles, boxsize) 
v(:,:) = rand_v(nparticles, int(1 * 1E5))
r(:,:,-1) = r(:,:,0) - v(:,:) * dt
!call print_rv(r(:,:,0:0), v, nparticles) 

! Energy minimization
!print '(/,A,/)', "** ENERGY MINIMIZATION **"
!r = minimize_energy(r, nparticles, properties, boxsize, 5E-11_dp, cutoff)

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


  ! Redo interaction: 
  ! Implement 3D nparticle x nparticle rij and fij
! cdef void pbc(real* v1, real* v2, real* box, real* boxh, real* v1v2) nogil:
!     """
!     Calculates the vector from v1 to v2 with periodic boundary conditions.
!     """
!     cdef int i
!     v1v2[0] = v2[0] - v1[0]
!     v1v2[1] = v2[1] - v1[1]
!     v1v2[2] = v2[2] - v1[2]

!     if v1v2[0] < -1 * boxh[0]:
!             v1v2[0] += box[0]
!     elif v1v2[0] > boxh[0]:
!             v1v2[0] -= box[0]

!     if v1v2[1] < -1*boxh[1]:
!             v1v2[1] += box[1]
!     elif v1v2[1] > boxh[1]:
!             v1v2[1] -= box[1]

!     if v1v2[2] < -1 * boxh[2]:
!             v1v2[2] += box[2]
!     elif v1v2[2] > boxh[2]:
!             v1v2[2] -= box[2]
!     return


  ! Calculation of particle interactions
  f = 0
  do particle1 = 1, nparticles   ! Interaction for each particle
    ! This will be parallelized
    do particle2 = 1, nparticles ! ..with each other particle
      if (particle1 /= particle2) then  

        d0 = r(:, particle2, 0) - r(:, particle1, 0)
        do dim1 = -1,1
          do dim2 = -1,1
            do dim3 = -1,1

              ! Define new position of image and for dim1=dim2=dim3=1 of real particle interaction
              d = d0
              ! abstandsmatrix implementieren (drekecksmatrix) !
              d(1) = d(1) + (dim1) * boxsize(1)
              d(2) = d(2) + (dim2) * boxsize(2)
              d(3) = d(3) + (dim3) * boxsize(3)

              l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)

              if (l < cutoff) then
                f(:,particle1,1) = f(:,particle1,1) + fitted_lj_f(d, l, ff_parameters(:,atom_type))
              end if

            end do
          end do
        end do

      end if
    end do

  end do

! Do Verlet integration step
  r(:,:,1) = 2 * r(:,:,0) - r(:,:,-1) + (f(:,:,1) / mass(atom_type)) * (dt ** 2)

! Calculate velocities by numerical derivation
  v = (r(:,:,1) - r(:,:,-1)) / (2 * dt) ! v(t+dt)

! Move array to t = t+dt for next integration step
  r = cshift(r,  1, DIM=3)

! Thermostat
  ! Calculate E and instantaneous temperature
  E = (mass(atom_type) / 2) * sum(v ** 2)
  T_inst = ((2 * E) / (3 * kB * nparticles))
  ! Weak coupling thermostat
  v = v * SQRT(1 + (dt / tau_T) * (T0 / T_inst - 1))
  ! Apply velocity to r(t-dt) for algorithm
  r(:,:,-1) = r(:,:,0) - (v(:,:) * dt)

! Barostat
  ! Calculate instantaneous volume and pressure
  V_inst = boxsize(1) * boxsize(2) * boxsize(3)
  !!! WRONG: CALCULATION OF INST PRESSURE !!!
  ! Doppelsumme des Virial muss Abstand der Teilchen mal absolute Kraft zwischen diesen Teilchen
  P_inst = ((2 * E) / (3 * V_inst)) + (1 / (3 * V_inst)) * (sum(r(:,:,-1) * abs(f(:,:,1))))
  ! Scaling factor
  mu = (1 - (dt / tau_P) * (P0 - P_inst)) ** (1./3.)
  ! Scale coordinates and boxsize by mu
  r = r * mu
  boxsize = boxsize * mu

! Console notifier during run
  if (mod(iter, 100)==0) then
    print *, " "
    print *, "Energy (J)",  E, "Temperature (K)", T_inst
    print *, "Volume (nm^3)",  (V_inst * 1E9 ** 3), "Pressure (Pa)", P_inst
    print *, (real(iter)/real(maxiter)*100), "% ", (iter * dt * 1E12), "ps", iter, "/", maxiter
    print *, " "
  end if

! Save trajectories, velocities, forces and thermodynamics data every .. iterations
  ! if (mod(iter, 1)==0) then
  !   ! do particle = 1, nparticles
  !   !   write (file_name,'(i0)') particle
  !   !   call save_trajectory("export/r_"//trim(file_name)//".txt", r(:,particle,0))
  !   !   call save_trajectory("export/v_"//trim(file_name)//".txt", v(:,particle))
  !   !   call save_trajectory("export/f_"//trim(file_name)//".txt", f(:,particle,1))
  !   ! end do
  !   call save_td_data("export/time.txt", (iter*dt))
  !   call save_td_data("export/T.txt", T_inst)
  !   call save_td_data("export/V.txt", V_inst)
  !   call save_td_data("export/P.txt", P_inst)
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

if (iter == 10) exit integrator

end do integrator

! -------------------------------------------------------------------------------------------- !
  
!call print_rv(r, v, nparticles) 

print *, "Mean velocity", sum(abs(v)) / nparticles

call system_clock(end)
print *, "Elapsed time: ", real(end - beginning) / real(rate), "seconds."

end program main
