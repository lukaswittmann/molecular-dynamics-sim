module energy_minimization
  
    use types, only: dp
    use forcefields
    use utils

    implicit none

    contains

function minimize_energy(r_in, nparticles, properties, boxsize, l0, cutoff)
    integer, intent(in):: nparticles
    real(dp), dimension(3, nparticles, -1:2), intent(in) :: r_in
    real(dp), dimension(3, nparticles, -1:2) :: minimize_energy
    real(dp), dimension(3, nparticles) :: v    
    real(dp), intent(in) :: properties(:,:,:)
    real(dp), dimension(3), intent(in) :: boxsize 
    real(dp), intent(in) :: l0, cutoff

    real(dp), parameter :: dt = 0.01
    integer :: iter = 0, particle1, particle2, particle, dim, dim1, dim2, dim3
    character (len=10) :: file_name
    
    real(dp), dimension(3, nparticles,1) :: f      ! ([f(x),f(y),f(z)], particle)
    real(dp), dimension(3) :: d, d0
    real(dp) :: l ! Needed for cutoff and calculation of force

    minimize_energy = r_in ! Its the positional array
    !call print_rv(minimize_energy(:,:,0), v, nparticles)

! -------------------------------------------------------------------------------------------- !
    
! Integration loop
integrator: do
  iter = iter +1



  f = 0
  do particle1 = 1, nparticles   ! Interaction for each particle
    do particle2 = 1, nparticles ! ..with each other particle
      if (particle1 /= particle2) then  

      d0 = minimize_energy(:, particle2, 0) - boxsize(:) - minimize_energy(:, particle1, 0)

        do dim1 = 1,3
          do dim2 = 1,3
            do dim3 = 1,3
  
              d = d0
              d(1) = d(1) + (dim1 - 1) * boxsize(1)
              d(2) = d(2) + (dim2 - 1) * boxsize(2)
              d(3) = d(3) + (dim3 - 1) * boxsize(3)
  
              l = SQRT(d(1) ** 2 + d(2) ** 2 + d(3) ** 2)
  
              if (l < cutoff) then
                f(:,particle1,1) = f(:,particle1,1) + min_f(d, l, l0)
              end if
  
            end do
          end do
        end do
  
      end if 
    end do
  end do
        
  ! Do Verlet integration step using
  minimize_energy(:,:,1) = 2 * minimize_energy(:,:,0) - minimize_energy(:,:,-1) + (f(:,:,1) / properties(:,:,1)) * (dt ** 2)
    
  ! Calculate velocities by numerical derivation
  v = (minimize_energy(:,:,0) - minimize_energy(:,:,-1)) / (1 * dt)

  ! Damping and implementation of velocity
  minimize_energy(:,:,0) = minimize_energy(:,:,1) - (v(:,:)/2 * dt)
    
  ! Move array to t = t+dt for next integration step
  minimize_energy = cshift(minimize_energy,  1, DIM=3)

  ! "Particle-teleporter" to fulfill periodic boundary condition
  do particle = 1, nparticles ! Iterate over all particles
      do dim = 1, 3 ! Iterate over 3d space
      
        if (minimize_energy(dim, particle, 0) > boxsize(dim)) then
            minimize_energy(dim, particle, :) = minimize_energy(dim, particle, :) - boxsize(dim)
        else if (minimize_energy(dim, particle, 0) < 0) then
            minimize_energy(dim, particle, :) = minimize_energy(dim, particle, :) + boxsize(dim)
        end if
      
      end do
    end do
    
  if (mod(iter, 500)==0) then
    ! Calculate temperature (kinetic energy) in [a.u]
    print *, "Energy =",  (properties(1,1,1) / 2 * sum(v ** 2)), "at step =", iter
    ! Quit energy minimization if done
    if ((properties(1,1,1) / 2 * sum(v ** 2)) < 1E-4) EXIT integrator
  end if
       
    ! Save trajectories, velocities and forces every .. iterations
    !if (mod(iter, 1)==0) then
    !  do particle = 1, nparticles
    !    write (file_name,'(i0)') particle
    !    call save_trajectory("export/r_"//trim(file_name)//".txt", minimize_energy(:,particle,0))
    !    call save_trajectory("export/v_"//trim(file_name)//".txt", v(:,particle))
    !    call save_trajectory("export/f_"//trim(file_name)//".txt", f(:,particle,1))
    !  end do
    !end if

  ! Hardexit after .. iterations
  if (iter == 20000) then
      print *, "Minimization aborted, maximum iterations reached!"
      exit integrator
  end if
   
end do integrator

! -------------------------------------------------------------------------------------------- !

print '(/,A,/)', "Minimization done!"

end function minimize_energy


end module energy_minimization
    