subroutine timing_calcs
!
! Basic calculations for number of current timesteps
!
! Called by: Background_calcs
! Calls:
!
use amrex_base_module
use amrex_amr_module
use amr_processes, only : amr_max_lvl
use amr_info_holder, only : timestep,max_timesteps,dt,elbm_epsilon,&
  time,time_prev,init_epsilon,startup_iterations,amr,regrid_times,regrid_timesteps
use precise
use constants
use freestream_values
use grid_data
use timing
use ggm_stuff
implicit none
integer :: lvl,i,num

!write(*,*) 'Starting timing calculations.'
startup_iterations = 100
dt = amrex_geom%dx(1)/speed_of_sound
!elbm_epsilon = (dt/characteristic_time)!*100D0
elbm_epsilon = (dt/1.0D0)
!init_epsilon = (1.0D0/real(startup_iterations,kind=dp))*dt/characteristic_time
init_epsilon = elbm_epsilon/startup_iterations!*10
timestep = 0
time = 0.0D0
time_prev = -dt
do lvl = 0,amrex_max_level
  max_timesteps(lvl) = CEILING(total_time/dt(lvl))
end do
max_timesteps = max_timesteps+1

do lvl = 0,amrex_max_level
  WRITE(11,101) dt(lvl),lvl
  write(11,100) max_timesteps(lvl),lvl
  write(11,102) elbm_epsilon(lvl),lvl
  write(11,103) rdx(lvl),lvl
!  WRITE(*,101) dt(lvl),lvl
!  write(*,100) max_timesteps(lvl),lvl
end do
!
!
! XXXX this will need to be updated if coarse time functions are eventually added
!
!
if (coarse_time > very_small_number .or. amr) then
  if (coarse_time > very_small_number) then
    num_regrids = floor((total_time-coarse_time)/ggm_dt)+1
  else
    num_regrids = floor((total_time-coarse_time)/ggm_dt)
  end if


  allocate(regrid_times(num_regrids))
  allocate(regrid_timesteps(num_regrids,0:amrex_max_level))

  do i = 1,num_regrids
    if (coarse_time > very_small_number) then
      regrid_times(i) = coarse_time + ggm_dt*(i-1)
    else
      regrid_times(i) = ggm_dt*i
    end if


  !  do num = 1,num_ggm_vars
      do lvl = 0,amrex_max_level
        regrid_timesteps(i,lvl) = floor(regrid_times(i)/dt(lvl))

        if (regrid_timesteps(i,lvl) < 1) then
          regrid_timesteps(i,lvl) = -1
        end if

  !    end do
    end do

  end do
!  write(*,*) 'regrid times',regrid_times
!  write(*,*) 'regrid timesteps',regrid_timesteps(:,1)

  regrid_basis = 1

  if (allocated(ggm_limit)) deallocate(ggm_limit)
  allocate(ggm_limit(num_ggm_vars,2))
  if (allocated(ggm_tag_limit)) deallocate(ggm_tag_limit)
  allocate(ggm_tag_limit(0:amr_max_lvl))

  ggm_limit(1:num_ggm_vars,1) = 0.04D0
  ggm_limit(1:num_ggm_vars,2) = -0.01D0

  ggm_tag_limit(0) = ggm_sensitivity
  do lvl = 1,amr_max_lvl
    ggm_tag_limit(lvl) = 0.9D0*ggm_tag_limit(lvl-1)

  end do

end if

!
!
!
!
!
 100  FORMAT ('Total number of timesteps  = ',I12,'at level ',I2)
 101  FORMAT ('Timestep size = ', F9.7,' at level ',I2)
 102  format ('Epsilon values for ELBM = ', F11.8,' at level ',I2)
 103 format ('Node length ', F11.5,' at level ',I2)
end subroutine

!  write(*,*) 'regrid times',regrid_times(i),regrid_timesteps(i,0:amrex_max_level)
