SUBROUTINE viscosity_and_timing(characteristic_length)
!
! Basic Sutherland's Law calculations, implemented for air only
! at this point.  Later it'll be for more
!
! Called by: background_calcs
! Calls:
!
!USE grid_data, only: rdx
USE freestream_values
USE constants
USE precise
use quick_calcs
USE timing
use amrex_base_module
use amrex_amr_module
!use amr_processes, only: beta_1_calculator,beta_2_calculator
!use amr_info_holder, only : dt,timestep,max_timesteps
IMPLICIT NONE
REAL(KIND=dp) :: sutherlands_const,beta_1_diff,beta_2_diff
REAL(KIND=dp) :: visc_ref,temp_ref
real(kind=dp) :: characteristic_velocity,kappa_base
REAL(KIND=dp) :: characteristic_length,beta_1,beta_2
logical :: startup
!integer :: lvl

!write(*,*) 'viscosity assignment'

sutherlands_const = 110.4
visc_ref = 1.716E-5
temp_ref = 273.15

viscosity = visc_ref*((t_ref/temp_ref)**1.5)*((temp_ref+&
         sutherlands_const)/(t_ref+sutherlands_const))
nu = viscosity/rho_ref
!mu_ref = viscosity!rho_ref*speed_of_sound*default_char_leng
mu_ref = rho_ref*u_ref*mach*1!characteristic_length
!write(*,*) 'references for visc',mu_ref,viscosity,rho_ref,u_ref*mach,characteristic_length
!mu_inf = 1.0D0

characteristic_velocity = v_inf*u_ref
!
! Get the reference thermal conductivity
!
!kappa_ref = gas_const_R*viscosity
startup = .true.

call quick_therm_cond(kappa_base,rho_inf,temperature,startup)
write(*,*) kappa_ref,rho_ref,gas_const_R,viscosity
!kappa_ref = rho_inf*gas_const_R*viscosity

kappa_ref = rho_inf*speed_of_sound**3/temperature

if (characteristic_velocity < 1.0D0) then
  characteristic_time = 1.0D0
else
!  characteristic_time = characteristic_length/characteristic_velocity
  characteristic_time = 1.0D0!/(characteristic_velocity)!speed_of_sound
end if

init_beta_1 = 0.5001D0!exp(-mach/6.0D0)-0.01D0
init_beta_2 = 0.5001D0!exp(-mach/6.0D0)-0.025D0

beta_1 = 1.0D0/((2.0D0*viscosity)/(rho_inf*temperature)+1.0D0)
beta_2 = 1.0D0/((2.0D0*kappa_base/kappa_ref)/(rho_inf*const_press*temperature)+1.0D0)

beta_1_scaling = (beta_1 - init_beta_1)/real(startup_timesteps,kind=dp)
beta_2_scaling = (beta_2 - init_beta_2)/real(startup_timesteps,kind=dp)

!write(*,*)
!write(*,*) 'beta values',init_beta_1,init_beta_2,beta_1,beta_2,beta_1_scaling,beta_2_scaling
!write(*,*) 'beta pre-calcs values',kappa_base,kappa_ref,viscosity,mu_ref,rho_ref,temperature
reynolds = rho_ref*u_ref*characteristic_length/mu_ref

write(11,108) reynolds
write(11,109) viscosity
WRITE(11,111) characteristic_time
write(11,112) mu_ref
write(11,113) kappa_ref


 108  FORMAT ('Reynolds number = ',F17.4)
! 111  FORMAT ('Lattice Reynolds number = ',F13.4)
 109  FORMAT ('Viscosity = ',F13.11)
 110  FORMAT ('Kinematic Viscosity =', F15.11)
 111  format ('Characteristic time =',F9.7)
 112  format ('Refernce Viscosity, mu = ',F13.11)
 113  format ('Reference Thermal Conductivity, kappa = ',F13.11)


END SUBROUTINE

!kappa_ref = 1.0D0
!CALL quick_therm_cond(dummy_kappa,rho_inf,temperature)

!WRITE(*,108) reynolds
!WRITE(*,111) lattice_reynolds
!WRITE(*,109) viscosity
!WRITE(*,110) nu

!WRITE(11,301) tau
!WRITE(11,100) timestep
!write(11,305)

! 200  FORMAT ('Value of tau is  ', F9.6)
! 301  FORMAT ('Tau = ',F9.7)
! 305  format ('Total number of timesteps at finest level = ',I12)

!do i = 1,amrex_max_level
!  dt(i) = dt(i-1)/2.0D0
!end do
!!      tau = nu/(cs**2.0D0*dt) + 0.5D0
!      tau = nu*dt/(dx**2.0D0 * cs**2.0D0) +.5
!
!      IF (tau < 0.50001) THEN
!
!        CALL tau_adjuster(nu)
!
!      ELSE IF (tau > 20.0) THEN
!
!        CALL tau_adjuster(nu)
!
!      ELSE
!        WRITE(11,200) tau
!      END IF

!max_ts = CEILING(total_time/dt(amrex_max_level))
