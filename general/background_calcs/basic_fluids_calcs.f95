SUBROUTINE basic_fluid_calcs
!
! This is to do basic fluid calculations.
!
! Called by: background_calcs
! Calls: sutherlands_calc
!
USE constants
USE freestream_values
USE timing
USE precise
use quick_calcs
IMPLICIT NONE
REAL(KIND=dp) :: dummy_kappa
!      REAL(KIND=dp) :: characteristic_length
! Define the gas constant first
gas_const_R = 287.052874D0
!
! Basic fluid calcs for use elsewhere
!
rho_ref = pressure/(gas_const_R*temperature)
! Define a few reference quantities
t_ref = temperature
p_ref = pressure
non_dim_R = 1.0D0
const_vol = 1.0D0/(gama-1.0D0)
const_press = gama/(gama-1.0D0)
cv_unit_air = const_vol*gas_const_R
cp_unit_air = const_press*gas_const_R

temperature = 1.0D0/gama
pressure = 1.0D0/gama
rho_inf = 1.0D0


!      speed_of_sound = SQRT(gama*gas_const_R*temperature)
!      iso_therm_sos = SQRT(const_R*temperature)
speed_of_sound = SQRT(gama*gas_const_R*t_ref)
u_ref = speed_of_sound
!
! Basic LBM simulations are isothermal, this changes the speed of sound so
! this is something to return to later. For now, the velocities have been low
! because of this fact. XXXX
!
write(*,*) v_inf
mach = v_inf/u_ref
!      mach = v_inf/iso_therm_sos
v_inf = mach
if (primary_flow_dir == 1) then
  fs_u = v_inf
  fs_v = 0.0D0
  fs_w = 0.0D0
else if (primary_flow_dir == 2) then
  fs_u = 0.0D0
  fs_v = v_inf
  fs_w = 0.0D0
else if (primary_flow_dir == 3) then
  fs_u = 0.0D0
  fs_v = 0.0D0
  fs_w = v_inf
else if (primary_flow_dir == 4) then
  fs_u = -v_inf
  fs_v = 0.0D0
  fs_w = 0.0D0
else if (primary_flow_dir == 5) then
  fs_u = 0.0D0
  fs_v = -v_inf
  fs_w = 0.0D0
else if (primary_flow_dir == 6) then
  fs_u = 0.0D0
  fs_v = 0.0D0
  fs_w = -v_inf
end if
!
! Stagnation calculations
!
stag_temperature = (1.0D0 + ((gama-1.0D0)/2.0D0)*&
  mach**2)*temperature
stag_density = ((1.0D0 + ((gama-1.0D0)/2.0D0)*&
  mach**2)**(1.0D0/(gama-1.0D0)))*rho_inf
stag_pressure = ((1.0D0 + ((gama-1.0D0)/2.0D0)*&
  mach**2)**(gama/(gama-1.0D0)))*pressure
!
! Find the mass of a single molecule of air (yes, I know, don't give me that look)
!
molec_weight = molec_mass_air/(avocados_number*1000)
!
! Viscosity calculation with a call to a Sutherland's Law subroutine
!
!      CALL sutherlands_calc
! Reynolds number calculation
!      reynolds = rho_inf*v_inf*characteristic_length/viscosity
!
WRITE(11,*) 'Freestream values'
WRITE(11,102) rho_inf
WRITE(11,103) speed_of_sound
WRITE(11,104) mach
WRITE(11,105) stag_temperature
WRITE(11,106) stag_density
WRITE(11,107) stag_pressure
!      WRITE(11,108) reynolds
!      WRITE(11,109) viscosity
WRITE(11,*) ''

! 101  FORMAT (F9.7)
 102  FORMAT ('Freestream density = ',F11.5)
 103  FORMAT ('Speed of sound = ',F11.5)
 104  FORMAT ('Mach number = ',F11.5)
 105  FORMAT ('Stagnation temperature = ',F11.5)
 106  FORMAT ('Stagnation density = ',F11.5)
 107  FORMAT ('Stagnation pressure = ',F13.5)
! 108  FORMAT ('Reynolds number = ',F13.4)
! 109  FORMAT ('Viscosity = ',F11.9)


END SUBROUTINE
