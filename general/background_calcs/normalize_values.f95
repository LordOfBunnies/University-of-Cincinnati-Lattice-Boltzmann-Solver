SUBROUTINE normalize_values(characteristic_length)
!
! This is intended to normalize all the flow variables so they can be used
! easily later on.  Most will simply be normalized by their reference value.
! This will also normalize the inlet and outlet variables of the same types
! as though in other parts of the input file.
!
! Called by: background_calcs
! Calls:
!
USE nml_inlet_outlet
USE constants
USE nml_init
USE freestream_values
USE precise
USE nml_output
IMPLICIT NONE
INTEGER :: i
REAL(KIND=dp) :: characteristic_length

!write(*,*) 'preparing to normalize inlets and outlets'

DO i = 1,num_inlets
!
! Density and velocity inlet
!
  IF (inlet_type(1,i) == 1) THEN
    IF (inlet_flow_val(1,i) <=0.0) THEN
      inlet_flow_val(1,i) = ABS(inlet_flow_val(1,i))/rho_ref
    ELSE
      inlet_flow_val(1,i) = inlet_flow_val(1,i)/rho_ref
    END IF
    IF (inlet_flow_val(2,i) <= 0.0) THEN
      inlet_flow_val(2,i) = ABS(inlet_flow_val(2,i))/speed_of_sound
    ELSE
      inlet_flow_val(2,i) = inlet_flow_val(2,i)/speed_of_sound
    END IF
    inlet_flow_val(3,i) = inlet_flow_val(3,i)/(t_ref*gama)

    WRITE(11,251) i
    WRITE(11,201) inlet_flow_val(1,i)
    WRITE(11,205) inlet_flow_val(2,i)
!
! Pressure and temperature inlet
!
  ELSE IF (inlet_type(1,i) == 2) THEN
    IF (inlet_flow_val(1,i) <=0.0) THEN
      inlet_flow_val(1,i) = ABS(inlet_flow_val(1,i))/(p_ref*gama)
    ELSE
      inlet_flow_val(1,i) = inlet_flow_val(1,i)/(p_ref*gama)
    END IF
    IF (inlet_flow_val(2,i) <= 0.0) THEN
      inlet_flow_val(2,i) = ABS(inlet_flow_val(2,i))/(t_ref*gama)
    ELSE
      inlet_flow_val(2,i) = inlet_flow_val(2,i)/(t_ref*gama)
    END IF
    WRITE(11,251) i
    WRITE(11,204) inlet_flow_val(1,i)
    WRITE(11,206) inlet_flow_val(2,i)
!
! Freestream boundary condition
!

  ELSE IF (inlet_type(1,i) == 3) THEN
    WRITE(11,251) i
    WRITE(11,201) rho_ref
    WRITE(11,205) mach
!
! Velocity and temperature inlet
!
  ELSE IF (inlet_type(1,i) == 4) THEN
    IF (inlet_flow_val(1,i) <=0.0) THEN
      inlet_flow_val(1,i) = ABS(inlet_flow_val(1,i))/speed_of_sound
    ELSE
      inlet_flow_val(1,i) = inlet_flow_val(1,i)/speed_of_sound
    END IF
    IF (inlet_flow_val(2,i) <= 0.0) THEN
      inlet_flow_val(2,i) = ABS(inlet_flow_val(2,i))/(t_ref*gama)
    ELSE
      inlet_flow_val(2,i) = inlet_flow_val(2,i)/(t_ref*gama)
    END IF

    WRITE(11,251) i
    WRITE(11,205) inlet_flow_val(1,i)
    WRITE(11,206) inlet_flow_val(2,i)

  else if (inlet_type(1,i) == 7) then
    inlet_flow_val(1,i) = ABS(inlet_flow_val(1,i))/rho_ref
    inlet_flow_val(2,i) = ABS(inlet_flow_val(2,i))/speed_of_sound
    inlet_flow_val(3,i) = ABS(inlet_flow_val(3,i))/speed_of_sound
    inlet_flow_val(4,i) = ABS(inlet_flow_val(4,i))/speed_of_sound
    inlet_flow_val(5,i) = inlet_flow_val(5,i)/(t_ref*gama)
  END IF
END DO

DO i = 1,num_outlets
!
! Density and velocity outlet
!
  IF (outlet_type(1,i) == 1) THEN
    IF (outlet_flow_val(1,i) <=0.0) THEN
      outlet_flow_val(1,i) = ABS(outlet_flow_val(1,i))/rho_ref
    ELSE
      outlet_flow_val(1,i) = outlet_flow_val(1,i)/rho_ref
    END IF
    IF (outlet_flow_val(2,i) <= 0.0) THEN
      outlet_flow_val(2,i) = ABS(outlet_flow_val(2,i))/speed_of_sound
    ELSE
      outlet_flow_val(2,i) = outlet_flow_val(2,i)/speed_of_sound
    END IF
    outlet_flow_val(3,i) = outlet_flow_val(3,i)/(t_ref*gama)

    WRITE(11,252) i
    WRITE(11,211) outlet_flow_val(1,i)
    WRITE(11,215) outlet_flow_val(2,i)
!
! Pressure and temperature outlet
!
  ELSE IF (outlet_type(1,i) == 2) THEN
    IF (outlet_flow_val(1,i) <=0.0) THEN
      outlet_flow_val(1,i) = ABS(outlet_flow_val(1,i))/(p_ref*gama)
    ELSE
      outlet_flow_val(1,i) = outlet_flow_val(1,i)/(p_ref*gama)
    END IF
    IF (outlet_flow_val(2,i) <= 0.0) THEN
      outlet_flow_val(2,i) = ABS(outlet_flow_val(2,i))/(t_ref*gama)
    ELSE
      outlet_flow_val(2,i) = outlet_flow_val(2,i)/(t_ref*gama)
    END IF

    WRITE(11,252) i
    WRITE(11,214) outlet_flow_val(1,i)
    WRITE(11,216) outlet_flow_val(2,i)
!
! Freestream boundary condition
!
  ELSE IF (outlet_type(1,i) == 3) THEN
    WRITE(11,252) i
    WRITE(11,211) rho_ref
    WRITE(11,215) mach
!
! Velocity and temperature outlet
!
  ELSE IF (outlet_type(1,i) == 4) THEN
    IF (outlet_flow_val(1,i) <=0.0) THEN
      outlet_flow_val(1,i) = ABS(outlet_flow_val(1,i))/speed_of_sound
    ELSE
      outlet_flow_val(1,i) = outlet_flow_val(1,i)/speed_of_sound
    END IF
    IF (outlet_flow_val(2,i) <= 0.0) THEN
      outlet_flow_val(2,i) = ABS(outlet_flow_val(2,i))/(t_ref*gama)
    ELSE
      outlet_flow_val(2,i) = outlet_flow_val(2,i)/(t_ref*gama)
    END IF

    WRITE(11,251) i
    WRITE(11,215) outlet_flow_val(1,i)
    WRITE(11,216) outlet_flow_val(2,i)

  else if (outlet_type(1,i) == 7) then
!    write(*,*) 'outlet pre-normalize',outlet_flow_val(1,i),outlet_flow_val(2,i),&
!      outlet_flow_val(3,i),outlet_flow_val(4,i),outlet_flow_val(5,i)
    outlet_flow_val(1,i) = ABS(outlet_flow_val(1,i))/rho_ref
    outlet_flow_val(2,i) = ABS(outlet_flow_val(2,i))/speed_of_sound
    outlet_flow_val(3,i) = ABS(outlet_flow_val(3,i))/speed_of_sound
    outlet_flow_val(4,i) = ABS(outlet_flow_val(4,i))/speed_of_sound
    outlet_flow_val(5,i) = outlet_flow_val(5,i)/(t_ref*gama)
!    write(*,*) 'outlet post-normalize',outlet_flow_val(1,i),outlet_flow_val(2,i),&
!      outlet_flow_val(3,i),outlet_flow_val(4,i),outlet_flow_val(5,i)
  END IF



END DO
!
!
!
DO i = 1,num_init_regions
  init_data(1,i) = init_data(1,i)/(t_ref*gama)
  init_data(2,i) = init_data(2,i)/(p_ref*gama)
  init_data(3,i) = init_data(3,i)/speed_of_sound
END DO
!
!
!
WRITE(11,101) stag_density
WRITE(11,102) stag_temperature
WRITE(11,103) v_inf
WRITE(11,104) rho_inf
WRITE(11,105) temperature
WRITE(11,106) gas_const_R
WRITE(11,107) speed_of_sound
WRITE(11,108) viscosity
!      WRITE(11,109) reynolds
!
!
 101  FORMAT ('Normalized stagnation density = ', F12.7)
 102  FORMAT ('Normalized stagnation temperature = ', F12.7)
 103  FORMAT ('Normalized freestream velocity = ', F12.7)
 104  FORMAT ('Normalized freestream density = ', F12.7)
 105  FORMAT ('Normalized freestream temperature = ', F12.7)
 106  FORMAT ('Normalized gas constant = ', F12.7)
 107  FORMAT ('Normalized speed of sound = ', F12.7)
 108  FORMAT ('Normalized viscosity = ', F13.11)
! 109  FORMAT ('Reynolds number = ', F9.7)


 201  FORMAT ('Inlet density = ', F11.7)
 204  FORMAT ('Inlet pressure = ',F17.5)
 205  FORMAT ('Inlet velocity ',F9.5)
 206  FORMAT ('Inlet temperature = ',F9.4)

 211  FORMAT ('Outlet density = ', F11.7)
 214  FORMAT ('Outlet pressure = ',F17.5)
 215  FORMAT ('Outlet velocity = ',F9.5)
 216  FORMAT ('Outlet temperature = ',F9.4)

 251  FORMAT ('For inlet number ',I8)
 252  FORMAT ('For outlet number ',I8)

END SUBROUTINE


!stag_density = stag_density/rho_ref
!stag_temperature = stag_temperature/(t_ref*gama)
!stag_pressure = stag_pressure/(p_ref*gama)

!v_inf = mach
!gas_const_R = 1.0D0
!rho_inf = 1.0D0
!temperature = 1.0D0/gama
!pressure = 1.0D0/gama
!speed_of_sound = 1.0D0
!viscosity = 1.0D0

!reynolds = rho_inf*(mach*u_ref)*characteristic_length/viscosity
