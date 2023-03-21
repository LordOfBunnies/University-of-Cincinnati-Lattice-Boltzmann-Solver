SUBROUTINE populate_q(lvl,tux_lo,tux_hi,g,mfi)
!
! Populate the q matrix with the needed variable
!
! Calls:
! Called by: ggm_calculator
!
USE precise
use constants
USE ggm_stuff
USE freestream_values
use amrex_amr_module
use amrex_base_module
use amr_info_holder
use amr_processes
IMPLICIT NONE
INTEGER :: lvl,g,tux_lo(3),tux_hi(3),choice

type(amrex_mfiter) :: mfi

!write(*,*) 'holdin McCoque'

SELECT CASE (ggm_variables(g))
  CASE(1) !Density
    rho => mfrho(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  rho(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,1)
  CASE(2) !u-velocity
    u_vel => mfu_vel(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  u_vel(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,1)
  CASE(3) !v-velocity
    v_vel => mfv_vel(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  v_vel(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,1)
  CASE(4) !w-velocity
    w_vel => mfw_vel(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  w_vel(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts:tux_hi(3)+nghosts_mid,1)
  CASE(5) !pressure
    press => mfpress(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  press(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,1)
  CASE(6) !temperature
    temp => mftemp(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  temp(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,1)
  CASE(7) !velocity magnitude
    choice = 1
    CALL calculate_v_mag(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(8) !kinetic energy
    choice = 2
    CALL calculate_v_mag(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(9) !Mach number
    choice = 3
    CALL calculate_v_mag(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(10) !x-vorticity
    choice = 1
    CALL calculate_vort(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(11) !y-vorticity
    choice = 2
    CALL calculate_vort(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(12) !z-vorticity
    choice = 3
    CALL calculate_vort(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(13) !vorticity magnitude
    choice = 4
    CALL calculate_vort(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(14) !xy-shear stress
    choice = 1
    CALL calculate_shear_stress(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(15) !xz-shear stress
    choice = 2
    CALL calculate_shear_stress(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(16) !yz-shear stress
    choice = 3
    CALL calculate_shear_stress(lvl,tux_lo,tux_hi,g,choice,mfi)
  CASE(17) !shear stress magnitude
    choice = 4
    CALL calculate_shear_stress(lvl,tux_lo,tux_hi,g,choice,mfi)
  case(20) !alpha
    alpha => mfalpha(lvl)%dataptr(mfi)
    holder(tux_lo(1)-nghosts:tux_hi(1)+nghosts_mid,tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
      tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,g)  =  alpha(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
      tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,1)
  CASE(27)
    WRITE(*,*) 'Turbulence is not yet implemented, thank you, come again!'
    CALL error_out
  CASE(28)
    WRITE(*,*) 'Turbulence is not yet implemented, thank you, come again!'
    CALL error_out
  CASE(29)
    CALL calculate_q_criterion
  END SELECT


END SUBROUTINE
