SUBROUTINE calculate_v_mag(lvl,tux_lo,tux_hi,g,choice,mfi)
!
! Calculate the velocity magnitude and apply it to the q array
!
! Calls:
! Called by: populate_q
!
USE ggm_stuff
use freestream_values
USE linkwise
USE precise
USE constants
use amr_info_holder
use amrex_base_module
use amrex_amr_module
IMPLICIT NONE
INTEGER :: a,b,c
integer, intent(in) :: tux_lo(3),tux_hi(3),lvl,choice,g
!real(kind=dp),intent(out) ::  holder(tux_lo(1)-nghosts:tux_hi(1)+nghosts,&
!  tux_lo(2)-nghosts:tux_hi(2)+nghosts,tux_lo(3)-nghosts:tux_hi(3)+nghosts,num_ggm_vars)
real(kind=dp) :: temp_rat
type(amrex_mfiter) :: mfi
!
! Assign the proper values to q before it's evaulated later.
!
! 2D first
!
IF (dimensions == 2) THEN


!
!CALL amrex_mfiter_build(mfi,mfu_vel(level),tiling=.false.)

ELSE IF (dimensions ==3) THEN
  u_vel => mfu_vel(lvl)%dataptr(mfi)
  v_vel => mfv_vel(lvl)%dataptr(mfi)
  w_vel => mfw_vel(lvl)%dataptr(mfi)
  rho => mfrho(lvl)%dataptr(mfi)
  temp => mftemp(lvl)%dataptr(mfi)

  do c = tux_lo(3)-nghosts_mid,tux_hi(3)+nghosts_mid
    do b = tux_lo(2)-nghosts_mid,tux_hi(2)+nghosts_mid
      do a = tux_lo(1)-nghosts_mid,tux_hi(1)+nghosts_mid
        IF (choice == 1) THEN
          holder(a,b,c,g) = SQRT(u_vel(a,b,c,1)**2+v_vel(a,b,c,1)**2+w_vel(a,b,c,1)**2)
        ELSE IF (choice == 2) THEN
          holder(a,b,c,g) = 0.5D0*rho(a,b,c,1)*SQRT(u_vel(a,b,c,1)**2 + v_vel(a,b,c,1)**2 + w_vel(a,b,c,1)**2)
        ELSE IF (choice == 3) THEN
          temp_rat = temp(a,b,c,1)/temperature
          holder(a,b,c,g) = SQRT(u_vel(a,b,c,1)**2+v_vel(a,b,c,1)**2+w_vel(a,b,c,1)**2)/temp_rat
        END IF
      end do
    end do
  end do
END IF


END SUBROUTINE
!          IF (bounding_method == -1) THEN
!  IF (choice == 1) THEN
!    q(a) = SQRT(u_vel(a)**2.0D0+v_vel(a)**2.0D0)
!  ELSE IF (choice == 2) THEN
!    q(a) = 0.5D0*rho(a)*SQRT(u_vel(a)**2.0D0 + v_vel(a)**2.0D0)
!  ELSE IF (choice == 3) THEN
!    q(a) = SQRT(u_vel(a)**2.0D0+v_vel(a)**2.0D0)
!  END IF


!real(kind=dp),intent(in) ::  rho(tux_lo(1)-nghosts:tux_hi(1)+nghosts,&
!  tux_lo(2)-nghosts:tux_hi(2)+nghosts,tux_lo(3)-nghosts:tux_hi(3)+nghosts,1),&
!  u_vel(tux_lo(1)-nghosts:tux_hi(1)+nghosts,&
!  tux_lo(2)-nghosts:tux_hi(2)+nghosts,tux_lo(3)-nghosts:tux_hi(3)+nghosts,1),&
!  v_vel(tux_lo(1)-nghosts:tux_hi(1)+nghosts,&
!  tux_lo(2)-nghosts:tux_hi(2)+nghosts,tux_lo(3)-nghosts:tux_hi(3)+nghosts,1),&
!  w_vel(tux_lo(1)-nghosts:tux_hi(1)+nghosts,&
!  tux_lo(2)-nghosts:tux_hi(2)+nghosts,tux_lo(3)-nghosts:tux_hi(3)+nghosts,1)
!      REAL(KIND=dp) :: vmag

!          ELSE IF (bounding_method == -2) THEN
!            IF (choice == 1) THEN
!              q(a) = SQRT(u_vel(a)**2.0D0+w_vel(a)**2.0D0)
!            ELSE IF (choice == 2) THEN
!              q(a) = 0.5D0*rho(a)*SQRT(u_vel(a)**2.0D0 + w_vel(a)**2.0D0)
!            ELSE IF (choice == 3) THEN
!              q(a) = SQRT(u_vel(a)**2.0D0+w_vel(a)**2.0D0)
!            END IF
!          ELSE IF (bounding_method == -1) THEN
!            IF (choice == 1) THEN
!              q(a) = SQRT(v_vel(a)**2.0D0+w_vel(a)**2.0D0)
!            ELSE IF (choice == 2) THEN
!              q(a) = 0.5D0*rho(a)*SQRT(v_vel(a)**2.0D0 + w_vel(a)**2.0D0)
!            ELSE IF (choice == 3) THEN
!              q(a) = SQRT(v_vel(a)**2.0D0+w_vel(a)**2.0D0)
!            END IF
!          END IF
!
! Three dimensional is actually easier to do.

