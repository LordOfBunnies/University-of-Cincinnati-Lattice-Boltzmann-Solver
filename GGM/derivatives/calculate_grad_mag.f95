SUBROUTINE calculate_grad_mag(lvl,tux_lo,tux_hi,g,mfi)
!
! Calculate the gradient magnitude of the matrix q
!
! Calls:
! Called by: ggm_calculator
!
USE precise
USE constants
use derivatives
USE ggm_stuff
use amr_info_holder
use amrex_amr_module
use amrex_base_module
IMPLICIT NONE
INTEGER :: a,b,c
integer,intent(in) :: tux_lo(3),tux_hi(3),g,lvl
integer :: reach
REAL(KIND=dp) :: x_der,y_der,z_der,max_vel,min_vel,wind_dir

type(amrex_mfiter) :: mfi

reach = nghosts-4

!write(*,*) 'giggle buster',lvl

gm => mfgm(lvl)%dataptr(mfi)
u_vel => mfu_vel(lvl)%dataptr(mfi)
v_vel => mfv_vel(lvl)%dataptr(mfi)
w_vel => mfw_vel(lvl)%dataptr(mfi)
state => mfstate(lvl)%dataptr(mfi)

do c = tux_lo(3)-reach,tux_hi(3)+reach
  do b = tux_lo(2)-reach,tux_hi(2)+reach
    do a = tux_lo(1)-reach,tux_hi(1)+reach

      if (state(a,b,c,1) < 0 .and. state(a,b,c,1) /= -666) then
        gm(a,b,c,1) = 0.0D0
        cycle
      end if

      x_der = derivative_windage_secondord((/holder(a,b,c,g),holder(a+1,b,c,g),holder(a+2,b,c,g),&
                holder(a-1,b,c,g),holder(a-2,b,c,g)/),(/state(a,b,c,1),state(a,b,c,2),&
                state(a,b,c,16),state(a,b,c,5),state(a,b,c,19)/),u_vel(a,b,c,1),&
                upwind,downwind,central)
!      write(*,*) 'bubble',a,b,c
      y_der = derivative_windage_secondord((/holder(a,b,c,g),holder(a,b+1,c,g),holder(a,b+2,c,g),&
                holder(a,b-1,c,g),holder(a,b-2,c,g)/),(/state(a,b,c,1),state(a,b,c,3),&
                state(a,b,c,17),state(a,b,c,6),state(a,b,c,20)/),v_vel(a,b,c,1),&
                upwind,downwind,central)
!      write(*,*) 'sparkle',a,b,c
      z_der = derivative_windage_secondord((/holder(a,b,c,g),holder(a,b,c+1,g),holder(a,b,c+2,g),&
                holder(a,b,c-1,g),holder(a,b,c-2,g)/),(/state(a,b,c,1),state(a,b,c,4),&
                state(a,b,c,18),state(a,b,c,7),state(a,b,c,21)/),w_vel(a,b,c,1),&
                upwind,downwind,central)
!write(*,*) 'after derivatives',x_der,y_der,z_der
      gm(a,b,c,g) = sqrt(x_der**2 + y_der**2 + z_der**2)
!write(*,*) 'after gm calc'
    end do
  end do
end do


END SUBROUTINE
!      DO a = 1,link_number
!        IF (dimensions == 2) THEN
!          IF (bounding_method == -1) THEN
!
!          ELSE IF (bounding_method == -2) THEN
!
!          ELSE IF (bounding_method == -3) THEN
!
!
!          END IF
!
!        ELSE IF (dimensions == 3) THEN
!
!          IF (deriv_form(1,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = q(link(5,a))
!            in_values(4) = q(link(2,a))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = q(a)
!            in_values(2) = q(link(5,a))
!            in_values(1) = q(link(5,link(5,a)))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = q(a)
!            in_values(4) = q(link(2,a))
!            in_values(5) = q(link(2,link(2,a)))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = q(a)
!            in_values(2) = q(link(5,a))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = q(a)
!            in_values(4) = q(link(2,a))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE
!            x_der = 0.0D0
!          END IF
!
!          IF (deriv_form(2,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = q(link(6,a))
!            in_values(4) = q(link(3,a))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = q(a)
!            in_values(2) = q(link(6,a))
!            in_values(1) = q(link(6,link(6,a)))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = q(a)
!            in_values(4) = q(link(3,a))
!            in_values(5) = q(link(3,link(3,a)))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = q(a)
!            in_values(2) = q(link(6,a))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = q(a)
!            in_values(4) = q(link(3,a))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE
!            y_der = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = q(link(7,a))
!            in_values(4) = q(link(4,a))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = q(a)
!            in_values(2) = q(link(7,a))
!            in_values(1) = q(link(7,link(7,a)))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = q(a)
!            in_values(4) = q(link(4,a))
!            in_values(5) = q(link(4,link(4,a)))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = q(a)
!            in_values(2) = q(link(7,a))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = q(a)
!            in_values(4) = q(link(4,a))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE
!            z_der = 0.0D0
!          END IF
!
!          gm(a) = SQRT(x_der**2.0D0 + y_der**2.0D0 + z_der**2.0D0)
!
!        END IF
!
!      END DO
