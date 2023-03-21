SUBROUTINE calculate_gradgrad_mag(lvl,settag,cleartag,tux_lo,tux_hi,&
  tag_lo,tag_hi,tagarr,g,mfi)
!
! Calculate the gradient magnitude of the matrix gm
!
! Calls:
! Called by: ggm_calculator
!
USE precise
USE constants
use iso_c_binding
use derivatives
USE ggm_stuff
use amr_info_holder
use amrex_amr_module
use amrex_base_module
IMPLICIT NONE
INTEGER :: a,b,c
integer,intent(in) :: tux_lo(3),tux_hi(3),g,lvl,tag_lo(4),tag_hi(4)
REAL(KIND=dp) :: x_der,y_der,z_der,wind_dir
real(kind=dp) :: max_vel,min_vel
character(kind=c_char),intent(in) :: settag,cleartag
character(kind=c_char),intent(inout) :: tagarr(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),&
  tag_lo(3):tag_hi(3))

type(amrex_mfiter) :: mfi

!write(*,*) 'sparkly',lvl

gm => mfgm(lvl)%dataptr(mfi)
ggm => mfggm(lvl)%dataptr(mfi)
u_vel => mfu_vel(lvl)%dataptr(mfi)
v_vel => mfv_vel(lvl)%dataptr(mfi)
w_vel => mfw_vel(lvl)%dataptr(mfi)
state => mfstate(lvl)%dataptr(mfi)

do c = tux_lo(3),tux_hi(3)
  do b = tux_lo(2),tux_hi(2)
    do a = tux_lo(1),tux_hi(1)

      if (state(a,b,c,1) < 0) then
        ggm(a,b,c,1) = 0.0D0
        cycle

      end if
!
!
!

      x_der = derivative_windage_secondord((/gm(a,b,c,g),gm(a+1,b,c,g),gm(a+2,b,c,g),&
                gm(a-1,b,c,g),gm(a-2,b,c,g)/),(/state(a,b,c,1),state(a,b,c,2),&
                state(a,b,c,16),state(a,b,c,5),state(a,b,c,19)/),u_vel(a,b,c,1),&
                upwind,downwind,central)
      y_der = derivative_windage_secondord((/gm(a,b,c,g),gm(a,b+1,c,g),gm(a,b+2,c,g),&
                gm(a,b-1,c,g),gm(a,b-2,c,g)/),(/state(a,b,c,1),state(a,b,c,3),&
                state(a,b,c,17),state(a,b,c,6),state(a,b,c,20)/),v_vel(a,b,c,1),&
                upwind,downwind,central)
      z_der = derivative_windage_secondord((/gm(a,b,c,g),gm(a,b,c+1,g),gm(a,b,c+2,g),&
                gm(a,b,c-1,g),gm(a,b,c-2,g)/),(/state(a,b,c,1),state(a,b,c,4),&
                state(a,b,c,18),state(a,b,c,7),state(a,b,c,21)/),w_vel(a,b,c,1),&
                upwind,downwind,central)

      ggm(a,b,c,g) = sqrt(x_der**2 + y_der**2 + z_der**2)

!      if (ggm(a,b,c,g) >= ggm_limit(g,1)) then
!      if (ggm(a,b,c,g) >= ggm_sensitivity) then
      if (ggm(a,b,c,g) >= ggm_tag_limit(lvl)) then
        tagarr(a,b,c) = settag
!      else if (ggm(a,b,c,g) <= ggm_limit(g,2)) then
!        tagarr(a,b,c) = cleartag
      end if

    end do
  end do
end do

if (use_pgs) then
  call predictive_splitting(lvl,settag,cleartag,tux_lo,tux_hi,&
           tag_lo,tag_hi,tagarr,g,mfi)
end if


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
!!            in_values(3) = gm(a)
!            in_values(2) = gm(link(5,a))
!            in_values(4) = gm(link(2,a))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = gm(a)
!            in_values(2) = gm(link(5,a))
!            in_values(1) = gm(link(5,link(5,a)))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = gm(a)
!            in_values(4) = gm(link(2,a))
!            in_values(5) = gm(link(2,link(2,a)))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = gm(a)
!            in_values(2) = gm(link(5,a))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = gm(a)
!            in_values(4) = gm(link(2,a))
!            CALL x_dir_derivative(a,x_der,in_values)
!          ELSE
!            x_der = 0.0D0
!          END IF
!
!          IF (deriv_form(2,1,a)) THEN
!!            in_values(3) = gm(a)
!            in_values(2) = gm(link(6,a))
!            in_values(4) = gm(link(3,a))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = gm(a)
!            in_values(2) = gm(link(6,a))
!            in_values(1) = gm(link(6,link(6,a)))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = gm(a)
!            in_values(4) = gm(link(3,a))
!            in_values(5) = gm(link(3,link(3,a)))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = gm(a)
!            in_values(2) = gm(link(6,a))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = gm(a)
!            in_values(4) = gm(link(3,a))
!            CALL y_dir_derivative(a,y_der,in_values)
!          ELSE
!            y_der = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!!            in_values(3) = gm(a)
!            in_values(2) = gm(link(7,a))
!            in_values(4) = gm(link(4,a))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = gm(a)
!            in_values(2) = gm(link(7,a))
!            in_values(1) = gm(link(7,link(7,a)))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = gm(a)
!            in_values(4) = gm(link(4,a))
!            in_values(5) = gm(link(4,link(4,a)))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = gm(a)
!            in_values(2) = gm(link(7,a))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = gm(a)
!            in_values(4) = gm(link(4,a))
!            CALL z_dir_derivative(a,z_der,in_values)
!          ELSE
!            z_der = 0.0D0
!          END IF
!
!          ggm(a) = SQRT(x_der**2.0D0 + y_der**2.0D0 + z_der**2.0D0)
!
!        END IF
!
!      END DO
