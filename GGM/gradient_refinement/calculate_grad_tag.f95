SUBROUTINE calculate_grad_tag(lvl,settag,cleartag,tux_lo,tux_hi,&
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
u_vel => mfu_vel(lvl)%dataptr(mfi)
v_vel => mfv_vel(lvl)%dataptr(mfi)
w_vel => mfw_vel(lvl)%dataptr(mfi)
state => mfstate(lvl)%dataptr(mfi)

do c = tux_lo(3),tux_hi(3)
  do b = tux_lo(2),tux_hi(2)
    do a = tux_lo(1),tux_hi(1)

      if (state(a,b,c,1) < 0) then
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

      if (gm(a,b,c,g) >= gradient_limit) then
        tagarr(a,b,c) = settag
      end if

    end do
  end do
end do

if (use_pgs) then
  call predictive_splitting(lvl,settag,cleartag,tux_lo,tux_hi,&
           tag_lo,tag_hi,tagarr,g,mfi)
end if


END SUBROUTINE
