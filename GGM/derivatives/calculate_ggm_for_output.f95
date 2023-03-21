subroutine calculate_ggm_for_output(lvl,tux_lo,tux_hi,g,mfi)
!
!
!
!
!
! Called by: ggm_output_shortcut
! Calls:
! External calls:
!
USE precise
USE constants
!use iso_c_binding
use derivatives
USE ggm_stuff
use amr_info_holder
use amrex_amr_module
use amrex_base_module
implicit none
integer,intent(in) :: tux_lo(3),tux_hi(3),g,lvl
integer :: a,b,c
REAL(KIND=dp) :: x_der,y_der,z_der,wind_dir
real(kind=dp) :: max_vel,min_vel

type(amrex_mfiter) :: mfi

!write(*,*) 'entering ggm for output'

gm => mfgm(lvl)%dataptr(mfi)
ggm => mfggm(lvl)%dataptr(mfi)
u_vel => mfu_vel(lvl)%dataptr(mfi)
v_vel => mfv_vel(lvl)%dataptr(mfi)
w_vel => mfw_vel(lvl)%dataptr(mfi)
state => mfstate(lvl)%dataptr(mfi)

do c = tux_lo(3),tux_hi(3)
  do b = tux_lo(2),tux_hi(2)
    do a = tux_lo(1),tux_hi(1)

      if (state(a,b,c,1) < 0 .and. state(a,b,c,1) /= -666) then
        ggm(a,b,c,1) = 0.0D0
        cycle

      end if

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
!      if (ggm(a,b,c,g) > 1) then
!        write(*,*) 'massive ggm values',ggm(a,b,c,g),gm(a,b,c,g),x_der,y_der,z_der,a,b,c
!      end if
!write(*,*) 'after derivatives',gm(a,b,c,g),x_der,y_der,z_der,a,b,c
!write(*,*) 'after gm calc'
    end do
  end do
end do

end subroutine
