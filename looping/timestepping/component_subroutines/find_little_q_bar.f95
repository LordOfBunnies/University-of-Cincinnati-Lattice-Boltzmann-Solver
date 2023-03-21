subroutine find_little_q_bar(loc_gi_diff,loc_u,loc_v,&
  loc_w,loc_temp,contr_heat_flux,lvl)
!
! Find the small q bar for the Pr < 1 (air)
!
!
! Called by: collision
! Calls:
!
use precise
use constants
use linkwise
use grid_data
use amr_info_holder, only: dt,timestep
use amr_processes, only: epsilon_calculator,amr_max_lvl
implicit none
real(kind=dp),intent(in) :: loc_gi_diff(dir),loc_u,loc_v,loc_w,loc_temp
real(kind=dp) ::little_q_bar(3),u_term,v_term,w_term, x_value,y_value,z_value,&
  q_lead,loc_epsilon
real(kind=dp),intent(out) :: contr_heat_flux(dir)
integer :: i,lvl!a,b,c,

little_q_bar = 0.0D0
contr_heat_flux = 0.0D0

loc_epsilon = epsilon_calculator(loc_u,loc_v,loc_w,loc_temp,&
  timestep(amr_max_lvl),-1,-2,-3,lvl)

q_lead = loc_epsilon/loc_temp

do i = 1,dir
  little_q_bar(1) = little_q_bar(1) + (rcx(i)-loc_u)*(loc_gi_diff(i))
  little_q_bar(2) = little_q_bar(2) + (rcy(i)-loc_v)*(loc_gi_diff(i))
  if (dimensions ==3) then
    little_q_bar(3) = little_q_bar(3) + (rcz(i)-loc_w)*(loc_gi_diff(i))
  end if
end do
!!
!!little_q_bar = little_q_bar/loc_temp
!
!!write(*,*) 'In heat flux sub',loc_u,loc_v,loc_w,loc_temp,little_q_bar
!
select case(dir)

case(39)

do i = 1,dir
  contr_heat_flux(i) = q_lead*(cx(i)*little_q_bar(1)+cy(i)*little_q_bar(2)+&
    cz(i)*little_q_bar(3))
end do

end select

end subroutine

!! First shell
!! dir 1
!contr_heat_flux(1) = q_lead*(cx(1)*little_q_bar(1)+cy(!0.0D0
!! Second shell
!contr_heat_flux(2) = init_epsilon(lvl)*cx(2)*little_q_bar(1)/loc_temp
!contr_heat_flux(3) = init_epsilon(lvl)*cy(3)*little_q_bar(2)/loc_temp
!contr_heat_flux(4) = init_epsilon(lvl)*cz(4)*little_q_bar(3)/loc_temp
!contr_heat_flux(5) = init_epsilon(lvl)*cx(5)*little_q_bar(1)/loc_temp
!contr_heat_flux(6) = init_epsilon(lvl)*cy(6)*little_q_bar(2)/loc_temp
!contr_heat_flux(7) = init_epsilon(lvl)*cz(7)*little_q_bar(3)/loc_temp
!! Third shell
!contr_heat_flux(8) = init_epsilon(lvl)*(little_q_bar(1)+little_q_bar(2)+little_q_bar(3))/loc_temp
!contr_heat_flux(9) = -init_epsilon(lvl)*(little_q_bar(1)+little_q_bar(2)+little_q_bar(3))/loc_temp
!contr_heat_flux(10) = -init_epsilon(lvl)*(little_q_bar(1)-little_q_bar(2)+little_q_bar(3))/loc_temp
!contr_heat_flux(11) = init_epsilon(lvl)*(little_q_bar(1)-little_q_bar(2)+little_q_bar(3))/loc_temp
!contr_heat_flux(12) = init_epsilon(lvl)*(little_q_bar(1)+little_q_bar(2)-little_q_bar(3))/loc_temp
!contr_heat_flux(13) = -init_epsilon(lvl)*(little_q_bar(1)+little_q_bar(2)-little_q_bar(3))/loc_temp
!contr_heat_flux(14) = -init_epsilon(lvl)*(little_q_bar(1)-little_q_bar(2)-little_q_bar(3))/loc_temp
!contr_heat_flux(15) = init_epsilon(lvl)*(little_q_bar(1)-little_q_bar(2)-little_q_bar(3))/loc_temp
!! Fourth shell
!contr_heat_flux(16) = init_epsilon(lvl)*little_q_bar(1)*2.0D0/loc_temp
!contr_heat_flux(17) = init_epsilon(lvl)*little_q_bar(2)*2.0D0/loc_temp
!contr_heat_flux(18) = init_epsilon(lvl)*little_q_bar(3)*2.0D0/loc_temp
!contr_heat_flux(19) = -init_epsilon(lvl)*little_q_bar(1)*2.0D0/loc_temp
!contr_heat_flux(20) = -init_epsilon(lvl)*little_q_bar(2)*2.0D0/loc_temp
!contr_heat_flux(21) = -init_epsilon(lvl)*little_q_bar(3)*2.0D0/loc_temp
!! Fifth shell
!contr_heat_flux(22) = 2.0D0*init_epsilon(lvl)*(little_q_bar(1) + little_q_bar(2))/loc_temp
!contr_heat_flux(23) = 2.0D0*init_epsilon(lvl)*(little_q_bar(1) + little_q_bar(3))/loc_temp
!contr_heat_flux(24) = 2.0D0*init_epsilon(lvl)*(little_q_bar(3) + little_q_bar(2))/loc_temp
!contr_heat_flux(25) = 2.0D0*init_epsilon(lvl)*(-little_q_bar(1) + little_q_bar(2))/loc_temp
!contr_heat_flux(26) = -2.0D0*init_epsilon(lvl)*(little_q_bar(1) + little_q_bar(2))/loc_temp
!contr_heat_flux(27) = 2.0D0*init_epsilon(lvl)*(little_q_bar(1) - little_q_bar(2))/loc_temp
!contr_heat_flux(28) = 2.0D0*init_epsilon(lvl)*(little_q_bar(1) - little_q_bar(3))/loc_temp
!contr_heat_flux(29) = -2.0D0*init_epsilon(lvl)*(little_q_bar(1) + little_q_bar(3))/loc_temp
!contr_heat_flux(30) = 2.0D0*init_epsilon(lvl)*(-little_q_bar(1) + little_q_bar(3))/loc_temp
!contr_heat_flux(31) = 2.0D0*init_epsilon(lvl)*(little_q_bar(2) - little_q_bar(3))/loc_temp
!contr_heat_flux(32) = -2.0D0*init_epsilon(lvl)*(little_q_bar(2) + little_q_bar(3))/loc_temp
!contr_heat_flux(33) = 2.0D0*init_epsilon(lvl)*(-little_q_bar(2) + little_q_bar(3))/loc_temp
!! Sixth shell
!contr_heat_flux(34) = init_epsilon(lvl)*little_q_bar(1)*3.0D0/loc_temp
!contr_heat_flux(35) = init_epsilon(lvl)*little_q_bar(2)*3.0D0/loc_temp
!contr_heat_flux(36) = init_epsilon(lvl)*little_q_bar(3)*3.0D0/loc_temp
!contr_heat_flux(37) = -init_epsilon(lvl)*little_q_bar(1)*3.0D0/loc_temp
!contr_heat_flux(38) = -init_epsilon(lvl)*little_q_bar(2)*3.0D0/loc_temp
!contr_heat_flux(39) = -init_epsilon(lvl)*little_q_bar(3)*3.0D0/loc_temp












!
!!do i = 1,dir
!!  if (contr_heat_flux(i) > 1.0D-7) then
!!    write(*,*) 'contracted heat fluxes large',contr_heat_flux(i),little_q_bar,loc_u,loc_v,loc_w,&
!!      loc_temp
!!  end if
!!end do
!
!end select
!contr_heat_flux(1) = 0.0D0
!!! Second shell
!!! dir 2
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(2) = init_epsilon(lvl)*x_value/loc_temp
!!if (lvl == 3) then
!!  write(*,*) 'contracted heat flux 2',x_value,u_term,loc_u,contr_heat_flux(2),loc_gi_diff
!!end if
!! dir 3
!y_value = 0.0D0
!v_term = 1.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(3) = init_epsilon(lvl)*y_value/loc_temp
!
!! dir 4
!z_value = 0.0D0
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(4) = init_epsilon(lvl)*z_value/loc_temp
!! dir 5
!x_value = 0.0D0
!u_term = -1.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(5) = -init_epsilon(lvl)*x_value/loc_temp
!! dir 6
!y_value = 0.0D0
!v_term = -1.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(6) = -init_epsilon(lvl)*y_value/loc_temp
!! dir 7
!z_value = 0.0D0
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(7) = -init_epsilon(lvl)*z_value/loc_temp
!! Third Shell
!! dir 8
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(8) = init_epsilon(lvl)*(x_value+y_value+z_value)/loc_temp
!! dir 9
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(9) = init_epsilon(lvl)*(-x_value + y_value + z_value)/loc_temp
!! dir 10
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(10) = init_epsilon(lvl)*(-x_value - y_value + z_value)/loc_temp
!! dir 11
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(11) = init_epsilon(lvl)*(x_value - y_value + z_value)/loc_temp
!! dir 12
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(12) = init_epsilon(lvl)*(x_value + y_value - z_value)/loc_temp
!! dir 13
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(13) = init_epsilon(lvl)*(-x_value + y_value - z_value)/loc_temp
!! dir 14
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(14) = -init_epsilon(lvl)*(x_value + y_value + z_value)/loc_temp
!! dir 15
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(15) = init_epsilon(lvl)*(x_value - y_value - z_value)/loc_temp
!! Fourth Shell
!! dir 16
!x_value = 0.0D0
!u_term = 2.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(16) = 2.0D0*init_epsilon(lvl)*x_value/loc_temp
!! dir 17
!y_value = 0.0D0
!v_term = 2.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(17) = 2.0D0*init_epsilon(lvl)*y_value/loc_temp
!! dir 18
!z_value = 0.0D0
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(18) = 2.0D0*init_epsilon(lvl)*z_value/loc_temp
!! dir 19
!x_value = 0.0D0
!u_term = -2.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(19) = -2.0D0*init_epsilon(lvl)*x_value/loc_temp
!! dir 20
!y_value = 0.0D0
!v_term = -2.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(20) = -2.0D0*init_epsilon(lvl)*y_value/loc_temp
!! dir 21
!z_value = 0.0D0
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(21) = -2.0D0*init_epsilon(lvl)*z_value/loc_temp
!! Fifth Shell
!! dir 22, +2,+2,0
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = 2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(22) = 2.0D0*init_epsilon(lvl)*(x_value + y_value)/loc_temp
!! dir 23, +2,0,+2
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = 2.0D0-loc_u
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(23) = 2.0D0*init_epsilon(lvl)*(x_value + z_value)/loc_temp
!! dir 24, 0,+2,+2
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = 2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(24) = 2.0D0*init_epsilon(lvl)*(y_value + z_value)/loc_temp
!! dir 25
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = 2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(25) = 2.0D0*init_epsilon(lvl)*(-x_value + y_value)/loc_temp
!! dir 26
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(26) = -2.0D0*init_epsilon(lvl)*(x_value + y_value)/loc_temp
!! dir 27
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(27) = 2.0D0*init_epsilon(lvl)*(x_value - y_value)/loc_temp
!! dir 28
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = 2.0D0-loc_u
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(28) = 2.0D0*init_epsilon(lvl)*(x_value - z_value)/loc_temp
!! dir 29
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = -2.0D0-loc_u
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(29) = -2.0D0*init_epsilon(lvl)*(x_value + z_value)/loc_temp
!! dir 30
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = -2.0D0-loc_u
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(30) = 2.0D0*init_epsilon(lvl)*(-x_value + z_value)/loc_temp
!! dir 31
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = 2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(31) = 2.0D0*init_epsilon(lvl)*(y_value - z_value)/loc_temp
!! dir 32
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = -2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(32) = -2.0D0*init_epsilon(lvl)*(y_value + z_value)/loc_temp
!! dir 33
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = -2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(33) = 2.0D0*init_epsilon(lvl)*(-y_value + z_value)/loc_temp
!! Sixth Shell
!! dir 34
!x_value = 0.0D0
!u_term = 3.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(34) = 3.0D0*init_epsilon(lvl)*x_value/loc_temp
!! dir 35
!y_value = 0.0D0
!v_term = 3.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(35) = 3.0D0*init_epsilon(lvl)*y_value/loc_temp
!! dir 36
!z_value = 0.0D0
!w_term = 3.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(36) = 3.0D0*init_epsilon(lvl)*z_value/loc_temp
!! dir 37
!x_value = 0.0D0
!u_term = -3.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(37) = -3.0D0*init_epsilon(lvl)*x_value/loc_temp
!! dir 38
!y_value = 0.0D0
!v_term = -3.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(38) = -3.0D0*init_epsilon(lvl)*y_value/loc_temp
!! dir 39
!z_value = 0.0D0
!w_term = -3.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(39) = -3.0D0*init_epsilon(lvl)*z_value/loc_temp

!
!do i = 1,dir
!  if (abs(contr_heat_flux(i)) > 10.0D0) then
!    write(*,*) 'In heat flux sub',loc_u,loc_v,loc_w,loc_temp,contr_heat_flux(i),&
!      loc_gi_diff(i),i,rcx(i),rcy(i),rcz(i)
!!    write(*,*) 'what the actual fuck',rcx(i)-loc_u,rcy(i)-loc_v,rcz(i)-loc_w,loc_gi(a,b,c,i)-loc_gieq(i)
!  end if
!end do





! First shell
! dir 1
!contr_heat_flux(1) = 0.0D0
!! Second shell
!! dir 2
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(2) = x_value/loc_temp
!! dir 3
!y_value = 0.0D0
!v_term = 1.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(3) = y_value/loc_temp
!
!! dir 4
!z_value = 0.0D0
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(4) = z_value/loc_temp
!! dir 5
!x_value = 0.0D0
!u_term = -1.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(5) = -x_value/loc_temp
!! dir 6
!y_value = 0.0D0
!v_term = -1.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(6) = -y_value/loc_temp
!! dir 7
!z_value = 0.0D0
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(7) = -z_value/loc_temp
!! Third Shell
!! dir 8
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(8) = (x_value+y_value+z_value)/loc_temp
!! dir 9
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(9) = (-x_value + y_value + z_value)/loc_temp
!! dir 10
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(10) = (-x_value - y_value + z_value)/loc_temp
!! dir 11
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(11) = (x_value - y_value + z_value)/loc_temp
!! dir 12
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(12) = (x_value + y_value - z_value)/loc_temp
!! dir 13
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(13) = (-x_value + y_value - z_value)/loc_temp
!! dir 14
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(14) = -(x_value + y_value + z_value)/loc_temp
!! dir 15
!x_value = 0.0D0
!y_value = 0.0D0
!z_value = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(15) = (x_value - y_value - z_value)/loc_temp
!! Fourth Shell
!! dir 16
!x_value = 0.0D0
!u_term = 2.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(16) = 2.0D0*x_value/loc_temp
!! dir 17
!y_value = 0.0D0
!v_term = 2.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(17) = 2.0D0*y_value/loc_temp
!! dir 18
!z_value = 0.0D0
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(18) = 2.0D0*z_value/loc_temp
!! dir 19
!x_value = 0.0D0
!u_term = -2.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(19) = -2.0D0*x_value/loc_temp
!! dir 20
!y_value = 0.0D0
!v_term = -2.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(20) = -2.0D0*y_value/loc_temp
!! dir 21
!z_value = 0.0D0
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(21) = -2.0D0*z_value/loc_temp
!! Fifth Shell
!! dir 22, +2,+2,0
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = 2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(22) = 2.0D0*(x_value + y_value)/loc_temp
!! dir 23, +2,0,+2
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = 2.0D0-loc_u
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(23) = 2.0D0*(x_value + z_value)/loc_temp
!! dir 24, 0,+2,+2
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = 2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(24) = 2.0D0*(y_value + z_value)/loc_temp
!! dir 25
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = 2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(25) = 2.0D0*(-x_value + y_value)/loc_temp
!! dir 26
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(26) = -2.0D0*(x_value + y_value)/loc_temp
!! dir 27
!x_value = 0.0D0
!y_value = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -2.0D0-loc_v
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(27) = 2.0D0*(x_value - y_value)/loc_temp
!! dir 28
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = 2.0D0-loc_u
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(28) = 2.0D0*(x_value - z_value)/loc_temp
!! dir 29
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = -2.0D0-loc_u
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(29) = -2.0D0*(x_value + z_value)/loc_temp
!! dir 30
!x_value = 0.0D0
!z_value = 0.0D0
!u_term = -2.0D0-loc_u
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(30) = 2.0D0*(-x_value + z_value)/loc_temp
!! dir 31
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = 2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(31) = 2.0D0*(y_value - z_value)/loc_temp
!! dir 32
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = -2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(32) = -2.0D0*(y_value + z_value)/loc_temp
!! dir 33
!y_value = 0.0D0
!z_value = 0.0D0
!v_term = -2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(33) = 2.0D0*(-y_value + z_value)/loc_temp
!! Sixth Shell
!! dir 34
!x_value = 0.0D0
!u_term = 3.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(34) = 3.0D0*x_value/loc_temp
!! dir 35
!y_value = 0.0D0
!v_term = 3.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(35) = 3.0D0*y_value/loc_temp
!! dir 36
!z_value = 0.0D0
!w_term = 3.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(36) = 3.0D0*z_value/loc_temp
!! dir 37
!x_value = 0.0D0
!u_term = -3.0D0-loc_u
!do i = 1,dir
!  x_value = x_value + u_term*loc_gi_diff(i)
!end do
!contr_heat_flux(37) = -3.0D0*x_value/loc_temp
!! dir 38
!y_value = 0.0D0
!v_term = -3.0D0-loc_v
!do i = 1,dir
!  y_value = y_value + v_term*loc_gi_diff(i)
!end do
!contr_heat_flux(38) = -3.0D0*y_value/loc_temp
!! dir 39
!z_value = 0.0D0
!w_term = -3.0D0-loc_w
!do i = 1,dir
!  z_value = z_value + w_term*loc_gi_diff(i)
!end do
!contr_heat_flux(39) = -3.0D0*z_value/loc_temp
!end select

