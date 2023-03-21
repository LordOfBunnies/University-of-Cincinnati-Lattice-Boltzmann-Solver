SUBROUTINE find_q_bar(loc_fidiff,loc_fieq,q_bar_dir,&
  loc_u,loc_v,loc_w,loc_temp,lvl)
!
! Find the heat flux and equilibrium heat flux
!
!
!
! Called by: collision
! Calls:
!
use precise
use constants
use amr_info_holder, only: dt,timestep
use amr_processes, only: epsilon_calculator,amr_max_lvl
use linkwise
use grid_data
IMPLICIT NONE
REAL(kind=dp),intent(in) :: loc_fidiff(dir),loc_fieq(dir)
real(kind=dp) :: u_term,v_term,w_term,loc_epsilon
real(kind=dp) :: q_array
real(kind=dp) :: q_bar_dir(dir),xxx_diff,yyy_diff,zzz_diff,xxy_diff,xxz_diff,&
  xyy_diff,xyz_diff,xzz_diff,yyz_diff,yzz_diff
real(kind=dp) :: qabc(27)

real(kind=dp),intent(in) :: loc_u,loc_v,loc_w,loc_temp
real(kind=dp) :: temp_coeff
!real(kind=dp) :: heat_dir_1,heat_dir_111,heat_dir_p11n1,heat_dir_2,heat_dir_p22,&
!  heat_dir_3
INTEGER :: i,lvl
!
! Set the sums to 0
!
q_bar_dir(1) = 0.0D0
temp_coeff = (6.0D0*loc_temp**3)

loc_epsilon = epsilon_calculator(loc_u,loc_v,loc_w,loc_temp,&
  timestep(amr_max_lvl),-4,-5,-6,lvl)

select case (dir)
case(19)

case(39)

!!
!
xxx_diff = 0.0D0
yyy_diff = 0.0D0
zzz_diff = 0.0D0
xxy_diff = 0.0D0
xxz_diff = 0.0D0
xyy_diff = 0.0D0
xyz_diff = 0.0D0
xzz_diff = 0.0D0
yyz_diff = 0.0D0
yzz_diff = 0.0D0

do i = 1,dir
  xxx_diff = xxx_diff + (rcx(i)-loc_u)**3*loc_fidiff(i)
  yyy_diff = yyy_diff + (rcy(i)-loc_v)**3*loc_fidiff(i)
  zzz_diff = zzz_diff + (rcz(i)-loc_w)**3*loc_fidiff(i)
!
  xxy_diff = xxy_diff + (rcx(i)-loc_u)**2 * (rcy(i)-loc_v) * loc_fidiff(i)
!
  xxz_diff = xxz_diff + (rcx(i)-loc_u)**2 * (rcz(i)-loc_w) * loc_fidiff(i)
!
  xyy_diff = xyy_diff + (rcy(i)-loc_v)**2 * (rcx(i)-loc_u) * loc_fidiff(i)
!
  xzz_diff = xzz_diff + (rcz(i)-loc_w)**2 * (rcx(i)-loc_u) * loc_fidiff(i)
!
  yyz_diff = yyz_diff + (rcy(i)-loc_v)**2 * (rcz(i)-loc_w) * loc_fidiff(i)
!
  yzz_diff = yzz_diff + (rcz(i)-loc_w)**2 * (rcy(i)-loc_v) * loc_fidiff(i)
!
  xyz_diff = xyz_diff + (rcx(i)-loc_u)*(rcy(i)-loc_v)* (rcz(i)-loc_w) * loc_fidiff(i)

end do

!if (.not. shifted) then
!!Second shell
!! dir 2
!q_bar_dir(2) = loc_epsilon*(xxx_diff*(cx(2)**2-3*cx(2)*loc_temp)-&
!                 3*cx(2)*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(3) = loc_epsilon*(yyy_diff*(cy(3)**2-3*cy(3)*loc_temp)-&
!                 3*cy(3)*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(4) = loc_epsilon*(zzz_diff*(cz(4)**2-3*cz(4)*loc_temp)-&
!                 3*cz(4)*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(5) = loc_epsilon*(xxx_diff*(cx(5)**2-3*cx(5)*loc_temp)-&
!                 3*cx(5)*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(6) = loc_epsilon*(yyy_diff*(cy(6)**2-3*cy(6)*loc_temp)-&
!                 3*cy(6)*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(7) = loc_epsilon*(zzz_diff*(cz(7)**2-3*cz(7)*loc_temp)-&
!                 3*cz(7)*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!!Third Shell
!q_bar_dir(8) = loc_epsilon*((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+yyz_diff+yzz_diff)+&
!                 6*xyz_diff+&
!                 2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(9) = loc_epsilon*((1-3*loc_temp)*(-xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff-xyy_diff-xzz_diff+yyz_diff+yzz_diff)-&
!                 6*xyz_diff+&
!                 2*(xxy_diff-xyy_diff+xxz_diff-xzz_diff+yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(10) = loc_epsilon*((1-3*loc_temp)*(-xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff-xyy_diff-xzz_diff+yyz_diff-yzz_diff)+&
!                 6*xyz_diff+&
!                 2*(-xxy_diff-xyy_diff+xxz_diff-xzz_diff+yyz_diff-yzz_diff))/temp_coeff
!q_bar_dir(11) = loc_epsilon*((1-3*loc_temp)*(xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff+xyy_diff+xzz_diff+yyz_diff-yzz_diff)-&
!                  6*xyz_diff+&
!                  2*(-xxy_diff+xyy_diff+xxz_diff+xzz_diff+yyz_diff-yzz_diff))/temp_coeff
!q_bar_dir(12) = loc_epsilon*((1-3*loc_temp)*(xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff+xyy_diff+xzz_diff-yyz_diff+yzz_diff)-&
!                  6*xyz_diff+&
!                  2*(xxy_diff+xyy_diff-xxz_diff+xzz_diff-yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(13) = loc_epsilon*((1-3*loc_temp)*(-xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff-xyy_diff-xzz_diff-yyz_diff+yzz_diff)+&
!                  6*xyz_diff+&
!                  2*(xxy_diff-xyy_diff-xxz_diff-xzz_diff-yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(14) = -loc_epsilon*((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+yyz_diff+yzz_diff)+&
!                   6*xyz_diff+&
!                   2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(15) = loc_epsilon*((1-3*loc_temp)*(xxx_diff-yyy_diff-zzz_diff-xxy_diff-xxz_diff+xyy_diff+xzz_diff- yyz_diff-yzz_diff)+&
!                  6*xyz_diff + &
!                  2*(-xxy_diff+xyy_diff-xxz_diff+xzz_diff-yyz_diff-yzz_diff))/temp_coeff
!! Fourth shell
!q_bar_dir(16) = loc_epsilon*(xxx_diff*(cx(16)**2-3*cx(16)*loc_temp)-&
!                 3*cx(16)*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(17) = loc_epsilon*(yyy_diff*(cy(17)**2-3*cy(17)*loc_temp)-&
!                 3*cy(17)*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(18) = loc_epsilon*(zzz_diff*(cz(18)**2-3*cz(18)*loc_temp)-&
!                 3*cz(18)*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(19) = loc_epsilon*(xxx_diff*(cx(19)**2-3*cx(19)*loc_temp)-&
!                 3*cx(19)*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(20) = loc_epsilon*(yyy_diff*(cy(20)**2-3*cy(20)*loc_temp)-&
!                 3*cy(20)*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(21) = loc_epsilon*(zzz_diff*(cz(21)**2-3*cz(21)*loc_temp)-&
!                 3*cz(21)*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!!q_bar_dir(16) = loc_epsilon*(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!!q_bar_dir(17) = loc_epsilon*(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!!q_bar_dir(18) = loc_epsilon*(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!!q_bar_dir(19) = -loc_epsilon*(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!!q_bar_dir(20) = -loc_epsilon*(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!!q_bar_dir(21) = -loc_epsilon*(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!! Fifth shell
!q_bar_dir(22) = loc_epsilon*((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!q_bar_dir(23) = loc_epsilon*((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!q_bar_dir(24) = loc_epsilon*((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!q_bar_dir(25) = loc_epsilon*((8-6*loc_temp)*(-xxx_diff+yyy_diff+xxy_diff-xyy_diff)+&
!                  16*(xxy_diff-xyy_diff)&
!                  -6*loc_temp*(-xzz_diff+yzz_diff))/temp_coeff
!q_bar_dir(26) = -loc_epsilon*((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!q_bar_dir(27) = loc_epsilon*((8-6*loc_temp)*(xxx_diff-yyy_diff-xxy_diff+xyy_diff)+&
!                  16*(-xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff-yzz_diff))/temp_coeff
!q_bar_dir(28) = loc_epsilon*((8-6*loc_temp)*(xxx_diff-zzz_diff-xxz_diff+xzz_diff)+&
!                 16*(-xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff-yyz_diff))/temp_coeff
!q_bar_dir(29) = -loc_epsilon*((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!q_bar_dir(30) = loc_epsilon*((8-6*loc_temp)*(-xxx_diff+zzz_diff+xxz_diff-xzz_diff)+&
!                 16*(xxz_diff-xzz_diff)-&
!                 6*loc_temp*(-xyy_diff+yyz_diff))/temp_coeff
!q_bar_dir(31) = loc_epsilon*((8-6*loc_temp)*(yyy_diff-zzz_diff-yyz_diff+yzz_diff)+&
!                 16*(-yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff-xxz_diff))/temp_coeff
!q_bar_dir(32) = -loc_epsilon*((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!q_bar_dir(33) = loc_epsilon*((8-6*loc_temp)*(-yyy_diff+zzz_diff+yyz_diff-yzz_diff)+&
!                 16*(yyz_diff-yzz_diff)-&
!                 6*loc_temp*(-xxy_diff+xxz_diff))/temp_coeff
!! Sixth shell
!q_bar_dir(34) = loc_epsilon*(xxx_diff*(cx(34)**2-3*cx(34)*loc_temp)-&
!                 3*cx(34)*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(35) = loc_epsilon*(yyy_diff*(cy(35)**2-3*cy(35)*loc_temp)-&
!                 3*cy(35)*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(36) = loc_epsilon*(zzz_diff*(cz(36)**2-3*cz(36)*loc_temp)-&
!                 3*cz(36)*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(37) = loc_epsilon*(xxx_diff*(cx(37)**2-3*cx(37)*loc_temp)-&
!                 3*cx(37)*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(38) = loc_epsilon*(yyy_diff*(cy(38)**2-3*cy(38)*loc_temp)-&
!                 3*cy(38)*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(39) = loc_epsilon*(zzz_diff*(cz(39)**2 - 3*cz(39)*loc_temp)-&
!                 3*cz(39)*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!else
qabc(1) = xxx_diff
qabc(2) = xxy_diff
qabc(3) = xxz_diff
qabc(4) = xxy_diff
qabc(5) = xyy_diff
qabc(6) = xyz_diff
qabc(7) = xxz_diff
qabc(8) = xyz_diff
qabc(9) = xzz_diff
qabc(10) = xxy_diff
qabc(11) = xyy_diff
qabc(12) = xyz_diff
qabc(13) = xyy_diff
qabc(14) = yyy_diff
qabc(15) = yyz_diff
qabc(16) = xyz_diff
qabc(17) = yyz_diff
qabc(18) = yzz_diff
qabc(19) = xxz_diff
qabc(20) = xyz_diff
qabc(21) = xzz_diff
qabc(22) = xyz_diff
qabc(23) = yyz_diff
qabc(24) = yzz_diff
qabc(25) = xzz_diff
qabc(26) = yzz_diff
qabc(27) = zzz_diff

do i = 1,39
  q_bar_dir(i) = loc_epsilon*(sum(qabc*yi(:,i)) - 3*loc_temp*sum(qabc*kroen_val(:,i)))/temp_coeff
end do

!end if

!q_bar_dir(34) = loc_epsilon*(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(35) = loc_epsilon*(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(36) = loc_epsilon*(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(37) = -loc_epsilon*(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(38) = -loc_epsilon*(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(39) = -loc_epsilon*(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
end select

!do i = 1,dir
!  if (abs(q_bar_dir(i)) >1.0D-7) then
!    write(*,*) 'q_bar_dir values',q_bar_dir(i),loc_u,loc_v,loc_w,loc_temp,i
!!    write(*,*) 'local values',loc_rho,loc_u,loc_v,loc_w,loc_temp
!    write(*,*) 'associated fi diffs',q_bar_dir(i),loc_fidiff,i
!    write(*,*) 'all together q bars',q_bar_dir
!  end if
!end do

!if (lvl == 3) then
!  write(*,*) 'q_bar_dir',q_bar_dir
!end if

!do i = 1,dir
!  if (abs(q_bar_dir(i)) > 100.0D0 .or. isnan(q_bar_dir(i))) then
!    write(*,*) 'In q_bar_dir sub',loc_u,loc_v,loc_w,loc_temp,loc_fidiff(i),loc_fieq(i),q_bar_dir(i)
!  end if
!end do

END SUBROUTINE
! dir 2
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!!do i = 1,dir
!!  xxx_diff = xxx_diff + (rcx(i)-loc_u)**3*loc_fidiff(i)
!!  xyy_diff = xyy_diff + (rcx(i)-loc_u)*(rcy(i)-loc_v)**2*loc_fidiff(i)
!!  xzz_diff = xzz_diff + (rcx(i)-loc_u)*(rcz(i)-loc_w)**2*loc_fidiff(i)
!!end do
!q_bar_dir(2) = init_epsilon(lvl)*(xxx_diff*(1-3*loc_temp)-3*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!!if (lvl == 3 .and. abs(q_bar_dir(2)) > 1.0D-8) then
!!  write(*,*) 'q_bar_dir 2',xxx_diff,xyy_diff,xzz_diff,u_term,v_term,w_term,loc_temp,q_bar_dir(2),loc_fidiff
!!end if
!! dir 3
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 1.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(3) = init_epsilon(lvl)*(yyy_diff*(1-3*loc_temp)-3*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!
!! dir 4
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(4) = init_epsilon(lvl)*(zzz_diff*(1-3*loc_temp)-3*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 5
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(5) = -init_epsilon(lvl)*(xxx_diff*(1-3*loc_temp)-3*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!!if (lvl == 3 .and. abs(q_bar_dir(5)) > 1.0D-8) then
!!  write(*,*) 'q_bar_dir 5',xxx_diff,xyy_diff,xzz_diff,u_term,v_term,w_term,loc_temp,q_bar_dir(5),loc_fidiff
!!end if
!! dir 6
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -1.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(6) = -init_epsilon(lvl)*(yyy_diff*(1-3*loc_temp)-3*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!! dir 7
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(7) = -init_epsilon(lvl)*(zzz_diff*(1-3*loc_temp)-3*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!! dir 8
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(8) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 9
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(9) = init_epsilon(lvl)*((1-3*loc_temp)*(-xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff-xyy_diff-xzz_diff+ &
!                   yyz_diff+yzz_diff)-6*xyz_diff+2*(xxy_diff-xyy_diff+xxz_diff-xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 10
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(10) = init_epsilon(lvl)*((1-3*loc_temp)*(-xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff-xyy_diff-xzz_diff+ &
!                   yyz_diff-yzz_diff)+6*xyz_diff+2*(-xxy_diff-xyy_diff+xxz_diff-xzz_diff+&
!                   yyz_diff-yzz_diff))/temp_coeff
!
!
!! dir 11
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(11) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff-yzz_diff)-6*xyz_diff+2*(-xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff-yzz_diff))/temp_coeff
!
!
!! dir 12
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(12) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff+xyy_diff+xzz_diff- &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff-xxz_diff+xzz_diff-&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 13
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(13) = init_epsilon(lvl)*((1-3*loc_temp)*(-xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff-xyy_diff-xzz_diff- &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff-xyy_diff-xxz_diff-xzz_diff-&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 14
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(14) = -init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 15
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(15) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff-yyy_diff-zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!
!! dir 16
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(16) = init_epsilon(lvl)*(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 17
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(17) = init_epsilon(lvl)*(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!
!! dir 18
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(18) = init_epsilon(lvl)*(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 19
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(19) = -init_epsilon(lvl)*(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 20
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(20) = -init_epsilon(lvl)*(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!! dir 21
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(21) = -init_epsilon(lvl)*(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 22, +2,+2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = 2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(22) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!
!! dir 23, +2,0,+2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(23) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!
!! dir 24, 0,+2,+2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(24) = init_epsilon(lvl)*((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!
!! dir 25, -2,+2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = 2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(25) = init_epsilon(lvl)*((8-6*loc_temp)*(-xxx_diff+yyy_diff+xxy_diff-xyy_diff)+&
!                  16*(xxy_diff-xyy_diff)&
!                  -6*loc_temp*(-xzz_diff+yzz_diff))/temp_coeff
!
!! dir 26, -2,-2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(26) = -init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!
!! dir 27,+2,-2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(27) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff-yyy_diff-xxy_diff+xyy_diff)+&
!                  16*(-xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff-yzz_diff))/temp_coeff
!
!! dir 28 +2,0,-2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(28) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff-zzz_diff-xxz_diff+xzz_diff)+&
!                 16*(-xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff-yyz_diff))/temp_coeff
!
!! dir 29, -2,0,-2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(29) = -init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!
!! dir 30, -2,0,+2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(30) = init_epsilon(lvl)*((8-6*loc_temp)*(-xxx_diff+zzz_diff+xxz_diff-xzz_diff)+&
!                 16*(xxz_diff-xzz_diff)-&
!                 6*loc_temp*(-xyy_diff+yyz_diff))/temp_coeff
!
!! dir 31, 0,+2,-2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(31) = init_epsilon(lvl)*((8-6*loc_temp)*(yyy_diff-zzz_diff-yyz_diff+yzz_diff)+&
!                 16*(-yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff-xxz_diff))/temp_coeff
!
!! dir 32, 0,-2,-2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(32) = -init_epsilon(lvl)*((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!
!! dir 33, 0,-2,+2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(33) = init_epsilon(lvl)*((8-6*loc_temp)*(-yyy_diff+zzz_diff+yyz_diff-yzz_diff)+&
!                 16*(yyz_diff-yzz_diff)-&
!                 6*loc_temp*(-xxy_diff+xxz_diff))/temp_coeff
!
!! dir 34
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = 3.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(34) = init_epsilon(lvl)*(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 35
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 3.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(35) = init_epsilon(lvl)*(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!
!! dir 36
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = 3.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(36) = init_epsilon(lvl)*(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 37
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = -3.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(37) = -init_epsilon(lvl)*(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!! dir 38
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -3.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(38) = -init_epsilon(lvl)*(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!! dir 39
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = -3.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(39) = -init_epsilon(lvl)*(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!



!  q_bar_diff(1) = q_bar_diff(1) + (rcx(i)-loc_u)*loc_fidiff(i)
!  q_bar_diff(2) = q_bar_diff(2) + (rcy(i)-loc_v)*loc_fidiff(i)
!  q_bar_diff(3) = q_bar_diff(3) + (rcz(i)-loc_w)*loc_fidiff(i)
!DO i = 1,dir
!!        big_q_bar = big_q_bar + ((rcx(i)-loc_u) * (rcy(i)-loc_v) * &
!!          (rcz(i)-loc_w))*loc_fi(i)
!!        big_q_bar_eq = big_q_bar_eq + ((rcx(i)-loc_u) * (rcy(i)-loc_v) * &
!!          (rcz(i)-loc_w))*loc_fieq(i)
!
!  q_bar_diff = q_bar_diff - ((rcx(i)-loc_u) * (rcy(i)-loc_v) * &
!    (rcz(i)-loc_w))*loc_fidiff(i)
!
!END DO
!
!!      q_bar_diff = q_bar_diff/(6.0D0*loc_temp**3)
!
!!      DO i = 1,dir
!
!select case (dir)
!!          case(1)
!  case (39)
!  heat_dir_1 = (1.0D0-9*loc_temp)/(6.0D0*loc_temp**3)
!  heat_dir_111 = (27.0D0-27.0D0*loc_temp)/(6.0D0*loc_temp**3)
!  heat_dir_p11n1 = (1.0D0-9.0D0*loc_temp)/(6.0D0*loc_temp**3)
!  heat_dir_2 = (8.0D0-18.0D0*loc_temp)/(6.0D0*loc_temp**3)
!  heat_dir_p22 = (64.0D0-36.0D0*loc_temp)/(6.0D0*loc_temp**3)
!  heat_dir_3 = (27.0D0-27.0D0*loc_temp)/(6.0D0*loc_temp**3)
!
!   q_bar_dir(1) = 0.0D0
!! case(2)
!   q_bar_dir(2) = q_bar_diff*heat_dir_1
!! case(3)
!   q_bar_dir(3) = q_bar_diff*heat_dir_1
!! case(4)
!   q_bar_dir(4) = q_bar_diff*heat_dir_1
!! case(5)
!   q_bar_dir(5) = -q_bar_diff*heat_dir_1
!! case(6)
!   q_bar_dir(6) = -q_bar_diff*heat_dir_1
!! case(7)
!   q_bar_dir(7) = -q_bar_diff*heat_dir_1
!! case(8)
!   q_bar_dir(8) = q_bar_diff*heat_dir_111
!! case(9)
!   q_bar_dir(9) = q_bar_diff*heat_dir_p11n1
!! case(10)
!   q_bar_dir(10) = -q_bar_diff*heat_dir_p11n1
!! case(11)
!   q_bar_dir(11) = q_bar_diff*heat_dir_p11n1
!! case(12)
!   q_bar_dir(12) = q_bar_diff*heat_dir_p11n1
!! case(13)
!   q_bar_dir(13) = -q_bar_diff*heat_dir_p11n1
!! case(14)
!   q_bar_dir(14) = -q_bar_diff*heat_dir_111
!! case(15)
!   q_bar_dir(15) = -q_bar_diff*heat_dir_p11n1
!! case(16)
!   q_bar_dir(16) = q_bar_diff*heat_dir_2
!! case(17)
!   q_bar_dir(17) = q_bar_diff*heat_dir_2
!! case(18)
!   q_bar_dir(18) = q_bar_diff*heat_dir_2
!! case(19)
!   q_bar_dir(19) = -q_bar_diff*heat_dir_2
!! case(20)
!   q_bar_dir(20) = -q_bar_diff*heat_dir_2
!! case(21)
!   q_bar_dir(21) = -q_bar_diff*heat_dir_2
!! case(22)
!   q_bar_dir(22) = q_bar_diff*heat_dir_p22
!! case(23)
!   q_bar_dir(23) = q_bar_diff*heat_dir_p22
!! case(24)
!   q_bar_dir(24) = q_bar_diff*heat_dir_p22
!! case(25)
!   q_bar_dir(25) = 0.0D0!-q_bar_diff*heat_dir_2
!! case(26)
!   q_bar_dir(26) = -q_bar_diff*heat_dir_p22
!! case(27)
!   q_bar_dir(27) = 0.0D0!-q_bar_diff*heat_dir_2
!! case(28)
!   q_bar_dir(28) = 0.0D0!-q_bar_diff*heat_dir_2
!! case(29)
!   q_bar_dir(29) = -q_bar_diff*heat_dir_p22
!! case(30)
!   q_bar_dir(30) = 0.0D0!-q_bar_diff*heat_dir_2
!! case(31)
!   q_bar_dir(31) = 0.0D0!-q_bar_diff*heat_dir_2
!! case(32)
!   q_bar_dir(32) = -q_bar_diff*heat_dir_p22
!! case(33)
!   q_bar_dir(33) = 0.0D0!-q_bar_diff*heat_dir_2
!! case(34)
!   q_bar_dir(34) = q_bar_diff*heat_dir_3
!! case(35)
!   q_bar_dir(35) = q_bar_diff*heat_dir_3
!! case(36)
!   q_bar_dir(36) = q_bar_diff*heat_dir_3
!! case(37)
!   q_bar_dir(37) = -q_bar_diff*heat_dir_3
!! case(37)
!   q_bar_dir(38) = -q_bar_diff*heat_dir_3
!! case(37)
!   q_bar_dir(39) = -q_bar_diff*heat_dir_3

!q_bar_dir(1) = 0.0D0
!!
!q_bar_dir(2) = q_bar_diff(1)*temp_coeff
!!
!q_bar_dir(3) = q_bar_diff(2)*temp_coeff
!!
!q_bar_dir(4) = q_bar_diff(3)*temp_coeff
!!
!q_bar_dir(5) = -q_bar_diff(1)*temp_coeff
!!
!q_bar_dir(6) = -q_bar_diff(2)*temp_coeff
!!
!q_bar_dir(7) = -q_bar_diff(3)*temp_coeff
!!
!q_bar_dir(8) = temp_coeff*(q_bar_diff(1)+q_bar_diff(2)+q_bar_diff(3))
!!
!q_bar_dir(9) = temp_coeff*(-q_bar_diff(1)+q_bar_diff(2)+q_bar_diff(3))
!!
!q_bar_dir(10) = temp_coeff*(-q_bar_diff(1)-q_bar_diff(2)+q_bar_diff(3))
!!
!q_bar_dir(11) = temp_coeff*(q_bar_diff(1)-q_bar_diff(2)+q_bar_diff(3))
!!
!q_bar_dir(12) = temp_coeff*(q_bar_diff(1)+q_bar_diff(2)-q_bar_diff(3))
!!
!q_bar_dir(13) = temp_coeff*(-q_bar_diff(1)+q_bar_diff(2)-q_bar_diff(3))
!!
!q_bar_dir(14) = -temp_coeff*(q_bar_diff(1)+q_bar_diff(2)+q_bar_diff(3))
!!
!q_bar_dir(15) = temp_coeff*(q_bar_diff(1)-q_bar_diff(2)-q_bar_diff(3))
!!
!q_bar_dir(16) = 2.0D0*q_bar_diff(1)*temp_coeff
!!
!q_bar_dir(17) = 2.0D0*q_bar_diff(2)*temp_coeff
!!
!q_bar_dir(18) = 2.0D0*q_bar_diff(3)*temp_coeff
!!
!q_bar_dir(19) = -2.0D0*q_bar_diff(1)*temp_coeff
!!
!q_bar_dir(20) = -2.0D0*q_bar_diff(2)*temp_coeff
!!
!q_bar_dir(21) = -2.0D0*q_bar_diff(3)*temp_coeff
!!
!q_bar_dir(22) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(23) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(24) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(25) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(26) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(27) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(28) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(29) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(30) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(31) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(32) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(33) = 2.0D0*temp_coeff*(q_bar_diff(1)+q_bar_diff(2))
!!
!q_bar_dir(16) = 3.0D0*q_bar_diff(1)*temp_coeff
!!
!q_bar_dir(17) = 3.0D0*q_bar_diff(2)*temp_coeff
!!
!q_bar_dir(18) = 3.0D0*q_bar_diff(3)*temp_coeff
!!
!q_bar_dir(19) = -3.0D0*q_bar_diff(1)*temp_coeff
!!
!q_bar_dir(38) = -3.0D0*q_bar_diff(2)*temp_coeff
!!
!q_bar_dir(39) = -3.0D0*q_bar_diff(3)*temp_coeff


!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!!u_term = 1.0D0-loc_u
!!v_term = -loc_v
!!w_term = -loc_w
!!do i = 1,dir
!!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!!end do
!do i = 1,dir
!  xxx_diff = xxx_diff + (rcx(i)-loc_u)**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + (rcx(i)-loc_u)*(rcy(i)-loc_v)**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + (rcx(i)-loc_u)*(rcz(i)-loc_w)**2*loc_fidiff(i)
!end do
!q_bar_dir(2) = (xxx_diff*(1-3*loc_temp)-3*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!if (lvl == 3) then
!  write(*,*) 'q_bar_dir 2',xxx_diff,xyy_diff,xzz_diff,u_term,v_term,w_term,loc_temp,q_bar_dir(2),loc_fidiff
!end if
!! dir 3
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 1.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(3) = (yyy_diff*(1-3*loc_temp)-3*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!
!! dir 4
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(4) = (zzz_diff*(1-3*loc_temp)-3*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 5
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(5) = -(xxx_diff*(1-3*loc_temp)-3*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 6
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -1.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(6) = -(yyy_diff*(1-3*loc_temp)-3*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!! dir 7
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(7) = -(zzz_diff*(1-3*loc_temp)-3*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!! dir 8
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(8) = ((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 9
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(9) = ((1-3*loc_temp)*(-xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff-xyy_diff-xzz_diff+ &
!                   yyz_diff+yzz_diff)-6*xyz_diff+2*(xxy_diff-xyy_diff+xxz_diff-xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 10
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(10) = ((1-3*loc_temp)*(-xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff-xyy_diff-xzz_diff+ &
!                   yyz_diff-yzz_diff)+6*xyz_diff+2*(-xxy_diff-xyy_diff+xxz_diff-xzz_diff+&
!                   yyz_diff-yzz_diff))/temp_coeff
!
!
!! dir 11
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = 1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(11) = ((1-3*loc_temp)*(xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff-yzz_diff)-6*xyz_diff+2*(-xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff-yzz_diff))/temp_coeff
!
!
!! dir 12
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(12) = ((1-3*loc_temp)*(xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff+xyy_diff+xzz_diff- &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff-xxz_diff+xzz_diff-&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 13
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = 1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(13) = ((1-3*loc_temp)*(-xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff-xyy_diff-xzz_diff- &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff-xyy_diff-xxz_diff-xzz_diff-&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 14
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(14) = -((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!! dir 15
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 1.0D0-loc_u
!v_term = -1.0D0-loc_v
!w_term = -1.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xyz_diff = xyz_diff + u_term*v_term*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!  q_bar_dir(15) = ((1-3*loc_temp)*(xxx_diff-yyy_diff-zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+ &
!                   yyz_diff+yzz_diff)+6*xyz_diff+2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+&
!                   yyz_diff+yzz_diff))/temp_coeff
!
!
!! dir 16
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(16) = (xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 17
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(17) = (yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!
!! dir 18
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(18) = (zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 19
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(19) = -(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 20
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(20) = -(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!! dir 21
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(21) = -(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 22, +2,+2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = 2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(22) = ((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!
!! dir 23, +2,0,+2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(23) = ((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!
!! dir 24, 0,+2,+2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(24) = ((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!
!! dir 25, -2,+2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = 2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(25) = ((8-6*loc_temp)*(-xxx_diff+yyy_diff+xxy_diff-xyy_diff)+&
!                  16*(xxy_diff-xyy_diff)&
!                  -6*loc_temp*(-xzz_diff+yzz_diff))/temp_coeff
!
!! dir 26, -2,-2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(26) = -((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!
!! dir 27,+2,-2,0
!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -2.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(27) = ((8-6*loc_temp)*(xxx_diff-yyy_diff-xxy_diff+xyy_diff)+&
!                  16*(-xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff-yzz_diff))/temp_coeff
!
!! dir 28 +2,0,-2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = 2.0D0-loc_u
!v_term = -loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(28) = ((8-6*loc_temp)*(xxx_diff-zzz_diff-xxz_diff+xzz_diff)+&
!                 16*(-xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff-yyz_diff))/temp_coeff
!
!! dir 29, -2,0,-2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(29) = -((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!
!! dir 30, -2,0,+2
!xxx_diff = 0.0D0
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -2.0D0-loc_u
!v_term = -loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + v_term**2*u_term*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!end do
!q_bar_dir(30) = ((8-6*loc_temp)*(-xxx_diff+zzz_diff+xxz_diff-xzz_diff)+&
!                 16*(xxz_diff-xzz_diff)-&
!                 6*loc_temp*(-xyy_diff+yyz_diff))/temp_coeff
!
!! dir 31, 0,+2,-2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(31) = ((8-6*loc_temp)*(yyy_diff-zzz_diff-yyz_diff+yzz_diff)+&
!                 16*(-yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff-xxz_diff))/temp_coeff
!
!! dir 32, 0,-2,-2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -2.0D0-loc_v
!w_term = -2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(32) = -((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!
!! dir 33, 0,-2,+2
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -2.0D0-loc_v
!w_term = 2.0D0-loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + u_term**2*w_term*loc_fidiff(i)
!  xxy_diff = xxy_diff + u_term**2*v_term*loc_fidiff(i)
!  yyz_diff = yyz_diff + v_term**2*w_term*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(33) = ((8-6*loc_temp)*(-yyy_diff+zzz_diff+yyz_diff-yzz_diff)+&
!                 16*(yyz_diff-yzz_diff)-&
!                 6*loc_temp*(-xxy_diff+xxz_diff))/temp_coeff
!
!! dir 34
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = 3.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(34) = (xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!
!! dir 35
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = 3.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(35) = (yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!
!! dir 36
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = 3.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(36) = (zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!
!! dir 37
!xxx_diff = 0.0D0
!xyy_diff = 0.0D0
!xzz_diff = 0.0D0
!u_term = -3.0D0-loc_u
!v_term = -loc_v
!w_term = -loc_w
!do i = 1,dir
!  xxx_diff = xxx_diff + u_term**3*loc_fidiff(i)
!  xyy_diff = xyy_diff + u_term*v_term**2*loc_fidiff(i)
!  xzz_diff = xzz_diff + u_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(37) = -(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!! dir 38
!yyy_diff = 0.0D0
!xxy_diff = 0.0D0
!yzz_diff = 0.0D0
!u_term = -loc_u
!v_term = -3.0D0-loc_v
!w_term = -loc_w
!do i = 1,dir
!  yyy_diff = yyy_diff + v_term**3*loc_fidiff(i)
!  xxy_diff = xxy_diff + v_term*u_term**2*loc_fidiff(i)
!  yzz_diff = yzz_diff + v_term*w_term**2*loc_fidiff(i)
!end do
!q_bar_dir(38) = -(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!! dir 39
!zzz_diff = 0.0D0
!xxz_diff = 0.0D0
!yyz_diff = 0.0D0
!u_term = -loc_u
!v_term = -loc_v
!w_term = -3.0D0-loc_w
!do i = 1,dir
!  zzz_diff = zzz_diff + w_term**3*loc_fidiff(i)
!  xxz_diff = xxz_diff + w_term*u_term**2*loc_fidiff(i)
!  yyz_diff = yyz_diff + w_term*v_term**2*loc_fidiff(i)
!end do
!q_bar_dir(39) = -(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff







!xxx_diff = 0.0D0
!yyy_diff = 0.0D0
!zzz_diff = 0.0D0
!xxy_diff = 0.0D0
!xxz_diff = 0.0D0
!xyy_diff = 0.0D0
!xyz_diff = 0.0D0
!xzz_diff = 0.0D0
!yyz_diff = 0.0D0
!yzz_diff = 0.0D0
!
!do i = 1,dir
!  xxx_diff = xxx_diff + (rcx(i)-loc_u)**3*loc_fidiff(i)
!  yyy_diff = yyy_diff + (rcy(i)-loc_v)**3*loc_fidiff(i)
!  zzz_diff = zzz_diff + (rcz(i)-loc_w)**3*loc_fidiff(i)
!!
!  xxy_diff = xxy_diff + (rcx(i)-loc_u)**2 * (rcy(i)-loc_v) * loc_fidiff(i)
!!
!  xxz_diff = xxz_diff + (rcx(i)-loc_u)**2 * (rcz(i)-loc_w) * loc_fidiff(i)
!!
!  xyy_diff = xyy_diff + (rcy(i)-loc_v)**2 * (rcx(i)-loc_u) * loc_fidiff(i)
!!
!  xzz_diff = xzz_diff + (rcz(i)-loc_w)**2 * (rcx(i)-loc_u) * loc_fidiff(i)
!!
!  yyz_diff = yyz_diff + (rcy(i)-loc_v)**2 * (rcz(i)-loc_w) * loc_fidiff(i)
!!
!  yzz_diff = yzz_diff + (rcz(i)-loc_w)**2 * (rcy(i)-loc_v) * loc_fidiff(i)
!!
!  xyz_diff = xyz_diff + (rcx(i)-loc_u)*(rcy(i)-loc_v)* (rcz(i)-loc_w) * loc_fidiff(i)
!
!end do
!
!!Second shell
!! dir 2
!q_bar_dir(2) = init_epsilon(lvl)*(xxx_diff*(1-3*loc_temp)-3*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(3) = init_epsilon(lvl)*(yyy_diff*(1-3*loc_temp)-3*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(4) = init_epsilon(lvl)*(zzz_diff*(1-3*loc_temp)-3*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(5) = -init_epsilon(lvl)*(xxx_diff*(1-3*loc_temp)-3*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(6) = -init_epsilon(lvl)*(yyy_diff*(1-3*loc_temp)-3*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(7) = -init_epsilon(lvl)*(zzz_diff*(1-3*loc_temp)-3*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!!Third Shell
!q_bar_dir(8) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+yyz_diff+yzz_diff)+&
!                 6*xyz_diff+&
!                 2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(9) = init_epsilon(lvl)*((1-3*loc_temp)*(-xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff-xyy_diff-xzz_diff+yyz_diff+yzz_diff)-&
!                 6*xyz_diff+&
!                 2*(xxy_diff-xyy_diff+xxz_diff-xzz_diff+yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(10) = init_epsilon(lvl)*((1-3*loc_temp)*(-xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff-xyy_diff-xzz_diff+yyz_diff-yzz_diff)+&
!                 6*xyz_diff+&
!                 2*(-xxy_diff-xyy_diff+xxz_diff-xzz_diff+yyz_diff-yzz_diff))/temp_coeff
!q_bar_dir(11) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff-yyy_diff+zzz_diff-xxy_diff+xxz_diff+xyy_diff+xzz_diff+yyz_diff-yzz_diff)-&
!                  6*xyz_diff+&
!                  2*(-xxy_diff+xyy_diff+xxz_diff+xzz_diff+yyz_diff-yzz_diff))/temp_coeff
!q_bar_dir(12) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff+xyy_diff+xzz_diff-yyz_diff+yzz_diff)-&
!                  6*xyz_diff+&
!                  2*(xxy_diff+xyy_diff-xxz_diff+xzz_diff-yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(13) = init_epsilon(lvl)*((1-3*loc_temp)*(-xxx_diff+yyy_diff-zzz_diff+xxy_diff-xxz_diff-xyy_diff-xzz_diff-yyz_diff+yzz_diff)+&
!                  6*xyz_diff+&
!                  2*(xxy_diff-xyy_diff-xxz_diff-xzz_diff-yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(14) = -init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff+yyy_diff+zzz_diff+xxy_diff+xxz_diff+xyy_diff+xzz_diff+yyz_diff+yzz_diff)+&
!                   6*xyz_diff+&
!                   2*(xxy_diff+xyy_diff+xxz_diff+xzz_diff+yyz_diff+yzz_diff))/temp_coeff
!q_bar_dir(15) = init_epsilon(lvl)*((1-3*loc_temp)*(xxx_diff-yyy_diff-zzz_diff-xxy_diff-xxz_diff+xyy_diff+xzz_diff- yyz_diff-yzz_diff)+&
!                  6*xyz_diff + &
!                  2*(-xxy_diff+xyy_diff-xxz_diff+xzz_diff-yyz_diff-yzz_diff))/temp_coeff
!! Fourth shell
!q_bar_dir(16) = init_epsilon(lvl)*(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(17) = init_epsilon(lvl)*(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(18) = init_epsilon(lvl)*(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(19) = -init_epsilon(lvl)*(xxx_diff*(8-6*loc_temp)-6*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(20) = -init_epsilon(lvl)*(yyy_diff*(8-6*loc_temp)-6*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(21) = -init_epsilon(lvl)*(zzz_diff*(8-6*loc_temp)-6*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!! Fifth shell
!q_bar_dir(22) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!q_bar_dir(23) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!q_bar_dir(24) = init_epsilon(lvl)*((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!q_bar_dir(25) = init_epsilon(lvl)*((8-6*loc_temp)*(-xxx_diff+yyy_diff+xxy_diff-xyy_diff)+&
!                  16*(xxy_diff-xyy_diff)&
!                  -6*loc_temp*(-xzz_diff+yzz_diff))/temp_coeff
!q_bar_dir(26) = -init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+yyy_diff+xxy_diff+xyy_diff)+&
!                  16*(xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff+yzz_diff))/temp_coeff
!q_bar_dir(27) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff-yyy_diff-xxy_diff+xyy_diff)+&
!                  16*(-xxy_diff+xyy_diff)&
!                  -6*loc_temp*(xzz_diff-yzz_diff))/temp_coeff
!q_bar_dir(28) = init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff-zzz_diff-xxz_diff+xzz_diff)+&
!                 16*(-xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff-yyz_diff))/temp_coeff
!q_bar_dir(29) = -init_epsilon(lvl)*((8-6*loc_temp)*(xxx_diff+zzz_diff+xxz_diff+xzz_diff)+&
!                 16*(xxz_diff+xzz_diff)-&
!                 6*loc_temp*(xyy_diff+yyz_diff))/temp_coeff
!q_bar_dir(30) = init_epsilon(lvl)*((8-6*loc_temp)*(-xxx_diff+zzz_diff+xxz_diff-xzz_diff)+&
!                 16*(xxz_diff-xzz_diff)-&
!                 6*loc_temp*(-xyy_diff+yyz_diff))/temp_coeff
!q_bar_dir(31) = init_epsilon(lvl)*((8-6*loc_temp)*(yyy_diff-zzz_diff-yyz_diff+yzz_diff)+&
!                 16*(-yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff-xxz_diff))/temp_coeff
!q_bar_dir(32) = -init_epsilon(lvl)*((8-6*loc_temp)*(yyy_diff+zzz_diff+yyz_diff+yzz_diff)+&
!                 16*(yyz_diff+yzz_diff)-&
!                 6*loc_temp*(xxy_diff+xxz_diff))/temp_coeff
!q_bar_dir(33) = init_epsilon(lvl)*((8-6*loc_temp)*(-yyy_diff+zzz_diff+yyz_diff-yzz_diff)+&
!                 16*(yyz_diff-yzz_diff)-&
!                 6*loc_temp*(-xxy_diff+xxz_diff))/temp_coeff
!! Sixth shell
!q_bar_dir(34) = init_epsilon(lvl)*(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(35) = init_epsilon(lvl)*(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(36) = init_epsilon(lvl)*(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff
!q_bar_dir(37) = -init_epsilon(lvl)*(xxx_diff*(27-9*loc_temp)-9*loc_temp*(xyy_diff+xzz_diff))/temp_coeff
!q_bar_dir(38) = -init_epsilon(lvl)*(yyy_diff*(27-9*loc_temp)-9*loc_temp*(xxy_diff+yzz_diff))/temp_coeff
!q_bar_dir(39) = -init_epsilon(lvl)*(zzz_diff*(27-9*loc_temp)-9*loc_temp*(xxz_diff+yyz_diff))/temp_coeff



