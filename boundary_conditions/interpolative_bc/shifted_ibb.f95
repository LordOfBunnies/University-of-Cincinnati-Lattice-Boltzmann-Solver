SUBROUTINE shifted_ibb(loc_fout, loc_fi,loc_fi_opp,loc_fieq,loc_fieq_opp,&
  fout_adj,gout_adj,loc_gout, loc_gi,loc_gi_opp,loc_gieq, loc_gieq_opp,&
  u_tgt, v_tgt, w_tgt, temp_tgt, i, lvl, rest_state, &
  loc_rho,loc_u,loc_v,loc_w, loc_temp,state_pass,fine,a,b,c)
!
! This does the interpolative bounceback operations for the boundary condition
!
! Note, the fieq and gieq values are for the final direction (after the fout becomes
!   the fi value).  Thus, the opp values are actually NOT the ones of final interest.
!
!
!
!
! Calls: quick_sutherlands,quick_therm_cond
! Called by: streaming
!
USE precise
use constants
use quick_calcs
USE grid_data
USE amrex_base_module
use linkwise
use freestream_values
!use amrex_amr_module
use amr_info_holder, only: elbm_epsilon,init_epsilon,timestep,startup
use amr_processes, only: self,beta_1_calculator,beta_2_calculator,amr_max_lvl,epsilon_calculator
use derivatives
IMPLICIT NONE
INTEGER,INTENT(IN) :: i,lvl,rest_state(a:a,b:b,c:c,1:dir),a,b,c,fine,state_pass(7)
real(kind=dp),intent(in) :: loc_u(7),loc_v(7),loc_w(7),loc_temp(7)
REAL(KIND=dp) :: loc_fout,loc_fi,loc_gout,loc_gi,loc_fieq,loc_gieq,loc_fieq_opp,loc_gieq_opp
real(kind=dp) :: dir_state,inc_fout,inc_gout,loc_fi_opp,loc_gi_opp,fout_adj,gout_adj
!REAL(KIND=dp),INTENT(IN) :: fout_opp
real(kind=dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
real(kind=dp) :: dTdx,dTdy,dTdz
real(kind=dp),intent(in) :: temp_tgt,u_tgt,v_tgt,w_tgt,loc_rho
real(kind=dp) :: loc_mu,loc_kappa,p_comp,q_comp
real(kind=dp) :: k_one,g_alpha,r_ab,loc_epsilon
real(kind=dp) :: fi_neq,gi_neq,cross_der,pri_der,temp_der_sum
real(kind=dp) :: k_fi,k_gi,beta_1,beta_2,cv_const,cv_const_2
real(kind=dp) :: p_term_1,p_term_2,p_term_3,q_term_1,q_term_2,q_term_3
real(kind=dp) :: pab(27),qabc(27),rab(27)
logical :: upwind,downwind

upwind = .true.
downwind = .false.

cv_const = 1/const_vol
cv_const_2 = 2/const_vol
k_fi = -(loc_rho*temp_tgt)/2.0D0
k_gi = -(loc_rho*temp_tgt*(2.0D0*const_vol-dimensions))/2.0D0

loc_epsilon = epsilon_calculator(u_tgt,v_tgt,w_tgt,temp_tgt,timestep(amr_max_lvl),a,b,c,lvl)

!write(*,*) loc_temp,temp_tgt,loc_rho,lvl
call quick_sutherlands(loc_mu,temp_tgt,.false.)
call quick_therm_cond(loc_kappa,loc_rho,temp_tgt)

!beta_1 = 1.0D0/((2.0D0*loc_mu)/(loc_rho*temp_tgt)+1.0D0)
!beta_2 = 1.0D0/((2.0D0*loc_kappa)/(loc_rho*const_press*temp_tgt)+1.0D0)

beta_1 = beta_1_calculator(loc_mu,loc_rho,temp_tgt,lvl,timestep(amr_max_lvl))
beta_2 = beta_2_calculator(loc_kappa,loc_rho,temp_tgt,lvl,timestep(amr_max_lvl))
!if (timestep(fine) < 10) then
!  beta_1 = 0.90D0 + 0.01D0*timestep(fine)
!end if

dir_state = (abs(real(rest_state(a,b,c,i),kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)

inc_fout = loc_fout
inc_gout = loc_gout

!
! for card, x = 1, y = 2, z = 3
!
! Fix it later if there's ever time
  if (primary_flow_dir == 1) then
    dudx = first_order_shifted_vel_cond_unknown(loc_u(1),u_tgt,loc_u(5),loc_u(2),i,1,state_pass(1),&
             state_pass(5),state_pass(2))
    dudy = first_order_shifted_vel_cond_unknown(loc_u(1),u_tgt,loc_u(6),loc_u(3),i,2,state_pass(1),&
             state_pass(6),state_pass(3))
    dudz = first_order_shifted_vel_cond_unknown(loc_u(1),u_tgt,loc_u(7),loc_u(4),i,3,state_pass(1),&
             state_pass(7),state_pass(4))
    dvdx = first_order_shifted_vel_cond_unknown(loc_v(1),v_tgt,loc_v(5),loc_v(2),i,1,state_pass(1),&
             state_pass(5),state_pass(2))
    dvdy = first_order_shifted_vel_cond_unknown(loc_v(1),v_tgt,loc_v(6),loc_v(3),i,2,state_pass(1),&
             state_pass(6),state_pass(3))
    dvdz = first_order_shifted_vel_cond_unknown(loc_v(1),v_tgt,loc_v(7),loc_v(4),i,3,state_pass(1),&
             state_pass(7),state_pass(4))
    dwdx = first_order_shifted_vel_cond_unknown(loc_w(1),w_tgt,loc_w(5),loc_w(2),i,1,state_pass(1),&
             state_pass(5),state_pass(2))
    dwdy = first_order_shifted_vel_cond_unknown(loc_w(1),w_tgt,loc_w(6),loc_w(3),i,2,state_pass(1),&
             state_pass(6),state_pass(3))
    dwdz = first_order_shifted_vel_cond_unknown(loc_w(1),w_tgt,loc_w(7),loc_w(4),i,3,state_pass(1),&
             state_pass(7),state_pass(4))
    dTdx = first_order_shifted_temp_cond_unknown(loc_temp(1),temp_tgt,loc_temp(5),loc_temp(2),i,1,state_pass(1),&
             state_pass(5),state_pass(2))
    dTdy = first_order_shifted_temp_cond_unknown(loc_temp(1),temp_tgt,loc_temp(6),loc_temp(3),i,2,state_pass(1),&
             state_pass(6),state_pass(3))
    dTdz = first_order_shifted_temp_cond_unknown(loc_temp(1),temp_tgt,loc_temp(7),loc_temp(4),i,3,state_pass(1),&
             state_pass(7),state_pass(4))
  else

  end if
!
    pab(1) = 2*dudx - dudx/const_vol !xxx
    pab(2) = 2*dudx - dvdy/const_vol !xxy
    pab(3) = 2*dudx - dwdz/const_vol !xxz
    pab(4) = dudy+dvdx !xyx
    pab(5) = dudy+dvdx !xyy
    pab(6) = dudy+dvdx !xyz
    pab(7) = dudz+dwdx !xzx
    pab(8) = dudz+dwdx !xzy
    pab(9) = dudz+dwdx !xzz
    pab(10) = dudy+dvdx
    pab(11) = dudy+dvdx
    pab(12) = dudy+dvdx
    pab(13) = 2*dvdy - dudx/const_vol
    pab(14) = 2*dvdy - dvdy/const_vol
    pab(15) = 2*dvdy - dwdz/const_vol
    pab(16) = dvdz+dwdy
    pab(17) = dvdz+dwdy
    pab(18) = dvdz+dwdy
    pab(19) = dudz+dwdx
    pab(20) = dudz+dwdx
    pab(21) = dudz+dwdx
    pab(22) = dvdz+dwdy
    pab(23) = dvdz+dwdy
    pab(24) = dvdz+dwdy
    pab(25) = 2*dwdz - dudx/const_vol
    pab(26) = 2*dwdz - dvdy/const_vol
    pab(27) = 2*dwdz - dwdz/const_vol

    pab = pab/beta_1

    qabc(1) = 3*(dTdx/beta_2 + u_tgt*pab(1))    !xxx
    qabc(2) = (dTdy/beta_2 + v_tgt*pab(2) + 2*u_tgt*pab(4)) !xxy
    qabc(3) = (dTdz/beta_2 + w_tgt*pab(3) + 2*u_tgt*pab(7)) !xxz
    qabc(4) = (dTdy/beta_2 + 2*u_tgt*pab(4) + v_tgt*pab(1)) !xyx
    qabc(5) = (dTdx/beta_2 + u_tgt*pab(13) + 2*v_tgt*pab(4))!xyy
    qabc(6) = (u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!xyz
    qabc(7) = (dTdz/beta_2 + w_tgt*pab(3) + 2*u_tgt*pab(7))!xzx
    qabc(8) = (u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!xzy
    qabc(9) = (dTdx/beta_2 + u_tgt*pab(25) + 2*w_tgt*pab(7))!xzz
    qabc(10) = (dTdy/beta_2 + v_tgt*pab(2) + 2*u_tgt*pab(4))!yxx
    qabc(11) = (dTdx/beta_2 + 2*v_tgt*pab(4) + u_tgt*pab(13))!yxy
    qabc(12) = (u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4))!yxz
    qabc(13) = (dTdx/beta_2 + u_tgt*pab(13) + 2*v_tgt*pab(4))!yyx
    qabc(14) = 3*(dTdy/beta_2 + v_tgt*pab(14))!yyy
    qabc(15) = (dTdz/beta_2 + w_tgt*pab(15) + 2*v_tgt*pab(22))!yyz
    qabc(16) = (u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!yzx
    qabc(17) = (dTdz/beta_2 + w_tgt*pab(15) + 2*v_tgt*pab(22))!yzy
    qabc(18) = (dTdy/beta_2 + v_tgt*pab(26) + 2*w_tgt*pab(22))!yzz
    qabc(19) = (dTdx/beta_2 + w_tgt*pab(3) + 2*u_tgt*pab(7))!zxx
    qabc(20) = (u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!zxy
    qabc(21) = (dTdx/beta_2 + u_tgt*pab(25) + 2*w_tgt*pab(7))!zxz
    qabc(22) = (u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!zyx
    qabc(23) = (dTdz/beta_2 + w_tgt*pab(15) + 2*v_tgt*pab(22))!zyy
    qabc(24) = (dTdy/beta_2 + v_tgt*pab(26) + 2*w_tgt*pab(22))!zyz
    qabc(25) = (dTdx/beta_2 + u_tgt*pab(25) + 2*w_tgt*pab(7))!zzx
    qabc(26) = (dTdy/beta_2 + v_tgt*pab(26) + 2*w_tgt*pab(22))!zzy
    qabc(27) = 3*(dTdz/beta_2 + w_tgt*pab(27))!zzz

!    fi_neq = sum(pab*(p_coeff(:,i)-temp_tgt*p_kroen))/(2*temp_tgt**2) + &
!             sum(qabc*(yi(:,i)-3*temp_tgt*kroen_val(:,i)))/(6*temp_tgt**3)
p_comp = (pab(1)*(cx(i)**2-temp_tgt) + pab(2)*(cx(i)**2-temp_tgt) + pab(3)*(cx(i)**2-temp_tgt) + &
  pab(4)*cx(i)*cy(i) + pab(5)*cx(i)*cy(i) + pab(6)*cx(i)*cy(i) + pab(7)*cx(i)*cz(i) + &
  pab(8)*cx(i)*cz(i) + pab(9)*cx(i)*cz(i) + pab(10)*cy(i)*cx(i) + pab(11)*cy(i)*cx(i) + &
  pab(12)*cy(i)*cx(i) + pab(13)*(cy(i)**2-temp_tgt) + pab(14)*(cy(i)**2-temp_tgt) + &
  pab(15)*(cy(i)**2-temp_tgt) + pab(16)*cy(i)*cz(i) + pab(17)*cy(i)*cz(i) + pab(18)*cy(i)*cz(i) + &
  pab(19)*cz(i)*cx(i) + pab(20)*cz(i)*cx(i) + pab(21)*cz(i)*cx(i) + pab(22)*cz(i)*cy(i) + &
  pab(23)*cz(i)*cy(i) + pab(24)*cz(i)*cy(i) + pab(25)*(cz(i)**2-temp_tgt) + pab(25)*(cz(i)**2-temp_tgt) +&
  pab(26)*(cz(i)**2-temp_tgt) + pab(27)*(cz(i)**2-temp_tgt))/(2*temp_tgt**2)
!
q_comp = (qabc(1)*(cx(i)**3-3*cx(i)*temp_tgt) + qabc(2)*(cx(i)**2*cy(i)-3*cy(i)*temp_tgt) + &
  qabc(3)*(cx(i)**2*cz(i)-3*cz(i)*temp_tgt) + qabc(4)*cx(i)**2*cy(i) + qabc(5)*cx(i)*cy(i)**2 + &
  qabc(6)*cx(i)*cy(i)*cz(i) + qabc(7)*cx(i)**2*cz(i) + qabc(8)*cx(i)*cy(i)*cz(i) + qabc(9)*cx(i)*cz(i)**2 + &
  qabc(10)*cx(i)**2*cy(i) + qabc(11)*cy(i)**2*cx(i) + qabc(12)*cx(i)*cy(i)*cz(i) + &
  qabc(13)*(cy(i)**2*cx(i)-3*cx(i)*temp_tgt) + qabc(14)*(cy(i)**3-3*cy(i)*temp_tgt)+ &
  qabc(15)*(cy(i)**2*cz(i)-3*cz(i)*temp_tgt) + qabc(16)*cx(i)*cy(i)*cz(i) + qabc(17)*cy(i)**2*cz(i) + &
  qabc(18)*cz(i)**2*cy(i) + qabc(19)*cx(i)**2*cz(i) + qabc(20)*cx(i)*cy(i)*cz(i) + qabc(21)*cz(i)**2*cx(i) +&
  qabc(22)*cx(i)*cy(i)*cz(i) + qabc(23)*cz(i)*cy(i)**2 + qabc(24)*cz(i)**2*cy(i) + &
  qabc(25)*(cz(i)**2*cx(i)-3*cx(i)*temp_tgt) + qabc(26)*(cz(i)**2*cy(i)-3*cy(i)*temp_tgt) + &
   qabc(27)*(cz(i)**3-3*cz(i)*temp_tgt))/(6*temp_tgt**3)


!    loc_fout = loc_fieq + loc_epsilon*fi_neq
!    loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
    loc_fout = loc_fieq_opp + loc_epsilon*k_fi*(p_comp + q_comp)
!
!
!
    k_one = 3*(dTdx+dTdy+dTdz)/(beta_1*const_vol) !as this is used multiple time, multiply later
    g_alpha = cx(i)/temp_tgt*(9*dTdx/beta_2 + u_tgt*k_one) + &
      cy(i)/temp_tgt*(9*dTdy/beta_2 + v_tgt*k_one) + &
      cz(i)/temp_tgt*(9*dTdz/beta_2 + w_tgt*k_one)

    rab(1) = 2*temp_tgt*dudx + 2*u_tgt*dTdx - dudx/const_vol*(2*temp_tgt + u_tgt**2)
    rab(2) = 2*temp_tgt*dudx + 2*u_tgt*dTdx - dvdy/const_vol*(2*temp_tgt + u_tgt**2)
    rab(3) = 2*temp_tgt*dudx + 2*u_tgt*dTdx - dwdz/const_vol*(2*temp_tgt + u_tgt**2)
    rab(4) = temp_tgt*(dudy+dvdx) + u_tgt*dTdy+v_tgt*dTdx - dudx/const_vol*u_tgt*v_tgt
    rab(5) = temp_tgt*(dudy+dvdx) + u_tgt*dTdy+v_tgt*dTdx - dvdy/const_vol*u_tgt*v_tgt
    rab(6) = temp_tgt*(dudy+dvdx) + u_tgt*dTdy+v_tgt*dTdx - dwdz/const_vol*u_tgt*v_tgt
    rab(7) = temp_tgt*(dudz+dwdx) + u_tgt*dTdz+w_tgt*dTdx - dudx/const_vol*u_tgt*w_tgt
    rab(8) = temp_tgt*(dudz+dwdx) + u_tgt*dTdz+w_tgt*dTdx - dvdy/const_vol*u_tgt*w_tgt
    rab(9) = temp_tgt*(dudz+dwdx) + u_tgt*dTdz+w_tgt*dTdx - dwdz/const_vol*u_tgt*w_tgt
    rab(10) = rab(4)
    rab(11) = rab(5)
    rab(12) = rab(6)
    rab(13) = 2*temp_tgt*dvdy + 2*v_tgt*dTdy - dudx/const_vol*(2*temp_tgt + v_tgt**2)
    rab(14) = 2*temp_tgt*dvdy + 2*v_tgt*dTdy - dvdy/const_vol*(2*temp_tgt + v_tgt**2)
    rab(15) = 2*temp_tgt*dvdy + 2*v_tgt*dTdy - dwdz/const_vol*(2*temp_tgt + v_tgt**2)
    rab(16) = temp_tgt*(dvdz+dwdy) + w_tgt*dTdy+v_tgt*dTdz - dudx/const_vol*v_tgt*w_tgt
    rab(17) = temp_tgt*(dvdz+dwdy) + w_tgt*dTdy+v_tgt*dTdz - dvdy/const_vol*v_tgt*w_tgt
    rab(18) = temp_tgt*(dvdz+dwdy) + w_tgt*dTdy+v_tgt*dTdz - dwdz/const_vol*v_tgt*w_tgt
    rab(19) = rab(7)
    rab(20) = rab(8)
    rab(21) = rab(9)
    rab(22) = rab(13)
    rab(23) = rab(14)
    rab(24) = rab(15)
    rab(25) = 2*temp_tgt*dwdz + 2*w_tgt*dTdz - dudx/const_vol*(2*temp_tgt + w_tgt**2)
    rab(26) = 2*temp_tgt*dwdz + 2*w_tgt*dTdz - dvdy/const_vol*(2*temp_tgt + w_tgt**2)
    rab(27) = 2*temp_tgt*dwdz + 2*w_tgt*dTdz - dwdz/const_vol*(2*temp_tgt + w_tgt**2)
!
gi_neq = k_gi*(3*k_one + g_alpha + (rab(1)*(cx(i)**2-temp_tgt) + rab(2)*(cx(i)**2-temp_tgt) + &
  rab(3)*(cx(i)**2-temp_tgt) + rab(4)*cx(i)*cy(i) + rab(5)*cx(i)*cy(i) + rab(6)*cx(i)*cy(i) + &
  rab(7)*cx(i)*cz(i) + rab(8)*cx(i)*cz(i) + rab(9)*cx(i)*cz(i) + rab(10)*cy(i)*cx(i) + &
  rab(11)*cy(i)*cx(i) + rab(12)*cy(i)*cx(i) + rab(13)*(cy(i)**2-temp_tgt) + rab(14)*(cy(i)**2-temp_tgt) + &
  rab(15)*(cy(i)**2-temp_tgt) + rab(16)*cy(i)*cz(i) + rab(17)*cy(i)*cz(i) + rab(18)*cy(i)*cz(i) + &
  rab(19)*cz(i)*cx(i) + rab(20)*cz(i)*cx(i) + rab(21)*cz(i)*cx(i) + rab(22)*cz(i)*cy(i) + &
  rab(23)*cz(i)*cy(i) +rab(24)*cz(i)*cy(i) + rab(25)*(cz(i)**2-temp_tgt) + rab(26)*(cz(i)**2-temp_tgt) +&
  rab(27)*(cz(i)**2-temp_tgt))/(2*beta_1*temp_tgt))

  loc_gout = loc_gieq_opp + loc_epsilon*gi_neq

!  write(*,*) 'outgoing IBB values',loc_fout,loc_fi,loc_fieq,loc_fi_opp,loc_fieq_opp,&
!    loc_gout,loc_gi,loc_gieq,loc_gi_opp,loc_gieq_opp,i,self
!  write(*,*) 'target values',loc_rho,u_tgt,v_tgt,w_tgt,temp_tgt,i,rest_state(a,b,c,i),dir_state,lvl
!  write(*,*) 'ibb states',state_pass,rest_state(a,b,c,1:dir),a,b,c,i,self
!  write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!  write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz
!  write(*,*) 'equilibrium values, constants',loc_fieq,loc_gieq,beta_1,beta_2,k_fi,k_gi
!  write(*,*) 'density values',loc_rho
!  write(*,*) 'local u values',loc_u
!  write(*,*) 'local v values',loc_v
!  write(*,*) 'local w values',loc_w
!  write(*,*) 'local temperature',loc_temp
!  write(*,*) 'local fi information',pri_der,p_comp,q_comp,fi_neq*loc_epsilon,loc_epsilon,loc_fout
!  write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*loc_epsilon,loc_epsilon,loc_gout
!
!    write(*,*) 'pab',pab
!    write(*,*) 'qabc',qabc
!    write(*,*) 'k_one',k_one
!    write(*,*) 'g_alpha',g_alpha
!    write(*,*) 'rab',rab

!    gi_neq = k_gi*(3*k_one + g_alpha + sum(rab*p_coeff(:,i) - &
!               temp_tgt*p_kroen)/(2*beta_1*temp_tgt**2))
!    loc_gout = loc_gieq + loc_epsilon*gi_neq

    if (loc_fout < 0.0D0) then
      loc_fout = abs(loc_fout)
    end if
    if (loc_gout < 0.0D0) then
      loc_gout = abs(loc_gout)
    end if

end subroutine
