SUBROUTINE interpolative_bounceback_bc(loc_fout, loc_fi,loc_fi_opp,loc_fieq,loc_fieq_opp,&
  fout_adj,gout_adj,loc_gout, loc_gi,loc_gi_opp,loc_gieq, loc_gieq_opp,&
  u_tgt, v_tgt, w_tgt, temp_tgt, i, lvl, rest_state, &
  loc_rho,loc_u,loc_v,loc_w, loc_temp,fine,a,b,c)
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
INTEGER,INTENT(IN) :: i,lvl,rest_state(a:a,b:b,c:c,1:dir),a,b,c,fine
real(kind=dp),intent(in) :: loc_u(7),loc_v(7),loc_w(7),loc_temp(7)
REAL(KIND=dp) :: loc_fout,loc_fi,loc_gout,loc_gi,loc_fieq,loc_gieq,loc_fieq_opp,loc_gieq_opp
real(kind=dp) :: dir_state,inc_fout,inc_gout,loc_fi_opp,loc_gi_opp,fout_adj,gout_adj
!REAL(KIND=dp),INTENT(IN) :: fout_opp
real(kind=dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
real(kind=dp) :: dTdx,dTdy,dTdz
real(kind=dp),intent(in) :: temp_tgt,u_tgt,v_tgt,w_tgt,loc_rho
real(kind=dp) :: loc_mu,loc_kappa,loc_epsilon
real(kind=dp) :: k_one,g_alpha,r_ab
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

loc_epsilon = epsilon_calculator(u_tgt,v_tgt,w_tgt,temp_tgt,timestep(amr_max_lvl),a,b,c,lvl)

!write(*,*) 'In IBB bc',k_fi,k_gi,beta_1,beta_2
!write(*,*) 'IBB inputs',
!if (a == 1 .and. b == 39 .and. c == 20) then
!  write(*,*) 'incoming ibb values near wall', loc_fout, loc_fi,loc_fi_opp,loc_fieq,loc_fieq_opp,&
!  fout_adj,gout_adj,loc_gout, loc_gi,loc_gi_opp,loc_gieq, loc_gieq_opp,i
!  write(*,*) 'incoming targets',loc_rho,u_tgt, v_tgt, w_tgt, temp_tgt,loc_u,loc_v,loc_w, loc_temp,a,b,c
!  write(*,*) 'incoming states',rest_state(a,b,c,1:dir),a,b,c
!end if
!
!write(*,*) 'incoming IBB values',loc_fieq,loc_gieq,rest_state,dir_state,&
!  loc_rho,u_tgt,v_tgt,w_tgt,temp_tgt,i,self
!
if (dir == 19) then
select case(i)
  case(2)
! 2 -  (+1,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)

  pri_der = dudx+dvdy+dwdz

! f_missing
  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) + &
    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) - &
    3/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

! g_missing
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(3)
! 3 -  (0,+1,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  pri_der = dudx+dvdy+dwdz
! f_missing
  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) + &
    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - &
    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) - &
    3/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  loc_fout = loc_fieq + loc_epsilon*fi_neq

!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(4)
! 4 -  (0,0,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f_missing
  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(dudx+dvdy+dwdz))) + &
    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(dudx+dvdy+dwdz)) - &
    3/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  loc_fout = loc_fieq + loc_epsilon*fi_neq


!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(5)
! 5 -  (-1,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  pri_der = dudx+dvdy+dwdz
! f_missing
  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
    3/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  loc_fout = loc_fieq + loc_epsilon*fi_neq


!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(6)
! 6 -  (0,-1,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  pri_der = dudx+dvdy+dwdz
! f_missing
  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) - &
    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - &
    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
    3/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  loc_fout = loc_fieq + loc_epsilon*fi_neq

!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(7)
! 7 -  (0,0,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  pri_der = dudx+dvdy+dwdz

!
  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(dudx+dvdy+dwdz))) - &
    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(dudx+dvdy+dwdz)) + &
    3/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  loc_fout = loc_fieq + loc_epsilon*fi_neq


!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(8)
! 8 -  (+1,+1,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dudy+dvdx)
  p_term_3 = 6*dwdz-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdx+dTdy)+(u_tgt*(6*dudx+2*(dvdx+dudy+dvdy)-4*cv_const*dudx)+&
    v_tgt*(6*dvdy+2*(dvdx+dudy+dudx)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdx+dTdy)+(u_tgt*(4*(dudy+dvdx+dvdy)-cv_const_2*dudx)+&
    v_tgt*(4*(dudy+dvdx+dudx)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdx+dTdy)/beta_2 + (2.0D0*w_tgt*(dudz+dwdx+dvdz+dwdy) + &
    u_tgt*(2*dwdz-cv_const*dudx) + v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
!
  case(9)
! 9 -  (0,+1,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dwdy+dvdz)
  p_term_3 = 6*dudx-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(-dTdz+dTdy)+(-w_tgt*(6*dwdz+2*(-dvdz-dwdy+dvdy)-4*cv_const*dwdz)+&
    v_tgt*(6*dvdy+2*(-dvdz-dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(-dTdz+dTdy)+(-w_tgt*(4*(-dwdy-dvdz+dvdy)-cv_const_2*dwdz)+&
    v_tgt*(4*(-dwdy-dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((-dTdz+dTdy)/beta_2 + (2.0D0*u_tgt*(-dudz-dwdx+dvdx+dudy) - &
    w_tgt*(2*dudx-cv_const*dwdz) + v_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(const_vol*beta_1)
  g_alpha = 9*(dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(v_tgt-w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
 !
  case(10)
! 10 - (-1,+1,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dudy+dvdx)
  p_term_3 = 6*dwdz-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(-dTdx+dTdy)+(-u_tgt*(6*dudx+2*(-dvdx-dudy+dvdy)-4*cv_const*dudx)+&
    v_tgt*(6*dvdy+2*(-dvdx-dudy+dudx)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(-dTdx+dTdy)+(-u_tgt*(4*(-dudy-dvdx+dvdy)-cv_const_2*dudx)+&
    v_tgt*(4*(-dudy-dvdx+dudx)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((-dTdx+dTdy)/beta_2 + (2.0D0*w_tgt*(-dudz-dwdx+dvdz+dwdy) - &
    u_tgt*(2*dwdz-cv_const*dudx) + v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdy-dTdx)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(v_tgt-u_tgt)
  r_ab = (1-temp_tgt)/(4*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(11)
! 11 - (0,+1,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dwdy+dvdz)
  p_term_3 = 6*dudx-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdz+dTdy)+(w_tgt*(6*dwdz+2*(dvdz+dwdy+dvdy)-4*cv_const*dwdz)+&
    v_tgt*(6*dvdy+2*(dvdz+dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdz+dTdy)+(w_tgt*(4*(dwdy+dvdz+dvdy)-cv_const_2*dwdz)+&
    v_tgt*(4*(dwdy+dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdz+dTdy)/beta_2 + (2.0D0*u_tgt*(dudz+dwdx+dvdx+dudy) + &
    w_tgt*(2*dudx-cv_const*dwdz) + v_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt+v_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(12)
! 12 - (+1,0,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = 6*(dudz+dwdx)
  p_term_3 = 6*dvdy-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdx+dTdz)+(u_tgt*(6*dudx+2*(dwdx+dudz+dwdz)-4*cv_const*dudx)+&
    w_tgt*(6*dwdz+2*(dwdx+dudz+dudx)-4*cv_const*dwdz))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdx+dTdz)+(u_tgt*(4*(dudz+dwdx+dwdz)-cv_const_2*dudx)+&
    w_tgt*(4*(dudz+dwdx+dudx)-cv_const_2*dwdz))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx+dvdz+dwdy) + &
    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(13)
! 13 - (+1,0,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = 6*(dudz+dwdx)
  p_term_3 = 6*dvdy-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdx-dTdz)+(u_tgt*(6*dudx+2*(-dwdx-dudz+dwdz)-4*cv_const*dudx) - &
    w_tgt*(6*dwdz+2*(-dwdx-dudz+dudx)-4*cv_const*dwdz))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdx-dTdz)+(u_tgt*(4*(-dudz-dwdx+dwdz)-cv_const_2*dudx) - &
    w_tgt*(4*(-dudz-dwdx+dudx)-cv_const_2*dwdz))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdx-dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx-dvdz-dwdy) + &
    u_tgt*(2*dvdy-cv_const*dudx) - w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(14)
! 14 - (-1,0,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)

!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = 6*(dudz+dwdx)
  p_term_3 = 6*dvdy-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdx+dTdz)+(u_tgt*(6*dudx+2*(dwdx+dudz+dwdz)-4*cv_const*dudx)+&
    w_tgt*(6*dwdz+2*(dwdx+dudz+dudx)-4*cv_const*dwdz))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdx+dTdz)+(u_tgt*(4*(dudz+dwdx+dwdz)-cv_const_2*dudx)+&
    w_tgt*(4*(dudz+dwdx+dudx)-cv_const_2*dwdz))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx+dvdz+dwdy) + &
    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(15)
! 15 - (-1,0,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = 6*(dudz+dwdx)
  p_term_3 = 6*dvdy-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(-dTdx+dTdz)+(-u_tgt*(6*dudx+2*(-dwdx-dudz+dwdz)-4*cv_const*dudx)+&
    w_tgt*(6*dwdz+2*(-dwdx-dudz+dudx)-4*cv_const*dwdz))/beta_1)
  q_term_2 = (2.0D0/beta_2*(-dTdx+dTdz)+(-u_tgt*(4*(-dudz-dwdx+dwdz)-cv_const_2*dudx)+&
    w_tgt*(4*(-dudz-dwdx+dudx)-cv_const_2*dwdz))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((-dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(-dudy-dvdx+dvdz+dwdy) - &
    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdz-dTdx)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt-u_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(16)
! 16 - (+1,-1,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dudy+dvdx)
  p_term_3 = 6*dwdz-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdx-dTdy)+(u_tgt*(6*dudx+2*(-dvdx-dudy+dvdy)-4*cv_const*dudx)-&
    v_tgt*(6*dvdy+2*(-dvdx-dudy+dudx)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdx-dTdy)+(u_tgt*(4*(-dudy-dvdx+dvdy)-cv_const_2*dudx)-&
    v_tgt*(4*(-dudy-dvdx+dudx)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdx-dTdy)/beta_2 + (2.0D0*w_tgt*(dudz+dwdx-dvdz-dwdy) + &
    u_tgt*(2*dwdz-cv_const*dudx) - v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx-dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt)
  r_ab = (1-temp_tgt)/(4*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(17)
! 17 - (0,-1,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dwdy+dvdz)
  p_term_3 = 6*dudx-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdz+dTdy)+(w_tgt*(6*dwdz+2*(dvdz+dwdy+dvdy)-4*cv_const*dwdz)+&
    v_tgt*(6*dvdy+2*(dvdz+dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdz+dTdy)+(w_tgt*(4*(dwdy+dvdz+dvdy)-cv_const_2*dwdz)+&
    v_tgt*(4*(dwdy+dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdz+dTdy)/beta_2 + (2.0D0*u_tgt*(dudz+dwdx+dvdx+dudy) + &
    w_tgt*(2*dudx-cv_const*dwdz) + v_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt+v_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(18)
! 18 - (-1,-1,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dudy+dvdx)
  p_term_3 = 6*dwdz-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdx+dTdy)+(u_tgt*(6*dudx+2*(dvdx+dudy+dvdy)-4*cv_const*dudx)+&
    v_tgt*(6*dvdy+2*(dvdx+dudy+dudx)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdx+dTdy)+(u_tgt*(4*(dudy+dvdx+dvdy)-cv_const_2*dudx)+&
    v_tgt*(4*(dudy+dvdx+dudx)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdx+dTdy)/beta_2 + (2.0D0*w_tgt*(dudz+dwdx+dvdz+dwdy) + &
    u_tgt*(2*dwdz-cv_const*dudx) + v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt)
  r_ab = (1-temp_tgt)/(4*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)
!
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(19)
! 19 - (0,-1,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
  p_term_1 = (1.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
  p_term_2 = 6*(dwdy+dvdz)
  p_term_3 = 6*dudx-cv_const*pri_der
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(4*(dTdz-dTdy)+(w_tgt*(6*dwdz+2*(-dvdz-dwdy+dvdy)-4*cv_const*dwdz) - &
    v_tgt*(6*dvdy+2*(-dvdz-dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
  q_term_2 = (2.0D0/beta_2*(dTdz-dTdy)+(w_tgt*(4*(-dwdy-dvdz+dvdy)-cv_const_2*dwdz) - &
    v_tgt*(4*(-dwdy-dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
  q_term_3 = (-3.0D0*temp_tgt)*((dTdz-dTdy)/beta_2 + (2.0D0*u_tgt*(dudz+dwdx-dvdx-dudy) + &
    w_tgt*(2*dudx-cv_const*dwdz) - v_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdz-dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt-v_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dudx+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
end select

else if (dir == 39) then

if (.not. shifted) then
select case(i)
  case(2)
! 2 - (+1,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-(pri_der)/const_vol)) + & ! term 1
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*dudx/const_vol)) - & ! term 2
!    1/(2*beta_1*temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*(pri_der)) - & ! term 3
!    3/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & ! term 4
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!
  p_term_1 = (1.0D0-temp_tgt)*(6.0D0*dudx-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dwdz)-2.0D0/const_vol*pri_der)
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(3.0D0/beta_2*dTdx+u_tgt/beta_1*dudx*(6.0D0-3.0D0/const_vol))
  q_term_2 = -3.0D0*temp_tgt*(2.0D0/beta_2*dTdx + 2.0D0/beta_1*(u_tgt*(dwdz+dvdy-dudx/const_vol) + &
    v_tgt*(dudy+dvdx) + w_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

! g_missing
  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
!  write(*,*) 'boom boom'

  case(3)
! 3 - (0,+1,0)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(pri_der))) + &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy + v_tgt/beta_1*(6*dvdy - 3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*(dudx+dwdz)-cv_const_2*(pri_der)) - &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  p_term_1 = (1.0D0-temp_tgt)*(6.0D0*dvdy - 1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dudx+dwdz) - 2.0D0/const_vol*pri_der)
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(3.0D0/beta_2*dTdy + v_tgt/beta_1*dvdy*(6.0D0-3.0D0/const_vol))
  q_term_2 = -3.0D0*temp_tgt*(2.0D0/beta_2*dTdy + 2.0D0/beta_1*(v_tgt*(dudx+dwdz-dvdy/const_vol) + &
    u_tgt*(dudy+dvdx) + w_tgt*(dvdz+dwdy)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(4)
! 4 - (0,0,+1)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(dudx+dvdy+dwdz))) + &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(dudx+dvdy+dwdz)) - &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (1.0D0-temp_tgt)*(6.0D0*dwdz-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dudx)-2.0D0/const_vol*pri_der)
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(3.0D0/beta_2*dTdz + w_tgt/beta_1*dwdz*(6.0D0-3.0D0/const_vol))
  q_term_2 = -3.0D0*temp_tgt*(2.0D0/beta_2*dTdz + 2.0D0/beta_1*(w_tgt*(dudx+dvdy-dwdz/const_vol) + &
    v_tgt*(dwdy+dvdz) + u_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(5)
! 5 - (-1,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (1.0D0-temp_tgt)*(6.0D0*dudx-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dwdz)-2.0D0/const_vol*pri_der)
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(3.0D0/beta_2*dTdx+u_tgt/beta_1*dudx*(6.0D0-3.0D0/const_vol))
  q_term_2 = -3.0D0*temp_tgt*(2.0D0/beta_2*dTdx + 2.0D0/beta_1*(u_tgt*(dwdz+dvdy-dudx/const_vol) + &
    v_tgt*(dudy+dvdx) + w_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(6)
! 6 - (0,-1,0)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) - &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  p_term_1 = (1.0D0-temp_tgt)*(6.0D0*dvdy - 1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dudx+dwdz) - 2.0D0/const_vol*pri_der)
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(3.0D0/beta_2*dTdy + v_tgt/beta_1*dvdy*(6.0D0-3.0D0/const_vol))
  q_term_2 = -3.0D0*temp_tgt*(2.0D0/beta_2*dTdy + 2.0D0/beta_1*(v_tgt*(dudx+dwdz-dvdy/const_vol) + &
    u_tgt*(dudy+dvdx) + w_tgt*(dvdz+dwdy)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(7)
! 7 - (0,0,-1)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
!
! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) - &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(pri_der)) + &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (1.0D0-temp_tgt)*(6.0D0*dwdz-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dudx)-2.0D0/const_vol*pri_der)
  q_term_1 = (1.0D0-3.0D0*temp_tgt)*(3.0D0/beta_2*dTdz + w_tgt/beta_1*dwdz*(6.0D0-3.0D0/const_vol))
  q_term_2 = -3.0D0*temp_tgt*(2.0D0/beta_2*dTdz + 2.0D0/beta_1*(w_tgt*(dudx+dvdy-dwdz/const_vol) + &
    v_tgt*(dwdy+dvdz) + u_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(8)
! +1,+1,+1
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)

  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx+dTdy+dTdz
  cross_der = dudy+dvdx+dwdx+dudz+dvdz+dwdy
! f_missing
!  fi_neq = k_fi((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(dudy+dudz+dvdx+dvdz+dwdx+dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + &
     v_tgt*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + &
     w_tgt*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
     v_tgt*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
     w_tgt*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(temp_der_sum)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy+dTdz)+ &!term 2, 1/2T^2
            v_tgt*(dTdx+dTdz)+ w_tgt*(dTdx+dTdy)) - &
            2*pri_der/const_vol*(u_tgt*v_tgt+u_tgt*w_tgt+v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(9)
! 9 - (-1,+1,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  cross_der = -dudy-dudz-dvdx+dvdz-dwdx+dwdy
  temp_der_sum = -dTdx+dTdy+dTdz
  pri_der = dudx+dvdy+dwdz

! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(-dudy-dudz-dvdx+dvdz-dwdx+dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdz+dTdy-dTdx) + &
!    u_tgt/beta_1*(-6*dudx+2*(dvdx+dudy+dwdx+dudz-cv_const*(dvdy+dwdz))+5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(-dvdx-dudy+dwdy+dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(-dudz-dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(-dvdy-dwdz+dudy+dvdx+dwdx+dudz+cv_const*dudx)-6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz-dvdx-dudy+dwdy+dvdz-cv_const*dvdy)-6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy-dwdx+dwdy-dudz+dvdz-cv_const*dwdz)-6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1-temp_tgt)*(6*pri_der-3/cv_const*pri_der)
  p_term_2 = 6*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(-6*dudx+2*( dvdx+dudy+dwdx+dudz-dvdy-dwdz)+5*cv_const*dudx) + &
     v_tgt*(6*dvdy + 2*( -dvdx-dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + &
     w_tgt*(6*dwdz + 2*( -dudz-dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(-dvdy-dwdz+dudy+dvdx+dwdx+dudz + cv_const*dudx) - 6*(dwdy+dvdz)) + &
     v_tgt*(4*(-dudx+dwdz-dvdx-dudy+dwdy+dvdz - cv_const*dvdy) - 6*(dudz+dwdx)) + &
     w_tgt*(4*(-dudx+dvdy-dwdx+dwdy-dudz+dvdz - cv_const*dwdz) - 6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq

!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(temp_der_sum)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(-u_tgt+v_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)-u_tgt*(dTdy+dTdz)+ & !term 2, 1/2T^2
            v_tgt*(-dTdx+dTdz)+ w_tgt*(-dTdx+dTdy)) - &
            2*pri_der/const_vol*(-u_tgt*v_tgt-u_tgt*w_tgt+v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

 !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
 loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(10)
! 10 - (-1,-1,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  cross_der = dudy-dudz+dvdx-dvdz-dwdx-dwdy
  temp_der_sum = -dTdx-dTdy+dTdz
  pri_der = dudx+dvdy+dwdz

! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(dudy-dudz+dvdx-dvdz-dwdx-dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(-dTdx-dTdy+dTdz) + &
!    u_tgt/beta_1*(-6*dudx+2*(dvdx+dudy+dwdx+dudz-cv_const*(dvdy+dwdz))+5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(-6*dvdy+2*(dvdx+dudy+dwdy+dvdz-cv_const*(dudx+dwdz))+5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz+cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz+cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(-6*dudx+2*(-dvdx-dudy+dwdx+dudz-dvdy-dwdz)+5*cv_const*dudx) + &
     v_tgt*(-6*dvdy+2*(-dvdx-dudy+dwdy+dvdz-dudx-dwdz)+5*cv_const*dvdy) + &
     w_tgt*(6*dwdz+2*(-dudz-dwdx-dvdz-dwdy+dudx+dvdy)-5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(-dvdy-dwdz-dudy-dvdx+dwdx+dudz+cv_const*dudx)+6*(dwdy+dvdz)) + &
     v_tgt*(4*(-dudx-dwdz-dvdx-dudy+dwdy+dvdz+cv_const*dvdy)+6*(dudz+dwdx)) + &
     w_tgt*(4*(dudx+dvdy-dwdx-dwdy-dudz-dvdz-cv_const*dwdz)+6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(-dTdx-dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(-u_tgt-v_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(-dTdy+dTdz) + &!term 2, 1/2T^2
            v_tgt*(-dTdx+dTdz) - w_tgt*(dTdx+dTdy)) - &
            2*pri_der/const_vol*(u_tgt*v_tgt-u_tgt*w_tgt-v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(11)
! 11 - (+1,-1,+1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  cross_der = -dudy+dudz-dvdx-dvdz+dwdx-dwdy
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx-dTdy+dTdz
! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(-dudy+dudz-dvdx-dvdz+dwdx-dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx-dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(-dvdx-dudy+dwdx+dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(-6*dvdy+2*(dvdx+dudy+dwdy+dvdz-cv_const*(dudx+dwdz))+5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx-dvdz-dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx-dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz-dudy-dvdx+dwdx+dudz-cv_const*dudx)-6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(-dudx-dwdz+dvdx+dudy+dwdy+dvdz+cv_const*dvdy)-6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy-dwdx+dwdy+dudz-dvdz-cv_const*dwdz)-6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(6*dudx+2*(-dvdx-dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + &
     v_tgt*(-6*dvdy+2*(dvdx+dudy+dwdy+dvdz-dudx-dwdz)+5*cv_const*dvdy) + &
     w_tgt*(6*dwdz+2*(dudz+dwdx-dvdz-dwdy+dudx+dvdy)-5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(-dvdy+dwdz-dudy-dvdx+dwdx+dudz-cv_const*dudx)-6*(dwdy+dvdz)) + &
     v_tgt*(4*(-dudx-dwdz+dvdx+dudy+dwdy+dvdz+cv_const*dvdy)-6*(dudz+dwdx)) + &
     w_tgt*(4*(dudx-dvdy+dwdx-dwdy+dudz-dvdz-cv_const*dwdz)-6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx-dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(-dTdy+dTdz) - &!term 2, 1/2T^2
            v_tgt*(dTdx+dTdz)+ w_tgt*(dTdx-dTdy)) - &
            2*pri_der/const_vol*(-u_tgt*v_tgt+u_tgt*w_tgt-v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(12)
! 12 - (+1,+1,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  cross_der = dudy-dudz+dvdx-dvdz-dwdx-dwdy
  temp_der_sum = dTdx+dTdy-dTdz
  pri_der = dudx+dvdy+dwdz

! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(dudy-dudz+dvdx-dvdz-dwdx-dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy-dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy-dwdx-dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy-dwdy-dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(-6*dwdz+2*(dudz+dwdx+dvdz+dwdy-cv_const*(dudx+dvdy))+5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy-dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx-dwdx-dudz-cv_const*dudx)-6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy-dwdy-dvdz-cv_const*dvdy)-6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(-dudx-dvdy+dwdx+dwdy+dudz+dvdz+cv_const*dwdz)-6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(6*dudx+2*(dvdx+dudy-dwdx-dudz+dvdy+dwdz)-5*cv_const*dudx) + &
     v_tgt*(6*dvdy+2*(dvdx+dudy-dwdy-dvdz+dudx+dwdz)-5*cv_const*dvdy) + &
     w_tgt*(-6*dwdz+2*(dudz+dwdx+dvdz+dwdy-dudx-dvdy)+5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(dvdy-dwdz+dudy+dvdx-dwdx-dudz-cv_const*dudx)-6*(dwdy+dvdz)) + &
     v_tgt*(4*(dudx-dwdz+dvdx+dudy-dwdy-dvdz-cv_const*dvdy)-6*(dudz+dwdx)) + &
     w_tgt*(4*(-dudx-dvdy+dwdx+dwdy+dudz+dvdz+cv_const*dwdz)-6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt-w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy-dTdz)+ &!term 2, 1/2T^2
            v_tgt*(dTdx-dTdz) - w_tgt*(dTdx+dTdy)) - &
            2*pri_der/const_vol*(u_tgt*v_tgt - u_tgt*w_tgt - v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(13)
! 13 - (-1,+1,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)

  cross_der = -dudy+dudz-dvdx-dvdz+dwdx-dwdy
  temp_der_sum = -dTdx+dTdy-dTdz
  pri_der = dudx+dvdy+dwdz

! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(cross_der) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(-6*dudx+2*(dvdx+dudy-dwdx-dudz-dvdy-dwdz)+5*cv_const*dudx) + &
     v_tgt*(6*dvdy+2*(-dvdx-dudy-dwdy-dvdz+dudx+dwdz)-5*cv_const*dvdy) + &
     w_tgt*(-6*dwdz+2*(-dudz-dwdx+dvdz+dwdy-dudx-dvdy)+5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(-dvdy-dwdz+dudy+dvdx-dwdx-dudz+cv_const*dudx)+6*(dwdy+dvdz)) + &
     v_tgt*(4*(dudx+dwdz-dvdx-dudy-dwdy-dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
     w_tgt*(4*(-dudx-dvdy-dwdx+dwdy-dudz+dvdz+cv_const*dwdz)+6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(-dTdx+dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(-u_tgt+v_tgt-w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy-dTdz) - &!term 2, 1/2T^2
            v_tgt*(dTdx+dTdz)+ w_tgt*(-dTdx+dTdy)) - &
            2*pri_der/const_vol*(-u_tgt*v_tgt+u_tgt*w_tgt-v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(14)
! 14 - (-1,-1,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  cross_der = dudy+dudz+dvdx+dvdz+dwdx+dwdy
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx+dTdy+dTdz
! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der)) + & !first term
!    6/(2*beta_1*temp_tgt**2)*(cross_der) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) - &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz) - &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + &
     v_tgt*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + &
     w_tgt*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
     v_tgt*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
     w_tgt*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt+w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy+dTdz)+ &!term 2, 1/2T^2
            v_tgt*(dTdx+dTdz)+ w_tgt*(dTdx+dTdy)) - &
            2*pri_der/const_vol*(u_tgt*v_tgt+u_tgt*w_tgt+v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(15)
! 15 - (+1,-1,-1)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  cross_der = -dudy-dudz-dvdx+dvdz-dwdx+dwdy
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx-dTdy-dTdz
! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(cross_der) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
  p_term_1 = (1.0D0-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
  p_term_2 = 6.0D0*cross_der
  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(temp_der_sum) + &
    (u_tgt*(6*dudx+2*(-dvdx-dudy-dwdx-dudz+dvdy+dwdz)-5*cv_const*dudx) + &
     v_tgt*(-6*dvdy+2*(dvdx+dudy-dwdy-dvdz-dudx-dwdz)+5*cv_const*dvdy) + &
     w_tgt*(-6*dwdz+2*(dudz+dwdx-dvdz-dwdy-dudx-dvdy)+5*cv_const*dwdz))/beta_1)
  q_term_2 = 4/beta_2*(temp_der_sum) + &
    (u_tgt*(4*(dvdy+dwdz-dudy-dvdx-dwdx-dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
     v_tgt*(4*(-dudx-dwdz+dvdx+dudy-dwdy-dvdz+cv_const*dvdy)+6*(dudz+dwdx)) + &
     w_tgt*(4*(-dudx-dvdy+dwdx+dwdy-dudz-dvdz+cv_const*dwdz)+6*(dudy+dvdx)))/beta_1

  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx-dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt-w_tgt)
  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der) - u_tgt*(dTdy+dTdz)+ &!term 2, 1/2T^2
            v_tgt*(-dTdx+dTdz)+ w_tgt*(-dTdx+dTdy)) - &
            2*pri_der/const_vol*(-u_tgt*v_tgt-u_tgt*w_tgt+v_tgt*w_tgt)) !
  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(16)
! 16 - (+2,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) + & !first term
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - & !third term
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) - & !second term
!    1/(temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & !fourth term
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (4.0D0-temp_tgt)*(6.0D0*dudx-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dwdz)-2.0D0/const_vol*pri_der)
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(3.0D0/beta_2*dTdx+u_tgt/beta_1*dudx*(6.0D0-3.0D0/const_vol))
  q_term_2 = -6.0D0*temp_tgt*(2.0D0/beta_2*dTdx + 2.0D0/beta_1*(u_tgt*(dwdz+dvdy-dudx/const_vol) + &
    v_tgt*(dudy+dvdx) + w_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(17)
! 17 - (0,+2,0)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) + & !first term
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - & !third term
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) - & !second term
!    1/(temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & !fourth term
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (4.0D0-temp_tgt)*(6.0D0*dvdy - 1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dudx+dwdz) - 2.0D0/const_vol*pri_der)
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(3.0D0/beta_2*dTdy + v_tgt/beta_1*dvdy*(6.0D0-3.0D0/const_vol))
  q_term_2 = -6.0D0*temp_tgt*(2.0D0/beta_2*dTdy + 2.0D0/beta_1*(v_tgt*(dudx+dwdz-dvdy/const_vol) + &
    u_tgt*(dudy+dvdx) + w_tgt*(dvdz+dwdy)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt + u_tgt**2 + w_tgt**2))
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(18)
! 18 - (0,0,+2)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(dudx+dvdy+dwdz))) + &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(dudx+dvdy+dwdz)) - &
!    1/(temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (4.0D0-temp_tgt)*(6.0D0*dwdz-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dudx)-2.0D0/const_vol*pri_der)
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(3.0D0/beta_2*dTdz + w_tgt/beta_1*dwdz*(6.0D0-3.0D0/const_vol))
  q_term_2 = -6.0D0*temp_tgt*(2.0D0/beta_2*dTdz + 2.0D0/beta_1*(w_tgt*(dudx+dvdy-dwdz/const_vol) + &
    v_tgt*(dwdy+dvdz) + u_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(19)
! 19 - (-2,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    1/(temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (4.0D0-temp_tgt)*(6.0D0*dudx-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dwdz)-2.0D0/const_vol*pri_der)
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(3.0D0/beta_2*dTdx+u_tgt/beta_1*dudx*(6.0D0-3.0D0/const_vol))
  q_term_2 = -6.0D0*temp_tgt*(2.0D0/beta_2*dTdx + 2.0D0/beta_1*(u_tgt*(dwdz+dvdy-dudx/const_vol) + &
    v_tgt*(dudy+dvdx) + w_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(20)
! 20 - (0,-2,0)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) - &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    1/(temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  p_term_1 = (4.0D0-temp_tgt)*(6.0D0*dvdy - 1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dudx+dwdz) - 2.0D0/const_vol*pri_der)
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(3.0D0/beta_2*dTdy + v_tgt/beta_1*dvdy*(6.0D0-3.0D0/const_vol))
  q_term_2 = -6.0D0*temp_tgt*(2.0D0/beta_2*dTdy + 2.0D0/beta_1*(v_tgt*(dudx+dwdz-dvdy/const_vol) + &
    u_tgt*(dudy+dvdx) + w_tgt*(dvdz+dwdy)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(21)
! 21 - (0,0,-2)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) - &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(pri_der)) + &
!    1/(temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (4.0D0-temp_tgt)*(6.0D0*dwdz-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dudx)-2.0D0/const_vol*pri_der)
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(3.0D0/beta_2*dTdz + w_tgt/beta_1*dwdz*(6.0D0-3.0D0/const_vol))
  q_term_2 = -6.0D0*temp_tgt*(2.0D0/beta_2*dTdz + 2.0D0/beta_1*(w_tgt*(dudx+dvdy-dwdz/const_vol) + &
    v_tgt*(dwdy+dvdz) + u_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(22)
! 22 - (+2,+2,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx+dTdy
  cross_der = dudy+dvdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = 24*cross_der
  p_term_3 = -temp_tgt*(6*dwdz-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (u_tgt*(6*dudx + 2*(cross_der+dvdy) - 4*cv_const*dudx)+ &
      v_tgt*(6*dvdy + 2*(cross_der+dudx) - 4*cv_const*dvdy))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (u_tgt*(4*(cross_der+dvdy) - cv_const_2*dudx)+&
      v_tgt*(4*(cross_der+dudx) - cv_const_2*dvdy))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*w_tgt*(dudz+dwdx+dvdz+dwdy) + &
      u_tgt*(2*dwdz-cv_const*dudx) + &
      v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) + &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(23)
! 23 - (+2,0,+2)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx+dTdz
  cross_der = dudz+dwdx
! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudz+dwdx)
!  p_term_3 = 6*dvdy-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx+dTdz)+(u_tgt*(6*dudx+2*(dwdx+dudz+dwdz)-4*cv_const*dudx)+&
!    w_tgt*(6*dwdz+2*(dwdx+dudz+dudx)-4*cv_const*dwdz))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx+dTdz)+(u_tgt*(4*(dudz+dwdx+dwdz)-cv_const_2*dudx)+&
!    w_tgt*(4*(dudz+dwdx+dudx)-cv_const_2*dwdz))/ubeta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx+dvdz+dwdy) + &
!    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = 24*cross_der
  p_term_3 = -6*temp_tgt*(dvdy-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
   (u_tgt*(6*dudx + 2*(cross_der+dwdz) - 4*cv_const*dudx)+ &
    w_tgt*(6*dwdz + 2*(cross_der+dudx) - 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
   (u_tgt*(4*(cross_der+dwdz) - cv_const_2*dudx)+&
    w_tgt*(4*(cross_der+dudx) - cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
   (2.0D0*v_tgt*(dvdx+dudy+dwdy+dvdz) + &
    u_tgt*(2*dvdy-cv_const*dudx) + &
    w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+w_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) + &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(24)
! 24 - (0,+2,+2)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdy+dTdz
  cross_der = dvdz+dwdy
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*pri_der)
  p_term_2 = 24*cross_der
  p_term_3 = -6*temp_tgt*(dudx-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (v_tgt*(6*dvdy + 2*(cross_der+dwdz) - 4*cv_const*dvdy)+ &
      w_tgt*(6*dwdz + 2*(cross_der+dvdy) - 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (v_tgt*(4*(cross_der+dwdz) - cv_const_2*dvdy)+&
      w_tgt*(4*(cross_der+dvdy) - cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*u_tgt*(dudy+dvdx+dwdx+dudz) + &
      v_tgt*(2*dudx-cv_const*dvdy) + &
      w_tgt*(2*dudx-cv_const*dwdz))/beta_1)
!
  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt+v_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) + &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(25)
! 25 - (-2,+2,0)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = -dTdx+dTdy
  cross_der = dudy+dvdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = -24*cross_der
  p_term_3 = -temp_tgt*(6*dwdz-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (u_tgt*(-6*dudx + 2*(cross_der-dvdy) + 4*cv_const*dudx)+ &
      v_tgt*(6*dvdy - 2*(cross_der-dudx) - 4*cv_const*dvdy))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (u_tgt*(4*(cross_der-dvdy) + cv_const_2*dudx)+&
      v_tgt*(-4*(cross_der-dudx) - cv_const_2*dvdy))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*w_tgt*(-dudz-dwdx+dvdz+dwdy) - &
      u_tgt*(2*dwdz-cv_const*dudx) + &
      v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)
  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdy-dTdx)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(v_tgt-u_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) - &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq
  case(26)
! 26 - (-2,-2,0)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx+dTdy
  cross_der = dudy+dvdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = 24*cross_der
  p_term_3 = -temp_tgt*(6*dwdz-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (u_tgt*(6*dudx + 2*(cross_der+dvdy) - 4*cv_const*dudx)+ &
      v_tgt*(6*dvdy + 2*(cross_der+dudx) - 4*cv_const*dvdy))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (u_tgt*(4*(cross_der+dvdy) - cv_const_2*dudx)+&
      v_tgt*(4*(cross_der+dudx) - cv_const_2*dvdy))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*w_tgt*(dudz+dwdx+dvdz+dwdy) + &
      u_tgt*(2*dwdz-cv_const*dudx) + &
      v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)
  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) + &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(27)
! 27 - (+2,-2,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx-dTdy
  cross_der = dudy+dvdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
  p_term_2 = -24*cross_der
  p_term_3 = -temp_tgt*(6*dwdz-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (u_tgt*(6*dudx - 2*(cross_der-dvdy) - 4*cv_const*dudx)+ &
      v_tgt*(-6*dvdy + 2*(cross_der-dudx) + 4*cv_const*dvdy))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (u_tgt*(-4*(cross_der-dvdy) - cv_const_2*dudx)+&
      v_tgt*(4*(cross_der-dudx) + cv_const_2*dvdy))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*w_tgt*(dudz+dwdx-dvdz-dwdy) + &
      u_tgt*(2*dwdz-cv_const*dudx) - &
      v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx-dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) - &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(28)
! 28 - (+2,0,-2)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx-dTdz
  cross_der = dudz+dwdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = -24*cross_der
  p_term_3 = -6*temp_tgt*(dvdy-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
   (u_tgt*(6*dudx - 2*(cross_der-dwdz) - 4*cv_const*dudx) + &
    w_tgt*(-6*dwdz + 2*(cross_der-dudx) + 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum) + &
   (u_tgt*(-4*(cross_der-dwdz) - cv_const_2*dudx) + &
    w_tgt*(4*(cross_der-dudx) + cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
   (2.0D0*v_tgt*(dvdx+dudy-dwdy-dvdz) + &
    u_tgt*(2*dvdy-cv_const*dudx) - &
    w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-w_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) - &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(29)
! 29 - (-2,0,-2)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdx+dTdz
  cross_der = dudz+dwdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = 24*cross_der
  p_term_3 = -6*temp_tgt*(dvdy-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
   (u_tgt*(6*dudx + 2*(cross_der+dwdz) - 4*cv_const*dudx)+ &
    w_tgt*(6*dwdz + 2*(cross_der+dudx) - 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
   (u_tgt*(4*(cross_der+dwdz) - cv_const_2*dudx)+&
    w_tgt*(4*(cross_der+dudx) - cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
   (2.0D0*v_tgt*(dvdx+dudy+dwdy+dvdz) + &
    u_tgt*(2*dvdy-cv_const*dudx) + &
    w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdx+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+w_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) + &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(30)
! 30 - (-2,0,+2)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = -dTdx+dTdz
  cross_der = dudz+dwdx
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
  p_term_2 = -24*cross_der
  p_term_3 = -6*temp_tgt*(dvdy-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
   (u_tgt*(-6*dudx + 2*(cross_der-dwdz) + 4*cv_const*dudx)+ &
    w_tgt*(6*dwdz - 2*(cross_der-dudx) - 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
   (u_tgt*(4*(cross_der-dwdz) + cv_const_2*dudx)+&
    w_tgt*(-4*(cross_der-dudx) - cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
   (2.0D0*v_tgt*(-dvdx-dudy+dwdy+dvdz) - &
    u_tgt*(2*dvdy-cv_const*dudx) + &
    w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 + p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdz-dTdx)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt-u_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) - &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(31)
! 31 - (0,+2,-2)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdy-dTdz
  cross_der = dvdz+dwdy
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*pri_der)
  p_term_2 = -24*cross_der
  p_term_3 = -6*temp_tgt*(dudx-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (v_tgt*(6*dvdy + 2*(cross_der+dwdz) - 4*cv_const*dvdy)+ &
      w_tgt*(-6*dwdz - 2*(cross_der-dvdy) + 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (v_tgt*(-4*(cross_der-dwdz) - cv_const_2*dvdy)+&
      w_tgt*(4*(cross_der-dvdy) + cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*u_tgt*(dudy+dvdx-dwdx-dudz) + &
      v_tgt*(2*dudx-cv_const*dvdy) - &
      w_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(const_vol*beta_1)
  g_alpha = 9*(dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(v_tgt-w_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) - &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(32)
! 32 - (0,-2,-2)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = dTdy+dTdz
  cross_der = dvdz+dwdy
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*pri_der)
  p_term_2 = 24*cross_der
  p_term_3 = -6*temp_tgt*(dudx-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (v_tgt*(6*dvdy + 2*(cross_der+dwdz) - 4*cv_const*dvdy)+ &
      w_tgt*(6*dwdz + 2*(cross_der+dvdy) - 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (v_tgt*(4*(cross_der+dwdz) - cv_const_2*dvdy)+&
      w_tgt*(4*(cross_der+dvdy) - cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*u_tgt*(dudy+dvdx+dwdx+dudz) + &
      v_tgt*(2*dudx-cv_const*dvdy) + &
      w_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))

  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt+v_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) + &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(33)
! 33 - (0,-2,+2)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
  temp_der_sum = -dTdy+dTdz
  cross_der = dvdz+dwdy
! f
  p_term_1 = (4.0D0-temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*pri_der)
  p_term_2 = -24*cross_der
  p_term_3 = -6*temp_tgt*(dudx-cv_const*pri_der)
!
  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*temp_der_sum/beta_2 + &
     (v_tgt*(-6*dvdy - 2*(cross_der-dwdz) + 4*cv_const*dvdy)+ &
      w_tgt*(6*dwdz + 2*(cross_der-dvdy) - 4*cv_const*dwdz))/beta_1)
  q_term_2 = 8.0D0*(2.0D0/beta_2*(temp_der_sum)+&
     (v_tgt*(4*(cross_der-dwdz) + cv_const_2*dvdy)+&
      w_tgt*(-4*(cross_der-dvdy) - cv_const_2*dwdz))/beta_1)
  q_term_3 = (-6.0D0*temp_tgt)*(temp_der_sum/beta_2 + &
     (2.0D0*u_tgt*(-dudy-dvdx+dwdx+dudz) - &
      v_tgt*(2*dudx-cv_const*dvdy) + &
      w_tgt*(2*dudx-cv_const*dwdz))/beta_1)

  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*(dTdz-dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt-v_tgt)
  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
            1/(2*temp_tgt)*(6*(temp_tgt*dudx+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) - &
            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(34)
! 34 - (+3,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-(pri_der)/const_vol)) + & ! term 1
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*dudx/const_vol)) - & ! term 2
!    1/(2*beta_1*temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*(pri_der)) - & ! term 3
!    9/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & ! term 4
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (9.0D0-temp_tgt)*(6.0D0*dudx-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dwdz)-2.0D0/const_vol*pri_der)
  q_term_1 = (27.0D0-9.0D0*temp_tgt)*(3.0D0/beta_2*dTdx+u_tgt/beta_1*dudx*(6.0D0-3.0D0/const_vol))
  q_term_2 = -9.0D0*temp_tgt*(2.0D0/beta_2*dTdx + 2.0D0/beta_1*(u_tgt*(dwdz+dvdy-dudx/const_vol) + &
    v_tgt*(dudy+dvdx) + w_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
  gi_neq = k_gi*(k_one + 3*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(35)
! 35 - (0,+3,0)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(pri_der))) + &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy + v_tgt/beta_1*(6*dvdy - 3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*(dudx+dwdz)-cv_const_2*(pri_der)) - &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  p_term_1 = (9.0D0-temp_tgt)*(6.0D0*dvdy - 1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dudx+dwdz) - 2.0D0/const_vol*pri_der)
  q_term_1 = (27.0D0-9.0D0*temp_tgt)*(3.0D0/beta_2*dTdy + v_tgt/beta_1*dvdy*(6.0D0-3.0D0/const_vol))
  q_term_2 = -9.0D0*temp_tgt*(2.0D0/beta_2*dTdy + 2.0D0/beta_1*(v_tgt*(dudx+dwdz-dvdy/const_vol) + &
    u_tgt*(dudy+dvdx) + w_tgt*(dvdz+dwdy)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq

! g
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one + 3*g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(36)
! 36 - (0,0,+3)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,downwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,upwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) + &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz + w_tgt/beta_1*(6*dwdz - 3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*(dudx+dvdy)-cv_const_2*(pri_der)) - &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dudx+dvdy-cv_const*dwdz)+2*u_tgt/beta_1* &
!    (dudz+dwdx)+2*v_tgt/beta_1*(dvdz+dwdy)))
  p_term_1 = (9.0D0-temp_tgt)*(6.0D0*dwdz-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dudx)-2.0D0/const_vol*pri_der)
  q_term_1 = (27.0D0-9.0D0*temp_tgt)*(3.0D0/beta_2*dTdz + w_tgt/beta_1*dwdz*(6.0D0-3.0D0/const_vol))
  q_term_2 = -9.0D0*temp_tgt*(2.0D0/beta_2*dTdz + 2.0D0/beta_1*(w_tgt*(dudx+dvdy-dwdz/const_vol) + &
    v_tgt*(dwdy+dvdz) + u_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) + &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
  gi_neq = k_gi*(k_one + 3*g_alpha + r_ab/beta_1)


  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(37)
! 37 - (-3,0,0)
  dudx = dudx_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdx = dTdx_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (9.0D0-temp_tgt)*(6.0D0*dudx-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dwdz)-2.0D0/const_vol*pri_der)
  q_term_1 = (27.0D0-9.0D0*temp_tgt)*(3.0D0/beta_2*dTdx+u_tgt/beta_1*dudx*(6.0D0-3.0D0/const_vol))
  q_term_2 = -9.0D0*temp_tgt*(2.0D0/beta_2*dTdx + 2.0D0/beta_1*(u_tgt*(dwdz+dvdy-dudx/const_vol) + &
    v_tgt*(dudy+dvdx) + w_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - 3*g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  case(38)
! 38 - (0,-3,0)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dudz = dudz_local(loc_u(1),loc_u(7),loc_u(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dvdz = dvdz_local(loc_v(1),loc_v(7),loc_v(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dwdz = dwdz_local(loc_w(1),loc_w(7),loc_w(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdy = dTdy_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
  dTdz = dTdz_local(loc_temp(1),loc_temp(7),loc_temp(4),rest_state(a,b,c,1),rest_state(a,b,c,7),rest_state(a,b,c,4))
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(pri_der))) - & !P term 1
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - & !Q term 1
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(pri_der)) + & !P term 2
!    9/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* & !Q term 2
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
  p_term_1 = (9.0D0-temp_tgt)*(6.0D0*dvdy - 1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dudx+dwdz) - 2.0D0/const_vol*pri_der)
  q_term_1 = (27.0D0-9.0D0*temp_tgt)*(3.0D0/beta_2*dTdy + v_tgt/beta_1*dvdy*(6.0D0-3.0D0/const_vol))
  q_term_2 = -9.0D0*temp_tgt*(2.0D0/beta_2*dTdy + 2.0D0/beta_1*(v_tgt*(dudx+dwdz-dvdy/const_vol) + &
    u_tgt*(dudy+dvdx) + w_tgt*(dvdz+dwdy)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq
!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
  gi_neq = k_gi*(k_one - 3*g_alpha + r_ab/beta_1)

  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  loc_gout = loc_gieq + loc_epsilon*gi_neq

  case(39)
! 39 - (0,0,-3)
  dudx = dudx_local(loc_u(1),loc_u(5),loc_u(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dudy = dudy_local(loc_u(1),loc_u(6),loc_u(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dudz = dudz_ibb(u_tgt,0.0D0,downwind,dir_state,i)
  dvdx = dvdx_local(loc_v(1),loc_v(5),loc_v(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dvdy = dvdy_local(loc_v(1),loc_v(6),loc_v(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dvdz = dvdz_ibb(v_tgt,0.0D0,downwind,dir_state,i)
  dwdx = dwdx_local(loc_w(1),loc_w(5),loc_w(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dwdy = dwdy_local(loc_w(1),loc_w(6),loc_w(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
  dwdz = dwdz_ibb(w_tgt,0.0D0,downwind,dir_state,i)
  dTdx = dTdx_local(loc_temp(1),loc_temp(5),loc_temp(2),rest_state(a,b,c,1),rest_state(a,b,c,5),rest_state(a,b,c,2))
  dTdy = dTdy_local(loc_temp(1),loc_temp(6),loc_temp(3),rest_state(a,b,c,1),rest_state(a,b,c,6),rest_state(a,b,c,3))
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
  dTdz = dTdz_ibb(temp_tgt,temp_tgt,downwind,dir_state,i)
!
  pri_der = dudx+dvdy+dwdz
! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) - &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(pri_der)) + &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
  p_term_1 = (9.0D0-temp_tgt)*(6.0D0*dwdz-1.0D0/const_vol*pri_der)
  p_term_2 = -temp_tgt*(6.0D0*(dvdy+dudx)-2.0D0/const_vol*pri_der)
  q_term_1 = (27.0D0-9.0D0*temp_tgt)*(3.0D0/beta_2*dTdz + w_tgt/beta_1*dwdz*(6.0D0-3.0D0/const_vol))
  q_term_2 = -9.0D0*temp_tgt*(2.0D0/beta_2*dTdz + 2.0D0/beta_1*(w_tgt*(dudx+dvdy-dwdz/const_vol) + &
    v_tgt*(dwdy+dvdz) + u_tgt*(dudz+dwdx)))

  fi_neq = k_fi*((p_term_1+p_term_2)/(2.0D0*beta_1*temp_tgt**2) - &
    (q_term_1+q_term_2)/(6.0D0*temp_tgt**3))
  !loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
  loc_fout = loc_fieq + loc_epsilon*fi_neq

!
! g missing
!
  k_one = 9*pri_der/(beta_1*const_vol)
  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
  gi_neq = k_gi*(k_one - 3*g_alpha + r_ab/beta_1)

  loc_gout = loc_gieq + loc_epsilon*gi_neq
  !loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
  end select

  else
    dudx = first_order_vel_cond_unknown(loc_u(1),u_tgt,loc_u(5),loc_u(2),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,5),rest_state(a,b,c,2))
    dudy = first_order_vel_cond_unknown(loc_u(1),u_tgt,loc_u(6),loc_u(3),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,6),rest_state(a,b,c,3))
    dudz = first_order_vel_cond_unknown(loc_u(1),u_tgt,loc_u(7),loc_u(4),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,7),rest_state(a,b,c,4))
    dvdx = first_order_vel_cond_unknown(loc_v(1),v_tgt,loc_v(5),loc_v(2),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,5),rest_state(a,b,c,2))
    dvdy = first_order_vel_cond_unknown(loc_v(1),v_tgt,loc_v(6),loc_v(3),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,6),rest_state(a,b,c,3))
    dvdz = first_order_vel_cond_unknown(loc_v(1),v_tgt,loc_v(7),loc_v(4),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,7),rest_state(a,b,c,4))
    dwdx = first_order_vel_cond_unknown(loc_w(1),w_tgt,loc_w(5),loc_w(2),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,5),rest_state(a,b,c,2))
    dwdy = first_order_vel_cond_unknown(loc_w(1),w_tgt,loc_w(6),loc_w(3),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,6),rest_state(a,b,c,3))
    dwdz = first_order_vel_cond_unknown(loc_w(1),w_tgt,loc_w(7),loc_w(4),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,7),rest_state(a,b,c,4))
    dTdx = first_order_temp_cond_unknown(loc_temp(1),temp_tgt,loc_temp(5),loc_temp(2),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,5),rest_state(a,b,c,2))
    dTdy = first_order_temp_cond_unknown(loc_temp(1),temp_tgt,loc_temp(6),loc_temp(3),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,6),rest_state(a,b,c,3))
    dTdz = first_order_temp_cond_unknown(loc_temp(1),temp_tgt,loc_temp(7),loc_temp(4),i,rest_state(a,b,c,1),&
             rest_state(a,b,c,7),rest_state(a,b,c,4))
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

    pab = k_fi/beta_1*pab

    qabc(1) = 3*k_fi*(dTdx/beta_2 + u_tgt*pab(1))    !xxx
    qabc(2) = k_fi*(dTdy/beta_2 + v_tgt*pab(2) + 2*u_tgt*pab(4)) !xxy
    qabc(3) = k_fi*(dTdz/beta_2 + w_tgt*pab(3) + 2*u_tgt*pab(7)) !xxz
    qabc(4) = k_fi*(dTdy/beta_2 + 2*u_tgt*pab(4) + v_tgt*pab(1)) !xyx
    qabc(5) = k_fi*(dTdx/beta_2 + u_tgt*pab(13) + 2*v_tgt*pab(4))!xyy
    qabc(6) = k_fi*(u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!xyz
    qabc(7) = k_fi*(dTdz/beta_2 + w_tgt*pab(3) + 2*u_tgt*pab(7))!xzx
    qabc(8) = k_fi*(u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!xzy
    qabc(9) = k_fi*(dTdx/beta_2 + u_tgt*pab(25) + 2*w_tgt*pab(7))!xzz
    qabc(10) = k_fi*(dTdy/beta_2 + v_tgt*pab(2) + 2*u_tgt*pab(4))!yxx
    qabc(11) = k_fi*(dTdx/beta_2 + 2*v_tgt*pab(4) + u_tgt*pab(13))!yxy
    qabc(12) = k_fi*(u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4))!yxz
    qabc(13) = k_fi*(dTdx/beta_2 + u_tgt*pab(13) + 2*v_tgt*pab(4))!yyx
    qabc(14) = 3*k_fi*(dTdy/beta_2 + v_tgt*pab(14))!yyy
    qabc(15) = k_fi*(dTdz/beta_2 + w_tgt*pab(15) + 2*v_tgt*pab(22))!yyz
    qabc(16) = k_fi*(u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!yzx
    qabc(17) = k_fi*(dTdz/beta_2 + w_tgt*pab(15) + 2*v_tgt*pab(22))!yzy
    qabc(18) = k_fi*(dTdy/beta_2 + v_tgt*pab(26) + 2*w_tgt*pab(22))!yzz
    qabc(19) = k_fi*(dTdx/beta_2 + w_tgt*pab(3) + 2*u_tgt*pab(7))!zxx
    qabc(20) = k_fi*(u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!zxy
    qabc(21) = k_fi*(dTdx/beta_2 + u_tgt*pab(25) + 2*w_tgt*pab(7))!zxz
    qabc(22) = k_fi*(u_tgt*pab(22) + v_tgt*pab(7) + w_tgt*pab(4) )!zyx
    qabc(23) = k_fi*(dTdz/beta_2 + w_tgt*pab(15) + 2*v_tgt*pab(22))!zyy
    qabc(24) = k_fi*(dTdy/beta_2 + v_tgt*pab(26) + 2*w_tgt*pab(22))!zyz
    qabc(25) = k_fi*(dTdx/beta_2 + u_tgt*pab(25) + 2*w_tgt*pab(7))!zzx
    qabc(26) = k_fi*(dTdy/beta_2 + v_tgt*pab(26) + 2*w_tgt*pab(22))!zzy
    qabc(27) = 3*k_fi*(dTdz/beta_2 + w_tgt*pab(27))!zzz

    fi_neq = sum(pab*(p_coeff(:,i)-temp_tgt*p_kroen))/(2*temp_tgt**2) + &
             sum(qabc*(yi(:,i)-3*temp_tgt*kroen_val(:,i)))/(6*temp_tgt**3)
!    loc_fout = loc_fieq + loc_epsilon*fi_neq
    loc_fout = loc_fieq_opp + loc_epsilon*fi_neq
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

    gi_neq = k_gi*(3*k_one + g_alpha + sum(rab*p_coeff(:,i) - &
               temp_tgt*p_kroen)/(2*beta_1*temp_tgt**2))
!    loc_gout = loc_gieq + loc_epsilon*gi_neq
    loc_gout = loc_gieq_opp + loc_epsilon*gi_neq
    if (loc_fout < 0.0D0) then
      loc_fout = -loc_fout
    end if
    if (loc_gout < 0.0D0) then
      loc_gout = -loc_gout
    end if
  end if


end if


if (loc_fout < 0.0D0 .or. loc_gout < 0.0D0 .or. loc_fout > 1.0D0 .or. loc_gout > 1.0D0 .or. &
    isnan(loc_fout) .or. isnan(loc_gout)) then
!    write(*,*) 'check after ibb',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),rest_state(a,b,c,opp(i)),&
!      rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl
    write(*,*) 'outgoing IBB values, fout bad',loc_fout,loc_fi,loc_fi_opp,loc_fieq,loc_fieq_opp,fout_adj,&
      loc_gout,loc_gi,loc_gi_opp,loc_gieq,loc_gieq_opp,gout_adj,rest_state(a,b,c,i),dir_state,loc_rho,u_tgt,&
      v_tgt,w_tgt,temp_tgt,i,lvl,a,b,c
    write(*,*) 'IBB bad states',rest_state(a,b,c,1:dir),a,b,c
    write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,i
    write(*,*) 'miscellaneous',beta_1,beta_2,const_vol,k_fi,k_gi,loc_epsilon
    write(*,*) 'ibb temperature derivatives',dTdx,dTdy,dTdz,i
    write(*,*) 'ibb flow values',loc_rho,loc_u,loc_v,loc_w,loc_temp,a,b,c

    write(*,*) 'ibb fi information',pri_der,p_term_1,p_term_2,p_term_3,q_term_1,q_term_2,q_term_3,&
      fi_neq*loc_epsilon,loc_fout,loc_fieq_opp,fi_neq,inc_fout
    write(*,*) 'ibb gi information',k_one,g_alpha,r_ab,gi_neq*loc_epsilon,loc_epsilon, loc_gout,&
      loc_gieq,gi_neq,inc_gout

    write(*,*) 'pab',pab
    write(*,*) 'qabc',qabc
    write(*,*) 'k_one',k_one
    write(*,*) 'g_alpha',g_alpha
    write(*,*) 'rab',rab
end if


!if (a == 191 .and. b == 1 .and. c == 43) then
!  write(*,*) 'outgoing IBB values',loc_fout,loc_fi,loc_fieq,loc_fi_opp,loc_fieq_opp,&
!    loc_gout,loc_gi,loc_gieq,loc_gi_opp,loc_gieq_opp,i,self
!  write(*,*) 'target values',loc_rho,u_tgt,v_tgt,w_tgt,temp_tgt,i,rest_state(a,b,c,i),dir_state,lvl
!  write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!  write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz
!  write(*,*) 'equilibrium values, constants',loc_fieq,loc_gieq,beta_1,beta_2
!  write(*,*) 'density values',loc_rho
!  write(*,*) 'local u values',loc_u
!  write(*,*) 'local v values',loc_v
!  write(*,*) 'local w values',loc_w
!  write(*,*) 'local temperature',loc_temp
!  write(*,*) 'local fi information',pri_der,p_term_1,p_term_2,p_term_3,q_term_1,&
!    q_term_2,q_term_3,fi_neq*loc_epsilon,loc_epsilon,loc_fout
!  write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*loc_epsilon,loc_epsilon,loc_gout
!
!    write(*,*) 'pab',pab
!    write(*,*) 'qabc',qabc
!    write(*,*) 'k_one',k_one
!    write(*,*) 'g_alpha',g_alpha
!    write(*,*) 'rab',rab
!end if
!
if (a == 349 .and. b == 253 .and. c == 3) then
  write(*,*) 'outgoing IBB values',loc_fout,loc_fi,loc_fieq,loc_fi_opp,loc_fieq_opp,&
    loc_gout,loc_gi,loc_gieq,loc_gi_opp,loc_gieq_opp,a,b,c,i,self
  write(*,*) 'target values',loc_rho,u_tgt,v_tgt,w_tgt,temp_tgt,i,rest_state(a,b,c,i),dir_state,lvl
  write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz
  write(*,*) 'equilibrium values, constants',loc_fieq,loc_gieq,beta_1,beta_2
  write(*,*) 'density values',loc_rho
  write(*,*) 'local u values',loc_u
  write(*,*) 'local v values',loc_v
  write(*,*) 'local w values',loc_w
  write(*,*) 'local temperature',loc_temp
  write(*,*) 'local fi information',pri_der,p_term_1,p_term_2,p_term_3,q_term_1,&
    q_term_2,q_term_3,fi_neq*loc_epsilon,loc_epsilon,loc_fout
  write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*loc_epsilon,loc_epsilon,loc_gout
end if

!if (i == 37) then
!  write(*,*) 'outgoing IBB values',loc_fout,loc_fi,loc_fieq,loc_fi_opp,loc_fieq_opp,&
!    loc_gout,loc_gi,loc_gieq,loc_gi_opp,loc_gieq_opp,a,b,c,i,self
!  write(*,*) 'target values',loc_rho,u_tgt,v_tgt,w_tgt,temp_tgt,i,rest_state(a,b,c,1:39),dir_state,lvl
!  write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!  write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz
!  write(*,*) 'equilibrium values, constants',loc_fieq,loc_gieq,beta_1,beta_2
!!  write(*,*) 'density values',loc_rho
!!  write(*,*) 'local u values',loc_u
!!  write(*,*) 'local v values',loc_v
!!  write(*,*) 'local w values',loc_w
!!  write(*,*) 'local temperature',loc_temp
!  write(*,*) 'local fi information',pri_der,p_term_1,p_term_2,p_term_3,q_term_1,&
!    q_term_2,q_term_3,fi_neq*loc_epsilon,loc_epsilon,loc_fout
!  write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*loc_epsilon,loc_epsilon,loc_gout
!end if

END SUBROUTINE

!if (startup) then
!  call startup_transient_damping(loc_fout, loc_fi,loc_fi_opp, loc_fieq_opp,fout_adj,&
!    loc_gout, loc_gi, loc_gi_opp, loc_gieq_opp,gout_adj,rest_state(a,b,c,i),rest_state(a,b,c,opp(i)),lvl)
!end if

!if (i == 34 .or. i == 37) then
!  write(*,*) 'looking at largest changers in fi',loc_fout,inc_fout,loc_fieq,fi_neq,fi_opp,u_tgt,i
!  write(*,*) 'looking at largest changers in gi',loc_gout,inc_gout,loc_gieq,gi_neq,gi_opp,u_tgt,i
!  write(*,*) ''
!end if

!if (abs(inc_fout/loc_fout) > 3.0D0 .or. abs(inc_gout/loc_gout) > 3.0D0) then
!    write(*,*) 'large out values change',loc_fout,loc_gout,wall_state,rest_state,loc_rho,u_tgt,&
!      v_tgt,w_tgt,temp_tgt,init_epsilon(lvl),init_epsilon(lvl),i,self
!    write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,i
!    write(*,*) 'miscellaneous',beta_1,beta_2,const_vol,k_fi,k_gi,init_epsilon(lvl),init_epsilon(lvl)
!    write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz,i
!    write(*,*) 'local fi information',pri_der,p_term_1,p_term_2,p_term_3,q_term_1,q_term_2,q_term_3,&
!      fi_neq*init_epsilon(lvl),loc_fout,loc_fieq,inc_fout
!    write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*init_epsilon(lvl),loc_gout,loc_gieq,inc_gout
!
!end if


!
!  if (loc_fout > 0.07 .and. i > 7) then
!!    write(*,*) 'fi check after ibb',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),rest_state(a,b,c,opp(i)),&
!!      rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl
!    write(*,*) 'outgoing IBB values, fout large',loc_fout,loc_gout,wall_state,rest_state,loc_rho,u_tgt,&
!      v_tgt,w_tgt,temp_tgt,i,self
!    write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!    write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz
!    write(*,*) 'local fi information',pri_der,fi_neq*init_epsilon(lvl),loc_fout,loc_fieq,inc_fout
!    write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*init_epsilon(lvl),loc_gout,loc_gieq,inc_gout
!!    write(*,*) 'all fis ',fi(a,b,c,1:dir)
!!    write(*,*) 'states after ibb',rest_state(a,b,c,1:dir)
!  end if
!  if (loc_gout > 0.085 .and. i > 7) then
!!    write(*,*) 'gi check after ibb',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),rest_state(a,b,c,opp(i)),&
!!      rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl
!    write(*,*) 'outgoing IBB values, gout large',loc_fout,loc_gout,wall_state,rest_state,loc_rho,u_tgt,&
!      v_tgt,w_tgt,temp_tgt,i,self
!    write(*,*) 'local velocity derivatives',dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!    write(*,*) 'local temperature derivatives',dTdx,dTdy,dTdz
!    write(*,*) 'local fi information',pri_der,fi_neq*init_epsilon(lvl),loc_fout,loc_fieq,inc_fout
!    write(*,*) 'local gi information',k_one,g_alpha,r_ab,gi_neq*init_epsilon(lvl),loc_gout,loc_gieq,inc_gout
!
!!    write(*,*) 'all gis ',gi(a,b,c,1:dir)
!!    write(*,*) 'states after ibb',rest_state(a,b,c,1:dir)
!  end if


!else if (dir == 39) then
!select case(i)
!  case(2)
!! 2 - (+1,0,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-(pri_der)/const_vol)) + & ! term 1
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*dudx/const_vol)) - & ! term 2
!    1/(2*beta_1*temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*(pri_der)) - & ! term 3
!    3/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & ! term 4
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
!            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!! g_missing
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!
!  case(3)
!! 3 - (0,+1,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(pri_der))) + &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy + v_tgt/beta_1*(6*dvdy - 3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*(dudx+dwdz)-cv_const_2*(pri_der)) - &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
!            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(4)
!! 4 - (0,0,+1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(dudx+dvdy+dwdz))) + &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(dudx+dvdy+dwdz)) - &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
!            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(5)
!! 5 - (-1,0,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
!            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(6)
!! 6 - (0,-1,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) - &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
!            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(7)
!! 7 - (0,0,-1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!!
!! f_missing
!  fi_neq = k_fi*(((1-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) - &
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(pri_der)) + &
!    3/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
!            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
!  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(8)
!! +1,+1,+1
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  pri_der = dudx+dvdy+dwdz
!  cross_der = dudy+dvdx+dwdx+dudz+dvdz+dwdy
!! f_missing
!!  fi_neq = k_fi((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!!    6/(2*beta_1*temp_tgt**2)*(dudy+dudz+dvdx+dvdz+dwdx+dwdy) +  & !second term
!!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) + &
!!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + & !u-term
!!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + & !v-term
!!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
!  p_term_1 = (1-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
!  p_term_2 = 6*(dudy+dvdx+dwdx+dudz+dvdz+dwdy)
!  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(dTdx+dTdy+dTdz) + &
!    (u_tgt*(6*dudx+2*(dvdx+dudy+dwdx+dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + &
!     v_tgt*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + &
!     w_tgt*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz))/beta_1)
!  q_term_2 = 4/beta_2*(dTdx+dTdy+dTdz) + &
!    (u_tgt*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!     v_tgt*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!     w_tgt*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))/beta_1
!
!  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2)/(6*temp_tgt**3))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt+w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy+dTdz)+ &!term 2, 1/2T^2
!            v_tgt*(dTdx+dTdz)+ w_tgt*(dTdx+dTdy)) - &
!            2*pri_der/const_vol*(u_tgt*v_tgt+u_tgt*w_tgt+v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(9)
!! 9 - (-1,+1,+1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = -dudy-dudz-dvdx+dvdz-dwdx+dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!!    6/(2*beta_1*temp_tgt**2)*(-dudy-dudz-dvdx+dvdz-dwdx+dwdy) +  & !second term
!!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdz+dTdy-dTdx) + &
!!    u_tgt/beta_1*(-6*dudx+2*(dvdx+dudy+dwdx+dudz-cv_const*(dvdy+dwdz))+5*cv_const*dudx) + & !u-term
!!    v_tgt/beta_1*(6*dvdy+2*(-dvdx-dudy+dwdy+dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + & !v-term
!!    w_tgt/beta_1*(6*dwdz+2*(-dudz-dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!!    u_tgt/beta_1*(4*(-dvdy-dwdz+dudy+dvdx+dwdx+dudz+cv_const*dudx)-6*(dwdy+dvdz)) + &
!!    v_tgt/beta_1*(4*(dudx+dwdz-dvdx-dudy+dwdy+dvdz-cv_const*dvdy)-6*(dudz+dwdx)) + &
!!    w_tgt/beta_1*(4*(dudx+dvdy-dwdx+dwdy-dudz+dvdz-cv_const*dwdz)-6*(dudy+dvdx)))) !fourth term
!  p_term_1 = (1-temp_tgt)*(6*(pri_der)-3/cv_const*pri_der)
!  p_term_2 = 6*cross_der
!  q_term_1 = (1-3*temp_tgt)*( 5/beta_2*(dTdx+dTdy+dTdz) + &
!    (u_tgt*(6*dudx+2*(dvdx+dudy+dwdx+dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + &
!     v_tgt*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + &
!     w_tgt*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz))/beta_1)
!  q_term_2 = 4/beta_2*(dTdx+dTdy+dTdz) + &
!    (u_tgt*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!     v_tgt*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!     w_tgt*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))/beta_1
!
!  fi_neq = k_fi*((p_term_1 + p_term_2)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2)/(6*temp_tgt**3))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(-dTdx+dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(-u_tgt+v_tgt+w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)-u_tgt*(dTdy+dTdz)+ & !term 2, 1/2T^2
!            v_tgt*(-dTdx+dTdz)+ w_tgt*(-dTdx+dTdy)) - &
!            2*pri_der/const_vol*(-u_tgt*v_tgt-u_tgt*w_tgt+v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(10)
!! 10 - (-1,-1,+1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = dudy-dudz+dvdx-dvdz-dwdx-dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(dudy-dudz+dvdx-dvdz-dwdx-dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(-dTdx-dTdy+dTdz) + &
!    u_tgt/beta_1*(-6*dudx+2*(dvdx+dudy+dwdx+dudz-cv_const*(dvdy+dwdz))+5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(-6*dvdy+2*(dvdx+dudy+dwdy+dvdz-cv_const*(dudx+dwdz))+5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz+cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz+cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(-dTdx-dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(-u_tgt-v_tgt+w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(-dTdy+dTdz) + &!term 2, 1/2T^2
!            v_tgt*(-dTdx+dTdz) - w_tgt*(dTdx+dTdy)) - &
!            2*pri_der/const_vol*(u_tgt*v_tgt-u_tgt*w_tgt-v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(11)
!! 11 - (+1,-1,+1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = -dudy+dudz-dvdx-dvdz+dwdx-dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(-dudy+dudz-dvdx-dvdz+dwdx-dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx-dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(-dvdx-dudy+dwdx+dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(-6*dvdy+2*(dvdx+dudy+dwdy+dvdz-cv_const*(dudx+dwdz))+5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx-dvdz-dwdy+cv_const*(dudx+dvdy))-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx-dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz-dudy-dvdx+dwdx+dudz-cv_const*dudx)-6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(-dudx-dwdz+dvdx+dudy+dwdy+dvdz+cv_const*dvdy)-6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy-dwdx+dwdy+dudz-dvdz-cv_const*dwdz)-6*(dudy+dvdx)))) !fourth term
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx-dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt+w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(-dTdy+dTdz) - &!term 2, 1/2T^2
!            v_tgt*(dTdx+dTdz)+ w_tgt*(dTdx-dTdy)) - &
!            2*pri_der/const_vol*(-u_tgt*v_tgt+u_tgt*w_tgt-v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(12)
!! 12 - (+1,+1,-1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = dudy-dudz+dvdx-dvdz-dwdx-dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(dudy-dudz+dvdx-dvdz-dwdx-dwdy) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy-dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy-dwdx-dudz+cv_const*(dvdy+dwdz))-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy-dwdy-dvdz+cv_const*(dudx+dwdz))-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(-6*dwdz+2*(dudz+dwdx+dvdz+dwdy-cv_const*(dudx+dvdy))+5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy-dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx-dwdx-dudz-cv_const*dudx)-6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy-dwdy-dvdz-cv_const*dvdy)-6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(-dudx-dvdy+dwdx+dwdy+dudz+dvdz+cv_const*dwdz)-6*(dudy+dvdx)))) !fourth term
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt-w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy-dTdz)+ &!term 2, 1/2T^2
!            v_tgt*(dTdx-dTdz) - w_tgt*(dTdx+dTdy)) - &
!            2*pri_der/const_vol*(u_tgt*v_tgt - u_tgt*w_tgt - v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(13)
!! 13 - (-1,+1,-1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = -dudy+dudz-dvdx-dvdz+dwdx-dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(cross_der) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(-dTdx+dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(-u_tgt+v_tgt-w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy-dTdz) - &!term 2, 1/2T^2
!            v_tgt*(dTdx+dTdz)+ w_tgt*(-dTdx+dTdy)) - &
!            2*pri_der/const_vol*(-u_tgt*v_tgt+u_tgt*w_tgt-v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(14)
!! 14 - (-1,-1,-1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = dudy+dudz+dvdx+dvdz+dwdx+dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der)) + & !first term
!    6/(2*beta_1*temp_tgt**2)*(cross_der) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) - &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz) - &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt+w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) - &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der)+u_tgt*(dTdy+dTdz)+ &!term 2, 1/2T^2
!            v_tgt*(dTdx+dTdz)+ w_tgt*(dTdx+dTdy)) - &
!            2*pri_der/const_vol*(u_tgt*v_tgt+u_tgt*w_tgt+v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one - g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(15)
!! 15 - (+1,-1,-1)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!
!  cross_der = -dudy-dudz-dvdx+dvdz-dwdx+dwdy
!  pri_der = dudx+dvdy+dwdz
!
!! f_missing
!  fi_neq = k_fi*((1-temp_tgt)/(2*beta_1*temp_tgt**2)*(6*(pri_der)-3/cv_const*(pri_der))+ & !first term
!    6/(2*beta_1*temp_tgt**2)*(cross_der) +  & !second term
!    (1-3*temp_tgt)/(6*temp_tgt**3)*(5/beta_2*(dTdx+dTdy+dTdz) + &
!    u_tgt/beta_1*(6*dudx+2*(dvdx+dudy+dwdx+dudz+dvdy+dwdz)-5*cv_const*dudx) + & !u-term
!    v_tgt/beta_1*(6*dvdy+2*(dvdx+dudy+dwdy+dvdz+dudx+dwdz)-5*cv_const*dvdy) + & !v-term
!    w_tgt/beta_1*(6*dwdz+2*(dudz+dwdx+dvdz+dwdy+dudx+dvdy)-5*cv_const*dwdz)) + & !w-term, third term
!    1/(6*temp_tgt**3)*(4/beta_2*(dTdx+dTdy+dTdz)+ &
!    u_tgt/beta_1*(4*(dvdy+dwdz+dudy+dvdx+dwdx+dudz-cv_const*dudx)+6*(dwdy+dvdz)) + &
!    v_tgt/beta_1*(4*(dudx+dwdz+dvdx+dudy+dwdy+dvdz-cv_const*dvdy)+6*(dudz+dwdx)) + &
!    w_tgt/beta_1*(4*(dudx+dvdy+dwdx+dwdy+dudz+dvdz-cv_const*dwdz)+6*(dudy+dvdx)))) !fourth term
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx-dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt-w_tgt)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy+dwdz)+u_tgt*dTdx+v_tgt*dTdy+w_tgt*dTdz) - & !term 1, 1-T/2T^2
!            pri_der*(6*temp_tgt+u_tgt**2+v_tgt**2+w_tgt**2)/const_vol) + &
!            1/(2*temp_tgt**2)*(6*(temp_tgt*(cross_der) - u_tgt*(dTdy+dTdz)+ &!term 2, 1/2T^2
!            v_tgt*(-dTdx+dTdz)+ w_tgt*(-dTdx+dTdy)) - &
!            2*pri_der/const_vol*(-u_tgt*v_tgt-u_tgt*w_tgt+v_tgt*w_tgt)) !
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(16)
!! 16 - (+2,0,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) + & !first term
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - & !third term
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) - & !second term
!    1/(temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & !fourth term
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
!            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(17)
!! 17 - (0,+2,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) + & !first term
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - & !third term
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) - & !second term
!    1/(temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & !fourth term
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (1-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
!            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt + u_tgt**2 + w_tgt**2))
!  gi_neq = k_gi*(k_one + g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(18)
!! 18 - (0,0,+2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(dudx+dvdy+dwdz))) + &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(dudx+dvdy+dwdz)) - &
!    1/(temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
!            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(19)
!! 19 - (-2,0,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    1/(temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
!            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(20)
!! 20 - (0,-2,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(dudx+dvdy+dwdz))) - &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    1/(temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
!            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(21)
!! 21 - (0,0,-2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((4-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) - &
!    (8-6*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(pri_der)) + &
!    1/(temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
!            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
!  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(22)
!! 22 - (+2,+2,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudy+dvdx)
!  p_term_3 = 6*dwdz-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx+dTdy)+(u_tgt*(6*dudx+2*(dvdx+dudy+dvdy)-4*cv_const*dudx)+&
!    v_tgt*(6*dvdy+2*(dvdx+dudy+dudx)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx+dTdy)+(u_tgt*(4*(dudy+dvdx+dvdy)-cv_const_2*dudx)+&
!    v_tgt*(4*(dudy+dvdx+dudx)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx+dTdy)/beta_2 + (2.0D0*w_tgt*(dudz+dwdx+dvdz+dwdy) + &
!    u_tgt*(2*dwdz-cv_const*dudx) + v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) + &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
!            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(23)
!! 23 - (+2,0,+2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudz+dwdx)
!  p_term_3 = 6*dvdy-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx+dTdz)+(u_tgt*(6*dudx+2*(dwdx+dudz+dwdz)-4*cv_const*dudx)+&
!    w_tgt*(6*dwdz+2*(dwdx+dudz+dudx)-4*cv_const*dwdz))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx+dTdz)+(u_tgt*(4*(dudz+dwdx+dwdz)-cv_const_2*dudx)+&
!    w_tgt*(4*(dudz+dwdx+dudx)-cv_const_2*dwdz))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx+dvdz+dwdy) + &
!    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+w_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
!            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) + &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
!            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(24)
!! 24 - (0,+2,+2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dwdy+dvdz)
!  p_term_3 = 6*dudx-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdz+dTdy)+(w_tgt*(6*dwdz+2*(dvdz+dwdy+dvdy)-4*cv_const*dwdz)+&
!    v_tgt*(6*dvdy+2*(dvdz+dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdz+dTdy)+(w_tgt*(4*(dwdy+dvdz+dvdy)-cv_const_2*dwdz)+&
!    v_tgt*(4*(dwdy+dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdz+dTdy)/beta_2 + (2.0D0*u_tgt*(dudz+dwdx+dvdx+dudy) + &
!    w_tgt*(2*dudx-cv_const*dwdz) + v_tgt*(2*dudx-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt+v_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) + &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
!            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(25)
!! 25 - (-2,+2,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudy+dvdx)
!  p_term_3 = 6*dwdz-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(-dTdx+dTdy)+(-u_tgt*(6*dudx+2*(-dvdx-dudy+dvdy)-4*cv_const*dudx)+&
!    v_tgt*(6*dvdy+2*(-dvdx-dudy+dudx)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(-dTdx+dTdy)+(-u_tgt*(4*(-dudy-dvdx+dvdy)-cv_const_2*dudx)+&
!    v_tgt*(4*(-dudy-dvdx+dudx)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((-dTdx+dTdy)/beta_2 + (2.0D0*w_tgt*(-dudz-dwdx+dvdz+dwdy) - &
!    u_tgt*(2*dwdz-cv_const*dudx) + v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdy-dTdx)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(v_tgt-u_tgt)
!  r_ab = (4-temp_tgt)/(4*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) - &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
!            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!!
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(26)
!! 26 - (-2,-2,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudy+dvdx)
!  p_term_3 = 6*dwdz-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx+dTdy)+(u_tgt*(6*dudx+2*(dvdx+dudy+dvdy)-4*cv_const*dudx)+&
!    v_tgt*(6*dvdy+2*(dvdx+dudy+dudx)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx+dTdy)+(u_tgt*(4*(dudy+dvdx+dvdy)-cv_const_2*dudx)+&
!    v_tgt*(4*(dudy+dvdx+dudx)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx+dTdy)/beta_2 + (2.0D0*w_tgt*(dudz+dwdx+dvdz+dwdy) + &
!    u_tgt*(2*dwdz-cv_const*dudx) + v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+v_tgt)
!  r_ab = (4-temp_tgt)/(4*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) + &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
!            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(27)
!! 27 - (+2,-2,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudy+dvdx)
!  p_term_3 = 6*dwdz-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx-dTdy)+(u_tgt*(6*dudx+2*(-dvdx-dudy+dvdy)-4*cv_const*dudx)-&
!    v_tgt*(6*dvdy+2*(-dvdx-dudy+dudx)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx-dTdy)+(u_tgt*(4*(-dudy-dvdx+dvdy)-cv_const_2*dudx)-&
!    v_tgt*(4*(-dudy-dvdx+dudx)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx-dTdy)/beta_2 + (2.0D0*w_tgt*(dudz+dwdx-dvdz-dwdy) + &
!    u_tgt*(2*dwdz-cv_const*dudx) - v_tgt*(2*dwdz-cv_const*dvdy))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx-dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-v_tgt)
!  r_ab = (4-temp_tgt)/(4*temp_tgt**2)*(6*(temp_tgt*(dudx+dvdy)+u_tgt*dTdx+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+u_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*dwdz+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+w_tgt**2)) - &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudy+dvdx)+u_tgt*dTdy+v_tgt*dTdx) - &
!            (2*u_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(28)
!! 28 - (+2,0,-2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudz+dwdx)
!  p_term_3 = 6*dvdy-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx-dTdz)+(u_tgt*(6*dudx+2*(-dwdx-dudz+dwdz)-4*cv_const*dudx) - &
!    w_tgt*(6*dwdz+2*(-dwdx-dudz+dudx)-4*cv_const*dwdz))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx-dTdz)+(u_tgt*(4*(-dudz-dwdx+dwdz)-cv_const_2*dudx) - &
!    w_tgt*(4*(-dudz-dwdx+dudx)-cv_const_2*dwdz))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx-dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx-dvdz-dwdy) + &
!    u_tgt*(2*dvdy-cv_const*dudx) - w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt-w_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
!            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) - &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
!            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(29)
!! 29 - (-2,0,-2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudz+dwdx)
!  p_term_3 = 6*dvdy-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdx+dTdz)+(u_tgt*(6*dudx+2*(dwdx+dudz+dwdz)-4*cv_const*dudx)+&
!    w_tgt*(6*dwdz+2*(dwdx+dudz+dudx)-4*cv_const*dwdz))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdx+dTdz)+(u_tgt*(4*(dudz+dwdx+dwdz)-cv_const_2*dudx)+&
!    w_tgt*(4*(dudz+dwdx+dudx)-cv_const_2*dwdz))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(dudy+dvdx+dvdz+dwdy) + &
!    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdx+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(u_tgt+w_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
!            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) + &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
!            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(30)
!! 30 - (-2,0,+2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dudx+dwdz)-cv_const_2*pri_der)
!  p_term_2 = 24*(dudz+dwdx)
!  p_term_3 = 6*dvdy-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(-dTdx+dTdz)+(-u_tgt*(6*dudx+2*(-dwdx-dudz+dwdz)-4*cv_const*dudx)+&
!    w_tgt*(6*dwdz+2*(-dwdx-dudz+dudx)-4*cv_const*dwdz))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(-dTdx+dTdz)+(-u_tgt*(4*(-dudz-dwdx+dwdz)-cv_const_2*dudx)+&
!    w_tgt*(4*(-dudz-dwdx+dudx)-cv_const_2*dwdz))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((-dTdx+dTdz)/beta_2 + (2.0D0*v_tgt*(-dudy-dvdx+dvdz+dwdy) - &
!    u_tgt*(2*dvdy-cv_const*dudx) + w_tgt*(2*dvdy-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdz-dTdx)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt-u_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 1
!            pri_der*(2*temp_tgt+u_tgt**2+w_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy)+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+v_tgt**2)) - &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dudz+dwdx)+u_tgt*dTdz+w_tgt*dTdx) - &
!            (2*u_tgt*w_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(31)
!! 31 - (0,+2,-2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dwdy+dvdz)
!  p_term_3 = 6*dudx-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(-dTdz+dTdy)+(-w_tgt*(6*dwdz+2*(-dvdz-dwdy+dvdy)-4*cv_const*dwdz)+&
!    v_tgt*(6*dvdy+2*(-dvdz-dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(-dTdz+dTdy)+(-w_tgt*(4*(-dwdy-dvdz+dvdy)-cv_const_2*dwdz)+&
!    v_tgt*(4*(-dwdy-dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((-dTdz+dTdy)/beta_2 + (2.0D0*u_tgt*(-dudz-dwdx+dvdx+dudy) - &
!    w_tgt*(2*dudx-cv_const*dwdz) + v_tgt*(2*dudx-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(const_vol*beta_1)
!  g_alpha = 9*(dTdy-dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(v_tgt-w_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) - &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
!            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(32)
!! 32 - (0,-2,-2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dwdy+dvdz)
!  p_term_3 = 6*dudx-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdz+dTdy)+(w_tgt*(6*dwdz+2*(dvdz+dwdy+dvdy)-4*cv_const*dwdz)+&
!    v_tgt*(6*dvdy+2*(dvdz+dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdz+dTdy)+(w_tgt*(4*(dwdy+dvdz+dvdy)-cv_const_2*dwdz)+&
!    v_tgt*(4*(dwdy+dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdz+dTdy)/beta_2 + (2.0D0*u_tgt*(dudz+dwdx+dvdx+dudy) + &
!    w_tgt*(2*dudx-cv_const*dwdz) + v_tgt*(2*dudx-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 + p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) - &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdy+dTdz)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt+v_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx)+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) + &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
!            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one - 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(33)
!! 33 - (0,-2,+2)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  p_term_1 = (4.0D0-temp_tgt)*(6*(dwdz+dvdy)-cv_const_2*pri_der)
!  p_term_2 = 24*(dwdy+dvdz)
!  p_term_3 = 6*dudx-cv_const*pri_der
!  q_term_1 = (8.0D0-6.0D0*temp_tgt)*(4*(dTdz-dTdy)+(w_tgt*(6*dwdz+2*(-dvdz-dwdy+dvdy)-4*cv_const*dwdz) - &
!    v_tgt*(6*dvdy+2*(-dvdz-dwdy+dwdz)-4*cv_const*dvdy))/beta_1)
!  q_term_2 = 8.0D0*(2.0D0/beta_2*(dTdz-dTdy)+(w_tgt*(4*(-dwdy-dvdz+dvdy)-cv_const_2*dwdz) - &
!    v_tgt*(4*(-dwdy-dvdz+dwdz)-cv_const_2*dvdy))/beta_1)
!  q_term_3 = (-6.0D0*temp_tgt)*((dTdz-dTdy)/beta_2 + (2.0D0*u_tgt*(dudz+dwdx-dvdx-dudy) + &
!    w_tgt*(2*dudx-cv_const*dwdz) - v_tgt*(2*dudx-cv_const*dwdz))/beta_1)
!
!  fi_neq = k_fi*((p_term_1 - p_term_2 - p_term_3)/(2*beta_1*temp_tgt**2) + &
!    (q_term_1 + q_term_2 + q_term_3)/(6*temp_tgt**3))
!
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*(dTdz-dTdy)/(beta_2*temp_tgt)+3*pri_der/(beta_1*temp_tgt*const_vol)*(w_tgt-v_tgt)
!  r_ab = (4-temp_tgt)/(2*temp_tgt**2)*(6*(temp_tgt*(dwdz+dvdy)+w_tgt*dTdz+v_tgt*dTdy) - & !term 1
!            pri_der*(4*temp_tgt+w_tgt**2+v_tgt**2)/const_vol) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*dudx+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(2*temp_tgt+u_tgt**2)) - &
!            4/(2*temp_tgt**2)*(6*(temp_tgt*(dwdy+dvdz)+w_tgt*dTdy+v_tgt*dTdz) - &
!            (2*w_tgt*v_tgt)*pri_der/const_vol) !term 3
!  gi_neq = k_gi*(k_one + 2*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(34)
!! 34 - (+3,0,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-(pri_der)/const_vol)) + & ! term 1
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*dudx/const_vol)) - & ! term 2
!    1/(2*beta_1*temp_tgt)*(6*(dvdy+dwdz)-cv_const_2*(pri_der)) - & ! term 3
!    9/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* & ! term 4
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
!            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt + v_tgt**2 + w_tgt**2))
!  gi_neq = k_gi*(k_one + 3*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(35)
!! 35 - (0,+3,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(pri_der))) + &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy + v_tgt/beta_1*(6*dvdy - 3*cv_const*dvdy)) - &
!    1/(2*beta_1*temp_tgt)*(6*(dudx+dwdz)-cv_const_2*(pri_der)) - &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!
!! g
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
!            (2*temp_tgt+v_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one + 3*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(36)
!! 36 - (0,0,+3)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) + &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz + w_tgt/beta_1*(6*dwdz - 3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*(dudx+dvdy)-cv_const_2*(pri_der)) - &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dudx+dvdy-cv_const*dwdz)+2*u_tgt/beta_1* &
!    (dudz+dwdx)+2*v_tgt/beta_1*(dvdz+dwdy)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
!            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+u_tgt*dTdx+v_tgt*dTdy) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2+u_tgt**2))
!  gi_neq = k_gi*(k_one + 3*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(37)
!! 37 - (-3,0,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dudx-cv_const*(dudx+dvdy+dwdz))) - &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdx+u_tgt/beta_1*(6*dudx-3*cv_const*dudx)) - &
!    1/(2*beta_1*temp_tgt)*(6*dvdy+6*dwdz-cv_const_2*(dudx+dvdy+dwdz)) + &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdx + 2*u_tgt/beta_1*(dvdy+dwdz-cv_const*dudx)+2*v_tgt/beta_1* &
!    (dudy+dvdx)+2*w_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdx/(beta_2*temp_tgt)+3*u_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dudx + 6*u_tgt*dTdx - & !term 1
!            (2*temp_tgt+u_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dwdz)+v_tgt*dTdy+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one - 3*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!  case(38)
!! 38 - (0,-3,0)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dvdy-cv_const*(pri_der))) - & !P term 1
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdy+v_tgt/beta_1*(6*dvdy-3*cv_const*dvdy)) - & !Q term 1
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dwdz-cv_const_2*(pri_der)) + & !P term 2
!    9/(6*temp_tgt**2)*(2/beta_2*dTdy + 2*v_tgt/beta_1*(dudx+dwdz-cv_const*dvdy)+2*u_tgt/beta_1* & !Q term 2
!    (dudy+dvdx)+2*w_tgt/beta_1*(dvdz+dwdy)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdy/(beta_2*temp_tgt)+3*v_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dvdy + 6*v_tgt*dTdy - & !term 1
!            (2*temp_tgt+v_tgt**2)/const_vol*3*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dudx+dwdz)+u_tgt*dTdx+w_tgt*dTdz) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +u_tgt**2+w_tgt**2))
!  gi_neq = k_gi*(k_one - 3*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  case(39)
!! 39 - (0,0,-3)
!  dudx = dudx_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudy = dudy_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dudz = dudz_ibb(u_tgt,0.0D0,upwind,dir_state,i)
!  dvdx = dvdx_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdy = dvdy_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dvdz = dvdz_ibb(v_tgt,0.0D0,upwind,dir_state,i)
!  dwdx = dwdx_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdy = dwdy_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dwdz = dwdz_ibb(w_tgt,0.0D0,upwind,dir_state,i)
!  dTdx = dTdx_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdy = dTdy_ibb(temp_tgt,temperature,upwind,dir_state,i)
!  dTdz = dTdz_ibb(temp_tgt,temperature,upwind,dir_state,i)
!!
!  pri_der = dudx+dvdy+dwdz
!! f
!  fi_neq = k_fi*(((9-temp_tgt)/(2*temp_tgt**2*beta_1)*(6*dwdz-cv_const*(pri_der))) - &
!    (27-9*temp_tgt)/(6*temp_tgt**3)*(3/beta_2*dTdz+w_tgt/beta_1*(6*dwdz-3*cv_const*dwdz)) - &
!    1/(2*beta_1*temp_tgt)*(6*dudx+6*dvdy-cv_const_2*(pri_der)) + &
!    9/(6*temp_tgt**2)*(2/beta_2*dTdz + 2*w_tgt/beta_1*(dvdy+dudx-cv_const*dwdz)+2*v_tgt/beta_1* &
!    (dvdz+dwdy)+2*u_tgt/beta_1*(dudz+dwdx)))
!  loc_fout = loc_fieq + init_epsilon(lvl)*fi_neq
!
!!
!! g missing
!!
!  k_one = 9*pri_der/(beta_1*const_vol)
!  g_alpha = 9*dTdz/(beta_2*temp_tgt)+3*w_tgt*pri_der/(beta_1*temp_tgt*const_vol)
!  r_ab = (9-temp_tgt)/(2*temp_tgt**2)*(6*temp_tgt*dwdz + 6*w_tgt*dTdz - & !term 1
!            (2*temp_tgt+w_tgt**2)/const_vol*pri_der) - &
!            1/(2*temp_tgt)*(6*(temp_tgt*(dvdy+dudx)+v_tgt*dTdy+u_tgt*dTdx) - & !term 2
!            pri_der/const_vol*(4*temp_tgt +v_tgt**2 + u_tgt**2))
!  gi_neq = k_gi*(k_one - 3*g_alpha + r_ab/beta_1)
!
!  loc_gout = loc_gieq + init_epsilon(lvl)*gi_neq
!
!  end select








!
! For the case where the distance to the wall is more than half a lattice length
!
!      IF (wall_dist(wall_node_number) >= 0.5D0) THEN
!        fout_bc= (1.0D0/(2.0D0*wall_dist(wall_node_number))) * fout_bc +&
!          ((2.0D0*wall_dist(wall_node_number)-1.0D0) / (2.0D0*&
!          wall_dist(wall_node_number)))*fout(opp_dir,loc_dir)
!
!!          WRITE(*,*) fout_bc,wall_dist(wall_node_number),loc_dir
!!
!! For the case where the distance to the wall is less than half a lattice length
!!
!      ELSE IF (wall_dist(wall_node_number) < 0.5D0) THEN
!!
!! Test if the node is in a corner. If it is, set it equal to its opposite fout
!! and treat it like a standard bounceback
!!
!        upstream_node = link(opp_dir,loc_node)
!        IF (upstream_node < 0) THEN
!!          WRITE(*,*) fout_bc,wall_dist(wall_node_number),loc_dir,upstream_node,fout_opp
!          fout_bc = fout_opp
!!          WRITE(*,*) fout_bc,wall_dist(wall_node_number),loc_dir,upstream_node,fout_opp
!
!        ELSE
!!          WRITE(*,*) fout_bc,fout(loc_dir,upstream_node),upstream_node
!          fout_bc = 2.0D0*wall_dist(wall_node_number)*fout_bc + (1.0D0-2.0D0 &
!          *wall_dist(wall_node_number))*fout(loc_dir,upstream_node)
!
!!          WRITE(*,*) fout_bc,wall_dist(wall_node_number),loc_dir,fout(loc_dir,upstream_node)
!        END IF
!
!      END IF
