module derivatives
!
! Derivatives of all stripes!!!
!
!
!
!
use precise
use constants
use linkwise
use grid_data
implicit none
save
private
public :: dudx_upwind_first,dudy_upwind_first,dudz_upwind_first,&
  dvdx_upwind_first,dvdy_upwind_first,dvdz_upwind_first,dwdx_upwind_first,&
  dwdy_upwind_first,dwdz_upwind_first,dTdx_upwind_first,dTdy_upwind_first,&
  dTdz_upwind_first
public :: dudx_downwind_first,dudy_downwind_first,dudz_downwind_first,&
  dvdx_downwind_first,dvdy_downwind_first,dvdz_downwind_first,dwdx_downwind_first,&
  dwdy_downwind_first,dwdz_downwind_first,dTdx_downwind_first,dTdy_downwind_first,&
  dTdz_downwind_first
public :: dudx_central_second,dudy_central_second,dudz_central_second,&
  dvdx_central_second,dvdy_central_second,dvdz_central_second,dwdx_central_second,&
  dwdy_central_second,dwdz_central_second,dTdx_central_second,dTdy_central_second,&
  dTdz_central_second
public :: dudx_ibb,dudy_ibb,dudz_ibb,&
  dvdx_ibb,dvdy_ibb,dvdz_ibb,dwdx_ibb,&
  dwdy_ibb,dwdz_ibb,dTdx_ibb,dTdy_ibb,&
  dTdz_ibb
public :: first_order_vel_cond_unknown,first_order_temp_cond_unknown
public :: first_order_shifted_vel_cond_unknown,first_order_shifted_temp_cond_unknown

public :: dudx_local,dudy_local,dudz_local,dvdx_local,dvdy_local,dvdz_local,&
  dwdx_local,dwdy_local,dwdz_local,dTdx_local,dTdy_local,dTdz_local
public :: derivative_windage_secondord,derivative_windage_firstord,&
  derivative_windage_highord
public :: high_order_known,second_order_known

integer :: order_of_acc


contains

function dudx_upwind_first(u_loc,u_upwind) result (dudx)
  real(kind=dp) :: u_loc,u_upwind,dudx
  dudx = (u_upwind-u_loc)
end function

function dudy_upwind_first(u_loc,u_upwind) result (dudy)
  real(kind=dp) :: u_loc,u_upwind,dudy
  dudy = (u_upwind-u_loc)
end function

function dudz_upwind_first(u_loc,u_upwind) result (dudz)
  real(kind=dp) :: u_loc,u_upwind,dudz
  dudz = (u_upwind-u_loc)
end function

function dvdx_upwind_first(v_loc,v_upwind) result (dvdx)
  real(kind=dp) :: v_loc,v_upwind,dvdx
  dvdx = (v_upwind-v_loc)
end function

function dvdy_upwind_first(v_loc,v_upwind) result (dvdy)
  real(kind=dp) :: v_loc,v_upwind,dvdy
  dvdy = (v_upwind-v_loc)
end function

function dvdz_upwind_first(v_loc,v_upwind) result (dvdz)
  real(kind=dp) :: v_loc,v_upwind,dvdz
  dvdz = (v_upwind-v_loc)
end function

function dwdx_upwind_first(w_loc,w_upwind) result (dwdx)
  real(kind=dp) :: w_loc,w_upwind,dwdx
  dwdx = (w_upwind-w_loc)
end function

function dwdy_upwind_first(w_loc,w_upwind) result (dwdy)
  real(kind=dp) :: w_loc,w_upwind,dwdy
  dwdy = (w_upwind-w_loc)
end function

function dwdz_upwind_first(w_loc,w_upwind) result (dwdz)
  real(kind=dp) :: w_loc,w_upwind,dwdz
  dwdz = (w_upwind-w_loc)
end function

function dTdx_upwind_first(T_loc,T_upwind) result (dTdx)
  real(kind=dp) :: T_loc,T_upwind,dTdx
  dTdx = (T_upwind-T_loc)
end function

function dTdy_upwind_first(T_loc,T_upwind) result (dTdy)
  real(kind=dp) :: T_loc,T_upwind,dTdy
  dTdy = (T_upwind-T_loc)
end function

function dTdz_upwind_first(T_loc,T_upwind) result (dTdz)
  real(kind=dp) :: T_loc,T_upwind,dTdz
  dTdz = (T_upwind-T_loc)
end function
!
!
! Downwind derivatives
!
!
function dudx_downwind_first(u_loc,u_downwind) result (dudx)
  real(kind=dp) :: u_loc,u_downwind,dudx
  dudx = (u_loc-u_downwind)
end function

function dudy_downwind_first(u_loc,u_downwind) result (dudy)
  real(kind=dp) :: u_loc,u_downwind,dudy
  dudy = (u_loc-u_downwind)
end function

function dudz_downwind_first(u_loc,u_downwind) result (dudz)
  real(kind=dp) :: u_loc,u_downwind,dudz
  dudz = (u_loc-u_downwind)
end function

function dvdx_downwind_first(v_loc,v_downwind) result (dvdx)
  real(kind=dp) :: v_loc,v_downwind,dvdx
  dvdx = (v_loc-v_downwind)
end function

function dvdy_downwind_first(v_loc,v_downwind) result (dvdy)
  real(kind=dp) :: v_loc,v_downwind,dvdy
  dvdy = (v_loc-v_downwind)
end function

function dvdz_downwind_first(v_loc,v_downwind) result (dvdz)
  real(kind=dp) :: v_loc,v_downwind,dvdz
  dvdz = (v_loc-v_downwind)
end function

function dwdx_downwind_first(w_loc,w_downwind) result (dwdx)
  real(kind=dp) :: w_loc,w_downwind,dwdx
  dwdx = (w_loc-w_downwind)
end function

function dwdy_downwind_first(w_loc,w_downwind) result (dwdy)
  real(kind=dp) :: w_loc,w_downwind,dwdy
  dwdy = (w_loc-w_downwind)
end function

function dwdz_downwind_first(w_loc,w_downwind) result (dwdz)
  real(kind=dp) :: w_loc,w_downwind,dwdz
  dwdz = (w_loc-w_downwind)
end function

function dTdx_downwind_first(T_loc,T_downwind) result (dTdx)
  real(kind=dp) :: T_loc,T_downwind,dTdx
  dTdx = (T_loc-T_downwind)
end function

function dTdy_downwind_first(T_loc,T_downwind) result (dTdy)
  real(kind=dp) :: T_loc,T_downwind,dTdy
  dTdy = (T_loc-T_downwind)
end function

function dTdz_downwind_first(T_loc,T_downwind) result (dTdz)
  real(kind=dp) :: T_loc,T_downwind,dTdz
  dTdz = (T_loc-T_downwind)
end function
!
!
! Central Difference
!
!
function dudx_central_second(u_upwind,u_downwind) result (dudx)
  real(kind=dp) :: u_upwind,u_downwind,dudx
  dudx = (u_upwind-u_downwind)/2
end function

function dudy_central_second(u_upwind,u_downwind) result (dudy)
  real(kind=dp) :: u_upwind,u_downwind,dudy
  dudy = (u_upwind-u_downwind)/2
end function

function dudz_central_second(u_upwind,u_downwind) result (dudz)
  real(kind=dp) :: u_upwind,u_downwind,dudz
  dudz = (u_upwind-u_downwind)/2
end function

function dvdx_central_second(v_upwind,v_downwind) result (dvdx)
  real(kind=dp) :: v_upwind,v_downwind,dvdx
  dvdx = (v_upwind-v_downwind)/2
end function

function dvdy_central_second(v_upwind,v_downwind) result (dvdy)
  real(kind=dp) :: v_upwind,v_downwind,dvdy
  dvdy = (v_upwind-v_downwind)/2
end function

function dvdz_central_second(v_upwind,v_downwind) result (dvdz)
  real(kind=dp) :: v_upwind,v_downwind,dvdz
  dvdz = (v_upwind-v_downwind)/2
end function

function dwdx_central_second(w_upwind,w_downwind) result (dwdx)
  real(kind=dp) :: w_upwind,w_downwind,dwdx
  dwdx = (w_upwind-w_downwind)/2
end function

function dwdy_central_second(w_upwind,w_downwind) result (dwdy)
  real(kind=dp) :: w_upwind,w_downwind,dwdy
  dwdy = (w_upwind-w_downwind)/2
end function

function dwdz_central_second(w_upwind,w_downwind) result (dwdz)
  real(kind=dp) :: w_upwind,w_downwind,dwdz
  dwdz = (w_upwind-w_downwind)/2
end function

function dTdx_central_second(T_upwind,T_downwind) result (dTdx)
  real(kind=dp) :: T_upwind,T_downwind,dTdx
  dTdx = (T_upwind-T_downwind)/2
end function

function dTdy_central_second(T_upwind,T_downwind) result (dTdy)
  real(kind=dp) :: T_upwind,T_downwind,dTdy
  dTdy = (T_upwind-T_downwind)/2
end function

function dTdz_central_second(T_upwind,T_downwind) result (dTdz)
  real(kind=dp) :: T_upwind,T_downwind,dTdz
  dTdz = (T_upwind-T_downwind)/2
end function
!
!
! IBB semi-central derivatives
!
!
function dudx_ibb(u_opp,u_wall,up_der,wall_dist,i) result (dudx)
  real(kind=dp) :: u_opp,u_wall,wall_dist,dudx
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dudx = (u_opp-u_wall)/(1.0D0+wall_dist)
  else
    dudx = (u_wall-u_opp)/(1.0D0+wall_dist)
  end if
end function

function dudy_ibb(u_opp,u_wall,up_der,wall_dist,i) result (dudy)
  real(kind=dp) :: u_opp,u_wall,wall_dist,dudy
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dudy = (u_opp-u_wall)/(1.0D0+wall_dist)
  else
    dudy = (u_wall-u_opp)/(1.0D0+wall_dist)
  end if
end function

function dudz_ibb(u_opp,u_wall,up_der,wall_dist,i) result (dudz)
  real(kind=dp) :: u_opp,u_wall,wall_dist,dudz
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dudz = (u_opp-u_wall)/(1.0D0+wall_dist)
  else
    dudz = (u_wall-u_opp)/(1.0D0+wall_dist)
  end if
end function

function dvdx_ibb(v_opp,v_wall,up_der,wall_dist,i) result (dvdx)
  real(kind=dp) :: v_opp,v_wall,wall_dist,dvdx
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dvdx = (v_opp-v_wall)/(1.0D0+wall_dist)
  else
    dvdx = (v_wall-v_opp)/(1.0D0+wall_dist)
  end if
end function

function dvdy_ibb(v_opp,v_wall,up_der,wall_dist,i) result (dvdy)
  real(kind=dp) :: v_opp,v_wall,wall_dist,dvdy
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dvdy = (v_opp-v_wall)/(1.0D0+wall_dist)
  else
    dvdy = (v_wall-v_opp)/(1.0D0+wall_dist)
  end if
end function

function dvdz_ibb(v_opp,v_wall,up_der,wall_dist,i) result (dvdz)
  real(kind=dp) :: v_opp,v_wall,wall_dist,dvdz
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dvdz = (v_opp-v_wall)/(1.0D0+wall_dist)
  else
    dvdz = (v_wall-v_opp)/(1.0D0+wall_dist)
  end if
end function

function dwdx_ibb(w_opp,w_wall,up_der,wall_dist,i) result (dwdx)
  real(kind=dp) :: w_opp,w_wall,wall_dist,dwdx
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dwdx = (w_opp-w_wall)/(1.0D0+wall_dist)
  else
    dwdx = (w_wall-w_opp)/(1.0D0+wall_dist)
  end if
end function

function dwdy_ibb(w_opp,w_wall,up_der,wall_dist,i) result (dwdy)
  real(kind=dp) :: w_opp,w_wall,wall_dist,dwdy
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dwdy = (w_opp-w_wall)/(1.0D0+wall_dist)
  else
    dwdy = (w_wall-w_opp)/(1.0D0+wall_dist)
  end if
end function

function dwdz_ibb(w_opp,w_wall,up_der,wall_dist,i) result (dwdz)
  real(kind=dp) :: w_opp,w_wall,wall_dist,dwdz
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dwdz = (w_opp-w_wall)/(1.0D0+wall_dist)
  else
    dwdz = (w_wall-w_opp)/(1.0D0+wall_dist)
  end if
end function

function dTdx_ibb(T_opp,T_wall,up_der,wall_dist,i) result (dTdx)
  real(kind=dp) :: T_opp,T_wall,wall_dist,dTdx
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dTdx = (T_opp-T_wall)/(1.0D0+wall_dist)
  else
    dTdx = (T_wall-T_opp)/(1.0D0+wall_dist)
  end if
end function

function dTdy_ibb(T_opp,T_wall,up_der,wall_dist,i) result (dTdy)
  real(kind=dp) :: T_opp,T_wall,wall_dist,dTdy
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dTdy = (T_opp-T_wall)/(1.0D0+wall_dist)
  else
    dTdy = (T_wall-T_opp)/(1.0D0+wall_dist)
  end if
end function

function dTdz_ibb(T_opp,T_wall,up_der,wall_dist,i) result (dTdz)
  real(kind=dp) :: T_opp,T_wall,wall_dist,dTdz
  integer :: i
  logical :: up_der

!  wall_dist = ABS(real((link_number-1001)/1000000000))

  if (up_der) then
    dTdz = (T_opp-T_wall)/(1.0D0+wall_dist)
  else
    dTdz = (T_wall-T_opp)/(1.0D0+wall_dist)
  end if
end function
!
function dudx_local(u_rest,u_up,u_down,state_rest,state_up,state_down) result (dudx)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: u_rest,u_up,u_down,dudx

  if (state_up >=0 .and. state_down >= 0) then
    dudx = (u_down - u_up)/2.0D0
  else if (state_up >= 0) then
    dudx = u_rest-u_up
  else if (state_down >= 0) then
    dudx = u_down - u_rest
  else
    dudx = 0.0D0
  end if
end function
!
function dudy_local(u_rest,u_up,u_down,state_rest,state_up,state_down) result (dudy)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: u_rest,u_up,u_down,dudy

  if (state_up >=0 .and. state_down >= 0) then
    dudy = (u_down - u_up)/2.0D0
  else if (state_up >= 0) then
    dudy = u_rest-u_up
  else if (state_down >= 0) then
    dudy = u_down - u_rest
  else
    dudy = 0.0D0
  end if
end function
!
function dudz_local(u_rest,u_up,u_down,state_rest,state_up,state_down) result (dudz)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: u_rest,u_up,u_down,dudz

  if (state_up >=0 .and. state_down >= 0) then
    dudz = (u_down - u_up)/2.0D0
  else if (state_up >= 0) then
    dudz = u_rest-u_up
  else if (state_down >= 0) then
    dudz = u_down - u_rest
  else
    dudz = 0.0D0
  end if
end function
!
function dvdx_local(v_rest,v_up,v_down,state_rest,state_up,state_down) result (dvdx)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: v_rest,v_up,v_down,dvdx

  if (state_up >=0 .and. state_down >= 0) then
    dvdx = (v_down - v_up)/2.0D0
  else if (state_up >= 0) then
    dvdx = v_rest-v_up
  else if (state_down >= 0) then
    dvdx = v_down - v_rest
  else
    dvdx = 0.0D0
  end if
end function
!
function dvdy_local(v_rest,v_up,v_down,state_rest,state_up,state_down) result (dvdy)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: v_rest,v_up,v_down,dvdy

  if (state_up >=0 .and. state_down >= 0) then
    dvdy = (v_down - v_up)/2.0D0
  else if (state_up >= 0) then
    dvdy = v_rest-v_up
  else if (state_down >= 0) then
    dvdy = v_down - v_rest
  else
    dvdy = 0.0D0
  end if
end function
!
function dvdz_local(v_rest,v_up,v_down,state_rest,state_up,state_down) result (dvdz)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: v_rest,v_up,v_down,dvdz

  if (state_up >=0 .and. state_down >= 0) then
    dvdz = (v_down - v_up)/2.0D0
  else if (state_up >= 0) then
    dvdz = v_rest-v_up
  else if (state_down >= 0) then
    dvdz = v_down - v_rest
  else
    dvdz = 0.0D0
  end if
end function
!
function dwdx_local(w_rest,w_up,w_down,state_rest,state_up,state_down) result (dwdx)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: w_rest,w_up,w_down,dwdx

  if (state_up >=0 .and. state_down >= 0) then
    dwdx = (w_down - w_up)/2.0D0
  else if (state_up >= 0) then
    dwdx = w_rest - w_up
  else if (state_down >= 0) then
    dwdx = w_down - w_rest
  else
    dwdx = 0.0D0
  end if
end function
!
function dwdy_local(w_rest,w_up,w_down,state_rest,state_up,state_down) result (dwdy)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: w_rest,w_up,w_down,dwdy

  if (state_up >=0 .and. state_down >= 0) then
    dwdy = (w_down - w_up)/2.0D0
  else if (state_up >= 0) then
    dwdy = w_rest - w_up
  else if (state_down >= 0) then
    dwdy = w_down - w_rest
  else
    dwdy = 0.0D0
  end if
end function
!
function dwdz_local(w_rest,w_up,w_down,state_rest,state_up,state_down) result (dwdz)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: w_rest,w_up,w_down,dwdz

  if (state_up >=0 .and. state_down >= 0) then
    dwdz = (w_down - w_up)/2.0D0
  else if (state_up >= 0) then
    dwdz = w_rest - w_up
  else if (state_down >= 0) then
    dwdz = w_down - w_rest
  else
    dwdz = 0.0D0
  end if
end function
!
function dTdx_local(temp_rest,temp_up,temp_down,state_rest,state_up,state_down) result (dTdx)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: temp_rest,temp_up,temp_down,dTdx

  if (state_up >=0 .and. state_down >= 0) then
    dTdx = (temp_down - temp_up)/2.0D0
  else if (state_up >= 0) then
    dTdx = temp_rest - temp_up
  else if (state_down >= 0) then
    dTdx = temp_down - temp_rest
  else
    dTdx = 0.0D0
  end if
end function
!
function dTdy_local(temp_rest,temp_up,temp_down,state_rest,state_up,state_down) result (dTdy)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: temp_rest,temp_up,temp_down,dTdy

  if (state_up >=0 .and. state_down >= 0) then
    dTdy = (temp_down - temp_up)/2.0D0
  else if (state_up >= 0) then
    dTdy = temp_rest - temp_up
  else if (state_down >= 0) then
    dTdy = temp_down - temp_rest
  else
    dTdy = 0.0D0
  end if
end function
!
function dTdz_local(temp_rest,temp_up,temp_down,state_rest,state_up,state_down) result (dTdz)
  integer :: state_rest,state_up,state_down
  real(kind=dp) :: temp_rest,temp_up,temp_down,dTdz

  if (state_up >=0 .and. state_down >= 0) then
    dTdz = (temp_down - temp_up)/2.0D0
  else if (state_up >= 0) then
    dTdz = temp_rest - temp_up
  else if (state_down >= 0) then
    dTdz = temp_down - temp_rest
  else
    dTdz = 0.0D0
  end if
end function
!
!
! Functions for GGM derivatives, twitchier than those above.
!
!
function derivative_windage_firstord(loc_q,loc_state,wind_dir,upwind,downwind,central) result (dqda)
  integer,intent(in) :: loc_state(3)
  real(kind=dp),intent(in) :: loc_q(3),wind_dir
  real(kind=dp) :: dqda
  logical,intent(in) :: upwind,downwind,central
!
  if (central) then
    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
      dqda = (loc_q(4)-loc_q(2))/2.0D0
    else
      if (loc_state(2) >= 0) then
        dqda = loc_q(2) - loc_q(1)
      else if (loc_state(4) >= 0) then
        dqda = loc_q(1) - loc_q(4)
      else
        dqda = 0.0D0
      end if
    end if
!
  else if (upwind) then
    if (wind_dir > 0.0D0) then
      if (loc_state(3) >= 0) then
        dqda = loc_q(3)-loc_q(1)
      else
        dqda = 0.0D0
      end if
    else
      if (loc_state(2) >= 0) then
        dqda = loc_q(1)-loc_q(2)
      else
        dqda = 0.0D0
      end if
    end if
!
  else if (downwind) then
    if (wind_dir < 0.0D0) then
      if (loc_state(3) >= 0) then
        dqda = loc_q(3)-loc_q(1)
      else
        dqda = 0.0D0
      end if
    else
      if (loc_state(2) >= 0) then
        dqda = loc_q(1)-loc_q(2)
      else
        dqda = 0.0D0
      end if
    end if
!
  else
    dqda = 0.0D0
  end if
!
end function
!
function derivative_windage_secondord(loc_q,loc_state,wind_dir,upwind,downwind,central) result (dqda)
  integer,intent(in) :: loc_state(5)
  real(kind=dp),intent(in) :: loc_q(5),wind_dir
  real(kind=dp) :: dqda
  logical,intent(in) :: upwind,downwind,central
!
  if (central) then
    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
      dqda = (loc_q(4)-loc_q(2))/2.0D0
    else
      if (loc_state(2) >= 0) then
        dqda = loc_q(2) - loc_q(1)
      else if (loc_state(4) >= 0) then
        dqda = loc_q(1) - loc_q(4)
      else
        dqda = 0.0D0
      end if
    end if
!
  else if (upwind) then
    if (wind_dir > 0.0D0) then
      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
        dqda = (3*loc_q(1)-4*loc_q(4)+loc_q(5))/2
      else
        dqda = derivative_windage_firstord((/loc_q(1),loc_q(2),loc_q(4)/),&
          (/loc_state(1),loc_state(2),loc_state(4)/),wind_dir,upwind,downwind,central)
      end if
    else
      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
        dqda = (-3*loc_q(1)+4*loc_q(4)-loc_q(5))/2
      else
        dqda = derivative_windage_firstord((/loc_q(1),loc_q(2),loc_q(4)/),&
          (/loc_state(1),loc_state(2),loc_state(4)/),wind_dir,upwind,downwind,central)
      end if
    end if
!
  else if (downwind) then
    if (wind_dir < 0.0D0) then
      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
        dqda = (3*loc_q(1)-4*loc_q(4)+loc_q(5))/2
      else
        dqda = derivative_windage_firstord((/loc_q(1),loc_q(2),loc_q(4)/),&
          (/loc_state(1),loc_state(2),loc_state(4)/),wind_dir,upwind,downwind,central)
      end if
    else
      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
        dqda = (-3*loc_q(1)+4*loc_q(4)-loc_q(5))/2
      else
        dqda = derivative_windage_firstord((/loc_q(1),loc_q(2),loc_q(4)/),&
          (/loc_state(1),loc_state(2),loc_state(4)/),wind_dir,upwind,downwind,central)
      end if
    end if
!
  else
    dqda = 0.0D0
  end if
end function
!
function derivative_windage_highord(loc_q,loc_state,wind_dir,upwind,downwind,central) result (dqda)
  integer,intent(in) :: loc_state(7)
  real(kind=dp),intent(in) :: loc_q(7),wind_dir
  real(kind=dp) :: dqda
  logical,intent(in) :: upwind,downwind,central

!  if (central) then
!    if
!  else if (upwind) then
!
!  else if(downwind) then
!
!  else
!    dqda = 0.0D0
!  end if

end function

function high_order_known(loc_q,central,upwind,downwind)  result (dqda)
  real(kind=dp),intent(in) :: loc_q(7)
  real(kind=dp) :: dqda
  logical,intent(in) :: central,upwind,downwind

  if (central) then
    dqda = (loc_q(6) - 8*loc_q(5) + 8*loc_q(2) - loc_q(3))/12.0D0
  else if (upwind) then
    dqda = (loc_q(6) - 4.0D0*loc_q(5) + 3.0D0*loc_q(1))/2.0D0
  else if (downwind) then
    dqda = (-3.0D0*loc_q(1) + 4.0D0*loc_q(2) - loc_q(3))/2.0D0
  else
    dqda = 0.0D0
  end if
end function

function second_order_known(loc_q,central,upwind,downwind)  result (dqda)
  real(kind=dp),intent(in) :: loc_q(5)
  real(kind=dp) :: dqda
  logical,intent(in) :: central,upwind,downwind

  if (central) then
    dqda = (loc_q(5) - 8*loc_q(4) + 8*loc_q(2) - loc_q(3))/12.0D0
  else if (upwind) then
    dqda = (loc_q(5) - 4.0D0*loc_q(4) + 3.0D0*loc_q(1))/2.0D0
  else if (downwind) then
    dqda = (-3.0D0*loc_q(1) + 4.0D0*loc_q(2) - loc_q(3))/2.0D0
  else
    dqda = 0.0D0
  end if
end function


function first_order_vel_cond_unknown(q_rest,q_target,q_up,q_down,i,&
  rest_state,up_state,down_state) result (dqda)
!
!
!
real(kind=dp) :: q_rest,q_up,q_down,q_target,dir_state,dqda
integer :: rest_state,up_state,down_state,i
logical :: winds

!
if (up_state >= 0 .and. down_state >= 0) then
dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
!
else if (up_state >= 0 .and. down_state < -1000) then
dir_state = (abs(real(down_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
winds = .false.
dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
!
else if (up_state < -1000 .and. down_state >= 0) then
dir_state = (abs(real(up_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
winds = .true.
dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
!
else if (up_state < -1000 .and. down_state < -1000) then
  dqda = 0.0D0
else
  dqda = 0.0D0
end if

end function

function first_order_temp_cond_unknown(q_rest,q_target,q_up,q_down,i,&
  rest_state,up_state,down_state) result (dqda)
use freestream_values
!
!
!
real(kind=dp) :: q_rest,q_target,q_up,q_down,dir_state,dqda
integer :: rest_state,up_state,down_state,i
logical :: winds

!
if (up_state >= 0 .and. down_state >= 0) then
dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
!
else if (up_state >= 0 .and. down_state < -1000) then
dir_state = (abs(real(down_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
winds = .false.
dqda = dudx_ibb(q_target,temperature,winds,dir_state,i)
!
else if (up_state < -1000 .and. down_state >= 0) then
dir_state = (abs(real(up_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
winds = .true.
dqda = dudx_ibb(q_target,temperature,winds,dir_state,i)
!
else if (up_state < -1000 .and. down_state < -1000) then
  dqda = 0.0D0
else
  dqda = 0.0D0
end if

end function

function first_order_shifted_vel_cond_unknown(q_rest,q_target,q_up,q_down,i,&
  card,rest_state,up_state,down_state) result (dqda)
!
! Currently only for flows shifted in the -x direction, i.e. M>1 in the x-dir
!
real(kind=dp) :: q_rest,q_up,q_down,q_target,dir_state,dqda
integer :: rest_state,up_state,down_state,solid_state,i,card
logical :: winds
!
dir_state = (abs(real(down_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
!
select case(i)
  case(1)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(2)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(3)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(4)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(5)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(6)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(7)
    if (card == 3) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(8)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(9)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(10)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(11)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(12)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(13)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(14)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(15)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    end if
  case(16)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(17)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(18)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(19)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(20)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(21)
    if (card == 3) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(22)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(23)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(24)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(25)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(26)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(27)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(28)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(29)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(30)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(31)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(32)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(33)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(34)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(35)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(36)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(37)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(38)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(39)
    if (card == 3) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if

end select
!if (up_state >= 0 .and. down_state >= 0) then
!dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
!!
!else if (up_state < -1000 .and. down_state == -1001) then
!
!else if (down_start < -1000 .and. up_state == -1001) then
!
!
!
!else if (up_state >= 0 .and. down_state < -1000) then
!dir_state = (abs(real(down_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
!winds = .false.
!dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
!!
!else if (up_state < -1000 .and. down_state >= 0) then
!dir_state = (abs(real(up_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
!winds = .true.
!dqda = dudx_ibb(q_target,0.0D0,winds,dir_state,i)
!!
!else if (up_state < -1000 .and. down_state < -1000) then
!  dqda = 0.0D0
!else
!  dqda = 0.0D0
!end if

end function

function first_order_shifted_temp_cond_unknown(q_rest,q_target,q_up,q_down,i,&
  card,rest_state,up_state,down_state) result (dqda)
use freestream_values
! Currently only for flows shifted in the -x direction, i.e. M>1 in the x-dir
!
real(kind=dp) :: q_rest,q_up,q_down,q_target,dir_state,dqda
integer :: rest_state,up_state,down_state,i,card
logical :: winds
!
dir_state = (abs(real(down_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
!
select case(i)
  case(1)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(2)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(3)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(4)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(5)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(6)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(7)
    if (card == 3) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(8)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(9)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(10)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(11)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(12)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(13)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(14)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(15)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    end if
  case(16)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(17)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(18)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(19)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(20)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(21)
    if (card == 3) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(22)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(23)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(24)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(25)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(26)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(27)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(28)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(29)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(30)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(31)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(32)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(33)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 3) then !dqdy
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(34)
    if (card == 1) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(35)
    if (card == 2) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(36)
    if (card == 3) then !dqdx
      winds = .false.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(37)
    if (card == 1) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(38)
    if (card == 2) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 1) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 3) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if
  case(39)
    if (card == 3) then !dqdx
      winds = .true.
      dqda = dudx_ibb(q_target,q_target,winds,dir_state,i)
    else if (card == 2) then !dqdy
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    else if (card == 1) then !dqdz
      dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
    end if

end select

!real(kind=dp) :: q_rest,q_target,q_up,q_down,dir_state,dqda
!integer :: rest_state,up_state,down_state,i
!logical :: winds
!
!!
!if (up_state >= 0 .and. down_state >= 0) then
!dqda = dudx_local(q_rest,q_up,q_down,rest_state,up_state,down_state)
!!
!else if (up_state >= 0 .and. down_state < -1000) then
!dir_state = (abs(real(down_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
!winds = .false.
!dqda = dudx_ibb(q_target,temperature,winds,dir_state,i)
!!
!else if (up_state < -1000 .and. down_state >= 0) then
!dir_state = (abs(real(up_state,kind=dp) + 1001.0D0)/1.0D9) * abs_latt_leng(i)
!winds = .true.
!dqda = dudx_ibb(q_target,temperature,winds,dir_state,i)
!!
!else if (up_state < -1000 .and. down_state < -1000) then
!  dqda = 0.0D0
!else
!  dqda = 0.0D0
!end if

end function

!

!
!
!
!
!
!
!
!
!
!
!
!
!
end module


!function dudx_unknown_secondord(loc_u,loc_state,upwind,downwind,central) result (dudx)
!  integer :: loc_state(5)
!  real(kind=dp) :: loc_u(5),dudx
!  logical :: upwind,downwind,central
!!
!  if (central) then
!    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
!      dudx = (loc_u(2)+loc_u(4))/2
!    else
!      dudx = dudx_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!    end if
!!
!  else if (upwind) then
!    if (loc_u(1) > 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dudx = (3*loc_u(1)-4*loc_u(4)+loc_u(5))/2
!      else
!        dudx = dudx_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dudx = (-3*loc_u(1)+4*loc_u(4)-loc_u(5))/2
!      else
!        dudx = dudx_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else if (downwind) then
!    if (loc_u(1) < 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dudx = (3*loc_u(1)-4*loc_u(4)+loc_u(5))/2
!      else
!        dudx = dudx_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dudx = (-3*loc_u(1)+4*loc_u(4)-loc_u(5))/2
!      else
!        dudx = dudx_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else
!    dudx = 0.0D0
!  end if
!end function
!!
!function dudy_unknown_secondord(loc_u,loc_state,wind_dir,upwind,downwind,central) result (dudy)
!  integer :: loc_state(5)
!  real(kind=dp) :: loc_u(5),dudy
!  logical :: upwind,downwind,central
!!
!  if (central) then
!    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
!      dudx = (loc_u(2)+loc_u(4))/2
!    else
!      dudx = dudx_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!    end if
!!
!  else if (upwind) then
!    if (wind_dir > 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dudy = (3*loc_u(1)-4*loc_u(4)+loc_u(5))/2
!      else
!        dudy = dudy_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dudy = (-3*loc_u(1)+4*loc_u(4)-loc_u(5))/2
!      else
!        dudy = dudy_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else if (downwind) then
!    if (wind_dir < 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dudy = (3*loc_u(1)-4*loc_u(4)+loc_u(5))/2
!      else
!        dudy = dudy_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dudy = (-3*loc_u(1)+4*loc_u(4)-loc_u(5))/2
!      else
!        dudy = dudy_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else
!    dudy = 0.0D0
!  end if
!end function
!!
!function dudz_unknown_secondord(loc_u,loc_state,wind_dir,upwind,downwind,central) result (dudz)
!  integer :: loc_state(5)
!  real(kind=dp) :: loc_u(5),dudz
!  logical :: upwind,downwind,central
!!
!  if (central) then
!    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
!      dudz = (loc_u(2)+loc_u(4))/2
!    else
!      dudz = dudz_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!    end if
!!
!  else if (upwind) then
!    if (wind_dir > 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dudz = (3*loc_u(1)-4*loc_u(4)+loc_u(5))/2
!      else
!        dudz = dudz_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dudz = (-3*loc_u(1)+4*loc_u(4)-loc_u(5))/2
!      else
!        dudz = dudz_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else if (downwind) then
!    if (wind_dir < 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dudz = (3*loc_u(1)-4*loc_u(4)+loc_u(5))/2
!      else
!        dudz = dudz_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dudz = (-3*loc_u(1)+4*loc_u(4)-loc_u(5))/2
!      else
!        dudz = dudz_local(loc_u(1),loc_u(2),loc_u(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else
!    dudz = 0.0D0
!  end if
!end function
!!
!function dvdx_unknown_secondord(loc_v,loc_state,wind_dir,upwind,downwind,central) result (dvdx)
!  integer :: loc_state(5)
!  real(kind=dp) :: loc_v(5),dvdx
!  logical :: upwind,downwind,central
!!
!  if (central) then
!    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
!      dvdx = (loc_v(2)+loc_v(4))/2
!    else
!      dvdx = dvdx_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!    end if
!!
!  else if (upwind) then
!    if (wind_dir > 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dvdx = (3*loc_v(1)-4*loc_v(4)+loc_v(5))/2
!      else
!        dvdx = dvdx_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dvdx = (-3*loc_v(1)+4*loc_v(2)-loc_v(3))/2
!      else
!        dvdx = dvdx_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else if (downwind) then
!    if (wind_dir < 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dvdx = (3*loc_v(1)-4*loc_v(4)+loc_v(5))/2
!      else
!        dvdx = dvdx_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dvdx = (-3*loc_v(1)+4*loc_v(2)-loc_v(3))/2
!      else
!        dvdx = dvdx_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else
!    dvdx = 0.0D0
!  end if
!end function
!!
!function dvdy_unknown_secondord(loc_u,loc_state,upwind,downwind,central) result (dvdy)
!  integer :: loc_state(5)
!  real(kind=dp) :: loc_v(5),dvdy
!  logical :: upwind,downwind,central
!!
!  if (central) then
!    if (loc_state(2) >= 0 .and. loc_state(4) >= 0) then
!      dudx = (loc_v(2)+loc_v(4))/2
!    else
!      dvdy = dvdy_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!    end if
!!
!  else if (upwind) then
!    if (loc_v(1) > 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dvdy = (3*loc_v(1)-4*loc_v(4)+loc_v(5))/2
!      else
!        dvdy = dvdy_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dvdy = (-3*loc_v(1)+4*loc_v(2)-loc_v(3))/2
!      else
!        dvdy = dvdy_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else if (downwind) then
!    if (loc_v(1) < 0.0D0) then
!      if (loc_state(4) >= 0 .and. loc_state(5) >= 0) then
!        dvdy = (3*loc_v(1)-4*loc_v(4)+loc_v(5))/2
!      else
!        dvdy = dvdy_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    else
!      if (loc_state(2) >= 0 .and. loc_state(3) >= 0) then
!        dvdy = (-3*loc_v(1)+4*loc_v(2)-loc_v(3))/2
!      else
!        dvdy = dvdy_local(loc_v(1),loc_v(2),loc_v(4),loc_state(1),loc_state(2),loc_state(4))
!      end if
!    end if
!!
!  else
!    dvdy = 0.0D0
!  end if
!end function
