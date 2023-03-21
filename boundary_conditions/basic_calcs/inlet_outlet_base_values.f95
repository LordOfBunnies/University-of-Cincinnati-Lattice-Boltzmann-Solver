subroutine inlet_outlet_base_values(io_num,is_it_inlet,loc_rho,&
  loc_u,loc_v,loc_w,loc_temp,loc_press)
!
! Calculate the rest of the required variables based on the information for the
! inlet or outlet type
!
!
! Called by:
! Calls:
!
use precise
use constants
use nml_inlet_outlet
use freestream_values
implicit none
integer,intent(IN) :: io_num
real(kind=dp),intent(out) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_press
real(kind=dp) :: vel_mag
logical,intent(in) :: is_it_inlet

write(*,*) 'giggle snort',io_num,is_it_inlet,inlet_type(1,io_num),outlet_type(1,io_num)

if (inlet_type(1,io_num) == 1 .and. is_it_inlet) then
!  if (is_it_inlet) then
    loc_rho = inlet_flow_val(1,io_num)
    loc_temp = inlet_flow_val(3,io_num)
    if (inlet_type(2,io_num) == 1) then
      loc_u = inlet_flow_val(2,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 2) then
      loc_u = 0.0D0
      loc_v = inlet_flow_val(2,io_num)
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 3) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = inlet_flow_val(2,io_num)
    else if (inlet_type(2,io_num) == 4) then
      loc_u = inlet_flow_val(2,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 5) then
      loc_u = 0.0D0
      loc_v = inlet_flow_val(2,io_num)
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 6) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = inlet_flow_val(2,io_num)
    end if
    loc_press = loc_rho*loc_temp+loc_rho*inlet_flow_val(2,io_num)**2
  else if (outlet_type(1,io_num) == 1 .and. .not. is_it_inlet) then
    loc_rho = outlet_flow_val(1,io_num)
    loc_temp = outlet_flow_val(3,io_num)
    if (outlet_type(2,io_num) == 1) then
      loc_u = outlet_flow_val(2,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 2) then
      loc_u = 0.0D0
      loc_v = outlet_flow_val(2,io_num)
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 3) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = outlet_flow_val(2,io_num)
    else if (outlet_type(2,io_num) == 4) then
      loc_u = outlet_flow_val(2,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 5) then
      loc_u = 0.0D0
      loc_v = outlet_flow_val(2,io_num)
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 6) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = outlet_flow_val(2,io_num)
    end if
    loc_press = loc_rho*loc_temp+loc_rho*outlet_flow_val(2,io_num)**2
 ! end if

else if (inlet_type(1,io_num) == 2 .and. is_it_inlet) then
!  if (is_it_inlet) then
    loc_press = inlet_flow_val(1,io_num)
    loc_temp = inlet_flow_val(2,io_num)
    loc_rho = inlet_flow_val(1,io_num)/(non_dim_R*inlet_flow_val(2,io_num))
    loc_u = 0.0D0
    loc_v = 0.0D0
    loc_w = 0.0D0
!    bernoulli = SQRT((2.0D0/loc_rho) * (loc_press + 0.5D0*loc_rho*&
!      vel_mag**2 - loc_press))
else if (outlet_type(1,io_num) == 2 .and. .not. is_it_inlet) then
    loc_press = outlet_flow_val(1,io_num)
    loc_temp = outlet_flow_val(2,io_num)
    loc_rho = outlet_flow_val(1,io_num)/(non_dim_R*outlet_flow_val(2,io_num))
    loc_u = 0.0D0
    loc_v = 0.0D0
    loc_w = 0.0D0
!  end if

else if (inlet_type(1,io_num) == 3 .and. is_it_inlet) then
! Freestream data
!  write(*,*) 'gruff'
  return
else if (outlet_type(1,io_num) ==3 .and. .not. is_it_inlet) then

  return
else if (inlet_type(1,io_num) == 4 .and. is_it_inlet) then
  loc_rho = rho_inf
!  if (is_it_inlet) then
    loc_temp = inlet_flow_val(2,io_num)
    if (inlet_type(2,io_num) == 1) then
      loc_u = inlet_flow_val(1,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 2) then
      loc_u = 0.0D0
      loc_v = inlet_flow_val(1,io_num)
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 3) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = inlet_flow_val(1,io_num)
    else if (inlet_type(2,io_num) == 4) then
      loc_u = inlet_flow_val(1,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 5) then
      loc_u = 0.0D0
      loc_v = inlet_flow_val(1,io_num)
      loc_w = 0.0D0
    else if (inlet_type(2,io_num) == 6) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = inlet_flow_val(1,io_num)
    end if
    loc_press = loc_rho*loc_temp+loc_rho*inlet_flow_val(2,io_num)**2

else if (outlet_type(1,io_num) == 4 .and. .not. is_it_inlet) then
    loc_temp = outlet_flow_val(2,io_num)
    if (outlet_type(2,io_num) == 1) then
      loc_u = outlet_flow_val(1,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 2) then
      loc_u = 0.0D0
      loc_v = outlet_flow_val(1,io_num)
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 3) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = outlet_flow_val(1,io_num)
    else if (outlet_type(2,io_num) == 4) then
      loc_u = outlet_flow_val(1,io_num)
      loc_v = 0.0D0
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 5) then
      loc_u = 0.0D0
      loc_v = outlet_flow_val(1,io_num)
      loc_w = 0.0D0
    else if (outlet_type(2,io_num) == 6) then
      loc_u = 0.0D0
      loc_v = 0.0D0
      loc_w = outlet_flow_val(1,io_num)
    end if
    loc_press = loc_rho*loc_temp+loc_rho*outlet_flow_val(2,io_num)**2
!  end if

else if (inlet_type(1,io_num) == 7 .and. is_it_inlet) then
!  write(*,*) 'eh?'
!  if (is_it_inlet) then
    loc_rho = inlet_flow_val(1,io_num)
    loc_u = inlet_flow_val(2,io_num)
    loc_v = inlet_flow_val(3,io_num)
    loc_w = inlet_flow_val(4,io_num)
    loc_temp = inlet_flow_val(5,io_num)
else if (outlet_type(1,io_num) == 7 .and. .not. is_it_inlet) then
!  write(*,*) 'outlet pre-assign',outlet_flow_val(1,io_num),outlet_flow_val(2,io_num),&
!      outlet_flow_val(3,io_num),outlet_flow_val(4,io_num),outlet_flow_val(5,io_num)
    loc_rho = outlet_flow_val(1,io_num)
    loc_u = outlet_flow_val(2,io_num)
    loc_v = outlet_flow_val(3,io_num)
    loc_w = outlet_flow_val(4,io_num)
    loc_temp = outlet_flow_val(5,io_num)
!  end if

else
  write(*,*) 'raffle'
  return
end if

end subroutine
