subroutine flow_direction(a,b,c,dudx,dudy,dudz,dvdx,dvdy,dvdz,&
  dwdx,dwdy,dwdz,dTdx,dTdy,dTdz,order_of_acc)
!
! Subroutine to determine the direction of the flow and apply the correct derivative
!
!
! Called by:
! Calls:
! External calls:
!
use amrex_base_module
use amrex_amr_module
use amrex_info_holder
use precise
use grid_data, only: cx,cy,cz,opp
implicit none
real(kind=dp) :: loc_u,loc_v,loc_w,loc_temp
real(kind=dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
real(kind=dp) :: dTdx,dTdy,dTdz
integer :: i,a,b,c,order_of_acc








end subroutine
