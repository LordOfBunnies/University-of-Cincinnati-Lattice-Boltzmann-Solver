subroutine coarse_derivative_checker(loc_state_px,loc_state_py,loc_state_pz,&
  loc_state_nx,loc_state_ny,loc_state_nz,loc_node,box_lo,box_hi,up,down,center)
!
! Determine if
!
!
! Called by: coarse_fine_node_id
! Calls:
! External calls:
!
use precise
use constants
use grid_data
use linkwise
implicit none
integer,intent(in) :: loc_state_px,loc_state_py,loc_state_pz,loc_state_nx,&
  loc_state_ny,loc_state_nz
integer,intent(in) :: box_lo(3),box_hi(3),loc_node(3)
integer :: up,down,center,i
logical :: upwind(3),downwind(3),central(3)

central = .true.
up = 0
down = 0
center = 0
!
if (dimensions == 2) then

else
!
! 3D cases
!
! Null nodes still contain valid fluid variable information from averaging down
! in the -x direction
  if (loc_state_nx == -666 .or. loc_state_nx >= 0) then
! If it exceeds the bounds of the box, it dies.
    if (loc_node(1)+cx(5) < box_lo(1)) then
      upwind(1) = .false.
      central(1) = .false.
    else
      upwind(1) = .true.
    end if

  else
    upwind(1) = .false.
    central(1) = .false.
  end if
! in the -y direction
  if (loc_state_ny == -666 .or. loc_state_ny >= 0) then
    if (loc_node(2)+cy(6) < box_lo(2)) then
      upwind(2) = .false.
      central(2) = .false.
    else
      upwind(2) = .true.
    end if
  else
    upwind(2) = .false.
    central(2) = .false.
  end if
! in the -z direction
  if (loc_state_nz == -666 .or. loc_state_nz >= 0) then
    if (loc_node(3)+cz(7) < box_lo(3)) then
      upwind(3) = .false.
      central(3) = .false.
    else
      upwind(3) = .true.
    end if
  else
    upwind(3) = .false.
    central(3) = .false.
  end if
! in the +x direction
  if (loc_state_px == -666 .or. loc_state_px >= 0) then
    if (loc_node(1)+cx(2) > box_hi(1)) then
      downwind(1) = .false.
      central(1) = .false.
    else
      downwind(1) = .true.
    end if
  else
    central(1) = .false.
    downwind(1) = .false.
  end if
! in the +y direction
  if (loc_state_py == -666 .or. loc_state_py >= 0) then

    if (loc_node(2)+cy(3) > box_hi(2)) then
      downwind(2) = .false.
      central(2) = .false.
    else
      downwind(2) = .true.
    end if
  else
    downwind(2) = .false.
    central(2) = .false.
  end if
! in the +z direction
  if (loc_state_pz == -666 .or. loc_state_pz >= 0) then

    if (loc_node(3)+cz(4) > box_hi(3)) then
      downwind(3) = .false.
      central(3) = .false.
    else
      downwind(3) = .true.
    end if
  else
    downwind(3) = .false.
    central(3) = .false.
  end if

  do i = 1,3
    if (upwind(i)) up = up + 2**(i-1)
    if (downwind(i)) down = down + 2**(i-1)
    if (central(i)) center = center + 2**(i-1)
  end do

end if

end subroutine
