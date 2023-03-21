subroutine ghost_derivative_checker(loc_state,state_far_mx,state_far_px,&
  state_far_my,state_far_py,state_far_mz,state_far_pz,&
  lo_box,hi_box,a,b,c,second_up,second_down,second_center,&
  up,down,center,overlap)
!
! Checks which derivatives are viable for a given ghost node
!
! The values are, in essence, binary, being 1 if true and 0 if false.
! In order, the digits are x, y, and z, thus, a value of 7 (4 + 2+ 1)
!   is that the derivative is good in all three directions
!
! The value saved depends on the its placement:
!   100s = upwind
!   10s = downwind
!   1s = central difference
!
! Called by: ghost_node_id
! Calls:
!
use precise
use constants
use grid_data, only: dir
use amr_info_holder, only : nghosts,nghosts_mid,state
use linkwise
implicit none
integer,intent(in) :: loc_state(a:a,b:b,c:c,1:dir),lo_box(3),hi_box(3),a,b,c
integer :: second_up,second_down,second_center
integer :: i,up,down,center
logical :: upwind(3,2),downwind(3,2),central(3,2),overlap
logical :: state_far_mx,state_far_px,&
  state_far_my,state_far_py,state_far_mz,state_far_pz


central = .true.
second_up = 0
second_down = 0
second_center = 0
up = 0
down = 0
center = 0

if (dimensions == 2) then

else
!
! Overlapped nodes, these will (hopefully) take a second order upwind/downwind derivative or
!   a fourth order central difference.
!
  if (.not. shifted) then
  if (overlap) then
    if (state_far_mx .and. loc_state(a,b,c,19) >= 0) then
      upwind(1,2) = .true.
      upwind(1,1) = .true.
    else if (loc_state(a,b,c,19) >= 0) then
      upwind(1,1) = .true.
      upwind(1,2) = .false.
      central(1,2) = .false.
    else
      upwind(1,1) = .false.
      upwind(1,2) = .false.
      central(1,1) = .false.
      central(1,2) = .false.
    end if

    if (state_far_px .and. loc_state(a,b,c,16) >= 0) then
      downwind(1,2) = .true.
      downwind(1,1) = .true.
    else if (loc_state(a,b,c,16) >= 0) then
      downwind(1,1) = .true.
      downwind(1,2) = .false.
      central(1,2) = .false.
    else
      downwind(1,1) = .false.
      downwind(1,2) = .false.
      central(1,1) = .false.
      central(1,2) = .false.
    end if

    if (state_far_my .and. loc_state(a,b,c,20) >= 0) then
      upwind(2,2) = .true.
      upwind(2,1) = .true.
    else if (loc_state(a,b,c,20) >= 0) then
      upwind(2,1) = .true.
      upwind(2,2) = .false.
      central(2,2) = .false.
    else
      upwind(2,1) = .false.
      upwind(2,2) = .false.
      central(2,1) = .false.
      central(2,2) = .false.
    end if

    if (state_far_py .and. loc_state(a,b,c,17) >= 0) then
      downwind(2,2) = .true.
      downwind(2,1) = .true.
    else if (loc_state(a,b,c,17) >= 0) then
      downwind(2,1) = .true.
      downwind(2,2) = .false.
      central(2,2) = .false.
    else
      downwind(2,1) = .false.
      downwind(2,2) = .false.
      central(2,1) = .false.
      central(2,2) = .false.
    end if

    if (state_far_mz .and. loc_state(a,b,c,21) >= 0) then
      upwind(3,2) = .true.
      upwind(3,1) = .true.
    else if (loc_state(a,b,c,21) >= 0) then
      upwind(3,1) = .true.
      upwind(3,2) = .false.
      central(3,2) = .false.
    else
      upwind(3,1) = .false.
      upwind(3,2) = .false.
      central(3,1) = .false.
      central(3,2) = .false.
    end if

    if (state_far_pz .and. loc_state(a,b,c,18) >= 0) then
      downwind(3,2) = .true.
      downwind(3,1) = .true.
    else if (loc_state(a,b,c,18) >= 0) then
      downwind(3,1) = .true.
      downwind(3,2) = .false.
      central(3,2) = .false.
    else
      downwind(3,1) = .false.
      downwind(3,2) = .false.
      central(3,1) = .false.
      central(3,2) = .false.
    end if
!
! For non-overlapped nodes, these will be interpolated from overlapped nodes, so need a slightly
!   different stencil (-3,-1,+1,+3)
!
  else
    if (state_far_mx .and. loc_state(a,b,c,5) > 0) then
      upwind(1,2) = .true.
      upwind(1,1) = .true.
    else if (loc_state(a,b,c,5) > 0) then
      upwind(1,1) = .true.
      upwind(1,2) = .false.
      central(1,2) = .false.
    else
      upwind(1,1) = .false.
      upwind(1,2) = .false.
      central(1,1) = .false.
      central(1,2) = .false.
    end if

    if (state_far_px .and. loc_state(a,b,c,2) > 0) then
      downwind(1,2) = .true.
      downwind(1,1) = .true.
    else if (loc_state(a,b,c,2) > 0) then
      downwind(1,1) = .true.
      downwind(1,2) = .false.
      central(1,2) = .false.
    else
      downwind(1,1) = .false.
      downwind(1,2) = .false.
      central(1,1) = .false.
      central(1,2) = .false.
    end if

    if (state_far_my .and. loc_state(a,b,c,6) > 0) then
      upwind(2,2) = .true.
      upwind(2,1) = .true.
    else if (loc_state(a,b,c,6) > 0) then
      upwind(2,1) = .true.
      upwind(2,2) = .false.
      central(2,2) = .false.
    else
      upwind(2,1) = .false.
      upwind(2,2) = .false.
      central(2,1) = .false.
      central(2,2) = .false.
    end if

    if (state_far_py .and. loc_state(a,b,c,3) > 0) then
      downwind(2,2) = .true.
      downwind(2,1) = .true.
    else if (loc_state(a,b,c,3) > 0) then
      downwind(2,1) = .true.
      downwind(2,2) = .false.
      central(2,2) = .false.
    else
      downwind(2,1) = .false.
      downwind(2,2) = .false.
      central(2,1) = .false.
      central(2,2) = .false.
    end if

    if (state_far_mz .and. loc_state(a,b,c,7) > 0) then
      upwind(3,2) = .true.
      upwind(3,1) = .true.
    else if (loc_state(a,b,c,7) > 0) then
      upwind(3,1) = .true.
      upwind(3,2) = .false.
      central(3,2) = .false.
    else
      upwind(3,1) = .false.
      upwind(3,2) = .false.
      central(3,1) = .false.
      central(3,2) = .false.
    end if

    if (state_far_pz .and. loc_state(a,b,c,4) > 0) then
      downwind(3,2) = .true.
      downwind(3,1) = .true.
    else if (loc_state(a,b,c,4) > 0) then
      downwind(3,1) = .true.
      downwind(3,2) = .false.
      central(3,2) = .false.
    else
      downwind(3,1) = .false.
      downwind(3,2) = .false.
      central(3,1) = .false.
      central(3,2) = .false.
    end if
  end if
!
!
!
!
!
!
!
!
!
  else
!
    if (overlap) then
      if (a-2 < lo_box(1)-nghosts_mid) then
          upwind(1,1) = .false.
          upwind(1,2) = .false.
          central(1,1) = .false.
          central(1,2) = .false.
      else
        if (state_far_mx .and. state(a-2,b,c,1) >= 0) then
          upwind(1,2) = .true.
          upwind(1,1) = .true.
        else if (state(a-2,b,c,1) >= 0) then
          upwind(1,1) = .true.
          upwind(1,2) = .false.
          central(1,2) = .false.
        else
          upwind(1,1) = .false.
          upwind(1,2) = .false.
          central(1,1) = .false.
          central(1,2) = .false.
        end if
      end if
!
      if (a+2 > hi_box(1)+nghosts_mid) then
        downwind(1,1) = .false.
        downwind(1,2) = .false.
        central(1,1) = .false.
        central(1,2) = .false.
      else
      if (state_far_px .and. state(a+2,b,c,1) >= 0) then
        downwind(1,2) = .true.
        downwind(1,1) = .true.
      else if (state(a+2,b,c,1) >= 0) then
        downwind(1,1) = .true.
        downwind(1,2) = .false.
        central(1,2) = .false.
      else
        downwind(1,1) = .false.
        downwind(1,2) = .false.
        central(1,1) = .false.
        central(1,2) = .false.
      end if
      end if

      if (b-2 < lo_box(2)-nghosts_mid) then
        upwind(2,1) = .false.
        upwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      else
      if (state_far_my .and. state(a,b-2,c,1) >= 0) then
        upwind(2,2) = .true.
        upwind(2,1) = .true.
      else if (state(a,b-2,c,1) >= 0) then
        upwind(2,1) = .true.
        upwind(2,2) = .false.
        central(2,2) = .false.
      else
        upwind(2,1) = .false.
        upwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      end if
      end if

      if (b+2 > hi_box(2)+nghosts_mid) then
        downwind(2,1) = .false.
        downwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      else
      if (state_far_py .and. state(a,b+2,c,1) >= 0) then
        downwind(2,2) = .true.
        downwind(2,1) = .true.
      else if (state(a,b+2,c,1) >= 0) then
        downwind(2,1) = .true.
        downwind(2,2) = .false.
        central(2,2) = .false.
      else
        downwind(2,1) = .false.
        downwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      end if
      end if

      if (c-2 < lo_box(3)) then
        upwind(3,1) = .false.
        upwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      else
      if (state_far_mz .and. state(a,b,c-2,1) >= 0) then
        upwind(3,2) = .true.
        upwind(3,1) = .true.
      else if (state(a,b,c-2,1) >= 0) then
        upwind(3,1) = .true.
        upwind(3,2) = .false.
        central(3,2) = .false.
      else
        upwind(3,1) = .false.
        upwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      end if
      end if

      if (c+2 > hi_box(3)+nghosts_mid) then
        downwind(3,1) = .false.
        downwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      else
      if (state_far_pz .and. state(a,b,c+2,1) >= 0) then
        downwind(3,2) = .true.
        downwind(3,1) = .true.
      else if (state(a,b,c+2,1) >= 0) then
        downwind(3,1) = .true.
        downwind(3,2) = .false.
        central(3,2) = .false.
      else
        downwind(3,1) = .false.
        downwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      end if
      end if
  !
  ! For non-overlapped nodes, these will be interpolated from overlapped nodes, so need a slightly
  !   different stencil (-3,-1,+1,+3)
  !
    else
      if (a-1 < lo_box(1)) then
        upwind(1,1) = .false.
        upwind(1,2) = .false.
        central(1,1) = .false.
        central(1,2) = .false.
      else
      if (state_far_mx .and. state(a-1,b,c,1) > 0) then
        upwind(1,2) = .true.
        upwind(1,1) = .true.
      else if (state(a-1,b,c,1) > 0) then
        upwind(1,1) = .true.
        upwind(1,2) = .false.
        central(1,2) = .false.
      else
        upwind(1,1) = .false.
        upwind(1,2) = .false.
        central(1,1) = .false.
        central(1,2) = .false.
      end if
      end if

      if (a+1 > hi_box(1)+nghosts_mid) then
        downwind(1,1) = .false.
        downwind(1,2) = .false.
        central(1,1) = .false.
        central(1,2) = .false.
      else
      if (state_far_px .and. state(a+1,b,c,1) > 0) then
        downwind(1,2) = .true.
        downwind(1,1) = .true.
      else if (state(a+1,b,c,1) > 0) then
        downwind(1,1) = .true.
        downwind(1,2) = .false.
        central(1,2) = .false.
      else
        downwind(1,1) = .false.
        downwind(1,2) = .false.
        central(1,1) = .false.
        central(1,2) = .false.
      end if
      end if

      if (b-1 < lo_box(2)-nghosts_mid) then
        upwind(2,1) = .false.
        upwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      else
      if (state_far_my .and. state(a,b-1,c,1) > 0) then
        upwind(2,2) = .true.
        upwind(2,1) = .true.
      else if (state(a,b-1,c,1) > 0) then
        upwind(2,1) = .true.
        upwind(2,2) = .false.
        central(2,2) = .false.
      else
        upwind(2,1) = .false.
        upwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      end if
      end if

      if (b+1 > hi_box(2)+nghosts_mid) then
        downwind(2,1) = .false.
        downwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      else
      if (state_far_py .and. state(a,b+1,c,1) > 0) then
        downwind(2,2) = .true.
        downwind(2,1) = .true.
      else if (state(a,b+1,c,1) > 0) then
        downwind(2,1) = .true.
        downwind(2,2) = .false.
        central(2,2) = .false.
      else
        downwind(2,1) = .false.
        downwind(2,2) = .false.
        central(2,1) = .false.
        central(2,2) = .false.
      end if
      end if

      if (c-1 < lo_box(3)-nghosts_mid) then
        upwind(3,1) = .false.
        upwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      else
      if (state_far_mz .and. state(a,b,c-1,1) > 0) then
        upwind(3,2) = .true.
        upwind(3,1) = .true.
      else if (state(a,b,c-1,1) > 0) then
        upwind(3,1) = .true.
        upwind(3,2) = .false.
        central(3,2) = .false.
      else
        upwind(3,1) = .false.
        upwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      end if
      end if

      if (c+1 > hi_box(3)+nghosts_mid) then
        downwind(3,1) = .false.
        downwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      else
      if (state_far_pz .and. state(a,b,c+1,1) > 0) then
        downwind(3,2) = .true.
        downwind(3,1) = .true.
      else if (state(a,b,c+1,1) > 0) then
        downwind(3,1) = .true.
        downwind(3,2) = .false.
        central(3,2) = .false.
      else
        downwind(3,1) = .false.
        downwind(3,2) = .false.
        central(3,1) = .false.
        central(3,2) = .false.
      end if
      end if


    end if
  end if

!if (abs(loc_node(1)-lo_box(1)) == nghosts .or. loc_state_nx < 0) then
!  upwind(1) = .false.
!  central(1) = .false.
!else
!  upwind(1) = .true.
!end if
!if (abs(loc_node(2)-lo_box(2)) == nghosts .or. loc_state_ny < 0) then
!  upwind(2) = .false.
!  central(2) = .false.
!else
!  upwind(2) = .true.
!end if
!if (abs(loc_node(3)-lo_box(3)) == nghosts .or. loc_state_nz < 0) then
!  upwind(3) = .false.
!  central(3) = .false.
!else
!  upwind(3) = .true.
!end if
!
!if (abs(loc_node(1)-hi_box(1)) == nghosts .or. loc_state_px < 0) then
!  downwind(1) = .false.
!  central(1) = .false.
!else
!  downwind(1) = .true.
!end if
!if (abs(loc_node(2)-hi_box(2)) == nghosts .or. loc_state_py < 0) then
!  downwind(2) = .false.
!  central(2) = .false.
!else
!  downwind(2) = .true.
!end if
!if (abs(loc_node(3)-hi_box(3)) == nghosts .or. loc_state_pz < 0) then
!  downwind(3) = .false.
!  central(3) = .false.
!else
!  downwind(3) = .true.
!end if
!write(*,*) upwind,downwind,central
do i = 1,3
  if (upwind(i,1)) up = up + 2**(i-1)
  if (downwind(i,1)) down = down + 2**(i-1)
  if (central(i,1)) center = center + 2**(i-1)

  if (upwind(i,2)) second_up = second_up + 2**(i-1)
  if (downwind(i,2)) second_down = second_down + 2**(i-1)
  if (central(i,2)) second_center = second_center + 2**(i-1)
end do

end if

end subroutine
