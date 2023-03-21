subroutine fine_to_coarse_node_derivatives(lvl,fine,state_value,&
  loc_state,a,b,c)
!
! This finds nodes which are at the boundary of the valid fine cells, but on
! the coarse level and labels them similarly to overlapped ghost nodes.
!
! Called by: nullify_nodes
! Calls:
!
!
!use precise
!use constants
!implicit none
!integer,intent(in) :: lvl,fine,a,b,c,loc_state(a:a,b:b,c:c,2:7)
!integer, intent(out) :: state_value
!integer :: up,down,center,i
!logical :: upwind(3),downwind(3),central(3)
!
!! in the -x direction
!  if (state(a,b,c,5) == -666 .or. state(a,b,c,5) >= 0 .or. &
!      state(a,b,c,5) == -667) then
!! If it exceeds the bounds of the box, it dies.
!    if (loc_node(1)+cx(5) < box_lo(1)) then
!      upwind(1) = .false.
!      central(1) = .false.
!    else
!      upwind(1) = .true.
!    end if
!
!  else
!    upwind(1) = .false.
!    central(1) = .false.
!  end if
!! in the -y direction
!  if (loc_state_ny == -666 .or. loc_state_ny >= 0) then
!    if (loc_node(2)+cy(6) < box_lo(2)) then
!      upwind(2) = .false.
!      central(2) = .false.
!    else
!      upwind(2) = .true.
!    end if
!  else
!    upwind(2) = .false.
!    central(2) = .false.
!  end if
!! in the -z direction
!  if (loc_state_nz == -666 .or. loc_state_nz >= 0) then
!    if (loc_node(3)+cz(7) < box_lo(3)) then
!      upwind(3) = .false.
!      central(3) = .false.
!    else
!      upwind(3) = .true.
!    end if
!  else
!    upwind(3) = .false.
!    central(3) = .false.
!  end if
!! in the +x direction
!  if (loc_state_px == -666 .or. loc_state_px >= 0) then
!    if (loc_node(1)+cx(2) > box_hi(1)) then
!      downwind(1) = .false.
!      central(1) = .false.
!    else
!      downwind(1) = .true.
!    end if
!  else
!    central(1) = .false.
!    downwind(1) = .false.
!  end if
!! in the +y direction
!  if (loc_state_py == -666 .or. loc_state_py >= 0) then
!
!    if (loc_node(2)+cy(3) > box_hi(2)) then
!      downwind(2) = .false.
!      central(2) = .false.
!    else
!      downwind(2) = .true.
!    end if
!  else
!    downwind(2) = .false.
!    central(2) = .false.
!  end if
!! in the +z direction
!  if (loc_state_pz == -666 .or. loc_state_pz >= 0) then
!
!    if (loc_node(3)+cz(4) > box_hi(3)) then
!      downwind(3) = .false.
!      central(3) = .false.
!    else
!      downwind(3) = .true.
!    end if
!  else
!    downwind(3) = .false.
!    central(3) = .false.
!  end if
!
!  do i = 1,3
!    if (upwind(i)) up = up + 2**(i-1)
!    if (downwind(i)) down = down + 2**(i-1)
!    if (central(i)) center = center + 2**(i-1)
!  end do



end subroutine
