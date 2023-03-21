subroutine predictive_splitting(lvl,settag,cleartag,tux_lo,tux_hi,&
  tag_lo,tag_hi,tagarr,g,mfi)
!
! This is to attempt to split the grid ahead of a moving disturbance so the it moves
!   refined grid and less information is lost
!
!
!
! Called by: calculate_gradgrad_mag
! Calls:
! External calls:
!
USE precise
USE constants
use iso_c_binding
use derivatives
USE ggm_stuff
use amr_info_holder
use amrex_amr_module
use amrex_base_module
implicit none
real(kind=dp) :: vel_dir,pgs_limit,pgs_upper
integer,intent(in) :: tux_lo(3),tux_hi(3),g,lvl,tag_lo(4),tag_hi(4)
integer :: a,b,c,x,y,z,v,reach(6)
character(kind=c_char),intent(in) :: settag,cleartag
character(kind=c_char),intent(inout) :: tagarr(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),&
  tag_lo(3):tag_hi(3))
logical :: x_upwind,y_upwind,z_upwind,u_dominant,v_dominant,w_dominant,victory

type(amrex_mfiter) :: mfi
!
! reach(1) = +x, 2 = +y, 3 = +z, 4 = -x
!
!
!
pgs_limit = 2*ggm_sensitivity
pgs_upper = 3*ggm_sensitivity

u_vel => mfu_vel(lvl)%dataptr(mfi)
v_vel => mfv_vel(lvl)%dataptr(mfi)
w_vel => mfw_vel(lvl)%dataptr(mfi)

call mfggm(lvl)% fill_boundary(amrex_geom(lvl),1,g)

do c = tux_lo(3),tux_hi(3)
  do b = tux_lo(2),tux_hi(2)
    do a = tux_lo(1),tux_hi(1)

      if (state(a,b,c,1) < 0) then
        cycle
      end if

      u_dominant = .false.
      v_dominant = .false.
      w_dominant = .false.
      victory = .false.
!
! Determine if you need to look upwind or downwind
!   upwind = true means the point is downwind
!
      if (u_vel(a,b,c,1) > 0.0D0) then
        x_upwind = .true.
      else
        x_upwind = .false.
      end if

      if (v_vel(a,b,c,1) > 0.0D0) then
        y_upwind = .true.
      else
        y_upwind = .false.
      end if

      if (w_vel(a,b,c,1) > 0.0D0) then
        z_upwind = .true.
      else
        z_upwind = .false.
      end if
!
! Determine if one velocity is far stronger than the others
!
      if (abs(u_vel(a,b,c,1)) > 5*abs(v_vel(a,b,c,1)) .and. &
          abs(u_vel(a,b,c,1)) > 5*abs(w_vel(a,b,c,1)) ) then

        u_dominant = .true.

!
! Determine which way we're going
!
        if (x_upwind) then
! Determine how far to scan
          if (abs(tux_lo(1)-nghosts - a) > pgs_max_scan) then
            reach(4) = a - pgs_max_scan
          else
            reach(4) = tux_lo(1)-nghosts
          end if
! Go from the edge of the bounding box with ghosts to the current node
          do x = reach(4),a
! Only scan valid fluid nodes and overlapped ghosts, and only above the limit
            if (state(x,b,c,1) >= 0 .and. state(x,b,c,1) < 100000 .and. &
                  ggm(x,b,c,g) >= pgs_limit) then
! Declare victory then walk away, no need to scan the rest because of the dominant
!   velocity here
              victory = .true.
              exit
            end if
          end do
! Tag all the nodes between the two points
          if (victory) then
            do v = x,a
              if (state(v,b,c,1) >= 0) tagarr(v,b,c) = settag
            end do
          end if
!
        else
!
          if (abs(tux_hi(1)+nghosts - a) > pgs_max_scan) then
            reach(1) = a + pgs_max_scan
          else
            reach(1) = tux_hi(1)+nghosts
          end if
! Go from the edge of the bounding box with ghosts to the current node
          do x = a,reach(1)
! Only scan valid fluid nodes and overlapped ghosts, and only above the limit
            if (state(x,b,c,1) >= 0 .and. state(x,b,c,1) < 100000 .and. &
                  ggm(x,b,c,g) >= pgs_limit) then
! Declare victory then walk away, no need to scan the rest because of the dominant
!   velocity here
              victory = .true.
              exit
            end if
          end do
! Tag all the nodes between the two points
          if (victory) then
            do v = a,x
              if (state(v,b,c,1) >= 0) tagarr(v,b,c) = settag
            end do
          end if
        end if
!
      else if (abs(v_vel(a,b,c,1)) > 5*abs(u_vel(a,b,c,1)) .and. &
          abs(v_vel(a,b,c,1)) > 5*abs(w_vel(a,b,c,1)) ) then

        v_dominant = .true.

!
! Determine which way we're going
!
        if (y_upwind) then
! Determine how far to scan
          if (abs(tux_lo(2)-nghosts - b) > pgs_max_scan) then
            reach(5) = b - pgs_max_scan
          else
            reach(5) = tux_lo(2)-nghosts
          end if
! Go from the edge of the bounding box with ghosts to the current node
          do y = reach(5),b
! Only scan valid fluid nodes and overlapped ghosts, and only above the limit
            if (state(a,y,c,1) >= 0 .and. state(a,y,c,1) < 100000 .and. &
                  ggm(a,y,c,g) >= pgs_limit) then
! Declare victory then walk away, no need to scan the rest because of the dominant
!   velocity here
              victory = .true.
              exit
            end if
          end do
! Tag all the nodes between the two points
          if (victory) then
            do v = y,b
              if (state(a,v,c,1) >= 0) tagarr(a,v,c) = settag
            end do
          end if
!
        else
!
          if (abs(tux_hi(2)+nghosts - b) > pgs_max_scan) then
            reach(2) = b + pgs_max_scan
          else
            reach(2) = tux_hi(2)+nghosts
          end if
! Go from the edge of the bounding box with ghosts to the current node
          do y = b,reach(2)
! Only scan valid fluid nodes and overlapped ghosts, and only above the limit
            if (state(a,y,c,1) >= 0 .and. state(a,y,c,1) < 100000 .and. &
                  ggm(a,y,c,g) >= pgs_limit) then
! Declare victory then walk away, no need to scan the rest because of the dominant
!   velocity here
              victory = .true.
              exit
            end if
          end do
! Tag all the nodes between the two points
          if (victory) then
            do v = b,y
              if (state(a,v,c,1) >= 0) tagarr(a,v,c) = settag
            end do
          end if
        end if
!
      else if (abs(w_vel(a,b,c,1)) > 5*abs(v_vel(a,b,c,1)) .and. &
          abs(w_vel(a,b,c,1)) > 5*abs(u_vel(a,b,c,1)) ) then
        w_dominant = .true.
!
! Determine which way we're going
!
        if (z_upwind) then
! Determine how far to scan
          if (abs(tux_lo(3)-nghosts - c) > pgs_max_scan) then
            reach(6) = c - pgs_max_scan
          else
            reach(6) = tux_lo(3)-nghosts
          end if
! Go from the edge of the bounding box with ghosts to the current node
          do z = reach(6),c
! Only scan valid fluid nodes and overlapped ghosts, and only above the limit
            if (state(a,b,z,1) >= 0 .and. state(a,b,z,1) < 100000 .and. &
                  ggm(a,b,z,g) >= pgs_limit) then
! Declare victory then walk away, no need to scan the rest because of the dominant
!   velocity here
              victory = .true.
              exit
            end if
          end do
! Tag all the nodes between the two points
          if (victory) then
            do v = z,c
              if (state(a,b,v,1) >= 0) tagarr(a,b,v) = settag
            end do
          end if
!
        else
!
          if (abs(tux_hi(3)+nghosts - c) > pgs_max_scan) then
            reach(3) = c + pgs_max_scan
          else
            reach(3) = tux_hi(3)+nghosts
          end if
! Go from the edge of the bounding box with ghosts to the current node
          do z = c,reach(3)
! Only scan valid fluid nodes and overlapped ghosts, and only above the limit
            if (state(a,b,z,1) >= 0 .and. state(a,b,z,1) < 100000 .and. &
                  ggm(a,b,z,g) >= pgs_limit) then
! Declare victory then walk away, no need to scan the rest because of the dominant
!   velocity here
              victory = .true.
              exit
            end if
          end do
! Tag all the nodes between the two points
          if (victory) then
            do v = c,z
              if (state(a,b,v,1) >= 0) tagarr(a,b,v) = settag
            end do
          end if
        end if
!
!
!
      else
        if (abs(tux_hi(1)+nghosts - a) > pgs_max_scan) then
          reach(1) = a + pgs_max_scan
        else
          reach(1) = tux_hi(1)+nghosts
        end if

        if (abs(tux_hi(2)+nghosts - b) > pgs_max_scan) then
          reach(2) = b + pgs_max_scan
        else
          reach(2) = tux_hi(2)+nghosts
        end if

        if (abs(tux_hi(3)+nghosts - c) > pgs_max_scan) then
          reach(3) = c + pgs_max_scan
        else
          reach(3) = tux_hi(3)+nghosts
        end if

        if (abs(tux_lo(1)-nghosts - a) > pgs_max_scan) then
          reach(4) = a - pgs_max_scan
        else
          reach(4) = tux_lo(1)-nghosts
        end if

        if (abs(tux_lo(2)-nghosts - b) > pgs_max_scan) then
          reach(5) = b - pgs_max_scan
        else
          reach(5) = tux_lo(2)-nghosts
        end if

        if (abs(tux_lo(3)-nghosts - c) > pgs_max_scan) then
          reach(6) = c - pgs_max_scan
        else
          reach(6) = tux_lo(3)-nghosts
        end if

        do x = reach(4),reach(1)
          if (state(x,b,c,1) >= 0 .and. state(x,b,c,1) < 100000 .and. &
                  ggm(x,b,c,g) >= pgs_upper) then
            tagarr(x,b,c) = settag
          end if
        end do

        do y = reach(5),reach(2)
          if (state(a,y,c,1) >= 0 .and. state(a,y,c,1) < 100000 .and. &
                  ggm(a,y,c,g) >= pgs_upper) then
            tagarr(a,y,c) = settag
          end if
        end do

        do z = reach(6),reach(3)
          if (state(a,b,z,1) >= 0 .and. state(a,b,z,1) < 100000 .and. &
                  ggm(a,b,z,g) >= pgs_upper) then
            tagarr(a,b,z) = settag
          end if
        end do

      end if
!
!  Scan for disturbances
!
    end do
  end do
end do


end subroutine
