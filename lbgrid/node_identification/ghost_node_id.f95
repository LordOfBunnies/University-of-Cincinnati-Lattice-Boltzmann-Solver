subroutine ghost_node_id(fine,num)
!
! Identifies ghost nodes and
!
! The ghost node id's work as follows:
!   If the node overlaps a box on the same level, it keeps its level designation
!
!   If a node does not overlap, it is sent to a derivativer checker to find which
!     derivatives are viable for the node.
!   The check returns positives in a binary form
!
!
!
! The final form of the ghost node id is:
!   level * 1k
!   upwind-derivatives (x = 1, y = 2, z = 3), so if all are valid, it's 7*100
!   downwind-derivatives (x = 1, y = 2, z = 3), so if all are valid, it's 7*10
!   central difference-derivatives (x = 1, y = 2, z = 3), so if all are valid, it's 7*1
!
! A final node number might look like 4755
!
!
! Called by: grid_gen
! Calls: ghost_derivative_checker
! External Calls:
!
use precise
use constants
use amr_info_holder, only: grid_info,state,mfstate,nghosts,nghosts_state,nghosts_mid
use amr_processes
use amrex_amr_module
use amrex_base_module
implicit none
integer :: lvl,fine,a,b,c,g,up,down,center,second_up,second_down,second_center,num
logical :: inside,overlap,also,state_far_mx,state_far_px,state_far_my,&
  state_far_py,state_far_mz,state_far_pz

type(amrex_mfiter) :: mfi
type(amrex_box) :: socks

!write(*,*) 'identifying ghost nodes'

do lvl = 0,fine
  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

  do while (mfi%next())
    state => mfstate(lvl)%dataptr(mfi)
    socks = mfi%validbox()
!    write(*,*) 'socks values',socks%lo,socks%hi
    do c = socks%lo(3)-nghosts_state,socks%hi(3)+nghosts_state
      do b = socks%lo(2)-nghosts_state,socks%hi(2)+nghosts_state
        do a = socks%lo(1)-nghosts_state,socks%hi(1)+nghosts_state
!
! check if node overlaps another box on the same level
!
          if (state(a,b,c,1) < 0) cycle
          if (state(a,b,c,1) >= 200 .and. state(a,b,c,1) <= 299) cycle
          if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
              c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
              b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then

            if (state(a,b,c,1) >= -299 .and. state(a,b,c,1) <-200) then
              cycle
            else
              state(a,b,c,1) = -666
              cycle
            end if
          end if
!          if (state(a,b,c,1) >= 51000 .and. state(a,b,c,1) < 100000 ) cycle

          inside = magic_box(socks%lo,socks%hi,(/a,b,c/))
          if (.not. inside) then
            if (state(a,b,c,1) >= 0) then
              also = .false.
              do g = 1,grid_info(lvl,num)%big_zones
! So it doesn't flag itself
                if (grid_info(lvl,num)%lower(1,g) == socks%lo(1) .and. &
                    grid_info(lvl,num)%lower(2,g) == socks%lo(2) .and. &
                    grid_info(lvl,num)%lower(3,g) == socks%lo(3) .and. &
                    grid_info(lvl,num)%higher(1,g) == socks%hi(1) .and. &
                    grid_info(lvl,num)%higher(2,g) == socks%hi(2) .and. &
                    grid_info(lvl,num)%higher(3,g) == socks%hi(3)) then
                  cycle
                end if

                also = magic_box(grid_info(lvl,num)%lower(1:3,g),&
                  grid_info(lvl,num)%higher(1:3,g),(/a,b,c/))
!
! If yes, declare it simply
!
                if (also) then
                  state(a,b,c,1) = (lvl+1)*1000
!                  if (a == 66 .and. b == 37 .and. c == 32) then
!                   write(*,*) 'ARGH!!!',state(a,b,c,1:39),self!,ducks%lo,ducks%hi,self
!                  end if
! gtfo
                  exit
                end if
              end do
! shortcut everything else
              if (also) then
                cycle
              else
!
! Check for the ability to do a second order derivative
!
                if (mod(a,2) == 0 .and. mod(b,2) == 0 .and. mod(c,2)==0) then
                  overlap = .true.

                  if (a-4 >= socks%lo(1)-nghosts_mid) then
                    if (state(a-4,b,c,1) > 0) then
                      state_far_mx = .true.
                    end if
                  else
                    state_far_mx = .false.
                  end if
                  if (a+4 <= socks%hi(1)+nghosts_mid) then
                    if (state(a+4,b,c,1) > 0) then
                      state_far_px = .true.
                    end if
                  else
                    state_far_px = .false.
                  end if
!
                  if (b-4 >= socks%lo(2)-nghosts_mid) then
                    if (state(a,b-4,c,1) > 0) then
                      state_far_my = .true.
                    end if
                  else
                    state_far_my = .false.
                  end if
                  if (b+4 <= socks%hi(2)+nghosts_mid) then
                    if (state(a,b+4,c,1) > 0) then
                      state_far_py = .true.
                    end if
                  else
                    state_far_py = .false.
                  end if
!
                  if (c-4 >= socks%lo(3)-nghosts_mid) then
                    if (state(a,b,c-4,1) > 0) then
                      state_far_mz = .true.
                    end if
                  else
                    state_far_mz = .false.
                  end if
                  if (c+4 <= socks%hi(3)+nghosts_mid) then
                    if (state(a,b,c+4,1) > 0) then
                      state_far_pz = .true.
                    end if
                  else
                    state_far_pz = .false.
                  end if
!
!
!
                else
                  overlap = .false.

                  if (a-3 >= socks%lo(1)-nghosts_mid) then
                    if (state(a-3,b,c,1) > 0) then
                      state_far_mx = .true.
                    end if
                  else
                    state_far_mx = .false.
                  end if
                  if (a+3 <= socks%hi(1)+nghosts_mid) then
                    if (state(a+3,b,c,1) > 0) then
                      state_far_px = .true.
                    end if
                  else
                    state_far_px = .false.
                  end if
!
                  if (b-3 >= socks%lo(2)-nghosts_mid) then
                    if (state(a,b-3,c,1) > 0) then
                      state_far_my = .true.
                    end if
                  else
                    state_far_my = .false.
                  end if
                  if (b+3 <= socks%hi(2)+nghosts_mid) then
                    if (state(a,b+3,c,1) > 0) then
                      state_far_py = .true.
                    end if
                  else
                    state_far_py = .false.
                  end if
!
                  if (c-3 >= socks%lo(3)-nghosts_mid) then
                    if (state(a,b,c-3,1) > 0) then
                      state_far_mz = .true.
                    end if
                  else
                    state_far_mz = .false.
                  end if
                  if (c+3 <= socks%hi(3)+nghosts_mid) then
                    if (state(a,b,c+3,1) > 0) then
                      state_far_pz = .true.
                    end if
                  else
                    state_far_pz = .false.
                  end if
                end if
!
!
!
                if (lvl == 0) then
                  cycle
                else
                  call ghost_derivative_checker(state(a,b,c,:),state_far_mx,state_far_px,&
                    state_far_my,state_far_py,state_far_mz,state_far_pz,socks%lo,socks%hi,&
                    a,b,c,second_up,second_down,second_center,up,down,center,overlap)
                  !write(*,*) 'snarf',up,down,center,state(a,b,c,1)
                  state(a,b,c,1) = (lvl+1)*10000000 + second_up*100000 + &
                    second_down*10000 + second_center*1000 + up*100 + down*10 + center
                end if
!                if (a == 66 .and. b == 37 .and. c == 32) then
!                 write(*,*) 'ARGH!!!',state(a,b,c,1:39),self!,ducks%lo,ducks%hi,self
!               end if
!                else
!                  state(a,b,c,1) = (lvl+51)*1000
!                end if
              end if
            end if
          end if
!                if (a == 66 .and. b == 37 .and. c == 32) then
!                 write(*,*) 'ARGH!!!',state(a,b,c,1:39),self!,ducks%lo,ducks%hi,self
!               end if
        end do
      end do
    end do


  end do

  call amrex_mfiter_destroy(mfi)

end do

end subroutine
                !write(*,*) 'barf'
!                if (mod(a,2) == 0 .and. mod(b,2) == 0 .and. mod(c,2) == 0) then
!          if (a == 51 .and. b == 37 .and. c == 2 .and. lvl == 1) then
!            write(*,*) 'ARGH!!!',state(a,b,c,1),state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!                  state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),socks%lo,socks%hi,self
!          end if
!            if (a == 113 .and. b == 49 .and. c == 49 .and. lvl == 2) then
!              write(*,*) 'node status, ghost',state(a,b,c,1),a,b,c,lvl
!            end if
