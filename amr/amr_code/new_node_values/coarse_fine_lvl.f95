subroutine coarse_fine_lvl(lvl,num)
!
!
!
!
!
!
!
!
! Called by:
! Calls:
!
use precise
use constants
use grid_data
use linkwise
use amrex_amr_module
use amrex_base_module
use amr_info_holder, only: state,mfstate,nghosts,grid_info,nghosts_state
use amr_processes, only: magic_box,total_big_zones
implicit none
integer,intent(in) :: lvl,num
integer :: a,b,c,i,gr,fine,overlap_counter,coarse_fine_counter
integer :: up,down,center,loc_zones,lvl_zones
logical :: valid
logical,allocatable :: inside(:),not_ghost(:)

type(amrex_mfiter) :: mfi
type(amrex_box) :: rox

lvl_zones = total_big_zones(lvl+1)

  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

  allocate(inside(lvl_zones))
  allocate(not_ghost(lvl_zones))
  coarse_fine_counter = 0

  do while (mfi%next())
    state => mfstate(lvl)%dataptr(mfi)
    rox = mfi%validbox()

    do c = rox%lo(3)-nghosts_state,rox%hi(3)+nghosts_state
      do b = rox%lo(2)-nghosts_state,rox%hi(2)+nghosts_state
        do a = rox%lo(1)-nghosts_state,rox%hi(1)+nghosts_state
!
! Skip unneeded nodes
!
          if (state(a,b,c,1) < 0) cycle
!
! check node against all boxes on level

          overlap_counter = 0
          inside = .false.
          not_ghost = .false.
          valid = .false.
!
!
!
          do i = 1,lvl_zones

!
! Inside check if the node is inside the valid box on a finer level
!
            inside(i) = magic_box((grid_info(lvl+1,num)%lower(1:3,i)-nghosts)/2,&
                (grid_info(lvl+1,num)%higher(1:3,i)+nghosts)/2,(/a,b,c/))
!
! not_ghost checks if the node is inside a different box on the same level
!
            not_ghost(i) = magic_box(grid_info(lvl+1,num)%lower(1:3,i)/2,&
              grid_info(lvl+1,num)%higher(1:3,i)/2,(/a,b,c/))
          end do
!
!
!
          if (any(not_ghost)) cycle

          if (any(inside) .and. .not. any(not_ghost)) then


!
              call coarse_derivative_checker(state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
                state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),(/a,b,c/),&
                rox%lo,rox%hi,up,down,center)
!
              if (any(inside) .and. .not. any(not_ghost)) then
                state(a,b,c,1) = (lvl+1)*100000 + up*100 + down*10 + center
              end if
!              coarse_fine_counter = coarse_fine_counter +1
          end if


        end do
      end do
    end do

  end do
!  write(*,*) 'overlapped nodes changed',coarse_fine_counter
  call amrex_mfiter_destroy(mfi)
!
! Deallocate before the next time around, the size won't change except based on level
!
  deallocate(inside)
  deallocate(not_ghost)

end subroutine

!            if (a == 113 .and. b == 49 .and. c == 49 .and. lvl == 2) then
!              write(*,*) 'node status, cfid',state(a,b,c,1),a,b,c,lvl
!            end if
!
!              write(*,*) 'overlap counter ',overlap_counter,state(a,b,c,1),state(a,b,c,2),&
!                state(a,b,c,3),state(a,b,c,4),&
!                state(a,b,c,5),state(a,b,c,6),state(a,b,c,7)
!
            !end if
!            do gr = 1,grid_info(lvl+1)%big_zones
!              if (inside(gr) .and. not_ghost(gr)) then
!                overlap_counter = overlap_counter + 1
!                valid = .true.
!              else if (inside(gr) .and. .not. not_ghost(gr)) then
!                state(a,b,c,1) = -666
!                valid = .false.
!                exit
!              end if
!            end do
            !if (valid) then
