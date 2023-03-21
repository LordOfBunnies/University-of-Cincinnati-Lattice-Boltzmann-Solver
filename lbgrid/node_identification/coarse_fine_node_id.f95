subroutine coarse_fine_node_id(fine)
!
! Fine which coarse nodes are overlapped by ghost nodes.
!
!
! Called by: grid_gen
! Calls:
! External calls:
!
use precise
use constants
use grid_data
use linkwise
use amrex_amr_module
use amrex_base_module
use amr_info_holder, only: state,mfstate,nghosts,grid_info,nghosts_state
use amr_processes, only: magic_box
implicit none
integer :: a,b,c,i,gr,lvl,fine,overlap_counter,coarse_fine_counter
integer :: up,down,center,num
logical :: valid
logical,allocatable :: inside(:),not_ghost(:)

type(amrex_mfiter) :: mfi
type(amrex_box) :: rox

num = 1

write(*,*) 'coarse/fine node id'

!
! The basic idea is to find nodes on a coarser grid which overlap the ghost nodes
!   on a finer grid.  The nodes must be exclusively fine ghost nodes, and not overlap
!   other fine valid nodes, i.e. overlapped ghost nodes.
!
! The nodes are
!
!

!
! Go through every level which include ghost node overlaps
!
do lvl = 0,fine-1

  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

  allocate(inside(grid_info(lvl+1,num)%big_zones))
  allocate(not_ghost(grid_info(lvl+1,num)%big_zones))
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
          do i = 1,grid_info(lvl+1,num)%big_zones

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
              state(a,b,c,1) = (lvl+1)*100000 + up*100 + down*10 + center
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

end do



end subroutine


!            if (a == 113 .and. b == 49 .and. c == 49 .and. lvl == 2) then
!              write(*,*) 'node status, cfid',state(a,b,c,1),a,b,c,lvl
!            end if
!          inside = .false.
!          inside = magic_box(grid_info(lvl+1)%lower(1:3,i)/2,grid_info(lvl+1)%higher(1:3,i)/2,(/a,b,c/))
!          if (.not. inside) then

!
!          inside = .false.
!
!          do i = 1,grid_info(lvl,num)%big_zones
!
!            inside = magic_box(grid_info(lvl+1)%lower(1:3,i)/2,grid_info(lvl+1)%higher(1:3,i)/2,(/a,b,c/))
!            if (inside) then
!              state(a,b,c,1) = -666
!              exit
!            end if
!
!          end do

!          if (inside .and. .not. not_ghost) then
!
!            call coarse_derivative_checker(state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!              state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),up,down,center)
!
!            state(a,b,c,1) = (lvl+1)*1E7 + up*100 + down*10 + center
!
!          else if (inside .and. not_ghost) then
!            do j = 2,7
!              if (state(a+cx(j),b+cy(j),c+cz(j),1) >=0) then
!                call coarse_derivative_checker(state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!                  state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),up,down,center)
!
!                state(a,b,c,1) = (lvl+1)*1E7 + up*100 + down*10 + center
!              else
!                state(a,b,c,1) = -666
!              end if
!            end do
!          !else
!
!          end if
!              coarse_fine_counter = coarse_fine_counter +1
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
!
!
!
