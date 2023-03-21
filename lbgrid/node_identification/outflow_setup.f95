subroutine outlfow_setup(fine)
!
! For outflow boundary conditions, valid nodes are tacked onto the end
!   of the flow region.  During streaming, these are updated with the old
!   fout values WITHOUT the streaming action.  Additionally, directions
!   pointing back into the valid fluid region are streamed as normal.
!
!
! Called by: grid_gen
! Calls:
! External calls:
!
use precise
use constants
use nml_inlet_outlet
use amr_info_holder, only: mfstate,state,nghosts_mid
use amr_processes, only: magic_box
use amrex_amr_module
use amrex_base_module
implicit none
integer,intent(in) :: fine
integer :: a,b,c,i,r,lvl
integer :: outflow_ghosts,out_num
logical :: inside

type(amrex_mfiter) :: mfi
type(amrex_box) :: fox


!outflow_ghosts = nghosts_mid/2
!
!! copy basic intro
!do lvl = 0,fine
!  !write(*,*) 'fluffy bunnies... ATTACK!!!'
!  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
!!
!  do while(mfi%next())
!    fox = mfi%validbox()
!    state => mfstate(lvl)%dataptr(mfi)
!!
!!
!!
!    if (outlet_side(1)) then
!      do r = 1,num_outlets
!        if (outlet_type(2,r) == 1) then
!          out_num = r
!          exit
!        end if
!      end do
!
!      do c =fox%lo(3),fox%hi(3)
!        do b = fox%lo(2),fox%hi(2)
!          do a = fox%lo(1)-outflow_ghosts,fox%lo(1)-1
!
!            inside = magic_box(amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi+1,(/a,b,c/))
!
!            if (state(a,b,c,1) == -666 .and. .not. inside) then
!
!              state(a,b,c,1) = 200+out_num
!
!            end if
!
!          end do
!        end do
!      end do
!
!    end if
!!
!!
!!
!    if (outlet_side(2) ) then
!      do r = 1,num_outlets
!        if (outlet_type(2,r) == 2) then
!          out_num = r
!          exit
!        end if
!      end do
!
!      do c =fox%lo(3),fox%hi(3)
!        do b = fox%lo(2)-outflow_ghosts,fox%lo(2)-1
!          do a = fox%lo(1),fox%hi(1)
!
!            inside = magic_box(amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi+1,(/a,b,c/))
!
!            if (state(a,b,c,1) == -666 .and. .not. inside) then
!
!              state(a,b,c,1) = 200+out_num
!
!            end if
!
!          end do
!        end do
!      end do
!    end if
!!
!!
!!
!    if (outlet_side(3) ) then
!      do r = 1,num_outlets
!        if (outlet_type(2,r) == 3) then
!          out_num = r
!          exit
!        end if
!      end do
!
!      do c =fox%lo(3)-outflow_ghosts,fox%lo(3)-1
!        do b = fox%lo(2),fox%hi(2)
!          do a = fox%lo(1),fox%hi(1)
!
!            inside = magic_box(amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi+1,(/a,b,c/))
!
!            if (state(a,b,c,1) == -666 .and. .not. inside) then
!
!              state(a,b,c,1) = 200+out_num
!
!            end if
!
!          end do
!        end do
!      end do
!    end if
!!
!!
!!
!    if (outlet_side(4) ) then
!      !write(*,*) 'Destroy the furry weasels!'
!      do r = 1,num_outlets
!        if (outlet_type(2,r) == 4) then
!          !write(*,*) 'Burn the teddy bears!'
!          out_num = r
!          exit
!        end if
!      end do
!
!      do c =fox%lo(3),fox%hi(3)
!        do b = fox%lo(2),fox%hi(2)
!          do a = fox%hi(1)+1,fox%hi(1)+outflow_ghosts
!            !write(*,*) 'Down with the tyrannous squirrels!',a,b,c
!!            if (a == 201 .and. b == 20 .and. c == 20) then
!!              write(*,*) 'Great googly moogly',state(a,b,c,1),a,b,c
!!            end if
!            inside = magic_box(amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi+1,(/a,b,c/))
!
!            if (state(a,b,c,1) == -666 .and. .not. inside) then
!
!              state(a,b,c,1) = 200+out_num
!
!            end if
!!            if (a == 201 .and. b == 20 .and. c == 20) then
!!              write(*,*) 'Great googly foogly',state(a,b,c,1),a,b,c
!!            end if
!
!          end do
!        end do
!      end do
!    end if
!!
!!
!!
!    if (outlet_side(5) ) then
!      do r = 1,num_outlets
!        if (outlet_type(2,r) == 5) then
!          out_num = r
!          exit
!        end if
!      end do
!
!      do c =fox%lo(3),fox%hi(3)
!        do b = fox%hi(2)+1,fox%hi(2)+outflow_ghosts
!          do a = fox%lo(1),fox%hi(1)
!
!            inside = magic_box(amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi+1,(/a,b,c/))
!
!            if (state(a,b,c,1) == -666 .and. .not. inside) then
!
!              state(a,b,c,1) = 200+out_num
!
!            end if
!
!          end do
!        end do
!      end do
!    end if
!!
!!
!!
!    if (outlet_side(6) ) then
!      do r = 1,num_outlets
!        if (outlet_type(2,r) == 6) then
!          out_num = r
!          exit
!        end if
!      end do
!
!      do c =fox%hi(3)+1,fox%hi(3)+outflow_ghosts
!        do b = fox%lo(2),fox%hi(2)
!          do a = fox%lo(1),fox%hi(1)
!
!            inside = magic_box(amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi+1,(/a,b,c/))
!
!            if (state(a,b,c,1) == -666 .and. .not. inside) then
!
!              state(a,b,c,1) = 200+out_num
!
!            end if
!
!          end do
!        end do
!      end do
!    end if
!
!
!  end do
!
!  call amrex_mfiter_destroy(mfi)
!
!end do

end subroutine
