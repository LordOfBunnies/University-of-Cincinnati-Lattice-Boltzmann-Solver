subroutine nullify_lvl(lvl,num)
!
! Same as the routine nullify_nodes, but only for one level and the update
!   of grid info is done elsewhere
!
!
! Called by: from_coarse,remake_lvl (amr_processes
! Calls: ml_overlapping_boxes
! External calls: amrex_mfiter_build/destroy
!
use precise
use constants
use grid_data, only: dir
use mpi
use linkwise
use amrex_amr_module
use amrex_base_module
use amr_info_holder, only: state,mfstate,grid_info,nghosts,nghosts_state,nghosts_mid
use amr_processes

implicit none
integer,intent(in) :: lvl,num
integer :: over_lo(3),over_hi(3),a,b,c,i,gr,loc_zones,lvl_zones
logical :: overlap,inside_lo,inside_hi

type(amrex_box) :: rox
type(amrex_mfiter) :: mfi

lvl_zones = total_big_zones(lvl+1)
!loc_zones = count_big_zones(lvl+1)
!
call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

!  convert_counter = 0
!write(*,*) 'lvl zones',lvl_zones,loc_zones,self
do while (mfi%next())
  state => mfstate(lvl)%dataptr(mfi)
  rox = mfi%validbox()

  do gr = 1,lvl_zones
    !See if the coarser box contains the finer box
    overlap = .false.
!    write(*,*) 'check check',gr,self
    call ml_overlapping_boxes(rox%lo-nghosts_state,&
      grid_info(lvl+1,num)%lower(1:3,gr),rox%hi+nghosts_state,&
      grid_info(lvl+1,num)%higher(1:3,gr),over_lo,over_hi,overlap,lvl,lvl+1)
!    write(*,*) 'all clear',gr,self
!    write(*,*) 'overlapping stuff',over_lo,over_hi,lvl,self

    if (overlap) then

!    write(*,*) 'box with overlap',rox%lo,rox%hi,over_lo,over_hi,lvl,self

!        write(*,*) 'coarse box corners',rox%lo-nghosts,rox%hi+nghosts
!        write(*,*) 'fine box corners',grid_info(lvl+1)%lower(1:3,gr)/2,grid_info(lvl+1)%higher(1:3,gr)/2
!        write(*,*) 'overlap region',over_lo,over_hi,lvl

      do c = over_lo(3),over_hi(3)
        do b = over_lo(2),over_hi(2)
          do a = over_lo(1),over_hi(1)
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
!            if (over_lo(1) < rox%lo(1)-nghosts .or. over_lo(2) < rox%lo(2)-nghosts .or. &
!                over_lo(3) < rox%lo(3)-nghosts .or. over_hi(1) > rox%hi(1)+nghosts .or. &
!                over_hi(2) > rox%hi(2)+nghosts .or. over_hi(3) > rox%hi(3)+nghosts) then
!              write(*,*) 'jiggle',rox%lo,rox%hi,over_lo,over_hi,a,b,c,lvl,self
!            end if
            if (state(a,b,c,1) >= 0) then
              state(a,b,c,1) = -666

!              if (a == 108 .and. b == 72 .and. c == 35) then
!                write(*,*) 'overlap checks',rox%lo,rox%hi,grid_info(lvl+1,num)%lower(1:3,gr),&
!                grid_info(lvl+1,num)%higher(1:3,gr),over_lo,over_hi,lvl
!              end if
              !write(*,*) 'jiggle',a,b,c,lvl
!                convert_counter = convert_counter+1
            end if

          end do
        end do
      end do

    else
      cycle
    end if !overlapping regions conditional

  end do
end do
! write(*,*) 'Number of nodes converted to null',convert_counter,lvl

call amrex_mfiter_destroy(mfi)
end subroutine




    !inside_lo = magic_box(rox_lo,rox_hi,grid_info(lvl+1)%lower(1,gr))
    !inside_hi = magic_box(rox_lo,rox_hi,grid_info(lvl+1)%lhigher(1,gr))
