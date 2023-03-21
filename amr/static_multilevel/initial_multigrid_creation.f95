subroutine initial_multigrid_creation(lvl,settag,cleartag,tag,&
  dox_lo,dox_hi,tag_lo,tag_hi,initial)
!
! Create the initial multigrid based on user input. This is done
! without grid adaptation, which comes once the solution is running.
!
!
! Called by: tag_nodes
! Calls:
!
!
use precise
use constants
use grid_data, only: rdx
!use nml_bounding_box
use amrex_base_module
use amrex_amr_module
use iso_c_binding
use amr_info_holder, only: box_lvl,boi_hi,boi_lo,num_boi,nodal_domain_limits
use amr_processes, only: self,get_real_coords
implicit none
integer :: i,a,b,c,x_lo_node_start,y_lo_node_start,&
  z_lo_node_start,x_hi_node_diff,y_hi_node_diff,z_hi_node_diff,lvl
integer :: x_hi_node_end,y_hi_node_end,z_hi_node_end
integer,intent(in) :: tag_lo(4),tag_hi(4),dox_lo(3),dox_hi(3)
integer :: x_max,x_min,y_max,y_min,z_max,z_min
integer :: tag_counter
logical,intent(in) :: initial
!integer :: x_min_mod,y_min_mod,z_min_mod,x_max_mod,y_max_mod,z_max_mod
real(kind=dp) :: box_lo_geo(dimensions),box_hi_geo(dimensions)

CHARACTER(kind = c_char),intent(in) :: settag,cleartag
character(kind=c_char),intent(inout) :: tag(tag_lo(1):tag_hi(1),&
  tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
!TYPE(amrex_tagboxarray) :: tag
!TYPE(amrex_mfiter) :: mfi
!type(amrex_box) :: bex

!call amrex_mfiter_build(mfi,state(lvl),tiling=.false.)
!bex = mfi%validbox()
!tagarr => tag%dataptr(mfi)

!write(*,*) 'starting initial multigrid creation',lvl,amrex_max_level
if (initial) then
box_lo_geo = get_real_coords(amrex_geom(lvl),&
  (/dox_lo(1),dox_lo(2),dox_lo(3)/))
box_hi_geo = get_real_coords(amrex_geom(lvl),&
  (/dox_hi(1),dox_hi(2),dox_hi(3)/))

tag_counter = 0
!level_diff = box_level - lvl
do i = 1,num_boi
!  write(*,*) 'initial multigrid',self,lvl,box_lvl(i),boi_lo,boi_hi
!  write(*,*) 'grid bounds',amrex_problo,amrex_probhi,rdx(lvl)
  if (lvl > box_lvl(i)) return

  x_lo_node_start = FLOOR((boi_lo(1,i)-amrex_problo(1))/rdx(lvl))
  y_lo_node_start = FLOOR((boi_lo(2,i)-amrex_problo(2))/rdx(lvl))
  z_lo_node_start = FLOOR((boi_lo(3,i)-amrex_problo(3))/rdx(lvl))
  x_hi_node_diff = CEILING((amrex_probhi(1)-boi_hi(1,i))/rdx(lvl))
  y_hi_node_diff = CEILING((amrex_probhi(2)-boi_hi(2,i))/rdx(lvl))
  z_hi_node_diff = CEILING((amrex_probhi(3)-boi_hi(3,i))/rdx(lvl))

  x_hi_node_end = nodal_domain_limits(1,2,lvl)-x_hi_node_diff
  y_hi_node_end = nodal_domain_limits(2,2,lvl)-y_hi_node_diff
  z_hi_node_end = nodal_domain_limits(3,2,lvl)-z_hi_node_diff

!    write(*,*) 'num node differences before adjustment',self,lvl,x_lo_node_start,y_lo_node_start,&
!    z_lo_node_start,x_hi_node_end,y_hi_node_end,z_hi_node_end
!
! For if the box level is finer than the current refinement level
!
  if (box_lvl(i) > lvl) then
    x_lo_node_start = x_lo_node_start - (box_lvl(i)-lvl)
    y_lo_node_start = y_lo_node_start - (box_lvl(i)-lvl)
    z_lo_node_start = z_lo_node_start - (box_lvl(i)-lvl)

    x_hi_node_end = x_hi_node_end + (box_lvl(i)-lvl)
    y_hi_node_end = y_hi_node_end + (box_lvl(i)-lvl)
    z_hi_node_end = z_hi_node_end + (box_lvl(i)-lvl)

  end if


  if (x_lo_node_start < nodal_domain_limits(1,1,lvl)) x_lo_node_start = nodal_domain_limits(1,1,lvl)
  if (y_lo_node_start < nodal_domain_limits(2,1,lvl)) y_lo_node_start = nodal_domain_limits(2,1,lvl)
  if (z_lo_node_start < nodal_domain_limits(3,1,lvl)) z_lo_node_start = nodal_domain_limits(3,1,lvl)

  if (x_hi_node_end > nodal_domain_limits(1,2,lvl)) x_hi_node_end = nodal_domain_limits(1,2,lvl)
  if (y_hi_node_end > nodal_domain_limits(2,2,lvl)) y_hi_node_end = nodal_domain_limits(2,2,lvl)
  if (z_hi_node_end > nodal_domain_limits(3,2,lvl)) z_hi_node_end = nodal_domain_limits(3,2,lvl)

!  write(*,*) 'num node differences after adjustment',self,lvl,x_lo_node_start,y_lo_node_start,&
!    z_lo_node_start,x_hi_node_end,y_hi_node_end,z_hi_node_end,self,lvl
!  write(*,*) 'tag box limits',self,lvl,tag_lo,tag_hi,self,lvl
!
! Keep it in the box!
!
  if (x_lo_node_start >= dox_lo(1)) then
    x_min = x_lo_node_start
  else
    x_min = dox_lo(1)
  end if
  if (y_lo_node_start >= dox_lo(2)) then
    y_min = y_lo_node_start
  else
    y_min = dox_lo(2)
  end if
  if (z_lo_node_start >= dox_lo(3)) then
    z_min = z_lo_node_start
  else
    z_min = dox_lo(3)
  end if
!
! A very high end box
!
  if (x_hi_node_end <= dox_hi(1)) then
    x_max = x_hi_node_end
  else
    x_max = dox_hi(1)
  end if
  if (y_hi_node_end <= dox_hi(2)) then
    y_max = y_hi_node_end
  else
    y_max = dox_hi(2)
  end if
  if (z_hi_node_end <= dox_hi(3)) then
    z_max = z_hi_node_end
  else
    z_max = dox_hi(3)
  end if

!write(*,*) 'max and min values',self,lvl,x_min,y_min,z_min,x_max,y_max,z_max,self,lvl
!
!
!
  if (x_min >= x_max .OR. y_min>=y_max .OR. z_min >= z_max) then
!    write(*,*) 'ranges bad, returning',x_min,x_max,y_min,y_max,z_min,z_max
    return
  end if

  do c = z_min,z_max
    do b = y_min,y_max
      do a = x_min,x_max

        tag(a,b,c) = settag
        tag_counter = tag_counter + 1
      end do
    end do
  end do

end do

else
  do c = dox_lo(3),dox_hi(3)
    do b = dox_lo(2),dox_hi(2)
      do a = dox_lo(1),dox_hi(1)
        tag(a,b,c) = settag
!        if () then
!
!        end if

      end do
    end do
  end do

end if

!write(*,*) 'number of nodes tagged ',tag_counter,self,lvl

end subroutine
!call amrex_mfiter_destroy(mfi)
!if (i == 0) then
!
!
!
!else
!
!end if

!  if (x_lo_node_diff <= 3) then
!    x_min = tag_lo(1)
!  !else if (x_lo_node_diff ==1) then
!  !  x_min = lo(1)
!  else
!    x_min = tag_lo(1)+x_lo_node_diff-2
!  end if
!
!  if (y_lo_node_diff <= 3) then
!    y_min = tag_lo(2)
!  !else if (y_lo_node_diff ==1) then
!  !  y_min = -1
!  else
!    y_min = tag_lo(2)+y_lo_node_diff-2
!  end if
!
!  if (z_lo_node_diff <= 3) then
!    z_min = tag_lo(3)
!  !else if (z_lo_node_diff ==1) then
!  !  z_min = -1
!  else
!    z_min = tag_lo(3)+y_lo_node_diff-2
!  end if
!
!  if (x_hi_node_diff < 3) then
!    x_max = tag_hi(1)
!  !else if (x_hi_node_diff ==1) then
!  !  x_max = -1
!  else
!    x_max = tag_hi(1)-x_hi_node_diff+2
!  end if
!
!  if (y_hi_node_diff < 3) then
!    y_max = tag_hi(2)
!  !else if (x_hi_node_diff ==1) then
!  !  x_max = -1
!  else
!    y_max = tag_hi(2)-y_hi_node_diff+2
!  end if
!
!  if (z_hi_node_diff < 3) then
!    z_max = tag_hi(3)
!  !else if (x_hi_node_diff ==1) then
!  !  x_max = -1
!  else
!    z_max = tag_hi(3)-z_hi_node_diff+2
!  end if
