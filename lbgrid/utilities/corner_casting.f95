subroutine corner_casting(ray_end1,ray_end2,corner_1,corner_2,&
  min_length,first_corner_fail)
!
! This casts a ray from the starting node to one of the corners of the grid
!   and beyond.  This occurs if the edge counter becomes extreme (maybe the
!   cross counter too, we'll see)
!
! Called by: fluid_or_solid
! Calls:
! External Calls:
!
use precise
use constants
implicit none
real(kind=dp),intent(in) :: ray_end1(3),min_length,corner_1(3),corner_2(3)
real(kind=dp),intent(out) :: ray_end2(3)
real(kind=dp) :: new_ray_length
logical :: first_corner_fail

!
! If we've been through this once, ray end at the opposite corner
!
if (first_corner_fail) then
!  write(*,*) 'hrm'
  new_ray_length = SQRT((ray_end1(1)-corner_2(1))**2 + (ray_end1(2)-corner_2(2))**2 + &
    (ray_end1(3)-corner_2(3))**2)
!
! If the ray is too short
!
!  if (new_ray_length < min_length) then
!    ray_end2 = corner_2*min_length/new_ray_length
!  else
!    ray_end2 = corner_2
!  end if
call ray_casting(new_ray_length,ray_end2)
!
! put the other end of the ray in the lowest corner of the grid
!
else
  new_ray_length = SQRT((ray_end1(1)-corner_1(1))**2 + (ray_end1(2)-corner_1(2))**2 + &
    (ray_end1(3)-corner_1(3))**2)
!call ray_casting(new_ray_length,ray_end2)
    if (new_ray_length < min_length) then
    ray_end2 = corner_1*min_length/new_ray_length
  else
    ray_end2 = corner_1
  end if
end if


end subroutine
