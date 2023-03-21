SUBROUTINE ray_casting_away(ray_length,ray_end1,ray_end2)
!
! This is a modification to the ray casting algorithim. Instead of
! casting the ray in a random direction, it is cast away from the grid
! center. This is to allow better adjusting to
!
! Calls:
! Called by: fluid_or_solid
!
USE precise
USE grid_data
IMPLICIT NONE
!      INTEGER,INTENT(IN) :: ray_end1(3),ray_length
REAL(KIND = dp),INTENT(IN) :: ray_length,ray_end1(3)
REAL(KIND = dp) :: dist_to_pt,ray_end2(3),&
  ray_delta(3)!,very_small_number

!      very_small_number = 1.0D-6

!      center(1) = (x_max+x_min)/2.0D0
!      center(2) = (y_max+y_min)/2.0D0
!      center(3) = (z_max+z_min)/2.0D0

ray_delta(1) = ray_end1(1) - grid_center(1)
ray_delta(2) = ray_end1(2) - grid_center(2)
ray_delta(3) = ray_end1(3) - grid_center(3)


dist_to_pt = SQRT(ray_delta(1)**2 + ray_delta(2)**2 &
  + ray_delta(3)**2)

ray_end2(1) = ray_end1(1) +ray_length*ray_delta(1)/dist_to_pt
ray_end2(2) = ray_end1(2) +ray_length*ray_delta(2)/dist_to_pt
ray_end2(3) = ray_end1(3) +ray_length*ray_delta(3)/dist_to_pt


END SUBROUTINE ray_casting_away


!      IF (ABS(dist_to_pt) <= very_small_number) THEN
!
!
!
!      ELSE
!        ray_end2(1) = ray_end1(1)*ray_length*ABS(ray_delta(1)/dist_to_pt
!        ray_end2(2) = ray_end1(2)*ray_length*ray_delta(2)/dist_to_pt
!        ray_end2(3) = ray_end1(3)*ray_length*ray_delta(3)/dist_to_pt
!      END IF
