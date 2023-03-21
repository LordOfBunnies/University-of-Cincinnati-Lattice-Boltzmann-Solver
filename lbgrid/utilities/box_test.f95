SUBROUTINE box_test(ray_point1,ray_point2,max_x,min_x,max_y,&
     min_y,max_z,min_z,outside)
!
! Tests whether a ray will pass through the region of the geometry
! (either whole geometry or just a triangle)
! This is meant to short circuit the need to test all the rays and
! points and save time and processing.
!
! Called by: fluid_or_solid
! Calls:
!
USE precise
IMPLICIT NONE
REAL(KIND=dp), INTENT(IN) :: max_x,max_y,max_z,min_x,min_y,min_z
INTEGER, INTENT(OUT) :: outside
REAL(KIND=dp), INTENT(IN) :: ray_point1(3), ray_point2(3)
!
! Repeated checks to find if the ray lies completely outside the geometry
!
IF ((ray_point1(1) >= max_x + plus_epsilon ) .AND. (ray_point2(1) >=&
   max_x + plus_epsilon)) THEN
  outside = 1
  RETURN
END IF
IF ((ray_point1(1) <= min_x + minus_epsilon) .AND. (ray_point2(1) <= &
  min_x + minus_epsilon)) THEN
  outside = 1
  RETURN
END IF
IF ((ray_point1(2) >= max_y + plus_epsilon) .AND. (ray_point2(2) >= &
  max_y + plus_epsilon)) THEN
  outside = 1
  RETURN
END IF
IF ((ray_point1(2) <= min_y + minus_epsilon) .AND. (ray_point2(2) <= &
  min_y + minus_epsilon)) THEN
  outside = 1
  RETURN
END IF
IF ((ray_point1(3) >= max_z + plus_epsilon) .AND. (ray_point2(3) >= &
  max_z + plus_epsilon)) THEN
  outside = 1
  RETURN
END IF
IF ((ray_point1(3) <= min_z + minus_epsilon) .AND. (ray_point2(3) <= &
  min_z + minus_epsilon)) THEN
  outside = 1
  RETURN
END IF

!IF ((ray_point1(1) >= max_x ) .AND. (ray_point2(1) >= max_x )) THEN
!  outside = 1
!  RETURN
!END IF
!IF ((ray_point1(1) <= min_x) .AND. (ray_point2(1) <= min_x)) THEN
!  outside = 1
!  RETURN
!END IF
!IF ((ray_point1(2) >= max_y) .AND. (ray_point2(2) >=  max_y)) THEN
!  outside = 1
!  RETURN
!END IF
!IF ((ray_point1(2) <= min_y) .AND. (ray_point2(2) <= min_y)) THEN
!  outside = 1
!  RETURN
!END IF
!IF ((ray_point1(3) >= max_z) .AND. (ray_point2(3) >= max_z)) THEN
!  outside = 1
!  RETURN
!END IF
!IF ((ray_point1(3) <= min_z) .AND. (ray_point2(3) <= min_z)) THEN
!  outside = 1
!  RETURN
!END IF
!
! If it makes it all the way to the end then it passes through the region
! of the geometry in question (whole geometry or triangle)
!
outside = 0

RETURN
END


