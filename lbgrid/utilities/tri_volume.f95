subroutine tri_volume(ray_crossing,ray_end,x1_in,y1_in,z1_in,&
     x2_in,y2_in,z2_in,x3_in,y3_in,z3_in,cross_tri)
!
! Determines if the ray crosses through the triangle.  It does this by
! effectively moving the node to the origin and the triangle with it by
! that amount.  This makes the calculation MUCH easier.  The volumes are
! then calculated, if they are either all positive or all negative then
! the ray passes through the triangle
!
! Called by: fluid_or_solid
! Calls:
!
USE precise
IMPLICIT NONE
REAL(KIND=dp),INTENT(IN) :: x1_in,y1_in,z1_in,x2_in,y2_in,&
  z2_in,x3_in,y3_in,z3_in
REAL(KIND=dp) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
REAL(KIND=dp) :: ray_pointx,ray_pointy,ray_pointz
INTEGER :: cross_tri,zero_counter,i

REAL(KIND=dp), DIMENSION(3) :: ray_crossing,ray_end,vol
!
! Move all the points to the origin
!
x1 = x1_in-ray_end(1)
x2 = x2_in-ray_end(1)
x3 = x3_in-ray_end(1)
y1 = y1_in-ray_end(2)
y2 = y2_in-ray_end(2)
y3 = y3_in-ray_end(2)
z1 = z1_in-ray_end(3)
z2 = z2_in-ray_end(3)
z3 = z3_in-ray_end(3)
ray_pointx = ray_crossing(1)-ray_end(1)
ray_pointy = ray_crossing(2)-ray_end(2)
ray_pointz = ray_crossing(3)-ray_end(3)
cross_tri=0
!
! This is much easier than doing a 4x4 matrix
!
! Cross product
!     | rp_x rp_y rp_z |
!     | p1_x p1_y p1_z |
!     | p2_x p2_y p2_z |
!
!   volume = rp_x*p1_y*p2_z + rp_y*p1_z*p2_x + rp_z*p1_x*p2_y
!     - rp_x*p1_z*p2_y - rp_y*p1_x*p2_Z - rp_z*p1_y*p2_x
!
! Beware pivots and going in order, i.e.
!
!
!
!      vol(1)= ray_pointx*y1*z2 + ray_pointy*z1*x2 + ray_pointz*x1*y2 - &
!         ray_pointx*z1*y2 - ray_pointy*x1*z2 - ray_pointz*y1*x2
!      vol(2)= ray_pointx*y2*z3 + ray_pointy*z2*x3 + ray_pointz*x2*y3 - &
!         ray_pointx*z2*y3 - ray_pointy*x2*z3 - ray_pointz*y2*x3
!      vol(3)= ray_pointx*y3*z1 + ray_pointy*z3*x1 + ray_pointz*x3*y1 - &
!         ray_pointx*z3*y1 - ray_pointy*x3*z1 - ray_pointz*y3*x1
!vol(1) = ray_pointx*(y2*z3 - z2*y3) - ray_pointy*(z2*x3 - x2*z3) + &
!         ray_pointz*(x2*y3 - x3*y2)
!
!vol(2) = -ray_pointx*(y1*z3 - z1*y3) + ray_pointy*(x3*z1-x1*z3) - &
!         ray_pointz*(x1*y3 - x3*y1)
!
!vol(3) = ray_pointx*(y1*z2 - z1*y2) - ray_pointy*(x2*z1 - z2*x1) + &
!         ray_pointz*(x1*y2 - y1*x2)
! volume 1, nodes 1 and 2, replace row 3, 2 row exchanges
vol(1) = ray_pointx*(y1*z2-z1*y2) + ray_pointy*(z1*x2-x1*z2) + &
  ray_pointz*(x1*y2-y1*x2)
! volume 2, nodes 1 and 3, replace row 2, 1 row exchanges
vol(2) = -ray_pointx*(y1*z3-z1*y3) - ray_pointy*(z1*x3-x1*z3) - &
  ray_pointz*(x1*y3-y1*x3)
! volume 3, nodes 2 and 3, replace row 1, 0 row exchanges
vol(3) = ray_pointx*(y2*z3-z2*y3) + ray_pointy*(z2*x3-x2*z3) + &
  ray_pointz*(x2*y3-y2*x3)
!
zero_counter = 0
! Count the 0 volumes
DO i = 1,3
  IF (vol(i) >= minus_epsilon .AND. vol(i) <= plus_epsilon) THEN
    zero_counter = zero_counter + 1
  END IF
END DO
!
! Determine if or how it crosses the triangle
!
IF ((vol(1) > 0.0D0) .AND. (vol(2) > 0.0D0) .AND. (vol(3) > 0.0D0)) THEN
!      IF ((vol(1) > minus_epsilon) .AND. (vol(2) > minus_epsilon) &
!        .AND. (vol(3) > minus_epsilon)) THEN
  cross_tri = 1
!        WRITE(*,*) 'Successful cross'
ELSE IF ((vol(1) <0.0D0) .AND. (vol(2)<0.0D0) .AND. (vol(3)<0.0D0)) THEN
!      ELSE IF ((vol(1) < plus_epsilon) .AND. (vol(2)<plus_epsilon) &
!        .AND. (vol(3)<plus_epsilon)) THEN
  cross_tri = 1
!        WRITE(*,*) 'Successful cross'
ELSE IF (zero_counter == 1) THEN
  cross_tri = -1
!        WRITE(*,*) 'Ray coincident with edge'
ELSE IF (zero_counter == 2) THEN
  cross_tri = -2
!        WRITE(*,*) 'Ray coincident with vertex'
ELSE IF (zero_counter == 3) THEN
  cross_tri = -3
!        WRITE(11,*) 'Degenerate ray'
ELSE
  cross_tri = 0
END IF

END SUBROUTINE
