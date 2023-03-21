SUBROUTINE surface_point_adjuster(ray_end1,ray_end2,normal_dir,i,source,&
  point_moved,lvl)
!
! This moves the point
!
! Calls:
! Called by: fluid_or_solid,state_directions_cgns
!
USE precise
USE grid_data
USE linkwise
IMPLICIT NONE
REAL(KIND = dp) :: ray_end1(3),ray_end2(3),direction_move(3)
REAL(KIND = dp) :: normalized_dist,ray_dist,normal_dir(3)
INTEGER :: source,i,lvl
LOGICAL :: point_moved
!
! Source 1 is fluid_or_solid
! Source 2 is state_directions_cgns
!
! 1) fluid_or_solid move it along the ray direction toward the geometric center
! 2) state_direction_cgns moves it along the lattice in the opposite direction
!   of the problem, inward in this case. It will move a fraction of the lattice
!   distance
!
IF (source == 1) THEN
!
! Determine the direction to move it, essentially walking back the ray
!
  direction_move(1) = ray_end1(1)-ray_end2(1)
  direction_move(2) = ray_end1(2)-ray_end2(2)
  direction_move(3) = ray_end1(3)-ray_end2(3)
!
! Determine the ray distance then make the movement a unit vector
!
  ray_dist = SQRT(direction_move(1)**2.0D0 + direction_move(2)**2.0D0 &
    + direction_move(3)**2.0D0)

  direction_move = direction_move/ray_dist

  ray_end1(1) = ray_end1(1) + direction_move(1)*rdx(lvl)/10.0D0
  ray_end1(2) = ray_end1(2) + direction_move(2)*rdx(lvl)/10.0D0
  ray_end1(3) = ray_end1(3) + direction_move(3)*rdx(lvl)/10.0D0


  point_moved = .TRUE.
ELSE IF (source == 2) THEN
!
! This deals with determining which triangles directions are passing through
! so deals with MUCH smaller distances. This will move the point the opposite
! direction to the direction in question, i.e. pull it towards the node itself
!
!  ray_end1(1) = ray_end1(1)+cx(opp(i))*rdx(lvl)/100.0D0
!  ray_end1(2) = ray_end1(2)+cy(opp(i))*rdx(lvl)/100.0D0
!  ray_end1(3) = ray_end1(3)+cz(opp(i))*rdx(lvl)/100.0D0
  ray_end1(1) = ray_end1(1) + sign(normal_dir(1),1.0D0)*rdx(lvl)/10.0D0
  ray_end1(2) = ray_end1(2) + sign(normal_dir(2),1.0D0)*rdx(lvl)/10.0D0
  ray_end1(3) = ray_end1(3) + sign(normal_dir(3),1.0D0)*rdx(lvl)/10.0D0

  point_moved = .TRUE.
END IF


END SUBROUTINE surface_point_adjuster
