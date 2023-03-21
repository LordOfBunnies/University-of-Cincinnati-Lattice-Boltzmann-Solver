      SUBROUTINE closeness_checker(ray_end1,x1_in,y1_in,z1_in,&
           x2_in,y2_in,z2_in,x3_in,y3_in,z3_in,closeness)
!
!
!
!
!
!
!
      USE precise
      IMPLICIT NONE
      REAL(KIND = dp),INTENT(IN) :: x1_in,y1_in,z1_in,&
           x2_in,y2_in,z2_in,x3_in,y3_in,z3_in,ray_end1(3)
      REAL(KIND = dp) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
      INTEGER :: i
      LOGICAL :: closeness

      x1 = ray_end1(1)-x1_in
      x2 = ray_end1(1)-x2_in
      x3 = ray_end1(1) - x3_in
      y1 = ray_end1(2) - y1_in
      y2 = ray_end1(2) - y2_in
      y3 = ray_end1(2) - y3_in
      z1 = ray_end1(3) - z1_in
      z2 = ray_end1(3) - z2_in
      z3 = ray_end1(3) - z3_in

!      DO i=1,3
      IF (ABS(x1) <= sp_plus_epsilon .AND. ABS(y1) <= sp_plus_epsilon &
        .AND. ABS(z1) <= sp_plus_epsilon) THEN

        closeness = .TRUE.
        RETURN
      END IF
      IF (ABS(x2) <= sp_plus_epsilon .AND. ABS(y2) <= sp_plus_epsilon &
        .AND. ABS(z2) <= sp_plus_epsilon) THEN

        closeness = .TRUE.
        RETURN
      END IF
      IF (ABS(x3) <= sp_plus_epsilon .AND. ABS(y3) <= sp_plus_epsilon &
        .AND. ABS(z3) <= sp_plus_epsilon) THEN

        closeness = .TRUE.
        RETURN
      END IF

!      END DO

      END SUBROUTINE closeness_checker
