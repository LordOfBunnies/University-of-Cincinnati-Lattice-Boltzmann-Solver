      SUBROUTINE tri_area(ray_crossing,x1,y1,z1,x2,y2,z2,&
          x3,y3,z3,triangle_area)
!
! Determines if the node is within the triangle, only called if ray
! is on the plane of the triangle
!
! Called by: fluid_or_solid
! Calls:
!
      USE precise
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(IN) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
      REAL(KIND=dp) :: ray_pointx,ray_pointy,ray_pointz
      INTEGER :: triangle_area,zero_counter=0,i
 
      REAL(KIND=dp), DIMENSION(3) :: ray_crossing,area
!
! Projects points onto (first) the x-y plane to make calculations
! easier.  It then calculates the areas, if all the areas are
! either positive or negative, the point is within the triangle
!
      area(1) = ((x1-ray_crossing(1))*(y2-ray_crossing(2))-(x2-&
        ray_crossing(1))*(y1-ray_crossing(2)))/2.0D0
      area(2) = ((x2-ray_crossing(1))*(y3-ray_crossing(2))-(x3-&
        ray_crossing(1))*(y2-ray_crossing(2)))/2.0D0
      area(3) = ((x3-ray_crossing(1))*(y1-ray_crossing(2))-(x1-&
        ray_crossing(1))*(y3-ray_crossing(2)))/2.0D0

! If all areas are equal to 0.0 then the triangle is perpindicular to the
! x-y plane, so it then projects onto the x-z plane
      IF (area(1) >= sp_minus_epsilon .AND. area(1) <= sp_plus_epsilon &
        .AND. area(2) >= sp_minus_epsilon .AND. area(2) <= sp_plus_epsilon &
        .AND. area(3) >= sp_minus_epsilon .AND. area(3) <= sp_plus_epsilon) THEN
        area(1) = ((x1-ray_crossing(1))*(z2-ray_crossing(3))-(x2-&
          ray_crossing(1))*(z1-ray_crossing(3)))/2.0D0
        area(2) = ((x2-ray_crossing(1))*(z3-ray_crossing(3))-(x3-&
          ray_crossing(1))*(z2-ray_crossing(3)))/2.0D0
        area(3) = ((x3-ray_crossing(1))*(z1-ray_crossing(3))-(x1-&
          ray_crossing(1))*(z3-ray_crossing(3)))/2.0D0
!
! If it is also perpindicular to the x-z plane it is projected onto the
! y-z plane as a last effort to solve the problem
!
        IF (area(1) >= sp_minus_epsilon .AND. area(1) <= sp_plus_epsilon &
          .AND. area(2) >= sp_minus_epsilon .AND. area(2) <= sp_plus_epsilon &
          .AND. area(3) >= sp_minus_epsilon .AND. area(3) <= sp_plus_epsilon) THEN
          area(1) = ((y1-ray_crossing(2))*(z2-ray_crossing(3))-(y2-&
            ray_crossing(2))*(z1-ray_crossing(3)))/2.0D0
          area(2) = ((y2-ray_crossing(2))*(z3-ray_crossing(3))-(y3-&
            ray_crossing(2))*(z2-ray_crossing(3)))/2.0D0
          area(3) = ((y3-ray_crossing(2))*(z1-ray_crossing(3))-(y1-&
            ray_crossing(2))*(z3-ray_crossing(3)))/2.0D0
!
! If is perpindicular to all three plane it is a degenerate triangle and
! the use needs new geometry.
!
          IF (area(1) >= sp_minus_epsilon .AND. area(1) <= sp_plus_epsilon &
            .AND. area(2) >= sp_minus_epsilon .AND. area(2) <= sp_plus_epsilon &
            .AND. area(3) >= sp_minus_epsilon .AND. area(3) <= sp_plus_epsilon) THEN
            WRITE(11,*) 'Warning, degenerate triangle.'
            CALL error_out
            triangle_area = -1000

          ELSE IF (area(1) > 0.0D0 .AND. area(2) > 0.0D0 .AND. area(3) > 0.0D0) THEN
            triangle_area = 1
!            WRITE(*,*) 'Triangle valid, 3 projections'
          ELSE IF (area(1) < 0.0D0 .AND. area(2) < 0.0D0 .AND. area(3) < 0.0D0) THEN
            triangle_area = 1
!            WRITE(*,*) 'Triangle valid, 3 projections'
          ELSE

            DO i = 1,3
              IF (area(i) >= sp_minus_epsilon .AND. area(i) <= sp_plus_epsilon) THEN
                zero_counter=zero_counter+1
              END IF
            END DO
!
! If one of the triangles' areas equals 0.0 then the point is coincident
! with an edge of the triangle
!
! If two of the triangles' areas equal 0.0 then the point is conincident
! with a vertex and bad things happen
!
            IF (zero_counter == 1) THEN
              triangle_area = -1
            ELSE IF (zero_counter == 2) THEN
              triangle_area = -2
            ELSE
              triangle_area = 0
            END IF
          END IF
!
! This process must be repeated for all planes
!
        ELSE IF (area(1) > 0.0D0 .AND. area(2) > 0.0D0 .AND. area(3) > 0.0D0) THEN
          triangle_area = 1
!          WRITE(*,*) 'Triangle valid, 2 projections'
        ELSE IF (area(1) < 0.0D0 .AND. area(2) < 0.0D0 .AND. area(3) < 0.0D0) THEN
          triangle_area = 1
!          WRITE(*,*) 'Triangle valid, 2 projections'
        ELSE

          DO i = 1,3
            IF (area(i) >= sp_minus_epsilon .AND. area(i) <= sp_plus_epsilon) THEN
              zero_counter=zero_counter+1
            END IF
          END DO
          IF (zero_counter == 1) THEN
            triangle_area = -1
          ELSE IF (zero_counter == 2) THEN
            triangle_area = -2
          ELSE
            triangle_area = 0
          END IF
        END IF
! Rinse, repeat
      ELSE IF (area(1) > 0.0D0 .AND. area(2) > 0.0D0 .AND. area(3) > 0.0D0) THEN
        triangle_area = 1
!        WRITE(*,*) 'Triangle valid, 1 projections'
      ELSE IF (area(1) < 0.0D0 .AND. area(2) < 0.0D0 .AND. area(3) < 0.0D0) THEN
        triangle_area = 1
!        WRITE(*,*) 'Triangle valid, 1 projections'
      ELSE 

        DO i = 1,3
          IF (area(i) >= sp_minus_epsilon .AND. area(i) <= sp_plus_epsilon) THEN
            zero_counter=zero_counter+1
          END IF
        END DO

        IF (zero_counter == 1) THEN
          triangle_area = -1
        ELSE IF (zero_counter == 2) THEN
          triangle_area = -2
        ELSE
          triangle_area = 0
        END IF

      END IF
      
      RETURN
      END SUBROUTINE
