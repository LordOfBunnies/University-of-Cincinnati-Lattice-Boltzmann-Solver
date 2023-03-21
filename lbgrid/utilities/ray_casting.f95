      SUBROUTINE ray_casting(ray_leng,ray_end)
!
! Creates random numbers to make a ray much longer than the largest
! direction of the bounding box.  This keeps it from inappropriately
! crossing things.
!
! Called by: fluid_or_solid
! Calls: DATE_AND_TIME,RANDOM_SEED,RANDOM_NUMBER
!
      USE precise
!      USE nml_bounding_box
      IMPLICIT NONE
      REAL(KIND=dp) :: ray_leng,d2
      REAL(KIND=dp) :: ray_end(3)
      INTEGER :: c,i
      
      INTEGER :: values(8),a(8)
      INTEGER,ALLOCATABLE :: b(:)


      CALL DATE_AND_TIME(VALUES=values)
      CALL RANDOM_SEED(SIZE = c)
      ALLOCATE(b(c))
      DO i = 1,c
        b(i) = values(6)+values(7)**i+ values(8)*i
      END DO
      CALL RANDOM_SEED(PUT=b)
      CALL RANDOM_NUMBER(ray_end)
      DO i = 1,3
        ray_end(i)=ray_end(i)-0.5D0
      END DO
      d2 = ray_end(1)**2+ray_end(2)**2+ray_end(3)**2
!      WRITE(*,*) d2,ray_leng,ray_end(1),ray_end(2),ray_end(3)
!      IF (bounding_method > 0) THEN
        ray_end(1) = ray_leng*ray_end(1)/SQRT(d2)
        ray_end(2) = ray_leng*ray_end(2)/SQRT(d2)
        ray_end(3) = ray_leng*ray_end(3)/SQRT(d2)
!      ELSE IF (bounding_method == -1) THEN
!        ray_end(1) = ray_leng*ray_end(1)/SQRT(d2)
!        ray_end(2) = ray_leng*ray_end(2)/SQRT(d2)
!        ray_end(3) = ray_end(3)
!      ELSE IF (bounding_method == -2) THEN
!        ray_end(1) = ray_leng*ray_end(1)/SQRT(d2)
!        ray_end(2) = ray_end(2)
!        ray_end(3) = ray_leng*ray_end(3)/SQRT(d2)
!      ELSE IF (bounding_method == -3) THEN
!        ray_end(1) = ray_end(1)
!        ray_end(2) = ray_leng*ray_end(2)/SQRT(d2)
!        ray_end(3) = ray_leng*ray_end(3)/SQRT(d2)
!      END IF
!      WRITE(*,*) ray_end(1),ray_end(2),ray_end(3)
      DEALLOCATE(b)
!
      END SUBROUTINE
