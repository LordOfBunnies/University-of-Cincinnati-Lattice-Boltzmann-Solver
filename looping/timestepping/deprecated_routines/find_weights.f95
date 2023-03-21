      SUBROUTINE find_weights(wi,temper)
!
! This finds the weights to be used in the collision step
!
!
! Called by: collision
! Calls:
!
      USE precise
      USE grid_data
      use linkwise
      USE constants
      IMPLICIT NONE
      REAL(kind=dp) :: wi(dir),temper,wx,wy,wz,w_sum
      REAL(kind=dp) :: w100,w110,w111,w200,w220,w222,w300,w0,w1,w2,w3
      INTEGER :: i

      w_sum = 0.0D0
      DO i = 1,dir
        IF (dimensions == 2) THEN
          IF (dir == 9) THEN

            DO i = 1,dir

            END DO
          ELSE IF (dir == 17) THEN
            DO i = 1,dir

            END DO
          ELSE IF (dir == 25) THEN
            DO i = 1,dir

            END DO
          END IF
        ELSE IF (dimensions == 3) THEN
          IF (dir == 19) THEN
            DO i = 1,dir
              IF (cx(i) == 0) THEN
                wx = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
                w_sum = w_sum + wx
              ELSE IF (cx(i) == 1 .OR. cx(i) == -1) THEN
                wx = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
                w_sum = w_sum + wx
              END IF
              IF (cx(i) == 0) THEN
                wy = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
                w_sum = w_sum + wy
              ELSE IF (cy(i) == 1 .OR. cy(i) == -1) THEN
                wy = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
                w_sum = w_sum + wy
              END IF
              IF (cz(i) == 0) THEN
                wz = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
              ELSE IF (cz(i) == 1 .OR. cz(i) == -1) THEN
                wz = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
              END IF
              wi(i) = wx*wy*wz
            END DO
          ELSE IF (dir == 31) THEN
            DO i = 1,dir

            END DO
          ELSE IF (dir == 39) THEN
            ! Cover the bases for 0, 1, 2, and 3 links
            w0 = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
            w1 = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
            w2 = temper*(-3.0D0 +10.0D0*temper - 5.0D0*temper**2)/40.0D0
            w3 = temper*(4.0D0 - 15.0D0*temper +15.0D0*temper**2)/720.0D0
            ! for ci = (0,0,0)
              !w000
            ! for ci = (+/-1, 0, 0)
              w100 = w1 * w0**2
            ! for ci = (+/-1, +/-1, +/-1)
              w111 = w1**3
            ! for ci = (+/-2, 0, 0)
              w200 = w2 * w0**2
            ! for ci = (+/-2, +/-2, 0)
              w220 = w0 * w2**2
            ! for ci = (+/-3, 0, 0)
              w300 = w0**2 * w3
            DO i = 1,dir
!              IF (cx(i) == 0) THEN
!                wx = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
!              ELSE IF (cx(i) == 1 .OR. cx(i) == -1) THEN
!                wx = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
!              ELSE IF (cx(i) == 2 .OR. cx(i) == -2) THEN
!                wx = temper*(-3.0D0 +10.0D0*temper - 5.0D0*temper**2)/40.0D0
!              ELSE IF (cx(i) == 3 .OR. cx(i) == -3) THEN
!                wx = temper*(4.0D0 - 15.0D0*temper +15.0D0*temper**2)/720.0D0
!              END IF
!              IF (cx(i) == 0) THEN
!                wy = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
!              ELSE IF (cy(i) == 1 .OR. cy(i) == -1) THEN
!                wy = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
!              ELSE IF (cy(i) == 2 .OR. cy(i) == -2) THEN
!                wy = temper*(-3.0D0 +10.0D0*temper - 5.0D0*temper**2)/40.0D0
!              ELSE IF (cy(i) == 3 .OR. cy(i) == -3) THEN
!                wy = temper*(4.0D0 - 15.0D0*temper +15.0D0*temper**2)/720.0D0
!              END IF
!              IF (cz(i) == 0) THEN
!                wz = (36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0
!              ELSE IF (cz(i) == 1 .OR. cz(i) == -1) THEN
!                wz = temper*(12.0D0 -13.0D0*temper +5.0D0*temper**2)/16.0D0
!              ELSE IF (cz(i) == 2 .OR. cz(i) == -2) THEN
!                wz = temper*(-3.0D0 +10.0D0*temper - 5.0D0*temper**2)/40.0D0
!              ELSE IF (cz(i) == 3 .OR. cz(i) == -3) THEN
!                wz = temper*(4.0D0 - 15.0D0*temper +15.0D0*temper**2)/720.0D0
!              END IF
!              wi(i) = wx*wy*wz
              IF (i == 1) THEN
                wi(i) = ((36.0D0 - 49.0D0*temper + 42.0D0*temper**2 - 15.0D0*temper**3)/36.0D0)**3
              ELSE IF (i > 1 .AND. i <= 7) THEN
                wi(i) = w100
              ELSE IF (i > 7 .AND. i <= 15) THEN
                wi(i) = w111
              ELSE IF (i > 15 .AND. i <= 21) THEN
                wi(i) = w200
              ELSE IF (i > 21 .AND. i <= 33) THEN
                wi(i) = w220
              ELSE IF (i > 33) THEN
                wi(i) = w300
              END IF

            END DO


          END IF
        END IF

      END DO

      END SUBROUTINE
