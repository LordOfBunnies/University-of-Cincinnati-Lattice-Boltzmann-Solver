      SUBROUTINE namelist_error_checker(characteristic_length)
!
! This checks for potential errors before the program gets further and
! bombs the program out if it detects bad inputs.  This will eventually
! be quite extensive and there will always be new stuff to add here.
!
! Called by: namelist_global
! Calls: error_out
!
!      USE nml_bounding_box
      USE nml_geometry
      USE nml_inlet_outlet
      USE nml_output
      USE nml_init
      USE precise
      USE freestream_values
      USE constants
      USE timing
      USE startup
      USE periodic_data
      IMPLICIT NONE
      INTEGER :: i,j,k,periodic_counter,periodic_mod,counter
      REAL(KIND=dp) :: characteristic_length
      LOGICAL :: periodic_pair
!
!
!
!
!
!      IF (min_vr > max_vr) THEN
!        WRITE(*,*) 'Minimum vr greater than maximum vr, exiting.'
!        WRITE(11,*) 'Minimum vr greater than maximum vr, exiting.'
!        CALL error_out
!      END IF

      IF (dimensions /= 2 .AND. dimensions /= 3) THEN
        WRITE(*,*) 'Run is neither 2 nor 3 dimensional, exiting.'
        WRITE(11,*) 'Run is neither 2 nor 3 dimensional, exiting.'
        CALL error_out
      END IF

      IF (cells_per_length < 1) THEN
        WRITE(*,*) 'Cells per characteristic length is less than 1.'
        WRITE(11,*) 'Cells per characteristic length is less than 1.'
        CALL error_out
      ELSE IF (cells_per_length > 10000) THEN
        WRITE(*,*) 'Cells per characteristic length is too large.'
        WRITE(11,*) 'Cells per characteristic length is too large.'
        CALL error_out
      END IF

      IF (primary_flow_dir < 1 .OR. primary_flow_dir > 6) THEN
        WRITE(11,*) 'Primary flow direction does not correspond to any',&
          ', primary flow direction must be between 1 and 6.'
        CALL error_out
      END IF


!      IF (bounding_method > 20) THEN
!        IF (primary_flow_dir == 1 .OR. primary_flow_dir == 4) THEN
!          IF (bounding_method /= 21 .OR. bounding_method /= 24) THEN
!            WRITE(*,*) 'Flow direction does not match cylinder bounding',&
!              ' box direction.'
!            WRITE(11,*) 'Flow direction does not match cylinder bounding',&
!              ' box direction.'
!            CALL error_out
!          END IF
!        ELSE IF (primary_flow_dir == 2 .OR. primary_flow_dir == 5) THEN
!          IF (bounding_method /= 22 .OR. bounding_method /= 25) THEN
!            WRITE(*,*) 'Flow direction does not match cylinder bounding',&
!              ' box direction.'
!            WRITE(11,*) 'Flow direction does not match cylinder bounding',&
!              ' box direction.'
!            CALL error_out
!          END IF
!        ELSE IF (primary_flow_dir == 3 .OR. primary_flow_dir == 6) THEN
!          IF (bounding_method /= 23 .OR. bounding_method /= 26) THEN
!            WRITE(*,*) 'Flow direction does not match cylinder bounding',&
!              ' box direction.'
!            WRITE(11,*) 'Flow direction does not match cylinder bounding',&
!              ' box direction.'
!            CALL error_out
!          END IF
!        END IF
!      END IF

      IF (prandtl < 0.0) THEN
        WRITE(*,*) 'Prandtl number is less than 0.'
        WRITE(11,*) 'Prandtl number is less than 0.'
        CALL error_out
      END IF

      IF (total_time < 0.0) THEN
        WRITE(*,*) 'Total run time is negative, exiting.'
        WRITE(11,*) 'Total run time is negative, exiting.'
        CALL error_out
      END IF

      IF (temperature < 0.0) THEN
        WRITE(*,*) 'Declared temperature below absolute zero, what''s',&
          ' wrong with you!!!'
        WRITE(11,*) 'Declared temperature below absolute zero, what''s',&
          ' wrong with you!!!'
        CALL error_out
      ELSE IF (temperature > 10000.0) THEN
        WRITE(*,*) 'Temperature would dissociate almost all fluids.'
        WRITE(11,*) 'Temperature would dissociate almost all fluids.'
        CALL error_out
      END IF

      IF (pressure < 10.0) THEN
        WRITE(*,*) 'Pressure is too low, this is non-continuum flow.'
        WRITE(11,*) 'Pressure is too low, this is non-continuum flow.'
        CALL error_out
      ELSE IF (pressure > 10130000.0) THEN
        WRITE(*,*) 'Pressure greater than 100 atmospheres,',&
          ' are you sure?'
        WRITE(11,*) 'Pressure greater than 100 atmospheres,',&
          ' are you sure?'
        CALL error_out
      END IF
!
!
!
      IF (num_inlets < 1) THEN
        WRITE(*,*) 'No inlets, exiting.'
        WRITE(11,*) 'No inlets, exiting.'
        CALL error_out
      ELSE IF (num_inlets > 5 .AND. .NOT. cgns_start) THEN
        WRITE(*,*) 'Every face cannot be an inlet.'
        WRITE(11,*) 'Every face cannot be an inlet.'
        CALL error_out
      END IF

      IF (cgns_start .AND. num_inlets > 98) THEN
        WRITE(*,*) 'More than 98 inlets not supported.'
        WRITE(11,*) 'More than 98 inlets not supoorted.'
        CALL error_out
      END IF

      IF (num_outlets <1) THEN
        WRITE(*,*) 'No outlets, exiting.'
        WRITE(11,*) 'No outlets, exiting.'
        CALL error_out
      ELSE IF (num_outlets > 5 .AND. .NOT. cgns_start) THEN
        WRITE(*,*) 'Every face cannot be an outlet, exiting.'
        WRITE(11,*) 'Every face cannot be an outlet, exiting.'
        CALL error_out
      END IF

      IF (cgns_start .AND. num_outlets > 98) THEN
        WRITE(*,*) 'More than 98 outlets not supported.'
        WRITE(11,*) 'More than 98 outlets not supoorted.'
        CALL error_out
      END IF
!
!
!
      IF (mixed_start .AND. .NOT. cgns_start) THEN
        WRITE(*,*) 'Cannot be a mixed start without a CGNS start.'
        CALL error_out
      END IF
!
!
!
      IF (num_mics == 0 .AND. num_saves == 0) THEN
        WRITE(*,*) 'You''re not saving anything? Why run this?'
        WRITE(11,*) 'You''re not saving anything? Why run this?'
        CALL error_out
      END IF

      IF (v_inf < 0.0) THEN
        WRITE(*,*) 'If you want a negative velocity, specify it in',&
          ' primary flow direction.'
        WRITE(11,*) 'If you want a negative velocity, specify it in',&
          ' primary flow direction.'
        CALL error_out
      END IF


!
! Check for an even number of periodic boundary conditions
!
      periodic_sums = 0
      counter = 1
      periodic_flow_dir = 0
      DO i = 1,num_outlets

        IF (outlet_type(1,i) == 5) THEN
          periodic_sums = periodic_sums + 1
          periodic_flow_dir(counter) = outlet_type(2,i)
          counter = counter+1
        END IF

      END DO
      periodic_mod =  MODULO(periodic_sums,2)
      IF (periodic_mod == 1) THEN
        WRITE(*,*) 'Uneven number of periodic boundary conditions, exiting.'
        WRITE(11,*) 'Uneven number of periodic boundary conditions, exiting.'
        CALL error_out
      END IF
!
! Check that the right boundary condition is applied
!
!      IF (periodic_sums > 0 .AND. bounding_method > 10) THEN
!        WRITE(*,*) 'Periodic boundary conditions being applied to wrong', &
!          ' bounding box shape, exiting.'
!        WRITE(11,*) 'Periodic boundary conditions being applied to wrong', &
!          ' bounding box shape, exiting.'
!        CALL error_out
!      END IF
!
!      IF (bounding_method < 0) THEN
!        IF (bounds(2) <= bounds(1)) THEN
!          WRITE(11,400) bounds(1),bounds(2)
!          WRITE(11,*) 'Exiting program'
!          CALL error_out
!        END IF
!        IF (bounds(4) <= bounds(3)) THEN
!          WRITE(11,400) bounds(3),bounds(4)
!          WRITE(11,*) 'Exiting program'
!          CALL error_out
!        END IF
!      END IF
! Next check the periodic BCs are on opposite sides
!
!      periodic_counter = 1

!      DO i = 1,num_outlets
!        IF (outlet_type(1,i) == 5) THEN
!          periodic_pair = .FALSE.
!!          DO j = 1,num_outlets
!!          IF (outlet_type(1,j) ==5 .AND. i /= j) THEN
!!              IF (outlet_type(2,i) == 1 .AND. &
!!                    outlet_type(2,j) == 4) THEN
!!                DO k = 1,3
!!                  IF (i == k .OR. j == k)
!!                    paired_periodic_bcs( )
!!                  ELSE IF
!
!!                  END IF
!!                END DO
!
!!              ELSE IF (outlet_type(2,i) == 2 .AND. &
!!                    outlet_type(2,j) == 5) THEN
!
!
!
!!              ELSE IF (outlet_type(2,i) == 3 .AND. &
!!                    outlet_type(2,j) == 6) THEN
!          SELECT CASE (outlet_type(2,i))
!          CASE (1)
!!            IF (outlet_type(1,i) == 5) THEN
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 4) THEN
!                  periodic_pair = .TRUE.
!                END IF
!              END DO
!!            END IF
!          CASE (2)
!            DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 5) THEN
!                  periodic_pair = .TRUE.
!                END IF
!              END DO
!          CASE (3)
!            DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 6) THEN
!                  periodic_pair = .TRUE.
!                END IF
!              END DO
!          CASE (4)
!            DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 1) THEN
!                  periodic_pair = .TRUE.
!                END IF
!              END DO
!
!          CASE (5)
!            DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 2) THEN
!                  periodic_pair = .TRUE.
!                END IF
!              END DO
!          CASE (6)
!            DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 3) THEN
!                  periodic_pair = .TRUE.
!                END IF
!              END DO
!          END SELECT
!
!          IF (.NOT. periodic_pair) THEN
!            WRITE(11,200) i
!            CALL error_out
!          END IF
!!          END IF
!
!!          END DO
!
!        END IF
!
!      END DO

      CLOSE(11)

 200  FORMAT ('Periodic boundary condition not paired with matching',&
              'boundary condition. Condition in error is number ',I8)
 400  FORMAT ('For planar bounding box, axis maximum must be greater than axis ',&
        'minimum. Axis minimum is ',F10.6,' and axis maximum is ',F10.6)
      END SUBROUTINE
