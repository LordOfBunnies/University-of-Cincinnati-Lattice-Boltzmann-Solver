      SUBROUTINE namelist_bound_20s
!
! Called to set up various cylidrical bounding boxes
!
! Called by: namelist_bounding_box
! Calls: error_out
!
      USE nml_bounding_box
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: istat
      REAL(KIND=dp) :: x_center,y_center,z_center,radius,&
            distance
! Declare the NAMELIST
      NAMELIST /BOUND_SHAPE/ x_center,y_center,z_center,radius,&
            distance
! Initialize the variables
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      radius = 1.0D0
      distance = 1.0D0
! Read the data
      READ(2,BOUND_SHAPE,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at bounding box read type ',&
          bounding_method
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign data to appropriate arrays
      bounds(1) = x_center
      bounds(2) = y_center
      bounds(3) = z_center
      bounds(4) = radius
      bounds(5) = distance
      WRITE(11,201) x_center
      WRITE(11,202) y_center
      WRITE(11,203) z_center
      WRITE(11,204) radius
!
! Write the bounding box data to the output file.
!
      IF (bounding_method == 21) THEN
        WRITE(11,101) distance
      ELSE IF (bounding_method == 22) THEN
        WRITE(11,102) distance
      ELSE IF (bounding_method == 23) THEN
        WRITE(11,103) distance
      ELSE IF (bounding_method == 24) THEN
        WRITE(11,104) distance
      ELSE IF (bounding_method == 25) THEN
        WRITE(11,105) distance
      ELSE IF (bounding_method == 26) THEN
        WRITE(11,106) distance
      END IF
      WRITE(11,*) ''
 101  FORMAT ('Distance in x-direction = ', F9.7)
 102  FORMAT ('Distance in y-direction = ', F9.7)
 103  FORMAT ('Distance in z-direction = ',F9.7)
 104  FORMAT ('Distance multiplier in x-direction = ', F9.7)
 105  FORMAT ('Distance multiplier in y-direction = ', F9.7)
 106  FORMAT ('Distance multiplier in z-direction = ', F9.7)

 201  FORMAT ('X center = ',F9.7)
 202  FORMAT ('Y center = ',F9.7)
 203  FORMAT ('Z center = ',F9.7)
 204  FORMAT ('Radius = ',F9.7)
      END SUBROUTINE 
