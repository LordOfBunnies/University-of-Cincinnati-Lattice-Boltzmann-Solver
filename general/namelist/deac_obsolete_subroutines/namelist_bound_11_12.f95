      SUBROUTINE namelist_bound_11_12
!
! Called to set up a spherical bounding box
!
! Called by: namelist_bounding_box
! Calls: error_out
!
      USE nml_bounding_box
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: istat
      REAL(KIND=dp) :: x_center,y_center,z_center,&
            radius
! Declare the NAMELIST
      NAMELIST /BOUND_SHAPE/ x_center,y_center,z_center,&
            radius
! Initialize the variables
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      radius = 1.0D0
! Read the data
      READ(2,BOUND_SHAPE,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at bounding box read type ',&
          bounding_method
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      bounds(1) = x_center
      bounds(2) = y_center
      bounds(3) = z_center
      bounds(4) = radius
!
! Write the bounding box data to the output file.
!
      WRITE(11,101) x_center
      WRITE(11,102) y_center
      WRITE(11,103) z_center
      IF (bounding_method == 11) THEN
        WRITE(11,104) radius
      ELSE
        WRITE(11,105) radius
      END IF
      WRITE(11,*) ''
 101  FORMAT ('X coordinate of center = ', F9.7)
 102  FORMAT ('Y coordinate of center = ', F9.7)
 103  FORMAT ('Z coordinate of center = ', F9.7)
 104  FORMAT ('Radius multiplier = ', F9.7)
 105  FORMAT ('Radius = ', F9.7)
      END SUBROUTINE
