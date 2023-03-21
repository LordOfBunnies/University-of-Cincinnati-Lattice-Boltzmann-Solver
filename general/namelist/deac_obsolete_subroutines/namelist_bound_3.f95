      SUBROUTINE namelist_bound_3
!
! Called to set up a rectangular bounding box
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
            x_multiplier,y_multiplier,z_multiplier
! Declare the NAMELIST
      NAMELIST /BOUND_SHAPE/ x_center,y_center,z_center,&
        x_multiplier,y_multiplier,z_multiplier
! Initialize the variables
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      x_multiplier = 2.0D0
      y_multiplier = 2.0D0
      z_multiplier = 2.0D0
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
      bounds(4) = x_multiplier
      bounds(5) = y_multiplier
      bounds(6) = z_multiplier
!
! Write the bounding box data to the output file.
!
      WRITE(11,101) x_center
      WRITE(11,102) y_center
      WRITE(11,103) z_center
      WRITE(11,104) x_multiplier
      WRITE(11,105) y_multiplier
      WRITE(11,106) z_multiplier
      WRITE(11,*) ''
 101  FORMAT ('X center = ', F9.7)
 102  FORMAT ('Y center = ', F9.7)
 103  FORMAT ('Z center = ', F9.7)
 104  FORMAT ('X multiplier = ', F9.7)
 105  FORMAT ('Y multiplier = ', F9.7)
 106  FORMAT ('Z multiplier = ', F9.7)
      END SUBROUTINE
