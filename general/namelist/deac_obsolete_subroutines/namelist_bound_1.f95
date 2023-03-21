      SUBROUTINE namelist_bound_1
!
! Call to set up a rectangular bounding box
!
! Called by: namelist_bounding_box
! Calls: error_out
!
      USE nml_bounding_box
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: istat
      REAL(KIND=dp) :: x_multiplier,y_multiplier,z_multiplier
! Declare the NAMELIST
      NAMELIST /BOUND_SHAPE/ x_multiplier,y_multiplier,z_multiplier
! Initialize the variables
      x_multiplier = 2.0D0
      y_multiplier = 2.0D0
      z_multiplier = 2.0D0
! Read the data
      READ(2,BOUND_SHAPE,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at bounding box read type ',&
          bounding_method
        WRITE(11,404) bounding_method
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      bounds(1) = x_multiplier
      bounds(2) = y_multiplier
      bounds(3) = z_multiplier
!
! Write out the data to the output file.
!
      WRITE(11,100) x_multiplier
      WRITE(11,101) y_multiplier
      WRITE(11,102) z_multiplier
      WRITE(11,*) ''

 100  FORMAT ('X multiplier = ', F11.5)
 101  FORMAT ('Y multiplier = ', F11.5)
 102  FORMAT ('Z multiplier = ', F11.5)
 404  FORMAT ('Error reading namelist at bounding box read type ',I4)
      END SUBROUTINE
