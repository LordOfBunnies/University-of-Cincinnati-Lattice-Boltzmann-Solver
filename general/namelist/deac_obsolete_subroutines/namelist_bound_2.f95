      SUBROUTINE namelist_bound_2
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
      REAL(KIND=dp) :: x_lower_lim,y_lower_lim,z_lower_lim,&
            x_upper_lim,y_upper_lim,z_upper_lim
! Declare the NAMELIST
      NAMELIST /BOUND_SHAPE/ x_lower_lim,y_lower_lim,z_lower_lim,&
            x_upper_lim,y_upper_lim,z_upper_lim
! Initialize the variables
      x_lower_lim = -1.0D0
      y_lower_lim = -1.0D0
      z_lower_lim = -1.0D0
      x_upper_lim = 1.0D0
      y_upper_lim = 1.0D0
      z_upper_lim = 1.0D0
! Read the data
      READ(2,BOUND_SHAPE,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at bounding box read type ',&
          bounding_method
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      bounds(1) = x_lower_lim
      bounds(2) = y_lower_lim
      bounds(3) = z_lower_lim
      bounds(4) = x_upper_lim
      bounds(5) = y_upper_lim
      bounds(6) = z_upper_lim
!
! Write the bounding box data to the output file.
!
      WRITE(11,101) x_lower_lim
      WRITE(11,102) y_lower_lim
      WRITE(11,103) z_lower_lim
      WRITE(11,104) x_upper_lim
      WRITE(11,105) y_upper_lim
      WRITE(11,106) z_upper_lim
      WRITE(11,*) ''
 101  FORMAT ('X lower limit = ', F9.6)
 102  FORMAT ('Y lower limit = ', F9.6)
 103  FORMAT ('Z lower limit = ', F9.6)
 104  FORMAT ('X upper limit = ', F9.6)
 105  FORMAT ('Y upper limit = ', F9.6)
 106  FORMAT ('Z upper limit = ', F9.6)
! 10  FORMAT (F9.7)

      END SUBROUTINE
