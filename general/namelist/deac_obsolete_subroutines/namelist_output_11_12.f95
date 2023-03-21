      SUBROUTINE namelist_output_11_12(shape_save,save_counter)
!
! Call to set up a spherical output region
!
! Called by: namelist_output
! Calls: error_out
!
      USE nml_output
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: shape_save,istat,save_counter
      REAL(KIND=dp) :: x_center,y_center,z_center,radius
! Declare the NAMELIST
      NAMELIST /SAVE_REGION/ x_center,y_center,z_center,&
            radius
! Initialize the data
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      radius = 1.0D0
! Read the data
      READ(2,SAVE_REGION,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at output read type ',&
          shape_save, ' at region ', save_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      save_regions(1,save_counter) = x_center
      save_regions(2,save_counter) = y_center
      save_regions(3,save_counter) = z_center
      save_regions(4,save_counter) = radius
      save_regions(5,save_counter) = 0.0D0
      save_regions(6,save_counter) = 0.0D0
!
! Write the output region data to the output file
!
      WRITE(11,101) x_center
      WRITE(11,102) y_center
      WRITE(11,103) z_center
      IF (shape_save == 11) THEN
        WRITE(11,104) radius
      ELSE
        WRITE(11,105)radius
      END IF
      WRITE(11,*) ''

 101  FORMAT ('X coordinate of center = ', F9.7)
 102  FORMAT ('Y coordinate of center = ', F9.7)
 103  FORMAT ('Z coordinate of center = ', F9.7)
 104  FORMAT ('Radius multiplier = ', F9.7)
 105  FORMAT ('Radius = ', F9.7)
      END SUBROUTINE
