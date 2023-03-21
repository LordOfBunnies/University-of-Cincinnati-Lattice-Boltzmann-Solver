      SUBROUTINE namelist_output_3(shape_save,save_counter)
!
! Call to set up a rectangular output region
!
! Called by: namelist_output
! Calls: error_out
!
      USE nml_output
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: shape_save,istat,save_counter
      REAL(KIND=dp) :: x_center,y_center,z_center,&
            x_multiplier,y_multiplier,z_multiplier
! Declare the NAMELIST
      NAMELIST /SAVE_REGION/ x_center,y_center,z_center,&
            x_multiplier,y_multiplier,z_multiplier
! Initialize the data
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      x_multiplier = 2.0D0
      y_multiplier = 2.0D0
      z_multiplier = 2.0D0
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
      save_regions(4,save_counter) = x_multiplier
      save_regions(5,save_counter) = y_multiplier
      save_regions(6,save_counter) = z_multiplier
!
! Write the output region data to the output file.
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
