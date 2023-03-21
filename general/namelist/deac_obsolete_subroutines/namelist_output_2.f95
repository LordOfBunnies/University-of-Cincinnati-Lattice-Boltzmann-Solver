      SUBROUTINE namelist_output_2(shape_save,save_counter)
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
      REAL(KIND=dp) :: x_lower_lim,y_lower_lim,z_lower_lim,&
            x_upper_lim,y_upper_lim,z_upper_lim
! Declare the NAMELIST
      NAMELIST /SAVE_REGION/ x_lower_lim,y_lower_lim,z_lower_lim,&
            x_upper_lim,y_upper_lim,z_upper_lim
! Initialize the data
      x_lower_lim = -1.0D0
      y_lower_lim = -1.0D0
      z_lower_lim = -1.0D0
      x_upper_lim = 1.0D0
      y_upper_lim = 1.0D0
      z_upper_lim = 1.0D0
! Read the data
      READ(2,SAVE_REGION,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at output read type ',&
          shape_save, ' at region ', save_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      save_regions(1,save_counter) = x_lower_lim
      save_regions(2,save_counter) = y_lower_lim
      save_regions(3,save_counter) = z_lower_lim
      save_regions(4,save_counter) = x_upper_lim
      save_regions(5,save_counter) = y_upper_lim
      save_regions(6,save_counter) = z_upper_lim
!
! Write the output region data to the output file
!
      WRITE(11,101) x_lower_lim
      WRITE(11,102) y_lower_lim
      WRITE(11,103) z_lower_lim
      WRITE(11,104) x_upper_lim
      WRITE(11,105) y_upper_lim
      WRITE(11,106) z_upper_lim
      WRITE(11,*) ''

 101  FORMAT ('X lower limit = ', F9.7)
 102  FORMAT ('Y lower limit = ', F9.7)
 103  FORMAT ('Z lower limit = ', F9.7)
 104  FORMAT ('X upper limit = ', F9.7)
 105  FORMAT ('Y upper limit = ', F9.7)
 106  FORMAT ('Z upper limit = ', F9.7)
      END SUBROUTINE
