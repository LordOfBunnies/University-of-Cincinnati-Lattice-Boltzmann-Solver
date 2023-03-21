      SUBROUTINE namelist_output_20s(shape_save,save_counter)
!
! Call to set up a cylindrical output region
!
! Called by: namelist_output
! Calls: error_out
!
      USE nml_output
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: shape_save,istat,save_counter
      REAL(KIND=dp) :: x_center,y_center,z_center,radius,&
            distance
! Declare the NAMELSIT
      NAMELIST /SAVE_REGION/ x_center,y_center,z_center,radius,&
            distance
! Initialize the data
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      radius = 1.0D0
      distance = 1.0D0
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
      save_regions(5,save_counter) = distance
      save_regions(6,save_counter) = 0.0D0
!
! Write the output region data to the output file
!
      WRITE(11,201) x_center
      WRITE(11,202) y_center
      WRITE(11,203) z_center
      WRITE(11,204) radius
!
!
      IF (shape_save == 21) THEN
        WRITE(11,101) distance
      ELSE IF (shape_save == 22) THEN
        WRITE(11,102) distance
      ELSE IF (shape_save == 23) THEN
        WRITE(11,103) distance
      ELSE IF (shape_save == 24) THEN
        WRITE(11,104) distance
      ELSE IF (shape_save == 25) THEN
        WRITE(11,105) distance
      ELSE IF (shape_save == 26) THEN
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
