      SUBROUTINE namelist_init_region_20s(region_shape,region_counter)
!
! Call to set up a cylidrical init region size
!
! Called by: namelist_init
! Calls: error_out
!
      USE nml_init
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: region_shape,istat,region_counter
      REAL(KIND=dp) :: x_center,y_center,z_center,radius,&
            distance
! Declare the NAMELIST
      NAMELIST /INIT_REGION/ x_center,y_center,z_center,radius,&
            distance
! Initialize the data
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      radius = 1.0D0
      distance = 1.0D0
! Read the data
      READ(2,INIT_REGION,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at init read type ',&
          region_shape, ' at region ', region_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      init_bounds(1,region_counter) = x_center
      init_bounds(2,region_counter) = y_center
      init_bounds(3,region_counter) = z_center
      init_bounds(4,region_counter) = radius
      init_bounds(5,region_counter) = distance
      init_bounds(6,region_counter) = 0.0D0
!
! Write the initialization region data to the output file
!
      WRITE(11,201) x_center
      WRITE(11,202) y_center
      WRITE(11,203) z_center
      WRITE(11,204) radius
!
!
      IF (region_shape == 21) THEN
        WRITE(11,101) distance
      ELSE IF (region_shape == 22) THEN
        WRITE(11,102) distance
      ELSE IF (region_shape == 23) THEN
        WRITE(11,103) distance
      ELSE IF (region_shape == 24) THEN
        WRITE(11,104) distance
      ELSE IF (region_shape == 25) THEN
        WRITE(11,105) distance
      ELSE IF (region_shape == 26) THEN
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
