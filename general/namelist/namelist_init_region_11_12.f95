      SUBROUTINE namelist_init_region_11_12(region_shape,region_counter)
!
! Call to set up a spherical init region
!
! Called by: namelist_init
! Calls: error_out
!
      USE nml_init
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: region_shape,istat,region_counter
      REAL(KIND=dp) :: x_center,y_center,z_center,&
            radius
! Declare the NAMELIST
      NAMELIST /INIT_REGION/ x_center,y_center,z_center,&
            radius
! Initialize the data
      x_center = 0.0D0
      y_center = 0.0D0
      z_center = 0.0D0
      radius = 1.0D0
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
      init_bounds(5,region_counter) = 0.0D0
      init_bounds(6,region_counter) = 0.0D0
!
! Write the initialization region data to the output file
!
      WRITE(11,101) x_center
      WRITE(11,102) y_center
      WRITE(11,103) z_center
      IF (region_shape == 11) THEN
        WRITE(11,104) radius
      ELSE
        WRITE(11,105) radius
      END IF
      WRITE(11,*) ''
 101  FORMAT ('X coordinate of center = ',F9.7)
 102  FORMAT ('Y coordinate of center = ',F9.7)
 103  FORMAT ('Z coordinate of center = ',F9.7)
 104  FORMAT ('Radius multiplier = ',F9.7)
 105  FORMAT ('Radius = ',F9.7)

      END SUBROUTINE
