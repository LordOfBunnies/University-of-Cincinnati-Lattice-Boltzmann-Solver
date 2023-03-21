      SUBROUTINE namelist_init_region_1(region_shape,region_counter)
!
! Call to set up a rectangular init region size
!
! Called by: namelist_init
! Calls:error_out
!
      USE nml_init
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      INTEGER :: region_shape,istat,region_counter
      REAL(KIND=dp) :: x_multiplier,y_multiplier,z_multiplier
! Declare the NAMELIST
      NAMELIST /INIT_REGION/ x_multiplier,y_multiplier,&
                  z_multiplier
! Initialize the variables
      x_multiplier = 2.0D0
      y_multiplier = 2.0D0
      z_multiplier = 2.0D0
! Read the data
      READ(2,INIT_REGION,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at init read type ',&
          region_shape, ' at region ', region_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign data to appropriate arrays
      init_bounds(1,region_counter) = x_multiplier
      init_bounds(2,region_counter) = y_multiplier
      init_bounds(3,region_counter) = z_multiplier
      init_bounds(4,region_counter) = 0.0D0
      init_bounds(5,region_counter) = 0.0D0
      init_bounds(6,region_counter) = 0.0D0
!
! Write the initialization regions data to the output file
!
      WRITE(11,101) x_multiplier
      WRITE(11,102) y_multiplier
      WRITE(11,103) z_multiplier
      WRITE(11,*) ''
 101  FORMAT ('X multiplier = ',F9.7)
 102  FORMAT ('Y multiplier = ',F9.7)
 103  FORMAT ('Z multiplier = ',F9.7)

      END SUBROUTINE
