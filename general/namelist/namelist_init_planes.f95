      SUBROUTINE namelist_init_planes(region_shape,region_counter)
!
! Read initialization bounds info for 2D cases, uses the bounding_box
! info to make sure everything stays in the same plane
!
! Called by: namelist_init
! Calls: error_out
!
      USE nml_init
!      USE nml_bounding_box
      USE precise
!      USE output_formatting
use amrex_base_module
      IMPLICIT NONE
      INTEGER :: region_shape,region_counter,istat
      REAL(KIND=dp) :: min_axis_1, max_axis_1,min_axis_2,max_axis_2,&
        constant_axis_loc
! NAMELIST setup, only need the 4 definitions
      NAMELIST /INIT_REGION/ min_axis_1,max_axis_1,min_axis_2,max_axis_2
! Initialize the variables
      min_axis_1 = 0.0D0
      max_axis_1 = 1.0D0
      min_axis_2 = 0.0D0
      max_axis_2 = 1.0D0
! Read the namelist
      READ(2,INIT_REGION,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at init read type ',&
          region_shape, ' at region ', region_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign data to appropriate matrices
      init_bounds(1,region_counter) = min_axis_1
      init_bounds(2,region_counter) = max_axis_1
      init_bounds(3,region_counter) = min_axis_2
      init_bounds(4,region_counter) = max_axis_2
      init_bounds(5,region_counter) = amrex_problo(3)
      init_bounds(6,region_counter) = 0.0D0
!
! Write out the initialization region bounds
!
      IF (region_shape == -1) THEN
        WRITE(11,*) 'X-Y plane'
        WRITE(11,101) min_axis_1
        WRITE(11,102) max_axis_1
        WRITE(11,103) min_axis_2
        WRITE(11,104) min_axis_2
        WRITE(11,109) constant_axis_loc

      ELSE IF (region_shape == -2) THEN
        WRITE(11,*) 'X-Z plane'
        WRITE(11,101) min_axis_1
        WRITE(11,102) max_axis_1
        WRITE(11,105) min_axis_2
        WRITE(11,106) min_axis_2
        WRITE(11,108) constant_axis_loc

      ELSE IF (region_shape == -3) THEN
        WRITE(11,*) 'Y-Z plane'
        WRITE(11,103) min_axis_1
        WRITE(11,104) max_axis_1
        WRITE(11,105) min_axis_2
        WRITE(11,106) min_axis_2
        WRITE(11,107) constant_axis_loc

      END IF

      WRITE(11,*) ''
 101  FORMAT ('X minimum = ', F9.7)
 102  FORMAT ('X maximum = ', F9.7)
 103  FORMAT ('Y minimum = ', F9.7)
 104  FORMAT ('Y maximum = ', F9.7)
 105  FORMAT ('Z minimum = ', F9.7)
 106  FORMAT ('Z maximum = ', F9.7)
 107  FORMAT ('Constant X location = ', F9.7)
 108  FORMAT ('Constant Y location = ', F9.7)
 109  FORMAT ('Constant Z location = ', F9.7)

      END SUBROUTINE
