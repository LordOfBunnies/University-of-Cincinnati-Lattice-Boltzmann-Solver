      SUBROUTINE namelist_output_planes(shape_save,save_counter)
!
! Deals with getting the dimensions of the output plane
!
! Called by: namelist_output
! Calls:
!
      USE nml_output
      USE precise
      IMPLICIT NONE
      REAL(KIND=dp) :: min_axis_1, max_axis_1,min_axis_2,max_axis_2
      REAL(KIND=dp) :: constant_axis_loc
      INTEGER :: shape_save,save_counter,istat
!
! Save shape values
! -1 = x-y plane with user defined min_x,max_x,min_y,max_y
! -2 = x-z plane with user defined min_x,max_x,min_z,max_z
! -3 = y-z plane with user defined min_y,max_y,min_z,max_z
!
! Declare the NAMELIST
      NAMELIST /SAVE_REGION_PLANES/ min_axis_1, max_axis_1,min_axis_2,&
        max_axis_2,constant_axis_loc
! Initialize the data
      min_axis_1 = 0.0D0
      max_axis_1 = 1.0D0
      min_axis_2 = 0.0D0
      max_axis_2 = 1.0D0
      constant_axis_loc = 1.0D0
 ! Read the data
      READ(2,SAVE_REGION_PLANES,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at output read type ',&
          shape_save, ' at region ', save_counter
        WRITE(11,404) shape_save,save_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign data to the appropriate arrays
      save_regions(1,save_counter) = min_axis_1
      save_regions(2,save_counter) = max_axis_1
      save_regions(3,save_counter) = min_axis_2
      save_regions(4,save_counter) = max_axis_2
      save_regions(5,save_counter) = constant_axis_loc
      save_regions(6,save_counter) = 0.0D0
!
! Write out the initialization region bounds
!
      IF (shape_save == -1) THEN
        WRITE(11,*) 'X-Y plane'
        WRITE(11,101) min_axis_1
        WRITE(11,102) max_axis_1
        WRITE(11,103) min_axis_2
        WRITE(11,104) min_axis_2
        WRITE(11,109) constant_axis_loc

      ELSE IF (shape_save == -2) THEN
        WRITE(11,*) 'X-Z plane'
        WRITE(11,101) min_axis_1
        WRITE(11,102) max_axis_1
        WRITE(11,105) min_axis_2
        WRITE(11,106) min_axis_2
        WRITE(11,108) constant_axis_loc

      ELSE IF (shape_save == -3) THEN
        WRITE(11) 'Y-Z plane'
        WRITE(11,103) min_axis_1
        WRITE(11,104) max_axis_1
        WRITE(11,105) min_axis_2
        WRITE(11,106) min_axis_2
        WRITE(11,107) constant_axis_loc

      END IF

      WRITE(11,*) ''
 101  FORMAT ('X minimum = ', F11.7)
 102  FORMAT ('X maximum = ', F11.7)
 103  FORMAT ('Y minimum = ', F11.7)
 104  FORMAT ('Y maximum = ', F11.7)
 105  FORMAT ('Z minimum = ', F11.7)
 106  FORMAT ('Z maximum = ', F11.7)
 107  FORMAT ('Constant X location = ', F11.7)
 108  FORMAT ('Constant Y location = ', F11.7)
 109  FORMAT ('Constant Z location = ', F11.7)
 404  FORMAT ('Error reading namelist at output read type ',&
          I4, ' at region ', I4)
      END SUBROUTINE
