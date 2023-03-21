      SUBROUTINE namelist_bound_planes
!
! Bounding box exclusive to 2D runs.
! Need to include more controls on 2D vs. 3D down the road.
!
! Called by: namelist_bounding_box
! Calls: error_out
!
      USE nml_bounding_box
      USE precise
!      USE output_formatting
      IMPLICIT NONE
      REAL(KIND=dp) :: min_axis_1, max_axis_1,min_axis_2,max_axis_2
      REAL(KIND=dp) :: constant_axis_loc
      INTEGER :: istat
!
! Bounding box shape values for 2D
! -1 = x-y plane with user defined min_x,max_x,min_y,max_y
! -2 = x-z plane with user defined min_x,max_x,min_z,max_z
! -3 = y-z plane with user defined min_y,max_y,min_z,max_z
!
! Declare the NAMELIST
      NAMELIST /BOUND_SHAPE/ min_axis_1, max_axis_1,min_axis_2,&
        max_axis_2,constant_axis_loc
! Initialize the variables, axis order always x -> y -> z
      min_axis_1 = 0.0D0
      max_axis_1 = 1.0D0
      min_axis_2 = 0.0D0
      max_axis_2 = 1.0D0
      constant_axis_loc = 1.0D0
! Read in the data
      READ(2,BOUND_SHAPE,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at bounding box read type ',&
          bounding_method
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      bounds(1) = min_axis_1
      bounds(2) = max_axis_1
      bounds(3) = min_axis_2
      bounds(4) = max_axis_2
      bounds(5) = constant_axis_loc
!
! Write out the bounding box data to the output file.
!
      IF (bounding_method == -1) THEN
        WRITE(11,*) 'X-Y plane'
        WRITE(11,101) min_axis_1
        WRITE(11,102) max_axis_1
        WRITE(11,103) min_axis_2
        WRITE(11,104) min_axis_2
        WRITE(11,109) constant_axis_loc

      ELSE IF (bounding_method == -2) THEN
        WRITE(11,*) 'X-Z plane'
        WRITE(11,101) min_axis_1
        WRITE(11,102) max_axis_1
        WRITE(11,105) min_axis_2
        WRITE(11,106) min_axis_2
        WRITE(11,108) constant_axis_loc

      ELSE IF (bounding_method == -3) THEN
        WRITE(11,*) 'Y-Z plane'
        WRITE(11,103) min_axis_1
        WRITE(11,104) max_axis_1
        WRITE(11,105) min_axis_2
        WRITE(11,106) min_axis_2
        WRITE(11,107) constant_axis_loc

      END IF

      WRITE(11,*) ''
 101  FORMAT ('X minimum = ', F11.5)
 102  FORMAT ('X maximum = ', F11.5)
 103  FORMAT ('Y minimum = ', F11.5)
 104  FORMAT ('Y maximum = ', F11.5)
 105  FORMAT ('Z minimum = ', F11.5)
 106  FORMAT ('Z maximum = ', F11.5)
 107  FORMAT ('Constant X location = ', F11.5)
 108  FORMAT ('Constant Y location = ', F11.5)
 109  FORMAT ('Constant Z location = ', F11.5)
! 110  FORMAT (F9.7)
      END SUBROUTINE
