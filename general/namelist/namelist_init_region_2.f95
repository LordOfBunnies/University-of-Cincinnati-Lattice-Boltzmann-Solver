      SUBROUTINE namelist_init_region_2(region_shape,region_counter)
!
! Call to set up rectangular init region size
!
! Called by: namelist_init
! Calls: error_out
!
      USE nml_init
      USE precise
      IMPLICIT NONE
      INTEGER :: region_shape,istat,region_counter
      REAL(KIND=dp) :: x_lower_lim,y_lower_lim,z_lower_lim,&
            x_upper_lim,y_upper_lim,z_upper_lim
! Declare the NAMELIST
      NAMELIST /INIT_REGION/ x_lower_lim,y_lower_lim,z_lower_lim,&
            x_upper_lim,y_upper_lim,z_upper_lim
! Initialize the variables
      x_lower_lim = -1.0D0
      y_lower_lim = -1.0D0
      z_lower_lim = -1.0D0
      x_upper_lim = 1.0D0
      y_upper_lim = 1.0D0
      z_upper_lim = 1.0D0
! Read the data
      READ(2,INIT_REGION,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at init read type ',&
          region_shape, ' at region ', region_counter
        WRITE(*,*) istat
        CALL error_out
      END IF
! Assign the data to appropriate arrays
      init_bounds(1,region_counter) = x_lower_lim
      init_bounds(2,region_counter) = y_lower_lim
      init_bounds(3,region_counter) = z_lower_lim
      init_bounds(4,region_counter) = x_upper_lim
      init_bounds(5,region_counter) = y_upper_lim
      init_bounds(6,region_counter) = z_upper_lim
!
! Write the initialization region data to the output file
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
