      SUBROUTINE namelist_init
!
! Namelist reader to set up initialization parameters in case those should
! be different from the basic freestream values
!
! Called by: namelist_global
! Calls: namelist_init_region_1,namelist_init_region_2,
!   namelist_init_region_3,namelist_init_region_11_12,
!   namelist_init_region_20s,namelist_init_planes,error_out
!
      USE nml_init
      USE precise
      USE constants
      USE freestream_values
      IMPLICIT NONE
      INTEGER :: i,region_shape,priority,istat
      REAL(KIND=dp) :: init_velocity,init_pressure, init_temperature
      REAL(KIND=dp) :: x_vel_comp,y_vel_comp,z_vel_comp
! Bounding Methods:
! 1 = box with x,y,z distance multipliers based on max geometry in 
!       that direction
! 2 = box with user defined max and min in x, y, and z
! 3 = box with user defined center and multiplier in x, y, and z
!       based on geometry dimensions
!
! 
! 11 = Sphere with user defined center and radius based on multiplier
!        from geometry
! 12 = Sphere with user defined center and radius
!
!
! 21 = Cylinder with user defined x-direction, starting point, radius, and
!       distance
! 22 = Cylinder with user defined y-direction, starting point, radius, and
!       distance
! 23 = Cylinder with user defined z-direction, starting point, radius, and
!       distance
! 24 = Cylinder with user defined starting point and x-direction, with 
!        radius and distance based on multiplier of geometry
! 25 = Cylinder with user defined starting point and y-direction, with 
!        radius and distance based on multiplier of geometry
! 26 = Cylinder with user defined starting point and z-direction, with 
!        radius and distance based on multiplier of geometry
!
! Region priority defines which data to use if two regions overlap, if 
!  priorities are equal for a node, the values will be averaged
!
! x,y,z components are only to tell which direction the flow will be going
!  in the initialized volume, they do not affect magnitude
!
! If nothing else is put for initialization values, freestream values will
!  be used.
!
      NAMELIST /INIT_VALUES/ region_shape,priority,init_temperature,&
        init_pressure,init_velocity,x_vel_comp,y_vel_comp,z_vel_comp
!
! Allocate all the necessary space for the init regions
!
! init_info(1,*) = region shape
! init_info(2,*) = region priority
! init_data(1,*) = temperature
! init_data(2,*) = pressure
! init_data(3,*) = velocity magnitude
! init_data(4,*) = x-component of velocity
! init_data(5,*) = y-component of velocity
! init_data(6,*) = z-component of velocity
! init_data 4-6 will be used to make a unit vector and multiplied by
! init_data 3 to get the components of velocity
!
      ALLOCATE(init_info(2,num_init_regions))
      ALLOCATE(init_data(6,num_init_regions))
      ALLOCATE(init_bounds(6,num_init_regions))
!
! Loop to pull the data from the NAMELIST
!
      DO i = 1,num_init_regions
        region_shape = 1
        priority = 1
        init_temperature = temperature
        init_pressure = pressure
        init_velocity = v_inf
        x_vel_comp = 1.0D0
        y_vel_comp = 0.0D0
        z_vel_comp = 0.0D0

! Read the data in
        READ(2,INIT_VALUES,IOSTAT=istat)
        IF (istat .NE. 0) THEN
          WRITE(*,*) 'Error reading namelist at initialization',&
            ' values.', i
          WRITE(*,*) istat
            CALL error_out
        END IF
! Apply the data to the appropriate array
        init_info(1,i) = region_shape
        init_info(2,i) = priority
        init_data(1,i) = init_temperature
        init_data(2,i) = init_pressure
        init_data(3,i) = init_velocity
        init_data(4,i) = x_vel_comp
        init_data(5,i) = y_vel_comp
        init_data(6,i) = z_vel_comp


      END DO
!
! Tells it where to go based on the region_shape to get its data
! Same basic shapes as in bounding_box
!
      DO i = 1,num_init_regions
!
!
!
        IF (init_info(1,i) == 1) THEN
!
          WRITE(11,207) i
          CALL namelist_init_region_1(init_info(1,i),i)
!
          WRITE(11,101) i,init_info(2,i)
          WRITE(11,102) i,init_data(1,i)
          WRITE(11,103) i,init_data(2,i)
          WRITE(11,104) i,init_data(3,i)
          WRITE(11,105) i,init_data(4,i)
          WRITE(11,106) i,init_data(5,i)
          WRITE(11,107) i,init_data(6,i)
          WRITE(11,*) ''
!
!
! 
        ELSE IF (init_info(1,i) == 2) THEN
!
          WRITE(11,206) i
          CALL namelist_init_region_2(init_info(1,i),i)
!
          WRITE(11,101) i,init_info(2,i)
          WRITE(11,102) i,init_data(1,i)
          WRITE(11,103) i,init_data(2,i)
          WRITE(11,104) i,init_data(3,i)
          WRITE(11,105) i,init_data(4,i)
          WRITE(11,106) i,init_data(5,i)
          WRITE(11,107) i,init_data(6,i)
          WRITE(11,*) ''
!
!
!           
        ELSE IF (init_info(1,i) == 3) THEN
!
          WRITE(11,205) i
          CALL namelist_init_region_3(init_info(1,i),i)
!
          WRITE(11,101) i,init_info(2,i)
          WRITE(11,102) i,init_data(1,i)
          WRITE(11,103) i,init_data(2,i)
          WRITE(11,104) i,init_data(3,i)
          WRITE(11,105) i,init_data(4,i)
          WRITE(11,106) i,init_data(5,i)
          WRITE(11,107) i,init_data(6,i)
          WRITE(11,*) ''
!
!
!
        ELSE IF (init_info(1,i) == 11 .OR. init_info(1,i) == 12) THEN
!
          IF (init_info(1,i) == 11) THEN
            WRITE(11,203) i
          ELSE
            WRITE(11,204) i
          END IF
          CALL namelist_init_region_11_12(init_info(1,i),i)
!
          WRITE(11,101) i,init_info(2,i)
          WRITE(11,102) i,init_data(1,i)
          WRITE(11,103) i,init_data(2,i)
          WRITE(11,104) i,init_data(3,i)
          WRITE(11,105) i,init_data(4,i)
          WRITE(11,106) i,init_data(5,i)
          WRITE(11,107) i,init_data(6,i)
          WRITE(11,*) ''
!
!
!
        ELSE IF (init_info(1,i) >= 21 .AND. init_info(1,i) <= 26) THEN
!
          IF (init_info(1,i) <=23) THEN
           WRITE(11,201) i
          ELSE
            WRITE(11,202) i
          END IF
          CALL namelist_init_region_20s(init_info(1,i),i)
!
          WRITE(11,101) i,init_info(2,i)
          WRITE(11,102) i,init_data(1,i)
          WRITE(11,103) i,init_data(2,i)
          WRITE(11,104) i,init_data(3,i)
          WRITE(11,105) i,init_data(4,i)
          WRITE(11,106) i,init_data(5,i)
          WRITE(11,107) i,init_data(6,i)
          WRITE(11,*) ''
!
! Call for 2D stuff only
!
        ELSE IF (init_info(1,i) <= -1 .AND. init_info(1,i) >= -3) THEN
!
          WRITE(11,*) 'Initialization region is a plane, this can',&
            ' only happen in 2D runs.'
          CALL namelist_init_planes(init_info(1,i),i)
!
          WRITE(11,101) i,init_info(2,i)
          WRITE(11,102) i,init_data(1,i)
          WRITE(11,103) i,init_data(2,i)
          WRITE(11,104) i,init_data(3,i)
          WRITE(11,105) i,init_data(4,i)
          WRITE(11,106) i,init_data(5,i)
          WRITE(11,107) i,init_data(6,i)
          WRITE(11,*) ''
!
        ELSE
          WRITE(*,*) 'Invalid initialization region chosen, exiting.'
          CALL error_out
        END IF
      END DO


 101  FORMAT ('Initialization region ',I8,', priority = ', I8)
 102  FORMAT ('Initialization region ',I8,', initial temperature = ',F11.7)
 103  FORMAT ('Initialization region ',I8,', initial pressure = ',F15.7)
 104  FORMAT ('Initialization region ',I8,', initial velocity = ',F11.7)
 105  FORMAT ('Initialization region ',I8,', x-velocity component = ',F9.7)
 106  FORMAT ('Initialization region ',I8,', y-velocity component = ',F9.7)
 107  FORMAT ('Initialization region ',I8,', z-velocity component = ',F9.7)


 201  FORMAT ('Initialization region ',I8,' is a cylinder with',&
            ' user defined center',&
              ', radius, and distance.')
 202  FORMAT ('Initialization region ',I8,' is a cylinder',&
            ' with user defined center and',&
            ' radius with distance based on multiplier from geometry.')
 203  FORMAT ('Initialization region ',I8,' is a sphere with',&
            ' user defined center',&
              ' and radius based on multiplier from geometry.')
 204  FORMAT ('Initialization region ',I8,' is a sphere with',&
            ' user defined center and radius.')
 205  FORMAT ('Initialization region ',I8,' is a rectangular',&
          ' prism with user defined x,y, and z center and multipliers.')
 206  FORMAT ('Initialization region ',I8,' is a rectangular',&
            ' prism with user defined min and max x,y, and z.')
 207  FORMAT ('Initialization region ',I8,' is a rectangular',&
            ' prism with x,y, and z multipliers.')

      END SUBROUTINE
