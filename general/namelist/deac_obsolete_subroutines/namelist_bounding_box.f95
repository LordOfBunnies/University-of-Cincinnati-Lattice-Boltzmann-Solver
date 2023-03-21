      SUBROUTINE namelist_bounding_box
!
! Determines the limits of the bounding box for the grid points, will 
! vary greatly depending on method and there will need to be a second 
! file in case of certain rare instances
!
! Called by: namelist_global
! Calls: namelist_bound_1,namelist_bound_2,namelist_bound_3,
!   namelist_bound_11_12,namelist_bound_20s,namelist_bound_planes,
!   error_out
!
      USE nml_bounding_box
      USE constants
      IMPLICIT NONE
      INTEGER :: istat
!
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
! 21 = Cylinder in x-direction with user defined starting point, radius, and
!       distance
! 22 = Cylinder in y-direction with user defined starting point, radius, and
!       distance
! 23 = Cylinder in z-direction with user defined starting point, radius, and
!       distance
! 24 = Cylinder with user defined starting point, x-direction,
!        radius, and with distance based on multiplier of geometry
! 25 = Cylinder with user defined starting point, y-direction,
!        radius, and with distance based on multiplier of geometry
! 26 = Cylinder with user defined starting point, z-direction,
!        radius, and with distance based on multiplier of geometry
!
!
! -1 = x-y plane with user defined min_x,max_x,min_y,max_y
! -2 = x-z plane with user defined min_x,max_x,min_z,max_z
! -3 = y-z plane with user defined min_y,max_y,min_z,max_z
!
! Declare the NAMELIST
      NAMELIST /BOUNDING_BOX/ bounding_method
      bounding_method = 1
! Read the bounding method in
      READ(2,BOUNDING_BOX,IOSTAT=istat)
      IF (istat .NE. 0) THEN
        WRITE(*,*) 'Error reading namelist at bounding box',&
            ' region read.'
        WRITE(11,*) 'Error reading namelist at bounding box',&
            ' region read.'
        WRITE(*,*) istat
        CALL error_out
      END IF
!
! Internal flow not implemented yet
! Allocates the array and sends the program to the appropriate
!
      IF (external_flow) THEN
!
        IF (bounding_method == 1) THEN
          ALLOCATE(bounds(3))
          WRITE(11,*) 'Bounding box is a rectangular prism with x,y,',&
            ' and z multipliers.'
          CALL namelist_bound_1
!
! 
        ELSE IF (bounding_method == 2) THEN
          ALLOCATE(bounds(6))
          WRITE(11,*) 'Bounding box is a rectangular prism with user',&
            ' defined mix and max x,y, and z.'
          CALL namelist_bound_2
!
!           
        ELSE IF (bounding_method == 3) THEN
          ALLOCATE(bounds(6))
          WRITE(11,*) 'Bounding box is a rectangular prism with user',&
            ' defined center x,y, and z multipliers.'
          CALL namelist_bound_3
!
!
        ELSE IF (bounding_method == 11 .OR. bounding_method == 12) THEN
          ALLOCATE(bounds(4))
          IF (bounding_method == 11) THEN
            WRITE(11,*) 'Bounding box is a sphere with user defined center',&
              ' and radius based on multiplier from geometry.'
          ELSE
            WRITE(11,*) 'Bounding box is a sphere with user defined center',&
            ' and radius.'
          END IF
          CALL namelist_bound_11_12
!
!
        ELSE IF (bounding_method >= 21 .AND. bounding_method <= 26) THEN
          ALLOCATE(bounds(5))
          IF (bounding_method <=23) THEN
            WRITE(11,*) 'Bounding box is a cylinder with user defined center',&
              ', radius, and distance.'
          ELSE
            WRITE(11,*) 'Bounding box is a cylinder with user defined center',&
              ' and radius with distance based on multiplier from geometry.'
          END IF
          CALL namelist_bound_20s
!
!
        ELSE IF (bounding_method >= -3 .AND. bounding_method <= -1) THEN
          ALLOCATE(bounds(5))
          WRITE(11,*) 'Bounding box is a plane, therefore this is a 2D run.'
          CALL namelist_bound_planes
!
!
        ELSE
          WRITE(*,*) 'Invalid bounding method chosen, exiting.'
          CALL error_out
        END IF
!
! Much more limited for internal flows
!
      ELSE
        WRITE(*,*) 'Internal flow not yet implemented.'
        CALL error_out
      END IF

      END SUBROUTINE
