subroutine namelist_output
!
! Namelist reader for dealing with outputs
!
! Called by: namelist_general
! Calls: namelist_output_1,namelist_output_2,namelist_output_3,
!   namelist_output_11_12,namelist_output_20s,namelist_output_planes,
!   error_out
!
use nml_output
use constants
use precise
use amrex_base_module
use amrex_amr_module
IMPLICIT NONE
REAL(KIND=dp) :: save_start_time,save_stop_time,restart_save_time,save_freq
INTEGER :: i,istat,min_save_level,max_save_level
LOGICAL :: density_save,u_vel_save,v_vel_save,w_vel_save,&
  pressure_save,temperature_save,alpha_save,lagrange_save
logical :: chi_save,zeta_x_save,zeta_y_save,zeta_z_save,pi_xx_save,&
  pi_yy_save,pi_zz_save,pi_xy_save,pi_xz_save,pi_yz_save,lambda_x_save,&
  lambda_y_save,lambda_z_save
! 
! Declare all the NAMELISTs to be used
!
namelist /SHAPE_OF_SAVE/ min_save_level,max_save_level
namelist /RESTART_INFO/ restart_save_time,restarting,write_restart_at_end
namelist /SAVE_TIMING/ save_start_time,save_stop_time,save_freq
namelist /SAVED_DATA/ density_save,u_vel_save,v_vel_save,&
  w_vel_save,pressure_save,temperature_save,alpha_save,lagrange_save
namelist /SAVED_LMS/ chi_save,zeta_x_save,zeta_y_save,zeta_z_save,pi_xx_save,&
  pi_yy_save,pi_zz_save,pi_xy_save,pi_xz_save,pi_yz_save,lambda_x_save,&
  lambda_y_save,lambda_z_save
!
! save_shapes defines the min and max level of output for a given output
! save_times is (start time,*) (stop_time,*)
! save_regions is defined similar to bounding_box info
! restart_saves is at what times to output a restartable files
! save_frequency says to output every X seconds.
!
allocate(save_shapes(2,num_saves))
allocate(save_times(2,num_saves))
!allocate(save_regions(6,num_saves))
allocate(restart_saves(num_restart_saves))
allocate(save_frequency(num_saves))
allocate(save_data(8,num_saves))
if (dimensions == 2) then
  allocate(lm_saves(8,num_saves))
else if (dimensions == 3) then
  allocate(lm_saves(13,num_saves))
end if
!
! Read the region shape so things can be assigned later
!
IF (num_saves >= 1) THEN
  DO i = 1, num_saves
    min_save_level = 0
    max_save_level = 1
    READ(2,SHAPE_OF_SAVE,IOSTAT=istat)
    IF (istat .NE. 0) THEN
      WRITE(*,*) 'Error reading namelist at output',&
          ' region ',i
      WRITE(11,*) 'Error reading namelist at output',&
          ' region ',i
      WRITE(*,*) istat
      CALL error_out
    END IF
    save_shapes(1,i) = min_save_level
    save_shapes(2,i) = max_save_level

    WRITE(11,*) 'Save number ',i,' will save all levels from ',save_shapes(1,i),' to ',&
      save_shapes(2,i)
!    IF (save_shapes(i) == 1) THEN
!      WRITE(11,*)'Box with x,y,z distance multipliers based on max geometry in',&
!        'that direction'
!!
!    ELSE IF (save_shapes(i) == 2) THEN
!      WRITE(11,*) 'Box with user defined max and min in x, y, and z'
!!
!    ELSE IF (save_shapes(i) == 3) THEN
!      WRITE(11,*) 'Box with user defined center and multiplier in x, y, and z',&
!        'based on geometry dimensions'
!!
!    ELSE IF (save_shapes(i) == 11) THEN
!      WRITE(11,*) 'Sphere with user defined center and radius based on multiplier',&
!        ' from geometry'
!!
!    ELSE IF (save_shapes(i) == 12) THEN
!      WRITE(11,*) 'Sphere with user defined center and radius'
!!
!    ELSE IF (save_shapes(i) >= 21 .AND. save_shapes(i) <= 23) THEN
!      WRITE(11,*) 'Cylinder in a direction with user defined starting point, radius, and',&
!        ' distance'
!!
!    ELSE IF (save_shapes(i) >= 24 .AND. save_shapes(i) <= 26) THEN
!      WRITE(11,*) 'Cylinder with user defined starting point, a direction,',&
!        ' radius, and with distance based on multiplier of geometry'
!!
!    ELSE IF (save_shapes(i) <= -1 .AND. save_shapes(i) >= -3) THEN
!      WRITE(11,*) 'Plane with user defined minimum and maximum axes.'
!
    if (save_shapes(2,i) < save_shapes(1,i)) then
      WRITE(*,*) 'Invalid output region chosen max is lower than min, exiting.'
      CALL error_out
!    else if (save_shapes(2,i) > amrex_max_level) then
!      WRITE(11,*) 'Upper output bound level finer than max level selected,',&
!        ' setting values to max level and continuing.'
!      WRITE(*,*) 'Upper output bound level finer than max level selected',&
!        ' setting values to max level and continuing.'
!      save_shapes(2,i) = amrex_max_level
    else if (save_shapes(1,i) < 0) then
       WRITE(11,*) 'Lower output bound level less than 0,',&
        ' setting min level to 0 and continuing.'
    END IF
!
! Initialize then read in the the save timing data
!
    save_start_time = 0.0D0
    save_stop_time = 1.0D0
    save_freq = 0.10D0
!
    READ(2,SAVE_TIMING,IOSTAT=istat)
    IF (istat .NE. 0) THEN
      WRITE(*,*) 'Error reading namelist at save',&
          ' timing ',i
      WRITE(11,*) 'Error reading namelist at save',&
          ' timing ',i
      WRITE(*,*) istat
      CALL error_out
    END IF
    save_times(1,i) = save_start_time
    save_times(2,i) = save_stop_time
    IF (save_freq < 0.0D0) THEN
      WRITE(11,*) 'Save frequency less than 0s, setting to 0.1s.'
      save_freq = 0.1D0
    END IF
    save_frequency(i) = save_freq

    WRITE(11,200) i
    WRITE(11,101) save_times(1,i),save_times(2,i)
    WRITE(11,201) save_frequency(i)

    density_save = .FALSE.
    u_vel_save = .FALSE.
    v_vel_save = .FALSE.
    w_vel_save = .FALSE.
    pressure_save = .FALSE.
    temperature_save = .FALSE.
    alpha_save = .FALSE.
    lagrange_save = .FALSE.
    READ(2,SAVED_DATA,IOSTAT=istat)
    IF (istat .NE. 0) THEN
      WRITE(*,*) 'Error reading namelist at output info.'
      WRITE(11,*) 'Error reading namelist at output info.'
      WRITE(*,*) istat
      CALL error_out
    END IF
    save_data(1,i) = density_save
    save_data(2,i) = u_vel_save
    save_data(3,i) = v_vel_save
    save_data(4,i) = w_vel_save
    save_data(5,i) = pressure_save
    save_data(6,i) = temperature_save
    save_data(7,i) = alpha_save
    save_data(8,i) = lagrange_save

!          WRITE(11,200) 'For save number ', i
!          WRITE(11,*) ''
    WRITE(11,405) save_data(1,i),save_data(2,i), &
      save_data(3,i),save_data(4,i),save_data(5,i),&
      save_data(6,i),save_data(7,i),save_data(8,i)
!    WRITE(*,405) save_data(1,i),save_data(2,i), &
!      save_data(3,i),save_data(4,i),save_data(5,i),&
!      save_data(6,i),save_data(7,i)

    chi_save = .false.
    zeta_x_save = .false.
    zeta_y_save = .false.
    zeta_z_save = .false.
    pi_xx_save = .false.
    pi_yy_save = .false.
    pi_zz_save = .false.
    pi_xy_save = .false.
    pi_xz_save = .false.
    pi_yz_save = .false.
    lambda_x_save = .false.
    lambda_y_save = .false.
    lambda_z_save = .false.

    if (save_data(7,i)) then
    READ(2,SAVED_LMS,IOSTAT=istat)
    IF (istat .NE. 0) THEN
      WRITE(*,*) 'Error reading namelist at Lagrange multiplier info.'
      WRITE(11,*) 'Error reading namelist at Lagrange multiplier info.'
      WRITE(*,*) istat
      CALL error_out
    END IF


    if (dimensions == 2) then
      lm_saves(1,i) = chi_save
      lm_saves(2,i) = zeta_x_save
      lm_saves(3,i) = zeta_y_save
      lm_saves(4,i) = pi_xx_save
      lm_saves(5,i) = pi_yy_save
      lm_saves(6,i) = pi_xy_save
      lm_saves(7,i) = lambda_x_save
      lm_saves(8,i) = lambda_y_save

      write(11,505) lm_saves(1,i),lm_saves(2,i),lm_saves(3,i),&
        lm_saves(4,i),lm_saves(5,i),lm_saves(6,i),lm_saves(7,i),&
        lm_saves(8,i)

    else if (dimensions == 3) then
      lm_saves(1,i) = chi_save
      lm_saves(2,i) = zeta_x_save
      lm_saves(3,i) = zeta_y_save
      lm_saves(4,i) = zeta_y_save
      lm_saves(5,i) = pi_xx_save
      lm_saves(6,i) = pi_yy_save
      lm_saves(7,i) = pi_zz_save
      lm_saves(8,i) = pi_xy_save
      lm_saves(9,i) = pi_xz_save
      lm_saves(10,i) = pi_yz_save
      lm_saves(11,i) = lambda_x_save
      lm_saves(12,i) = lambda_y_save
      lm_saves(13,i) = lambda_z_save

      write(11,506) lm_saves(1,i),lm_saves(2,i),lm_saves(3,i),&
        lm_saves(4,i),lm_saves(5,i),lm_saves(6,i),lm_saves(7,i),&
        lm_saves(8,i),lm_saves(9,i),lm_saves(10,i),lm_saves(11,i),&
        lm_saves(12,i),lm_saves(13,i)

    end if
    end if
!          WRITE(11,*) ''
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
! -1 = x-y plane with user defined min_x,max_x,min_y,max_y
! -2 = x-z plane with user defined min_x,max_x,min_z,max_z
! -3 = y-z plane with user defined min_y,max_y,min_z,max_z
!
! Send the program to an appropriate subroutine to read in bounding data
!
!      IF (num_saves >= 1) THEN
!        DO i = 1,num_saves

!    IF (save_shapes(i) == 1) THEN
!      CALL namelist_output_1(save_shapes(i),i)
!!
!    ELSE IF (save_shapes(i) == 2) THEN
!      CALL namelist_output_2(save_shapes(i),i)
!!
!    ELSE IF (save_shapes(i) == 3) THEN
!      CALL namelist_output_3(save_shapes(i),i)
!!
!    ELSE IF (save_shapes(i) == 11) THEN
!      CALL namelist_output_11_12(save_shapes(i),i)
!!
!    ELSE IF (save_shapes(i) == 12) THEN
!      CALL namelist_output_11_12(save_shapes(i),i)
!!
!    ELSE IF (save_shapes(i) >= 21 .AND. save_shapes(i) <= 26) THEN
!      CALL namelist_output_20s(save_shapes(i),i)
!!
!    ELSE IF (save_shapes(i) <= -1 .AND. save_shapes(i) >= -3) THEN
!      CALL namelist_output_planes(save_shapes(i),i)
!!
!    ELSE
!      WRITE(*,*) 'Invalid output region chosen, exiting.'
!      CALL error_out
!    END IF

  END DO
END IF
!
!
! Next we figure out what we are saving
!
!      IF (num_saves >= 1) THEN
!        DO i = 1,num_saves
!          restart_save_time = 999.0


!        END DO
!      END IF
!
! Restart saves aren't necessary, but you can have them written out at
! specific times.  If the time declared doesn't fall within the total
! time, the only restart file will be the one at the end.
!
IF (num_restart_saves>=1) THEN

  DO i = 1,num_restart_saves
    restart_save_time = 999.0
    restarting = .FALSE.
    write_restart_at_end = .FALSE.
    READ(2,RESTART_INFO,IOSTAT=istat)
    IF (istat .NE. 0) THEN
      write(*,*) 'Error reading namelist at restart info.'
      write(*,*) istat
      CALL error_out
    END IF
    restart_saves(i) = restart_save_time
    WRITE(11,151) i,restart_saves(i)

    IF (restarting) THEN
      WRITE(11,*) 'This is a restarted run.'
    ELSE
      WRITE(11,*) 'This is a new run.'
    END IF

    IF (write_restart_at_end)  THEN
      WRITE(11,*) 'This will write a restart file at the end ',&
        'of the run.'
    END IF

  END DO
END IF

!        END DO
!      END IF
  
 101  format ('Saves will start at time ',F8.4, &
            ' and end at time ',F8.4)
 151  format ('Restart save ', I8, ' to occuar at time ', &
        F9.3)
 200  format ('For save number ', I4)
 201  format ('It saves every ',F9.7,'  timesteps.')
 251  format (I8, A32)
 400  format (L5)
 405  format ('Density saved? ', L3/, &
            'U velocity saved? ', L3/,'V velocity saved? ',&
            L3/, 'W velocity saved? ',L3/,&
            'Pressure saved? ',L3/,'Temperature saved? ',&
            L3/,'Lagrange Multipliers saved?',L3/)

 505  format ('Chi saved? ',L3/,'Zeta x saved? ',L3/,'Zeta y saved? ', &
   L3/,'Pi xx saved?',L3/,'Pi yy saved?',L3/,'Pi xy saved?',L3/,'Lambda x saved?', &
   L3/,'Lambda y saved?',L3/)

 506  format ('Chi saved? ',L3/,'Zeta x saved?',L3/,'Zeta y saved?',L3/,&
   'Zeta z saved?',L3/,'Pi xx saved?',L3/,'Pi yy saved?',L3/,'Pi zz saved?',L3/, &
   'Pi xy saved?',L3/,'Pi xz saved?',L3/,'Pi yz saved?',L3/,'Lambda x saved?', &
   L3/,'Lambda y saved?',L3/,'Lambda z saved?',L3/)
end subroutine
