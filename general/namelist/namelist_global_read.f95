SUBROUTINE namelist_global_read(characteristic_length,mic_file)
!
! This reads the information for the global variables because it's
! being a pain in the ass in the main directory
!
! Called by: namelist_global
! Calls: error_out
!
USE nml_geometry
USE nml_init
USE nml_inlet_outlet
USE grid_data
USE startup
USE nml_output
use linkwise
USE constants
USE freestream_values
USE timing
USE precise
USE nml_inputs
use amr_info_holder, only: amr,multigrid
USE ggm_stuff
IMPLICIT NONE
REAL(KIND=dp) :: characteristic_length
INTEGER :: istat
CHARACTER(len=32) :: mic_file
!
! Declare the global namelist to be read in.
!
NAMELIST /GLOBAL/ primary_flow_dir,num_init_regions,&
  num_geom_files,num_restart_saves,num_saves,num_mics,num_cgns_files,&
  dimensions,characteristic_length,cells_per_length,&
  prandtl,total_time,coarse_time,startup_timesteps,gama,&
  temperature,pressure,v_inf,viscous,turbulent,&
  bomb_file,body_forces,restart_start,cgns_start,&
  mixed_start,external_flow,gas,turb_type,mic_file,amr,multigrid,use_ggm,&
  use_shifted,gradient_refine_only

! Integer initializations
primary_flow_dir = 1 ! 1=+x, 2=+y, 3=+z, 4=-x, 5=-y, 6=-z
num_init_regions = 0
num_geom_files = 1
num_restart_saves = 0
num_saves = 1
num_mics = 0
num_cgns_files = 0
dimensions = 3 !Can be either 2 or 3
cells_per_length = 20
startup_timesteps = 25
! Real initializations
prandtl = 0.72D0
total_time = 0.0D0
coarse_time = 0.0D0
gama = 1.4D0
temperature = 273.0D0
pressure = 101300.0D0
v_inf = 0.0D0
!      tau = .96D0
characteristic_length = 1.0D0
! Logical Initializations
viscous = .TRUE.
turbulent = .FALSE.
amr = .TRUE.
external_flow = .TRUE.
bomb_file = .FALSE.
use_shifted = .FALSE.
body_forces = .FALSE.
restart_start = .FALSE.
cgns_start = .FALSE.
mixed_start = .FALSE.
multigrid = .TRUE.
use_ggm = .FALSE.
gradient_refine_only = .FALSE.
! String initializations
gas = 'AIR'
turb_type = 'KE'
mic_file = 'mic_file.dat'

READ(2,NML=GLOBAL,IOSTAT=istat)
IF (istat .NE. 0) THEN
  WRITE(*,*) 'Error reading namelist file, exiting'
  WRITE(11,*) 'Error reading namelist file, exiting'
  WRITE(*,*) istat
  CALL error_out
END IF
!
! Begin writing the output file with all the data regurgitated
!
!      WRITE(11,201) min_vr
!      WRITE(11,202) max_vr
IF (primary_flow_dir == 1) THEN
  WRITE(11,*) 'Primary flow direction is +x direction.'
ELSE IF (primary_flow_dir == 2) THEN
  WRITE(11,*) 'Primary flow direction is +y direction.'
ELSE IF (primary_flow_dir == 3) THEN
  WRITE(11,*) 'Primary flow direction is +z direction.'
ELSE IF (primary_flow_dir == 4) THEN
  WRITE(11,*) 'Primary flow direction is -x direction.'
ELSE IF (primary_flow_dir == 5) THEN
  WRITE(11,*) 'Primary flow direction is -y direction.'
ELSE IF (primary_flow_dir == 6) THEN
  WRITE(11,*) 'Primary flow direction is -z direction.'
END IF

 201  FORMAT ('Minimum VR = ', I8)
 202  FORMAT ('Maximum VR = ', I8)
!
! Begin writing out the integer quantities to the output file.
!
WRITE(11,*) ''
WRITE(11,*) 'Number of saves and region inputs'
WRITE(11,101) num_init_regions
WRITE(11,102) num_geom_files
WRITE(11,103) num_restart_saves
WRITE(11,122) num_cgns_files
WRITE(11,104) num_saves
WRITE(11,105) num_mics
WRITE(11,106) dimensions
WRITE(11,107) cells_per_length
WRITE(11,*) ''
WRITE(11,*) 'Constants and ambient flow variables'
WRITE(11,108) prandtl
WRITE(11,109) total_time
WRITE(11,110) coarse_time
WRITE(11,111) gama
WRITE(11,112) temperature
WRITE(11,113) pressure
WRITE(11,114) v_inf
WRITE(11,116) characteristic_length
IF (viscous) THEN
  WRITE(11,*) 'Simulation is viscous'
ELSE
  WRITE(11,*) 'Simulation is inviscid'
END IF
      !
!      IF (vles) THEN
!        WRITE(11,*) 'Simulation will include VLES method, this is',&
!          ' currently unused.'
!      ELSE
!        WRITE(11,*) 'Simulation will not include VLES method, ',&
!          'this is currently unused.'
!      END IF
      !
!      IF (energy) THEN
!        WRITE(11,*) 'Simulation will include calculations of energy.'
!      ELSE
!        WRITE(11,*) 'Simulation will not include calculations of energy.'
!      END IF
      !
IF (turbulent) THEN
  WRITE(11,*) 'Simulation will be NOT turbulent, this is',&
    ' currently unused.'
  WRITE(11,117) turb_type
ELSE
  WRITE(11,*) 'Simulation will not be turbulent, ',&
    'this is currently unused.'
END IF
!
IF (external_flow) THEN
  WRITE(11,*) 'Simulation will be an external flow.'
ELSE
  WRITE(11,*) 'Simulation will be an internal flow.'
END IF
!
IF (bomb_file) THEN
  WRITE(11,*) 'Simulation will write a bomb file if it fails,',&
    ' this is currently unused.'
ELSE
  WRITE(11,*) 'Simulation will not write a bomb file if it fails,',&
    ' this is currently unused.'
END IF
!
IF (body_forces) THEN
  WRITE(11,*) 'Simulation will include body forces,',&
    'this is currently unused.'
ELSE
  WRITE(11,*) 'Simulation will not include body forces,',&
    ' this is currently unused.'
END IF
!
write(11,*) 'This run will be shifted',use_shifted
!
IF (restart_start) THEN
  WRITE(11,*) 'This will be a restarted run.'
ELSE
  IF (cgns_start) THEN
    IF (mixed_start) THEN
      WRITE(11,*) 'This will be a mixed start, using a CGNS',&
        ' file and imported geometry together.'
    ELSE
      WRITE(11,*) 'This will be a CGNS started run.'
    END IF
  ELSE
    WRITE(11,*) 'This will be a standard run of UCLBS.'
  END IF
END IF
WRITE(11,118) gas
WRITE(11,*) ''

!      dx = characteristic_length/DBLE(cells_per_length)
!      WRITE(11,*) ''
!      WRITE(11,208) dx
!      WRITE(11,*) ''


 101  FORMAT ('Number of initialization regions = ', I8)
 102  FORMAT ('Number of geometry files = ',I8)
 103  FORMAT ('Number of full restart saves = ', I8)
 104  FORMAT ('Number of outputs = ', I8)
 105  FORMAT ('Number of micorphones = ', I8)
 106  FORMAT ('Dimensions of the solutions = ', I8)
 107  FORMAT ('Number of cells along the characteristic length = ',I8)
 108  FORMAT ('Prandtl number = ', F5.3)
 109  FORMAT ('Total run time = ', F8.4)
 110  FORMAT ('Total amount of coarse grid time = ', F6.4)
 111  FORMAT ('Ratio of specific heats, gamma = ', F6.4)
 112  FORMAT ('Freestream temperature = ', F12.3)
 113  FORMAT ('Freestream pressure = ', F12.2)
 114  FORMAT ('Freestream velocity = ', F8.4)
 115  FORMAT ('BGK time constant, tau = ', F4.2)
 116  FORMAT ('Characteristic length = ', F11.5)
 117  FORMAT ('Turbulence type is ', A8)
 118  FORMAT ('The simulation fluid is ', A8)
 122  FORMAT ('This simulation will use ',I3,' CGNS files for input.')


 208  FORMAT ('Voxel size = ',F9.7)

END SUBROUTINE
