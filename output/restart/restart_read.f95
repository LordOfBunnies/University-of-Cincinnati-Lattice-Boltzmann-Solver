SUBROUTINE restart_read
!
! This is to read the restart file of a previously run simulation.
! It's essentially the inverse of the restart_write code, though this
! requires doing allocations as you go.
!
! Calls: error_out
! Called by: uclbs_main
!
!USE linkwise
!USE loop_data
!USE nml_inlet_outlet
!USE nml_output
!USE nml_geometry
!USE geom_data
!USE output_data
!USE precise
!USE timing
!USE startup
!USE constants
!USE freestream_values
!USE grid_data
!IMPLICIT NONE
!INTEGER :: ier,t,a,i,istat
!CHARACTER(len=80) :: filename_out
!filename_out = 'uclbs.out'
!
!OPEN(FILE=filename_out,UNIT = 11,STATUS='OLD',FORM = 'FORMATTED',&
!  ACTION='WRITE',POSITION='APPEND',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
!IF (istat /= 0) THEN
!  WRITE(*,*) 'Problems opening output file, closing.'
!  WRITE(*,*) istat
!  CALL error_out
!END IF
!
!OPEN(UNIT = 254,FILE='uclbs.restart',STATUS='OLD',&
!  ACTION='READ',FORM='UNFORMATTED',IOSTAT=ier)
!IF (ier /= 0) THEN
!  CALL error_out
!END IF
!
!!
!!
!!
!READ(254) a,t
!READ(254) total_tris
!READ(254)
!READ(254) link_number
!READ(254)
!READ(254) num_inlets,num_outlets
!READ(254) timesteps,t
!READ(254) dt,cs
!READ(254)
!READ(254) num_mics,num_restart_saves,num_saves
!READ(254)
!READ(254) total_time,coarse_time,comp_time
!READ(254)
!READ(254) x_geom_max,y_geom_max,z_geom_max,x_geom_min, &
!          y_geom_min,z_geom_min
!READ(254) x_grid_max,y_grid_max,z_grid_max,x_grid_min, &
!          y_grid_min,z_grid_min
!READ(254) dx,longest_dist
!READ(254) num_x,num_y,num_z,dir
!READ(254)
!
!WRITE(11,101) a
!WRITE(11,102) total_tris
!WRITE(11,103) link_number
!WRITE(11,104) num_inlets,num_outlets
!WRITE(11,105) timesteps,t
!WRITE(11,106) dt,cs
!WRITE(11,107) num_mics,num_restart_saves,num_saves
!WRITE(11,108) total_time,coarse_time,comp_time
!WRITE(11,109) x_geom_max,y_geom_max,z_geom_max,x_geom_min, &
!          y_geom_min,z_geom_min
!WRITE(11,110) x_grid_max,y_grid_max,z_grid_max,x_grid_min, &
!          y_grid_min,z_grid_min
!WRITE(11,111) dx,longest_dist
!WRITE(11,112) num_x,num_y,num_z,dir
!!
!! Allocation break!!! Feel free to open up into a guitar solo during this time
!!
!! Start with link information
!ALLOCATE(link(dir,link_number))
!ALLOCATE(link_grid(3,link_number))
!ALLOCATE(boundary_node(link_number))
!! Geometry info
!ALLOCATE(normal(3,total_tris))
!ALLOCATE(tri_vertices(3,3,total_tris))
!! Flow information
!ALLOCATE(u_vel(link_number))
!ALLOCATE(v_vel(link_number))
!ALLOCATE(w_vel(link_number))
!ALLOCATE(rho(link_number))
!ALLOCATE(fi(dir,link_number))
!ALLOCATE(fout(dir,link_number))
!! Allocate directions and weights
!ALLOCATE(wi(dir))
!ALLOCATE(cx(dir))
!ALLOCATE(cy(dir))
!ALLOCATE(cz(dir))
!ALLOCATE(opp(dir))
!! Inlet/outlet info
!ALLOCATE(inlet_type(2,num_inlets))
!ALLOCATE(outlet_type(2,num_outlets))
!ALLOCATE(in_name(num_inlets))
!ALLOCATE(out_name(num_outlets))
!ALLOCATE(inlet_flow_val(2,num_inlets))
!ALLOCATE(outlet_flow_val(2,num_outlets))
!! Microphone info
!ALLOCATE(mic_locations(3,num_mics))
!ALLOCATE(mic_times(2,num_mics))
!ALLOCATE(mic_type(6,num_mics))
!ALLOCATE(mic_frequency(num_mics))
!! Save data
!ALLOCATE(save_shapes(num_saves))
!ALLOCATE(save_times(2,num_saves))
!ALLOCATE(save_regions(6,num_saves))
!ALLOCATE(restart_saves(num_restart_saves))
!ALLOCATE(save_frequency(num_saves))
!ALLOCATE(save_data(6,num_saves))
!ALLOCATE(num_nodes_in_save(num_saves))
!!
!! Begin the read of the data
!!
!!
!!
!READ(254) tri_vertices,normal
!READ(254)
!READ(254) link
!READ(254) link_grid
!READ(254) boundary_node
!READ(254) wi,cx,cy,cz,opp
!READ(254)
!READ(254) fi
!READ(254) rho
!READ(254) u_vel
!READ(254) v_vel
!READ(254) w_vel
!READ(254)
!READ(254) inlet_type
!READ(254) outlet_type
!READ(254) inlet_flow_val
!READ(254) outlet_flow_val
!READ(254) in_name
!READ(254) out_name
!READ(254)
!READ(254) mic_locations
!READ(254) mic_times
!READ(254) mic_type
!READ(254) mic_frequency
!READ(254) save_times
!READ(254) save_regions
!!      READ(254) restart_saves
!READ(254) save_shapes
!READ(254) save_frequency
!READ(254) restart_timesteps
!READ(254) save_data
!READ(254)
!DO i = 1,num_saves
!
!  SELECT CASE (i)
!  CASE(1)
!    ALLOCATE(vol_info_2(num_nodes_in_save(i)))
!    READ(254) vol_info_2
!  CASE(2)
!    ALLOCATE(vol_info_3(num_nodes_in_save(i)))
!    READ(254) vol_info_3
!  CASE(3)
!    ALLOCATE(vol_info_4(num_nodes_in_save(i)))
!    READ(254) vol_info_4
!  CASE(4)
!    ALLOCATE(vol_info_4(num_nodes_in_save(i)))
!    READ(254) vol_info_4
!  CASE(5)
!    ALLOCATE(vol_info_5(num_nodes_in_save(i)))
!    READ(254) vol_info_5
!  CASE(6)
!    ALLOCATE(vol_info_6(num_nodes_in_save(i)))
!    READ(254) vol_info_6
!  CASE(7)
!    ALLOCATE(vol_info_7(num_nodes_in_save(i)))
!    READ(254) vol_info_7
!  CASE(8)
!    ALLOCATE(vol_info_8(num_nodes_in_save(i)))
!    READ(254) vol_info_8
!  CASE(9)
!    ALLOCATE(vol_info_9(num_nodes_in_save(i)))
!    READ(254) vol_info_9
!  CASE(10)
!    ALLOCATE(vol_info_10(num_nodes_in_save(i)))
!    READ(254) vol_info_10
!  CASE(11)
!    ALLOCATE(vol_info_11(num_nodes_in_save(i)))
!    READ(254) vol_info_11
!  CASE(12)
!    ALLOCATE(vol_info_12(num_nodes_in_save(i)))
!    READ(254) vol_info_12
!  CASE(13)
!    ALLOCATE(vol_info_13(num_nodes_in_save(i)))
!    READ(254) vol_info_13
!  CASE(14)
!    ALLOCATE(vol_info_14(num_nodes_in_save(i)))
!    READ(254) vol_info_14
!  CASE(15)
!    ALLOCATE(vol_info_15(num_nodes_in_save(i)))
!    READ(254) vol_info_15
!  CASE(16)
!    ALLOCATE(vol_info_16(num_nodes_in_save(i)))
!    READ(254) vol_info_16
!  CASE(17)
!    ALLOCATE(vol_info_17(num_nodes_in_save(i)))
!    READ(254) vol_info_17
!  CASE(18)
!    ALLOCATE(vol_info_18(num_nodes_in_save(i)))
!    READ(254) vol_info_18
!  CASE(19)
!    ALLOCATE(vol_info_19(num_nodes_in_save(i)))
!    READ(254) vol_info_19
!  CASE(20)
!    ALLOCATE(vol_info_20(num_nodes_in_save(i)))
!    READ(254) vol_info_20
!  CASE(21)
!    ALLOCATE(vol_info_21(num_nodes_in_save(i)))
!    READ(254) vol_info_21
!  CASE(22)
!    ALLOCATE(vol_info_22(num_nodes_in_save(i)))
!    READ(254) vol_info_22
!  CASE(23)
!    ALLOCATE(vol_info_23(num_nodes_in_save(i)))
!    READ(254) vol_info_23
!  CASE(24)
!    ALLOCATE(vol_info_24(num_nodes_in_save(i)))
!    READ(254) vol_info_24
!  CASE(25)
!    ALLOCATE(vol_info_25(num_nodes_in_save(i)))
!    READ(254) vol_info_25
!  END SELECT
!END DO
!CLOSE(254)
!CLOSE(11)
!
!t_start = t
!
! 101  FORMAT ('Starting from restart ',I2)
! 102  FORMAT ('Total number of triangles in geometry =',I8)
! 103  FORMAT ('Total number of active fluid nodes =', I12)
! 104  FORMAT ('Number of inlets =',I4,' and number of outlets =',I4)
! 105  FORMAT ('Total number of timesteps =',I8,&
!        'and restarting at timestep ',I8)
! 106  FORMAT ('Timestep is approximately ',F10.8,'seconds. Lattice',&
!        ' speed of sound is =',F9.7)
! 107  FORMAT ('Number of microphones =',I4/,'Number of restarts ',I2/,&
!        'Number of saves ',I3)
! 108  FORMAT ('Total time of simulation ',F8.5,'. Coarse time ',&
!        F8.5,'. Computational time ',F11.5)
! 109  FORMAT ('Geometric max, X-direction ',F12.7/,&
!      'Geometric max, Y-direction ',F12.7/,'Geometric max, Z-direction ',&
!      F12.7/,'Geometric min, X-direction ',F12.7/,'Geometric min, Y-direction ',&
!      F12.7/,'Geometric min, Z-direction ',F12.7)
! 110  FORMAT ('Grid max, X-direction ',F12.7/,&
!        'Grid max, Y-direction ',F12.7/,'Grid max, Z-direction ',&
!        F12.7/,'Grid min, X-direction ',F12.7/,'Grid min, Y-direction ',&
!        F12.7/,'Grid min, Z-direction ',F12.7)
! 111  FORMAT ('Base voxel size =',F15.12/,'Longest geometric distance =',F14.10)
! 112  FORMAT ('Number of cells in the x-direction =', I7/,&
!        'Number of cells in the y-direction =', I7/,&
!        'Number of cells in the z-direction =', I7/,&
!        'Number of directions of the lattice =', I7)
!! 113
! 201  FORMAT ('Beginning read of restart file ', I8, &
!        ' at timestep ', I8,'.')

      END SUBROUTINE restart_read
