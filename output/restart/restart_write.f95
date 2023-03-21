SUBROUTINE restart_write(a)
!
! This writes all data necessary to restart a case to a restart file.
!
! Called by: lb_output
! Calls: error_out
! XXXX
USE linkwise
USE nml_inlet_outlet
USE nml_output
USE nml_geometry
USE geom_data
USE output_data
USE precise
USE timing
USE constants
USE freestream_values
USE grid_data
IMPLICIT NONE
INTEGER :: ier,t,a,i

!OPEN(UNIT = 254,FILE='uclbs.restart',STATUS='REPLACE',&
!  ACTION='WRITE',FORM='UNFORMATTED',IOSTAT=ier)
!IF (ier /= 0) THEN
!  CALL error_out
!END IF
!!
!! All the individual numbers used for the run. The numbers may or may
!! not actually be used but they COULD be important. Many are needed
!! for allocation before reading the second part of the file.
!!
!WRITE(11,*) ''
!WRITE(11,201) a,t
!!
!!
!!
!WRITE(254) a
!WRITE(254) total_tris
!WRITE(254)
!WRITE(254) link_number
!WRITE(254)
!WRITE(254) num_inlets,num_outlets
!WRITE(254) timesteps,t
!WRITE(254) dt,cs
!WRITE(254)
!WRITE(254) num_mics,num_restart_saves,num_saves
!WRITE(254)
!WRITE(254) total_time,coarse_time,comp_time
!WRITE(254)
!WRITE(254) x_geom_max,y_geom_max,z_geom_max,x_geom_min, &
!          y_geom_min,z_geom_min
!WRITE(254) x_grid_max,y_grid_max,z_grid_max,x_grid_min, &
!          y_grid_min,z_grid_min
!WRITE(254) dx,longest_dist
!WRITE(254) num_x,num_y,num_z,dir
!WRITE(254)
!!
!! This is all the allocated stuff that needs to be read. This is the meat
!! of the data.
!!
!      WRITE(254) tri_vertices,normal
!      WRITE(254)
!      WRITE(254) link
!      WRITE(254) link_grid
!      WRITE(254) boundary_node
!      WRITE(254) wi,cx,cy,cz,opp
!      WRITE(254)
!      WRITE(254) fi
!      WRITE(254) rho
!      WRITE(254) u_vel
!      WRITE(254) v_vel
!      WRITE(254) w_vel
!      WRITE(254)
!      WRITE(254) inlet_type
!      WRITE(254) outlet_type
!      WRITE(254) inlet_flow_val
!      WRITE(254) outlet_flow_val
!      WRITE(254) in_name
!      WRITE(254) out_name
!      WRITE(254)
!      WRITE(254) mic_locations
!      WRITE(254) mic_times
!      WRITE(254) mic_type
!      WRITE(254) mic_frequency
!!      WRITE(254) num_nodes_in_save
!      WRITE(254) save_times
!      WRITE(254) save_regions
!!      WRITE(254) restart_saves
!      WRITE(254) save_shapes
!      WRITE(254) save_frequency
!      WRITE(254) restart_timesteps
!      WRITE(254) save_data
!      WRITE(254)
!
!!      DO i = 1,num_saves
!!        SELECT CASE (i)
!!        CASE(1)
!!          WRITE(254) vol_info_1
!!        CASE(2)
!!          WRITE(254) vol_info_2
!!        CASE(3)
!!          WRITE(254) vol_info_3
!!        CASE(4)
!!          WRITE(254) vol_info_4
!!        CASE(5)
!!          WRITE(254) vol_info_5
!!        CASE(6)
!!          WRITE(254) vol_info_6
!!        CASE(7)
!!          WRITE(254) vol_info_7
!!        CASE(8)
!!          WRITE(254) vol_info_8
!!        CASE(9)
!!          WRITE(254) vol_info_9
!!        CASE(10)
!!          WRITE(254) vol_info_10
!!        CASE(11)
!!          WRITE(254) vol_info_11
!!        CASE(12)
!!          WRITE(254) vol_info_12
!!        CASE(13)
!!          WRITE(254) vol_info_13
!!        CASE(14)
!!          WRITE(254) vol_info_14
!!        CASE(15)
!!          WRITE(254) vol_info_15
!!        CASE(16)
!!          WRITE(254) vol_info_16
!!        CASE(17)
!!          WRITE(254) vol_info_17
!!        CASE(18)
!!          WRITE(254) vol_info_18
!!        CASE(19)
!!          WRITE(254) vol_info_19
!!        CASE(20)
!!          WRITE(254) vol_info_20
!!        CASE(21)
!!          WRITE(254) vol_info_21
!!        CASE(22)
!!          WRITE(254) vol_info_22
!!        CASE(23)
!!          WRITE(254) vol_info_23
!!        CASE(24)
!!          WRITE(254) vol_info_24
!!        CASE(25)
!!          WRITE(254) vol_info_25
!!        END SELECT
!!
!!      END DO
!
!
!
!!      WRITE(254)
!!      WRITE(254)
!      WRITE(11,*) 'Restart write done, continuing timestepping.'
!      WRITE(11,*) ''
!
 201  FORMAT ('Beginning write of restart file ', I8, &
        ' at timestep ', I8,'.')
      END SUBROUTINE
