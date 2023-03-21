      SUBROUTINE inside_init(link_num,init_number,highest_priority,&
        inside,overlap_counter)
!!
!! This determines if the point is inside or outside the initialization
!! region and what to do if it is or isn't.  If it is, set its data to
!! that of the region, if not, move on.  If multiple regions overlap for
!! a point, do a running average with the other points.
!!
!! Called by: init_fi
!! Calls: error_out
!!
!      USE nml_init
!      USE linkwise
!      USE grid_data
!      USE loop_data
!      USE constants
!      USE precise
!      IMPLICIT NONE
!      INTEGER :: init_number,highest_priority,overlap_counter,link_num
!      REAL(KIND = dp) :: radius,x_dist,y_dist,z_dist,u_comp,v_comp,&
!        w_comp
!      REAL(KIND = dp) :: x_center,y_center,z_center,rho_local
!      REAL(KIND = dp) :: x_max,y_max,z_max,x_min,y_min,z_min,dist_to_pt
!      LOGICAL :: inside
!
!!      gas_const_R = 287.1
!
!      IF (highest_priority < init_info(2,init_number)) THEN
!
!        highest_priority = init_info(2,init_number)
!        x_dist = x_geom_max - x_geom_min
!        y_dist = y_geom_max - y_geom_min
!        z_dist = z_geom_max - z_geom_min
!        x_center = (x_geom_max-x_geom_min)/2.0D0
!        y_center = (y_geom_max-y_geom_min)/2.0D0
!        z_center = (z_geom_max-z_geom_min)/2.0D0
!
!        IF (init_info(1,init_number) == 1) THEN
!!
!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                  <= y_max) THEN
!              IF (link_grid(3,link_num) >= z_min .AND. &
!                    link_grid(3,link_num) <= z_max) THEN
!
!                rho(link_num) = init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!              END IF
!            END IF
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 2) THEN
!!
!          x_min = init_bounds(1,init_number)
!          y_min = init_bounds(2,init_number)
!          z_min = init_bounds(3,init_number)
!          x_max = init_bounds(4,init_number)
!          y_max = init_bounds(5,init_number)
!          z_max = init_bounds(6,init_number)
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                  <= y_max) THEN
!              IF (link_grid(3,link_num) >= z_min .AND. &
!                    link_grid(3,link_num) <= z_max) THEN
!
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!              END IF
!            END IF
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 3) THEN
!!
!          x_max = init_bounds(1,init_number) + (init_bounds(4,init_number)&
!                    /2.0D0)*x_dist
!          y_max = init_bounds(2,init_number) + (init_bounds(5,init_number)&
!                    /2.0D0)*y_dist
!          z_max = init_bounds(3,init_number) + (init_bounds(6,init_number)&
!                    /2.0D0)*z_dist
!          x_min = init_bounds(1,init_number) - (init_bounds(4,init_number)&
!                    /2.0D0)*x_dist
!          y_min = init_bounds(2,init_number) - (init_bounds(5,init_number)&
!                    /2.0D0)*y_dist
!          z_min = init_bounds(3,init_number) - (init_bounds(6,init_number)&
!                    /2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                  <= y_max) THEN
!              IF (link_grid(3,link_num) >= z_min .AND. &
!                    link_grid(3,link_num) <= z_max) THEN
!
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!              END IF
!            END IF
!          END IF
!!
!! Spherical bounds
!!
!        ELSE IF (init_info(1,init_number) == 11) THEN
!!
!!          x_grid_max = bounds(1) + bounds(4)*longest_dist
!!          y_grid_max = bounds(2) + bounds(4)*longest_dist
!!          z_grid_max = bounds(3) + bounds(4)*longest_dist
!!          x_grid_min = bounds(1) - bounds(4)*longest_dist
!!          y_grid_min = bounds(2) - bounds(4)*longest_dist
!!          z_grid_min = bounds(3) - bounds(4)*longest_dist
!          x_center = init_bounds(1,init_number)
!          y_center = init_bounds(2,init_number)
!          z_center = init_bounds(3,init_number)
!          radius = init_bounds(4,init_number)*longest_dist
!          dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!            init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!            init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!            init_number))**2.0D0)
!          IF (dist_to_pt <= radius) THEN
!            rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!            u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!            v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!            w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 12) THEN
!!
!!          x_grid_max = bounds(1) + bounds(4)
!!          y_grid_max = bounds(2) + bounds(4)
!!          z_grid_max = bounds(3) + bounds(4)
!!          x_grid_min = bounds(1) - bounds(4)
!!          y_grid_min = bounds(2) - bounds(4)
!!          z_grid_min = bounds(3) - bounds(4)
!          x_center = init_bounds(1,init_number)
!          y_center = init_bounds(2,init_number)
!          z_center = init_bounds(3,init_number)
!          radius = init_bounds(4,init_number)
!          dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!            init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!            init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!            init_number))**2.0D0)
!          IF (dist_to_pt <= radius) THEN
!            rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!            u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!            v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!            w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!          END IF
!!
!! Cylindrical bounds
!!
!        ELSE IF (init_info(1,init_number) == 21) THEN
!!
!! If the cylinder moves in the positive x-direction
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!! Test to see if it exceeds the x-bounds on the cylinder
!            IF (link_grid(1,link_num) >= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) <= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!!            x_grid_min = bounds(1)
!!            x_grid_max = bounds(1) + bounds(5)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1)
!!            x_grid_min = bounds(1) + bounds(5)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!            IF (link_grid(1,link_num) <= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) >= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!
!        ELSE IF (init_info(1,init_number) == 22) THEN
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) + bounds(5)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) >= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) <= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) - bounds(5)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) <= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) >= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!
!
!        ELSE IF (init_info(1,init_number) == 23) THEN
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3)
!!            z_grid_max = bounds(3) + bounds(5)
!
!            IF (link_grid(3,link_num) >= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) <= (init_bounds(3,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/&
!                     (gas_const_R*init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) + bounds(5)
!!            z_grid_max = bounds(3)
!
!            IF (link_grid(3,link_num) <= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) >= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!
!        ELSE IF (init_info(1,init_number) == 24) THEN
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1)
!!            x_grid_max = bounds(1) + bounds(5)*x_dist
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(1,link_num) >= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) <= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number)*x_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1)
!!            x_grid_min = bounds(1) + bounds(5)*x_dist
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(1,link_num) <= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) >= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number)*x_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 25) THEN
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) + bounds(5)*y_dist
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) >= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) <= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number)*y_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) - bounds(5)*y_dist
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) <= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) >= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number)*y_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 26) THEN
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3)
!!            z_grid_max = bounds(3) + bounds(5)*z_dist
!
!            IF (link_grid(3,link_num) >= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) <= (init_bounds(3,init_number)&
!                  +init_bounds(5,init_number)*z_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) + bounds(5)*z_dist
!!            z_grid_max = bounds(3)
!
!            IF (link_grid(3,link_num) <= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) >= (init_bounds(3,init_number)&
!                  +init_bounds(5,init_number)*z_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!!
!! Planar bounds
!!
!        ELSE IF (init_info(1,init_number) == -1) THEN
!! For x-y plane for 2D
!!          x_grid_min = bounds(1)
!!          x_grid_max = bounds(2)
!!          y_grid_min = bounds(3)
!!          y_grid_max = bounds(4)
!!          z_grid_min = bounds(5)
!!          z_grid_max = bounds(5)
!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND.link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(1,link_num) >= y_min .AND. link_grid(1,&
!                  link_num) <= y_max) THEN
!
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = 0.0D0
!
!            END IF
!          END IF
!
!!
!        ELSE IF (init_info(1,init_number) == -2) THEN
!! For x-z plane for 2D
!!          x_grid_min = bounds(1)
!!          x_grid_max = bounds(2)
!!          z_grid_min = bounds(3)
!!          z_grid_max = bounds(4)
!!          y_grid_min = bounds(5)
!!          y_grid_max = bounds(5)
!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(3,link_num) >= z_min .AND. link_grid(3,link_num)&
!                  <= z_max) THEN
!
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = init_data(3,init_number)*init_data(4,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                v_vel(link_num) = 0.0D0
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!            END IF
!          END IF
!
!!
!        ELSE IF (init_info(1,init_number) == -3) THEN
!! For y-z plane for 2D
!!          y_grid_min = bounds(1)
!!          y_grid_max = bounds(2)
!!          z_grid_min = bounds(3)
!!          z_grid_max = bounds(4)
!!          x_grid_min = bounds(5)
!!          x_grid_max = bounds(5)
!!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                <= y_max) THEN
!            IF (link_grid(3,link_num) >= z_min .AND. link_grid(3,link_num)&
!                  <= z_max) THEN
!
!                rho(link_num) = init_data(2,init_number)/(gas_const_R*&
!                    init_data(1,init_number))
!                u_vel(link_num) = 0.0D0
!                v_vel(link_num) = init_data(3,init_number)*init_data(5,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!                w_vel(link_num) = init_data(3,init_number)*init_data(6,&
!                      init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                      init_data(5,init_number)**2.0D0+init_data(6,&
!                      init_number)**2.0D0)
!
!            END IF
!          END IF
!!
!        ELSE
!          WRITE(*,*) 'Failure applying boundaries, exiting.'
!          CALL error_out
!        END IF
!!
!! In the case where there are two overlapping regions with the same
!! priority
!!
!      ELSE IF (highest_priority == init_info(2,init_number) .AND. &
!        init_info(2,init_number) /= 0) THEN
!!
!!
!! This section will be the same as above but with weighted averages if
!! the priority is the same for two nodes.  The format will go like this
!! for the weighted average where n=number of overlapping regions
!!
!! value = old_value*(n-1/n)+ new_value*(1/n)
!!
!!
!        x_dist = x_geom_max - x_geom_min
!        y_dist = y_geom_max - y_geom_min
!        z_dist = z_geom_max - z_geom_min
!        x_center = (x_geom_max-x_geom_min)/2.0D0
!        y_center = (y_geom_max-y_geom_min)/2.0D0
!        z_center = (z_geom_max-z_geom_min)/2.0D0
!
!        IF (init_info(1,init_number) == 1) THEN
!!
!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND.link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                  <= y_max) THEN
!              IF (link_grid(3,link_num) >= z_min .AND. &
!                    link_grid(3,link_num) <= z_max) THEN
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!              END IF
!            END IF
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 2) THEN
!!
!          x_min = init_bounds(1,init_number)
!          y_min = init_bounds(2,init_number)
!          z_min = init_bounds(3,init_number)
!          x_max = init_bounds(4,init_number)
!          y_max = init_bounds(5,init_number)
!          z_max = init_bounds(6,init_number)
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                  <= y_max) THEN
!              IF (link_grid(3,link_num) >= z_min .AND. &
!                    link_grid(3,link_num) <= z_max) THEN
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 3) THEN
!!
!          x_max = init_bounds(1,init_number) + (init_bounds(4,init_number)&
!                    /2.0D0)*x_dist
!          y_max = init_bounds(2,init_number) + (init_bounds(5,init_number)&
!                    /2.0D0)*y_dist
!          z_max = init_bounds(3,init_number) + (init_bounds(6,init_number)&
!                    /2.0D0)*z_dist
!          x_min = init_bounds(1,init_number) - (init_bounds(4,init_number)&
!                    /2.0D0)*x_dist
!          y_min = init_bounds(2,init_number) - (init_bounds(5,init_number)&
!                    /2.0D0)*y_dist
!          z_min = init_bounds(3,init_number) - (init_bounds(6,init_number)&
!                    /2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                  <= y_max) THEN
!              IF (link_grid(3,link_num) >= z_min .AND. &
!                    link_grid(3,link_num) <= z_max) THEN
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!          END IF
!!
!! Spherical bounds
!!
!        ELSE IF (init_info(1,init_number) == 11) THEN
!!
!!          x_grid_max = bounds(1) + bounds(4)*longest_dist
!!          y_grid_max = bounds(2) + bounds(4)*longest_dist
!!          z_grid_max = bounds(3) + bounds(4)*longest_dist
!!          x_grid_min = bounds(1) - bounds(4)*longest_dist
!!          y_grid_min = bounds(2) - bounds(4)*longest_dist
!!          z_grid_min = bounds(3) - bounds(4)*longest_dist
!          x_center = init_bounds(1,init_number)
!          y_center = init_bounds(2,init_number)
!          z_center = init_bounds(3,init_number)
!          radius = init_bounds(4,init_number)*longest_dist
!          dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!            init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!            init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!            init_number))**2.0D0)
!
!          IF (dist_to_pt <= radius) THEN
!            overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 12) THEN
!!
!!          x_grid_max = bounds(1) + bounds(4)
!!          y_grid_max = bounds(2) + bounds(4)
!!          z_grid_max = bounds(3) + bounds(4)
!!          x_grid_min = bounds(1) - bounds(4)
!!          y_grid_min = bounds(2) - bounds(4)
!!          z_grid_min = bounds(3) - bounds(4)
!          x_center = init_bounds(1,init_number)
!          y_center = init_bounds(2,init_number)
!          z_center = init_bounds(3,init_number)
!          radius = init_bounds(4,init_number)
!          dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!            init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!            init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!            init_number))**2.0D0)
!          IF (dist_to_pt <= radius) THEN
!
!           overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!          END IF
!!
!! Cylindrical bounds
!!
!        ELSE IF (init_info(1,init_number) == 21) THEN
!!
!! If the cylinder moves in the positive x-direction
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!! Test to see if it exceeds the x-bounds on the cylinder
!            IF (link_grid(1,link_num) >= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) <= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!!            x_grid_min = bounds(1)
!!            x_grid_max = bounds(1) + bounds(5)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1)
!!            x_grid_min = bounds(1) + bounds(5)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!            IF (link_grid(1,link_num) <= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) >= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!
!        ELSE IF (init_info(1,init_number) == 22) THEN
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) + bounds(5)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) >= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) <= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) - bounds(5)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) <= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) >= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!
!
!        ELSE IF (init_info(1,init_number) == 23) THEN
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3)
!!            z_grid_max = bounds(3) + bounds(5)
!
!            IF (link_grid(3,link_num) >= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) <= (init_bounds(3,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) + bounds(5)
!!            z_grid_max = bounds(3)
!
!            IF (link_grid(3,link_num) <= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) >= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number))) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!
!        ELSE IF (init_info(1,init_number) == 24) THEN
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1)
!!            x_grid_max = bounds(1) + bounds(5)*x_dist
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(1,link_num) >= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) <= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number)*x_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1)
!!            x_grid_min = bounds(1) + bounds(5)*x_dist
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(1,link_num) <= init_bounds(1,init_number) .AND.&
!                  link_grid(1,link_num) >= (init_bounds(1,init_number)&
!                  +init_bounds(5,init_number)*x_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 25) THEN
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) + bounds(5)*y_dist
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) >= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) <= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number)*y_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2)
!!            y_grid_min = bounds(2) - bounds(5)*y_dist
!!            z_grid_min = bounds(3) - bounds(4)
!!            z_grid_max = bounds(3) + bounds(4)
!
!            IF (link_grid(2,link_num) <= init_bounds(2,init_number) .AND.&
!                  link_grid(2,link_num) >= (init_bounds(2,init_number)&
!                  +init_bounds(5,init_number)*y_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(3,link_num)-init_bounds(3,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!!
!        ELSE IF (init_info(1,init_number) == 26) THEN
!!
!          IF (init_bounds(5,init_number) > 0.0) THEN
!!            x_grid_min = bounds(1) - bounds(4)
!!            x_grid_max = bounds(1) + bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3)
!!            z_grid_max = bounds(3) + bounds(5)*z_dist
!
!            IF (link_grid(3,link_num) >= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) <= (init_bounds(3,init_number)&
!                  +init_bounds(5,init_number)*z_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE IF (init_bounds(5,init_number) < 0.0) THEN
!!            x_grid_max = bounds(1) + bounds(4)
!!            x_grid_min = bounds(1) - bounds(4)
!!            y_grid_max = bounds(2) + bounds(4)
!!            y_grid_min = bounds(2) - bounds(4)
!!            z_grid_min = bounds(3) + bounds(5)*z_dist
!!            z_grid_max = bounds(3)
!
!            IF (link_grid(3,link_num) <= init_bounds(3,init_number) .AND.&
!                  link_grid(3,link_num) >= (init_bounds(3,init_number)&
!                  +init_bounds(5,init_number)*z_dist)) THEN
!! Test to see if it's inside the radius of the cylinder
!              dist_to_pt = SQRT((link_grid(1,link_num)-init_bounds(1,&
!                init_number))**2.0D0+(link_grid(2,link_num)-init_bounds(2,&
!                init_number))**2.0D0)
!              radius = init_bounds(4,init_number)
!              IF (dist_to_pt <= radius) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!              END IF
!            END IF
!
!          ELSE
!            WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!            CALL error_out
!          END IF
!!
!! Planar bounds
!!
!        ELSE IF (init_info(1,init_number) == -1) THEN
!! For x-y plane for 2D
!!          x_grid_min = bounds(1)
!!          x_grid_max = bounds(2)
!!          y_grid_min = bounds(3)
!!          y_grid_max = bounds(4)
!!          z_grid_min = bounds(5)
!!          z_grid_max = bounds(5)
!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(1,link_num) >= y_min .AND. link_grid(1,link_num)&
!                  <= y_max) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = 0.0D0
!
!            END IF
!          END IF
!
!!
!        ELSE IF (init_info(1,init_number) == -2) THEN
!! For x-z plane for 2D
!!          x_grid_min = bounds(1)
!!          x_grid_max = bounds(2)
!!          z_grid_min = bounds(3)
!!          z_grid_max = bounds(4)
!!          y_grid_min = bounds(5)
!!          y_grid_max = bounds(5)
!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(1,link_num) >= x_min .AND. link_grid(1,link_num)&
!                <= x_max) THEN
!            IF (link_grid(3,link_num) >= z_min .AND. link_grid(3,link_num)&
!                  <= z_max) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = u_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                  (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(4,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                v_vel(link_num) = 0.0D0
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!            END IF
!          END IF
!
!!
!        ELSE IF (init_info(1,init_number) == -3) THEN
!! For y-z plane for 2D
!!          y_grid_min = bounds(1)
!!          y_grid_max = bounds(2)
!!          z_grid_min = bounds(3)
!!          z_grid_max = bounds(4)
!!          x_grid_min = bounds(5)
!!          x_grid_max = bounds(5)
!!          x_max = x_center + (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_max = y_center + (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_max = z_center + (init_bounds(3,init_number)/2.0D0)*z_dist
!!          x_min = x_center - (init_bounds(1,init_number)/2.0D0)*x_dist
!          y_min = y_center - (init_bounds(2,init_number)/2.0D0)*y_dist
!          z_min = z_center - (init_bounds(3,init_number)/2.0D0)*z_dist
!
!          IF (link_grid(2,link_num) >= y_min .AND. link_grid(2,link_num)&
!                <= y_max) THEN
!            IF (link_grid(3,link_num) >= z_min .AND. link_grid(3,link_num)&
!                  <= z_max) THEN
!
!                overlap_counter = overlap_counter + 1
!
!                rho(link_num) = rho(link_num)*((DBLE(overlap_counter)-&
!                 1.0D0)/DBLE(overlap_counter))+(1.0D0/DBLE(&
!                  overlap_counter))*(init_data(2,init_number)/&
!                  (gas_const_R*init_data(1,init_number)))
!
!                u_vel(link_num) = 0.0D0
!
!                v_vel(link_num) = v_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))+&
!                 (1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(5,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!                w_vel(link_num) = w_vel(link_num)*((DBLE(&
!                  overlap_counter)-1.0D0)/DBLE(overlap_counter))&
!                  +(1.0D0/DBLE(overlap_counter))*&
!                  (init_data(3,init_number)*init_data(6,&
!                  init_number)/SQRT(init_data(4,init_number)**2.0D0+&
!                  init_data(5,init_number)**2.0D0+init_data(6,&
!                  init_number)**2.0D0))
!
!            END IF
!          END IF
!!
!        ELSE
!          WRITE(*,*) 'Failure applying boundaries, exiting.'
!          CALL error_out
!        END IF
!
!
!      END IF
!
END SUBROUTINE
