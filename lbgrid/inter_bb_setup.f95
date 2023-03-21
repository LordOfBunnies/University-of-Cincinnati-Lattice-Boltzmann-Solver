subroutine inter_bb_setup(dir_state,rest_state,t,a,b,c,lvl)
!
! This is to set up the interpolative bounceback to be used later. It will
! assign the the grid coordinates to wall nodes and determine the distance
! to the wall, then find that as a percentage of dx.
!
! Based on fluid_or_solid and state_directions_cgns
!
! Called by: state_directions
! Calls: !!!!wall_lattice
!
USE precise
USE grid_data, only : rdx
USE geom_data
USE linkwise
USE nml_inlet_outlet
use amr_processes, only : get_real_coords,self
use amrex_amr_module
use amrex_base_module
IMPLICIT NONE
INTEGER :: a,b,c,i,t,triangle_area,volume_tri,lvl
REAL(KIND=dp) :: x_tri_max,x_tri_min,y_tri_max,y_tri_min,z_tri_max
REAL(KIND=dp) :: z_tri_min,upper,lower,dot,planar_coeff,closest_diff_lo,closest_diff_hi
!      REAL(KIND=dp) :: ray_intersect_length
INTEGER :: cross_counter,edge_counter,outside,final_counter
INTEGER :: fluid_target,solid_target,null_target,source
integer :: dir_state,rest_state,planar_coeff_counter
LOGICAL :: closeness,point_moved,dir_assigned
REAL(KIND=dp) :: ray_end1(3),ray_end2(3),ray_crossing(3)
!      REAL(KIND = dp) :: loc_dist

null_target = 0
fluid_target = 1
solid_target = 0
source = 2
!
! Array to tell whether you should skip the long procedure to determine
! the link quality and go straight to assigning the link number. The skip
! is used if the direction exceeds the boundaries of the domain or the
! node in that direction is not the same type as the current one
!
point_moved = .FALSE.

ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))

if (.not. shifted) then
ray_end2(1) = ray_end1(1) + cx(t)*amrex_geom(lvl)%dx(1)
ray_end2(2) = ray_end1(2) + cy(t)*amrex_geom(lvl)%dx(2)
ray_end2(3) = ray_end1(3) + cz(t)*amrex_geom(lvl)%dx(3)
else
ray_end2(1) = ray_end1(1) - cx(t)*amrex_geom(lvl)%dx(1)
ray_end2(2) = ray_end1(2) - cy(t)*amrex_geom(lvl)%dx(2)
ray_end2(3) = ray_end1(3) - cz(t)*amrex_geom(lvl)%dx(3)
end if

dir_assigned = .FALSE.

!              WRITE(*,*) 'Entering triangle loop.'
closest_diff_lo = 2.0D0
closest_diff_hi = 2.0D0
planar_coeff_counter = 0
cross_counter = 0
edge_counter = 0
final_counter = 0
volume_tri = 0
closeness = .FALSE.

triangle_loop: DO i = 1,total_tris
!
! For cgns starts, we must determine if the triangle is viable
!
  x_tri_max = MAX(tri_vertices(1,1,i),tri_vertices(1,2,i),&
    tri_vertices(1,3,i))
  x_tri_min = MIN(tri_vertices(1,1,i),tri_vertices(1,2,i),&
    tri_vertices(1,3,i))
  y_tri_max = MAX(tri_vertices(2,1,i),tri_vertices(2,2,i),&
    tri_vertices(2,3,i))
  y_tri_min = MIN(tri_vertices(2,1,i),tri_vertices(2,2,i),&
    tri_vertices(2,3,i))
  z_tri_max = MAX(tri_vertices(3,1,i),tri_vertices(3,2,i),&
    tri_vertices(3,3,i))
  z_tri_min = MIN(tri_vertices(3,1,i),tri_vertices(3,2,i),&
    tri_vertices(3,3,i))
! Determines whether the ray passes through the volume of the triangle
 1000 CALL box_test(ray_end1,ray_end2,x_tri_max,&
        x_tri_min,y_tri_max,y_tri_min,z_tri_max,z_tri_min,&
        outside)
  IF (outside == 0) THEN
!
! Start of the point-in-volume problem
!
    dot = normal(1,i)*tri_vertices(1,1,i)+normal(2,i)*&
      tri_vertices(2,1,i)+normal(3,i)*tri_vertices(3,1,i)
    upper = dot - (ray_end1(1)*normal(1,i)+ray_end1(2)*&
      normal(2,i)+ray_end1(3)*normal(3,i))
    lower = (ray_end2(1)-ray_end1(1))*normal(1,i)+&
      (ray_end2(2)-ray_end1(2))*normal(2,i)+&
      (ray_end2(3)-ray_end1(3))*normal(3,i)
! Take care of degenerate cases
!
! if lower is 0, the ray is parallel to the plane
!
    IF (lower >= minus_epsilon .AND. lower <= plus_epsilon) THEN
!
! if upper is 0, the point is on the plane
! if both are 0, the point is on the plane and ray is parallel to the plane
! this should not happen in ibb as the node should be flagged as a solid before getting here
!
        !write(*,*) 'grumble grubmle',a,b,c,t,lvl
      IF (upper >= minus_epsilon .AND. upper <= plus_epsilon) THEN
!
! Cycling because
!
          !write(*,*) 'weird snooting noise'
        CYCLE triangle_loop

      ELSE
!
! Parallel, but not on the plane means the ray doesn't pass through the triangle.
!
!          write(*,*) 'special ibb case, parallel ray',upper,lower,dot,normal(1:3,i),ray_end1,ray_end2
!          dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
        CYCLE triangle_loop
      END IF
    ELSE
! If you won't get a divide by zero error, can continue with the rest of
! the analysis
      planar_coeff = upper/lower
      if (abs(planar_coeff-1.0D0) < closest_diff_hi) then
        closest_diff_hi = planar_coeff
      end if
      if (abs(planar_coeff) < closest_diff_lo) then
        closest_diff_lo = planar_coeff
      end if

      !closest_diff =
!
!        IF ((planar_coeff > minus_epsilon).AND.(planar_coeff < 1.0D0+plus_epsilon)) THEN
      IF ((planar_coeff > 0.0D0).AND.(planar_coeff < 1.0D0+plus_epsilon)) THEN
        ray_crossing(1) = ray_end1(1)+planar_coeff*(ray_end2(1)-ray_end1(1))
        ray_crossing(2) = ray_end1(2)+planar_coeff*(ray_end2(2)-ray_end1(2))
        ray_crossing(3) = ray_end1(3)+planar_coeff*(ray_end2(3)-ray_end1(3))

        planar_coeff_counter = planar_coeff_counter+1
          !write(*,*) 'crosses the plane',planar_coeff
!                      WRITE(11,*) planar_coeff
!
! tri_volume is the workhorse, it determines if the ray passes through the
! triangle
!                      tri_vol_calls = tri_vol_calls +1
        IF (upper >= lower -very_small_number .AND. upper <= lower + very_small_number) THEN
!
! End point lies directly on the plane
!
          write(*,*) 'calling tri_area'
          CALL tri_area(ray_end1,tri_vertices(1,1,i),tri_vertices(2,1,i&
            ),tri_vertices(3,1,i),tri_vertices(1,2,i),tri_vertices(2,2,i),&
            tri_vertices(3,2,i),tri_vertices(1,3,i),tri_vertices(2,3,i),&
            tri_vertices(3,3,i),triangle_area)
!
! Deal with the output
!
          IF (triangle_area == 1) THEN
!
! Normal function, all areas either positive or negative
!
!              write(*,*) 'node on plane ibb values',planar_coeff,upper,lower,dot,normal(1:3,i)
            dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
            return
          ELSE IF (triangle_area == -1) THEN
!
! End point lies directly on edge
!
!              write(*,*) 'ibb, node on edge',planar_coeff,upper,lower,ray_end1,ray_crossing
            dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
            return
          ELSE IF (triangle_area == -2) THEN
!
! End point lies directly on vertex
!
!              write(*,*) 'ibb, node on vertex',planar_coeff,upper,lower,ray_end1,ray_crossing
            dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
            return
!              CYCLE dir_loop_2
          END IF
        ELSE
          CALL tri_volume(ray_crossing, ray_end1,&
            tri_vertices(1,1,i),tri_vertices(2,1,i),tri_vertices(3,1,i),&
            tri_vertices(1,2,i),tri_vertices(2,2,i),tri_vertices(3,2,i),&
            tri_vertices(1,3,i),tri_vertices(2,3,i),&
            tri_vertices(3,3,i),volume_tri)
!
! Deal with the output of the call to tri_volume
!
          IF (volume_tri == 1) THEN
!
! Normal operation.
!
!              write(*,*) 'normal ibb values',planar_coeff,upper,lower,dot,normal(1:3,i)
            dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)

            return
!              CYCLE dir_loop_2
          ELSE IF (volume_tri == -1) THEN
!
! ray intersects an edge
!
              !write(*,*) 'ray intersects edge',planar_coeff,upper,lower,dot,normal(1:3,i)
          dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
          return
        ELSE IF (volume_tri == -2) THEN
!
! Ray intersects a vertex
!
!              write(*,*) 'ray intersects vertex',planar_coeff,upper,lower,dot,normal(1:3,i)
          dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
          return
        ELSE IF (volume_tri == -3 ) THEN
!
! Ray on a vertex
!
          dir_state = crossing_point_and_dist(ray_crossing,ray_end1,t,rdx(lvl),lvl)
          return
        END IF
      END IF
!
! Ray parallel to plane, need to go to tri_area
!
      ELSE IF (upper >= -very_small_number .AND. upper <= very_small_number) THEN
! Call the triangle area code if the ray lies completely on the plane
                    WRITE(*,*) 'Oook',t,a,b,c,ray_end1,ray_end2,upper,lower
        CALL surface_point_adjuster(ray_end1,ray_end2,t,source,point_moved)
!                      WRITE(*,*) 'Oook',t,a,b,c,ray_end1,ray_end2,upper,lower
        GO TO 1000

      END IF


    END IF !valid planar coeff

  END IF !inside the bounding box

END DO triangle_loop
!write(*,*) 'ibb setup end reached without finding an intersecting triangle',ray_end1,&
!  ray_end2,planar_coeff_counter,closest_diff_lo,closest_diff_hi,a,b,c,t,lvl,self
dir_state = lvl
!            END DO dir_loop_2

!          END DO x_loop
!        END DO y_loop
!      END DO z_loop

!WRITE(11,102)

!CALL wall_lattice

 101  FORMAT('Beginning setup of interpolative bounceback')
 102  FORMAT('Calling wall lattice to create a projected surface mesh')
 103  FORMAT('Surface node idices range from ',I9,' to ',I9)
END SUBROUTINE

!mod_total_wall_nodes = -1000-num_solid_dirs

!WRITE(11,101)
!WRITE(11,103) mod_total_wall_nodes,-1001
!WRITE(*,*) 'Beginning interpolative bb setup.'

!ALLOCATE(wall_nodes(3,mod_total_wall_nodes:-1001))
!ALLOCATE(wall_dist(mod_total_wall_nodes:-1001))
!ALLOCATE(associated_tri(mod_total_wall_nodes:-1001))

!wall_num_counter = -1001
!wall_diag_counter = -1001-(num_solid_dirs-num_card_solid_dirs)

!z_loop: DO c = 1,num_z
!  y_loop: DO b = 1,num_y
!    x_loop: DO a = 1,num_x

!IF (state(1,a,b,c) /=0 ) THEN
!              WRITE(*,*) 'Hokey smokes Bullwinkle.'
!        CYCLE x_loop
!      END IF
!
!      dir_loop_2: DO t = 2,dir
!
! Only need to deal with bounceback nodes.
!
!              IF (state(t,a,b,c) /= -1) THEN

!                CYCLE dir_loop_2
!              END IF

!ray_end1(1)=grid(1,a,b,c)
!ray_end1(2)=grid(2,a,b,c)
!ray_end1(3)=grid(3,a,b,c)
