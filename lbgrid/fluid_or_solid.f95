subroutine fluid_or_solid(fine)
!
! Determines whether a node is fluid or solid
!
! Called by: grid_gen
! Calls: ray_casting,box_test,tri_area,tri_volume,error_out,inter_bb_setup
! External calls: amrex_mfiter_build,amrex_mfiter_destroy
!
!use nml_bounding_box
use constants
use grid_data
use geom_data
use precise
use startup
use linkwise
use nml_inlet_outlet
use amrex_amr_module
use amrex_base_module
use amr_processes, only: self,commune,get_real_coords!,nprocs
use amr_info_holder, only: mfstate, state,nghosts,nghosts_state!,box_lvl
use mpi
IMPLICIT NONE
real(KIND=dp) :: x_tri_max,x_tri_min,y_tri_max,y_tri_min,z_tri_max
real(KIND=dp) :: z_tri_min,upper,lower,dot,planar_coeff,ray_length
real(KIND=dp) :: x_grid_delta,y_grid_delta,z_grid_delta
integer :: a,b,c,i,triangle_area,volume_tri,fluid_node_counter!,bix_lo(3),bix_hi(3)
integer :: cross_counter,edge_counter,outside,final_counter,bot_bound(4),top_bound(4)
integer :: fluid_target,solid_target,null_target,useless,source,fine,mx_lv,lvl
integer :: ier,total_fluid_nodes, solid_counter
!integer :: tri_vol_calls,tri_area_calls,big_bt_calls,lil_bt_calls
logical :: closeness,point_moved,first_corner_fail,exterior
real(KIND=dp) :: ray_end1(3),ray_end2(3),ray_crossing(3)

type(amrex_mfiter) :: mfi
type(amrex_box) :: bix
!
! Slightly different targets for number of crosses for normal startup vs.
! cgns starts and internal flows
!

!write(*,*) 'grid bounds',amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi

call mpi_barrier(commune,ier)

null_target = 0
source = 1
useless = 0
mx_lv = amrex_max_level

IF (cgns_start .OR. mixed_start .or. .not. external_flow) THEN
  fluid_target = 1
  solid_target = 0
ELSE
  fluid_target = 0
  solid_target = 1
END IF
! Set up the maximums
x_grid_delta = amrex_probhi(1)-amrex_problo(1)
y_grid_delta = amrex_probhi(2)-amrex_problo(2)
z_grid_delta = amrex_probhi(3)-amrex_problo(3)

grid_center(1) = (amrex_problo(1)+amrex_probhi(1))/2.0D0
grid_center(2) = (amrex_problo(2)+amrex_probhi(2))/2.0D0
grid_center(3) = (amrex_problo(3)+amrex_probhi(3))/2.0D0

!tri_vol_calls = 0
!tri_area_calls = 0
!big_bt_calls = 0
!lil_bt_calls = 0
!
! First we need to make sure the ray is long enough to go completely out
! of the bounds of the simulation volume.
!
!if (mx_lv > maxval(box_lvl)) then
!  refine_lvl = maxval(box_lvl)
!else
!  refine_lvl = mx_lv
!end if

ray_length = ABS(100.0D0*MAX(x_grid_delta,y_grid_delta,&
  z_grid_delta))
fluid_node_counter = 0

!write(*,*) 'Begininning determination of whether node is fluid or solid.',fine

do lvl = 0,fine
  call amrex_mfiter_build(mfi,mfstate(lvl), tiling=.false.)
  do while (mfi%next())

    state => mfstate(lvl)%dataptr(mfi)
    bix = mfi%validbox()

    solid_counter = 0

!  WRITE(*,*) 'Hokey smokes Bullwinkle.',self,bix%lo,bix%hi,amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi
!  bot_bound = lboun
!  top_bound = ubound(state)
! Loop through all the nodes
    z_loop: DO c = bix%lo(3)-nghosts_state,bix%hi(3)+nghosts_state
      y_loop: DO b = bix%lo(2)-nghosts_state,bix%hi(2)+nghosts_state
        x_loop: DO a = bix%lo(1)-nghosts_state,bix%hi(1)+nghosts_state
          if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
              c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
              b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then
            if (external_flow) then
              !state(a,b,c,1) = -666
              exterior = .true.
              !cycle
            else
              if (a < amrex_geom(lvl)%domain%lo(1) .and. inlet_side(1)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (b < amrex_geom(lvl)%domain%lo(2) .and. inlet_side(2)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (c < amrex_geom(lvl)%domain%lo(3) .and. inlet_side(3)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (a > amrex_geom(lvl)%domain%hi(1) .and. outlet_side(4)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (b > amrex_geom(lvl)%domain%hi(2) .and. outlet_side(5)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (c > amrex_geom(lvl)%domain%hi(3) .and. outlet_side(6)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (a < amrex_geom(lvl)%domain%lo(1) .and. outlet_side(1)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (b < amrex_geom(lvl)%domain%lo(2) .and. outlet_side(2)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (c < amrex_geom(lvl)%domain%lo(3) .and. outlet_side(3)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (a > amrex_geom(lvl)%domain%hi(1) .and. inlet_side(4)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (b > amrex_geom(lvl)%domain%hi(2) .and. inlet_side(5)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else if (c > amrex_geom(lvl)%domain%hi(3) .and. inlet_side(6)) then
                !state(a,b,c,1) = -666
                exterior = .true.
                !cycle
              else
                !state(a,b,c,1) = -1001
                exterior = .true.
                !cycle
!                state(a,b,c,1) = -666
!                exterior = .true.
!                cycle
              end if
            end if
          else
            exterior = .false.
          end if

!        IF (state(a,b,c,1) == -1) THEN

!          CYCLE x_loop
!        END IF
!        write(*,*) self,lvl,a,b,c
          ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))
          !write(*,*) 'ray start',ray_end1,a,b,c
!          write(*,*) 'geom maxes and mins',x_geom_max,&
!            x_geom_min,y_geom_max,y_geom_min,z_geom_max,z_geom_min
!
! ray_casting generates a directly away from the geometric center
!
          CALL ray_casting_away(ray_length,ray_end1,ray_end2)
!
! box_test tests whether the ray might pass through the geometry
!
          first_corner_fail = .false.
 2000     CALL box_test(ray_end1,ray_end2,x_geom_max,&
            x_geom_min,y_geom_max,y_geom_min,z_geom_max,z_geom_min,&
            outside)

          IF (outside == 1 .and. external_flow) THEN
            if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
              c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
              b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then
              state(a,b,c,1) = -666
              cycle
!            else
!              if (external_flow) then
!                state(a,b,c,1) = lvl
!                fluid_node_counter = fluid_node_counter + 1
!              else
!                state(a,b,c,1) = -1001
!              end if
            end if
          end if
!          ELSE
            cross_counter = 0
            edge_counter = 0
            final_counter = 0
            volume_tri = 0
            closeness = .FALSE.
!
! Loop through all the triangles
!
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

!            write(*,*) 'max and min',x_tri_max,&
!              x_tri_min,y_tri_max,y_tri_min,z_tri_max,z_tri_min
!            write(*,*) 'tris',tri_vertices(1:3,1:3,i),i
! Determines whether the ray passes through the volume of the triangle
!                lil_bt_calls = lil_bt_calls +1
 1000       CALL box_test(ray_end1,ray_end2,x_tri_max,&
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
!              write(*,*) 'googly moogly',ray_end1,ray_end2,normal(1:3,i),dot,upper,lower
! Take care of degenerate cases

              IF (lower >= minus_epsilon .AND. lower <=  plus_epsilon) THEN
!                    WRITE(11,*) dot,lower,upper
                IF (upper >= minus_epsilon .AND. upper <=  plus_epsilon) THEN
! tri_area is if the ray is on the plane of the triangle
!
              !CALL ray_parallel_adjuster(ray_end2,normal(1,i),&
              !  normal(2,i),normal(3,i))
                  write(*,*) 'barf barf!',a,b,c,ray_end1,ray_end2,normal(1:3,i),i
!                  ray_end2 = ray_parallel_adjuster(ray_end2(1),ray_end2(2),&
!                    ray_end2(3),normal(1,i),normal(2,i),normal(3,i),link_dist(2,lvl))
!                  edge_counter = 0
!                  cross_counter = 0
       !           call corner_casting(ray_end1,ray_end2,amrex_problo,&
       !             amrex_probhi,ray_length,first_corner_fail)
       !             first_corner_fail = .true.
       !           go to 2000
       !           ray_end1 = different_parallel_adjuster(normal(1,i),normal(2,i),normal(3,i),&
       !               lvl)
                  !call fluff_it_all(ray_end1,ray_end2, normal(1,i),normal(2,i),normal(3,i),&
                  !    lvl)
!                  GO TO 1000
                      CALL tri_area(ray_end1,tri_vertices(1,1,i),tri_vertices(2,1,i&
                        ),tri_vertices(3,1,i),tri_vertices(1,2,i),tri_vertices(2,2,i),&
                        tri_vertices(3,2,i),tri_vertices(1,3,i),tri_vertices(2,3,i),&
                        tri_vertices(3,3,i),triangle_area)
                      IF (triangle_area == 0) THEN
                        cross_counter = cross_counter +1
                      ELSE IF (triangle_area == -1) THEN
                        edge_counter = edge_counter +1
                      END IF
                ELSE
!                      state(a,b,c,1) = 0
                  CALL ray_casting(ray_length,ray_end2)
                  GO TO 1000
!                      planar_coeff = upper/lower
                END IF

              ELSE
!
! If you won't get a divide by zero error, can continue with the rest of
! the analysis
!
                planar_coeff = upper/lower
!
                IF (planar_coeff > 0.0D0 .and. planar_coeff < 1.0D0) THEN
!                IF ((planar_coeff > 0.0D0 + minus_epsilon).AND.(planar_coeff < 1.0D0+ plus_epsilon)) THEN
                  ray_crossing(1) = ray_end1(1)+planar_coeff*(ray_end2(1)-ray_end1(1))
                  ray_crossing(2) = ray_end1(2)+planar_coeff*(ray_end2(2)-ray_end1(2))
                  ray_crossing(3) = ray_end1(3)+planar_coeff*(ray_end2(3)-ray_end1(3))
!
! tri_volume is the workhorse, it determines if the ray passes through the
! triangle
!                      tri_vol_calls = tri_vol_calls +1

                  IF (upper >= lower + minus_epsilon .AND. upper <= lower + plus_epsilon) THEN
! Same issue as earlier
!                       WRITE(*,*) 'Ook'
                      WRITE(*,*) 'Calling tri_area',upper,lower
!                      tri_area_calls = tri_area_calls +1
                    CALL tri_area(ray_end1,tri_vertices(1,1,i),tri_vertices(2,1,i),&
                      tri_vertices(3,1,i),tri_vertices(1,2,i),tri_vertices(2,2,i),&
                      tri_vertices(3,2,i),tri_vertices(1,3,i),tri_vertices(2,3,i),&
                      tri_vertices(3,3,i),triangle_area)
!
! Deal with the output
!
                    IF (triangle_area == 1) THEN
                      cross_counter = cross_counter +1
                    ELSE IF (triangle_area == -1) THEN
                      edge_counter = edge_counter + 1
                      !write(*,*) 'spooky pants'
                    ELSE IF (triangle_area == -2) THEN
                      !WRITE(*,*) 'Degenerate ray, ray on vertex.'
!                      WRITE(*,*) 'Degenerate ray, ray on vertex.',a,b,c
!              cross_counter = cross_counter +1
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
                      cross_counter = cross_counter + 1
!              WRITE(*,*) cross_counter
                    ELSE IF (volume_tri == -1) THEN
                      edge_counter = edge_counter + 1
                      !WRITE(*,*) 'Edge encountered.'
                    ELSE IF (volume_tri == -2) THEN
!              WRITE(11,*) 'Ray passes through vertex,', &
!                ' checking ray.',upper,lower,dot
                     ! WRITE(*,*) 'Ray passes through vertex, checking ray.'
                      CALL closeness_checker(ray_end1,&
                        tri_vertices(1,1,i),tri_vertices(2,1,i),tri_vertices(3,1,i),&
                        tri_vertices(1,2,i),tri_vertices(2,2,i),tri_vertices(3,2,i),&
                        tri_vertices(1,3,i),tri_vertices(2,3,i),&
                        tri_vertices(3,3,i),closeness)
                      IF (closeness) THEN
                        write(*,*) 'close!',a,b,c
                        state(a,b,c,1) = lvl
                        CYCLE x_loop
                      ELSE
!                        WRITE(11,*) 'Ray far from vertex, regenerating ray.'
!                  WRITE(*,*) 'Ray far from vertex, regenerating ray.',a,b,c,self,lvl
                !WRITE(*,*) 'Ray far from vertex, regenerating ray.',a,b,c
                        CALL ray_casting(ray_length,ray_end2)
!                        edge_counter = 0
!                        cross_counter = 0
!                        first_corner_fail = .false.
!                        call corner_casting(ray_end1,ray_end2,amrex_problo,&
!                          amrex_probhi,ray_length,first_corner_fail)
                        GO TO 2000 !Don't judge me
                      END IF
!                        WRITE(*,*) 'Ing!',volume_tri
!                        GO TO 100 !Don't judge me
                    ELSE IF (volume_tri == -3 ) THEN
                      WRITE(11,*) 'Degenerate ray, node on vertex.'
  !                WRITE(*,*) 'Degenerate ray, node on vertex.',a,b,c
                      state(a,b,c,1) = -1001
                      solid_counter = solid_counter + 1
                      WRITE(*,*) 'SOLIDS!',a,b,c,lvl,self
                      write(*,*) 'nodal overlap',ray_end1,&
                        tri_vertices(1,1,i),tri_vertices(2,1,i),tri_vertices(3,1,i),&
                        tri_vertices(1,2,i),tri_vertices(2,2,i),tri_vertices(3,2,i),&
                        tri_vertices(1,3,i),tri_vertices(2,3,i),&
                        tri_vertices(3,3,i)
                      CYCLE x_loop
!              cross_counter = cross_counter + 1
                    END IF
                  END IF
! Ray parallel to plane, need to go to tri_area
!
                ELSE IF (upper >= -very_small_number .AND. upper <=&
                very_small_number ) THEN
!
!                  write(*,*) 'Adjusting point near surface.',a,b,c,self
                  CALL surface_point_adjuster(ray_end1,ray_end2,normal(1:3,i),&
                    useless,source,point_moved,lvl)
                  GO TO 1000
!
                END IF
              END IF
            END IF
          END DO triangle_loop

          if (edge_counter >= 2) then
            !write(*,*) 'noodle?',a,b,c,edge_counter
            edge_counter = 0
            cross_counter = 0
            call corner_casting(ray_end1,ray_end2,amrex_problo,&
              amrex_probhi,ray_length,first_corner_fail)
            first_corner_fail = .true.
            go to 2000
          end if
!
! Once all the triangles are done, we have to compute whether the point
! is inside or outside.  A ray that passes through geometry an even
! number of times is outside the geometry (if it's external flow) and
! inside if it's an odd number. A few exceptions occur like if a ray
! passes through a vertex or is on a vertex.
!
            !write(*,*) 'total crosses and edges',cross_counter,edge_counter
!            IF (cross_counter == 0 .AND. edge_counter == 0) THEN
!              state(a,b,c,1) = -666
!              cycle
!            END IF

!            if (a == 268 .and. b == 20 .and. c == 176) then
!              write(*,*) 'fluid/solid values',edge_counter,cross_counter,ray_end1,ray_end2
!            end if

            IF (modulo(edge_counter,2) == 1 .AND. MODULO(cross_counter,2)==0) THEN
                edge_counter = 0
                !WRITE(*,*) 'Single edge encountered.'
                cross_counter = cross_counter+1
            END IF
!
            edge_counter = edge_counter/2
            final_counter = cross_counter + edge_counter
            final_counter = MODULO(final_counter,2)
! Final ruling on if node is fluid or solid
            IF (final_counter == solid_target) THEN
              state(a,b,c,1) = -1001
              solid_counter = solid_counter + 1
!              if (a == -1 .and. b == 32 .and. c == 32) then
!                write(*,*) 'exterior node values',state(a,b,c,1),edge_counter,cross_counter,&
!                  a,b,c
!              end if
              !WRITE(*,*) 'SOLIDS!',a,b,c,self
!
            ELSE
!
! This is to ensure that fluid nodes that are outside the domain, are made null to prevent whoopsies.
!   It should not affect state_directions.
!
!              if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
!                c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
!                b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then
!
!                !write(*,*) 'boing?',a,b,c
!                state(a,b,c,1) = -666
!              else
                if (exterior) then
                  if (external_flow) then
                    state(a,b,c,1) = -666
                  else
                    if (inlet_side(1) .and. a < amrex_geom(lvl)%domain%lo(1)) then
                      state(a,b,c,1) = -666
                      !write(*,*) 'wiggles'
                    else if (inlet_side(2) .and. b < amrex_geom(lvl)%domain%lo(2)) then
                      state(a,b,c,1) = -666
                    else if (inlet_side(3) .and. c < amrex_geom(lvl)%domain%lo(3)) then
                      state(a,b,c,1) = -666
                    else if (inlet_side(4) .and. a > amrex_geom(lvl)%domain%hi(1)+1) then
                      state(a,b,c,1) = -666
                    else if (inlet_side(5) .and. b > amrex_geom(lvl)%domain%hi(2)+1) then
                      state(a,b,c,1) = -666
                    else if (inlet_side(6) .and. c > amrex_geom(lvl)%domain%hi(3)+1) then
                      state(a,b,c,1) = -666
                    else if (outlet_side(1) .and. a < amrex_geom(lvl)%domain%lo(1)) then
                      state(a,b,c,1) = -666
                    else if (outlet_side(2) .and. b < amrex_geom(lvl)%domain%lo(2)) then
                      state(a,b,c,1) = -666
                    else if (outlet_side(3) .and. c < amrex_geom(lvl)%domain%lo(3)) then
                      state(a,b,c,1) = -666
                    else if (outlet_side(4) .and. a > amrex_geom(lvl)%domain%hi(1)+1) then
                      state(a,b,c,1) = -666
                    else if (outlet_side(5) .and. b > amrex_geom(lvl)%domain%hi(2)+1) then
                      state(a,b,c,1) = -666
                    else if (outlet_side(6) .and. c > amrex_geom(lvl)%domain%hi(3)+1) then
                      state(a,b,c,1) = -666
                    else
                      state(a,b,c,1) = -1001
                      if (a == 268 .and. b == 20 .and. c == 176) then
                        write(*,*) 'exterior node values',state(a,b,c,1),edge_counter,cross_counter,&
                         a,b,c
                       end if
                      !write(*,*) 'snarls'
                    end if
                  end if
                else
                  !if (external_flow) then
                    state(a,b,c,1) = lvl
                    fluid_node_counter = fluid_node_counter + 1
                  !else
                  !  state(a,b,c,1) = -666
                  !end if
                end if
!              end if
            END IF

!          END IF

           if (a == 268 .and. b == 20 .and. c == 176) then
             write(*,*) 'last ditch node values',state(a,b,c,1),a,b,c
           end if

!          if (a == 95 .and. b == 38 .and. c == 0) then
!            write(*,*) 'in-bounds node check',state(a,b,c,1),a,b,c
!          end if
!          if (a == 95 .and. b == 39 .and. c == 1) then
!            write(*,*) 'out of bounds node check',state(a,b,c,1),a,b,c
!          end if

        END DO x_loop
      END DO y_loop
    END DO z_loop
  end do

!call mfstate(lvl)%fill_boundary(amrex_geom(lvl),.true.)

call amrex_mfiter_destroy(mfi)
call mpi_allreduce(fluid_node_counter,total_fluid_nodes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ier)
!write(*,*) 'current fluid nodes',fluid_node_counter,total_fluid_nodes,self,lvl
!
!write(*,*) 'solid nodes ',solid_counter,'on level ',lvl,'for proc',self

call mpi_barrier(commune,ier)
end do
!write(*,*) 'Calling cleanup subroutines.',self
CALL fluid_solid_cleanup(fine)
!CALL fluid_solid_cleanup(fine)
!CALL fluid_solid_cleanup(fine)
!CALL fluid_solid_cleanup(fine)

if (self ==0) WRITE(11,200) total_fluid_nodes
if (self ==0) write(*,200) total_fluid_nodes
!      WRITE(11,*) big_bt_calls,lil_bt_calls,tri_vol_calls,tri_area_calls
      WRITE(11,*) ''
 200  FORMAT ('Total valid fluid nodes = ',I8)
END SUBROUTINE
