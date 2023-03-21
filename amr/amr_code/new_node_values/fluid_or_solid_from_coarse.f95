subroutine fluid_or_solid_from_coarse(lvl)
!
! Same as fluid_solid before, but needs to be done on a new level.  A shortcut
!   needs to be found as imultifabs cannot be copied or accessed together.
!
!
! Called by: from_coarse (amr_processes)
! Calls:
! External calls: amrex_mfiter_build
!
use constants
use precise
use geom_data
use startup
use linkwise
use nml_inlet_outlet
use grid_data
use amr_info_holder, only: state,mfstate,nghosts_state
use amr_processes
use amrex_base_module
use amrex_amr_module
use mpi
implicit none
integer,intent(in) :: lvl
real(KIND=dp) :: x_tri_max,x_tri_min,y_tri_max,y_tri_min,z_tri_max
real(KIND=dp) :: z_tri_min,upper,lower,dot,planar_coeff,ray_length
real(KIND=dp) :: x_grid_delta,y_grid_delta,z_grid_delta
integer :: a,b,c,i,triangle_area,volume_tri,fluid_node_counter!,bix_lo(3),bix_hi(3)
integer :: cross_counter,edge_counter,outside,final_counter,bot_bound(4),top_bound(4)
integer :: fluid_target,solid_target,null_target,useless,source,fine,mx_lv
integer :: ier,total_fluid_nodes, solid_counter
!integer :: tri_vol_calls,tri_area_calls,big_bt_calls,lil_bt_calls
logical :: closeness,point_moved,first_corner_fail,exterior
real(KIND=dp) :: ray_end1(3),ray_end2(3),ray_crossing(3),dummy_ray(3)

type(amrex_mfiter) :: mfi
type(amrex_box) :: mux
!
!
!
!
!
x_grid_delta = amrex_probhi(1)-amrex_problo(1)
y_grid_delta = amrex_probhi(2)-amrex_problo(2)
z_grid_delta = amrex_probhi(3)-amrex_problo(3)
fluid_node_counter = 0
ray_length = ABS(1000.0D0*MAX(x_grid_delta,y_grid_delta,&
  z_grid_delta))

IF (cgns_start .OR. mixed_start .or. .not. external_flow) THEN
  fluid_target = 1
  solid_target = 0
ELSE
  fluid_target = 0
  solid_target = 1
END IF
!write(*,*) 'hi, Im ',self,x_grid_delta,fluid_node_counter,ray_length
!
!
!
call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

do while (mfi%next())
  mux = mfi%validbox()

  state => mfstate(lvl)%dataptr(mfi)

!  write(*,*) 'um, what?',self,mux%lo,mux%hi,lvl

  do c = mux%lo(3)-nghosts_state,mux%hi(3)+nghosts_state
    do b = mux%lo(2)-nghosts_state,mux%hi(2)+nghosts_state
      x_loop: do a = mux%lo(1)-nghosts_state,mux%hi(1)+nghosts_state
          if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
              c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
              b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then
            if (external_flow) then
              state(a,b,c,1) = -666
              exterior = .true.
              cycle
            else
              if (a < amrex_geom(lvl)%domain%lo(1) .and. inlet_side(1)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (b < amrex_geom(lvl)%domain%lo(2) .and. inlet_side(2)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (c < amrex_geom(lvl)%domain%lo(3) .and. inlet_side(3)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (a > amrex_geom(lvl)%domain%hi(1) .and. outlet_side(4)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (b > amrex_geom(lvl)%domain%hi(2) .and. outlet_side(5)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (c > amrex_geom(lvl)%domain%hi(3) .and. outlet_side(6)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (a < amrex_geom(lvl)%domain%lo(1) .and. outlet_side(1)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (b < amrex_geom(lvl)%domain%lo(2) .and. outlet_side(2)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (c < amrex_geom(lvl)%domain%lo(3) .and. outlet_side(3)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (a > amrex_geom(lvl)%domain%hi(1) .and. inlet_side(4)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (b > amrex_geom(lvl)%domain%hi(2) .and. inlet_side(5)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else if (c > amrex_geom(lvl)%domain%hi(3) .and. inlet_side(6)) then
                state(a,b,c,1) = -666
                exterior = .true.
                cycle
              else
                state(a,b,c,1) = -1001
                exterior = .true.
                cycle
              end if
            end if
          else
            exterior = .false.
          end if

!            state(a,b,c,1) = -666
!            exterior = .true.
!            cycle
!          else
!            exterior = .false.
!          end if
!          if (a == 301 .and. b == 221 .and. c == -3) then
!            write(*,*) ' ok, what the actual fluff',state(a,b,c,1),a,b,c,amrex_geom(lvl)%domain%lo,&
!              amrex_geom(lvl)%domain%hi
!          end if


          ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))
!      ray_end1(1)=grid(a,b,c,1)
!      ray_end1(2)=grid(2,a,b,c)
!      ray_end1(3)=grid(3,a,b,c)
! ray_casting generates a directly away from the geometric center
          CALL ray_casting_away(ray_length,ray_end1,ray_end2)
!          write(*,*) 'ray info',ray_length,ray_end1,ray_end2
! box_test tests whether the ray might pass through the geometry
          first_corner_fail = .false.

!          if (self == 2 .or. self == 3) then
!            write(*,*) 'f/s proc info',ray_end1,ray_end2,ray_length,mux%lo,mux%hi,a,b,c,self
!          end if

 2000     CALL box_test(ray_end1,ray_end2,x_geom_max,&
            x_geom_min,y_geom_max,y_geom_min,z_geom_max,z_geom_min,&
            outside)
!          if (a == 124 .and. b == 64 .and. c == 64) then
!            write(*,*) 'should be solid',ray_end1,ray_end2,outside,a,b,c,self,lvl
!          end if

          IF (outside == 1) THEN
            state(a,b,c,1) = lvl
            fluid_node_counter = fluid_node_counter + 1
            cycle
          ELSE
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
! Take care of degenerate cases

              IF (lower >= minus_epsilon .AND. lower <=  plus_epsilon) THEN
!                    WRITE(11,*) dot,lower,upper
                IF (upper >= minus_epsilon .AND. upper <=  plus_epsilon) THEN
! tri_area is if the ray is on the plane of the triangle
!
              !CALL ray_parallel_adjuster(ray_end2,normal(1,i),&
              !  normal(2,i),normal(3,i))
!                  ray_end2 = ray_parallel_adjuster(ray_end2(1),ray_end2(2),&
!                    ray_end2(3),normal(1,i),normal(2,i),normal(3,i),link_dist(2,lvl))
!!                  write(*,*) 'wrecker',self
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
!                  write(*,*) 'blecker',upper,lower,ray_end1,ray_end2
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
!
                IF ((planar_coeff > 0.0D0 + minus_epsilon).AND.(planar_coeff < 1.0D0+ plus_epsilon)) THEN
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
!                      WRITE(*,*) 'Calling tri_area',upper,lower
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
                    ELSE IF (triangle_area == -2) THEN
                      WRITE(11,*) 'Degenerate ray, ray on vertex.'
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
                      !WRITE(*,*) 'Ray passes through vertex, checking ray.'
                      CALL closeness_checker(ray_end1,&
                        tri_vertices(1,1,i),tri_vertices(2,1,i),tri_vertices(3,1,i),&
                        tri_vertices(1,2,i),tri_vertices(2,2,i),tri_vertices(3,2,i),&
                        tri_vertices(1,3,i),tri_vertices(2,3,i),&
                        tri_vertices(3,3,i),closeness)
                      IF (closeness) THEN
                        state(a,b,c,1) = lvl
                        CYCLE x_loop
                      ELSE
                        WRITE(11,*) 'Ray far from vertex, regenerating ray.'
!                  WRITE(*,*) 'Ray far from vertex, regenerating ray.',a,b,c,self,lvl
!                WRITE(*,*) 'Ray far from vertex, regenerating ray.'
                        CALL ray_casting(ray_length,ray_end2)
!                        write(*,*) 'blargh',self
                        GO TO 1000 !Don't judge me
                      END IF
!                        WRITE(*,*) 'Ing!',volume_tri
!                        GO TO 100 !Don't judge me
                    ELSE IF (volume_tri == -3 ) THEN
                      WRITE(11,*) 'Degenerate ray, node on vertex.'
  !                WRITE(*,*) 'Degenerate ray, node on vertex.',a,b,c
                      state(a,b,c,1) = -1001
                      solid_counter = solid_counter + 1
                      !WRITE(*,*) 'SOLIDS!',a,b,c,lvl,self
                      CYCLE x_loop
!              cross_counter = cross_counter + 1
                    END IF
                  END IF
! Ray parallel to plane, need to go to tri_area
!
                ELSE IF (upper >= minus_epsilon .AND. upper <= plus_epsilon) THEN
!
!                  write(*,*) 'Adjusting point near surface.',a,b,c,self
                  CALL surface_point_adjuster(ray_end1,ray_end2,&
                    useless,source,point_moved,lvl)
!                  write(*,*) 'doodle snoot',self
                  GO TO 1000
!
                END IF
              END IF
            END IF
          END DO triangle_loop

          if (edge_counter >= 2) then
!            write(*,*) 'harbormaster',edge_counter,self
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
            IF (cross_counter == 0 .AND. edge_counter == 0) THEN
              state(a,b,c,1) = -666
            END IF

            IF (edge_counter == 1 .AND. MODULO(cross_counter,2)==0) THEN
                edge_counter = 0
                !WRITE(*,*) 'Single edge encountered.'
                cross_counter = cross_counter+1
            END IF
!
!
            edge_counter = edge_counter/2
            final_counter = cross_counter + edge_counter
            final_counter = MODULO(final_counter,2)
! Final ruling on if node is fluid or solid
            IF (final_counter == solid_target) THEN
              state(a,b,c,1) = -1001
              solid_counter = solid_counter + 1
!
            ELSE
!
! This is to ensure that fluid nodes that are outside the domain, are made null to prevent whoopsies.
!   It should not affect state_directions.
!
              if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
                c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
                b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then

                state(a,b,c,1) = -666
              else
                state(a,b,c,1) = lvl
                fluid_node_counter = fluid_node_counter + 1
              end if
            END IF

          END IF

!        if (a == 301 .and. b == 221 .and. c == -3) then
!          write(*,*) 'exterior node values',state(a,b,c,1),a,b,c
!        end if
!
!        if (a == 132 .and. b == 60 .and. c == 67) then
!          dummy_ray = get_real_coords(amrex_geom(lvl),(/132,60,67/))
!          write(*,*) 'fluid/solid value',dummy_ray,cross_counter,final_counter,&
!            state(a,b,c,1),self,lvl,a,b,c
!!          write(*,*) 'triangle limits',x_geom_max,&
!!            x_geom_min,y_geom_max,y_geom_min,z_geom_max,z_geom_min,a,b,c,self,lvl
!        end if

      end do x_loop
    end do
  end do




end do
write(*,*) 'number of fluid nodes after regrid on lvl',fluid_node_counter,lvl,self
call mpi_barrier(commune,ier)
!write(*,*) 'dingleberry',self

CALL fluid_solid_cleanup_lvl(lvl)


end subroutine
