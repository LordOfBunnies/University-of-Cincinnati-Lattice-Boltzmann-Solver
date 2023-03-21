subroutine state_directions_shifted_cgns_lvl(lvl)
!
!
!
!
!
!
use grid_data
use amrex_amr_module
use amrex_base_module
use geom_data
use linkwise
use nml_inlet_outlet
!use nml_bounding_box
use amr_info_holder!, only:mfstate,state,box_lvl,boi_lo,boi_hi
use amr_processes, only: get_real_coords
use precise
IMPLICIT NONE
INTEGER :: a,b,c,i,t,triangle_area,volume_tri,bot_bound(3),top_bound(3)
REAL(KIND=dp) :: x_tri_max,x_tri_min,y_tri_max,y_tri_min,z_tri_max,bix_lo(3),bix_hi(3)
REAL(KIND=dp) :: z_tri_min,upper,lower,dot,planar_coeff,ray_length
INTEGER :: cross_counter,edge_counter,outside,final_counter,loop_counter
INTEGER :: fluid_target,solid_target,null_target,source,fine,lvl
LOGICAL :: closeness,point_moved,dir_assigned
!LOGICAL,ALLOCATABLE :: skippy(:)
REAL(KIND=dp) :: ray_end1(3),ray_end2(3),ray_crossing(3)

type(amrex_mfiter) :: mfi
type(amrex_box) :: bix

null_target = 0
fluid_target = 1
solid_target = 0
source = 2


!do lvl = 0,fine
call amrex_mfiter_build(mfi,mfstate(lvl), tiling=.false.)

  do while (mfi%next())
  state => mfstate(lvl)%dataptr(mfi)
  bix = mfi%tilebox()
  !bix_lo = lbound(bix)
  !bix_hi = ubound(bix)
  bot_bound = bix%lo - nghosts_mid
  top_bound = bix%hi + nghosts_mid

  z_loop: DO c = bix%lo(3)-nghosts_mid,bix%hi(3)+nghosts_mid
    y_loop: DO b = bix%lo(2)-nghosts_mid,bix%hi(2)+nghosts_mid
      x_loop: DO a = bix%lo(1)-nghosts_mid,bix%hi(1)+nghosts_mid

      point_moved = .FALSE.

      IF (state(a,b,c,1) == -1001 .OR. state(a,b,c,1) == -666) THEN
!              WRITE(*,*) 'Hokey smokes Bullwinkle.'
        CYCLE x_loop
      END IF



      dir_loop_2: DO t = 2,dir

        if (a+cx(t) < bot_bound(1) .OR. a+cx(t) > top_bound(1) .OR. &
            b + cy(t) < bot_bound(2) .OR. b+cy(t) > top_bound(2) .OR. &
            c + cz(t) < bot_bound(3) .OR. c+cz(t)  > top_bound(3)) then
          state(a,b,c,t) = -666
          cycle

!          ray_end1(1)=grid(1,a,b,c)  ! - DBLE(cx_27(t))*dx/10.0D0
!          ray_end1(2)=grid(2,a,b,c)  ! - DBLE(cy_27(t))*dx/10.0D0
!          ray_end1(3)=grid(3,a,b,c)  ! - DBLE(cz_27(t))*dx/10.0D0

!                CONTINUE
        else if (a-cx(t) < 0 .OR. a-cx(t) > amrex_geom(lvl)%domain%hi(1)+1 .OR. &
            b - cy(t) < 0 .OR. b-cy(t) > amrex_geom(lvl)%domain%hi(2)+1 .OR. &
            c - cz(t) < 0 .OR. c-cz(t)  > amrex_geom(lvl)%domain%hi(3)+1 ) THEN

!            write(*,*) 'boopin de snoot!',a,b,c
          ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))
          ray_end2(1) = ray_end1(1) - cx(t)*rdx(lvl)
          ray_end2(2) = ray_end1(2) - cy(t)*rdx(lvl)
          ray_end2(3) = ray_end1(3) - cz(t)*rdx(lvl)
        else if (a+cx(t) < 0 .OR. a+cx(t) > amrex_geom(lvl)%domain%hi(1)+1 .OR. &
            b + cy(t) < 0 .OR. b+cy(t) > amrex_geom(lvl)%domain%hi(2)+1 .OR. &
            c + cz(t) < 0 .OR. c+cz(t)  > amrex_geom(lvl)%domain%hi(3)+1 ) THEN

!         skippy(t) = .FALSE.
          !ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))

!        else if (  state(a+cx(t),b+cy(t),c+cz(t),1) == -1001 .and. shifted) then
          state(a,b,c,t) = -666
          cycle
!        else if (state(a+cx(t),b+cy(t),c+cz(t),1) < 0) then
!
!          if (state(a-cx(t),b-cy(t),c-cz(t),1) >= 0) then
!            state(a,b,c,t) = -666
!            cycle
!          else if (state(a-cx(t),b-cy(t),c-cz(t),1) < 0) then
!            ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))
!            ray_end2(1) = ray_end1(1) - cx(t)*rdx(lvl)
!            ray_end2(2) = ray_end1(2) - cy(t)*rdx(lvl)
!            ray_end2(3) = ray_end1(3) - cz(t)*rdx(lvl)
!          end if

        else if (state(a-cx(t),b-cy(t),c-cz(t),1) < 0) then
          ray_end1 = get_real_coords(amrex_geom(lvl),(/a,b,c/))
          ray_end2(1) = ray_end1(1) - cx(t)*rdx(lvl)
          ray_end2(2) = ray_end1(2) - cy(t)*rdx(lvl)
          ray_end2(3) = ray_end1(3) - cz(t)*rdx(lvl)
          !state(a,b,c,
        else if (state(a,b,c,1) /= state(a+cx(t),b+cy(t),c+cz(t),1)) then
          state(a,b,c,t) = -666
          cycle dir_loop_2
        ELSE
          state(a,b,c,t) = state(a+cx(t),b+cy(t),c+cz(t),1)
          CYCLE dir_loop_2
        END IF

!        if (state(a+cx(t),b+cy(t),c+cz(t),1) < -1000 .and. .not. shifted) then
!          ray_end2(1) = ray_end1(1) + cx(t)*rdx(lvl)
!          ray_end2(2) = ray_end1(2) + cy(t)*rdx(lvl)
!          ray_end2(3) = ray_end1(3) + cz(t)*rdx(lvl)
!        if (state(a-cx(t),b-cy(t),c-cz(t),1) < -1000) then
!          ray_end2(1) = ray_end1(1) - cx(t)*rdx(lvl)
!          ray_end2(2) = ray_end1(2) - cy(t)*rdx(lvl)
!          ray_end2(3) = ray_end1(3) - cz(t)*rdx(lvl)
!        else
!          ray_end2(1) = ray_end1(1) + cx(t)*rdx(lvl)
!          ray_end2(2) = ray_end1(2) + cy(t)*rdx(lvl)
!          ray_end2(3) = ray_end1(3) + cz(t)*rdx(lvl)
!        end if

        dir_assigned = .FALSE.

!              WRITE(*,*) 'Entering triangle loop.'

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

          loop_counter = 0
! Determines whether the ray passes through the volume of the triangle
 1000  CALL box_test(ray_end1,ray_end2,x_tri_max,&
         x_tri_min,y_tri_max,y_tri_min,z_tri_max,z_tri_min,&
         outside)

       loop_counter = loop_counter+1

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

            IF (lower >= minus_epsilon .AND. lower <= plus_epsilon) THEN
!               WRITE(11,*) dot,lower,upper
              IF (upper >= minus_epsilon .AND. upper <= plus_epsilon) THEN
! tri_area is if the ray is on the plane of the triangle
!                      tri_area_calls = tri_area_calls +1
                !WRITE(*,*) 'Spook'
                ray_end2 = cgns_parallel_adjuster(lvl,loop_counter)
                go to 1000
                !cycle triangle_loop

              ELSE
                state(a,b,c,t) = tri_meaning(i)
                dir_assigned = .TRUE.
                CYCLE dir_loop_2
!                      planar_coeff = upper/lower
             END IF
           ELSE
! If you won't get a divide by zero error, can continue with the rest of
! the analysis
              planar_coeff = upper/lower
              ray_crossing(1) = ray_end2(1)+planar_coeff*(&
                ray_end1(1)-ray_end2(1))
              ray_crossing(2) = ray_end2(2)+planar_coeff*(&
                ray_end1(2)-ray_end2(2))
              ray_crossing(3) = ray_end2(3)+planar_coeff*(&
                ray_end1(3)-ray_end2(3))
!
!
              IF ((planar_coeff > 0.0D0).AND.(planar_coeff < 1.0D0)) THEN
!
! tri_volume is the workhorse, it determines if the ray passes through the
! triangle
                IF (upper >= lower - plus_epsilon .AND. upper &
                <= lower + plus_epsilon) THEN
! Same issue as earlier
                      !WRITE(*,*) 'Calling tri_area',upper,lower
!                      tri_area_calls = tri_area_calls +1
                  CALL tri_area(ray_end1,tri_vertices(1,1,i),tri_vertices(2,1,i&
                    ),tri_vertices(3,1,i),tri_vertices(1,2,i),tri_vertices(2,2,i),&
                    tri_vertices(3,2,i),tri_vertices(1,3,i),tri_vertices(2,3,i),&
                    tri_vertices(3,3,i),triangle_area)
!
! Deal with the output
!
                  IF (triangle_area == 1) THEN
                    cross_counter = cross_counter +1
                    if  (tri_meaning(i) == -1001) then
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                    end if
                    cycle dir_loop_2
!
                  ELSE IF (triangle_area == -1) THEN
                    edge_counter = edge_counter + 1
                    if  (tri_meaning(i) == -1001) then
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                    end if
                    cycle dir_loop_2
!
                  ELSE IF (triangle_area == -2) THEN
                    WRITE(11,*) 'Degenerate ray, ray on vertex.'
                    !WRITE(*,*) 'Degenerate ray, ray on vertex.',a,b,c
!                          cross_counter = cross_counter +1
                    if (tri_meaning(i) == -1001) then
                      !call inter_bb_setup(state(a,b,c,t),state(a,b,c,1),t,a,b,c,lvl)
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                      cycle dir_loop_2
                      !state(a,b,c,t) = tri_meaning(i)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                      cycle dir_loop_2
                    end if
                    !state(a,b,c,t) = tri_meaning(i)
                    dir_assigned = .TRUE.
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
                    if  (tri_meaning(i) == -1001) then
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                    else
                      state(a,b,c,t) = tri_meaning(i)
!                      if (state(a,b,c,t) <= -100 .and. state(a,b,c,t) >= -199) then
!                        write(*,*) 'JIGGLE JIGGLE'
!                      end if
                    end if
                    cycle dir_loop_2
                  ELSE IF (volume_tri == -1) THEN
                    edge_counter = edge_counter + 1
                    if  (tri_meaning(i) == -1001) then
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                    else
                      state(a,b,c,t) = tri_meaning(i)
!                      if (state(a,b,c,t) <= -100 .and. state(a,b,c,t) >= -199) then
!                        write(*,*) 'JIGGLE JIGGLE',state(a,b,c,t),a,b,c,t
!                      end if
                    end if
                    cycle dir_loop_2
                  ELSE IF (volume_tri == -2) THEN
                    !state(a,b,c,t) = tri_meaning(i)
                    if (tri_meaning(i) == -1001) then
                      !call inter_bb_setup(state(a,b,c,t),state(a,b,c,1),t,a,b,c,lvl)
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                      cycle dir_loop_2
                      !state(a,b,c,t) = tri_meaning(i)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                      cycle dir_loop_2
                    end if
                    dir_assigned = .TRUE.
                    CYCLE dir_loop_2
                  ELSE IF (volume_tri == -3 ) THEN
!                        WRITE(11,*) 'Degenerate ray, node on vertex.'

                    if (tri_meaning(i) == -1001) then
                      !call inter_bb_setup(state(a,b,c,t),state(a,b,c,1),t,a,b,c,lvl)
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                      cycle dir_loop_2
                    else
                      state(a,b,c,t) = tri_meaning(i)
                      cycle dir_loop_2
                    end if
                    dir_assigned = .TRUE.
                    CYCLE dir_loop_2
!                        cross_counter = cross_counter + 1
                 END IF
               END IF
!
! Ray parallel to plane, need to go to tri_area
!
              ELSE IF (upper >= minus_epsilon .AND. upper <= plus_epsilon ) THEN
! Call the triangle area code if the ray lies completely on the plane
!                      tri_area_calls = tri_area_calls +1
!                WRITE(*,*) 'Oook',t,a,b,c,ray_end1,ray_end2,upper,lower
!                CALL surface_point_adjuster(ray_end1,ray_end2,normal(1:3,i),&
!                  t,source,point_moved,lvl)
!                WRITE(*,*) 'Oook',t,a,b,c,ray_end1,ray_end2,upper,lower
!!                      CALL error_out
!                GO TO 1000

                if (loop_counter <= 5) then
                  WRITE(*,*) 'Oook',t,a,b,c,ray_end1,ray_end2,upper,lower

                  CALL tri_area(ray_end1,tri_vertices(1,1,i),tri_vertices(2,1,i&
                    ),tri_vertices(3,1,i),tri_vertices(1,2,i),tri_vertices(2,2,i),&
                    tri_vertices(3,2,i),tri_vertices(1,3,i),tri_vertices(2,3,i),&
                    tri_vertices(3,3,i),triangle_area)
!
! Deal with the output
!
                  IF (triangle_area == 1) THEN
                    cross_counter = cross_counter +1
                    if  (tri_meaning(i) == -1001) then
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                    end if
                    cycle dir_loop_2
!
                  ELSE IF (triangle_area == -1) THEN
                    edge_counter = edge_counter + 1
                    if  (tri_meaning(i) == -1001) then
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                    end if
                    cycle dir_loop_2
!
                  ELSE IF (triangle_area == -2) THEN
                    WRITE(11,*) 'Degenerate ray, ray on vertex.'
                    WRITE(*,*) 'Degenerate ray, ray on vertex.',a,b,c
!                          cross_counter = cross_counter +1
                    if (tri_meaning(i) == -1001) then
                      !call inter_bb_setup(state(a,b,c,t),state(a,b,c,1),t,a,b,c,lvl)
                      state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                        rdx(lvl),lvl)
                      cycle dir_loop_2
                      !state(a,b,c,t) = tri_meaning(i)
                    else
                      state(a,b,c,t) = tri_meaning(i)
                      cycle dir_loop_2
                    end if
                    !state(a,b,c,t) = tri_meaning(i)
                    dir_assigned = .TRUE.
                  END IF

                  WRITE(*,*) 'Oook',t,a,b,c,ray_end1,ray_end2,upper,lower

                  GO TO 1000
                else if (loop_counter <=10) then
                  write(*,*) 'Eek',t,a,b,c,ray_end1,ray_end2,upper,lower
                  ray_end1 = cgns_node_push(lvl,t)

                  WRITE(*,*) 'Aack',t,a,b,c,ray_end1,ray_end2,upper,lower

                  go to 1000
                else
                  state(a,b,c,t) = tri_meaning(i)
                  cycle triangle_loop
                end if
              END IF
            END IF

          END IF
!
! This is different in regard to the final outcome. Because the ray is from
! one node to the next, it should only cross one or two triangles (worst case,
! pass through a vertex)
!
          IF (cross_counter == 1 .OR. edge_counter == 1) THEN
            !state(a,b,c,t) = tri_meaning(i)
!            if (tri_meaning(i) == -1001) then
            !call inter_bb_setup(state(a,b,c,t),state(a,b,c,1),t,a,b,c,lvl)
            if (tri_meaning(i) == -1001) then
                      !call inter_bb_setup(state(a,b,c,t),state(a,b,c,1),t,a,b,c,lvl)
              state(a,b,c,t) = crossing_point_and_dist(ray_crossing,ray_end1,t,&
                rdx(lvl),lvl)
              cycle dir_loop_2
            else
              state(a,b,c,t) = tri_meaning(i)
              cycle dir_loop_2
            end if
                     !state(a,b,c,t) = tri_meaning(i)
!            else
!              state(a,b,c,t) = tri_meaning(i)
!            end if
            dir_assigned = .TRUE.
!                  WRITE(*,*) 'Thingy assigned to boundary node dir',tri_meaning(i)
            CYCLE dir_loop_2
          END IF

          IF (.NOT. dir_assigned) THEN
            state(a,b,c,t) = -666

          END IF
        END DO triangle_loop

      END DO dir_loop_2

    END DO x_loop
  END DO y_loop
END DO z_loop

end do

call amrex_mfiter_destroy(mfi)

!end do

 505  FORMAT ('Starting the process of finding boundary conditions.')

END SUBROUTINE
