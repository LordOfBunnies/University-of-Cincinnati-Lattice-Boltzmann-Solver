      SUBROUTINE ggm_output_setup
!
! Set up the data output for GGM
!
! Called by: links
! Calls:
!
      USE output_data
      USE cgns
      USE grid_data
      USE linkwise
      USE ggm_stuff
      USE precise
      USE constants
      USE derivative_form
      IMPLICIT NONE
      INTEGER :: a,i,j,k,num,skips(7)
      INTEGER*8 :: isize(1,3),bound_elems,node_start,counter,nub,node_number
      INTEGER :: index_base,index_zone,index_coord,index_file,icelldim,&
        iphysdim,ier,index_flow,index_field,&
        index_section
      CHARACTER(len=32) :: filename,basename,zonename,&
        format_junk,char_save_number

        format_junk = '(I3.3)'

        index_base = 1
!        index_file = 1
        index_coord = 1
        index_zone = 1
        index_section = 1
        index_flow = 1
        index_field = 1
        node_start = 1

        bound_elems = 0

        IF (dimensions == 3) THEN
          icelldim = 3
          iphysdim = 3
        ELSE IF (dimensions == 2) THEN
          icelldim = 2
          iphysdim = 3
        END IF

!      WRITE(*,*) 'Preparing to allocate for output prep arrays.'
      ALLOCATE(box(num_x,num_y,num_z))
      ALLOCATE(done(num_x,num_y,num_z))
      ALLOCATE(nodes(8,link_number))
!      WRITE(*,*) 'Output prep arrays successfully allocated.'


      DO k = 1,num_z
      DO j = 1,num_y
      DO i = 1,num_x
        IF (state(dir+1,i,j,k) > 0) THEN
          box(i,j,k) = state(dir+1,i,j,k)
        ELSE
          box(i,j,k) = 0
        END IF
      END DO
      END DO
      END DO

      nodes = 0
      done = .FALSE.
      node_number = 1
      counter = 0
      DO k = 1,num_z
        DO j = 1,num_y
          inner_loop: DO i = 1,num_x
!            WRITE(*,*) 'Current i,j,k', i,j,k
            IF (state(dir+1,i,j,k) > 0) THEN
!            WRITE(*,*) 'Shooby doo!'
              IF (done(i,j,k)) THEN
                WRITE(*,*) 'Already done, skipping.'
                CYCLE inner_loop
!              WRITE(*,*) 'Already done, skipping.'
              ELSE
!
! Check for out of bounds points
! This will cause the next loop to skip out of bounds points
!
                skips = 0
!
                CALL card_dir_testing(i,j,k,skips,node_number)
!
              END IF
            END IF
!              CALL error_out
            IF (state(1,i,j,k) == 0 .AND. state(dir+1,i,j,k) >0) THEN
              node_loop: DO nub = 1,8
!
! Hit all the unassigned nodes
!                  WRITE(*,*) inside_counter, node_number
                IF (nodes(nub,node_number) == 0) THEN
!                    WRITE(*,*) 'Apply a new node number.'
                  counter = counter +1
                  nodes(nub,node_number) = counter
!                    counter = counter +1
!                    WRITE(11,*) 'Applied node ',counter, ' to element ',&
!                        node_number
                END IF
              END DO node_loop

              done(i,j,k) = .TRUE.
              node_number = node_number + 1
!                WRITE(*,*) inside_counter,node_number
!
            END IF
          END DO inner_loop

        END DO
      END DO
      WRITE(*,*) counter
      ALLOCATE(x_vertices(counter))
      ALLOCATE(y_vertices(counter))
      ALLOCATE(z_vertices(counter))
      x_vertices = -9999.9
      y_vertices = -9999.9
      z_vertices = -9999.9

      counter = 1
      node_number = 1
      DO k = 1, num_z
        DO j = 1, num_y
          DO i = 1, num_x

            IF (state(dir+1,i,j,k) > 0) THEN
              DO nub = 1,8
!              WRITE(*,*) nodes(nub,node_number),counter
                IF (nodes(nub,node_number)==counter) THEN

                  IF (x_vertices(counter) <= -9998.1) THEN
                    x_vertices(counter)=grid(1,i,j,k)+dir_node_x(nub)*dx
                  END IF

                  IF (y_vertices(counter) <= -9998.1) THEN
                    y_vertices(counter)=grid(2,i,j,k)+dir_node_y(nub)*dx
                  END IF

                  IF (z_vertices(counter) <= -9998.1) THEN
                    z_vertices(counter)=grid(3,i,j,k)+dir_node_z(nub)*dx
                  END IF

!                  WRITE(11,*) x_vertices(counter),y_vertices(counter),&
!                  z_vertices(counter),node_number,counter
!                  WRITE(11,*) x_vertices(counter),y_vertices(counter),&
!                  z_vertices(counter)

                counter = counter + 1
                END IF
!              counter = counter + 1
              END DO

            node_number = node_number + 1
            END IF
!            WRITE(*,*) x_vertices(counter),y_vertices(counter),&
!              z_vertices(counter),node_number,counter
          END DO
        END DO
      END DO
      counter = counter -1
      node_number = node_number - 1

      DO a = 1,num_ggm_vars

        WRITE(char_save_number,format_junk) a

!        WRITE(*,*) SIZE(x_vertices),SIZE(y_vertices),SIZE(z_vertices)

        filename = 'ggm_output_'//TRIM(char_save_number)//'.cgns'
        basename = 'outputbase_'//TRIM(char_save_number)
        zonename = 'outputzone_'//TRIM(char_save_number)

!        WRITE(*,*) filename,basename,zonename

        isize(1,1) = counter
        isize(1,2) = node_number
        isize(1,3) = 0
!        WRITE(*,*) 'Counter and node number ',counter,node_number
!        WRITE(*,*) isize,SIZE(x_vertices),SIZE(y_vertices),SIZE(z_vertices)
!        WRITE(*,*) bound_elems,solid_counter,&
!        node_start,counter,nub,node_number
!        WRITE(*,*) isize(1,1),isize(1,2),isize(1,3)
!
! Open a new file for all the data
!
        WRITE(11,*) ''
        WRITE(11,200) a
!        WRITE(11) ''
!        filename = 'grid.cgns'
        CALL cg_open_f(TRIM(filename),CG_MODE_WRITE,index_file,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'Open CGNS file done'
!
! Write the base with the dimensions
!
        WRITE(*,*) index_file,basename,icelldim,iphysdim,&
          index_base,ier
!        WRITE(*,*) ''
!        basename = 'basic'
        CALL cg_base_write_f(index_file,basename,icelldim,iphysdim,&
          index_base,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'Base write done'
!
! Write the zone with the size information
! It is unstructured because of the nature of LB grids. The cells are
! perfect cubes, but can be ignored because of geometry or other reasons
!
!        WRITE(*,*) index_file,index_base,zonename,isize,&
!         Unstructured,index_zone,ier
        CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
          Unstructured,index_zone,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'Zone write done'
!
! Write the x-coordinates
!
!        WRITE(*,*) index_file,index_base,index_zone,RealDouble,&
!          'CoordianteX',index_coord,ier
        CALL cg_coord_write_f(index_file,index_base,index_zone,RealDouble,&
          'CoordianteX',x_vertices,index_coord,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'X coord write done'
!
! Write the y-coordinates
!
!        WRITE(*,*) index_file,index_base,index_zone,RealDouble,&
!          'CoordianteY',index_coord,ier
        CALL cg_coord_write_f(index_file,index_base,index_zone,RealDouble,&
          'CoordianteY',y_vertices,index_coord,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'Y coord write done'
!
! Write the z-coordinates
!
!        WRITE(*,*) index_file,index_base,index_zone,RealDouble,&
!         'CoordianteZ',index_coord,ier
        CALL cg_coord_write_f(index_file,index_base,index_zone,RealDouble,&
          'CoordianteZ',z_vertices,index_coord,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'Z coord write done'
!
! Write the section with all information required to write the data
!
!       WRITE(*,*) index_file,index_base,index_zone,'Elem',&
!          node_start,node_number,bound_elems,&
 !          index_section,ier
        CALL cg_section_write_f(index_file,index_base,index_zone,'Elem',&
          HEXA_8,node_start,node_number,bound_elems,nodes,&
          index_section,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,*) 'Connectivity write done'

        CALL cg_close_f(index_file,ier)
        IF (ier /= CG_OK) CALL cg_error_exit_f
        WRITE(11,201) num

!        WRITE(*,*) SIZE(x_vertices),SIZE(y_vertices),SIZE(z_vertices),&
!          SIZE(nodes)

      END DO


      DEALLOCATE(box)
      DEALLOCATE(nodes)
      DEALLOCATE(done)
      DEALLOCATE(x_vertices)
      DEALLOCATE(y_vertices)
      DEALLOCATE(z_vertices)

      ALLOCATE(deriv_form(3,5,link_number))
      CALL deriv_form_checker

 200  FORMAT ('Starting write of CGNS file for save ',I8)
 201  FORMAT ('Close initial write of save ',I8,' done')

      END SUBROUTINE
