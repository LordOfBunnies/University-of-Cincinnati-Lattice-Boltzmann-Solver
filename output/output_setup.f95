SUBROUTINE output_setup
!
! This is to create the link this of all the output save regions
!
! Called by: links
! Calls: zone_box_setup
!
! Outside Calls: cg_open_f,cg_base_write_f,cg_zone_write_f,cg_coord_write_f,
!   cg_section_write_f,cg_close_f
!
USE output_data, only: save_steps,base_info
use constants
use zone_writing
!USE grid_data
!USE linkwise
USE nml_output
USE precise
USE cgns
USE amr_info_holder
use amr_processes
use amrex_amr_module
use amrex_base_module
use timing
use mpi
use ggm_stuff, only: use_ggm,ggm_output
IMPLICIT NONE
INTEGER :: i,j,num,fine,out_max,out_min,lvl,max_lvl,dummy
!INTEGER*8 :: bound_elems,solid_counter,node_start
!real(kind=dp),allocatable :: time_vals(:)
real(kind=dp) :: t_diff,time_start,time_stop,ts
!      REAL(KIND = dp) :: radius,x_dist,y_dist,z_dist
!      REAL(KIND = dp) :: x_center,y_center,z_center
!      REAL(KIND = dp) :: x_max,y_max,z_max,x_min,y_min,z_min

INTEGER :: index_base,index_zone,index_coord,index_file,icelldim,&
  iphysdim,ier,index_flow,index_field,index_array,&
  index_section,output_min_bound,output_max_bound,array_div
CHARACTER(len=32) :: filename,basename,sol_num,&
  format_junk,char_save_number,format_junk_ts,format_lvl_and_timestep


!write(*,*) 'output setup',self
fine = amrex_get_finest_level()
dummy = 1
!write(*,*) 'meh?',self,fine

WRITE(11,*) ''
WRITE(11,*) 'Begin setup of output files'
!WRITE(*,*) 'Begin setup of output files',self
format_junk = '(I3.3)'
format_junk_ts = '(I9.9)'
format_lvl_and_timestep = '(I2.2,a1,I9.9)'

if (allocated(save_steps)) deallocate(save_steps)

allocate(save_steps(num_saves))
allocate(base_info(num_saves))

!if (allocated(multigrid_save)) deallocate(multigrid_save)
!allocate(multigrid_save(num)
DO num = 1,num_saves
  index_file = num+1
  index_base = 1
  index_coord = 1
  index_zone = 1
  index_section = 1
  index_flow = 1
  index_field = 1
  index_array = 1
  ier =0

  IF (dimensions ==3) THEN
    icelldim = 3
    iphysdim = 3
  ELSE
    icelldim = 2
    iphysdim = 3
  END IF
!
!  Check output information
!
!  write(*,*) 'output info, start time ',save_times(1,num),' total_time ',total_time,&
!    'stop time',save_times(2,num),'local_timestep',dt(save_shapes(1,num)),'save frequency',&
!    save_frequency(num)
!
! Make sure the start time is positive
!
  if (save_times(1,num) < very_small_number) then
    time_start = 0.0D0
  else
    time_start = save_times(1,num)
  end if
!
! Catchall for late starts
!
  if (save_times(1,num) > total_time) then
    write(*,*) 'writing no data for save ',num
    write(11,*) 'writing no data for save',num
  end if
!
! Making sure of the total recording of the time
!
  if (save_times(2,num) > total_time) then
    time_stop = total_time
  else
    time_stop = save_times(2,num)
  end if
  t_diff = time_stop-time_start
!
! Number of saves
!
!  write(*,*) 'boop de snoot',self

  save_steps(num)%n_outsteps = ceiling(t_diff/save_frequency(num))
  save_steps(num)%n_outsteps = save_steps(num)%n_outsteps+1
! write(*,*) 'number of output steps ',save_steps(num)%n_outsteps

  if (allocated(save_steps(num)%saved_times)) deallocate(save_steps(num)%saved_times)
  allocate(save_steps(num)%saved_times(save_steps(num)%n_outsteps))

  if (allocated(save_steps(num)%saved_timestep)) deallocate(save_steps(num)%saved_timestep)
  allocate(save_steps(num)%saved_timestep(save_steps(num)%n_outsteps))

  if (allocated(save_steps(num)%sol_names)) deallocate(save_steps(num)%sol_names)
  allocate(save_steps(num)%sol_names(save_steps(num)%n_outsteps))

  if (allocated(save_steps(num)%write_done)) deallocate(save_steps(num)%write_done)
  allocate(save_steps(num)%write_done(save_steps(num)%n_outsteps))

  save_steps(num)%write_done = .false.
!
!  write(*,*) 'bounce'
  ts = 0.0D0
  i = 1
!  do while (ts < save_times(2,num))
  do while (ts < total_time)
    ts = time_start+(i-1)*save_frequency(num)
!      write(*,*) 'ts ', ts
    save_steps(num)%saved_times(i) = ts
    save_steps(num)%saved_timestep(i) = FLOOR(ts/dt(save_shapes(2,num)))

    write(sol_num,format_junk_ts) save_steps(num)%saved_timestep(i)
    save_steps(num)%sol_names(i) = 'Timestep_'//TRIM(sol_num)
!  write(*,*) 'jounce'
!      write(*,*) 'save_times, save_timesteps, solution names',&
!        save_steps(num)%saved_times(i),save_steps(num)%saved_timestep(i),&
!        save_steps(num)%sol_names(i),ts
    i = i+1
  end do

  WRITE(char_save_number,format_junk) num
!    write(*,*) char_save_number
!          WRITE(*,*) SIZE(x_vertices),SIZE(y_vertices),SIZE(z_vertices)

  filename = 'output_'//TRIM(char_save_number)//'.cgns'
  basename = 'outputbase_'//TRIM(char_save_number)
!    zonename = 'outputzone_'//TRIM(char_save_number)

!
!   Open a new file for all the data
!
!    write(*,*) 'boogey'
  WRITE(11,*) ''
  WRITE(11,200) num

  CALL cgp_open_f(trim(filename),CG_MODE_WRITE,index_file,ier)
    !write(*,*) index_file,ier
  IF (ier /= CG_OK) CALL cgp_error_exit_f
  WRITE(11,*) 'Open CGNS file done'
!    write(*,*) 'Open CGNS file done',self,index_file
!
!   Write the base with the dimensions
!
  CALL cg_base_write_f(index_file,basename,icelldim,iphysdim,index_base,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f

  CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!  write(*,*) 'CGNS simulation type written'
! XXXX big region to be changed for future!!!!!

  CALL cg_biter_write_f(index_file,index_base,'NumberOfSteps',&
    save_steps(num)%n_outsteps,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!!     write(*,*) 'Writing time values at start of file.',self,save_steps(num)%n_outsteps
!!
!!   Write the zone with the size information
!!   It is unstructured because of the nature of LB grids. The cells are
!!   perfect cubes, but can be ignored because of geometry or other reasons
!!
    CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',&
      1,'end')
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!   Both parts of the array write are needed and the cgsize_t argument is important as CGNS freaks
!   out if the wrong size is passed, but it doesn't always call it out.
!
!
!   This writes the number of timesteps to the CGNS file.  The actual writing of the timestep numbers
!   is done in the cgp_array_data_write_f below.
!
    CALL cgp_array_write_f('TimeValues',RealDouble,INT(1,cgsize_t),&
      INT(save_steps(num)%n_outsteps,cgsize_t),index_array,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f

!    write(*,*) 'squish'

    if (i-1 >= nprocs) then

      array_div = (i-1)/nprocs

      if (self == 0) then
        output_min_bound = 1
        output_max_bound = array_div
      else if (self ==nprocs-1) then
        output_min_bound = self*array_div+1
        output_max_bound = i-1
      else
        output_min_bound = self*array_div+1
        output_max_bound = output_min_bound+array_div-1
      end if

!      write(*,*) 'bonk',self,output_min_bound,output_max_bound,array_div
!      write(*,*) 'useful data',self,save_steps(num)%saved_times(output_min_bound:output_max_bound)
!      call cgp_array_write_data_f(index_array,INT(output_min_bound,cgsize_t),&
!      INT(output_max_bound,cgsize_t),save_steps(num)%saved_times(output_min_bound:output_max_bound),ier)
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!      write(*,*) 'boof',self
    else
!      write(*,*) 'doing the other thing',self
      do j = 0,i-1
        if (self == j) then
          output_min_bound = self+1
          output_max_bound = output_min_bound
!
!        else
!          output_min_bound = 0
!          output_max_bound = 0
        end if
      end do

!      if (self < j-1) then
!!        write(*,*) 'output junk',output_min_bound,output_max_bound,self,j,&
!!          save_steps(num)%saved_times
!        call cgp_array_write_data_f(index_array,INT(output_min_bound,cgsize_t),&
!        INT(output_max_bound,cgsize_t),save_steps(num)%saved_times(output_min_bound:output_max_bound),ier)
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!!        write(*,*) 'doing the after thing',self
!      end if
    end if

!    write(*,*) 'wiggles!',self
!
!
! If this is just a static multi-level run, we can output all the grid and zone information now.
!
!
  if (coarse_time < very_small_number .and. .not. amr) then
    if (.not. amr) then
      if (save_shapes(2,num) > fine) then
        out_max = fine
      else
        out_max = save_shapes(2,num)
      end if

      out_min = save_shapes(1,num)
!    write(*,*) 'zone storage bounds',self,out_min,out_max,index_file,num
      if (allocated(zone_storage)) deallocate(zone_storage)
      allocate(zone_storage(out_min:out_max))
!
!    do lvl = out_min,out_max
!      call zone_box_setup(lvl)
!    end do

      call zone_writer(out_min,out_max,num)
      call zone_coord_data_write(index_file,index_base,out_min,out_max,num)
    else
      call zone_dealloc_safety(fine,dummy,save_shapes(1,num),save_shapes(2,num))
    end if

    save_steps(num)%output_counter = 1
!
!
! For anything with either valid coarse time or adaptive mesh, it gets more complicated.
!
!
  else
 !   write(*,*) 'gibber',self
!    if(allocated(multigrid_save)) deallocate(multigrid_save)
!    if (save_shapes(1,num) < fine) return

    do lvl = save_shapes(1,num),save_shapes(2,num)
      allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%n_outsteps))
    end do
    allocate(base_info(num)%total_zones_at_ts(save_steps(num)%n_outsteps))
    allocate(base_info(num)%zones_per_lvl(save_shapes(1,num):save_shapes(2,num),&
      save_steps(num)%n_outsteps))

    do lvl = save_shapes(1,num),save_shapes(2,num)
      grid_info(lvl,num)%multigrid_save(1:save_steps(num)%n_outsteps)%regrid_between_outputs = .false.
    end do
!
! Regrid between outputs math
!
    do i = 1,save_steps(num)%n_outsteps-1
      do j = 1,size(regrid_times)
!        write(*,*) 'doof',size(regrid_times),save_steps(num)%saved_times(i),regrid_times(j)

!          write(*,*) 'woof',coarse_time,save_steps(num)%saved_times(i),&
!            save_steps(num)%saved_times(i+1),regrid_times(j)

        if (i == 1) then
!
          do lvl = save_shapes(1,num),save_shapes(2,num)
!            if (lvl > fixed_lvl) then
              grid_info(lvl,num)%multigrid_save(i)%regrid_between_outputs = .true.
!            end if
          end do

          exit
        else if (regrid_times(j) >= save_steps(num)%saved_times(i) .and. &
          regrid_times(j) <= save_steps(num)%saved_times(i+1)) then
!
          do lvl = save_shapes(1,num),save_shapes(2,num)
            if (lvl > fixed_lvl) then
              grid_info(lvl,num)%multigrid_save(i+1)%regrid_between_outputs = .true.
            end if
          end do
!          write(*,*) 'gong!',regrid_times(j),save_steps(num)%saved_times(i),&
!            save_steps(num)%saved_times(i+1),lvl,fixed_lvl

          exit
        else if (coarse_time >= save_steps(num)%saved_times(i) .and. &
          coarse_time <= save_steps(num)%saved_times(i+1)) then
!          write(*,*) 'woof',coarse_time,save_steps(num)%saved_times(i),&
!            save_steps(num)%saved_times(i+1)
!
          do lvl = save_shapes(1,num),save_shapes(2,num)
            if (lvl > fixed_lvl) then
             grid_info(lvl,num)%multigrid_save(i+1)%regrid_between_outputs = .true.
            end if
          end do
          exit
        else
          !write(*,*) 'no go',
          cycle
        end if
      end do
    end do

    base_info(num)%adapted_zone_start_number = 1
!
!
!
!    write(*,*) 'normal',self,grid_info(save_shapes(2,num),num)%&
!      multigrid_save(1:save_steps(num)%n_outsteps)%regrid_between_outputs
    save_steps(num)%output_counter = 1
  end if

!  write(*,*) 'after action report',self
!  CALL cgp_close_f(index_file,ier)
!  IF (ier /= CG_OK) CALL cgp_error_exit_f
  WRITE(11,201) num
!  call zone_regrid_checks

END DO

!if (amr) call zone_dealloc_safety(fine,save_shapes(1,num),save_shapes(2,num))

!if (ggm_output) then
!  call ggm_output_setup(fine)
!end if
WRITE(11,*) ''
!
! This section deals with getting the restart saves setup
!

 200  FORMAT ('Starting write of CGNS file for save ',I8)
 201  FORMAT ('Close initial write of save ',I8,' done')



END SUBROUTINE

!
! 1 = -dx,-dy,-dz
! 2 = -dx,+dy,-dx
! 3 = +dx,+dy,-dz
! 4 = +dx,-dy,-dz
! 5 = -dx,-dy,+dz
! 6 = -dx,+dy,+dz
! 7 = +dx,+dy,+dz
! 8 = +dx,-dy,+dz
!
!        WRITE(11,*) inside_counter
!        WRITE(11,*) x_start,x_end,y_start,y_end,z_start,z_end
!        WRITE(*,*) x_start,x_end,y_start,y_end,z_start,z_end
!        WRITE(*,*) 'Inside counter ',inside_counter
!        WRITE(11,*) ''
!        ALLOCATE(done(x_start:x_end,y_start:y_end,z_start:z_end))
!        inside_counter = inside_counter*2
!        ALLOCATE(nodes(8,inside_counter))
!        inside_counter = inside_counter/2
!        nodes = 0
!        ALLOCATE(box(x_start:x_end,y_start:y_end,z_start:z_end))
!        done = .FALSE.
!        counter = 0
!        node_number = 1
!  DO k = z_start,z_end
!    DO j = y_start,y_end
!      inner_loop: DO i = x_start,x_end
!
!!            WRITE(*,*) 'Current i,j,k', i,j,k
!
!        IF (box(i,j,k) /= 0) THEN
!!            WRITE(*,*) 'Shooby doo!'
!          IF (done(i,j,k)) THEN
!            WRITE(*,*) 'Already done, skipping.'
!            CYCLE inner_loop
!!              WRITE(*,*) 'Already done, skipping.'
!          ELSE
!!
!! Check for out of bounds points
!! This will cause the next loop to skip out of bounds points
!!
!                  skips = 0
!!
!!                  WRITE(*,*) x_start,x_end,y_start,y_end,z_start,z_end
!!                  WRITE(*,*) i,j,k,node_number
!            CALL card_dir_testing(i,j,k,skips,node_number)
!!
!          END IF
!        END IF
!!              CALL error_out
!        IF (state(1,i,j,k) == 0 .AND. box(i,j,k) >0) THEN
!          node_loop: DO nub = 1,8
!!
!! Hit all the unassigned nodes
!!                  WRITE(*,*) inside_counter, node_number
!            IF (nodes(nub,node_number) == 0) THEN
!!                    WRITE(*,*) 'Apply a new node number.'
!              counter = counter +1
!              nodes(nub,node_number) = counter
!!                    counter = counter +1
!!                    WRITE(11,*) 'Applied node ',counter, ' to element ',&
!!                        node_number
!            END IF
!          END DO node_loop
!!
!!                WRITE(11,*) nodes(1,node_number),nodes(2,node_number), &
!!           nodes(3,node_number),nodes(4,node_number),nodes(5,node_number),&
!!           nodes(6,node_number),nodes(7,node_number),nodes(8,node_number)
!!                done(i,j,k) = .TRUE.
!!                node_number = node_number + 1
!!                WRITE(*,*) inside_counter,node_number
!!
!        END IF
!      END DO inner_loop
!
!    END DO
!  END DO
!
!
!  counter = 1
!  node_number = 1
!  DO k = z_start, z_end
!    DO j = y_start, y_end
!      DO i = x_start, x_end

!        IF (box(i,j,k) > 0) THEN
!          DO nub = 1,8
!              WRITE(*,*) nodes(nub,node_number),counter
!            IF (nodes(nub,node_number)==counter) THEN

!              IF (x_vertices(counter) <= -9998.1) THEN
!          x_verts(counter)=grid(1,i,j,k)+dir_node_x(nub)*dx
!              END IF

!              IF (y_vertices(counter) <= -9998.1) THEN
!          y_verts(counter)=grid(2,i,j,k)+dir_node_y(nub)*dx
!              END IF

!              IF (z_vertices(counter) <= -9998.1) THEN
!          z_verts(counter)=grid(3,i,j,k)+dir_node_z(nub)*dx
!              END IF
!                  counter = counter + 1
!            END IF
!              counter = counter + 1
!          END DO

!          CALL output_link_num_save(num,node_number,i,j,k)

!          node_number = node_number + 1
!        END IF
!            WRITE(*,*) x_vertices(counter),y_vertices(counter),&
!              z_vertices(counter),node_number,counter
!      END DO
!    END DO
!  END DO
!  counter = counter -1
!  node_number = node_number - 1

!        DO i = 1,node_number
!          WRITE(*,*) vol_info_1(i)
!        END DO

!        filename = 'Output_file_'//CHAR(num)
!        basename = 'Base_'//CHAR(num)
!        zonename = 'Zone_'//CHAR(num)


!    allocate(save_steps(num)%multigrid_save(save_steps(num)%n_outsteps))
!! Allocate, for every possible output level and output timestep, whether things on that level
!!   can be regridded
!    do n = 1,save_steps(num)%n_outsteps
!      allocate(save_steps(num)%multigrid_save(n)%regrid_between_outputs(saves_shapes(1,num):&
!        save_shapes(2,num)))
!      allocate(save_steps(num)%multigrid_save(n)%zones_per_level(saves_shapes(1,num):&
!        save_shapes(2,num)))
!    end do
!
!    if (n == 1) then
!      save_steps(num)%multigrid_save(n)%regrid_between_outputs = .FALSE.
!    end if
!



!  write(*,*) 'number of outputs',save_steps(num)%n_outsteps
!  write(*,*) 'Saved times',i-1,save_steps(num)%saved_times
!  write(*,*) 'Saved timesteps',i-1,save_steps(num)%saved_timestep


!    WRITE(11,*) 'Base write done'
!    write(*,*) 'woogey'
!    WRITE(*,200) num
    !filename = 'output.cgns'
!    write(*,*) filename,index_file,ier

!   call mpi_barrier(commune,ier)
!          WRITE(11) ''
!          filename = 'grid.cgns'
