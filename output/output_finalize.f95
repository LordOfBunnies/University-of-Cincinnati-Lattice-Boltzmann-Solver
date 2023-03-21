SUBROUTINE output_finalize
!
! This finalizes the output files by putting in things like time
! dependecy information and
!
! Called by: uclbs_main
! Calls: error_out,external_code_finalizer,cpu_time
! External calls: cgp_open_f,cgp_error_exit_f,cg_goto_f,cg_biter_write_f,
!   cgp_close_f,cg_ziter_write_f,cgp_array_write_f,cgp_array_write_data_f
!
USE precise
USE grid_data
USE output_data
use zone_writing
USE linkwise
USE timing
USE constants
USE freestream_values
USE nml_output
USE cgns
use mpi
use amr_info_holder
use amr_processes, only: self,nprocs,commune
use amrex_base_module
use amrex_amr_module
IMPLICIT NONE
INTEGER :: ier,a,lvl,sol_max_lvl,fine,i,j,max_lvl,z,ts,l,num
integer :: array_div
!integer :: output_min_bound,output_max_bound,dummy_bound
INTEGER(cgsize_t) :: idata(2),min_arr(2),output_min_bound,output_max_bound,&
  dummy_bound,min_dummy(2),size_array(3),self_min_bounds(3),self_max_bounds(3)
INTEGER :: index_base,index_file,index_coord,index_zone,index_section,&
  index_flow,index_field,index_array,istat
!      REAL(KIND=dp) :: ts_start,ts_end,ts_save,
CHARACTER(len=32),ALLOCATABLE :: all_sol_names(:)
CHARACTER(len=32) :: filename,char_save_number,sol_num,zonename,&
  format_junk,format_junk_ts,filename_out,basename,format_nm_lvl_and_zone
logical :: marked

REAL(KIND=dp),ALLOCATABLE :: timing_data(:)
!      REAL,ALLOCATABLE :: timing_data(:)
REAL(KIND=dp) :: end_times,elapsed_cpu_time
!
!
write(*,*) 'Finalizing output files.'
!
!
!
filename_out = 'uclbs.out'
format_junk = '(I3.3)'
format_junk_ts = '(I9.9)'
format_nm_lvl_and_zone = '(I2.2,a1,I3.3)'
fine = amrex_get_finest_level()
min_arr = 1
dummy_bound = 1
min_dummy = 1
index_array = 1

OPEN(FILE=filename_out,UNIT = 11,STATUS='OLD',FORM = 'FORMATTED',&
  ACTION='WRITE',POSITION='APPEND',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Problems opening output file, closing.'
  WRITE(*,*) istat
  CALL error_out
END IF
!
! For all the output files using CGNS, they need to be finished up.
! This does that by writing the unsteady data into them and telling
! them what that data is.
!
DO num = 1,num_saves

  WRITE(11,201) num

!        filename = 'Output_file_'//CHAR(num)
  WRITE(char_save_number,format_junk) num

  filename = 'output_'//TRIM(char_save_number)//'.cgns'
  basename = 'outputbase_'//TRIM(char_save_number)
!        WRITE(*,*) filename
!        basename = 'output_'//TRIM(char_save_num)//'.cgns'
!        zonename = 'output_'//TRIM(char_save_num)//'.cgns'

  index_base = 1
  index_file = num+1
  index_coord = 1
  index_zone = 1
  index_section = 1
  index_flow = 1
  index_field = 1

  idata(1) = 32
  idata(2) = save_steps(num)%n_outsteps
!
  if (fine > max_lvl) then
    max_lvl = save_shapes(2,num)
  else
    max_lvl = fine
  end if
!
!  if (coarse_time < very_small_number .and. .not. amr) then !adapted mesh or not
!
!    if (allocated(zone_storage)) deallocate(zone_storage)
!    allocate(zone_storage(save_shapes(1,num):max_lvl))
!
!    do lvl = save_shapes(1,num),max_lvl
!    !  write(*,*) 'box level bounds',self,save_shapes(1,num),max_lvl
!      call zone_box_setup(lvl)
!    end do
!    call zone_output_num_setup(num,save_shapes(1,num),max_lvl)
!!
!!
!!
!    IF (total_time < save_times(1,num)) THEN !short time
!
!!      CALL cgp_open_f(filename,CG_MODE_MODIFY,index_file,ier)
!!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!      WRITE(11,*) 'Output file reopened for finalization.'
!!            WRITE(*,*) 'Output file reopened for finalization.'
!!
!!
!!
!!            WRITE(*,*) index_file,index_base,'TimeIterValues',&
!!              ts_save,ier
!      CALL cg_biter_write_f(index_file,index_base,'TimeIterValues',&
!        save_steps(num)%n_outsteps,ier)
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!      WRITE(11,*) 'Writing time iterations values to finalize file.'
!      WRITE(*,*) 'Writing time iterations values to finalize file.',idata
!!
!!
!!
!!            WRITE(*,*) index_file,index_base,TimeAccurate,&
!!              ier
!!      CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
!!        ier)
!!      IF (ier /= CG_OK) CALL cg_error_exit_f
!!      WRITE(11,*) 'Writing solution type to finalize file.'
!!      WRITE(*,*) 'Writing solution type to finalize file.'
!
!      CALL cgp_close_f(index_file,ier)
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!!
!    ELSE !normal time
!
!      !if (coarse_time > very_small_number .or. amr) then
!
!      if (fine >= save_shapes(2,num)) then
!        sol_max_lvl = save_shapes(2,num)
!      else
!        sol_max_lvl = fine
!      end if
!!
!      CALL cgp_open_f(filename,CG_MODE_MODIFY,index_file,ier)
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!      WRITE(11,*) 'Output file reopened for finalization.'
!!      WRITE(*,*) 'Output file reopened for finalization.',index_file
!!
!      CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',&
!        1,'end')
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!!      write(*,*) lbound(all_sol_names),ubound(all_sol_names)
!      index_array = 2
!!      write(*,*) 'step names',self,idata,save_steps(num)%sol_names
!!      write(*,*) save_shapes(1,a),sol_max_lvl,zone_storage(lvl)%output_number(1,a),&
!!          zone_storage(lvl)%output_number(zone_storage(lvl)%big_zones,a)
!      do lvl = save_shapes(1,num),sol_max_lvl
!        do i = zone_storage(lvl)%output_number(1,num),&
!          zone_storage(lvl)%output_number(zone_storage(lvl)%big_zones,num)
!
!!          write(*,*) 'cross check before output',index_file,index_base,i
!
!          CALL cg_ziter_write_f(index_file,index_base,i,&
!            'ZoneIterativeData',ier)
!          IF (ier /= CG_OK) CALL cgp_error_exit_f
!      !
!      ! Go to the zone data
!      !
!          CALL cg_goto_f(index_file,index_base,ier,'Zone_t',i,&
!            'ZoneIterativeData_t',1,'end')
!          IF (ier /= CG_OK) CALL cgp_error_exit_f
!      !
!      ! Tell it the solution names and their sizes
!      !
!!          write(*,*) 'solution names, making sure they match',self,all_sol_names
!!          CALL cgp_array_write_f('FlowSolutionPointers',Character,int(2,cgsize_t),idata,index_array,ier)
!          CALL cgp_array_write_f('FlowSolutionPointers',Character,2,idata,index_array,ier)
!
!          IF (ier /= CG_OK) CALL cgp_error_exit_f
!!          write(*,*) 'squishy face!!!',index_array,idata
!
!          if (save_steps(num)%n_outsteps >= nprocs) then
!
!            array_div = (save_steps(num)%n_outsteps-1)/nprocs
!
!            if (self == 0) then
!              output_min_bound = 1
!              output_max_bound = array_div
!            else if (self ==nprocs-1) then
!              output_min_bound = self*array_div+1
!              output_max_bound = save_steps(num)%n_outsteps
!            else
!              output_min_bound = self*array_div+1
!              output_max_bound = output_min_bound+array_div-1
!            end if
!
!              call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!                (/idata(1),output_max_bound/),save_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!              IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!          else
!            do j = 1,save_steps(num)%n_outsteps
!              if (self <= j) then
!                output_min_bound = self+1
!                output_max_bound = output_min_bound+1
!              end if
!            end do
!            if (self <= save_steps(num)%n_outsteps) then
!              call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!                (/idata(1),output_max_bound/), save_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!              IF (ier /= CG_OK) CALL cgp_error_exit_f
!            end if
!          end if
!
!!        call cgp_array_write_data_f(index_array,int((/1,1/),cgsize_t),idata,&
!!          all_sol_names,ier)
!!        call cg_error_print_f
!!        write(*,*) self,ier
!!        IF (ier /= CG_OK) CALL cg_error_exit_f
!        end do
!      end do
!
!      CALL cgp_close_f(index_file,ier)
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!      WRITE(11,200) num
!      WRITE(11,*) ''
!    END IF
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!  else ! Controls adapted mesh or not
!!    CALL cgp_open_f(filename,CG_MODE_MODIFY,index_file,ier)
!!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    WRITE(11,*) 'Output file reopened for finalization.'
!    WRITE(*,*) 'Output file reopened for finalization.'
!!
!!
!!
!    if (save_steps(num)%n_outsteps >= nprocs) then
!
!      array_div = (save_steps(num)%n_outsteps-1)/nprocs
!      if (self == 0) then
!        output_min_bound = 1
!        output_max_bound = array_div
!      else if (self ==nprocs-1) then
!        output_min_bound = self*array_div+1
!        output_max_bound = save_steps(num)%n_outsteps
!      else
!        output_min_bound = self*array_div+1
!        output_max_bound = output_min_bound+array_div-1
!      end if
!
!!      call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!!        (/idata(1),output_max_bound/),save_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!    else
!      do j = 1,save_steps(num)%n_outsteps
!        if (self <= j) then
!          output_min_bound = self+1
!          output_max_bound = output_min_bound+1
!        end if
!      end do
!!      if (self <= save_steps(num)%n_outsteps) then
!!!        call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!!!          (/idata(1),output_max_bound/), save_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!!        call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!!          (/idata(1),output_max_bound/), save_steps(num)%sol_names,ier)
!!
!!
!!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!!      end if
!    end if
!
!    size_array(1) = 3
!    size_array(2) = save_steps(num)%n_outsteps
!    size_array(3) = base_info(num)%max_zones
!!
!! zone_pointer_names contains zone names for the BaseIterativeData_t
!! amr_solution_names contains names for solutions
!!
!    allocate(base_info(num)%zone_pointer_names(base_info(num)%max_zones,save_steps(num)%n_outsteps))
!    allocate(base_info(num)%amr_solution_names(base_info(num)%max_zones,save_steps(num)%n_outsteps))
!
!    base_info(num)%zone_pointer_names = "Null"
!    base_info(num)%amr_solution_names = "Null"
!!
!! Get to the BaseIterativeData_t to write a bunch of stuff
!!
!! Write the number of zones at each timestep
!!
!    CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',1,'end')
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!    index_array = 2
!!    call cgp_array_write_f('NumberOfZones',INTEGER,dummy_bound,&
!!      int(save_steps(num)%n_outsteps,cgsize_t),index_array,ier)
!    call cgp_array_write_f('NumberOfZones',INTEGER,dummy_bound,&
!      (/save_steps(num)%n_outsteps/),index_array,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!    call cgp_array_write_data_f(index_array,INT(output_min_bound,cgsize_t),&
!      INT(output_max_bound,cgsize_t),base_info(num)%total_zones_at_ts(output_min_bound:&
!      output_max_bound),ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!
!      write(*,*) 'zones per ts',self,output_min_bound,output_max_bound,&
!      base_info(num)%total_zones_at_ts(output_min_bound:output_max_bound)
!
!!
!! Assign all the information to zone_pointer_names and amr_solution_names
!!
!!
!! For every timestep
!    do ts = 1,save_steps(num)%n_outsteps
!! For every zone
!      do z = 1,base_info(num)%max_zones
!
!        marked = .false.
!! Loop through all zones at a given timestep
!        do l = 1,base_info(num)%total_zones_at_ts(ts)
!! Loop through all the levels to satisfy grid_info
!          do lvl = save_shapes(1,num),save_shapes(2,num)
!
!          if (allocated(grid_info(lvl,num)%multigrid_save(ts)%zone_number)) then
!          !write(*,*) 'zones at ts',grid_info(lvl,num)%multigrid_save(ts)%zone_number
!! Check if any zones at the timestep match ones recorded during the run, if so
!!   put the zone name into zone_pointer_names.  Also put the solution name into
!!   amr_solution_names which goes into FlowSolutionPointers_t at the end
!            if (grid_info(lvl,num)%multigrid_save(ts)%zone_number(l) == z) then
!
!              base_info(num)%zone_pointer_names(z,ts) = save_steps(num)%sol_names(ts)
!
!              write(char_save_number,format_nm_lvl_and_zone) lvl,'_',z
!              zonename = 'zone_lvl_'//TRIM(char_save_number)
!
!              base_info(num)%amr_solution_names(z,ts) = zonename
!
!              marked = .true.
!            end if
!
!            if (marked) exit
!            end if
!          end do
!
!        end do
!
!      end do
!    end do
!!
!!   Write the zone with the size information
!!   It is unstructured because of the nature of LB grids. The cells are
!!   perfect cubes, but can be ignored because of geometry or other reasons
!!
!    CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',1,'end')
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!    index_array = 3
!    call cgp_array_write_f('ZonePointers',CHARACTER,int(3,cgsize_t),size_array,index_array,ier)
!
!    write(*,*) 'gigglestan',self,idata,index_array,size_array
!
!    self_min_bounds(1) = min_dummy(1)
!    self_min_bounds(2) = int(output_min_bound,cgsize_t)
!    self_min_bounds(3) = min_dummy(1)
!    self_max_bounds(1) = idata(1)
!    self_max_bounds(2) = int(output_max_bound,cgsize_t)
!    self_max_bounds(3) = int(save_steps(num)%n_outsteps,cgsize_t)
!
!    call cgp_array_write_data_f(index_array,self_min_bounds,self_max_bounds,&
!      base_info(num)%amr_solution_names(output_min_bound:output_max_bound,1:&
!      save_steps(num)%n_outsteps),ier)
!!
!! Divide the arrays up roughly so the zone_names and solution_names
!!   can be written to the CGNS file
!!
!
!    do z = 1,base_info(num)%max_zones
!      index_array = 4
!      CALL cg_ziter_write_f(index_file,index_base,z,'ZoneIterativeData' ,ier)
!      IF (ier /= CG_OK) CALL cgp_error_exit_f()
!      !
!      ! Go to the zone data
!      !
!      CALL cg_goto_f(index_file,index_base,ier,'Zone_t',z,'ZoneIterativeData_t',1,'end')
!      IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!      call cgp_array_write_f('FlowSolutionPointers',Character,int(2,cgsize_t),idata,index_array,ier)
!
!      write(*,*) 'vert de flerf',self,z,index_array,base_info(num)%zone_pointer_names(z,output_min_bound:output_max_bound)
!
!      call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/), &
!          (/idata(1),output_max_bound/), &
!          base_info(num)%zone_pointer_names(z,output_min_bound:output_max_bound),ier)
!    end do
!
!
!  end if

  !DEALLOCATE(timing_data)
  !DEALLOCATE(all_sol_names)

call cgp_close_f(index_file,ier)

END DO



!      CALL ggm_output_finalize !XXXX

CALL CPU_TIME(cpu_stop_time)
elapsed_cpu_time = cpu_stop_time - cpu_start_time
WRITE(11,501) elapsed_cpu_time
write(*,501) elapsed_cpu_time


WRITE(11,*) 'End of UCLBS. Program completed successfully.'
WRITE(11,*) 'Thank you for using UCLBS.'

CLOSE(11)


WRITE(*,*) 'End of UCLBS. Program completed successfully.'
WRITE(*,*) 'Thank you for using UCLBS.'

 200  FORMAT ('Output of save ',I4,' done.')
 201  FORMAT ('Beginning finalization of output file for save ',I4)
 501  FORMAT ('Total elapsed time running UCLBS = ',F14.6)



END SUBROUTINE


      !
! Find the first and last timestep recorded
!
!  ts_start = CEILING(save_times(1,a)/dt)
!  winkin = MODULO(ts_start,save_frequency(num))
!  ts_first_save = ts_start + (save_frequency(num) - winkin)
!
!  end_times = MIN(total_time,save_times(2,a))
!
!  ts_end = FLOOR(end_times/dt)
!  blinkin = MODULO(ts_end,save_frequency(num))
!  ts_last_save = ts_end - blinkin
!
!  max_time_end = FLOOR(total_time/dt)
!  noggin = MODULO(max_time_end,save_frequency(num))
!  time_max_save = max_time_end - noggin
!
!  enders_timer = MIN(time_max_save,ts_last_save)
!
! Find the number of timesteps saved over the course of the program
!
!        ts_save = (enders_timer - ts_first_save)/save_frequency(num) + 1
!        WRITE(*,*) ts_start,ts_end,winkin,blinkin,end_times,enders_timer
!        WRITE(*,*) ts_save,ts_last_save,ts_first_save,save_frequency(num)
!        WRITE(*,*) total_time,ts_first_save*dt

!        idata(1) = 32
!        idata(2) = ts_save
!
! Allocate timing data and solution names
!
!  ALLOCATE(all_sol_names(save_steps(num)%n_outsteps))
!        nod = 1
!  DO i = 1,ts_save
!
!    timing_data(i) = DBLE(i)*dt*DBLE(save_frequency(num))
!!          timing_data(i) = REAL(i)*dt*REAL(save_frequency(num))
!  WRITE(solut_save_number,format_junk_ts) i*save_frequency(num)
!  solut_name(i) = 'Timestep_'//TRIM(solut_save_number)
!!          WRITE(*,*) timing_data(i),solut_name(i)
!!          solut_name(nod) = 'Timestep_'//CHAR(i)
!!          nod = nod+1
!  END DO
!  WRITE(*,*) ts_save,idata



!          WRITE(*,*) filename
!          WRITE(*,*) ''
!          WRITE(*,*) ts_start,ts_end,ts_first_save,&
!        ts_last_save,winkin,blinkin,enders_timer,time_max_save,&
!        max_time_end,noggin,i
!          WRITE(*,*) ''
!          WRITE(*,*) ier,a,ts_save,istat
!          WRITE(*,*) ''
!          WRITE(*,*) idata(2),funk,end_times,timing_data
!          WRITE(*,*) ''
!          WRITE(*,*) index_base,index_file,index_coord,index_zone,index_section,&
!        index_flow,index_field
!          WRITE(*,*) ''
!          WRITE(*,*) filename,char_save_number,solut_save_number,&
!        format_junk,format_junk_ts,filename_out,basename
!          WRITE(*,*) ''
!          WRITE(*,*) solut_name(:)
!          WRITE(*,*) ''
!          WRITE(*,*) x_start,x_end,y_start,y_end,z_start,z_end
!          WRITE(*,*) num_nodes_in_save(:)



!          WRITE(*,*) 'Writing zone header to finalize file.'
!
! Go to the zone data
!
!          CALL cg_goto_f(index_file,index_base,ier,'Zone_t',index_zone,&
!            'ZoneIterativeData_t',1,'end')
!          IF (ier /= CG_OK) CALL cg_error_exit_f
!!
!! Tell it the solution names and their sizes
!!
!          CALL cg_array_write_f('FlowSolutionPoints',Character,INT(2,cgsize_t),idata,&
!            all_sol_names,ier)
!          IF (ier /= CG_OK) CALL cg_error_exit_f
!          WRITE(11,*) 'Writing solution names to finalize file.'
!          WRITE(*,*) 'Writing solution names to finalize file.'
!
! Tell it the solution is time accurate (no other way to do LBM)
!
!          CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
!            ier)
!          IF (ier /= CG_OK) CALL cgp_error_exit_f
!          WRITE(11,*) 'Writing solution type to finalize file.'
!          WRITE(*,*) 'Writing solution type to finalize file.'
!
! Close out the file
!
!          if (any(grid_info(save_shapes(1,num),save_shapes(2,num),num)%&
!              multigrid_save(ts)%zone_number) == z) then
!
! Open the file to be appended to
!
!         WRITE(*,*) ''
!          WRITE(*,*) 'Begin finalization of output save ', a
!  WRITE(11,*) 'Begin finalization of output save ', a
!  CALL cg_open_f(filename,CG_MODE_MODIFY,index_file,ier)
!  IF (ier /= CG_OK) CALL cg_error_exit_f
!  WRITE(11,*) 'Output file reopened for finalization.'
!          WRITE(*,*) 'Output file reopened for finalization.'
!
! Call the base iterative write to tell it how many timesteps there are
!
!          CALL cg_goto_f(index_file,index_base,ier,basename,&
!            1,'end')
!          IF (ier /= CG_OK) CALL cg_error_exit_f
!          WRITE(*,*) 'Going to base iterative data.'

!          WRITE(*,*) index_file,index_base,'TimeIterValues',&
!            ts_save,ier
!    do i = 1,save_steps(num)%n_outsteps
!      write(sol_num,format_junk) save_steps(num)%saved_timestep(i)
!      all_sol_names(i) = 'Timestep_'//trim(sol_num)
!      write(*,*) self,all_sol_names(i)
!    end do
!    CALL cg_biter_write_f(index_file,index_base,'NumberOfZones',&
!      base_info(num)%total_zones_at_ts(1:save_steps(num)%n_outsteps),ier)
!    CALL cg_biter_write_f(index_file,index_base,'NumberOfSteps',&
!      save_steps(num)%n_outsteps,ier)

!    call cgp_array_write_f('NumberOfZones',INTEGER,(/dummy_bound,output_min_bound/),&
!      (/idata(1),output_max_bound/),index_array,ier)
!
!    if (save_steps(num)%n_outsteps >= nprocs) then
!
!      array_div = (save_steps(num)%n_outsteps-1)/nprocs
!
!      if (self == 0) then
!        output_min_bound = 1
!        output_max_bound = array_div
!      else if (self ==nprocs-1) then
!        output_min_bound = self*array_div+1
!        output_max_bound = save_steps(num)%n_outsteps
!      else
!        output_min_bound = self*array_div+1
!        output_max_bound = output_min_bound+array_div-1
!      end if
!
!      call cgp_array_write_data_f(index_array,(/dummy_bound,dummy_bound,output_min_bound/), &
!          (/idata(1),int(base_info(num)% max_zones,kind=cgsize_t),output_max_bound/), &
!          save_steps(num)% sol_names(output_min_bound:output_max_bound),ier)
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!    else
!      do j = 1,save_steps(num)%n_outsteps
!        if (self <= j) then
!          output_min_bound = self+1
!          output_max_bound = output_min_bound+1
!        end if
!      end do
!      if (self <= save_steps(num)%n_outsteps) then
!        call cgp_array_write_data_f(index_array,(/dummy_bound,dummy_bound,&
!          output_min_bound/), (/idata(1),int(base_info(num)% max_zones,kind=cgsize_t),output_max_bound/), &
!          save_steps(num)% sol_names(output_min_bound:output_max_bound),ier)
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!      end if
!    end if
!    call cg_
