subroutine ggm_output_setup(fine)
!
! Set up the output files for studying GGM and PGS values
!
!
!
! Called by:
! Calls:
! External calls:
!
USE output_data, only: ggm_steps,ggm_base_info
use precise
use constants
use nml_output, only: num_saves
use mpi
use timing
use cgns
use amr_info_holder
use amr_processes
use zone_writing
use ggm_stuff
use amrex_amr_module
use amrex_base_module
use ggm_amr_writing
implicit none
INTEGER :: i,j,num,fine,out_max,out_min,lvl,max_lvl,pass_num,valid_ggm_outs,dummy
!integer,allocatable :: temp_output_timesteps(:,:)
real(kind=dp) :: t_diff,ts
!real(kind=dp),allocatable :: temp_output_times(:)
!
INTEGER :: index_base,index_zone,index_coord,index_file,icelldim,&
  iphysdim,ier,index_flow,index_field,index_array,&
  index_section,output_min_bound,output_max_bound,array_div,num_ggm_outputs
CHARACTER(len=32) :: filename,basename,sol_num,&
  format_junk,char_save_number,format_junk_ts
!
!
!
!
!
format_junk = '(I3.3)'
format_junk_ts = '(I9.9)'
dummy = 1
!

allocate(ggm_steps(num_ggm_vars))
allocate(ggm_base_info(num_ggm_vars))
!temp_output_times = -1.0D0
!temp_output_timesteps = -1
t_diff = total_time

if (allocated(ggm_steps)) deallocate(ggm_steps)
allocate(ggm_steps(num_ggm_vars))

do num = 1,num_ggm_vars
!  index_file = num +1
  index_file = num + num_saves + 1
  index_base = 1
  index_coord = 1
  index_zone = 1
  index_section = 1
  index_flow = 1
  index_field = 1
  ier =0

  IF (dimensions == 3) THEN
    icelldim = 3
    iphysdim = 3
  ELSE
    icelldim = 2
    iphysdim = 3
  END IF

  num_ggm_outputs = floor(total_time/ggm_output_freq(num)) + 1
!  write(*,*) 'ggm outs',num_ggm_outputs
!  allocate(temp_output_times(num_ggm_outputs))
!  allocate(temp_output_timesteps(num_ggm_outputs,0:amrex_max_level))
  if (allocated(ggm_steps(num)%saved_times)) deallocate(ggm_steps(num)%saved_times)
  allocate(ggm_steps(num)%saved_times(num_ggm_outputs))

  if (allocated(ggm_steps(num)%saved_timestep)) deallocate(ggm_steps(num)%saved_timestep)
  allocate(ggm_steps(num)%saved_timestep(num_ggm_outputs))

  if (allocated(ggm_steps(num)%sol_names)) deallocate(ggm_steps(num)%sol_names)
  allocate(ggm_steps(num)%sol_names(num_ggm_outputs))

  if (allocated(ggm_steps(num)%write_done)) deallocate(ggm_steps(num)%write_done)
  allocate(ggm_steps(num)%write_done(num_ggm_outputs))

    ggm_steps(num)%write_done = .false.
    ggm_steps(num)%n_outsteps = num_ggm_outputs

  WRITE(char_save_number,format_junk) num
  filename = 'ggm_out_'//TRIM(char_save_number)//'.cgns'
  basename = 'ggm_base_'//TRIM(char_save_number)

  ggm_steps(num)%output_counter = 1
  valid_ggm_outs = 0
!
! If the output frequency is greater than the regrid frequency, bind the output freq to the
!   regrid freq as you don't want to calculate extra ggms.
!

!
! If the ggm output frequency is less than the regrid frequency (regrid more often than output), then
!   find which regrids get output attached.
!
! As the regrid timesteps are known, it's possible to loop through them and find if an output would occur
!   before the next regrid, if so, the output is bound to that regrid.
!
    ts = ggm_output_freq(num)
    i = 1
    do while (ts < total_time + ggm_output_freq(num)/2.0D0)
!      do j = i,num_regrids
      ts = ggm_output_freq(num)*i
      ggm_steps(num)%saved_times(i) = ts
      ggm_steps(num)%saved_timestep(i) = floor(ts/dt(ggm_max_lvl(num)))
!      write(*,*) 'saved_times',ggm_steps(num)%saved_times(i),ggm_steps(num)%saved_timestep(i)
      write(sol_num,format_junk_ts) ggm_steps(num)%saved_timestep(i)
      ggm_steps(num)%sol_names(i) = 'Timestep_'//TRIM(sol_num)
      i = i+1
    end do

!  write(*,*) 'spinkle',filename,index_file,ier,ggm_steps(num)%n_outsteps

  CALL cgp_open_f(trim(filename),CG_MODE_WRITE,index_file,ier)
!  call cg_error_print_f
!  write(*,*) index_file,ier
  IF (ier /= CG_OK) CALL cgp_error_exit_f
  WRITE(11,*) 'Open CGNS file done'
!  write(*,*) 'Open CGNS file done',self,index_file
!
! Write the base with the dimensions
!
  CALL cg_base_write_f(index_file,basename,icelldim,iphysdim,&
    index_base,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!  WRITE(*,*) 'Base write done'
!
! Write the zone with the size information
! It is unstructured because of the nature of LB grids. The cells are
! perfect cubes, but can be ignored because of geometry or other reasons
!
  CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
    ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f

  CALL cg_biter_write_f(index_file,index_base,'NumberOfSteps',&
    ggm_steps(num)%n_outsteps,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f


  CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',&
    1,'end')
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!!  write(*,*) 'rawr',self
!!
!! Both parts of the array write are needed and the cgsize_t argument is important as CGNS freaks
!! out if the wrong size is passed, but it doesn't always call it out.
!!
!!
!! This writes the number of timesteps to the CGNS file.  The actual writing of the timestep numbers
!! is done in the cgp_array_data_write_f below.
!!
  CALL cgp_array_write_f('TimeValues',RealDouble,INT(1,cgsize_t),&
    INT(ggm_steps(num)%n_outsteps,cgsize_t),index_array,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!!  write(*,*) 'squish',self
!
  if (i-1 >= nprocs) then

    array_div = (i-1)/nprocs

    if (self == 0) then
      output_min_bound = 1
      output_max_bound = array_div
    else if (self == nprocs-1) then
      output_min_bound = self*array_div+1
      output_max_bound = i-1
    else
      output_min_bound = self*array_div+1
      output_max_bound = output_min_bound+array_div-1
    end if

!    write(*,*) 'bonk',self,output_min_bound,output_max_bound,array_div
!    write(*,*) 'useful data',self,ggm_steps(num)%saved_times(output_min_bound:output_max_bound)
    call cgp_array_write_data_f(index_array,INT(output_min_bound,cgsize_t),&
    INT(output_max_bound,cgsize_t),ggm_steps(num)%saved_times(output_min_bound:output_max_bound),ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'boof',self
  else
!    write(*,*) 'doing the other thing',self
    do j = 0,i-1
      if (self == j) then
        output_min_bound = self+1
        output_max_bound = output_min_bound
      end if
    end do
    if (self < j-1) then
      call cgp_array_write_data_f(index_array,INT(output_min_bound,cgsize_t),&
      INT(output_max_bound,cgsize_t),ggm_steps(num)%saved_times(output_min_bound:output_max_bound),ier)
      IF (ier /= CG_OK) CALL cgp_error_exit_f
!      write(*,*) 'doing the other thing',self
    end if
  end if

!  write(*,*) 'jingle',self
!
! If this is just a static multi-level run, we can output all the grid and zone information now.
!
  if (coarse_time < very_small_number .and. .not. amr) then
    if (.not. amr) then
      if (ggm_max_lvl(num) > fine) then
        out_max = fine
      else
        out_max = ggm_max_lvl(num)
      end if

      out_min = ggm_min_lvl(num)
!      write(*,*) 'zone storage bounds',self,out_min,out_max,index_file,num
      if (allocated(zone_storage)) deallocate(zone_storage)
      allocate(zone_storage(out_min:out_max))
!
!      do lvl = out_min,out_max
!        call zone_box_setup(lvl)
!      end do
      pass_num = num_saves+1
      !call zone_writer(out_min,out_max,num,pass_num)
      call zone_writer(out_min,out_max,num)
      call zone_coord_data_write(index_file,index_base,out_min,out_max,num)

    else
      call zone_dealloc_safety(fine,dummy,ggm_min_lvl(num),ggm_max_lvl(num))
    end if

  else
!    write(*,*) 'darm',self
!    if (ggm_min_lvl(num) < fine) then
!
!      return
!    end if

    do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
      allocate(ggm_save_info(lvl,num)% multigrid_save(ggm_steps(num)%n_outsteps))
    end do
    allocate(ggm_base_info(num)%total_zones_at_ts(ggm_steps(num)%n_outsteps))
    allocate(ggm_base_info(num)%zones_per_lvl(ggm_min_lvl(num):ggm_max_lvl(num),&
      ggm_steps(num)%n_outsteps))

    do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
      ggm_save_info(lvl,num)%multigrid_save(1:ggm_steps(num)%n_outsteps)%&
        regrid_between_outputs = .false.
    end do
!
! Regrid between ggm saves math
!
    do i = 1,ggm_steps(num)%n_outsteps-1
      do j = 1,size(regrid_times)
!        write(*,*) 'doof',size(regrid_times),ggm_steps(num)%saved_times(i),regrid_times(j)
! First iteration must be true to have output
        if (i == 1) then
          do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
!            if (lvl > fixed_lvl) then
              ggm_save_info(lvl,num)%multigrid_save(i)%regrid_between_outputs = .true.
!            end if
          end do

          exit
        else if (regrid_times(j) >= ggm_steps(num)%saved_times(i) .and. &
          regrid_times(j) <= ggm_steps(num)%saved_times(i+1)) then
! Loop through the levels to
          do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
            if (lvl > fixed_lvl) then
              ggm_save_info(lvl,num)%multigrid_save(i+1)%regrid_between_outputs = .true.
            end if
          end do

!          write(*,*) 'ggm regrid time flags',self,regrid_times(j),ggm_steps(num)%saved_times(i),&
!            ggm_steps(num)%saved_times(i+1),lvl,fixed_lvl
          exit

       else if (coarse_time >= ggm_steps(num)%saved_times(i) .and. &
          coarse_time <= ggm_steps(num)%saved_times(i+1)) then
!          write(*,*) 'goof'
          do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
            if (lvl > fixed_lvl) then
              ggm_save_info(lvl,num)%multigrid_save(i+1)%regrid_between_outputs = .true.
            end if
          end do
          exit
        else
          !write(*,*) 'bad timings',r
          cycle
        end if
      end do
    end do

    ggm_base_info(num)%adapted_zone_start_number = 1

!    write(*,*) 'ggm',self,ggm_save_info(ggm_max_lvl(num),num)%&
!      multigrid_save(1:ggm_steps(num)%n_outsteps)%regrid_between_outputs



  end if
!
!  write(*,*) 'after action report',self
!
!  CALL cgp_close_f(index_file,ier)
!  IF (ier /= CG_OK) CALL cgp_error_exit_f
!  WRITE(11,201) num
    ggm_steps(num)%output_counter = 1
!deallocate(temp_output_times)
!deallocate(temp_output_timesteps)

end do

end subroutine
!  if (save_times(1,num) < 0.0D0) then
!    time_start = 0.0D0
!  else
!    time_start = save_times(1,num)
!  end if

!!
!! Catchall for late starts
!!
!  if (save_times(1,num) > total_time) then
!    write(*,*) 'writing no data for save ',num
!    write(11,*) 'writing no data for save',num
!  end if
!
! Making sure of the total recording of the time
!
!  if (save_times(2,num) > total_time) then
!    time_stop = total_time
!  else
!    time_stop = save_times(2,num)
!  end if

!
! Number of saves
!
!  write(*,*) 'boop de snoot',self

!  ggm_steps(num)%n_outsteps = ceiling(t_diff/ggm_output_freq(num))
!  ggm_steps(num)%n_outsteps = ggm_steps(num)%n_outsteps+1
!  !write(*,*) 'number of output steps ',save_steps(num)%n_outsteps
!
!  if (allocated(ggm_steps(num)%saved_times)) deallocate(ggm_steps(num)%saved_times)
!  allocate(ggm_steps(num)%saved_times(ggm_steps(num)%n_outsteps))
!
!  if (allocated(ggm_steps(num)%saved_timestep)) deallocate(ggm_steps(num)%saved_timestep)
!  allocate(ggm_steps(num)%saved_timestep(ggm_steps(num)%n_outsteps))
!
!  if (allocated(ggm_steps(num)%sol_names)) deallocate(ggm_steps(num)%sol_names)
!  allocate(ggm_steps(num)%sol_names(ggm_steps(num)%n_outsteps))
!
!  if (allocated(ggm_steps(num)%write_done)) deallocate(ggm_steps(num)%write_done)
!  allocate(ggm_steps(num)%write_done(ggm_steps(num)%n_outsteps))


!        ts = (i-1)*dt(amrex_max_level)
!!    write(*,*) 'ts ', ts
!        if (ts > regrid_time(outly) .and. regrid_time(outly)-ts < dt(amrex_max_level) ) then
!          ggm_steps(num)%saved_times(i) = ts
!          ggm_steps(num)%saved_timestep(i) = FLOOR(ts/dt(ggm_max_lvl(num)))
!
!        end if
!
!        write(sol_num,format_junk_ts) ggm_steps(num)%saved_timestep(i)
!        ggm_steps(num)%sol_names(i) = 'Timestep_'//TRIM(sol_num)
!
!!    write(*,*) 'save_times, save_timesteps, solution names',&
!!      save_steps(num)%saved_times(i),save_steps(num)%saved_timestep(i),&
!!      save_steps(num)%sol_names(i),ts



!    do i = 1,num_regrids
!      write(sol_num,format_junk_ts) regrid_timesteps(i,ggm_max_lvl(num))
!      ggm_steps(num)%sol_names(i) = 'Timestep_'//TRIM(sol_num)
!
!      do lvl = 0,amrex_max_level
!
!        ggm_steps(num)%saved_timestep(i) = regrid_timesteps(i,ggm_max_lvl(num))
!
!      end do
!    end do

!        i = regrid_timesteps(j,amrex_max_level)
!        outly = 1

!        if (abs(regrid_times(j) - ts) < ggm_dt/2.0D0) then
!          temp_output_times(i) = regrid_times(j)
!          do lvl = amrex_max_level,0,-1
!
!            temp_output_timesteps(i,lvl) = regrid_timesteps(j,lvl)
!          end do
!          valid_ggm_outs = valid_ggm_outs+1
!          exit
!        else
!          cycle
!        end if

!      end do

!  if (ggm_output_freq(num) < ggm_dt) then
!    ggm_steps(num)%n_outsteps = num_regrids
!    if (allocated(ggm_steps(num)%saved_times)) deallocate(ggm_steps(num)%saved_times)
!    allocate(ggm_steps(num)%saved_times(num_regrids))
!    if (allocated(ggm_steps(num)%saved_timestep)) deallocate(ggm_steps(num)%saved_timestep)
!    allocate(ggm_steps(num)%saved_timestep(num_regrids))
!    if (allocated(ggm_steps(num)%sol_names)) deallocate(ggm_steps(num)%sol_names)
!    allocate(ggm_steps(num)%sol_names(num_regrids))
!    if (allocated(ggm_steps(num)%write_done)) deallocate(ggm_steps(num)%write_done)
!    allocate(ggm_steps(num)%write_done(num_regrids))
!
!    ggm_output_freq(num) = ggm_dt
!    ggm_steps(num)%saved_times = regrid_times
!    ggm_steps(num)%write_done = .false.
!
!    valid_ggm_outs = num_regrids
!    ggm_steps(num)%n_outsteps = valid_ggm_outs

!  else


!      do lvl = amrex_max_level,0,-1
!        ggm_steps(num)%saved_timestep( = floor(ts/ggm_max_lvl(lvl))
!      end do

!      valid_ggm_outs = valid_ggm_outs + 1
!      i = i+1
!      write(*,*) 'timesteps out',ts
!      ts = i*ggm_output_freq(num)

!    do i = 1,valid_ggm_outs
!      ggm_steps(num)%saved_times(i) = temp_output_times(i)
!!      do lvl = 0,amrex_max_level
!        ggm_steps(num)%saved_timestep(i) = temp_output_timesteps(i,ggm_max_lvl(num))
!!      end do
!!
!      write(sol_num,format_junk_ts) ggm_steps(num)%saved_timestep(i)
!      ggm_steps(num)%sol_names(i) = 'Timestep_'//TRIM(sol_num)
!!
!    end do

!    ggm_steps(num)%n_outsteps = valid_ggm_outs
!
!  end if
!  WRITE(char_save_number,format_junk) num
