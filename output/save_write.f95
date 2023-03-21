subroutine save_write(num,fine)
!
! Write all required information to the output file.
!
!
! Called by: lb_output
! Calls:
! External calls:
!
use precise
use amrex_amr_module
use amrex_base_module
use zone_writing
use output_data
use nml_output
use grid_data
use linkwise
use amr_grid_writing
use timing
use constants
use freestream_values
use cgns
use mpi
use amr_info_holder
use amr_processes, only: self,commune,count_big_zones,total_big_zones
implicit none
INTEGER :: index_base,index_zone,index_coord,index_file,icelldim,&
  iphysdim,ier,index_flow,index_field,index_section,dummy
integer :: fine,max_lvl,lvl,num,zee_count,total_zones,zer,lvl_zones
integer,allocatable :: level_zones(:)
integer*8 :: loc_lo(4),loc_hi(4)
CHARACTER(len=32) :: format_junk,char_save_number,format_junk_ts,filename
logical :: regridded

type(amrex_mfiter) :: mfi
type(amrex_box)  :: bux

if (time(fine) < save_times(1,num) .or. time(fine) > (save_times(2,num))) RETURN

write(*,*) 'Preparing to write output step',time(fine),self

index_base = 1
index_file = num+1
index_coord = 1
index_zone = 1
index_section = 1
index_flow = 1
index_field = 1
dummy = 1

format_junk = '(I3.3)'
format_junk_ts = '(I9.9)'

if (fine < save_shapes(1,num)) then
  WRITE(11,*) 'No appropriate levels to output for save ',num, &
    'at time ',time(fine)
  write(11,*) 'The finest level currently in the mesh ',fine, &
    'is smaller than the desired level of ',save_shapes(1,num)
  return
end if

IF (dimensions == 3) THEN
  icelldim = 3
  iphysdim = 3
ELSE
  icelldim = 2
  iphysdim = 3
END IF
!
!
!
!
!
WRITE(char_save_number,format_junk) num
filename = 'output_'//TRIM(char_save_number)//'.cgns'
!basename = 'outputbase_'//TRIM(char_save_number)
!
!  if (amr) then
if (fine > save_shapes(2,num)) then
  max_lvl = save_shapes(2,num)
else
  max_lvl = fine
end if

!write(*,*) 'max/fine',max_lvl,fine


!
!
! Set up the zone outputs
!
!
if (coarse_time < very_small_number .and. .not. amr) then
!

  if (allocated(zone_storage)) deallocate(zone_storage)
  allocate(zone_storage(save_shapes(1,num):max_lvl))
  do lvl = save_shapes(1,num),max_lvl
!    write(*,*) 'box level bounds',self,save_shapes(1,num),max_lvl
    call zone_box_setup(lvl)
  end do
  call zone_output_num_setup(num,save_shapes(1,num),max_lvl)
!
!
!
!  CALL cgp_open_f(TRIM(filename),CG_MODE_MODIFY,index_file,ier)
  !
  ! Write the solution name out
  !
  do lvl = save_shapes(1,num),max_lvl
    do index_zone = 1,zone_storage(lvl)%big_zones
  !    write(*,*) 'Before solution write',self,index_file,index_base,&
  !      zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%sol_names(save_steps(num)%output_counter),&
  !      save_steps(num)%output_counter,index_flow

  !    call cg_error_print_f
      call cg_sol_write_f(index_file,index_base,zone_storage(lvl)%output_number(index_zone,num),&
        save_steps(num)%sol_names(save_steps(num)%output_counter),Vertex,&
        save_steps(num)%output_counter,ier)

  !    call cg_error_print_f
  !    write(*,*) 'after solution',self,index_file,index_base,&
  !      zone_storage(lvl)%output_number(index_zone,1),index_flow
      IF (ier /= CG_OK) CALL cgp_error_exit_f
    end do
  end do
  zee_count = 0
  do lvl = save_shapes(1,num),max_lvl

    do index_zone = 1,zone_storage(lvl)%big_zones

      index_field = 1
  !    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
  !      zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(1,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'Density',index_field,ier)
        index_field = index_field+1
      end if
  !    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
  !      zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(2,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'U-velocity',index_field,ier)
        index_field = index_field+1
      end if
  !    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
  !      zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(3,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'V-velocity',index_field,ier)
        index_field = index_field+1
      end if
  !    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
  !      zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(4,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'W-velocity',index_field,ier)
        index_field = index_field+1
      end if
  !    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
  !      zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(5,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'Pressure',index_field,ier)
        index_field = index_field+1
      end if
      !write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !  zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(6,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'Temperature',index_field,ier)
        index_field = index_field+1
      end if
      !write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !  zone_storage(lvl)%output_number(index_zone,num),self
      if (save_data(7,num)) then
        call cgp_field_write_f(index_file,index_base,&
          zone_storage(lvl)%output_number(index_zone,num),save_steps(num)%output_counter,&
          RealDouble,'Alpha',index_field,ier)
        index_field = index_field+1
      end if

      IF (ier /= CG_OK) CALL cgp_error_exit_f
    end do
  end do

  do lvl = save_shapes(1,num),max_lvl
  zee_count = zone_storage(lvl)%zone_start(self)
    call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
    do while (mfi%next())
      bux = mfi%validbox()
  !    if (zone_storage(lvl)%full(b)) then
  !        write(*,*) self,zone_storage(lvl)%full,c
  !        write(*,*) self,zone_storage(lvl)%zones_owned
  !        write(*,*) 'bounds of rho',lbound(rho),ubound(rho)
  !        write(*,*) 'bounds of the box',bux%lo,bux%hi
  !        write(*,*) 'rho at low',rho(bux%lo(1),bux%lo(2),bux%lo(3),1)
  !        write(*,*) 'rho at hi',rho(bux%hi(1),bux%hi(2),bux%hi(3),1)
      loc_lo = 1
      loc_hi(1) = bux%hi(1)-bux%lo(1)+1
      loc_hi(2) = bux%hi(2)-bux%lo(2)+1
      loc_hi(3) = bux%hi(3)-bux%lo(3)+1
      loc_hi(4) = 1
  !
  !
  !
  !
      index_field = 1
  !        call cgp_field_write_data_f(index_file,index_base,&
  !          zone_storage(lvl)%zones_owned(c),&
  !          index_flow,index_field,loc_lo,loc_hi,rho(bux%lo(1):bux%hi(1),&
  !          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      if (save_data(1,num)) then
        rho => mfrho(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, rho',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,rho(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if
      if (save_data(2,num)) then
        u_vel => mfu_vel(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, u',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,u_vel(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if
      if (save_data(3,num)) then
        v_vel => mfv_vel(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, v',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,v_vel(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if
      if (save_data(4,num)) then
        w_vel => mfw_vel(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, w',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,w_vel(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if
      if (save_data(5,num)) then
        press => mfpress(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, pressure',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,press(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if
      if (save_data(6,num)) then
        temp => mftemp(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, temp',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,temp(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if

      if (save_data(7,num)) then
        alpha => mfalpha(lvl)%dataptr(mfi)
        !write(*,*) 'field data write, temp',index_file,index_base,save_steps(num)%output_counter,index_field,&
        !zone_storage(lvl)%output_number(zee_count,num),self
        call cgp_field_write_data_f(index_file,index_base,&
          zone_storage(lvl)%output_number(zee_count,num),&
          save_steps(num)%output_counter,index_field,loc_lo,loc_hi,alpha(bux%lo(1):bux%hi(1),&
          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
        index_field = index_field+1
      end if
  !      write(*,*) 'field data',self,index_file,index_base,b,&
  !        index_flow,index_field,loc_lo,loc_hi
      IF (ier /= CG_OK) CALL cgp_error_exit_f
      zee_count = zee_count+1
    end do
  end do

!  call cgp_close_f(index_file,ier)

  call zone_dealloc_safety(fine,num,save_shapes(1,num),max_lvl)

  save_steps(num)%write_done(save_steps(num)%output_counter) = .true.
  save_steps(num)%output_counter = save_steps(num)%output_counter + 1
!
!
!
!
!
! For any kind of adapting mesh
!
!
!
!
!
else
!  CALL cgp_open_f(TRIM(filename),CG_MODE_MODIFY,index_file,ier)
  !
  ! Write the solution name out
  !
  allocate(level_zones(save_shapes(1,num):save_shapes(2,num)))
!  write(*,*) 'cgns file opened',filename,index_file

!
! First timestep to establish some things
!
  if (save_steps(num)%output_counter == 1) then
!    write(*,*) 'wedge antilles'
!    call cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',1,'TimeValues',&
!      save_steps(num)%n_outsteps,'end')

    total_zones = 0

    base_info(num)%adapted_zone_start_number = 1

    do lvl = save_shapes(1,num),max_lvl
! If the requested output level is finer than the current finest, skip it
      if (fine < save_shapes(1,num)) then
        cycle
! Else, copy copy from zone_box_setup to find total number of zones for this iteration
      else
        lvl_zones = total_big_zones(lvl)
        !write(*,*) 'lvl_zones',lvl_zones

        total_zones = total_zones + lvl_zones
        base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) = lvl_zones!grid_info(lvl,num)%big_zones

        allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(lvl_zones))
        allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_lvls(lvl_zones))

      end if


    end do
    base_info(num)%max_zones = total_zones
!
! If nothing has changed, use the previous value
!
  else
    regridded = .false.
    total_zones = 0
    do lvl = save_shapes(1,num),max_lvl
!
      if (grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%regrid_between_outputs) then
!
        lvl_zones = total_big_zones(lvl)
!
        total_zones = total_zones + lvl_zones
!        write(*,*) 'inside box counting',lvl_zones,total_zones,lvl,self
        base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) = lvl_zones
! Otherwise, the number of zones on a level from a previous timestep can be used.

        regridded = .true.

        allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(lvl_zones))
        allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_lvls(lvl_zones))

      else
        total_zones = base_info(num)%total_zones_at_ts(save_steps(num)%output_counter - 1)

        base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) = &
          base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter -1)

        allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%&
          zone_number(base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)))
!
        allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%&
          zone_lvls(base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)))

      end if

    end do

    if (total_zones > base_info(num)%max_zones) then
      base_info(num)%max_zones = total_zones
    end if
  end if
!  write(*,*) 'jangle',total_zones
!
!
!
  if (total_zones > 0) then
    base_info(num)%total_zones_at_ts(save_steps(num)%output_counter) = total_zones
!
!    do lvl = save_shapes(1,num),save_shapes(2,num)
!      allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(total_zones))
!      allocate(grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_lvls(total_zones))
!    end do
!    write(*,*) 'crinkle',self
!
! Go through all the levels
!
    do lvl = save_shapes(1,num),max_lvl
! If there's nothing at this level, skip
      if (base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) == 0) cycle
!
! If regrid_between_outputs is true at the given level, output new zones for affected levels
!
      if (grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%regrid_between_outputs) then
!        write(*,*) 'BIFF!!!'
!        call cg_goto_f(index_file,index_base,'Zone_t',base_info(num)%adapted_zone_start_number,&
!          'end',ier)
!        write(*,*) 'zone_starts',base_info(num)%adapted_zone_start_number
        call amr_grid_zones(lvl,num)
        base_info(num)%adapted_zone_start_number = base_info(num)%adapted_zone_start_number + &
             grid_info(lvl,num)%big_zones
!        write(*,*) 'base info stuff',base_info(num)%adapted_zone_start_number
        call amr_grid_zone_writer(lvl,num,index_file)
        call amr_grid_coord_data_write(lvl,num,index_file)
!    do lvl = save_shapes(1,num),max_lvl
!      allocate(grid_info(lvl)%multigrid_save(save_steps(num)%output_counter)% &
!        zones_used(zone_storage(lvl)%little_zones))

!    end do
!
! else use the same information as last time
!
      else
        do zer = 1,base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)
!          write(*,*) 'BOW!!!'

          if (save_steps(num)%output_counter == 1) then
            call amr_grid_zones(lvl,num)
!            base_info(num)%adapted_zone_start_number = base_info(num)%adapted_zone_start_number + &
!              grid_info(lvl,num)%big_zones

            call amr_grid_zone_writer(lvl,num,index_file)
            call amr_grid_coord_data_write(lvl,num,index_file)
!          else
!            write(*,*) 'giggle'
!            grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(zer) = &
!              grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter - 1)%zone_number(zer)
!            grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_lvls(zer) = &
!              grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter - 1)%zone_lvls(zer)
!
!            write(*,*) 'zone nums',self,lvl,zer,grid_info(lvl,num)%&
!              multigrid_save(save_steps(num)%output_counter)%zone_number(zer)
          end if
        end do
      end if

    end do

!
! Do the actual save writes
!
    do lvl = save_shapes(1,num),max_lvl
! If there's nothing at this level, skip anything
      if (base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) == 0) cycle
!
!      write(*,*) 'Thor, puppy of thunder.'
      call amr_save_writer(lvl,num,fine,index_file,index_base)
      call mpi_barrier(commune,ier)
    end do



  end if
!  call mpi_barrier(commune,ier)
!  write(*,*) 'effigy',self
!  call cgp_close_f(index_file,ier)

  write(*,*) 'blargeist',self

  save_steps(num)%write_done(save_steps(num)%output_counter) = .true.
  save_steps(num)%output_counter = save_steps(num)%output_counter + 1
end if

if (allocated(level_zones)) deallocate(level_zones)
!call mpi_barrier(commune,ier)

end subroutine
!write(*,*) 'Save successfully written.'
!write(11,*) 'Save successfully written.'

!write(*,*) 'CGNS file opened to write',num,index_file,index_base,index_flow,filename
!write(*,*) 'CGNS file opened to write',num,index_file,index_base,index_flow

!
!write(sol_num,format_junk) save_steps(num)%saved_timestep(save_steps(num)%output_counter)
!solut_name = 'Timestep_'//trim(sol_num)







!
! If any levels have been regridded
!
!  else if (any(grid_info(0:amrex_max_level,num)%multigrid_save(save_steps(num)%output_counter)%&
!    regrid_between_outputs)) then
!
!    total_zones = 0
!
!    do lvl = save_shapes(1,num),save_shapes(2,num)
!! If a level has been regridded, we need to count the zones again
!      if (grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%&
!          regrid_between_outputs) then
!
!        lvl_zones = count_big_zones(lvl)
!
!        total_zones = total_zones + lvl_zones
!        base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) = lvl_zones
!! Otherwise, the number of zones on a level from a previous timestep can be used.
!      else
!        total_zones = total_zones + base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter - 1)
!        base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter) = &
!               base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter - 1)
!      end if
!
!      if (total_zones > base_info(num)%max_zones) then
!        base_info(num)%max_zones = total_zones
!      end if
!
!    end do
