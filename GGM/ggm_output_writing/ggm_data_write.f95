subroutine ggm_data_write(num,fine)
!
! Output GGM data to a CGNS file
!
!
!
!
!
!
!
!
use ggm_stuff
use amr_info_holder
use amr_processes
use cgns
use mpi
use nml_output, only: num_saves
use precise
use timing
use constants
use output_data, only: ggm_steps,ggm_base_info
use nml_output
use ggm_amr_writing
use amrex_amr_module
use amrex_base_module
use zone_writing
implicit none
!real(kind=dp) ::
INTEGER :: index_base,index_zone,index_coord,index_file,icelldim,&
  iphysdim,ier,index_flow,index_field,index_section
integer,allocatable :: level_zones(:)
integer :: fine,max_lvl,lvl,ggm_num,num,zee_count,total_zones,lvl_zones,zer
integer*8 :: loc_lo(4),loc_hi(4)
logical :: regridded
CHARACTER(len=32) :: format_junk,char_save_number,format_junk_ts,ggm_field_value,filename
!real(kind=dp) :: cur_t

type(amrex_mfiter) :: mfi
type(amrex_box)  :: flux
!
!
!
!
!
ggm_num = num_saves + num

index_base = 1
index_file = ggm_num + 1
index_coord = 1
index_zone = 1
index_section = 1
index_flow = 1
index_field = 1

format_junk = '(I3.3)'
format_junk_ts = '(I9.9)'
!
!
!
IF (dimensions == 3) THEN
  icelldim = 3
  iphysdim = 3
ELSE
  icelldim = 2
  iphysdim = 3
END IF
!
!write(*,*) 'CGNS file opened to write',num,index_file,index_base,index_flow
!
WRITE(char_save_number,format_junk) num
filename = 'ggm_out_'//TRIM(char_save_number)//'.cgns'
!basename = 'outputbase_'//TRIM(char_save_number)
!
!  if (amr) then
if (fine > ggm_max_lvl(num)) then
  max_lvl = ggm_max_lvl(num)
else
  max_lvl = fine
end if
!

!
if (coarse_time < very_small_number .and. .not. amr) then
!
  if (allocated(zone_storage)) deallocate(zone_storage)
  allocate(zone_storage(ggm_min_lvl(num):max_lvl))

  do lvl = ggm_min_lvl(num),max_lvl
  !      write(*,*) 'box level bounds',self,save_shapes(1,num),max_lvl
    call zone_box_setup(lvl)
  end do
  call zone_output_num_setup(ggm_num,ggm_min_lvl(num),max_lvl)
  !  end if
  !
  ! XXXX check on this to be sure
  !
  !
  !
  !write(*,*) 'CGNS file opened to write',num,index_file,index_base,index_flow
!  CALL cgp_open_f(TRIM(filename),CG_MODE_MODIFY,index_file,ier)
  !
  ! Write the solution name out
  !
  !write(sol_num,format_junk) ggm_steps(num)%saved_timestep(ggm_steps(num)%output_counter)
  !solut_name = 'Timestep_'//trim(sol_num)
  zee_count = 0
  do lvl = ggm_min_lvl(num),max_lvl
    do index_zone = 1,zone_storage(lvl)%big_zones

  !    write(*,*) 'Before solution write',self,index_file,index_base,&
  !      zone_storage(lvl)%output_number(index_zone,num),&
  !      ggm_steps(num)%sol_names(ggm_steps(num)%output_counter),&
  !      ggm_steps(num)%output_counter,index_flow
  !    call cg_error_print_f
      call cg_sol_write_f(index_file,index_base,zone_storage(lvl)%output_number(index_zone,num),&
        ggm_steps(num)%sol_names(ggm_steps(num)%output_counter),Vertex,&
        ggm_steps(num)%output_counter,ier)
  !    call cg_error_print_f
  !    write(*,*) 'after solution',self,index_file,index_base,&
  !      zone_storage(lvl)%output_number(index_zone,1),index_flow
      IF (ier /= CG_OK) CALL cgp_error_exit_f
    end do
  end do
  !
  !
  !
  do lvl = ggm_min_lvl(num),max_lvl
    do index_zone = 1,zone_storage(lvl)%big_zones
  !
  ! Because there are a lot of possibilities, just have one select instead of 20
  !
      select case (ggm_variables(num))
        case (1)
          ggm_field_value = 'GGM Density'
        case (2)
          ggm_field_value = 'GGM u-Velocity'
        case (3)
          ggm_field_value = 'GGM v-Velocity'
        case (4)
          ggm_field_value = 'GGM w-Velocity'
        case (5)
          ggm_field_value = 'GGM Temperature'
        case (6)
          ggm_field_value = 'GGM Pressure'
        case (7)
          ggm_field_value = 'GGM Velocity Magnitude'
        case (8)
          ggm_field_value = 'GGM Dynamic Pressure'
        case (9)
          ggm_field_value = 'GGM Mach Number'
        case (10)
          ggm_field_value = 'GGM xy-Vorticity'
        case (11)
          ggm_field_value = 'GGM xz-Vorticity'
        case (12)
          ggm_field_value = 'GGM yz-Vorticity'
        case (13)
          ggm_field_value = 'GGM Vorticity Magnitude'
        case (14)
          ggm_field_value = 'GGM xy-Shear Stress'
        case (15)
          ggm_field_value = 'GGM xz-Shear Stress'
        case (16)
          ggm_field_value = 'GGM yz-Shear Stress'
        case (17)
          ggm_field_value = 'GGM Shear Stress Magnitude'
        case (20)
          ggm_field_value = 'GGM Alpha'
        case (29)
          ggm_field_value = 'GGM Q-criterion'

      end select

      index_field = 1

      call cgp_field_write_f(index_file,index_base,&
        zone_storage(lvl)%output_number(index_zone,num),ggm_steps(num)%output_counter,&
        RealDouble,ggm_field_value,index_field,ier)

  !    write(*,*) 'after field write',index_file,index_base,ggm_field_value,index_field

    end do
  end do
  !
  !
  !
  do lvl = ggm_min_lvl(num),max_lvl
    zee_count = zone_storage(lvl)%zone_start(self)
    call amrex_mfiter_build(mfi,mfggm(lvl),tiling=.false.)
    do while (mfi%next())
      flux = mfi%validbox()
  !    if (zone_storage(lvl)%full(b)) then
  !        write(*,*) self,zone_storage(lvl)%full,c
  !        write(*,*) self,zone_storage(lvl)%zones_owned
  !        write(*,*) 'bounds of rho',lbound(rho),ubound(rho)
  !        write(*,*) 'bounds of the box',bux%lo,bux%hi
  !        write(*,*) 'rho at low',rho(bux%lo(1),bux%lo(2),bux%lo(3),1)
  !        write(*,*) 'rho at hi',rho(bux%hi(1),bux%hi(2),bux%hi(3),1)
      loc_lo = 1
      loc_hi(1) = flux%hi(1)-flux%lo(1)+1
      loc_hi(2) = flux%hi(2)-flux%lo(2)+1
      loc_hi(3) = flux%hi(3)-flux%lo(3)+1
      loc_hi(4) = 1

      index_field = 1

  !    ggm => mfggm(lvl)%dataptr(mfi)
      ggm => mfggm(lvl)%dataptr(mfi)
      call cgp_field_write_data_f(index_file,index_base,&
        zone_storage(lvl)%output_number(zee_count,num),&
        ggm_steps(num)%output_counter,index_field,loc_lo,loc_hi,ggm(flux%lo(1):flux%hi(1),&
        flux%lo(2):flux%hi(2),flux%lo(3):flux%hi(3),num:num),ier)

  !    write(*,*) 'after field data write',index_file,index_base,index_field

      IF (ier /= CG_OK) CALL cgp_error_exit_f
      zee_count = zee_count+1
    end do
  end do

!  call cgp_close_f(index_file,ier)

  ggm_steps(num)%write_done(ggm_steps(num)%output_counter) = .true.
  ggm_steps(num)%output_counter = ggm_steps(num)%output_counter + 1

  write(11,*) 'GGM save successfully written.'
  write(*,*) 'GGM save successfully written.'

  call zone_dealloc_safety(fine,num,ggm_min_lvl(num),max_lvl)
!
!
!
!
! Adapting grid stuff
!
!
!
!
else
!  CALL cgp_open_f(TRIM(filename),CG_MODE_MODIFY,index_file,ier)
  !
  ! Write the solution name out
  !
!  write(*,*) 'opening for ggm out'
  allocate(level_zones(ggm_min_lvl(num):ggm_max_lvl(num)))


  if (ggm_steps(num)%output_counter == 1) then

!    call cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',1,'TimeValues',&
!      ggm_steps(num)%n_outsteps,'end')

    total_zones = 0

    ggm_base_info(num)%adapted_zone_start_number = 1

    do lvl = ggm_min_lvl(num),max_lvl
! If the requested output level is finer than the current finest, skip it
      if (fine < ggm_min_lvl(num)) then
        cycle
! Else, copy copy from zone_box_setup to find total number of zones for this iteration
      else
        lvl_zones = total_big_zones(lvl)

        total_zones = total_zones + lvl_zones
        ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter) = lvl_zones

      end if


    end do
    ggm_base_info(num)%max_zones = total_zones
!
! If nothing has changed, use the previous value
!
  else
    regridded = .false.
    total_zones = 0
    do lvl = ggm_min_lvl(num),max_lvl

      if (ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%regrid_between_outputs) then

        lvl_zones = total_big_zones(lvl)

        total_zones = total_zones + lvl_zones
        ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter) = lvl_zones
! Otherwise, the number of zones on a level from a previous timestep can be used.

        regridded = .true.

      else
        total_zones = ggm_base_info(num)%total_zones_at_ts(ggm_steps(num)%output_counter - 1)

        ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter) = &
          ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter -1)

      end if
    end do

    if (total_zones > ggm_base_info(num)%max_zones) then
      ggm_base_info(num)%max_zones = total_zones
    end if
  end if
!
!
!
  if (total_zones > 0) then
    ggm_base_info(num)%total_zones_at_ts(ggm_steps(num)%output_counter) = total_zones
!
!
!
    write(*,*) 'jingle',total_zones,ggm_steps(num)%output_counter
    do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
      allocate(ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_number(total_zones))
      allocate(ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_lvls(total_zones))
    end do
!
! Go through all the levels
!
  do lvl = ggm_min_lvl(num),max_lvl
! If there's nothing at this level, skip
    if (ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter) == 0) cycle
!
! If regrid_between_outputs is true at the given level, output new zones for affected levels
!
    if (ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%regrid_between_outputs) then

      call ggm_amr_grid_zones(lvl,num)
      ggm_base_info(num)%adapted_zone_start_number = ggm_base_info(num)%adapted_zone_start_number + &
          ggm_save_info(lvl,num)%big_zones

      call ggm_amr_grid_zone_writer(lvl,num,index_file)
      call ggm_amr_grid_coord_data_write(lvl,num,index_file)
!    do lvl = save_shapes(1,num),max_lvl
!      allocate(ggm_save_info(lvl)%multigrid_save(ggm_steps(num)%output_counter)% &
!        zones_used(zone_storage(lvl)%little_zones))

!    end do
!
! else use the same information as last time
!
    else
      do zer = 1,ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter)


        if (ggm_steps(num)%output_counter == 1) then
          call ggm_amr_grid_zones(lvl,num)
!          ggm_base_info(num)%adapted_zone_start_number = ggm_base_info(num)%adapted_zone_start_number + &
!            ggm_save_info(lvl,num)%big_zones

          call ggm_amr_grid_zone_writer(lvl,num,index_file)
          call ggm_amr_grid_coord_data_write(lvl,num,index_file)
        else
          ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_number(zer) = &
            ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter - 1)%zone_number(zer)
          ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_lvls(zer) = &
            ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter - 1)%zone_lvls(zer)
        end if
      end do
    end if

  end do
!
! Do the actual save writes
!
    do lvl = ggm_min_lvl(num),max_lvl
! If there's nothing at this level, skip anything
      if (ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter) == 0) cycle
!
      call ggm_amr_save_writer(lvl,num,fine,index_file,index_base)
    end do

  end if

!  write(*,*) 'closed ggm data file',index_file,ier,ggm_steps(num)%write_done(ggm_steps(num)%output_counter),&
!    ggm_steps(num)%output_counter
!
!  call mpi_barrier(commune,ier)

!  call cgp_close_f(index_file,ier)
!  call cg_error_print_f()
  write(*,*) 'scooby doo!'

  ggm_steps(num)%write_done(ggm_steps(num)%output_counter) = .true.
  ggm_steps(num)%output_counter = ggm_steps(num)%output_counter + 1

end if

end subroutine
!write(*,*) 'CGNS file opened to write',num,index_file,index_base,index_flow,filename
!
!
!
!if (fine < ggm_max_lvl(num)) then
!  max_lvl = fine
!else
!  max_lvl = ggm_max_lvl(num)
!end if

!if (fine < save_shapes(1,num)) then
!  WRITE(11,*) 'No appropriate levels to output for save ',num, &
!    'at time ',time(fine)
!  write(11,*) 'The finest level currently in the mesh ',fine, &
!    'is smaller than the desired level of ',save_shapes(1,num)
!  return
!end if
