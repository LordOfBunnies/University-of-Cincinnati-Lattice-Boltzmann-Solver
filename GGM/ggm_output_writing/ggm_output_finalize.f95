subroutine ggm_output_finalize
!
! finalize the GGM output
!
!
!
!
! Called by: output_finalize
! Calls: error_out,external_code_finalizer,cpu_time
! External calls: cgp_open_f,cgp_error_exit_f,cg_goto_f,cg_biter_write_f,
!   cgp_close_f,cg_ziter_write_f,cgp_array_write_f,cgp_array_write_data_f
!
use precise
use constants
use ggm_stuff
USE grid_data
USE output_data, only: ggm_steps,ggm_base_info
use amr_info_holder
use amr_processes, only: self,nprocs,commune
use cgns
use mpi
use timing
USE linkwise
use nml_output, only: num_saves
use amrex_amr_module
use amrex_base_module
use zone_writing
use ggm_amr_writing
implicit none
INTEGER :: index_base,index_zone,index_coord,index_file,icelldim,&
  iphysdim,ier,index_flow,index_field,index_section,index_array
integer :: l,fine,max_lvl,lvl,num,zee_count,i,j,array_div,z,ts
INTEGER(cgsize_t) :: idata(2),output_min_bound,output_max_bound,dummy_bound,min_dummy
logical :: marked
CHARACTER(len=32) :: format_junk,char_save_number,format_junk_ts,&
  filename,ggm_field_value,basename,format_nm_lvl_and_zone,zonename

!real(kind=dp) :: cur_t

!type(amrex_mfiter) :: mfi
!type(amrex_box)  :: flux
!
!
!
!
!
format_junk = '(I3.3)'
format_junk_ts = '(I9.9)'
fine = amrex_get_finest_level()

!do num = 1,num_ggm_vars
!  index_file = num + num_saves + 1
!  index_base = 1
!  index_coord = 1
!  index_zone = 1
!  index_section = 1
!  index_flow = 1
!  index_field = 1
!  ier =0
!
!  idata(1) = 32
!  idata(2) = ggm_steps(num)%n_outsteps
!
!  IF (dimensions == 3) THEN
!    icelldim = 3
!    iphysdim = 3
!  ELSE
!    icelldim = 2
!    iphysdim = 3
!  END IF
!
!  WRITE(char_save_number,format_junk) num
!  filename = 'ggm_out_'//TRIM(char_save_number)//'.cgns'
!  !basename = 'outputbase_'//TRIM(char_save_number)
!  !
!  !  if (amr) then
!  if (fine > max_lvl) then
!    max_lvl = ggm_max_lvl(num)
!  else
!    max_lvl = fine
!  end if
!!
!!
!!
!!
!!
!!
! if (coarse_time < very_small_number .and. .not. amr) then
!  !
!  if (allocated(zone_storage)) deallocate(zone_storage)
!  allocate(zone_storage(ggm_min_lvl(num):max_lvl))
!  !
!  do lvl = ggm_min_lvl(num),max_lvl
!  !      write(*,*) 'box level bounds',self,save_shapes(1,num),max_lvl
!    call zone_box_setup(lvl)
!  end do
!  call zone_output_num_setup(num,ggm_min_lvl(num),max_lvl)
!
!
!  WRITE(char_save_number,format_junk) num
!  filename = 'ggm_out_'//TRIM(char_save_number)//'.cgns'
!  basename = 'ggm_base_'//TRIM(char_save_number)
!
!
!!  if (fine >= save_shapes(2,a)) then
!!    sol_max_lvl = save_shapes(2,a)
!!  else
!!    sol_max_lvl = fine
!!  end if
!
!!  CALL cgp_open_f(filename,CG_MODE_MODIFY,index_file,ier)
!!  IF (ier /= CG_OK) CALL cgp_error_exit_f
!  WRITE(11,*) 'Output file reopened for finalization.'
!
!  CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',&
!    1,'end')
!  IF (ier /= CG_OK) CALL cgp_error_exit_f
!!    write(*,*) lbound(all_sol_names),ubound(all_sol_names)
!  index_array = 2
!
!    do lvl = ggm_min_lvl(num),ggm_max_lvl(num)
!      do i = zone_storage(lvl)%output_number(1,num),&
!        zone_storage(lvl)%output_number(zone_storage(lvl)%big_zones,num)
!
!!        write(*,*) 'cross check before output',index_file,index_base,i
!
!        CALL cg_ziter_write_f(index_file,index_base,i,&
!          'ZoneIterativeData',ier)
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!! Go to the zone data
!!
!        CALL cg_goto_f(index_file,index_base,ier,'Zone_t',i,&
!          'ZoneIterativeData_t',1,'end')
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!! Tell it the solution names and their sizes
!!
!!        write(*,*) 'solution names, making sure they match',self,all_sol_names
!        CALL cgp_array_write_f('FlowSolutionPointers',Character,int(2,cgsize_t),idata,index_array,ier)
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!!        write(*,*) 'squishy face!!!',index_array,idata
!
!        if (ggm_steps(num)%n_outsteps >= nprocs) then
!
!          array_div = (ggm_steps(num)%n_outsteps-1)/nprocs
!
!          if (self == 0) then
!            output_min_bound = 1
!            output_max_bound = array_div
!          else if (self ==nprocs-1) then
!            output_min_bound = self*array_div+1
!            output_max_bound = ggm_steps(num)%n_outsteps
!          else
!            output_min_bound = self*array_div+1
!            output_max_bound = output_min_bound+array_div-1
!          end if
!
!            call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!              (/idata(1),output_max_bound/),ggm_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!            IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!        else
!          do j = 1,ggm_steps(num)%n_outsteps
!            if (self <= j) then
!              output_min_bound = self+1
!              output_max_bound = output_min_bound+1
!            end if
!          end do
!          if (self <= ggm_steps(num)%n_outsteps) then
!            call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!              (/idata(1),output_max_bound/), ggm_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!            IF (ier /= CG_OK) CALL cgp_error_exit_f
!          end if
!        end if
!
!      end do
!    end do
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
!  else
!!    call cgp_open_f(filename,CG_MODE_MODIFY,index_file,ier)
!!    if (ier /= CG_OK) call cgp_error_exit_f
!    write(11,*) 'Output file reopened for finalization.'
!!
!!
!!
!    if (ggm_steps(num)%n_outsteps >= nprocs) then
!
!      array_div = (ggm_steps(num)%n_outsteps-1)/nprocs
!      if (self == 0) then
!        output_min_bound = 1
!        output_max_bound = array_div
!      else if (self ==nprocs-1) then
!        output_min_bound = self*array_div+1
!        output_max_bound = ggm_steps(num)%n_outsteps
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
!      do j = 1,ggm_steps(num)%n_outsteps
!        if (self <= j) then
!          output_min_bound = self+1
!          output_max_bound = output_min_bound+1
!        end if
!      end do
!      if (self <= ggm_steps(num)%n_outsteps) then
!        call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/),&
!          (/idata(1),output_max_bound/), ggm_steps(num)%sol_names(output_min_bound:output_max_bound),ier)
!        IF (ier /= CG_OK) CALL cgp_error_exit_f
!      end if
!    end if
!!
!! zone_pointer_names contains zone names for the BaseIterativeData_t
!! amr_solution_names contains names for solutions
!!
!    allocate(ggm_base_info(num)%zone_pointer_names(ggm_base_info(num)%max_zones,&
!      ggm_steps(num)%n_outsteps))
!    allocate(ggm_base_info(num)%amr_solution_names(ggm_base_info(num)%max_zones,&
!      ggm_steps(num)%n_outsteps))
!
!    ggm_base_info(num)%zone_pointer_names = "Null"
!    ggm_base_info(num)%amr_solution_names = "Null"
!
!! Get to the BaseIterativeData_t to write a bunch of stuff
!
!    CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',1,'end')
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!    call cgp_array_write_f('NumberOfZones',INTEGER,dummy_bound,&
!      int(ggm_steps(num)%n_outsteps,cgsize_t),index_array,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!    call cgp_array_write_data_f(index_array,INT(output_min_bound,cgsize_t),&
!      INT(output_max_bound,cgsize_t),ggm_base_info(num)%total_zones_at_ts(output_min_bound:&
!      output_max_bound),ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!! Assign all the information to zone_pointer_names and amr_solution_names
!!
!!
!! For every timestep
!    do ts = 1,ggm_steps(num)%n_outsteps
!! For every zone
!      do z = 1,ggm_base_info(num)%max_zones
!
!        marked = .false.
!! Loop through all zones at a given timestep
!        do l = 1,ggm_base_info(num)%total_zones_at_ts(ts)
!! Loop through all the levels to satisfy grid_info
!          do lvl = ggm_min_lvl(num),max_lvl
!! Check if any zones at the timestep match ones recorded during the run, if so
!!   put the zone name into zone_pointer_names.  Also put the solution name into
!!   amr_solution_names which goes into FlowSolutionPointers_t at the end
!            if (ggm_save_info(lvl,num)%multigrid_save(ts)%zone_number(l) == z) then
!
!              ggm_base_info(num)%zone_pointer_names(z,ts) = ggm_steps(num)%sol_names(ts)
!
!              write(char_save_number,format_nm_lvl_and_zone) lvl,'_',z
!              zonename = 'zone_lvl_'//TRIM(char_save_number)
!
!              ggm_base_info(num)%amr_solution_names(z,ts) = zonename
!
!              marked = .true.
!            end if
!
!            if (marked) exit
!          end do !levels
!
!        end do !zones at a given timestep
!      end do !zones
!    end do !timesteps
!!
!!   Write the zone with the size information
!!   It is unstructured because of the nature of LB grids. The cells are
!!   perfect cubes, but can be ignored because of geometry or other reasons
!!
!    CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',1,'end')
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!    call cgp_array_write_f('ZonePointers',CHARACTER,int(3,cgsize_t),&
!      (/idata(1),int(ggm_base_info(num)%max_zones,cgsize_t),&
!      INT(ggm_steps(num)%n_outsteps,cgsize_t)/),index_array,ier)
!
!    call cgp_array_write_data_f(index_array,(/min_dummy,int(output_min_bound,cgsize_t),min_dummy /),&
!      (/idata(1),int(output_max_bound,cgsize_t),int(ggm_steps(num)%n_outsteps,cgsize_t)/),&
!      ggm_base_info(num)%amr_solution_names(output_min_bound:output_max_bound,1:&
!      ggm_steps(num)%n_outsteps),ier)
!!
!! Divide the arrays up roughly so the zone_names and solution_names
!!   can be written to the CGNS file
!!
!    do z = 1,ggm_base_info(num)%max_zones
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
!      call cgp_array_write_data_f(index_array,(/dummy_bound,output_min_bound/), &
!          (/idata(1),output_max_bound/), &
!          ggm_base_info(num)%amr_solution_names(z,output_min_bound:output_max_bound),ier)
!    end do
!
!
!
!  end if
!
!
!  CALL cgp_close_f(index_file,ier)
!  IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!  WRITE(11,200) num
!  WRITE(11,*) ''
!
!END DO

 200  FORMAT ('Finalizing GGM data file ',I4,' done.')

end subroutine

!  END IF

  !DEALLOCATE(timing_data)
  !DEALLOCATE(all_sol_names)

!        call cgp_array_write_data_f(index_array,int((/1,1/),cgsize_t),idata,&
!          all_sol_names,ier)
!        call cg_error_print_f
!        write(*,*) self,ier
!        IF (ier /= CG_OK) CALL cg_error_exit_f
