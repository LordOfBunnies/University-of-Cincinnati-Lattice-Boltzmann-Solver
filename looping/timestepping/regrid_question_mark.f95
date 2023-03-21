subroutine regrid_question_mark(fine)
!
! This tells asks Amrex to do a regrid operation.  It should only work if the AMR
!   option is turned on during initial data read.
!
!
! Called by: loop
! Calls:
! External calls: amrex_regrid
!
use amrex_amr_module
use amrex_base_module
use precise
use linkwise
use timing
use ggm_stuff
use constants
use amr_processes
use amr_info_holder
use startup
use amr_grid_writing, only: amr_grid_zones
implicit none
integer :: lvl,i,min_coarse,dummy,old_fine
integer :: fine,new_fine,ier
logical :: regrid_lvls(0:amrex_max_level),valid_coarse

regrid_lvls = .false.
!
! Loop through the levels to see if the timestep matches any of the regrid_timesteps
!

!
! regrid_basis records the last regrid number, so we don't have to loop through all of them
!   every time.
!

! Bounce back if new grid refinement is going to happen
if (.not. multigrid) return

if (.not. amr .and. coarse_time < very_small_number) return

if (.not. use_ggm .and. coarse_time < very_small_number) return

if (coarse_time > very_small_number) then
  min_coarse = minval(coarse_box)
end if



!write(*,*) 'minimum coarse thingy ',min_coarse,coarse_time,coarse_regrid_done
dummy = 0
!
! Refining the grid to the appropriate box levels once the given coarse time is reached.
!
if (coarse_time > very_small_number) then
  if (time(fine) >= coarse_time - dt(fine)/2 .and. time(fine) <= coarse_time + dt(fine)/2 .and. &
      .not. coarse_regrid_done) then
!    write(*,*) 'fuzzy wuzzy'
    valid_coarse = .false.

    do i = 1,num_boi
!      write(*,*) 'epple'
      if (box_lvl(i) <= fixed_lvl .and. .not. valid_coarse) then

      else
        valid_coarse = .true.
      end if

      do lvl = min_coarse,box_lvl(i)-1
        write(*,*) 'calling regrid',lvl,box_lvl(i)
        call amrex_regrid(lvl,time(lvl))
        regrid_lvls(lvl) = .true.
      end do
    end do

    new_fine = amrex_get_finest_level()

  do lvl = 0,new_fine
    call update_grid_info(lvl,dummy)
  end do


!    write(*,*) 'new finery',new_fine
    if (new_fine > min_coarse) then
      do lvl = min_coarse+1,new_fine
!      write(*,*)
      if (lvl == 0) cycle
        call fluid_or_solid_from_coarse(lvl)
        if (cgns_start .or. mixed_start) then
          call state_directions_cgns_lvl(lvl)
        else if (shifted) then
          call shifted_state_directions_lvl(lvl)
        else
          call state_directions_lvl(lvl)
        end if

        call unfuck_coarse_interp(lvl)

      end do
!
    call mpi_barrier(commune,ier)
!        write(*,*) 'after fluid/solid'
! Eliminate unnecessary nodes
      do lvl = min_coarse,new_fine-1
        call nullify_lvl(lvl,dummy)
! Reinstate required coarse nodes at the edge of the fine grids
        call find_overlapped_coarse(lvl)
        call simple_state_directions(lvl)
      end do
!
    call mpi_barrier(commune,ier)
!
      do lvl = min_coarse,new_fine
        call ghost_lvl(lvl,dummy)
        call simple_state_directions(lvl)
      end do
!
!
!
      call mpi_barrier(commune,ier)
    end if

    coarse_regrid_done = .true.
    regrid_basis = regrid_basis+1
!    call initial_output

    return
  end if

end if
!
! Regular AMR operations based on timestep
!
if (timestep(fine) == regrid_timesteps(regrid_basis,fine) .and. amr .and. &
  time(fine) > coarse_time+ggm_dt/2.0D0) then

!  write(*,*) 'kinky boots!'
  if (fine < amr_max_lvl) then
    old_fine = fine
  else
    old_fine = amr_max_lvl - 1
  end if


!  do lvl = fixed_lvl,old_fine
!    write(*,*) 'sparks',lvl,fixed_lvl,old_fine
!    call amrex_multifab_build(mfgm(lvl),mffi(lvl)%ba,mffi(lvl)%dm,num_ggm_vars,nghosts,node_based_3d)
!    call amrex_multifab_build(mfggm(lvl),mffi(lvl)%ba,mffi(lvl)%dm,num_ggm_vars,0,node_based_3d)
!    write(*,*) 'fire',lvl
!  end do

  do lvl = fixed_lvl,old_fine
!    call amrex_regrid(fixed_lvl+1,time(fixed_lvl+1))
    call amrex_regrid(lvl,time(lvl))
!    regrid_lvls(lvl) = .true.
  end do

!  do lvl = fixed_lvl,old_fine
!    write(*,*) 'thunder',lvl
!    call amrex_multifab_destroy(mfgm(lvl))
!    call amrex_multifab_destroy(mfggm(lvl))
!    write(*,*) 'lightning',lvl
    call mpi_barrier(commune,ier)
!  end do
!
!  write(*,*) 'asscrack'

  regrid_basis = regrid_basis+1
  new_fine = amrex_get_finest_level()
!
  do lvl = 0,new_fine
    call update_grid_info(lvl,dummy)
    call mpi_barrier(commune,ier)
  end do
!
!write(*,*) 'honky mufuggas'
  do lvl = fixed_lvl,new_fine
    if (lvl == 0) cycle

    call fluid_or_solid_from_coarse(lvl)
    if (cgns_start .or. mixed_start) then
      call state_directions_cgns_lvl(lvl)
    else if (shifted) then
      call shifted_state_directions_lvl(lvl)
    else
      call state_directions_lvl(lvl)
    end if

    call unfuck_coarse_interp(lvl)
  end do

!write(*,*) 'schizo'
! Eliminate unnecessary nodes
  do lvl = fixed_lvl,new_fine-1
    call nullify_lvl(lvl,dummy)
! Reinstate required coarse nodes at the edge of the fine grids
    call find_overlapped_coarse(lvl)
    call simple_state_directions(lvl)
  end do
!write(*,*) 'soft, what light through yonder window breaks'
  do lvl = fixed_lvl+1,new_fine
    call ghost_lvl(lvl,dummy)
    call simple_state_directions(lvl)
  end do
!
!
!write(*,*) 'tis the moon, and Juliet is the Sun!'
!
!  do lvl = min_coarse,new_fine
!    call state_directions_lvl(lvl)
!  end do
end if
!end do
!write(*,*) 'oh happy dagger, this is thy sheath'
end subroutine
!    if (regrid_lvls(lvl)) then
!      call update_grid_info(lvl,dummy)
!
!      call fluid_or_solid_from_coarse(lvl)
!
!      call nullify_lvl(lvl-1,dummy)
!
!      call ghost_lvl(lvl,dummy)
!
!      call coarse_fine_lvl(lvl-1,dummy)
!!    end if
!    end if


!      if (regrid_lvls(lvl)) then
!

!      call amr_grid_zones(lvl,num)
        !if (lvl > min_coarse)
! Find what's a fluid and solid node

    !call update_grid_info(lvl,num)

!        write(*,*) 'after nullify lvl',lvl,self

!        write(*,*) 'after ghost lvl',lvl,self

!        call coarse_fine_lvl(lvl-1,dummy)
!        call simple_state_directions(lvl)
!        write(*,*) 'after c/f id lvl',lvl,self
!      end if
