program uclbs_main
!     
! Main part of the program for the Lattice Boltzmann Solver.  UCLBS
! is meant to be a fully functional LBS that can be used beyond the scope
! of just my PhD dissertation.  To that end it's more complicated than I
! would normally code for just myself. This is an ongoing development and
! there are many, MANY parts I'd like to add eventually.
!
! Written by: Sean Duncan
! Date: 9/1/18
! Current version: alpha v0.02
!
! All rights reserved by Sean Duncan
!
use precise
use mpi
use amrex_base_module
!      USE nml_inputs
use startup
implicit none
real(KIND=dp) :: characteristic_length
!
!
! First step is to call the input reader.  Most of the important data is
! stored in modules that are USEd as needed
!
!
call namelist_global(characteristic_length)
!
call calcs_before_init
! This imports the geometry and sets up the arrays to be used later in
! dealing with the geometry
!
if (.not. restart_start) then
  if (cgns_start) then
!
! This is simply to calculate the basic flow quantities not fed to the program
! in the input file.
!

! This starts Amrex and MPI
    call external_code_starter()
    call uclbs_cgns_read

    call background_calcs(characteristic_length)
  else
    call geom_import
! This starts Amrex and MPI


    call external_code_starter()
    call background_calcs(characteristic_length)

  end if
else
!
! This is simply to calculate the basic flow quantities not fed to the program
! in the input file.
!
  call background_calcs(characteristic_length)
! This starts Amrex and MPI
  call external_code_starter()
  !call restart_read() !XXXX reimplement eventually once restart files with Amrex are figured out
end if

!if (.not. cgns_start) then
!  call external_code_starter()
!end if

!
! This creates the grid and determines whether the nodes are fluid or solid
!
call grid_gen(characteristic_length)
!
! This changes the grid from an x,y,z grid to a 2 dimensional list of links
! that are useful when doing LB work.  This will be especially important
! later when mesh refinement is used since basic x,y,z grids fly out the
! window.
!
!write(*,*) 'plow'
call links
!write(*,*) 'mrow'
!
! This initializes the data as specified in the namelist input
!
!
!      WRITE(*,*) 'Boop!'
! XXXX call temporarily suppressed in order to test
call init_fi
!write(*,*) 'cow'
!
!
!     Call the main looping subroutine
!      WRITE(*,*) 'Boop!'
! XXXX call temporarily suppressed in order to test
!else
!  call restart_read
!end if
! XXXX call temporarily suppressed in order to test
call loop
!
! Call the output finalization subroutine to close out all relevant
! files
!
!call output_finalize
!
!
call external_code_finalizer

end program uclbs_main
