subroutine external_code_starter
!
!
!
!
! Called by: uclbs_main
! Calls: amr_transfer_data_to_amrex
! External calls: mpi_init,amrex_init
!
use mpi
!use cgns
use amrex_base_module
use amrex_amr_module
use amr_processes

IMPLICIT NONE
integer :: ierr,mpi_ierr

!call mpi_init( mpi_ierr )
!if (mpi_ierr /= MPI_SUCCESS) then
!  WRITE(*,*) ' Unexpected return from MPI_INIT', ierr
!  call error_out
!endif
!write(*,*) 'huh?'
call amrex_init()
!write(*,*) 'bonk?'
call amrex_amrcore_init()
!write(*,*) 'wiggles'
call init_amrex()
call external_data_read

!WRITE(*,*) 'Initializing from scratch'
CALL amrex_init_from_scratch(0.0D0)
!WRITE(*,*) 'Initialized'

!call amr_transfer_data_to_amrex

!call init_things

end subroutine
