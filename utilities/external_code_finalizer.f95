subroutine external_code_finalizer
!
! Frees all the memory and terminates Amrex so there are no errors
! at the end of the program
!
! Called by: output_finalize
! Calls:
! External calls: amrex_multifab_destroy,amrex_imultifab_destroy,
!   amrex_amrcore_finalize,amrex_finalize
!
use amrex_amr_module
use amrex_base_module
use mpi
use amr_info_holder
use constants
implicit none
integer :: level,ier

!
! Destroy all the mutlifabs at every level
!
!do level = 0,amrex_max_level
!  call amrex_multifab_destroy(mfrho(level))
!  call amrex_multifab_destroy(mfu_vel(level))
!  call amrex_multifab_destroy(mfv_vel(level))
!  call amrex_multifab_destroy(mftemp(level))
!  call amrex_multifab_destroy(mfpress(level))
!
!  call amrex_multifab_destroy(mfchi(level))
!  call amrex_multifab_destroy(mfzeta_x(level))
!  call amrex_multifab_destroy(mfzeta_y(level))
!  call amrex_multifab_destroy(mfpi_xx(level))
!  call amrex_multifab_destroy(mfpi_xy(level))
!  call amrex_multifab_destroy(mfpi_yy(level))
!  call amrex_multifab_destroy(mflambda_x(level))
!  call amrex_multifab_destroy(mflambda_y(level))
!
!  if (dimensions == 3) then
!    call amrex_multifab_destroy(mfzeta_z(level) )
!    call amrex_multifab_destroy(mfpi_xz(level) )
!    call amrex_multifab_destroy(mfpi_yz(level) )
!    call amrex_multifab_destroy(mfpi_zz(level) )
!    call amrex_multifab_destroy(mflambda_z(level) )
!    call amrex_multifab_destroy(mfw_vel(level))
!  end if
!
!  call amrex_multifab_destroy(mfalpha(level)  )
!  call amrex_imultifab_destroy(mfstate(level) )
!
!  call amrex_multifab_destroy(mffi(level))
!  call amrex_multifab_destroy(mffout(level))
!  call amrex_multifab_destroy(mfgi(level))
!  call amrex_multifab_destroy(mfgout(level))
!end do
!
! Terminate Amrex and MPI along with it.
!
call amrex_amrcore_finalize()
call amrex_finalize()
!call mpi_finalize(ier)
end subroutine
