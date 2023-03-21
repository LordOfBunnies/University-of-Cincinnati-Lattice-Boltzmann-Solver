subroutine mic_on_regrid(fine)
!
! Adjust the mic lvl and associated nodes as needed when a regrid is performed
!
!
!
!
!
use precise
use amrex_base_module
use output_data
use nml_output
use constants

implicit none
integer :: fine,i

do i = 1,num_mics
  if (fine >= mic_lvl(i)) then
    cycle
  else

  end if


end do

end subroutine
