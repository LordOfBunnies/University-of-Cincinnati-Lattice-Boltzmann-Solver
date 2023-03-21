subroutine pure_outflow_bc(fout_outl,fout_opp,gout_outl,gout_opp,loc_dir,outl_number)
!
! This is meant to allow pure outflow, or as close as possible without
!   saving silly amounts of data. It takes the fout value opposite the
!   outflow direction and sets that as the fout at the boundary.  The
!   streaming subroutine then makes this fout the opposite fi.
!
! Called by: outlet_bc
! Calls:
! External calls:
!
use precise
use constants
implicit none
integer,intent(in) :: loc_dir,outl_number
real(kind=dp),intent(in) :: fout_opp,gout_opp
real(kind=dp) ::fout_outl,gout_outl

fout_outl = fout_opp
gout_outl = gout_opp

end subroutine
