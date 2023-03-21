SUBROUTINE bounceback(fout_1,fout_2,gout_1,gout_2)
!
! Swaps two values, used at walls for non-equilibrium bounce back
!
! Called by: loop
! Calls:
!
USE precise
IMPLICIT NONE
REAL(KIND = dp),INTENT(OUT) :: fout_1,fout_2,gout_1,gout_2
REAL(KIND = dp) :: dumdum

dumdum = fout_1
fout_1 = fout_2
fout_2 = dumdum

dumdum = gout_1
gout_1 = gout_2
gout_2 = dumdum

END SUBROUTINE
