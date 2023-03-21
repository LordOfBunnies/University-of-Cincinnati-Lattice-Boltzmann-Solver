MODULE precise
!
! This defines double precision for all REAL values for the entire program.
! This module is USEd by every other module and most subroutines.
!
! The epsilons are for use with closeness to machine zero. If something is
! between those two, it is effectively zero.
!
!
IMPLICIT NONE
SAVE
INTEGER,PARAMETER :: dp = (SELECTED_REAL_KIND(15,307))
REAL(KIND=dp) :: plus_epsilon = 1.24D-16,minus_epsilon = &
  -1.24D-16, grid_move_delta = 1.0D-7,sp_plus_epsilon = 1.25D-7,&
  sp_minus_epsilon = -1.25D-7,semi_epsilon = 1.0D-8
REAL(KIND=dp),PARAMETER :: very_small_number = 1.0D-14
REAL(KIND=dp) :: very_large_number = 1.0D16


END MODULE
