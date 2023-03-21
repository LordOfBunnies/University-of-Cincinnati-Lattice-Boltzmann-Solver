SUBROUTINE null_bc
!
! This shouldn't be needed, but
!
! Called by: inlet_bc
! Calls:
!
USE precise
USE constants
USE freestream_values
IMPLICIT NONE

WRITE(*,*) 'Using null BC, something has gone wrong.'


END SUBROUTINE
