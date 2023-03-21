SUBROUTINE error_out
!
! Called if there was a problem loading something and a need to shut down
!
!
! Called by: damn near everything
! Calls: MPI stuff when I get around to putting that in
!
IMPLICIT NONE

WRITE(11,*) 'Error occured, exiting.'
WRITE(*,*) 'Error in coding, shutting down gracefully.'

STOP
END
