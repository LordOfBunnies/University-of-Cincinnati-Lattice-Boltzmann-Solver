MODULE startup
!
! This defines various startup possibilities including
! a restart, a cold start, a CGNS start, or a mixed CGNS
! start with geometry import
!
!
!
IMPLICIT NONE
SAVE
LOGICAL :: restart_start,cgns_start,mixed_start
INTEGER :: t_start


END MODULE startup
