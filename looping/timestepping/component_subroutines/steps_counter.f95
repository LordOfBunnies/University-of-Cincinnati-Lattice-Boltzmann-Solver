SUBROUTINE steps_counter(level_step,fine)
!
! This increments the local timesteps and local time for every level.
! If a level is not doing an interation this time, nothing changes. If
! a level is doing an iteration it's time and step are increased. For levels
! above the current finest level, their steps and time are increased proportionately
!
! Called by: loop
! Calls:
!
use precise
use amr_info_holder, only: timestep,time,time_prev,dt
use amr_processes, only: amr_max_lvl
IMPLICIT NONE
INTEGER :: fine,lvl
LOGICAL :: level_step(0:amr_max_lvl)

DO lvl = 0,amr_max_lvl
  IF (level_step(lvl)) THEN
    timestep(lvl) = timestep(lvl)+1
    time_prev(lvl) = time(lvl)
    time(lvl) = time(lvl)+dt(lvl)


  ELSE IF (lvl > fine) THEN
    timestep(lvl) = timestep(lvl)+2**(lvl-fine)
    time_prev(lvl) = time(lvl)
    time(lvl) = time(lvl)+dt(lvl)*2**(lvl-fine)
  END IF

END DO

!if (all(level_step)) then
!  write(*,*) 'current time',time
!end if

END SUBROUTINE
