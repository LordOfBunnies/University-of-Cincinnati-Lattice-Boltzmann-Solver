SUBROUTINE valid_levels(level_step,fine,coarsest_lvl_run)
!
! This calculates which levels will be timestepping during the current
! timestep.
!
!
! Called by: loop
! Calls:
!
USE precise
use constants
use timing
use freestream_values
use amr_info_holder, only: timestep,elbm_epsilon,init_epsilon,startup_iterations,startup,time
use amr_processes, only: amr_max_lvl,self
!use amr_processes, only: amr_max_lvl
IMPLICIT NONE
INTEGER :: fine,lvl,coarsest_lvl_run
LOGICAL :: level_step(0:amr_max_lvl)

level_step = .false.

DO lvl = fine,0,-1

!  if (timestep(amr_max_lvl) < startup_timesteps) then
!    !elbm_epsilon(lvl) = (timestep(fine)+1)*init_epsilon(lvl)
!    init_epsilon(lvl) = (timestep(amr_max_lvl)+1)/real(startup_timesteps,kind=dp)*elbm_epsilon(lvl)
!!   init_epsilon(lvl) = elbm_epsilon(lvl)
!    write(*,*) 'Init elbm epsilon',elbm_epsilon(lvl),init_epsilon(lvl),startup_timesteps,lvl,self
!    startup = .true.
!  else
!    startup = .false.
!    init_epsilon(lvl) = elbm_epsilon(lvl)
!    write(*,*) 'Elbm epsilon',elbm_epsilon(lvl),init_epsilon(lvl),startup_timesteps,lvl,self
!  end if
!
! The fine level always runs
!
  IF (lvl == fine) THEN
    level_step(lvl) = .TRUE.
    coarsest_lvl_run = fine
!    write(*,*) 'valid level',lvl,level_step
  ELSE
!
! For all the coarser levels
!
    IF (MODULO(timestep(lvl + 1),2)==1 .and. level_step(lvl+1)) THEN
      level_step(lvl) = .TRUE.
      coarsest_lvl_run = lvl
!      write(*,*) 'valid level',lvl,level_step
    ELSE
      level_step(lvl) = .FALSE.
      !write(*,*) 'not running',lvl
    END IF

  END IF !fine level logical


END DO

if (self == 0) then
  write(*,*) 'time check',fine,timestep,time,level_step
end if


END SUBROUTINE
