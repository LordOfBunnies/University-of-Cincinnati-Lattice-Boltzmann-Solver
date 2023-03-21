SUBROUTINE lb_output(fine)
!
! This is called every timestep, though it should just shunt back out
! if it doesn't meet the requirements of
!
! Called by: loop
! Calls: save_write,restart_write,mic_write
!
USE linkwise
USE nml_output
USE precise
use output_data
use ggm_stuff
use zone_writing
use amr_info_holder, only: dt,time,timestep,amr
!USE ggm_stuff !XXXX
IMPLICIT NONE
INTEGER :: a,fine,max_lvl
REAL(KIND = dp) :: time_passed
!
! Deal with writing output
!
IF (num_saves > 0) THEN

!        CALL assign_output_dirs
!
! Need to redo the setup since it was based on the timestep and is now based on time.
!
! Basic premise is call the appropriate save routine when the time at the finest level is
! between t and t-dt, thus when it steps on it isn't called.
!
  DO a = 1, num_saves

    if (save_shapes(2,a) > fine) then
      max_lvl = fine
    else
      max_lvl = save_shapes(2,a)
    end if

    if (abs(time(max_lvl) - save_steps(a)%saved_times(save_steps(a)%output_counter)) &
        <= dt(max_lvl) .and. .not. save_steps(a)%write_done(save_steps(a)%output_counter) ) then
!        CALL save_write(a,fine)
         !call hdf5_output_test(fine)
         call plot3d_output(fine)
         save_steps(a)%output_counter = save_steps(a)%output_counter + 1
!
! Remove the zone information
!
!        if (amr) then
!          call zone_dealloc_safety(fine,a,0,fine)
!        end if
      END IF
!    END IF
  END DO

!  write(*,*) 'spank you very much'
END IF
!
! Deal with writing things out to mic files
!
IF (num_mics > 0) THEN
  mic_write_loop: DO a = 1,num_mics

    IF (.NOT. mic_valid(a)) THEN
      CYCLE mic_write_loop
    END IF
!    IF (time_passed > mic_times(1,a) .AND. time_passed < &
!          mic_times(2,a)) THEN
      IF (abs(timestep(fine)-mic_frequency(a)) <= dt(fine)) THEN

        CALL mic_write(a)

      END IF
!    END IF
  END DO mic_write_loop

END IF
!
!
!
!if (ggm_output) then
!!  if (num_ggm_vars > 0) then
!
!  DO a = 1, num_ggm_vars
!    if (ggm_max_lvl(a) > fine) then
!      max_lvl = fine
!    else
!      max_lvl = ggm_max_lvl(a)
!    end if
!!    write(*,*) 'ggm output checking',time(max_lvl),dt(max_lvl),a,fine,&
!!      ggm_steps(a)%saved_times(ggm_steps(a)%output_counter),ggm_steps(a)%output_counter,&
!!      ggm_steps(a)%write_done(ggm_steps(a)%output_counter)
!    if (abs(time(max_lvl) - ggm_steps(a)%saved_times(ggm_steps(a)%output_counter)) &
!      <= dt(max_lvl) .and. .not. ggm_steps(a)%write_done(ggm_steps(a)%output_counter) ) then
!
!!         CALL save_write(a,fine)
!!       call ggm_output_shortcut(a,fine)
!!
!! Remove the zone information
!!
!!      if (amr) then
!!        call zone_dealloc_safety(fine,a,ggm_min_lvl(a),ggm_max_lvl(a))
!!      end if
!    END IF
!!    END IF
!  END DO
!
!end if

!
! Deal with writing restart saves
!
DO a = 1,num_restart_saves
  IF (MOD(time(fine),restart_saves(a)) <= dt(fine) .AND. &
        MOD(time(fine),restart_saves(a)) > 0.) THEN
    CALL restart_write(a)
  END IF
END DO

END SUBROUTINE

!    IF (time_passed > save_times(1,a) .AND. time_passed < &
!          save_times(2,a)) THEN
!      IF (MOD(t,save_frequency(a)) == 0) THEN
!      if ( mod(time(save_shapes(2,a)),save_frequency(a)) <= dt(save_shapes(2,a)) .AND. &
!        abs(time(save_shapes(2,a))-save_frequency(a) > 0. ) then

!IF (use_ggm) THEN !XXXX
!
!!        CALL assign_output_dirs
!
!!        DO a = 1, num_ggm_vars
!    IF (MOD(t,ggm_output_freq) == 0) THEN
!
!      CALL ggm_calculator(t)
!
!    END IF
!!        END DO
!END IF

!      time_passed = dt*REAL(t)
