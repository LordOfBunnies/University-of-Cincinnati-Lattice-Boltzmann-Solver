      SUBROUTINE ggm_output_finalize
!
! Finalize the cgns outputs for the ggm
!
! Called by: output_finalize
! Calls:
! Outside Calls: cg_open_f,cg_biter_write_f,cg_simulation_type_write_f,cg_close_f
!
USE precise
USE cgns
USE loop_data
USE constants
USE ggm_stuff
USE timing
IMPLICIT NONE
INTEGER :: ier,a,ts_save
INTEGER*8 :: idata(2),funk,dummy1,dummy2
INTEGER :: index_base,index_file,index_coord,index_zone,index_section,&
  index_flow,index_field
INTEGER*8 :: i
!      REAL(KIND=dp) :: ts_start,ts_end,ts_save,
CHARACTER(len=32),ALLOCATABLE :: solut_name(:)
CHARACTER(len=32) :: filename,char_save_number,solut_save_number,&
  format_junk,format_junk_ts,basename

REAL(KIND=dp),ALLOCATABLE :: timing_data(:)

index_base = 1
!      index_file = 1
index_coord = 1
index_zone = 1
index_section = 1
index_flow = 1
index_field = 1
format_junk_ts = '(I9.9)'
format_junk = '(I3.3)'

dummy1 = 1
dummy2 = 2

DO a = 1,num_ggm_vars


  WRITE(char_save_number,format_junk) a

!        WRITE(*,*) SIZE(x_vertices),SIZE(y_vertices),SIZE(z_vertices)

  filename = 'ggm_output_'//TRIM(char_save_number)//'.cgns'
  basename = 'outputbase_'//TRIM(char_save_number)
!        zonename = 'outputzone_'//TRIM(char_save_number)

  ts_save = timesteps/ggm_output_freq

!        WRITE(*,*) ts_start,ts_end,winkin,blinkin,end_times,enders_timer
!        WRITE(*,*) ts_save,ts_last_save,ts_first_save,save_frequency(a)
!        WRITE(*,*) total_time,ts_first_save*dt

  idata(1) = 32
  idata(2) = ts_save
!
! Allocate timing data and solution names
!
  ALLOCATE(timing_data(ts_save))
  ALLOCATE(solut_name(ts_save))

!        nod = 1
  DO i = 1,ts_save

    timing_data(i) = DBLE(i)*dt*DBLE(ggm_output_freq)
!          timing_data(i) = REAL(i)*dt*REAL(save_frequency(a))
    WRITE(solut_save_number,format_junk_ts) i*ggm_output_freq
    solut_name(i) = 'Timestep_'//TRIM(solut_save_number)

  END DO
  WRITE(*,*) ts_save,idata

!
!
!
  IF (total_time < DBLE(ggm_output_freq)*dt ) THEN

    CALL cg_open_f(filename,CG_MODE_MODIFY,index_file,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'GGM output file reopened for finalization.'
!          WRITE(*,*) 'Output file reopened for finalization.'
!
!
!
!          WRITE(*,*) index_file,index_base,'TimeIterValues',&
!            ts_save,ier
    CALL cg_biter_write_f(index_file,index_base,'TimeIterValues',&
      ts_save,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing time iterations values to finalize GGM file.'
!          WRITE(*,*) 'Writing time iterations values to finalize file.'
!
!
!
!          WRITE(*,*) index_file,index_base,TimeAccurate,&
!            ier
    CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
      ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing solution type to finalize file.'
    WRITE(*,*) 'Writing solution type to finalize file.'

    CALL cg_close_f(index_file,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
!
!
  ELSE

!
! Open the file to be appended to
!
!         WRITE(*,*) ''
!          WRITE(*,*) 'Begin finalization of output save ', a
    WRITE(11,*) 'Begin finalization of ggm file ', a
    CALL cg_open_f(filename,CG_MODE_MODIFY,index_file,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Output file reopened for finalization.'
!          WRITE(*,*) 'Output file reopened for finalization.'
!
! Call the base iterative write to tell it how many timesteps there are
!
    CALL cg_biter_write_f(index_file,index_base,'TimeIterValues',&
      ts_save,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing time iterations values to finalize ggm file.'
!          WRITE(*,*) 'Writing time iterations values to finalize file.'
!
! Make it go to the end of the file for more writing
!

    CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',&
      dummy1,'end')
    IF (ier /= CG_OK) CALL cg_error_exit_f
!
! Write the actual times of the timesteps into the CGNS file
!
    funk = ts_save
    CALL cg_array_write_f('TimeValues',RealDouble,dummy1,funk,&
      timing_data,ier)

    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing time values to finalize file.'
!
! Tell it zone iterative data is present
!
    CALL cg_ziter_write_f(index_file,index_base,index_zone,&
      'ZoneIterativeData',ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing zone header to finalize file.'
!
! Go to the zone data
!
    CALL cg_goto_f(index_file,index_base,ier,'Zone_t',index_zone,&
      'ZoneIterativeData_t',1,'end')
    IF (ier /= CG_OK) CALL cg_error_exit_f
!
! Tell it the solution names and their sizes
!
    CALL cg_array_write_f('FlowSolutionPoints',Character,dummy2,idata,&
      solut_name,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing solution names to finalize file.'
!
! Tell it the solution is time accurate (no other way to do LBM)
!
    CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
      ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f
    WRITE(11,*) 'Writing solution type to finalize file.'
!
! Close out the file
!
    CALL cg_close_f(index_file,ier)
    IF (ier /= CG_OK) CALL cg_error_exit_f

    WRITE(11,200) a
    WRITE(11,*) ''


  END IF

  DEALLOCATE(timing_data)
  DEALLOCATE(solut_name)

END DO

 200  FORMAT ('Output of ggm file ',I4,' done.')
 201  FORMAT ('Beginning ggm finalization of output file for save ',I4)
! 501  FORMAT ('Total elapsed time running UCLBS = ',F11.5)

END SUBROUTINE
