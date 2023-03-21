subroutine namelist_mic_open(mic_file)
!
! Reads in microphone data from a file because other ways were too
! annoying to deal with
!
! Called by: namelist_global
! Calls: error_out
!
USE nml_output
USE precise
use amrex_base_module
IMPLICIT NONE
INTEGER :: istat,i,type_mic
CHARACTER(len=32) :: mic_file
! Allocate the memory for the arrays to be stored in the module
! mic_locations is (x/y/z,number of the mic)
! mic_times is (start recording time, stop recording time)
! mic_type tells it what to record
!
! mic_frequency is how often to write data to the file,
! if the number is 10, it will write every 10 timesteps
allocate(mic_locations(3,num_mics))
allocate(mic_times(2,num_mics))
allocate(mic_type(6,num_mics))
allocate(mic_frequency(num_mics))
allocate(mic_lvl(num_mics))
! Open a different number than everything else
OPEN(UNIT=14,FILE=mic_file,STATUS='OLD',ACTION='READ',FORM='FORMATTED',&
  IOSTAT=istat)
IF (istat .NE. 0) THEN
  WRITE(*,*) 'Error opening microphone file, exiting'
  WRITE(11,*) 'Error opening microphone file, exiting'
  WRITE(*,*) istat
  CALL error_out
END IF
!
! Loop through all the mics declared
!
mic_type = .FALSE.
WRITE(11,*) 'Microphone data'
DO i = 1,num_mics
!
! THIS IS NOT A NAMELIST READ
! This is just data from a file, possibly ported from Excel or Matlab
!
! mic_type will simply be a number, and it will be based on binary
! various components will be turned on and off based on bit placement
! i.e. 101101
! pressure,temperature,density,u,v,w, total up to 63
!
!
!
  READ(14,*,IOSTAT=istat) mic_locations(1,i),mic_locations(2,i),&
    mic_locations(3,i),type_mic,mic_times(1,i),mic_times(2,i)&
    ,mic_frequency(i),mic_lvl(i)
  WRITE(11,106) mic_locations(1,i),mic_locations(2,i)&
    ,mic_locations(3,i),&
    type_mic,mic_times(1,i),&
    mic_times(2,i),mic_frequency(i),mic_lvl(i)
  IF (istat .NE. 0) THEN
    WRITE(*,*) 'Error reading microphone file, exiting'
    WRITE(11,*) 'Error reading microphone file, exiting'
    WRITE(*,*) istat
    CALL error_out
  END IF

  IF (type_mic >=32) THEN
    mic_type(1,i) = .TRUE.
    type_mic = type_mic - 32
  END IF
  IF (type_mic >=16) THEN
    mic_type(2,i) = .TRUE.
    type_mic = type_mic - 16
  END IF
  IF (type_mic >=8) THEN
    mic_type(3,i) = .TRUE.
    type_mic = type_mic - 8
  END IF
  IF (type_mic >=4) THEN
    mic_type(4,i) = .TRUE.
    type_mic = type_mic - 4
  END IF
  IF (type_mic >=2) THEN
    mic_type(5,i) = .TRUE.
    type_mic = type_mic - 2
  END IF
  IF (type_mic >=1) THEN
    mic_type(6,i) = .TRUE.
    type_mic = type_mic - 1
  END IF

END DO

CLOSE(3)
WRITE(11,*) ''

! 100  FORMAT (F9.7)
 106  FORMAT ('x = ',F11.5,' y= ',F11.5,' z= ',F11.5,' microphone type = ',&
          I4,' start time = ',F11.5,' stop time = ',&
          F11.5,' microphone frequency = ',I8,'microphone max level',I2)
end subroutine
