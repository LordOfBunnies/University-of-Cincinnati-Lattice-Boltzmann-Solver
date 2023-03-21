SUBROUTINE mic_initialize
!
! This initializes the microphone directory and creates all the needed
! microphone files in order to be able to access them later. It will
! overwrite any mic files already extant so be careful
!
! Called by: mic_setup
! Calls: error_out
!
USE nml_output
USE output_data
USE precise
use mpi
use amr_processes, only: self,commune,nprocs
IMPLICIT NONE
INTEGER :: a,ier,mic_unit
!      CHARACTER(len=64) :: sys_func,mic_directory,
CHARACTER(len=64) :: mic_file,file_path
CHARACTER(len=96) :: header,header_adder,junk

!      sys_func = 'mkdir -p'
!      mic_directory = ' microphones'

junk = '(I5.5)'

!      WRITE(*,*) TRIM(sys_func)//TRIM(mic_directory)

!      CALL SYSTEM(TRIM(sys_func)//TRIM(mic_directory),ier)
!      IF (ier /=0) THEN
!        WRITE(*,*) 'Something went wrong creating microphone directory.'
!        WRITE(11,*) 'Something went wrong creating microphone directory.'

!        CALL error_out
!      END IF

!      WRITE(11,111) sys_func//TRIM(mic_directory)
!      mic_directory = 'microphones/'

WRITE(11,*) ''
WRITE(11,*) 'Beginning creation of microphone files.'
!WRITE(*,*) 'Beginning creation of microphone files.'
mic_initialize_loop: DO a = 1,num_mics

  IF (.NOT. mic_valid(a)) THEN
    CYCLE mic_initialize_loop
  END IF

  WRITE(mic_file,junk) a
!  write(*,*) 'processor doing the writing',self,mic_owner(a)
!
!        mic_unit = MOD(a,100)+100
  mic_unit = 10
!        file_path = TRIM(mic_directory)//'microphone'//TRIM(mic_file)//'.txt'
  file_path = 'microphone'//TRIM(mic_file)//'.txt'
!
!  WRITE(*,109) file_path
!        WRITE(11,109) file_path
!
  OPEN(UNIT = mic_unit,FILE = file_path,FORM = 'FORMATTED',&
    STATUS ='REPLACE',ACTION = 'WRITE',IOSTAT=ier)
  IF (ier /=0) THEN
    WRITE(*,*) 'Something went wrong creating microphone file ',&
      mic_unit, ' for microphone ', a
  CALL error_out
  END IF

!  if (self == mic_owner(a)) then
    !cycle mic_initialize_loop


!
!        mic_file = 'Microphone_'//CHAR(a)//'.txt'
!        WRITE(11,108) mic_file
!

!        WRITE(11,*) mic_unit,ier
!        WRITE(*,*) mic_unit,ier


    WRITE(11,200) a
!        mic_type_num = mic
        header = ''
! 8 spaces
    IF (mic_type(1,a)) THEN
      header_adder = 'P,          ;'
      header = TRIM(header)//TRIM(header_adder)
    END IF
    IF (mic_type(2,a)) THEN
      header_adder = 'Temp,       ;'
      header = TRIM(header)//TRIM(header_adder)
    END IF
    IF (mic_type(3,a)) THEN
      header_adder = 'Density,    ;'
      header = TRIM(header)//TRIM(header_adder)
    END IF
    IF (mic_type(4,a)) THEN
      header_adder = 'U-vel,      ;'
      header = TRIM(header)//TRIM(header_adder)
    END IF
    IF (mic_type(5,a)) THEN
      header_adder = 'V-vel,      ;'
      header = TRIM(header)//TRIM(header_adder)
    END IF
    IF (mic_type(6,a)) THEN
      header_adder = 'W-vel,      ;'
      header = TRIM(header)//TRIM(header_adder)
    END IF

!      CALL SYSTEM('pwd',ier)
!      IF (ier /=0) THEN
!!        WRITE(*,*) 'Something went wrong  microphone directory.'
!       WRITE(11,*) 'Something went wrong  microphone directory.'
!       CALL error_out
!      END IF

    WRITE(mic_unit,*) TRIM(header)
    WRITE(11,202) a
!    WRITE(*,202) a
!  end if
!  call mpi_barrier(commune,ier)
  CLOSE(mic_unit)

END DO mic_initialize_loop
!write(*,*) 'boople',self
!call mpi_barrier(commune,ier)
!      WRITE(11,*) 'Finished creating microphone files.'
!WRITE(*,*) 'What is going on here?',self

 200  FORMAT ('Initializing microphone ',I8)
 201  FORMAT ('Microphone file number ',A8)
 202  FORMAT ('Header written for microphone ',I8)
 108  FORMAT ('Creating microphone file ',A64)
 109  FORMAT ('Putting file in path ',A32,' !')
! 300  FORMAT (A32)
! 302  FORMAT (A64)

END SUBROUTINE
