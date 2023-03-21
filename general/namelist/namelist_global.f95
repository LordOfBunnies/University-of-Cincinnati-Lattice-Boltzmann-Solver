subroutine namelist_global(characteristic_length)
!
!     Reads the general information from the input file. 
!     Initialized here, read elsewhere and brought back in.
!
!     Called by: uclbs_main
!     Calls: namelist_bounding_box,namelist_init,namelist_turbulent
!       namelist_output,namelist_amr_data
!
use nml_geometry
use nml_init
use nml_inlet_outlet
use nml_output
use nml_inputs
use constants
use freestream_values
use timing
use precise
use startup
use amr_info_holder, only: amr,multigrid
use ggm_stuff !XXXX
IMPLICIT NONE
INTEGER :: istat
REAL(KIND=dp) :: characteristic_length

CHARACTER(len=80) :: filename_in,filename_out,mic_file
! Start the CPU timer
CALL CPU_TIME(cpu_start_time)
!
! Declare the namelist for all the variables to come from the global list
! Some of these are not implemented anywhere yet and won't be for a long
! time.
!
!
! Basic filename for input, this won't change unless we get fancy
filename_in = 'uclbs.in'
filename_out = 'uclbs.out'
!
!
! This is to write the user input data to uclbs.out.  Essentially it will
! reiterate everything the user put in so that the user can catch errors they
! made during input if that is needed. It will also include important data down
! the line that the user may want to see. This will be a catchall location for
! problems down the line.
!
OPEN(FILE=filename_out,UNIT = 11,STATUS='REPLACE',FORM = 'FORMATTED',&
  ACTION='WRITE',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Problems opening output file, closing.'
  WRITE(*,*) istat
  CALL error_out
END IF
WRITE(11,*) 'UCLBS output file'
WRITE(11,*) 'Beginning with basic user inputs.'
WRITE(11,*) 'VR starts at 10 and will change based on the data.'
!
! Open the input file
!
OPEN(UNIT=2,FILE=filename_in,STATUS='OLD',ACTION='READ',&
  FORM='FORMATTED',IOSTAT=istat)
IF (istat .NE. 0) THEN
  WRITE(*,*) 'Error opening namelist file, exiting'
  WRITE(11,*) 'Error opening namelist file, exiting'
  WRITE(*,*) istat
  CALL error_out
END IF
! Read the global data from the namelist file
!      WRITE(*,*) 'Things happening'
!
CALL namelist_global_read(characteristic_length,mic_file)
! Start up the namelist read calls
! 
! Call geometry namelist read subroutine
!
IF (num_geom_files >= 1) THEN
  CALL namelist_geometry
END IF
!
! Call the namelist read subroutine to read inlet and outlets
!
CALL namelist_inlet_outlet
!
! Call the namelist reader for info on the boudning box
!
!IF (.NOT. cgns_start .AND. .NOT. restart_start) THEN
!  CALL namelist_bounding_box
!END IF
!
! Call the namelist reader for info on initialization
!
IF (num_init_regions >=1) THEN
  CALL namelist_init
END IF
!
! Call the cgns reader if this is being used 
!
IF (num_cgns_files >= 1) THEN
  CALL namelist_cgns_files
!        CALL uclbs_cgns_read
END IF
!
!
!

!
!
! Call the namelist reader for output info
IF (num_saves >= 1) THEN
  CALL namelist_output
END IF
!
IF (num_mics >= 1) THEN
  CALL namelist_mic_open(mic_file)
END IF
!
if (amr .or. multigrid) then
  CALL namelist_amr_data
end if

IF (use_ggm) THEN !XXXX
  CALL namelist_ggm_data
END IF
!
CALL namelist_error_checker(characteristic_length)
!
! 100  FORMAT (F9.7)
! 101  FORMAT (2F9.7)
! 102  FORMAT (3F9.7)
! 103  FORMAT (4F9.7)
! 104  FORMAT (5F9.7)
! 105  FORMAT (6F9.7)
! 106  FORMAT (3F9.7,I8,2F9.7,I8)
! 150  FORMAT (F9.7,I8)
! 151  FORMAT (I8,F9.7)
!
! The 200s are primarily integer output formats
!
! 200  FORMAT (I8)
! 201  FORMAT (2I8)
! 202  FORMAT (3I8)
! 251  FORMAT (I8, A32)
!
! The 300s are primarily character output formats
!
! 300  FORMAT (A32)
!
! The 400s are primarily for logical output formats
!
END SUBROUTINE
