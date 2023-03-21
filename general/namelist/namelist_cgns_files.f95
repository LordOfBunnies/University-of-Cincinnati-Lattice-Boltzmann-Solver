SUBROUTINE namelist_cgns_files
!
! This subroutine gets all the CGNS file names out of the namelist file
!
! Calls: error_out
! Called by:
!
USE nml_inputs
IMPLICIT NONE
INTEGER :: i,istat,units
CHARACTER(len=64) :: cgns_file

NAMELIST /CGNS_FILES/ cgns_file,units

ALLOCATE(cgns_input_files(num_cgns_files))
ALLOCATE(cgns_units(num_cgns_files))

cgns_file = 'grid.cgns'
units = 1

DO i = 1,num_cgns_files

  READ(2,NML=CGNS_FILES,IOSTAT=istat)
  IF (istat .NE. 0) THEN
    WRITE(*,*) 'Error reading namelist file at CGNS file ',&
      i,', exiting'
    WRITE(11,*) 'Error reading namelist file at CGNS file ',&
      i,', exiting'
    CALL error_out
  END IF

  cgns_input_files(i) = cgns_file
  cgns_units(i) = units

END DO

END SUBROUTINE namelist_cgns_files
