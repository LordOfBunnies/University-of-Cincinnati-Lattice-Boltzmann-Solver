      SUBROUTINE namelist_geometry
!
! Reads in the filenames of geometry files from the namelist
!
!
! Called by: namelist_global
! Calls: error_out
!
      USE nml_geometry
      IMPLICIT NONE
      INTEGER :: i,istat,geom_units
      CHARACTER(len=80) :: filename
!
! Declare the NAMELIST for the geometry names
      NAMELIST /GEOMETRY/ filename,geom_units
! Allocate the array of names
      ALLOCATE(geom_filenames(num_geom_files))
      ALLOCATE(units(num_geom_files))
      ALLOCATE(num_tris(num_geom_files))
! Initialize the names
      filename = 'geom_file.stl'
      geom_units = 1
!
! Units meaning
! 1 = meters, basic unit for UCLBS
! 2 = millimeters
! 3 = centimeters
! 4 = inches
! 5 = feet
!
! Run through all the geometry files to load
      DO i = 1,num_geom_files
        READ(2,GEOMETRY,IOSTAT=istat)
        IF (istat .NE. 0) THEN
          WRITE(*,*) 'Error reading namelist at geometry ', i
          WRITE(11,404) i
          WRITE(*,*) istat
          CALL error_out
        END IF
! Assign the names to the array
        geom_filenames(i) = filename
        units(i) = geom_units
        WRITE(11,252) i, filename
        IF (units(i) == 1) THEN
          WRITE(11,*) 'Geometry is in meters, no conversion required.'
        ELSE IF (units(i) == 2) THEN
          WRITE(11,*) 'Geometry is in millimeters, conversion to meters',&
            ' will be automatic.'
        ELSE IF (units(i) == 3) THEN
          WRITE(11,*) 'Geometry is in centimeters, conversion to meters',&
            ' will be automatic.'
        ELSE IF (units(i) == 4) THEN
          WRITE(11,*) 'Geometry is in inches, conversion to meters',&
            ' will be automatic.'
        ELSE IF (units(i) == 5) THEN
          WRITE(11,*) 'Geometry is in feet, conversion to meters',&
            ' will be automatic.'
        ELSE
          WRITE(11,*) 'Invalid units selected.'
        END IF
      END DO
      WRITE(11,*) ''

 252  FORMAT ('Geometry file number ',I8,' is ', A80)
 404  FORMAT ('Error reading namelist at geometry ', I4)
      END SUBROUTINE
