SUBROUTINE geom_import
!
! This is a import all the geometry files.  It needs to do
! an rapid open and close on them in order to read the number
! of triangles in each file.
!
! Called by: uclbs_main
! Calls:
!
USE nml_geometry
USE geom_data
USE precise
USE startup
IMPLICIT NONE
INTEGER(kind=2) :: bytes
INTEGER :: tri_num,junk2
INTEGER :: istat, i,j,tri_counter
!      INTEGER ::
REAL :: normal1,normal2,normal3
REAL :: vertex1x1,vertex1x2,vertex1x3
REAL :: vertex2x1,vertex2x2,vertex2x3
REAL :: vertex3x1,vertex3x2,vertex3x3
REAL(KIND = dp) :: adjusted_units
CHARACTER(len=80) :: header,junk1

CHARACTER(len=80) :: filename_out

filename_out = 'uclbs.out'
!
OPEN(FILE=filename_out,UNIT = 11,STATUS='OLD',FORM = 'FORMATTED',&
  ACTION='WRITE',POSITION='APPEND',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Problems opening output file, closing.'
  WRITE(*,*) istat
  CALL error_out
END IF
!
total_tris = 0
stl_tris = 0
IF (.NOT. cgns_start) THEN
  cgns_tris = 0
END IF
!
! Rapid open and close of the geometry files to read info at the top.
! Can't think of a better way to do this so may be revised later.
!
DO i = 1,num_geom_files
  OPEN(4,FILE=geom_filenames(i),STATUS='OLD',ACTION='READ', &
     FORM='UNFORMATTED',ACCESS='STREAM', IOSTAT=istat)
  IF (istat /= 0) THEN
    WRITE(*,*) 'Problem loading geometry file to read header info.'
!          WRITE(*,*) 'Problem loading geometry file to read header info.'
    WRITE(*,*) istat
    CALL error_out
  END IF
!      WRITE(*,*) 'Checking this is working.'

  READ(4) header
  READ(4) tri_num

  num_tris(i) = tri_num

  stl_tris = stl_tris + tri_num
  WRITE(11,301) geom_filenames(i)
  WRITE(11,201) i,num_tris(i)

  CLOSE(4)
END DO


total_tris = stl_tris + cgns_tris

WRITE(11,200) total_tris
WRITE(11,*) ''
!
! Now open them again to read everything
!
!
      IF (mixed_start) THEN
! Must get the total number of triangles
!        total_tris = cgns_tris + stl_tris
!
! Transfer the cgns values into a holding array, reallocate the array
! then transfer them back
!
  ALLOCATE(vertex_holder(3,3,cgns_tris))
  vertex_holder = tri_vertices
  DEALLOCATE(tri_vertices)
  ALLOCATE(tri_vertices(3,3,total_tris))

  tri_vertices(1:3,1:3,cgns_tris+1:total_tris) = &
    vertex_holder(1:3,1:3,1:cgns_tris)
  DEALLOCATE(vertex_holder)
!
! Do the same with the normal values
!
  ALLOCATE(normal_holder(3,cgns_tris))
  normal_holder = normal
  DEALLOCATE(normal)
  ALLOCATE(normal(3,total_tris))

  normal(1:3,cgns_tris+1:total_tris) = &
    normal_holder(1:3,1:cgns_tris)
  DEALLOCATE(normal_holder)
!
! Do the same with the triangle meaning
!
  ALLOCATE(tri_meaning(cgns_tris))
  meaning_holder = tri_meaning
  DEALLOCATE(tri_meaning)
  ALLOCATE(tri_meaning(total_tris))

  tri_meaning(cgns_tris+1:total_tris) = &
    meaning_holder(1:cgns_tris)
  DEALLOCATE(meaning_holder)
!
! If it's only a cgns start and not a mixed start, set the cgns
! values to the total_values
!
ELSE IF (cgns_start .AND. .NOT. mixed_start) THEN
!!        total_tris = cgns_tris
  RETURN
!!
!! If geometry file start, set the stl triangle total to the total
!! number of triangles
!!
!      ELSE
!        total_tris = stl_tris
END IF
!
! Allocate all the important stuff for the geometry
!
ALLOCATE(normal(3,total_tris))
ALLOCATE(tri_vertices(3,3,total_tris))
ALLOCATE(tri_meaning(total_tris))
!      total_tris = 1

tri_counter = 0

geom_file_open: DO i=1,num_geom_files
  OPEN(4,FILE=geom_filenames(i),STATUS='OLD',ACTION='READ', &
     FORM='UNFORMATTED',ACCESS='STREAM', IOSTAT=istat)
  IF (istat /= 0) THEN
    WRITE(*,*) 'Problem loading geometry file to read triangles.'
    WRITE(11,*) 'Problem loading geometry file to read triangles.'
    WRITE(*,*) istat
    CALL error_out
  END IF

  IF (units(i) == 1) THEN
    WRITE(11,*) 'Units in meters, no adjustment needed.'
    WRITE(11,*) ''
    adjusted_units = 1.0D0
  ELSE IF (units(i) == 2) THEN
    WRITE(11,*) 'Units in millimeters, need to divide vertex array ',&
      'by 1000.'
    adjusted_units = 1000.0D0
  ELSE IF (units(i) == 3) THEN
    WRITE(11,*) 'Units in centimeters, need to divide vertex array ',&
      'by 100.'
    adjusted_units = 100.0D0
  ELSE IF (units(i) == 4) THEN
    WRITE(11,*) 'Units in inches, need to divide vertex array ',&
      'by 39.37'
    adjusted_units = 39.37D0
  ELSE IF (units(i) == 5) THEN
    WRITE(11,*) 'Units in feet, need to divide vertex array ',&
      'by 3.28'
    adjusted_units = 3.28D0
  END IF


!
! Normal vector order, (x/y/z,triangle number)
! Vertex order, (x/y/z,vertex of triangle,triangle number)
! Bytes is useless for now
!
READ(4) junk1
READ(4) junk2

triangle_read: DO j = 1,num_tris(i)
  tri_counter = tri_counter +1

  READ(4) normal1,normal2,normal3
  READ(4) vertex1x1,vertex1x2,vertex1x3 !vertex 1, x,y, and z
  READ(4) vertex2x1,vertex2x2,vertex2x3 !vertex 2, x,y, and z
  READ(4) vertex3x1,vertex3x2,vertex3x3 !vertex 3, x,y, and z

!          READ(4) normal(1,tri_counter),normal(2,tri_counter),&
!            normal(3,tri_counter)
!          READ(4) vertex(1,1,tri_counter),vertex(2,1,tri_counter),&
!            vertex(3,1,tri_counter)
!          READ(4) vertex(1,2,tri_counter),vertex(2,2,tri_counter),&
!            vertex(3,2,tri_counter)
!          READ(4) vertex(1,3,tri_counter),vertex(2,3,tri_counter),&
!            vertex(3,3,tri_counter)
    READ(4) bytes
!          WRITE(*,*) tri_counter

    normal(1,tri_counter) = DBLE(normal1)
    normal(2,tri_counter) = DBLE(normal2)
    normal(3,tri_counter) = DBLE(normal3)
    tri_vertices(1,1,tri_counter) = DBLE(vertex1x1)/adjusted_units
    tri_vertices(2,1,tri_counter) = DBLE(vertex1x2)/adjusted_units
    tri_vertices(3,1,tri_counter) = DBLE(vertex1x3)/adjusted_units
    tri_vertices(1,2,tri_counter) = DBLE(vertex2x1)/adjusted_units
    tri_vertices(2,2,tri_counter) = DBLE(vertex2x2)/adjusted_units
    tri_vertices(3,2,tri_counter) = DBLE(vertex2x3)/adjusted_units
    tri_vertices(1,3,tri_counter) = DBLE(vertex3x1)/adjusted_units
    tri_vertices(2,3,tri_counter) = DBLE(vertex3x2)/adjusted_units
    tri_vertices(3,3,tri_counter) = DBLE(vertex3x3)/adjusted_units

    tri_meaning(tri_counter) = 1

  END DO triangle_read

  CLOSE(4)
!      WRITE(*,*) vertex1x1,vertex1x2,vertex1x3
!      WRITE(*,*) vertex(1,1,tri_counter),vertex(2,1,tri_counter),vertex(3,1,tri_counter)

WRITE(11,*) ''
END DO geom_file_open

!      IF (mixed_start) THEN
!        CALL geometry_unify
!      END IF

CLOSE(11)

 200  FORMAT ('Total number of tris in all geometry file = ',&
        I8)
 201  FORMAT ('Total number of tris in geometry file ',I8,' = ',&
          I8)
 301  FORMAT ('Reading in file ',A80)
END SUBROUTINE
