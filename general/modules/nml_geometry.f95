      MODULE nml_geometry
!
! Saved information about the geometry
!
!
      IMPLICIT NONE
      SAVE
      INTEGER :: num_geom_files
      CHARACTER(LEN=80),ALLOCATABLE :: geom_filenames(:)
      INTEGER,ALLOCATABLE :: units(:),num_tris(:)


      END MODULE
