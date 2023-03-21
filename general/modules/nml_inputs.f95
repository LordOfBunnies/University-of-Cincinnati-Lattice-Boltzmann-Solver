      MODULE nml_inputs
!
! This stores information like cgns load data
!
!
!
!
      USE precise
      IMPLICIT NONE
      SAVE
      INTEGER :: num_cgns_files
      INTEGER,ALLOCATABLE :: cgns_units(:)
      CHARACTER(len=64),ALLOCATABLE :: cgns_input_files(:)




      END MODULE 
