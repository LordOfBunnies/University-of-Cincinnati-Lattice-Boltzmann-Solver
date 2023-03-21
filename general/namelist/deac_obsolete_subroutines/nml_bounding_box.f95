      MODULE nml_bounding_box
!
! Saved information about the bounding box.
!
! bounding_method defines the shape
! bounds defines the limits of the box
!
      USE precise
      IMPLICIT NONE
      SAVE
      INTEGER :: bounding_method
      REAL(KIND=dp),ALLOCATABLE :: bounds(:)
   
      END MODULE
