MODULE nml_init
!
! Saved information about initializations the user wants to do.
!
! init_bounds define the limits of the box
! init_info defines the shape and priority of the init region
! init_data defines what date the region will be starting with
! 
USE precise
IMPLICIT NONE
SAVE
INTEGER :: num_init_regions
REAL(KIND=dp),ALLOCATABLE :: init_bounds(:,:),init_data(:,:)
INTEGER,ALLOCATABLE :: init_info(:,:)

END MODULE
