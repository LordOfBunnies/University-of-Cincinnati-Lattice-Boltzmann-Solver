      MODULE periodic_data
!
! Stores data for use with the periodic boundary condition
!
! The periodic pairs will be like
!
!
!
!
      USE precise
      USE constants
      IMPLICIT NONE
      SAVE
      INTEGER :: paired_periodic_bcs(2,3),periodic_sums,&
        periodic_flow_dir(6)
      INTEGER,ALLOCATABLE :: periodic_pairs(:,:,:)

      END MODULE
