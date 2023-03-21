      MODULE geom_data
!
! Holds the pertinent information for the geometry
!
!
      USE precise
      IMPLICIT NONE
      SAVE
      REAL(KIND=dp), ALLOCATABLE :: tri_vertices(:,:,:),normal(:,:),nodes_x(:),&
        nodes_y(:),nodes_z(:),normal_holder(:,:),vertex_holder(:,:,:)
      INTEGER,ALLOCATABLE :: tri_meaning(:),meaning_holder(:)
      INTEGER :: total_tris,cgns_tris,stl_tris

      END MODULE
