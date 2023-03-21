MODULE nml_inlet_outlet
!
! Saves information about the inlets and outlets of the run
!
! inlet_type says what kind of inlet it is
! outlet_type says what kind of outlet it is
! inlet_flow_val defines the values for the inlet
! outlet_flow_val defines the values for the outlet
!
USE precise
IMPLICIT NONE
SAVE
INTEGER :: num_inlets,num_outlets

INTEGER,ALLOCATABLE :: inlet_type(:,:),outlet_type(:,:)
REAL(KIND=dp),ALLOCATABLE :: inlet_flow_val(:,:),outlet_flow_val(:,:)
CHARACTER(len=80),ALLOCATABLE :: in_name(:),out_name(:)
LOGICAL :: inlet_side(6),outlet_side(6)

END MODULE
