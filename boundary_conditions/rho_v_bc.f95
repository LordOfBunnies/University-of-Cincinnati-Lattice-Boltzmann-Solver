SUBROUTINE rho_v_bc(fout_bc,gout_bc,loc_dir,bc_number,inlet)
!
! This deals with a momentum inlet where rho and V are specified
!
! Called by: inlet_bc,outlet_bc
! Calls:
!
USE nml_inlet_outlet
use inlet_outlet_elbm_data, only: inl_fieq,inl_gieq,outl_fieq,outl_gieq
USE linkwise
USE freestream_values
USE constants
USE precise
IMPLICIT NONE
REAL(KIND=dp),INTENT(INOUT) :: fout_bc,gout_bc
REAL(KIND=dp) :: cu,vel_mag,u_velocity,v_velocity,&
  w_velocity,loc_rho
INTEGER :: opp_dir
INTEGER,INTENT(IN) :: loc_dir,bc_number
LOGICAL,INTENT(IN) :: inlet

!opp_dir = opp(loc_dir)
if (inlet) then
  fout_bc = inl_fieq(opp(loc_dir),bc_number)
  gout_bc = inl_gieq(opp(loc_dir),bc_number)
else
  fout_bc = outl_fieq(opp(loc_dir),bc_number)
  gout_bc = outl_gieq(opp(loc_dir),bc_number)
end if

END SUBROUTINE
!IF (inlet) THEN
!  loc_rho = inlet_flow_val(2,bc_number)
!  IF (inlet_type(2,bc_number) == 1) THEN
!    u_velocity = inlet_flow_val(2,bc_number)
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 2) THEN
!    u_velocity = 0.0D0
!    v_velocity = inlet_flow_val(2,bc_number)
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 3) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = inlet_flow_val(2,bc_number)
!  ELSE IF (inlet_type(2,bc_number) == 4) THEN
!    u_velocity = inlet_flow_val(2,bc_number)
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 5) THEN
!    u_velocity = 0.0D0
!    v_velocity = inlet_flow_val(2,bc_number)
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 6) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = inlet_flow_val(2,bc_number)
!  END IF
!!
!! Outlets
!!
!ELSE
!  loc_rho = outlet_flow_val(2,bc_number)
!  IF (outlet_type(2,bc_number) == 1) THEN
!    u_velocity = outlet_flow_val(2,bc_number)
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 2) THEN
!    u_velocity = 0.0D0
!    v_velocity = outlet_flow_val(2,bc_number)
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 3) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = outlet_flow_val(2,bc_number)
!  ELSE IF (outlet_type(2,bc_number) == 4) THEN
!    u_velocity = outlet_flow_val(2,bc_number)
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 5) THEN
!    u_velocity = 0.0D0
!    v_velocity = outlet_flow_val(2,bc_number)
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 6) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = outlet_flow_val(2,bc_number)
!  END IF
!END IF

!
!cu = 3.0D0*(cx(opp_dir)*u_velocity + cy(opp_dir)*v_velocity +&
!       cz(opp_dir)*w_velocity )
!vel_mag = (u_velocity**2.0D0 + v_velocity**2.0D0 + w_velocity**2.0D0)
!fout_bc = wi(opp_dir)*loc_rho*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)

!      fout_bc = fi - tau/dt*(fi-feq_loc)


