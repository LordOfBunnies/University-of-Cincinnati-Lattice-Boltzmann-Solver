SUBROUTINE press_temp_bc(fout_bc,gout_bc,loc_dir,bc_number,inlet)
!
! This handles the pressure/temperature boundary condition. The user
! specifies the pressure and temperature and this calculates the density
! and the velocity from the pressure drop into the flow volume
!
! Called by: inlet_bc,outlet_bc
! Calls: error_out
!
USE nml_inlet_outlet
use inlet_outlet_elbm_data
USE linkwise
USE freestream_values
USE constants
USE precise
IMPLICIT NONE
INTEGER,INTENT(IN) :: loc_dir,bc_number!,link_num
INTEGER :: opp_dir
REAL(KIND=dp) :: u_velocity,v_velocity,w_velocity,&
  vel_mag,cu,bc_rho,bernoulli,node_press
real(kind=dp) :: fout_bc,gout_bc
!      REAL(KIND = dp),INTENT(OUT) ::
LOGICAL :: inlet

!
fout_bc = inl_fieq(opp(loc_dir),bc_number)
gout_bc = inl_gieq(opp(loc_dir),bc_number)
!
END SUBROUTINE

!IF (inlet) THEN
!  bc_rho = inlet_flow_val(2,bc_number)/(gas_const_R*&
!    inlet_flow_val(1,bc_number))
!
!  node_press = rho(link_num)*gas_const_R*temperature
!
!  vel_mag = (u_vel(link_num)**2.0D0 + v_vel(link_num)**2.0D0+&
!    w_vel(link_num)**2.0D0)
!
!  bernoulli = SQRT((2.0D0/bc_rho) * (node_press + 0.5D0*rho(link_num)*&
!    vel_mag**2.0D0 - node_press))
!
!  IF (inlet_type(2,bc_number) == 1) THEN
!    u_velocity = bernoulli
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 2) THEN
!    u_velocity = 0.0D0
!    v_velocity = bernoulli
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 3) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = bernoulli
!  ELSE IF (inlet_type(2,bc_number) == 4) THEN
!    u_velocity = -bernoulli
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 5) THEN
!    u_velocity = 0.0D0
!    v_velocity = -bernoulli
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 6) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = -bernoulli
!  END IF
!
!  cu = 3.0D0*(cx(loc_dir)*u_velocity+cy(loc_dir)*v_velocity+&
!       cz(loc_dir)*w_velocity)
!
!  fout_bc = wi(loc_dir)*bc_rho*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)
!!
!! For outlets
!!
!ELSE
!  bc_rho = outlet_flow_val(2,bc_number)/(gas_const_R*&
!    outlet_flow_val(1,bc_number))
!
!  node_press = rho(link_num)*gas_const_R*temperature
!
!  vel_mag = (u_vel(link_num)**2.0D0 + v_vel(link_num)**2.0D0 + &
!    w_vel(link_num)**2.0D0)
!
!  bernoulli = SQRT((2.0D0/bc_rho) * (node_press + 0.5D0*rho(link_num)*&
!    vel_mag**2.0D0 - node_press))
!
!  IF (inlet_type(2,bc_number) == 1) THEN
!    u_velocity = bernoulli
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 2) THEN
!    u_velocity = 0.0D0
!    v_velocity = bernoulli
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 3) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = bernoulli
!  ELSE IF (inlet_type(2,bc_number) == 4) THEN
!    u_velocity = -bernoulli
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 5) THEN
!    u_velocity = 0.0D0
!    v_velocity = -bernoulli
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 6) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = -bernoulli
!  END IF
!
!  cu = 3.0D0*(cx(opp_dir)*u_velocity+cy(opp_dir)*v_velocity+&
!       cz(opp_dir)*w_velocity)
!
!  fout_bc = wi(opp_dir)*bc_rho*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)
!
!END IF

! This will use Bernoulli's Equation for now because it's effectively
! incompressible and I'll come back later to add more detailed stuff
! since eventually this will be running compressible flow. XXXX
!
! This part is for inlet
!
!opp_dir = opp(loc_dir)
! XXXXXX not the simplest thing, needs to be solve because pressure boundaries are common
!fout_bc =
!gout_bc =
!
