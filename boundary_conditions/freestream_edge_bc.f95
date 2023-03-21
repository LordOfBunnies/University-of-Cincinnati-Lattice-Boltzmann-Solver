SUBROUTINE freestream_edge_bc(fout_bc,gout_bc,loc_dir,loc_state)
!
! This deals with the freestream boundaries slightly differently than
! freestream_bc does.  This is because it
!
! Called by: loop
! Calls:
!
USE freestream_values
USE constants
USE linkwise
!USE loop_data
USE precise
IMPLICIT NONE
INTEGER,INTENT(IN) :: loc_dir,loc_state
REAL(KIND = dp),INTENT(OUT) :: fout_bc,gout_bc

fout_bc = fs_fieq(opp(loc_dir))
gout_bc = fs_gieq(opp(loc_dir))

END SUBROUTINE

!opp_dir = opp(loc_dir)
!
!loc_rho = rho_inf
!IF (primary_flow_dir == 1) THEN
!  u_velocity = v_inf
!  v_velocity = 0.0D0
!  w_velocity = 0.0D0
!ELSE IF (primary_flow_dir == 2) THEN
!  u_velocity = 0.0D0
!  v_velocity = v_inf
!  w_velocity = 0.0D0
!ELSE IF (primary_flow_dir == 3) THEN
!  u_velocity = 0.0D0
!  v_velocity = 0.0D0
!  w_velocity = v_inf
!ELSE IF (primary_flow_dir == 4) THEN
!  u_velocity = -v_inf
!  v_velocity = 0.0D0
!  w_velocity = 0.0D0
!ELSE IF (primary_flow_dir == 5) THEN
!  u_velocity = 0.0D0
!  v_velocity = -v_inf
!  w_velocity = 0.0D0
!ELSE IF (primary_flow_dir == 6) THEN
!  u_velocity = 0.0D0
!  v_velocity = 0.0D0
!  w_velocity = -v_inf
!
!END IF
!
!cu = 3.0D0*(cx(opp_dir)*u_velocity+cy(opp_dir)*v_velocity+&
!       cz(opp_dir)*w_velocity)
!vel_mag = (u_velocity**2.0D0+v_velocity**2.0D0+w_velocity**2.0D0)
!fout_bc = wi(opp_dir)*loc_rho*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)
!      WRITE(*,*) cu,vel_mag,fout_bc,u_velocity,opp_dir
