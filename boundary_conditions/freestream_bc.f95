SUBROUTINE freestream_bc(fout_bc,gout_bc,inlet,loc_dir,a,b,c)
!
! Deals with making the fout values appropriate for those coming from
! the freestream
!
! Called by: loop,inlet_bc,outlet_bc
! Calls:
!
USE freestream_values
USE constants
USE precise
USE nml_inlet_outlet
USE linkwise
use amr_info_holder
IMPLICIT NONE
INTEGER :: loc_dir
real(kind=dp),intent(out) :: fout_bc,gout_bc
integer, intent(in) :: a,b,c
logical :: inlet

!if (inlet) then
fout_bc = fs_fieq(opp(loc_dir))
gout_bc = fs_gieq(opp(loc_dir))
!write(*,*) 'freestream stuff',fs_fieq
!write(*,*) 'grog',fs_gieq
!else
!
!  if (state(a+cx(opp(loc_dir)),b+cy(opp(loc_dir)),c+cz(opp(loc_dir)),1) >= 0) then
!    fout_bc = fi(a+cx(opp(loc_dir)),b+cy(opp(loc_dir)),c+cz(opp(loc_dir)),opp(loc_dir))
!    gout_bc = gi(a+cx(opp(loc_dir)),b+cy(opp(loc_dir)),c+cz(opp(loc_dir)),opp(loc_dir))
!  else
!    fout_bc = fi(a,b,c,opp(loc_dir))
!    gout_bc = gi(a,b,c,opp(loc_dir))
!  end if
!
!end if

END SUBROUTINE


!REAL(KIND=dp) :: loc_rho,u_velocity,v_velocity,w_velocity,&
!  vel_mag,cu
!      opp_dir = opp(loc_dir)
!      loc_rho = rho_inf
!
!
!
! Inlet
!
!
!IF (inlet) THEN
!
!  IF (inlet_type(2,bc_number) == 1) THEN
!    u_velocity = v_inf
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 2) THEN
!    u_velocity = 0.0D0
!    v_velocity = v_inf
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 3) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = v_inf
!  ELSE IF (inlet_type(2,bc_number) == 4) THEN
!    u_velocity = -v_inf
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 5) THEN
!    u_velocity = 0.0D0
!    v_velocity = -v_inf
!    w_velocity = 0.0D0
!  ELSE IF (inlet_type(2,bc_number) == 6) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = -v_inf
!  END IF
!!
!!
!! Outlet
!!
!!
!ELSE
!  IF (outlet_type(2,bc_number) == 1) THEN
!    u_velocity = -v_inf
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 2) THEN
!    u_velocity = 0.0D0
!    v_velocity = -v_inf
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 3) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = -v_inf
!  ELSE IF (outlet_type(2,bc_number) == 4) THEN
!    u_velocity = v_inf
!    v_velocity = 0.0D0
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 5) THEN
!    u_velocity = 0.0D0
!    v_velocity = v_inf
!    w_velocity = 0.0D0
!  ELSE IF (outlet_type(2,bc_number) == 6) THEN
!    u_velocity = 0.0D0
!    v_velocity = 0.0D0
!    w_velocity = v_inf
!  END IF
!END IF
!cu = 3.0D0*(cx(opp_dir)*u_velocity+cy(opp_dir)*v_velocity+&
!       cz(opp_dir)*w_velocity)
!vel_mag = (u_velocity**2.0D0+v_velocity**2.0D0+w_velocity**2.0D0)
!fout_bc = wi(opp_dir)*loc_rho*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)

