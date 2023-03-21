SUBROUTINE outlet_bc(fout_outl,fout_opp,gout_outl,gout_opp,loc_dir,&
  outl_link,opp_state,away_fi,away_gi,a,b,c)
!
! Deal with outlet boundary conditions
!
! Called by: loop
! Calls: rho_v_bc,press_temp_bc,freestream_bc,v_temp_bc,periodic_bc
!
USE nml_inlet_outlet
USE precise
use linkwise
IMPLICIT NONE
INTEGER :: loc_dir,link_num,outl_link,outl_number
integer,intent(in) :: a,b,c,opp_state
integer :: pass_dir
REAL(KIND=dp),INTENT(INOUT) :: fout_outl,gout_outl,fout_opp,gout_opp !,outlet_stuff
real(kind=dp),intent(in) :: away_fi,away_gi
LOGICAL :: inlet

inlet = .FALSE.
!outl_number = ABS(outl_link)-200

if (outl_link >= -299 .and. outl_link <= -200) then
  outl_number = ABS(outl_link)-200
  pass_dir = loc_dir
else if (outl_link >= -499 .and. outl_link <= -400) then
  outl_number = ABS(outl_link)-400
  pass_dir = opp(loc_dir)
end if
!
! Call for a rho-V outlet
!
IF (outlet_type(1,outl_number) == 1) THEN
  CALL rho_v_bc(fout_outl,gout_outl,pass_dir,outl_number,inlet)

!
! Call for a pressure-temperature outlet
!
ELSE IF (outlet_type(1,outl_number) == 2) THEN
  CALL press_temp_bc(fout_outl,gout_outl,pass_dir,outl_number,inlet)

!
! Call for a freestream outlet
!
ELSE IF (outlet_type(1,outl_number) == 3) THEN
  CALL freestream_bc(fout_outl,gout_outl,inlet,pass_dir,a,b,c)

!
! Call for a velocity-temperature outlet
!
ELSE IF (outlet_type(1,outl_number) == 4) THEN
  CALL v_temp_bc(fout_outl,gout_outl,pass_dir,outl_number,inlet)

else if (outlet_type(1,outl_number) == 6) THEN
  CALL pure_outflow_bc(fout_outl,fout_opp,gout_outl,gout_opp,pass_dir,outl_number)
!      ELSE IF (outlet_type(1,outl_number) == 5) THEN
!        CALL periodic_bc(fout_outl,loc_dir,outl_number,inlet)
END IF

END SUBROUTINE
