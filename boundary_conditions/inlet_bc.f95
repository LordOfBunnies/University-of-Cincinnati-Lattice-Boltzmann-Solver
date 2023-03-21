SUBROUTINE inlet_bc(fout_inl,gout_inl,loc_dir,inl_link)
!
! This is for dealing with inlet boundary conditions.  Essentially,
! this takes the opposite direction from what was being used in
! loop and calculates the values for it from the freestream
!
! Called by: loop
! Calls: rho_v_bc,press_temp_bc,freestream_bc,v_temp_bc,periodic_bc
!
USE nml_inlet_outlet
USE linkwise
USE precise
IMPLICIT NONE
INTEGER :: loc_dir,link_num,inl_link,inl_number
integer :: pass_dir
REAL(KIND=dp),INTENT(OUT) :: fout_inl,gout_inl
LOGICAL :: inlet

inlet = .TRUE.

if (inl_link >= -199 .and. inl_link <= -100) then
  inl_number = ABS(inl_link)-100
  pass_dir = loc_dir
else if (inl_link >= -399 .and. inl_link <= -300) then
  inl_number = ABS(inl_link)-300
  pass_dir = opp(loc_dir)
end if
!      WRITE(*,*) inl_number

!
! Call for a rho-V outlet
!
IF (inlet_type(1,inl_number) == 1) THEN
  CALL rho_v_bc(fout_inl,gout_inl,pass_dir,inl_number,inlet)
!
! Call for a pressure-temperature outlet
!
ELSE IF (inlet_type(1,inl_number) == 2) THEN
  CALL press_temp_bc(fout_inl,gout_inl,pass_dir,inl_number,inlet)
!
! Call for a freestream outlet
!
ELSE IF (inlet_type(1,inl_number) == 3) THEN
  CALL freestream_bc(fout_inl,gout_inl,inlet,pass_dir)
!        WRITE(*,*) 'Goof!'
!
! Call for a velocity-temperature outlet
!
ELSE IF (inlet_type(1,inl_number) == 4) THEN
  CALL v_temp_bc(fout_inl,gout_inl,pass_dir,inl_number,inlet)
!
!
END IF

END SUBROUTINE
