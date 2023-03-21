SUBROUTINE calculate_lagrangian_mults(loc_rho,loc_u,loc_v,&
  loc_w, loc_temp, loc_pxx, loc_pyy, loc_pzz, loc_pxy, loc_pxz, loc_pyz,&
  loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
  loc_pi_xy,loc_pi_xz,loc_pi_yz,&
  loc_lambda_x,loc_lambda_y,loc_lambda_z,loc_fieq,a,b,c)
!
! Calculates the Langrangian multipliers for every location in the grid
!
!
! Called by: collision
! Calls:
!
use freestream_values
USE amrex_base_module
USE amrex_amr_module
USE precise
USE amr_info_holder
use amr_processes!, only: fieq_comp_19,fieq_comp_39,self,fieq_comp_39_shifted,self
use constants
use grid_data
use linkwise
IMPLICIT NONE
INTEGER :: a,b,c,i,en,em,loop_counter,smooth_counter
REAL(KIND=dp),INTENT(IN) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,&
  loc_pxy,loc_pxz,loc_pyz
real(kind=dp) ::  loc_q_x,loc_q_y,loc_q_z
REAL(KIND=dp) :: vel_mag,max_diff
real(kind=dp) :: orig_rho,orig_u,orig_v,orig_w,diff_rho,diff_u,diff_v,diff_w,diff_pxx,&
  diff_pyy,diff_pzz,diff_pxy,diff_pxz,diff_pyz,diff_q_x,diff_q_y,diff_q_z,diff_vel_mag
REAL(kind=dp),INTENT(OUT) :: loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_lambda_x,loc_lambda_y,&
  loc_lambda_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz
!
real(kind=dp) :: sp_cx,sp_cy,sp_cz,sp_cxcy,sp_cxcz,sp_cycz,sp_cx2,sp_cy2,sp_cz2,sp_8dir_pos
real(kind=dp) :: sp_c2cx,sp_c2cy,sp_c2cz,sp_cx2cy,sp_cxcy2,sp_cx2cz,sp_cxcz2,sp_cy2cz,sp_cycz2
real(kind=dp) :: loc_fieq(dir)
REAL(kind=dp),ALLOCATABLE :: expo(:)
REAL(kind=dp),ALLOCATABLE :: jacobian(:,:)!,b(:)
!
! sp = save point
!
!TYPE(amrex_mfiter) :: mfi
!write(*,*) 'Beginning calculation of Lagrangian multipliers',self

max_diff = 1.0D0
vel_mag = (loc_u**2+loc_v**2+loc_w**2)

loc_q_x = (loc_rho*vel_mag + 5*loc_rho*loc_temp)*loc_u
loc_q_y = (loc_rho*vel_mag + 5*loc_rho*loc_temp)*loc_v
loc_q_z = (loc_rho*vel_mag + 5*loc_rho*loc_temp)*loc_w
!write(*,*) 'Values for solutions ',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,&
!  loc_pzz,loc_pxy,loc_pxz,loc_pyz,loc_q_x,loc_q_y,loc_q_z
loop_counter = 0

!if (a == 191 .and. b == 1 .and. c == 43) then
!  write(*,*) 'incoming values',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,&
!  loc_pxy,loc_pxz,loc_pyz,loc_q_x,loc_q_y,loc_q_z
!  write(*,*) 'incoming LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!  loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z
!end if

IF (dimensions == 2) THEN
  en = 8
  ALLOCATE(jacobian(en,en))
  ALLOCATE(expo(dir))

  jacobian = 0.0D0



  DEALLOCATE(jacobian)
  deallocate(expo)
!
!
! 3D section
!
!
  ELSE IF (dimensions == 3) THEN
    em = 13
    ALLOCATE(jacobian(em,em))
    ALLOCATE(expo(dir))



  if (dir == 19) then
    do while (max_diff > 1.0D-7)
!      write(*,*) 'local values',loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
!        loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
!        loc_lambda_z,self

      do i = 1,dir
        expo(i) = fieq_comp_19(loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
        loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
        loc_lambda_z,i)

        !write(*,*) 'fieq value',expo(i),' for direction ',i
      end do
      write(*,*) 'dummy/infinite loop check',expo(1),self
! Row 1, density derivatives, coefficient = 1
    jacobian(1,1) = SUM(expo) !1
    jacobian(1,2) = expo(2)-expo(5)+expo(8)-expo(10)+expo(12)+expo(13)-expo(14)-expo(15)+&
      expo(16)-expo(18)  !drho_dzeta_x, cx
    jacobian(1,3) = expo(3)-expo(6)+expo(8)+expo(9)+expo(10)+expo(11)-expo(16)-expo(17)-&
      expo(18)-expo(19)  !drho_dzeta_y, cy
    jacobian(1,4) = expo(4)-expo(7)-expo(9)+expo(11)+expo(12)-expo(13)-expo(14)+expo(15)-&
      expo(17)+expo(19)  !drho_dpi_z, cz
    jacobian(1,5) = expo(2)+expo(5)+expo(8)+expo(10)+expo(12)+expo(13)+expo(14)+expo(15)+&
      expo(16)+expo(18)  !drho_dpi_xx, cx^2
    jacobian(1,6) = expo(3)+expo(6)+expo(8)+expo(9)+expo(10)+expo(11)+expo(16)+expo(17)+&
      expo(18)+expo(19)  !drho_dpi_yy, cy^2
    jacobian(1,7) = expo(4)+expo(7)+expo(9)+expo(11)+expo(12)+expo(13)+expo(14)+expo(15)+&
      expo(17)+expo(19) !drho_dpi_zz, cz^2
    jacobian(1,8) = expo(8)-expo(10)-expo(16)+expo(18) !drho_dpi_xy, cx*cy
    jacobian(1,9) = expo(12)-expo(13)+expo(14)-expo(15) !drho_dpi_xz, cx*cz
    jacobian(1,10) = expo(11)-expo(9)+expo(17)-expo(19) !drho_dpi_yz, cy*cz 
    jacobian(1,11) = expo(2)-expo(5)+2.0D0*(expo(8)-expo(10)+expo(12)+expo(13)-expo(14)-expo(15)+&
      expo(16)-expo(18)) !drho_dlambda_x, cx*c^2
    jacobian(1,12) = expo(3)-expo(6)+2.0D0*(expo(8)+expo(9)+expo(10)+expo(11)-expo(16)-expo(17)-&
      expo(18)-expo(19)) !drho_dlambda_z, cy*c^2
    jacobian(1,13) = expo(4)-expo(7)+2.0D0*(expo(11)-expo(9)+expo(12)-expo(13)-expo(14)+expo(15)-&
      expo(17)+expo(19))!drho_dlambda_z, cz*c^2
!
! Row 2, x-momentum derivatives, cx
!
    jacobian(2,1) = jacobian(1,2) !cx
    jacobian(2,2) = jacobian(1,5) !cx^2
    jacobian(2,3) = jacobian(1,8) !cxcy
    jacobian(2,4) = jacobian(1,9) !cxcz
    jacobian(2,5) = jacobian(1,2) !cx^3
!
    jacobian(2,6) = expo(8)-expo(10)+expo(16)-expo(18) !cx*cy^2
!
    jacobian(2,7) = expo(12)+expo(13)-expo(14)-expo(15)!cx*cz^2
!
    jacobian(2,8) = expo(8)+expo(10)-expo(16)-expo(18)!cx^2*cy
!
    jacobian(2,9) = expo(12)-expo(13)-expo(14)+expo(15)!cx^2*cz
!
    jacobian(2,10) = 0.0D0 !cx*cy*cz
!
    jacobian(2,11) =2.0D0*jacobian(1,5)-expo(2)-expo(5)!cx^2*c^2

    jacobian(2,12) = 2.0D0*jacobian(1,8) !cx*cy*c^2

    jacobian(2,13) = 2.0D0*jacobian(1,9) !cx*cz*c^2
! Row 3, y-momentum derivatives, cy
    jacobian(3,1) = jacobian(1,3) !cy
    jacobian(3,2) = jacobian(1,8) !cx*cy
    jacobian(3,3) = jacobian(1,6) !cy^2
    jacobian(3,4) = jacobian(1,10) !cy*cz
    jacobian(3,5) = jacobian(2,8) !cx^2*cy
    jacobian(3,6) = jacobian(1,3) !cy^3
    jacobian(3,7) = expo(9)+expo(11)-expo(17)-expo(19) !cy*cz^2
    jacobian(3,8) = jacobian(2,6) !cx*cy^2
    jacobian(3,9) = jacobian(2,10) !cx*cy*cz
    jacobian(3,10) = expo(11)-expo(9)-expo(17)+expo(19) !cy^2*cz
    jacobian(3,11) = jacobian(2,12) !cx*cy*c^2
    jacobian(3,12) = 2.0D0*jacobian(1,6)-expo(3)-expo(6) !cy^2*c^2
    jacobian(3,13) = 2.0D0*jacobian(1,10) !cy*cz*c^2
! Row 4, z-momentum derivatives, cz
    jacobian(4,1) = jacobian(1,4) !cz
    jacobian(4,2) = jacobian(2,4) !cx*cz
    jacobian(4,3) = jacobian(3,4) !cy*cz
    jacobian(4,4) = jacobian(1,7) !cz^2
    jacobian(4,5) = jacobian(2,9) !cx^2*cz
    jacobian(4,6) = jacobian(3,10) !cy^2*cz
    jacobian(4,7) = jacobian(1,4) !cz^3
    jacobian(4,8) = jacobian(3,9) !cx*cy*cz
    jacobian(4,9) = jacobian(2,7) !cx*cz^2
    jacobian(4,10) = jacobian(3,7) !cy*cz^2
    jacobian(4,11) = jacobian(2,13) !cx*cz*c^2
    jacobian(4,12) = jacobian(3,13) !cy*cz*c^2
    jacobian(4,13) = 2.0D0*jacobian(1,7)-expo(4)-expo(7) !cz^2*c^2
! Row 5, x-directional pressure, cx^2
    jacobian(5,1) = jacobian(1,5) !cx^2
    jacobian(5,2) = jacobian(2,5) !cx^3
    jacobian(5,3) = jacobian(3,5) !cx^2*cy or (2,8)
    jacobian(5,4) = jacobian(4,5) !cx^2*cz
    jacobian(5,5) = jacobian(1,5) !cx^4
    jacobian(5,6) = expo(8)+expo(10)+expo(16)+expo(18) !cx^2*cy^2
    jacobian(5,7) = expo(12)+expo(13)+expo(14)+expo(15) !cx^2*cz^2
    jacobian(5,8) = jacobian(3,4) !cx^3*cy
    jacobian(5,9) = jacobian(2,4) ! cx^3*cz
    jacobian(5,10) = 0.0D0 !cx^2*cy*cz
    jacobian(5,11) = jacobian(1,11) !cx^3*c^2
    jacobian(5,12) = 2.0D0*jacobian(2,8) !cx^2*cy*c^2
    jacobian(5,13) = 2.0D0*jacobian(4,5) !cx^2*cz*c^2
! Row 6, y-directional pressure, cy^2
    jacobian(6,1) = jacobian(1,6) !cy^2
    jacobian(6,2) = jacobian(2,6) !cx*cy^2
    jacobian(6,3) = jacobian(3,6) !cy^3
    jacobian(6,4) = jacobian(4,6) !cy^2*cz
    jacobian(6,5) = jacobian(5,6) !cx^2*cy^2
    jacobian(6,6) = jacobian(1,6) !cy^4
    jacobian(6,7) = expo(9)+expo(11)+expo(17)+expo(19)  !cy^2*cz^2
    jacobian(6,8) = jacobian(1,8) !cy^3*cx
    jacobian(6,9) = 0.0D0 !cx*cy^2*cz
    jacobian(6,10) = jacobian(1,10) !cy^3*cz
    jacobian(6,11) = 2.0D0*jacobian(2,6) !cx*cy^2*c^2
    jacobian(6,12) = jacobian(1,12) !cy^3*c^2
    jacobian(6,13) = 2.0D0*jacobian(4,6) !cz*cy^2*c^2
! Row 7, z-directional pressure, cz^2
    jacobian(7,1) = jacobian(1,7) !cz^2
    jacobian(7,2) = jacobian(2,7) !cx*cz^2
    jacobian(7,3) = jacobian(3,7) !cy*cz^2
    jacobian(7,4) = jacobian(4,7) !cz^3
    jacobian(7,5) = jacobian(5,7) !cx^2*cz^2
    jacobian(7,6) = jacobian(6,7) !cy^2*cz^2
    jacobian(7,7) = jacobian(1,7) !cz^4
    jacobian(7,8) = 0.0D0  !cx*cy*cz^2
    jacobian(7,9) = jacobian(2,4) !cx*cz^3
    jacobian(7,10) = jacobian(3,4) !cy*cz^3
    jacobian(7,11) = 2.0D0*jacobian(2,7) !cx*cz^2*c^2
    jacobian(7,12) = 2.0D0*jacobian(3,7)  !cy*cz^2*c^2
    jacobian(7,13) = jacobian(1,13) !cz^3*c^2
! Row 8, xy-shear, cx*cy
    jacobian(8,1) = jacobian(1,8) !cx*cy
    jacobian(8,2) = jacobian(2,8) !cx^2*cy
    jacobian(8,3) = jacobian(3,8) !cx*cy^2
    jacobian(8,4) = jacobian(4,8) !cx*cy*cz
    jacobian(8,5) = jacobian(5,8) !cx^3*cy
    jacobian(8,6) = jacobian(6,8) !cx*cy^3
    jacobian(8,7) = jacobian(7,8) !cx^2*cy^2
    jacobian(8,8) = jacobian(5,6)  !cx^2*cy^2
    jacobian(8,9) = jacobian(5,10) !cx^2*cy*cz
    jacobian(8,10) = jacobian(6,9) !cx*cy^2*cz
    jacobian(8,11) = jacobian(5,12) !cx^2*cy*c^2
    jacobian(8,12) = jacobian(6,11) !cx*cy^2*c^2
    jacobian(8,13) = 0.0D0 !cx*cy*cz*c^2
! Row 9, xz-shear, cx*cz
    jacobian(9,1) = jacobian(1,9) !cx*cz
    jacobian(9,2) = jacobian(2,9) !cx^2*cz
    jacobian(9,3) = jacobian(8,4) !cx*cy*cz
    jacobian(9,4) = jacobian(4,9) !cx*cz^2
    jacobian(9,5) = jacobian(5,9) !cx^3*cz
    jacobian(9,6) = jacobian(6,9) !cx*cy^2*cz
    jacobian(9,7) = jacobian(7,9) !cx*cz^3
    jacobian(9,8) = jacobian(8,9) !cx^2*cy*cz
    jacobian(9,9) = jacobian(5,7) !cx^2*cz^2
    jacobian(9,10) = jacobian(8,7) !cx*cy*cz^2
    jacobian(9,11) = jacobian(5,13) !cx^2*cz*c^2
    jacobian(9,12) = jacobian(8,13) !cx*cy*cz*c^2
    jacobian(9,13) = jacobian(7,11) !cx*cz^2*c^2
! Row 10, yz-shear, cy*cz
    jacobian(10,1) = jacobian(1,10) !cy*cz
    jacobian(10,2) = jacobian(2,10) !cx*cy*cz
    jacobian(10,3) = jacobian(3,10) !cy^2*cz
    jacobian(10,4) = jacobian(4,10) !cy*cz^2
    jacobian(10,5) = jacobian(5,10) !cx^2*cy*cz
    jacobian(10,6) = jacobian(6,10) !cy^3*cz
    jacobian(10,7) = jacobian(7,10) !cy*cz^3
    jacobian(10,8) = jacobian(8,10) !cx*cy^2*cz
    jacobian(10,9) = jacobian(9,10) !cx*cy*cz^2
    jacobian(10,10) = jacobian(6,7) !cy^2*cz^2
    jacobian(10,11) = jacobian(9,12) !cx*cy*cz*c^2
    jacobian(10,12) = jacobian(6,13) !cy^2*cz*c^2
    jacobian(10,13) = jacobian(7,12) !cy*cz^2*c^2
! Row 11, x-direction contracted heat flux, cx*c^2
    jacobian(11,1) = jacobian(1,11) !cx*c^2
    jacobian(11,2) = jacobian(2,11) !cx^2*c^2
    jacobian(11,3) = jacobian(3,11) !cx*cy*c^2
    jacobian(11,4) = jacobian(4,11) !cx*cz*c^2
    jacobian(11,5) = jacobian(5,11) !cx^3*cy*c^2
    jacobian(11,6) = jacobian(6,11) !cx*cy^2*c^2
    jacobian(11,7) = jacobian(7,11) !cx*cz^2*c^2
    jacobian(11,8) = jacobian(8,11) !cx^2*cy*c^2
    jacobian(11,9) = jacobian(9,11) !cx^2*cz*c^2
    jacobian(11,10) = jacobian(10,11) !cx*cy*cz*c^2
    jacobian(11,11) = 2.0D0*jacobian(2,11)-expo(2)-expo(5) !cx^2*c^4
    jacobian(11,12) = 2.0D0*jacobian(2,12) !cx*cy*c^4
    jacobian(11,13) = 2.0D0*jacobian(2,13) !cx*cz*c^4
! Row 12, y-direction contracted heat flux, cy*c^2
    jacobian(12,1) = jacobian(1,12) !cy*c^2
    jacobian(12,2) = jacobian(2,12) !cx*cy*c^2
    jacobian(12,3) = jacobian(3,12) !cy^2*c^2
    jacobian(12,4) = jacobian(4,12) !cy*cz*c^2
    jacobian(12,5) = jacobian(5,12) !cx^2*cy*c^2
    jacobian(12,6) = jacobian(6,12) !cy^3*c^2
    jacobian(12,7) = jacobian(7,12) !cy*cz^2*c^2
    jacobian(12,8) = jacobian(8,12) !cx*cy^2*c^2
    jacobian(12,9) = jacobian(9,12) !cx*cy*cz*c^2
    jacobian(12,10) = jacobian(10,12) !cy^2*cz*c^2
    jacobian(12,11) = jacobian(11,12) !cx*cy*c^4
    jacobian(12,12) = 2.0D0*jacobian(3,12)-expo(3)-expo(6) !cy^2*c^4
    jacobian(12,13) = 2.0D0*jacobian(3,13) !cy*cz*c^4
! Row 13, z-direction contracted heat flux, cz*c^2
    jacobian(13,1) = jacobian(1,13) !cz*c^2
    jacobian(13,2) = jacobian(2,13) !cx*cz*c^2
    jacobian(13,3) = jacobian(3,13) !cy*cz*c^2
    jacobian(13,4) = jacobian(4,13) !cx^2*c^2
    jacobian(13,5) = jacobian(5,13) !cx^2*cz*c^2
    jacobian(13,6) = jacobian(6,13) !cy^2*cz*c^2
    jacobian(13,7) = jacobian(7,13) !cz^3*c^2
    jacobian(13,8) = jacobian(8,13) !cx*cz*cz*c^2
    jacobian(13,9) = jacobian(9,13) !cx*cz^2*c^2
    jacobian(13,10) = jacobian(10,13) !cy*cz^2*c^2
    jacobian(13,11) = jacobian(11,13) !cx*cz*c^4
    jacobian(13,12) = jacobian(12,13) !cy*cz*c^4
    jacobian(13,13) = 2.0D0*jacobian(4,13)-expo(4)-expo(7)  !cz^2*c^4
!
! Output to check if things are going nuts.  Comment out when no longer needed
!
!write(*,202) jacobian(1,1),jacobian(1,2),jacobian(1,3),jacobian(1,4),jacobian(1,5),jacobian(1,6),jacobian(1,7),&
!  jacobian(1,8),jacobian(1,9),jacobian(1,10),jacobian(1,11),jacobian(1,12),jacobian(1,13)
!write(*,202) jacobian(2,1),jacobian(2,2),jacobian(2,3),jacobian(2,4),jacobian(2,5),jacobian(2,6),jacobian(2,7),&
!  jacobian(2,8),jacobian(2,9),jacobian(2,10),jacobian(2,11),jacobian(2,12),jacobian(2,13)
!write(*,202) jacobian(3,1),jacobian(3,2),jacobian(3,3),jacobian(3,4),jacobian(3,5),jacobian(3,6),jacobian(3,7),&
!  jacobian(3,8),jacobian(3,9),jacobian(3,10),jacobian(3,11),jacobian(3,12),jacobian(3,13)
!write(*,202) jacobian(4,1),jacobian(4,2),jacobian(4,3),jacobian(4,4),jacobian(4,5),jacobian(4,6),jacobian(4,7),&
!  jacobian(4,8),jacobian(4,9),jacobian(4,10),jacobian(4,11),jacobian(4,12),jacobian(4,13)
!write(*,202) jacobian(5,1),jacobian(5,2),jacobian(5,3),jacobian(5,4),jacobian(5,5),jacobian(5,6),jacobian(5,7),&
!  jacobian(5,8),jacobian(5,9),jacobian(5,10),jacobian(5,11),jacobian(5,12),jacobian(5,13)
!write(*,202) jacobian(6,1),jacobian(6,2),jacobian(6,3),jacobian(6,4),jacobian(6,5),jacobian(6,6),jacobian(6,7),&
!  jacobian(6,8),jacobian(6,9),jacobian(6,10),jacobian(6,11),jacobian(6,12),jacobian(6,13)
!write(*,202) jacobian(7,1),jacobian(7,2),jacobian(7,3),jacobian(7,4),jacobian(7,5),jacobian(7,6),jacobian(7,7),&
!  jacobian(7,8),jacobian(7,9),jacobian(7,10),jacobian(7,11),jacobian(7,12),jacobian(7,13)
!write(*,202) jacobian(8,1),jacobian(8,2),jacobian(8,3),jacobian(8,4),jacobian(8,5),jacobian(8,6),jacobian(8,7),&
!  jacobian(8,8),jacobian(8,9),jacobian(8,10),jacobian(8,11),jacobian(8,12),jacobian(8,13)
!write(*,202) jacobian(9,1),jacobian(9,2),jacobian(9,3),jacobian(9,4),jacobian(9,5),jacobian(9,6),jacobian(9,7),&
!  jacobian(9,8),jacobian(9,9),jacobian(9,10),jacobian(9,11),jacobian(9,12),jacobian(9,13)
!write(*,202) jacobian(10,1),jacobian(10,2),jacobian(10,3),jacobian(10,4),jacobian(10,5),jacobian(10,6),jacobian(10,7),&
!  jacobian(10,8),jacobian(10,9),jacobian(10,10),jacobian(10,11),jacobian(10,12),jacobian(10,13)
!write(*,202) jacobian(11,1),jacobian(11,2),jacobian(11,3),jacobian(11,4),jacobian(11,5),jacobian(11,6),jacobian(11,7),&
!  jacobian(11,8),jacobian(11,9),jacobian(11,10),jacobian(11,11),jacobian(11,12),jacobian(11,13)
!write(*,202) jacobian(12,1),jacobian(12,2),jacobian(12,3),jacobian(12,4),jacobian(12,5),jacobian(12,6),jacobian(12,7),&
!  jacobian(12,8),jacobian(12,9),jacobian(12,10),jacobian(12,11),jacobian(12,12),jacobian(12,13)
!write(*,202) jacobian(13,1),jacobian(13,2),jacobian(13,3),jacobian(13,4),jacobian(13,5),jacobian(13,6),jacobian(13,7),&
!  jacobian(13,8),jacobian(13,9),jacobian(13,10),jacobian(13,11),jacobian(13,12),jacobian(13,13)
!
!
!
    CALL lu_decomp(jacobian,em,loc_rho,loc_u,loc_v,loc_w,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
      loc_pyz,loc_q_x,loc_q_y,loc_q_z,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
      loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff)

    end do

  else if (dir == 39) then

  if (.not. shifted) then

  do while (max_diff > 1.0D-7)
!    write(*,*) 'local values',loc_rho,loc_u,loc_v,loc_w,loc_temp,&
!        loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
!        loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
!        loc_lambda_z,self
  loop_counter = loop_counter + 1

!  if (loop_counter >= 15) then
!    write(*,*) loop_counter,max_diff,loc_rho,loc_u,loc_v,&
!      loc_w, loc_temp, loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
!      loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
!      loc_lambda_z
!  end if

  do i = 1,dir
    expo(i) = fieq_comp_39(loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z,i)

!    write(*,*) 'fieq value',expo(i),' for direction ',i
  end do
!  write(*,*) 'dummy/infinite loop check',expo(1),max_diff,self
! Row 1, density derivatives, 1
    jacobian(1,1) = SUM(expo) !drho_dchi
    !drho_dzeta_x
    sp_cx = expo(2)-expo(5)+expo(8)-expo(9)-expo(10) +&
      expo(11)+expo(12)-expo(13)-expo(14)+expo(15)
    jacobian(1,2) = sp_cx+2.0D0*(expo(16)-expo(19)+expo(22)+expo(23)-expo(25)-&
      expo(26)+expo(27)+expo(28)-expo(29)-expo(30))+&
      3.0D0*(expo(34)-expo(37))
    !drho_dzeta_y
    sp_cy =expo(3)-expo(6)+expo(8)+expo(9)-expo(10)-&
      expo(11)+expo(12)+expo(13)-expo(14)-expo(15)
    jacobian(1,3) = sp_cy+2.0D0*(expo(17)-expo(20)+expo(22)+expo(24)+expo(25)-expo(26)-expo(27)+&
      expo(31)-expo(32)-expo(33))+3.0D0*(expo(35)-expo(38))
    !drho_dzeta_z
    sp_cz = expo(4)-expo(7)+expo(8)+expo(9)+expo(10)+&
      expo(11)-expo(12)-expo(13)-expo(14)-expo(15)
    jacobian(1,4) = sp_cz+2.0D0*(expo(18)-expo(21)+expo(23)+expo(24)-expo(28)-&
      expo(29)+expo(30)-expo(31)-expo(32)+expo(33))+&
      3.0D0*(expo(36)-expo(39))
    !drho_dpi_xx
    sp_8dir_pos = SUM(expo(8:15))
    sp_cx2 = expo(2)+expo(5)+sp_8dir_pos
    jacobian(1,5) = sp_cx2 + 4.0D0*(expo(16)+expo(19)+expo(22)+&
      expo(23)+expo(25)+&
      expo(26)+expo(27)+expo(28)+expo(29)+expo(30))+&
      9.0D0*(expo(34)+expo(37))
    !drho_dpi_yy
    sp_cy2 = expo(3)+expo(6)+sp_8dir_pos
    jacobian(1,6) = sp_cy2+4.0D0*(expo(17)+&
      expo(20)+expo(22)+expo(24)+expo(25)+expo(26)&
      +expo(27)+expo(31)+expo(32)+expo(33))+9.0D0*(expo(35)+expo(38))
    !drho_dpi_zz
    sp_cz2 = expo(4)+expo(7)+sp_8dir_pos
    jacobian(1,7) = sp_cz2+4.0D0*(expo(18)+expo(21)+expo(23)+expo(24)+&
      expo(28)+expo(29)+expo(30)+expo(31)+&
      expo(32)+expo(33))+9.0D0*(expo(36)+expo(39))
    !drho_dpi_xy
    sp_cxcy = expo(8)-expo(9)+expo(10)-expo(11) &
      +expo(12)-expo(13)+expo(14)-expo(15)
    jacobian(1,8) = sp_cxcy+4.0D0*(expo(22)-expo(25)+expo(26) &
      -expo(27))
    !drho_dpi_xz
    sp_cxcz = expo(8)-expo(9)-expo(10)+expo(11) &
      -expo(12)+expo(13)+expo(14)-expo(15)
    jacobian(1,9) = sp_cxcz+4.0D0*(expo(23)-expo(28)+expo(29)-expo(30))
    !drho_dpi_yz
    sp_cycz = expo(8)+expo(9)-expo(10)-expo(11) &
      -expo(12)-expo(13)+expo(14)+expo(15)
    jacobian(1,10) = sp_cycz+4.0D0*(expo(24)-expo(31)+expo(32) &
      -expo(33))
    !drho_dlambda_x
    sp_c2cx = expo(2)-expo(5)+3.0D0*(expo(8)-expo(9)-expo(10) +&
      expo(11)+expo(12)-expo(13)-&
      expo(14)+expo(15))
    jacobian(1,11) = sp_c2cx +8.0D0*(expo(16)-expo(19))+16.0D0*(expo(22) &
      +expo(23)-expo(25)-expo(26)+expo(27) &
      +expo(28)-expo(29)-expo(30))+27.0D0*(expo(34)-expo(37)) !cx*c^2
    !drho_dlambda_y
    sp_c2cy = expo(3)-expo(6)+3.0D0*(expo(8)+expo(9)-expo(10) -&
      expo(11)+expo(12)+expo(13)-&
      expo(14)-expo(15))
    jacobian(1,12) = sp_c2cy+8.0D0*(expo(17)-expo(20))+16.0D0*(expo(22) &
      +expo(24)+expo(25)-expo(26)-expo(27) &
      +expo(31)-expo(32)-expo(33)) &
      +27.0D0*(expo(35)-expo(38)) !cy*c^2
    !drho_dlambda_z
    sp_c2cz = expo(4)-expo(7)+3.0D0*(expo(8)+expo(9)+expo(10) +&
      expo(11)-expo(12)-expo(13) &
      -expo(14)-expo(15))
    jacobian(1,13) = sp_c2cz +8.0D0*(expo(18)-expo(21))+16.0D0*(expo(23) &
      +expo(24)-expo(28)-expo(29)+expo(30) &
      -expo(31)-expo(32)+expo(33))+27.0D0*(expo(36) -expo(39)) !cz*c^2
!
! Row 2, x-momentum derivatives, cx
!
    jacobian(2,1) = jacobian(1,2) !cx
    jacobian(2,2) = jacobian(1,5) !cx^2
    jacobian(2,3) = jacobian(1,8) !cxcy
    jacobian(2,4) = jacobian(1,9) !cxcz
    jacobian(2,5) = sp_cx+8.0D0*(expo(16)-expo(19)+expo(22)+expo(23)-&
      expo(25)-expo(26)+expo(27)+expo(28)-expo(29)-&
      expo(30))+27.0D0*(expo(34)-expo(37)) !cx^3
!
    jacobian(2,6) = expo(8)-expo(9)-expo(10)+expo(11)+expo(12)-expo(13)-expo(14)+expo(15)+&
      8.0D0*(expo(22)-expo(25)-expo(26)+expo(27)) !cx*cy^2
!
    jacobian(2,7) = expo(8)-expo(9)-expo(10)+expo(11)+expo(12)-expo(13)-expo(14)+expo(15)+&
      8.0D0*(expo(23)+expo(28)-expo(29)-expo(30))!cx*cz^2
!
    jacobian(2,8) = expo(8)+expo(9)-expo(10)-expo(11)+expo(12)+expo(13)-expo(14)-expo(15)+&
      8.0D0*(expo(22)+expo(25)-expo(26)-expo(27))!cx^2*cy
!
    jacobian(2,9) = expo(8)+expo(9)+expo(10)+expo(11)-expo(12)-expo(13)-expo(14)-expo(15)+&
      8.0D0*(expo(23)-expo(28)-expo(29)+expo(30))!cx^2*cz
!
    jacobian(2,10) = expo(8)-expo(9)+expo(10) -expo(11)-expo(12)+expo(13)&
      -expo(14)+expo(15) !cx*cy*cz
!
    jacobian(2,11) =expo(2)+expo(5)+3.0D0*sp_8dir_pos+&
      16.0D0*(expo(16)+expo(19))+32.0D0*(expo(22)+expo(23) &
      +expo(25)+expo(26)+expo(27)+expo(28)+expo(29)+expo(30)) + &
      81.0D0*(expo(34)+expo(37)) !cx^2*c^2

    jacobian(2,12) = sp_cxcy*3.0D0 +32.0D0*(expo(22)-expo(25)+expo(26)-expo(27)) !cx*cy*c^2

    jacobian(2,13) = sp_cxcz*3.0D0 + 32.0D0*(expo(23)-expo(28)+expo(29)-expo(30)) !cx*cz*c^2
! Row 3, y-momentum derivatives, cy
    jacobian(3,1) = jacobian(1,3) !cy
    jacobian(3,2) = jacobian(1,8) !cx*cy
    jacobian(3,3) = jacobian(1,6) !cy^2
    jacobian(3,4) = jacobian(1,10) !cy*cz
    jacobian(3,5) = jacobian(2,8) !cx^2*cy
    jacobian(3,6) = sp_cy + 8.0D0*(expo(17) - expo(20)+ expo(22)+ expo(24) &
      + expo(25)- expo(26) - expo(27) + expo(31) - expo(32)- expo(33))&
      + 27.0D0*(expo(35)- expo(38)) !cy^3
    jacobian(3,7) = expo(8)+expo(9)-expo(10)-expo(11)+expo(12)+expo(13)-&
      expo(14)-expo(15)+8.0D0*(expo(24)+expo(31)-expo(32)-expo(33)) !cycz^2
    jacobian(3,8) = jacobian(2,6) !cx*cy^2
    jacobian(3,9) = jacobian(2,10) !cx*cy*cz
    jacobian(3,10) = expo(8)+expo(9)+expo(10)+expo(11)-expo(12)-expo(13)-expo(14)-expo(15)+&
      8.0D0*(expo(24)-expo(31)-expo(32)+expo(33)) !cy^2*cz
    jacobian(3,11) = jacobian(2,12) !cx*cy*c^2
    jacobian(3,12) = expo(3)+expo(6)+3.0D0*sp_8dir_pos+16.0D0*(expo(17)+expo(20))+32.0D0*(expo(22)+expo(24)+expo(25)+expo(26)+&
      expo(27)+expo(31)+expo(32)+expo(33))+81.0D0*(expo(35)+expo(38)) !cy^2*c^2
    jacobian(3,13) = 3.0D0*sp_cycz + 32.0D0*(expo(24)-expo(31)+expo(32)-expo(33)) !cy*cz*c^2
! Row 4, z-momentum derivatives, cz
    jacobian(4,1) = jacobian(1,4) !cz
    jacobian(4,2) = jacobian(2,4) !cx*cz
    jacobian(4,3) = jacobian(3,4) !cy*cz
    jacobian(4,4) = jacobian(1,7) !cx^2
    jacobian(4,5) = jacobian(2,9) !cx^2*cz
    jacobian(4,6) = jacobian(3,10) !cy^2*cz
    jacobian(4,7) = sp_cz + 8.0D0*(expo(18)-expo(21)+expo(23)+expo(24)-expo(28)&
      -expo(29)+expo(30)-expo(31)-expo(32)+expo(33))+27.0D0*(expo(36)-expo(39)) !cz^3
    jacobian(4,8) = jacobian(3,9) !cx*cy*cz
    jacobian(4,9) = jacobian(2,7) !cx*cz^2
    jacobian(4,10) = jacobian(3,7) !cy*cz^2
    jacobian(4,11) = jacobian(2,13) !cx*cz*c^2
    jacobian(4,12) = jacobian(3,13) !cy*cz*c^2
    jacobian(4,13) = expo(4)+expo(7) +3.0D0*sp_8dir_pos+16.0D0*(expo(18)+expo(21))+32.0D0*(expo(23)+&
      expo(24)+sum(expo(28:33)))+81.0D0*(expo(36)+expo(39)) !cz^2*c^2
! Row 5, x-directional pressure, cx^2
    jacobian(5,1) = jacobian(1,5) !cx^2
    jacobian(5,2) = jacobian(2,5) !cx^3
    jacobian(5,3) = jacobian(3,5) !cx^2*cy or (2,8)
    jacobian(5,4) = jacobian(4,5) !cx^2*cz
    jacobian(5,5) = expo(2)+expo(5)+sp_8dir_pos+16.0D0*(expo(16)+expo(19)+expo(22)+&
      expo(23)+expo(25)+expo(26)+expo(27)+expo(28)+expo(29)+expo(30))+81.0D0*(expo(34)+expo(37)) !cx^4
    jacobian(5,6) = sp_8dir_pos + 16.0D0*(expo(22)+expo(25)+expo(26)+expo(27)) !cx^2*cy^2
    jacobian(5,7) = sp_8dir_pos + 16.0D0*(expo(23)+expo(28)+expo(29)+expo(30)) !cx^2*cz^2
    jacobian(5,8) = sp_cxcy+16.0D0*(expo(22)-expo(25)+expo(26)-expo(27)) !cx^3*cy
    jacobian(5,9) = sp_cxcz+16.0D0*(expo(23)-expo(28)+expo(29)-expo(30)) ! cx^3*cz
    jacobian(5,10) = expo(8)+expo(9)-expo(10)-expo(11)-expo(12)-expo(13)+expo(14)+expo(15) !cx^2*cy*cz
    jacobian(5,11) = sp_c2cx + 32.0D0*(expo(16)-expo(19))+64.0D0*(expo(22)+expo(23)-expo(25)-expo(26)+&
      expo(27)+expo(28)-expo(29)-expo(30))+243.0D0*(expo(34)-expo(37)) !cx^3*c^2
    jacobian(5,12) = 3.0D0*(expo(8)+expo(9)-expo(10)-expo(11)+expo(12)+expo(13)-expo(14)-expo(15))+&
      64.0D0*(expo(22)+expo(25)-expo(26)-expo(27))!cx^2*cy*c^2
    jacobian(5,13) = 3.0D0*(expo(8)+expo(9)+expo(10)+expo(11)-expo(12)-expo(13)-expo(14)-expo(15))+&
      64.0D0*(expo(23)-expo(28)-expo(29)+expo(30))!cx^2*cz*c^2
! Row 6, y-directional pressure, cy^2
    jacobian(6,1) = jacobian(1,6) !cy^2
    jacobian(6,2) = jacobian(2,6) !cx*cy^2
    jacobian(6,3) = jacobian(3,6) !cy^3
    jacobian(6,4) = jacobian(4,6) !cy^2*cz
    jacobian(6,5) = jacobian(5,6) !cx^2*cy^2
    jacobian(6,6) = expo(3)+expo(6)+sp_8dir_pos+16.0D0*(expo(17)+expo(20)+expo(22)+&
      expo(24)+expo(25)+expo(26)+expo(27)+expo(31)+expo(32)+expo(33))+81.0D0*(expo(35)+expo(38)) !cy^4
    jacobian(6,7) = sp_8dir_pos+16.0D0*(expo(24)+expo(31)+expo(32)+expo(33)) !cy^2*cz^2
    jacobian(6,8) = sp_cxcy+16.0D0*(expo(22)-expo(25)+expo(26)-expo(27) )!cy^3*cx
    jacobian(6,9) = expo(8)-expo(9)-expo(10)+expo(11)-expo(12)+expo(13)+expo(14)-expo(15) !cx*cy^2*cz
    jacobian(6,10) = sp_cycz+16.0D0*(expo(24)-expo(31)+expo(32)-expo(33))!cy^3*cz
    jacobian(6,11) = 3.0D0*(expo(8)-expo(9)-expo(10)+expo(11)+expo(12)-expo(13)-expo(14)+expo(15))+&
      64.0D0*(expo(22)-expo(25)-expo(26)+expo(27)) !cx*cy^2*c^2
    jacobian(6,12) = sp_c2cy + 32.0D0*(expo(17)-expo(20))+64.0D0*(expo(22)+expo(24) &
      + expo(25)- expo(26) - expo(27) + expo(31)-expo(32)-expo(33))&
      + 243.0D0*(expo(35)- expo(38)) !cy^3*c^2
    jacobian(6,13) = 3.0D0*(expo(8)+expo(9)+expo(10)+expo(11)-expo(12)-expo(13)-expo(14)-expo(15))+&
      64.0D0*(expo(24)-expo(31)-expo(32)+expo(33))!cz*cy^2*c^2
! Row 7, z-directional pressure, cz^2
    jacobian(7,1) = jacobian(1,7) !cz^2
    jacobian(7,2) = jacobian(2,7) !cx*cz^2
    jacobian(7,3) = jacobian(3,7) !cy*cz^2
    jacobian(7,4) = jacobian(4,7) !cz^3
    jacobian(7,5) = jacobian(5,7) !cx^2*cz^2
    jacobian(7,6) = jacobian(6,7) !cy^2*cz^2
    jacobian(7,7) = expo(4)+expo(7)+sp_8dir_pos+16.0D0*(expo(18)+expo(21)+expo(23)+expo(24)+expo(28)+expo(29)+&
      expo(30)+expo(31)+expo(32)+expo(33))+81.0D0*(expo(36)+expo(39)) !cz^4
    jacobian(7,8) = expo(8)-expo(9)+expo(10)-expo(11)+expo(12)-expo(13)+expo(14)-expo(15)!cx*cy*cz^2
    jacobian(7,9) = sp_cxcz + 16.0D0*(expo(23)-expo(28)+expo(29)-expo(30)) !cx*cz^3
    jacobian(7,10) = sp_cycz + 16.0D0*(expo(24)-expo(31)+expo(32)-expo(33))!cy*cz^3
    jacobian(7,11) = 3.0D0*(expo(8)-expo(9)-expo(10)+expo(11)+expo(12)-expo(13)-expo(14)+expo(15))+&
      64.0D0*(expo(23)+expo(28)-expo(29)-expo(30))
    !expo(4)+expo(7)+3.0D0*sp_8dir_pos+16.0D0*(expo(18)+expo(21))+32.0D0*&
    !  (expo(18)+expo(21)+expo(23)+expo(24)+expo(28)+expo(29)+&
    !  expo(30)+expo(31)+expo(32)+expo(33))+81.0D0*(expo(36)+expo(39))!cx*cz^2*c^2
    jacobian(7,12) = 3.0D0*(expo(8)+expo(9)-expo(10)-expo(11)+expo(12)+expo(13)-expo(14)-expo(15))+&
      64.0D0*(expo(24)+expo(31)-expo(32)-expo(33))!cy*cz^2*c^2
    jacobian(7,13) = sp_c2cz + 32.0D0*(expo(18)-expo(21))+64.0D0*(expo(23)+&
      expo(24)-expo(28)-expo(29)+expo(30)-expo(31)-expo(32)+expo(33))+243.0D0*(expo(36)-expo(39))!cz^3*c^2
! Row 8, xy-shear, cx*cy
    jacobian(8,1) = jacobian(1,8)
    jacobian(8,2) = jacobian(2,8)
    jacobian(8,3) = jacobian(3,8)
    jacobian(8,4) = jacobian(4,8)
    jacobian(8,5) = jacobian(5,8)
    jacobian(8,6) = jacobian(6,8)
    jacobian(8,7) = jacobian(7,8)
    jacobian(8,8) = jacobian(5,6)
    jacobian(8,9) = jacobian(5,10)
    jacobian(8,10) = jacobian(6,9)
    jacobian(8,11) = jacobian(5,12)
    jacobian(8,12) = jacobian(6,11)
    jacobian(8,13) = 3.0D0*jacobian(2,10)
! Row 9, xz-shear, cx*cz
    jacobian(9,1) = jacobian(1,9)
    jacobian(9,2) = jacobian(2,9)
    jacobian(9,3) = jacobian(8,4)
    jacobian(9,4) = jacobian(4,9)
    jacobian(9,5) = jacobian(5,9)
    jacobian(9,6) = jacobian(6,9)
    jacobian(9,7) = jacobian(7,9)
    jacobian(9,8) = jacobian(8,9)
    jacobian(9,9) = jacobian(5,7)
    jacobian(9,10) = jacobian(8,7)
    jacobian(9,11) = jacobian(5,13)
    jacobian(9,12) = jacobian(8,13)
    jacobian(9,13) = jacobian(7,11)
! Row 10, yz-shear, cy*cz
    jacobian(10,1) = jacobian(1,10)
    jacobian(10,2) = jacobian(2,10)
    jacobian(10,3) = jacobian(3,10)
    jacobian(10,4) = jacobian(4,10)
    jacobian(10,5) = jacobian(5,10)
    jacobian(10,6) = jacobian(6,10)
    jacobian(10,7) = jacobian(7,10)
    jacobian(10,8) = jacobian(8,10)
    jacobian(10,9) = jacobian(9,10)
    jacobian(10,10) = jacobian(6,7)
    jacobian(10,11) = jacobian(9,12)
    jacobian(10,12) = jacobian(6,13)
    jacobian(10,13) = jacobian(7,12)
! Row 11, x-direction contracted heat flux
    jacobian(11,1) = jacobian(1,11)
    jacobian(11,2) = jacobian(2,11)
    jacobian(11,3) = jacobian(3,11)
    jacobian(11,4) = jacobian(4,11)
    jacobian(11,5) = jacobian(5,11)
    jacobian(11,6) = jacobian(6,11)
    jacobian(11,7) = jacobian(7,11)
    jacobian(11,8) = jacobian(8,11)
    jacobian(11,9) = jacobian(9,11)
    jacobian(11,10) = jacobian(10,11)
    jacobian(11,11) = expo(2)+expo(5)+9.0D0*sp_8dir_pos+64.0D0*(expo(16)+expo(19))+256.0D0*(expo(22)+expo(23)+&
      expo(25)+expo(26)+expo(27)+expo(28)+expo(29)+expo(30))+729.0D0*(expo(34)+expo(37)) !cx^2*c^4
    jacobian(11,12) = 9.0D0*sp_cxcy+256.0D0*(expo(22)-expo(25)+expo(26)-expo(27)) !cx*cy*c^4
    jacobian(11,13) = 9.0D0*sp_cxcz+256.0D0*(expo(23)-&
      expo(28)+expo(29)-expo(30)) !cx*cz*c^4
! Row 12, y-direction contracted heat flux
    jacobian(12,1) = jacobian(1,12)
    jacobian(12,2) = jacobian(2,12)
    jacobian(12,3) = jacobian(3,12)
    jacobian(12,4) = jacobian(4,12)
    jacobian(12,5) = jacobian(5,12)
    jacobian(12,6) = jacobian(6,12)
    jacobian(12,7) = jacobian(7,12)
    jacobian(12,8) = jacobian(8,12)
    jacobian(12,9) = jacobian(9,12)
    jacobian(12,10) = jacobian(10,12)
    jacobian(12,11) = jacobian(11,12)
    jacobian(12,12) = expo(3)+expo(6)+9.0D0*sp_8dir_pos+64.0D0*(expo(17)+expo(20))+256.0D0*(expo(22)+expo(24)+&
      expo(25)+expo(26)+expo(27)+expo(31)+expo(32)+expo(33))+729.0D0*(expo(35)+expo(38)) !cy^2*c^4
    jacobian(12,13) = 9.0D0*sp_cycz+256.0D0*(expo(24)-expo(31)+expo(32)-expo(33)) !cy*cz*c^4
! Row 13, z-direction contracted heat flux
    jacobian(13,1) = jacobian(1,13)
    jacobian(13,2) = jacobian(2,13)
    jacobian(13,3) = jacobian(3,13)
    jacobian(13,4) = jacobian(4,13)
    jacobian(13,5) = jacobian(5,13)
    jacobian(13,6) = jacobian(6,13)
    jacobian(13,7) = jacobian(7,13)
    jacobian(13,8) = jacobian(8,13)
    jacobian(13,9) = jacobian(9,13)
    jacobian(13,10) = jacobian(10,13)
    jacobian(13,11) = jacobian(11,13)
    jacobian(13,12) = jacobian(12,13)
    jacobian(13,13) = expo(4)+expo(7)+9.0D0*sp_8dir_pos+64.0D0*(expo(18)+expo(21))+256.0D0*(expo(23)+expo(24)+&
      expo(28)+expo(29)+expo(30)+expo(31)+expo(32)+expo(33))+729.0D0*(expo(36)+expo(39))
!
! Output to check if things are going nuts.  Comment out when no longer needed
!
!if (self == 0) then
!write(*,*) 'local values',loc_rho,loc_rho*loc_u,loc_rho*loc_v,loc_rho*loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,&
!  loc_pxy,loc_pxz,loc_pyz,loc_q_x,loc_q_y,loc_q_z,loc_u,loc_v,loc_w,vel_mag
!write(*,*) 'local LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!    loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff
!write(*,202) jacobian(1,1),jacobian(1,2),jacobian(1,3),jacobian(1,4),jacobian(1,5),jacobian(1,6),jacobian(1,7),&
!  jacobian(1,8),jacobian(1,9),jacobian(1,10),jacobian(1,11),jacobian(1,12),jacobian(1,13)
!write(*,202) jacobian(2,1),jacobian(2,2),jacobian(2,3),jacobian(2,4),jacobian(2,5),jacobian(2,6),jacobian(2,7),&
!  jacobian(2,8),jacobian(2,9),jacobian(2,10),jacobian(2,11),jacobian(2,12),jacobian(2,13)
!write(*,202) jacobian(3,1),jacobian(3,2),jacobian(3,3),jacobian(3,4),jacobian(3,5),jacobian(3,6),jacobian(3,7),&
!  jacobian(3,8),jacobian(3,9),jacobian(3,10),jacobian(3,11),jacobian(3,12),jacobian(3,13)
!write(*,202) jacobian(4,1),jacobian(4,2),jacobian(4,3),jacobian(4,4),jacobian(4,5),jacobian(4,6),jacobian(4,7),&
!  jacobian(4,8),jacobian(4,9),jacobian(4,10),jacobian(4,11),jacobian(4,12),jacobian(4,13)
!write(*,202) jacobian(5,1),jacobian(5,2),jacobian(5,3),jacobian(5,4),jacobian(5,5),jacobian(5,6),jacobian(5,7),&
!  jacobian(5,8),jacobian(5,9),jacobian(5,10),jacobian(5,11),jacobian(5,12),jacobian(5,13)
!write(*,202) jacobian(6,1),jacobian(6,2),jacobian(6,3),jacobian(6,4),jacobian(6,5),jacobian(6,6),jacobian(6,7),&
!  jacobian(6,8),jacobian(6,9),jacobian(6,10),jacobian(6,11),jacobian(6,12),jacobian(6,13)
!write(*,202) jacobian(7,1),jacobian(7,2),jacobian(7,3),jacobian(7,4),jacobian(7,5),jacobian(7,6),jacobian(7,7),&
!  jacobian(7,8),jacobian(7,9),jacobian(7,10),jacobian(7,11),jacobian(7,12),jacobian(7,13)
!write(*,202) jacobian(8,1),jacobian(8,2),jacobian(8,3),jacobian(8,4),jacobian(8,5),jacobian(8,6),jacobian(8,7),&
!  jacobian(8,8),jacobian(8,9),jacobian(8,10),jacobian(8,11),jacobian(8,12),jacobian(8,13)
!write(*,202) jacobian(9,1),jacobian(9,2),jacobian(9,3),jacobian(9,4),jacobian(9,5),jacobian(9,6),jacobian(9,7),&
!  jacobian(9,8),jacobian(9,9),jacobian(9,10),jacobian(9,11),jacobian(9,12),jacobian(9,13)
!write(*,202) jacobian(10,1),jacobian(10,2),jacobian(10,3),jacobian(10,4),jacobian(10,5),jacobian(10,6),jacobian(10,7),&
!  jacobian(10,8),jacobian(10,9),jacobian(10,10),jacobian(10,11),jacobian(10,12),jacobian(10,13)
!write(*,202) jacobian(11,1),jacobian(11,2),jacobian(11,3),jacobian(11,4),jacobian(11,5),jacobian(11,6),jacobian(11,7),&
!  jacobian(11,8),jacobian(11,9),jacobian(11,10),jacobian(11,11),jacobian(11,12),jacobian(11,13)
!write(*,202) jacobian(12,1),jacobian(12,2),jacobian(12,3),jacobian(12,4),jacobian(12,5),jacobian(12,6),jacobian(12,7),&
!  jacobian(12,8),jacobian(12,9),jacobian(12,10),jacobian(12,11),jacobian(12,12),jacobian(12,13)
!write(*,202) jacobian(13,1),jacobian(13,2),jacobian(13,3),jacobian(13,4),jacobian(13,5),jacobian(13,6),jacobian(13,7),&
!  jacobian(13,8),jacobian(13,9),jacobian(13,10),jacobian(13,11),jacobian(13,12),jacobian(13,13)
!
!write(*,*) expo
!end if
!
!
!
  CALL lu_decomp(jacobian,em,loc_rho,loc_u,loc_v,loc_w,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
    loc_pyz,loc_q_x,loc_q_y,loc_q_z,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
    loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff)

  end do

  else
    do while (max_diff > 1.0D-7)
      jacobian = 0.0D0
      loop_counter = loop_counter + 1

      do i = 1,dir
        expo(i) = fieq_comp_39_shifted(loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
        loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
        loc_lambda_z,i)

!    write(*,*) 'fieq value',expo(i),' for direction ',i
      end do

!      if (a == 139 .and. b == 127 .and. c == 120) then
!        write(*,*) 'expo',expo,a,b,c
!        write(*,*) 'input values',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
!    loc_pyz,loc_q_x,loc_q_y,loc_q_z
!        write(*,*) 'local LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!        loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff,self,loop_counter
!      end if
!        orig_rho = sum(expo)
!        orig_u = sum(cx*expo)
!        orig_v = sum(cy*expo)
!        orig_w = sum(cz*expo)

!        write(*,*) 'orig',orig_rho,orig_u,orig_v,orig_w,loop_counter
!      end if

!      if (abs(loc_rho-orig_rho) > 0.5D0 .or. abs(loc_u-orig_u) > 0.8D0 .or. &
!          abs(loc_v-orig_v) > 0.8D0 .or. abs(loc_w-orig_w) > 0.8D0 .and. &
!           loop_counter < 5) then
!write(*,*) 'pre- LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!          loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff
!        loc_chi = chi_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_zeta_x = zeta_x_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_zeta_y = zeta_y_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_zeta_z = zeta_z_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_pi_xx = pi_xx_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_pi_yy = pi_yy_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_pi_zz = pi_zz_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_pi_xy = pi_xy_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_pi_xz = pi_xz_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_pi_yz = pi_yz_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_lambda_x = lambda_x_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_lambda_y = lambda_y_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        loc_lambda_z = lambda_z_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp)
!        !loc_chi = ((loc_rho-stag_density)/(rho_inf-stag_density))*(loc_chi-rest_chi)+rest_chi
!        !loc_zeta_x = (loc_u/fs_u)*(loc_zeta_x-rest_zeta_x) + rest_zeta_x
!        !loc_zeta_y = 0.0D0!(loc_v/fs_v)*(loc_zeta_y-rest_zeta_y) + rest_zeta_y
!        !loc_zeta_z = 0.0D0!(loc_w/fs_w)*(loc_zeta_z-rest_zeta_z) + rest_zeta_z
!        !loc_pi_xx = (loc_pxx-rest_press)/(fs_pxx-rest_press)*(loc_pi_xx - rest_pixx) + rest_pixx
!        !loc_pi_yy = -0.20D0!(loc_pyy-rest_press)/(fs_pyy-rest_press)*(loc_pi_yy - rest_piyy) + rest_piyy
!        !loc_pi_zz = -0.20D0! (loc_pzz-rest_press)/(fs_pzz-rest_press)*(loc_pi_zz - rest_pizz) + rest_pizz
!        !loc_pi_xy = 0.0D0!loc_pxy/fs_pxy*(loc_pi_xy - rest_pixy) + rest_pixy
!        !loc_pi_xz = 0.0D0!loc_pxz/fs_pxz*(loc_pi_xz - rest_pixz) + rest_pixz
!        !loc_pi_yz = 0.0D0!loc_pyz/fs_pyz*(loc_pi_yz - rest_piyz) + rest_piyz
!        !loc_lambda_x = loc_q_x/fs_qx*(loc_lambda_x - rest_lambda_x) + rest_lambda_x
!        !loc_lambda_y = 0.0D0!loc_q_y/fs_qy*(loc_lambda_y - rest_lambda_y) + rest_lambda_y
!        !loc_lambda_z = 0.0D0!loc_q_z/fs_qz*(loc_lambda_z - rest_lambda_z) + rest_lambda_z
!        write(*,*) 'changed- LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!          loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff,self,loop_counter
!!        write(*,*) 'local values',loc_rho,loc_u,loc_v,loc_w,loc_temp,orig_rho,orig_u,orig_v,orig_w,self,loop_counter
!!        write(*,*) 'reference terms',stag_density,rest_chi,rest_zeta_x,rest_pixx,rest_lambda_x
!!        write(*,*) 'freestream values',rho_inf,fs_u,fs_v,fs_w,temperature,fs_pxx,fs_qx
!
!        do i = 1,dir
!          expo(i) = fieq_comp_39_shifted(loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
!          loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
!          loc_lambda_z,i)
!        end do
!      end if

! Row 1, density derivatives, 1
    jacobian(1,1) = SUM(expo) !drho_dchi
    !drho_dzeta_x
    jacobian(1,2) = sum(cx*expo)
    !drho_dzeta_y
    jacobian(1,3) = sum(cy*expo)
    !drho_dzeta_z
    jacobian(1,4) = sum(cz*expo)
    !drho_dpi_xx
    jacobian(1,5) = sum(cx2*expo)
    !drho_dpi_yy
    jacobian(1,6) = sum(cy2*expo)
    !drho_dpi_zz
    jacobian(1,7) = sum(cz2*expo)
    !drho_dpi_xy
    jacobian(1,8) = sum(cx*cy*expo)
    !drho_dpi_xz
    jacobian(1,9) = sum(cx*cz*expo)
    !drho_dpi_yz
    jacobian(1,10) = sum(cy*cz*expo)
    !drho_dlambda_x
    jacobian(1,11) =  sum(cxc2*expo) !cx*c^2
    !drho_dlambda_y
    jacobian(1,12) =  sum(cyc2*expo)!cy*c^2
    !drho_dlambda_z
    jacobian(1,13) =  sum(czc2*expo)!cz*c^2

!
! Row 2, x-momentum derivatives, cx
!
    jacobian(2,1) = jacobian(1,2) !cx
    jacobian(2,2) = jacobian(1,5) !cx^2
    jacobian(2,3) = jacobian(1,8) !cxcy
    jacobian(2,4) = jacobian(1,9) !cxcz
    jacobian(2,5) = sum(cxxx*expo) !cx^3
    jacobian(2,6) = sum(cx*cy2*expo) !cx*cy^2
    jacobian(2,7) = sum(cx*cz2*expo) !cx*cz^2
    jacobian(2,8) = sum(cx2*cy*expo) !cx^2*cy
    jacobian(2,9) = sum(cx2*cz*expo) !cx^2*cz
    jacobian(2,10) = sum(cxyz*expo) !cx*cy*cz
    jacobian(2,11) = sum(cx*cxc2*expo) !cx^2*c^2
    jacobian(2,12) = sum(cxc2*cy*expo) !cx*cy*c^2
    jacobian(2,13) = sum(cxc2*cz*expo) !cx*cz*c^2
! Row 3, y-momentum derivatives, cy
    jacobian(3,1) = jacobian(1,3) !cy
    jacobian(3,2) = jacobian(1,8) !cx*cy
    jacobian(3,3) = jacobian(1,6) !cy^2
    jacobian(3,4) = jacobian(1,10) !cy*cz
    jacobian(3,5) = jacobian(2,8) !cx^2*cy
    jacobian(3,6) = sum(cyyy*expo) !cy^3
    jacobian(3,7) = sum(cycz*cz*expo) !cycz^2
    jacobian(3,8) = jacobian(2,6) !cx*cy^2
    jacobian(3,9) = jacobian(2,10) !cx*cy*cz
    jacobian(3,10) = sum(cycz*cy*expo) !cy^2*cz
    jacobian(3,11) = jacobian(2,12) !cx*cy*c^2
    jacobian(3,12) = sum(cyc2*cy*expo) !cy^2*c^2
    jacobian(3,13) = sum(cyc2*cz*expo)!cy*cz*c^2
! Row 4, z-momentum derivatives, cz
    jacobian(4,1) = jacobian(1,4) !cz
    jacobian(4,2) = jacobian(2,4) !cx*cz
    jacobian(4,3) = jacobian(3,4) !cy*cz
    jacobian(4,4) = jacobian(1,7) !cx^2
    jacobian(4,5) = jacobian(2,9) !cx^2*cz
    jacobian(4,6) = jacobian(3,10) !cy^2*cz
    jacobian(4,7) = sum(czzz*expo) !cz^3
    jacobian(4,8) = jacobian(3,9) !cx*cy*cz
    jacobian(4,9) = jacobian(2,7) !cx*cz^2
    jacobian(4,10) = jacobian(3,7) !cy*cz^2
    jacobian(4,11) = jacobian(2,13) !cx*cz*c^2
    jacobian(4,12) = jacobian(3,13) !cy*cz*c^2
    jacobian(4,13) = sum(czc2*cz*expo) !cz^2*c^2
! Row 5, x-directional pressure, cx^2
    jacobian(5,1) = jacobian(1,5) !cx^2
    jacobian(5,2) = jacobian(2,5) !cx^3
    jacobian(5,3) = jacobian(3,5) !cx^2*cy or (2,8)
    jacobian(5,4) = jacobian(4,5) !cx^2*cz
    jacobian(5,5) = sum(cx2*cx2*expo) !cx^4
    jacobian(5,6) = sum(cx2*cy2*expo) !cx^2*cy^2
    jacobian(5,7) = sum(cx2*cz2*expo) !cx^2*cz^2
    jacobian(5,8) = sum(cx2*cxcy*expo) !cx^3*cy
    jacobian(5,9) = sum(cxcz*cx2*expo) ! cx^3*cz
    jacobian(5,10) = sum(cx2*cycz*expo) !cx^2*cy*cz
    jacobian(5,11) = sum(cx2*cxc2*expo) !cx^3*c^2
    jacobian(5,12) = sum(cx2*cyc2*expo) !cx^2*cy*c^2
    jacobian(5,13) = sum(cx2*czc2*expo) !cx^2*cz*c^2
! Row 6, y-directional pressure, cy^2
    jacobian(6,1) = jacobian(1,6) !cy^2
    jacobian(6,2) = jacobian(2,6) !cx*cy^2
    jacobian(6,3) = jacobian(3,6) !cy^3
    jacobian(6,4) = jacobian(4,6) !cy^2*cz
    jacobian(6,5) = jacobian(5,6) !cx^2*cy^2
    jacobian(6,6) = sum(cy2*cy2*expo) !cy^4
    jacobian(6,7) = sum(cy2*cz2*expo) !cy^2*cz^2
    jacobian(6,8) = sum(cxcy*cy2*expo) !cy^3*cx
    jacobian(6,9) = sum(cxcz*cy2*expo) !cx*cy^2*cz
    jacobian(6,10) = sum(cy2*cycz*expo) !cy^3*cz
    jacobian(6,11) = sum(cxc2*cy2*expo) !cx*cy^2*c^2
    jacobian(6,12) = sum(cy2*cyc2*expo) !cy^3*c^2
    jacobian(6,13) = sum(cy2*czc2*expo) !cz*cy^2*c^2
! Row 7, z-directional pressure, cz^2
    jacobian(7,1) = jacobian(1,7) !cz^2
    jacobian(7,2) = jacobian(2,7) !cx*cz^2
    jacobian(7,3) = jacobian(3,7) !cy*cz^2
    jacobian(7,4) = jacobian(4,7) !cz^3
    jacobian(7,5) = jacobian(5,7) !cx^2*cz^2
    jacobian(7,6) = jacobian(6,7) !cy^2*cz^2
    jacobian(7,7) = sum(cz2*cz2*expo) !cz^4
    jacobian(7,8) = sum(cxcy*cz2*expo) !cx*cy*cz^2
    jacobian(7,9) = sum(cxcz*cz2*expo) !cx*cz^3
    jacobian(7,10) = sum(cycz*cz2*expo) !cy*cz^3
    jacobian(7,11) = sum(cxc2*cz2*expo)
    !expo(4)+expo(7)+3.0D0*sp_8dir_pos+16.0D0*(expo(18)+expo(21))+32.0D0*&
    !  (expo(18)+expo(21)+expo(23)+expo(24)+expo(28)+expo(29)+&
    !  expo(30)+expo(31)+expo(32)+expo(33))+81.0D0*(expo(36)+expo(39))!cx*cz^2*c^2
    jacobian(7,12) = sum(cyc2*cz2*expo) !cy*cz^2*c^2
    jacobian(7,13) = sum(czc2*cz2*expo) !cz^3*c^2
! Row 8, xy-shear, cx*cy
    jacobian(8,1) = jacobian(1,8)
    jacobian(8,2) = jacobian(2,8)
    jacobian(8,3) = jacobian(3,8)
    jacobian(8,4) = jacobian(4,8)
    jacobian(8,5) = jacobian(5,8)
    jacobian(8,6) = jacobian(6,8)
    jacobian(8,7) = jacobian(7,8)
    jacobian(8,8) = jacobian(5,6)
    jacobian(8,9) = jacobian(5,10)
    jacobian(8,10) = jacobian(6,9)
    jacobian(8,11) = jacobian(5,12)
    jacobian(8,12) = jacobian(6,11)
    jacobian(8,13) = sum(cxcy*czc2*expo) !cx*cy*cz*c^2
! Row 9, xz-shear, cx*cz
    jacobian(9,1) = jacobian(1,9)
    jacobian(9,2) = jacobian(2,9)
    jacobian(9,3) = jacobian(8,4)
    jacobian(9,4) = jacobian(4,9)
    jacobian(9,5) = jacobian(5,9)
    jacobian(9,6) = jacobian(6,9)
    jacobian(9,7) = jacobian(7,9)
    jacobian(9,8) = jacobian(8,9)
    jacobian(9,9) = jacobian(5,7)
    jacobian(9,10) = jacobian(8,7)
    jacobian(9,11) = jacobian(5,13)
    jacobian(9,12) = jacobian(8,13)
    jacobian(9,13) = jacobian(7,11)
! Row 10, yz-shear, cy*cz
    jacobian(10,1) = jacobian(1,10)
    jacobian(10,2) = jacobian(2,10)
    jacobian(10,3) = jacobian(3,10)
    jacobian(10,4) = jacobian(4,10)
    jacobian(10,5) = jacobian(5,10)
    jacobian(10,6) = jacobian(6,10)
    jacobian(10,7) = jacobian(7,10)
    jacobian(10,8) = jacobian(8,10)
    jacobian(10,9) = jacobian(9,10)
    jacobian(10,10) = jacobian(6,7)
    jacobian(10,11) = jacobian(9,12)
    jacobian(10,12) = jacobian(6,13)
    jacobian(10,13) = jacobian(7,12)
! Row 11, x-direction contracted heat flux
    jacobian(11,1) = jacobian(1,11)
    jacobian(11,2) = jacobian(2,11)
    jacobian(11,3) = jacobian(3,11)
    jacobian(11,4) = jacobian(4,11)
    jacobian(11,5) = jacobian(5,11)
    jacobian(11,6) = jacobian(6,11)
    jacobian(11,7) = jacobian(7,11)
    jacobian(11,8) = jacobian(8,11)
    jacobian(11,9) = jacobian(9,11)
    jacobian(11,10) = jacobian(10,11)
    jacobian(11,11) = sum(cxc2**2*expo) !cx^2*c^4
    jacobian(11,12) = sum(cxc2*cyc2*expo) !cx*cy*c^4
    jacobian(11,13) = sum(cxc2*czc2*expo) !cx*cz*c^4
! Row 12, y-direction contracted heat flux
    jacobian(12,1) = jacobian(1,12)
    jacobian(12,2) = jacobian(2,12)
    jacobian(12,3) = jacobian(3,12)
    jacobian(12,4) = jacobian(4,12)
    jacobian(12,5) = jacobian(5,12)
    jacobian(12,6) = jacobian(6,12)
    jacobian(12,7) = jacobian(7,12)
    jacobian(12,8) = jacobian(8,12)
    jacobian(12,9) = jacobian(9,12)
    jacobian(12,10) = jacobian(10,12)
    jacobian(12,11) = jacobian(11,12)
    jacobian(12,12) = sum(cyc2**2*expo) !cy^2*c^4
    jacobian(12,13) = sum(cyc2*czc2*expo) !cy*cz*c^4
! Row 13, z-direction contracted heat flux
    jacobian(13,1) = jacobian(1,13)
    jacobian(13,2) = jacobian(2,13)
    jacobian(13,3) = jacobian(3,13)
    jacobian(13,4) = jacobian(4,13)
    jacobian(13,5) = jacobian(5,13)
    jacobian(13,6) = jacobian(6,13)
    jacobian(13,7) = jacobian(7,13)
    jacobian(13,8) = jacobian(8,13)
    jacobian(13,9) = jacobian(9,13)
    jacobian(13,10) = jacobian(10,13)
    jacobian(13,11) = jacobian(11,13)
    jacobian(13,12) = jacobian(12,13)
    jacobian(13,13) = sum(czc2**2*expo)
!
!
!
!      if (loop_counter == 1) then


!        CALL lu_decomp(jacobian,em,loc_rho,loc_u,loc_v,loc_w,loc_pxx,loc_pyy,loc_pzz,&
!          loc_pxy,loc_pxz,loc_pyz,loc_q_x,loc_q_y,loc_q_z,loc_chi,loc_zeta_x,&
!          loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!          loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff)

!if (self == 0) then
!  write(*,*) 'passed values',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
!    loc_pyz,loc_q_x,loc_q_y,loc_q_z,self,loop_counter
!  write(*,*) 'local LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!    loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff,self,loop_counter

!write(*,202) jacobian(1,1),jacobian(1,2),jacobian(1,3),jacobian(1,4),jacobian(1,5),jacobian(1,6),jacobian(1,7),&
!  jacobian(1,8),jacobian(1,9),jacobian(1,10),jacobian(1,11),jacobian(1,12),jacobian(1,13)
!write(*,202) jacobian(2,1),jacobian(2,2),jacobian(2,3),jacobian(2,4),jacobian(2,5),jacobian(2,6),jacobian(2,7),&
!  jacobian(2,8),jacobian(2,9),jacobian(2,10),jacobian(2,11),jacobian(2,12),jacobian(2,13)
!write(*,202) jacobian(3,1),jacobian(3,2),jacobian(3,3),jacobian(3,4),jacobian(3,5),jacobian(3,6),jacobian(3,7),&
!  jacobian(3,8),jacobian(3,9),jacobian(3,10),jacobian(3,11),jacobian(3,12),jacobian(3,13)
!write(*,202) jacobian(4,1),jacobian(4,2),jacobian(4,3),jacobian(4,4),jacobian(4,5),jacobian(4,6),jacobian(4,7),&
!  jacobian(4,8),jacobian(4,9),jacobian(4,10),jacobian(4,11),jacobian(4,12),jacobian(4,13)
!write(*,202) jacobian(5,1),jacobian(5,2),jacobian(5,3),jacobian(5,4),jacobian(5,5),jacobian(5,6),jacobian(5,7),&
!  jacobian(5,8),jacobian(5,9),jacobian(5,10),jacobian(5,11),jacobian(5,12),jacobian(5,13)
!write(*,202) jacobian(6,1),jacobian(6,2),jacobian(6,3),jacobian(6,4),jacobian(6,5),jacobian(6,6),jacobian(6,7),&
!  jacobian(6,8),jacobian(6,9),jacobian(6,10),jacobian(6,11),jacobian(6,12),jacobian(6,13)
!write(*,202) jacobian(7,1),jacobian(7,2),jacobian(7,3),jacobian(7,4),jacobian(7,5),jacobian(7,6),jacobian(7,7),&
!  jacobian(7,8),jacobian(7,9),jacobian(7,10),jacobian(7,11),jacobian(7,12),jacobian(7,13)
!write(*,202) jacobian(8,1),jacobian(8,2),jacobian(8,3),jacobian(8,4),jacobian(8,5),jacobian(8,6),jacobian(8,7),&
!  jacobian(8,8),jacobian(8,9),jacobian(8,10),jacobian(8,11),jacobian(8,12),jacobian(8,13)
!write(*,202) jacobian(9,1),jacobian(9,2),jacobian(9,3),jacobian(9,4),jacobian(9,5),jacobian(9,6),jacobian(9,7),&
!  jacobian(9,8),jacobian(9,9),jacobian(9,10),jacobian(9,11),jacobian(9,12),jacobian(9,13)
!write(*,202) jacobian(10,1),jacobian(10,2),jacobian(10,3),jacobian(10,4),jacobian(10,5),jacobian(10,6),jacobian(10,7),&
!  jacobian(10,8),jacobian(10,9),jacobian(10,10),jacobian(10,11),jacobian(10,12),jacobian(10,13)
!write(*,202) jacobian(11,1),jacobian(11,2),jacobian(11,3),jacobian(11,4),jacobian(11,5),jacobian(11,6),jacobian(11,7),&
!  jacobian(11,8),jacobian(11,9),jacobian(11,10),jacobian(11,11),jacobian(11,12),jacobian(11,13)
!write(*,202) jacobian(12,1),jacobian(12,2),jacobian(12,3),jacobian(12,4),jacobian(12,5),jacobian(12,6),jacobian(12,7),&
!  jacobian(12,8),jacobian(12,9),jacobian(12,10),jacobian(12,11),jacobian(12,12),jacobian(12,13)
!write(*,202) jacobian(13,1),jacobian(13,2),jacobian(13,3),jacobian(13,4),jacobian(13,5),jacobian(13,6),jacobian(13,7),&
!  jacobian(13,8),jacobian(13,9),jacobian(13,10),jacobian(13,11),jacobian(13,12),jacobian(13,13)

!end if

!        write(*,*) 'sparkle blast!!!'
!        max_diff = 1.0D0

 !     else
        CALL lu_decomp(jacobian,em,loc_rho,loc_u,loc_v,loc_w,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
          loc_pyz,loc_q_x,loc_q_y,loc_q_z,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
          loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff)


!write(*,202) jacobian(1,1),jacobian(1,2),jacobian(1,3),jacobian(1,4),jacobian(1,5),jacobian(1,6),jacobian(1,7),&
!  jacobian(1,8),jacobian(1,9),jacobian(1,10),jacobian(1,11),jacobian(1,12),jacobian(1,13)
!write(*,202) jacobian(2,1),jacobian(2,2),jacobian(2,3),jacobian(2,4),jacobian(2,5),jacobian(2,6),jacobian(2,7),&
!  jacobian(2,8),jacobian(2,9),jacobian(2,10),jacobian(2,11),jacobian(2,12),jacobian(2,13)
!write(*,202) jacobian(3,1),jacobian(3,2),jacobian(3,3),jacobian(3,4),jacobian(3,5),jacobian(3,6),jacobian(3,7),&
!  jacobian(3,8),jacobian(3,9),jacobian(3,10),jacobian(3,11),jacobian(3,12),jacobian(3,13)
!write(*,202) jacobian(4,1),jacobian(4,2),jacobian(4,3),jacobian(4,4),jacobian(4,5),jacobian(4,6),jacobian(4,7),&
!  jacobian(4,8),jacobian(4,9),jacobian(4,10),jacobian(4,11),jacobian(4,12),jacobian(4,13)
!write(*,202) jacobian(5,1),jacobian(5,2),jacobian(5,3),jacobian(5,4),jacobian(5,5),jacobian(5,6),jacobian(5,7),&
!  jacobian(5,8),jacobian(5,9),jacobian(5,10),jacobian(5,11),jacobian(5,12),jacobian(5,13)
!write(*,202) jacobian(6,1),jacobian(6,2),jacobian(6,3),jacobian(6,4),jacobian(6,5),jacobian(6,6),jacobian(6,7),&
!  jacobian(6,8),jacobian(6,9),jacobian(6,10),jacobian(6,11),jacobian(6,12),jacobian(6,13)
!write(*,202) jacobian(7,1),jacobian(7,2),jacobian(7,3),jacobian(7,4),jacobian(7,5),jacobian(7,6),jacobian(7,7),&
!  jacobian(7,8),jacobian(7,9),jacobian(7,10),jacobian(7,11),jacobian(7,12),jacobian(7,13)
!write(*,202) jacobian(8,1),jacobian(8,2),jacobian(8,3),jacobian(8,4),jacobian(8,5),jacobian(8,6),jacobian(8,7),&
!  jacobian(8,8),jacobian(8,9),jacobian(8,10),jacobian(8,11),jacobian(8,12),jacobian(8,13)
!write(*,202) jacobian(9,1),jacobian(9,2),jacobian(9,3),jacobian(9,4),jacobian(9,5),jacobian(9,6),jacobian(9,7),&
!  jacobian(9,8),jacobian(9,9),jacobian(9,10),jacobian(9,11),jacobian(9,12),jacobian(9,13)
!write(*,202) jacobian(10,1),jacobian(10,2),jacobian(10,3),jacobian(10,4),jacobian(10,5),jacobian(10,6),jacobian(10,7),&
!  jacobian(10,8),jacobian(10,9),jacobian(10,10),jacobian(10,11),jacobian(10,12),jacobian(10,13)
!write(*,202) jacobian(11,1),jacobian(11,2),jacobian(11,3),jacobian(11,4),jacobian(11,5),jacobian(11,6),jacobian(11,7),&
!  jacobian(11,8),jacobian(11,9),jacobian(11,10),jacobian(11,11),jacobian(11,12),jacobian(11,13)
!write(*,202) jacobian(12,1),jacobian(12,2),jacobian(12,3),jacobian(12,4),jacobian(12,5),jacobian(12,6),jacobian(12,7),&
!  jacobian(12,8),jacobian(12,9),jacobian(12,10),jacobian(12,11),jacobian(12,12),jacobian(12,13)
!write(*,202) jacobian(13,1),jacobian(13,2),jacobian(13,3),jacobian(13,4),jacobian(13,5),jacobian(13,6),jacobian(13,7),&
!  jacobian(13,8),jacobian(13,9),jacobian(13,10),jacobian(13,11),jacobian(13,12),jacobian(13,13)


!      end if
!write(*,*) 'post-LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
!          loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff,self,loop_counter
!      loop_counter = loop_counter + 1


      if (isnan(loc_chi) .or. isnan(loc_zeta_x) .or. isnan(loc_zeta_y) .or. isnan(loc_zeta_z) .or. &
        isnan(loc_pi_xx) .or. isnan(loc_pi_yy) .or. isnan(loc_pi_zz) .or. isnan(loc_pi_xy) .or. &
        isnan(loc_pi_xz) .or. isnan(loc_pi_yz) .or. isnan(loc_lambda_x) .or. isnan(loc_lambda_x) .or. &
        isnan(loc_lambda_y) .or. isnan(loc_lambda_z)) then
        write(*,*) 'post-LMs',loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
          loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff,self,loop_counter
        write(*,*) 'incoming values',loc_rho,loc_u,loc_v,loc_w,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
          loc_pyz,loc_q_x,loc_q_y,loc_q_z,a,b,c
write(*,202) jacobian(1,1),jacobian(1,2),jacobian(1,3),jacobian(1,4),jacobian(1,5),jacobian(1,6),jacobian(1,7),&
  jacobian(1,8),jacobian(1,9),jacobian(1,10),jacobian(1,11),jacobian(1,12),jacobian(1,13)
write(*,202) jacobian(2,1),jacobian(2,2),jacobian(2,3),jacobian(2,4),jacobian(2,5),jacobian(2,6),jacobian(2,7),&
  jacobian(2,8),jacobian(2,9),jacobian(2,10),jacobian(2,11),jacobian(2,12),jacobian(2,13)
write(*,202) jacobian(3,1),jacobian(3,2),jacobian(3,3),jacobian(3,4),jacobian(3,5),jacobian(3,6),jacobian(3,7),&
  jacobian(3,8),jacobian(3,9),jacobian(3,10),jacobian(3,11),jacobian(3,12),jacobian(3,13)
write(*,202) jacobian(4,1),jacobian(4,2),jacobian(4,3),jacobian(4,4),jacobian(4,5),jacobian(4,6),jacobian(4,7),&
  jacobian(4,8),jacobian(4,9),jacobian(4,10),jacobian(4,11),jacobian(4,12),jacobian(4,13)
write(*,202) jacobian(5,1),jacobian(5,2),jacobian(5,3),jacobian(5,4),jacobian(5,5),jacobian(5,6),jacobian(5,7),&
  jacobian(5,8),jacobian(5,9),jacobian(5,10),jacobian(5,11),jacobian(5,12),jacobian(5,13)
write(*,202) jacobian(6,1),jacobian(6,2),jacobian(6,3),jacobian(6,4),jacobian(6,5),jacobian(6,6),jacobian(6,7),&
  jacobian(6,8),jacobian(6,9),jacobian(6,10),jacobian(6,11),jacobian(6,12),jacobian(6,13)
write(*,202) jacobian(7,1),jacobian(7,2),jacobian(7,3),jacobian(7,4),jacobian(7,5),jacobian(7,6),jacobian(7,7),&
  jacobian(7,8),jacobian(7,9),jacobian(7,10),jacobian(7,11),jacobian(7,12),jacobian(7,13)
write(*,202) jacobian(8,1),jacobian(8,2),jacobian(8,3),jacobian(8,4),jacobian(8,5),jacobian(8,6),jacobian(8,7),&
  jacobian(8,8),jacobian(8,9),jacobian(8,10),jacobian(8,11),jacobian(8,12),jacobian(8,13)
write(*,202) jacobian(9,1),jacobian(9,2),jacobian(9,3),jacobian(9,4),jacobian(9,5),jacobian(9,6),jacobian(9,7),&
  jacobian(9,8),jacobian(9,9),jacobian(9,10),jacobian(9,11),jacobian(9,12),jacobian(9,13)
write(*,202) jacobian(10,1),jacobian(10,2),jacobian(10,3),jacobian(10,4),jacobian(10,5),jacobian(10,6),jacobian(10,7),&
  jacobian(10,8),jacobian(10,9),jacobian(10,10),jacobian(10,11),jacobian(10,12),jacobian(10,13)
write(*,202) jacobian(11,1),jacobian(11,2),jacobian(11,3),jacobian(11,4),jacobian(11,5),jacobian(11,6),jacobian(11,7),&
  jacobian(11,8),jacobian(11,9),jacobian(11,10),jacobian(11,11),jacobian(11,12),jacobian(11,13)
write(*,202) jacobian(12,1),jacobian(12,2),jacobian(12,3),jacobian(12,4),jacobian(12,5),jacobian(12,6),jacobian(12,7),&
  jacobian(12,8),jacobian(12,9),jacobian(12,10),jacobian(12,11),jacobian(12,12),jacobian(12,13)
write(*,202) jacobian(13,1),jacobian(13,2),jacobian(13,3),jacobian(13,4),jacobian(13,5),jacobian(13,6),jacobian(13,7),&
  jacobian(13,8),jacobian(13,9),jacobian(13,10),jacobian(13,11),jacobian(13,12),jacobian(13,13)
      end if

    end do

  end if

end if
  loc_fieq = expo
!  write(*,*) '39 direction fieq output',loc_fieq

  deallocate(expo)
  DEALLOCATE(jacobian)
END IF




202  FORMAT (13F15.9)

END SUBROUTINE









!    write(*,*) 'local values',loc_rho,loc_u,loc_v,loc_w,loc_temp,&
!        loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
!        loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
!        loc_lambda_z,self

!  if (loop_counter >= 15) then
!    write(*,*) loop_counter,max_diff,loc_rho,loc_u,loc_v,&
!      loc_w, loc_temp, loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
!      loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
!      loc_lambda_z
!  end if

!        jacobian = 0.0D0

        !etr = REAL(dimensions,kind=dp)/2.0D0 * loc_temp
!        DO i = 1,dir

!          cmag = (REAL(cx(i),kind=dp)**2+REAL(cy(i),kind=dp)**2 + &
!            REAL(cz(i),kind=dp)**2)
! Rest velocity 0,0,0
!  expo(1) = loc_rho*EXP(loc_chi)
!! First shell 1,0,0
!  expo(2) = loc_rho*EXP(loc_chi + loc_zeta_x + loc_pi_xx + loc_lambda_x)
!  expo(3) = loc_rho*EXP(loc_chi + loc_zeta_y + loc_pi_yy + loc_lambda_y)
!  expo(4) = loc_rho*EXP(loc_chi + loc_zeta_z + loc_pi_zz + loc_lambda_z)
!  expo(5) = loc_rho*EXP(loc_chi - loc_zeta_x + loc_pi_xx - loc_lambda_x)
!  expo(6) = loc_rho*EXP(loc_chi - loc_zeta_y + loc_pi_yy - loc_lambda_y)
!  expo(7) = loc_rho*EXP(loc_chi - loc_zeta_z + loc_pi_zz - loc_lambda_z)
!! Second shell 1,1,1
!  expo(8) = loc_rho*EXP(loc_chi + loc_zeta_x + loc_zeta_y + loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy + loc_pi_xz + loc_pi_yz +&
!    3.0D0*loc_lambda_x + 3.0D0*loc_lambda_y + 3.0D0*loc_lambda_z)
!  expo(9) = loc_rho*EXP(loc_chi - loc_zeta_x + loc_zeta_y + loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy - loc_pi_xz + loc_pi_yz +&
!    3.0D0*(-loc_lambda_x + loc_lambda_y + loc_lambda_z))
!  expo(10) = loc_rho*EXP(loc_chi - loc_zeta_x - loc_zeta_y + loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy - loc_pi_xz - loc_pi_yz -&
!    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y + 3.0D0*loc_lambda_z)
!  expo(11) = loc_rho*EXP(loc_chi + loc_zeta_x - loc_zeta_y + loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy + loc_pi_xz - loc_pi_yz +&
!    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y + 3.0D0*loc_lambda_z)
!  expo(12) = loc_rho*EXP(loc_chi + loc_zeta_x + loc_zeta_y - loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy - loc_pi_xz - loc_pi_yz + &
!    3.0D0*loc_lambda_x + 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
!  expo(13) = loc_rho*EXP(loc_chi - loc_zeta_x + loc_zeta_y - loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy + loc_pi_xz - loc_pi_yz - &
!    3.0D0*loc_lambda_x + 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
!  expo(14) = loc_rho*EXP(loc_chi - loc_zeta_x - loc_zeta_y - loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy + loc_pi_xz + loc_pi_yz -&
!    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
!  expo(15) = loc_rho*EXP(loc_chi + loc_zeta_x - loc_zeta_y - loc_zeta_z + &
!    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy - loc_pi_xz + loc_pi_yz +&
!    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
!! Third shell 2,0,0
!  expo(16) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x + 4.0D0*loc_pi_xx + 8.0D0*loc_lambda_x)
!  expo(17) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_y + 4.0D0*loc_pi_yy + 8.0D0*loc_lambda_y)
!  expo(18) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_z + 4.0D0*loc_pi_zz + 8.0D0*loc_lambda_z)
!  expo(19) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x + 4.0D0*loc_pi_xx - 8.0D0*loc_lambda_x)
!  expo(20) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_y + 4.0D0*loc_pi_yy - 8.0D0*loc_lambda_y)
!  expo(21) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_z + 4.0D0*loc_pi_zz - 8.0D0*loc_lambda_z)
!! Fourth shell 2,2,0
!  expo(22) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_y + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy + 4.0D0*loc_pi_xy + &
!    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_y)
!  ! 2,0,2
!  expo(23) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_z + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_xz +&
!    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_z)
!  ! 0,2,2
!  expo(24) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_y + 2.0D0*loc_zeta_z + &
!    4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_yz +&
!    16.0D0*loc_lambda_y + 16.0D0*loc_lambda_z)
!  ! -2,2,0
!  expo(25) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_y + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy - 4.0D0*loc_pi_xy - &
!    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_y)
!  ! -2,-2,0
!  expo(26) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_y + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy + 4.0D0*loc_pi_xy  - &
!    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_y)
!  ! 2,-2,0
!  expo(27) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_y + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy - 4.0D0*loc_pi_xy +&
!    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_y)
!  ! 2,0,-2
!  expo(28) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_z + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_xz +&
!    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_z)
!  ! -2,0,-2
!  expo(29) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_z + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_xz - &
!    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_z)
!  ! -2,0,2
!  expo(30) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_z + &
!    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_xz - &
!    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_z)
!  ! 0,2,-2
!  expo(31) = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_y - 2.0D0*loc_zeta_z + &
!     4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_yz + &
!     16.0D0*loc_lambda_y - 16.0D0*loc_lambda_z)
!  ! 0,-2,-2
!  expo(32) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_y - 2.0D0*loc_zeta_z + &
!    4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_yz - &
!    16.0D0*loc_lambda_y - 16.0D0*loc_lambda_z)
!  ! 0,-2,2
!  expo(33) = loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_y + 2.0D0*loc_zeta_z + &
!     4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_yz - &
!    16.0D0*loc_lambda_y + 16.0D0*loc_lambda_z)
!! Fifth shell
!  expo(34) = loc_rho*EXP(loc_chi + 3.0D0*loc_zeta_x + 9.0D0*loc_pi_xx + 27.0D0*loc_lambda_x)
!  expo(35) = loc_rho*EXP(loc_chi + 3.0D0*loc_zeta_y + 9.0D0*loc_pi_yy + 27.0D0*loc_lambda_y)
!  expo(36) = loc_rho*EXP(loc_chi + 3.0D0*loc_zeta_z + 9.0D0*loc_pi_zz + 27.0D0*loc_lambda_z)
!  expo(37) = loc_rho*EXP(loc_chi - 3.0D0*loc_zeta_x + 9.0D0*loc_pi_xx - 27.0D0*loc_lambda_x)
!  expo(38) = loc_rho*EXP(loc_chi - 3.0D0*loc_zeta_y + 9.0D0*loc_pi_yy - 27.0D0*loc_lambda_y)
!  expo(39) = loc_rho*EXP(loc_chi - 3.0D0*loc_zeta_z + 9.0D0*loc_pi_zz - 27.0D0*loc_lambda_z)
!          expo(i) = -loc_rho*EXP(loc_chi + rcx(i)*loc_zeta_x &
!            +rcy(i)*loc_zeta_y + rcz(i)*loc_zeta_z + rcx(i)**2*loc_pi_xx +&
!            +rcy(i)**2*loc_pi_yy + rcz(i)**2*loc_pi_zz + rcx(i)*rcy(i)*loc_pi_xy+ &
!            rcx(i)*rcz(i)*loc_pi_xz + rcy(i)*rcz(i)*loc_pi_yz + &
!            rcx(i)*cmag(i)*loc_lambda_x +rcy(i)*cmag(i)*loc_lambda_y +&
!            rcx(i)*cmag(i)*loc_lambda_z)
!        END DO


!
! Shells
! rest = 1
! 1-dir,speed-1 = 2-7
! 3-dir,speed-1 = 8-15
! 1-dir,speed-2 = 16-21
! 2-dir,speed-2 = 22-33
! 1-dir,speed-3 = 34-39
!
! x = 2,5,8,9,10,11,12,13,14,15,16,19,22,23,25,26,27,28,29,30,34,37
! y = 3,6,8,9,10,11,12,13,14,15,17,20,22,24,25,26,27,31,32,33,35,38
! z = 4,7,8,9,10,11,12,13,14,15,18,21,23,24,28,29,30,31,32,33,36,39
!
! xy = 8,9,10,11,12,13,14,15,22,25,26,27
! xz = 8,9,10,11,12,13,14,15,23,28,29,30
! yz = 8,9,10,11,12,13,14,15,24,31,32,33
!
! xyz = 8,9,10,11,12,13,14,15
!
!

!          jacobian(1,1) = jacobian(1,1) + expo
!          jacobian(2,1) = jacobian(2,1) + REAL(cx(i),kind=dp)*expo
!          jacobian(3,1) = jacobian(3,1) + REAL(cy(i),kind=dp)*expo
!          jacobian(4,1) = jacobian(4,1) + REAL(cz(i),kind=dp)*expo
!          jacobian(5,1) = jacobian(5,1) + cmag(i)/2.0D0*expo
!
!          jacobian(1,2) = jacobian(1,2) + REAL(cx(i),kind=dp)*expo
!          jacobian(2,2) = jacobian(2,2) + REAL(cx(i),kind=dp)**2*expo
!          jacobian(3,2) = jacobian(3,2) + REAL(cx(i),kind=dp)*REAL(cy(i),kind=dp)*expo
!          jacobian(4,2) = jacobian(4,2) + REAL(cx(i),kind=dp)*REAL(cz(i),kind=dp)*expo
!          jacobian(5,2) = jacobian(5,2) + REAL(cx(i),kind=dp)*cmag(i)/2.0D0*expo
!
!          jacobian(1,3) = jacobian(1,3) + REAL(cy(i),kind=dp) * expo
!          jacobian(2,3) = jacobian(2,3) + REAL(cy(i),kind=dp) * REAL(cx(i),kind=dp)*expo
!          jacobian(3,3) = jacobian(3,3) + REAL(cy(i),kind=dp)**2 * expo
!          jacobian(4,3) = jacobian(4,3) + REAL(cy(i),kind=dp) * REAL(cz(i),kind=dp)*expo
!          jacobian(5,3) = jacobian(5,3) + REAL(cy(i),kind=dp) * cmag(i)/2.0D0 * expo
!
!          jacobian(1,4) = jacobian(1,4) + REAL(cz(i),kind=dp) * expo
!          jacobian(2,4) = jacobian(2,4) + REAL(cz(i),kind=dp)*REAL(cx(i),kind=dp)*expo
!          jacobian(3,4) = jacobian(3,4) + REAL(cz(i),kind=dp)*REAL(cy(i),kind=dp)*expo
!          jacobian(4,4) = jacobian(4,4) + REAL(cz(i),kind=dp)**2 * expo
!          jacobian(5,4) = jacobian(5,4) + REAL(cz(i),kind=dp)*cmag(i)/2.0D0*expo
!
!          jacobian(1,5) = jacobian(1,5) + cmag(i)*expo
!          jacobian(2,5) = jacobian(2,5) + REAL(cx(i),kind=dp)*cmag(i)*expo
!          jacobian(3,5) = jacobian(3,5) + REAL(cy(i),kind=dp)*cmag(i)*expo
!          jacobian(4,5) = jacobian(4,5) + REAL(cz(i),kind=dp)*cmag(i)*expo
!          jacobian(5,5) = jacobian(5,5) + cmag(i)**2*expo/2.0D0




!         end do


!loc_pxx = loc_rho*(loc_temp+loc_u**2)
!loc_pyy = loc_rho*(loc_temp+loc_v**2)
!loc_pzz = loc_rho*(loc_temp+loc_u**2)
!loc_pxx = loc_rho*loc_u*loc_v
!loc_pxx = loc_rho*loc_u*loc_w
!loc_pxx = loc_rho*loc_v*loc_w

!write(*,*) 'values entering LM calcs',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,&
!  loc_pzz,loc_pxy,loc_pxz,loc_pyz,&
!  loc_q_x,loc_q_y,loc_q_z


        !etr = REAL(dimensions,kind=dp)/2.0D0 * loc_temp
!        DO i = 1,dir

!          cmag = (REAL(cx(i),kind=dp)**2+REAL(cy(i),kind=dp)**2 )

!          expo(i) = -loc_rho*wi(i)*EXP(loc_chi + REAL(cx(i),kind=dp)*loc_zeta_x &
!            +REAL(cy(i),kind=dp)*loc_zeta_y + cmag(i)*loc_lambda)
!        END DO


!          jacobian(1,1) = jacobian(1,1) + expo
!          jacobian(2,1) = jacobian(2,1) + REAL(cx(i),kind=dp)*expo
!          jacobian(3,1) = jacobian(3,1) + REAL(cy(i),kind=dp)*expo
!          jacobian(4,1) = jacobian(5,1) + cmag(i)/2.0D0*expo
!
!          jacobian(1,2) = jacobian(1,2) + REAL(cx(i),kind=dp)*expo
!          jacobian(2,2) = jacobian(2,2) + REAL(cx(i),kind=dp)**2*expo
!          jacobian(3,2) = jacobian(3,2) + REAL(cx(i),kind=dp)*REAL(cy(i),kind=dp)*expo
!          jacobian(5,2) = jacobian(4,2) + REAL(cx(i),kind=dp)*cmag(i)/2.0D0*expo
!
!          jacobian(1,3) = jacobian(1,3) + REAL(cy(i),kind=dp) * expo
!          jacobian(2,3) = jacobian(2,3) + REAL(cy(i),kind=dp) * REAL(cx(i),kind=dp)*expo
!          jacobian(3,3) = jacobian(3,3) + REAL(cy(i),kind=dp)**2 * expo
!          jacobian(4,3) = jacobian(4,3) + REAL(cy(i),kind=dp)*cmag(i)/2.0D0*expo
!
!          jacobian(1,4) = jacobian(1,4) + cmag(i)*expo
!          jacobian(2,4) = jacobian(2,4) + REAL(cx(i),kind=dp)*cmag(i)*expo
!          jacobian(3,4) = jacobian(3,4) + REAL(cy(i),kind=dp)*cmag(i)*expo
!          jacobian(4,4) = jacobian(4,4) + cmag(i)**2*expo/2.0D0


!        CALL gauss_elim_and_newton_raph(jacobian,deter,cramer_1,cramer_2,cramer_3,cramer_4,cramer_5)
