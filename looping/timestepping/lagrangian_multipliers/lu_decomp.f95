SUBROUTINE lu_decomp(jacobian,em,loc_rho,loc_u,loc_v,loc_w,&
  loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,loc_q_x,loc_q_y,loc_q_z,&
  loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
  loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z,max_diff)
!
! Use an LU decomposition strategey and cramers rule (working right to left) to
! find the values for chi, zeta, and lambda
!
! Called by: calculate_lagrangian_mults
! Calls:
!
USE precise
USE constants

IMPLICIT NONE
REAL(kind=dp),intent(in) :: loc_rho,loc_u,loc_v,&
  loc_w,loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,&
    loc_pyz,loc_q_x,loc_q_y,loc_q_z
REAL(kind=dp) :: loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,loc_pi_yy,loc_pi_zz,&
   loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,loc_lambda_z
REAL(kind=dp) :: max_diff
REAL(kind=dp) :: delta_chi,delta_zeta_x,delta_zeta_y,delta_zeta_z,delta_pixx,delta_piyy,delta_pizza,&
  delta_pixy,delta_pixz,delta_piyz,delta_lambda_x,delta_lambda_y,delta_lambda_z
real(kind=dp) :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13

INTEGER :: em,i,a

REAL(kind=dp) :: jacobian(em,em),lower(em,em),upper(em,em)



lower(1,1) = 1.0D0
lower(2,2) = 1.0D0
lower(3,3) = 1.0D0
lower(4,4) = 1.0D0
lower(5,5) = 1.0D0
lower(6,6) = 1.0D0
lower(7,7) = 1.0D0
lower(8,8) = 1.0D0
lower(9,9) = 1.0D0
lower(10,10) = 1.0D0
lower(11,11) = 1.0D0
lower(12,12) = 1.0D0
lower(13,13) = 1.0D0

upper(1,1) = jacobian(1,1)
upper(1,2) = jacobian(1,2)
upper(1,3) = jacobian(1,3)
upper(1,4) = jacobian(1,4)
upper(1,5) = jacobian(1,5)
upper(1,6) = jacobian(1,6)
upper(1,7) = jacobian(1,7)
upper(1,8) = jacobian(1,8)
upper(1,9) = jacobian(1,9)
upper(1,10) = jacobian(1,10)
upper(1,11) = jacobian(1,11)
upper(1,12) = jacobian(1,12)
upper(1,13) = jacobian(1,13)

lower(2,1) = jacobian(2,1)/upper(1,1)
upper(2,2) = jacobian(2,2)-lower(2,1)*upper(1,2)
upper(2,3) = jacobian(2,3)-lower(2,1)*upper(1,3)
upper(2,4) = jacobian(2,4)-lower(2,1)*upper(1,4)
upper(2,5) = jacobian(2,5)-lower(2,1)*upper(1,5)
upper(2,6) = jacobian(2,6)-lower(2,1)*upper(1,6)
upper(2,7) = jacobian(2,7)-lower(2,1)*upper(1,7)
upper(2,8) = jacobian(2,8)-lower(2,1)*upper(1,8)
upper(2,9) = jacobian(2,9)-lower(2,1)*upper(1,9)
upper(2,10) = jacobian(2,10)-lower(2,1)*upper(1,10)
upper(2,11) = jacobian(2,11)-lower(2,1)*upper(1,11)
upper(2,12) = jacobian(2,12)-lower(2,1)*upper(1,12)
upper(2,13) = jacobian(2,13)-lower(2,1)*upper(1,13)

lower(3,1) = jacobian(3,1)/upper(1,1)
lower(3,2) = (jacobian(3,2)-lower(3,1)*upper(1,2))/upper(2,2)
upper(3,3) = jacobian(3,3)-lower(3,1)*upper(1,3)-lower(3,2)*upper(2,3)
upper(3,4) = jacobian(3,4)-lower(3,1)*upper(1,4)-lower(3,2)*upper(2,4)
upper(3,5) = jacobian(3,5)-lower(3,1)*upper(1,5)-lower(3,2)*upper(2,5)
upper(3,6) = jacobian(3,6)-lower(3,1)*upper(1,6)-lower(3,2)*upper(2,6)
upper(3,7) = jacobian(3,7)-lower(3,1)*upper(1,7)-lower(3,2)*upper(2,7)
upper(3,8) = jacobian(3,8)-lower(3,1)*upper(1,8)-lower(3,2)*upper(2,8)
upper(3,9) = jacobian(3,9)-lower(3,1)*upper(1,9)-lower(3,2)*upper(2,9)
upper(3,10) = jacobian(3,10)-lower(3,1)*upper(1,10)-lower(3,2)*upper(2,10)
upper(3,11) = jacobian(3,11)-lower(3,1)*upper(1,11)-lower(3,2)*upper(2,11)
upper(3,12) = jacobian(3,12)-lower(3,1)*upper(1,12)-lower(3,2)*upper(2,12)
upper(3,13) = jacobian(3,13)-lower(3,1)*upper(1,13)-lower(3,2)*upper(2,13)

lower(4,1) = jacobian(4,1)/upper(1,1)
lower(4,2) = (jacobian(4,2)-lower(4,1)*upper(1,2))/upper(2,2)
lower(4,3) = (jacobian(4,3)-lower(4,1)*upper(1,3)-lower(4,2)*upper(2,3))/upper(3,3)
upper(4,4) = jacobian(4,4)-lower(4,1)*upper(1,4)-lower(4,2)*upper(2,4)-lower(4,3)*upper(3,4)
upper(4,5) = jacobian(4,5)-lower(4,1)*upper(1,5)-lower(4,2)*upper(2,5)-lower(4,3)*upper(3,5)
upper(4,6) = jacobian(4,6)-lower(4,1)*upper(1,6)-lower(4,2)*upper(2,6)-lower(4,3)*upper(3,6)
upper(4,7) = jacobian(4,7)-lower(4,1)*upper(1,7)-lower(4,2)*upper(2,7)-lower(4,3)*upper(3,7)
upper(4,8) = jacobian(4,8)-lower(4,1)*upper(1,8)-lower(4,2)*upper(2,8)-lower(4,3)*upper(3,8)
upper(4,9) = jacobian(4,9)-lower(4,1)*upper(1,9)-lower(4,2)*upper(2,9)-lower(4,3)*upper(3,9)
upper(4,10) = jacobian(4,10)-lower(4,1)*upper(1,10)-lower(4,2)*upper(2,10)-lower(4,3)*upper(3,10)
upper(4,11) = jacobian(4,11)-lower(4,1)*upper(1,11)-lower(4,2)*upper(2,11)-lower(4,3)*upper(3,11)
upper(4,12) = jacobian(4,12)-lower(4,1)*upper(1,12)-lower(4,2)*upper(2,12)-lower(4,3)*upper(3,12)
upper(4,13) = jacobian(4,13)-lower(4,1)*upper(1,13)-lower(4,2)*upper(2,13)-lower(4,3)*upper(3,13)
!
lower(5,1) = jacobian(5,1)/upper(1,1)
lower(5,2) = (jacobian(5,2)-lower(5,1)*upper(1,2))/upper(2,2)
lower(5,3) = (jacobian(5,3)-lower(5,1)*upper(1,3)-lower(5,2)*upper(2,3))/upper(3,3)
lower(5,4) = (jacobian(5,4)-lower(5,1)*upper(1,4)-lower(5,2)*upper(2,4)-lower(5,2)*upper(3,4))/upper(4,4)
upper(5,5) = jacobian(5,5)-lower(5,1)*upper(1,5)-lower(5,2)*upper(2,5)-lower(5,3)*upper(3,5)-lower(5,4)*upper(4,5)
upper(5,6) = jacobian(5,6)-lower(5,1)*upper(1,6)-lower(5,2)*upper(2,6)-lower(5,3)*upper(3,6)-lower(5,4)*upper(4,6)
upper(5,7) = jacobian(5,7)-lower(5,1)*upper(1,7)-lower(5,2)*upper(2,7)-lower(5,3)*upper(3,7)-lower(5,4)*upper(4,7)
upper(5,8) = jacobian(5,8)-lower(5,1)*upper(1,8)-lower(5,2)*upper(2,8)-lower(5,3)*upper(3,8)-lower(5,4)*upper(4,8)
upper(5,9) = jacobian(5,9)-lower(5,1)*upper(1,9)-lower(5,2)*upper(2,9)-lower(5,3)*upper(3,9)-lower(5,4)*upper(4,9)
upper(5,10) = jacobian(5,10)-lower(5,1)*upper(1,10)-lower(5,2)*upper(2,10)-lower(5,3)*upper(3,10)-lower(5,4)*upper(4,10)
upper(5,11) = jacobian(5,11)-lower(5,1)*upper(1,11)-lower(5,2)*upper(2,11)-lower(5,3)*upper(3,11)-lower(5,4)*upper(4,11)
upper(5,12) = jacobian(5,12)-lower(5,1)*upper(1,12)-lower(5,2)*upper(2,12)-lower(5,3)*upper(3,12)-lower(5,4)*upper(4,12)
upper(5,13) = jacobian(5,13)-lower(5,1)*upper(1,13)-lower(5,2)*upper(2,13)-lower(5,3)*upper(3,13)-lower(5,4)*upper(4,13)
!
lower(6,1) = jacobian(6,1)/upper(1,1)
lower(6,2) = (jacobian(6,2)-lower(6,1)*upper(1,2))/upper(2,2)
lower(6,3) = (jacobian(6,3)-lower(6,1)*upper(1,3)-lower(6,2)*upper(2,3))/upper(3,3)
lower(6,4) = (jacobian(6,4)-lower(6,1)*upper(1,4)-lower(6,2)*upper(2,4)-lower(6,3)*upper(3,4))/upper(4,4)
lower(6,5) = (jacobian(6,5)-lower(6,1)*upper(1,5)-lower(6,2)*upper(2,5)-lower(6,3)*upper(3,5)-lower(6,4)*upper(4,5))/upper(5,5)
upper(6,6) = jacobian(6,6)-lower(6,1)*upper(1,6)-lower(6,2)*upper(2,6)-lower(6,3)*upper(3,6)-lower(6,4)*upper(4,6)-&
  lower(6,5)*upper(5,6)
upper(6,7) = jacobian(6,7)-lower(6,1)*upper(1,7)-lower(6,2)*upper(2,7)-lower(6,3)*upper(3,7)-lower(6,4)*upper(4,7)-&
  lower(6,5)*upper(5,7)
upper(6,8) = jacobian(6,8)-lower(6,1)*upper(1,8)-lower(6,2)*upper(2,8)-lower(6,3)*upper(3,8)-lower(6,4)*upper(4,8)-&
  lower(6,5)*upper(5,8)
upper(6,9) = jacobian(6,9)-lower(6,1)*upper(1,9)-lower(6,2)*upper(2,9)-lower(6,3)*upper(3,9)-lower(6,4)*upper(4,9)-&
  lower(6,5)*upper(5,9)
upper(6,10) = jacobian(6,10)-lower(6,1)*upper(1,10)-lower(6,2)*upper(2,10)-lower(6,3)*upper(3,10)-lower(6,4)*upper(4,10)-&
  lower(6,5)*upper(5,10)
upper(6,11) = jacobian(6,11)-lower(6,1)*upper(1,11)-lower(6,2)*upper(2,11)-lower(6,3)*upper(3,11)-lower(6,4)*upper(4,11)-&
  lower(6,5)*upper(5,11)
upper(6,12) = jacobian(6,12)-lower(6,1)*upper(1,12)-lower(6,2)*upper(2,12)-lower(6,3)*upper(3,12)-lower(6,4)*upper(4,12)-&
  lower(6,5)*upper(5,12)
upper(6,13) = jacobian(6,13)-lower(6,1)*upper(1,13)-lower(6,2)*upper(2,13)-lower(6,3)*upper(3,13)-lower(6,4)*upper(4,13)-&
  lower(6,5)*upper(5,13)
!
!
!
lower(7,1) = jacobian(7,1)/upper(1,1)
lower(7,2) = (jacobian(7,2)-lower(7,1)*upper(1,2))/upper(2,2)
lower(7,3) = (jacobian(7,3)-lower(7,1)*upper(1,3)-lower(7,2)*upper(2,3))/upper(3,3)
lower(7,4) = (jacobian(7,4)-lower(7,1)*upper(1,4)-lower(7,2)*upper(2,4)-lower(7,3)*upper(3,4))/upper(4,4)
lower(7,5) = (jacobian(7,5)-lower(7,1)*upper(1,5)-lower(7,2)*upper(2,5)-lower(7,3)*upper(3,5)-lower(7,4)*upper(4,5))/upper(5,5)
lower(7,6) = (jacobian(7,6)-lower(7,1)*upper(1,6)-lower(7,2)*upper(2,6)-lower(7,3)*upper(3,6)-lower(7,4)*upper(4,6)-&
  lower(7,5)*upper(5,6))/upper(6,6)
!
upper(7,7) = jacobian(7,7)-lower(7,1)*upper(1,7)-lower(7,2)*upper(2,7)-lower(7,3)*upper(3,7)-lower(7,4)*upper(4,7) &
  -lower(7,5)*upper(5,7)-lower(7,6)*upper(6,7)
!
upper(7,8) = jacobian(7,8)-lower(7,1)*upper(1,8)-lower(7,2)*upper(2,8)-lower(7,3)*upper(3,8)-lower(7,4)*upper(4,8) &
  -lower(7,5)*upper(5,8)-lower(7,6)*upper(6,8)
!
upper(7,9) = jacobian(7,9)-lower(7,1)*upper(1,9)-lower(7,2)*upper(2,9)-lower(7,3)*upper(3,9)-lower(7,4)*upper(4,9) &
  -lower(7,5)*upper(5,9)-lower(7,6)*upper(6,9)
!
upper(7,10) = jacobian(7,10)-lower(7,1)*upper(1,10)-lower(7,2)*upper(2,10)-lower(7,3)*upper(3,10)-lower(7,4)*upper(4,10) &
  -lower(7,5)*upper(5,10)-lower(7,6)*upper(6,10)
!
upper(7,11) = jacobian(7,11)-lower(7,1)*upper(1,11)-lower(7,2)*upper(2,11)-lower(7,3)*upper(3,11)-lower(7,4)*upper(4,11) &
  -lower(7,5)*upper(5,11)-lower(7,6)*upper(6,11)
!
upper(7,12) = jacobian(7,12)-lower(7,1)*upper(1,12)-lower(7,2)*upper(2,12)-lower(7,3)*upper(3,12)-lower(7,4)*upper(4,12) &
  -lower(7,5)*upper(5,12)-lower(7,6)*upper(6,12)
!
upper(7,13) = jacobian(7,13)-lower(7,1)*upper(1,13)-lower(7,2)*upper(2,13)-lower(7,3)*upper(3,13)-lower(7,4)*upper(4,13) &
  -lower(7,5)*upper(5,13)-lower(7,6)*upper(6,13)
!
lower(8,1) = jacobian(8,1)/upper(1,1)
lower(8,2) = (jacobian(8,2)-lower(8,1)*upper(1,2))/upper(2,2)
lower(8,3) = (jacobian(8,3)-lower(8,1)*upper(1,3)-lower(8,2)*upper(2,3))/upper(3,3)
lower(8,4) = (jacobian(8,4)-lower(8,1)*upper(1,4)-lower(8,2)*upper(2,4)-lower(8,3)*upper(3,4))/upper(4,4)
lower(8,5) = (jacobian(8,5)-lower(8,1)*upper(1,5)-lower(8,2)*upper(2,5)-lower(8,3)*upper(3,5)-lower(8,4)*upper(4,5))/upper(5,5)
lower(8,6) = (jacobian(8,6)-lower(8,1)*upper(1,6)-lower(8,2)*upper(2,6)-lower(8,3)*upper(3,6)-lower(8,4)*upper(4,6) &
  -lower(8,5)*upper(5,6))/upper(6,6)
!
lower(8,7) = (jacobian(8,7)-lower(8,1)*upper(1,7)-lower(8,2)*upper(2,7)-lower(8,3)*upper(3,7)-lower(8,4)*upper(4,7) &
  -lower(8,5)*upper(5,7)-lower(8,6)*upper(6,7))/upper(7,7)
!
upper(8,8) = jacobian(8,8)-lower(8,1)*upper(1,8)-lower(8,2)*upper(2,8)-lower(8,3)*upper(3,8)-lower(8,4)*upper(4,8) &
  -lower(8,5)*upper(5,8)-lower(8,6)*upper(6,8)-lower(8,7)*upper(7,8)
!
upper(8,9) = jacobian(8,9)-lower(8,1)*upper(1,9)-lower(8,2)*upper(2,9)-lower(8,3)*upper(3,9)-lower(8,4)*upper(4,9) &
  -lower(8,5)*upper(5,9)-lower(8,6)*upper(6,9)-lower(8,7)*upper(7,9)
!
upper(8,10) = jacobian(8,10)-lower(8,1)*upper(1,10)-lower(8,2)*upper(2,10)-lower(8,3)*upper(3,10)-lower(8,4)*upper(4,10) &
  -lower(8,5)*upper(5,10)-lower(8,6)*upper(6,10)-lower(8,7)*upper(7,10)
!
upper(8,11) = jacobian(8,11)-lower(8,1)*upper(1,11)-lower(8,2)*upper(2,11)-lower(8,3)*upper(3,11)-lower(8,4)*upper(4,11) &
  -lower(8,5)*upper(5,11)-lower(8,6)*upper(6,11)-lower(8,7)*upper(7,11)
!
upper(8,12) = jacobian(8,12)-lower(8,1)*upper(1,12)-lower(8,2)*upper(2,12)-lower(8,3)*upper(3,12)-lower(8,4)*upper(4,12) &
  -lower(8,5)*upper(5,12)-lower(8,6)*upper(6,12)-lower(8,7)*upper(7,12)
!
upper(8,13) = jacobian(8,13)-lower(8,1)*upper(1,13)-lower(8,2)*upper(2,13)-lower(8,3)*upper(3,13)-lower(8,4)*upper(4,13) &
  -lower(8,5)*upper(5,13)-lower(8,6)*upper(6,13)-lower(8,7)*upper(7,13)
!
!
!
lower(9,1) = jacobian(9,1)/upper(1,1)
lower(9,2) = (jacobian(9,2)-lower(9,1)*upper(1,2))/upper(2,2)
lower(9,3) = (jacobian(9,3)-lower(9,1)*upper(1,3)-lower(9,2)*upper(2,3))/upper(3,3)
lower(9,4) = (jacobian(9,4)-lower(9,1)*upper(1,4)-lower(9,2)*upper(2,4)-lower(9,3)*upper(3,4))/upper(4,4)
lower(9,5) = (jacobian(9,5)-lower(9,1)*upper(1,5)-lower(9,2)*upper(2,5)-lower(9,3)*upper(3,5)-lower(9,4)*upper(4,5))/upper(5,5)
!
lower(9,6) = (jacobian(9,6)-lower(9,1)*upper(1,6)-lower(9,2)*upper(2,6)-lower(9,3)*upper(3,6)-lower(9,4)*upper(4,6) &
  -lower(9,5)*upper(5,6))/upper(6,6)
!
lower(9,7) = (jacobian(9,7)-lower(9,1)*upper(1,7)-lower(9,2)*upper(2,7)-lower(9,3)*upper(3,7)-lower(9,4)*upper(4,7) &
  -lower(9,5)*upper(5,7)-lower(9,6)*upper(6,7))/upper(7,7)
!
lower(9,8) = (jacobian(9,8)-lower(9,1)*upper(1,8)-lower(9,2)*upper(2,8)-lower(9,3)*upper(3,8)-lower(9,4)*upper(4,8) &
  -lower(9,5)*upper(5,8)-lower(9,6)*upper(6,8)-lower(9,7)*upper(7,8))/upper(8,8)
!
upper(9,9) = jacobian(9,9)-lower(9,1)*upper(1,9)-lower(9,2)*upper(2,9)-lower(9,3)*upper(3,9)-lower(9,4)*upper(4,9) &
  -lower(9,5)*upper(5,9)-lower(9,6)*upper(6,9)-lower(9,7)*upper(7,9)-lower(9,8)*upper(8,9)
!
upper(9,10) = jacobian(9,10)-lower(9,1)*upper(1,10)-lower(9,2)*upper(2,10)-lower(9,3)*upper(3,10)-lower(9,4)*upper(4,10) &
  -lower(9,5)*upper(5,10)-lower(9,6)*upper(6,10)-lower(9,7)*upper(7,10)-lower(9,8)*upper(8,10)
!
upper(9,11) = jacobian(9,11)-lower(9,1)*upper(1,11)-lower(9,2)*upper(2,11)-lower(9,3)*upper(3,11)-lower(9,4)*upper(4,11) &
  -lower(9,5)*upper(5,11)-lower(9,6)*upper(6,11)-lower(9,7)*upper(7,11)-lower(9,8)*upper(8,11)
!
upper(9,12) = jacobian(9,12)-lower(9,1)*upper(1,12)-lower(9,2)*upper(2,12)-lower(9,3)*upper(3,12)-lower(9,4)*upper(4,12) &
  -lower(9,5)*upper(5,12)-lower(9,6)*upper(6,12)-lower(9,7)*upper(7,12)-lower(9,8)*upper(8,12)
!
upper(9,13) = jacobian(9,13)-lower(9,1)*upper(1,13)-lower(9,2)*upper(2,13)-lower(9,3)*upper(3,13)-lower(9,4)*upper(4,13) &
  -lower(9,5)*upper(5,13)-lower(9,6)*upper(6,13)-lower(9,7)*upper(7,13)-lower(9,8)*upper(8,13)
!
!
!
lower(10,1) = jacobian(10,1)/upper(1,1)
lower(10,2) = (jacobian(10,2)-lower(10,1)*upper(1,2))/upper(2,2)
lower(10,3) = (jacobian(10,3)-lower(10,1)*upper(1,3)-lower(10,2)*upper(2,3))/upper(3,3)
lower(10,4) = (jacobian(10,4)-lower(10,1)*upper(1,4)-lower(10,2)*upper(2,4)-lower(10,3)*upper(3,4))/upper(4,4)
lower(10,5) = (jacobian(10,5)-lower(10,1)*upper(1,5)-lower(10,2)*upper(2,5)-lower(10,3)*upper(3,5)-&
  lower(10,4)*upper(4,5))/upper(5,5)
!
lower(10,6) = (jacobian(10,6)-lower(10,1)*upper(1,6)-lower(10,2)*upper(2,6)-lower(10,3)*upper(3,6)-lower(10,4)*upper(4,6) &
  -lower(10,5)*upper(5,6))/upper(6,6)
!
lower(10,7) = (jacobian(10,7)-lower(10,1)*upper(1,7)-lower(10,2)*upper(2,7)-lower(10,3)*upper(3,7)-lower(10,4)*upper(4,7) &
  -lower(10,5)*upper(5,7)-lower(10,6)*upper(6,7))/upper(7,7)
!
lower(10,8) = (jacobian(10,8)-lower(10,1)*upper(1,8)-lower(10,2)*upper(2,8)-lower(10,3)*upper(3,8)-lower(10,4)*upper(4,8) &
  -lower(10,5)*upper(5,8)-lower(10,6)*upper(6,8)-lower(10,7)*upper(7,8))/upper(8,8)
!
lower(10,9) = (jacobian(10,9)-lower(10,1)*upper(1,9)-lower(10,2)*upper(2,9)-lower(10,3)*upper(3,9)-lower(10,4)*upper(4,9) &
  -lower(10,5)*upper(5,9)-lower(10,6)*upper(6,9)-lower(10,7)*upper(7,9)-lower(10,8)*upper(8,9))/upper(9,9)
!
upper(10,10) = jacobian(10,10)-lower(10,1)*upper(1,10)-lower(10,2)*upper(2,10)-lower(10,3)*upper(3,10)-lower(10,4)*upper(4,10) &
  -lower(10,5)*upper(5,10)-lower(10,6)*upper(6,10)-lower(10,7)*upper(7,10)-lower(10,8)*upper(8,10)-lower(10,9)*upper(9,10)
!
upper(10,11) = jacobian(10,11)-lower(10,1)*upper(1,11)-lower(10,2)*upper(2,11)-lower(10,3)*upper(3,11)-lower(10,4)*upper(4,11) &
  -lower(10,5)*upper(5,11)-lower(10,6)*upper(6,11)-lower(10,7)*upper(7,11)-lower(10,8)*upper(8,11)-lower(10,9)*upper(9,11)
!
upper(10,12) = jacobian(10,12)-lower(10,1)*upper(1,12)-lower(10,2)*upper(2,12)-lower(10,3)*upper(3,12)-lower(10,4)*upper(4,12) &
  -lower(10,5)*upper(5,12)-lower(10,6)*upper(6,12)-lower(10,7)*upper(7,12)-lower(10,8)*upper(8,12)-lower(10,9)*upper(9,12)
!
upper(10,13) = jacobian(10,13)-lower(10,1)*upper(1,13)-lower(10,2)*upper(2,13)-lower(10,3)*upper(3,13)-lower(10,4)*upper(4,13) &
  -lower(10,5)*upper(5,13)-lower(10,6)*upper(6,13)-lower(10,7)*upper(7,13)-lower(10,8)*upper(8,13)-lower(10,9)*upper(9,13)
!
!
!
lower(11,1) = jacobian(11,1)/upper(1,1)
lower(11,2) = (jacobian(11,2)-lower(11,1)*upper(1,2))/upper(2,2)
lower(11,3) = (jacobian(11,3)-lower(11,1)*upper(1,3)-lower(11,2)*upper(2,3))/upper(3,3)
lower(11,4) = (jacobian(11,4)-lower(11,1)*upper(1,4)-lower(11,2)*upper(2,4)-lower(11,3)*upper(3,4))/upper(4,4)
lower(11,5) = (jacobian(11,5)-lower(11,1)*upper(1,5)-lower(11,2)*upper(2,5)-lower(11,3)*upper(3,5)-&
  lower(11,4)*upper(4,5))/upper(5,5)
!
lower(11,6) = (jacobian(11,6)-lower(11,1)*upper(1,6)-lower(11,2)*upper(2,6)-lower(11,3)*upper(3,6)-lower(11,4)*upper(4,6) &
  -lower(11,5)*upper(5,6))/upper(6,6)
!
lower(11,7) = (jacobian(11,7)-lower(11,1)*upper(1,7)-lower(11,2)*upper(2,7)-lower(11,3)*upper(3,7)-lower(11,4)*upper(4,7) &
  -lower(11,5)*upper(5,7)-lower(11,6)*upper(6,7))/upper(7,7)
!
lower(11,8) = (jacobian(11,8)-lower(11,1)*upper(1,8)-lower(11,2)*upper(2,8)-lower(11,3)*upper(3,8)-lower(11,4)*upper(4,8) &
  -lower(11,5)*upper(5,8)-lower(11,6)*upper(6,8)-lower(11,7)*upper(7,8))/upper(8,8)
!
lower(11,9) = (jacobian(11,9)-lower(11,1)*upper(1,9)-lower(11,2)*upper(2,9)-lower(11,3)*upper(3,9)-lower(11,4)*upper(4,9) &
  -lower(11,5)*upper(5,9)-lower(11,6)*upper(6,9)-lower(11,7)*upper(7,9)-lower(11,8)*upper(8,9))/upper(9,9)
!
lower(11,10) = (jacobian(11,10)-lower(11,1)*upper(1,10)-lower(11,2)*upper(2,10)-lower(11,3)*upper(3,10)-lower(11,4)*upper(4,10) &
  -lower(11,5)*upper(5,10)-lower(11,6)*upper(6,10)-lower(11,7)*upper(7,10)-lower(11,8)*upper(8,10)-&
  lower(11,9)*upper(9,10))/upper(10,10)
!
upper(11,11) = jacobian(11,11)-lower(11,1)*upper(1,11)-lower(11,2)*upper(2,11)-lower(11,3)*upper(3,11)-lower(11,4)*upper(4,11) &
  -lower(11,5)*upper(5,11)-lower(11,6)*upper(6,11)-lower(11,7)*upper(7,11)-lower(11,8)*upper(8,11)-lower(11,9)*upper(9,11) &
  -lower(11,10)*upper(10,11)
!
upper(11,12) = jacobian(11,12)-lower(11,1)*upper(1,12)-lower(11,2)*upper(2,12)-lower(11,3)*upper(3,12)-lower(11,4)*upper(4,12) &
  -lower(11,5)*upper(5,12)-lower(11,6)*upper(6,12)-lower(11,7)*upper(7,12)-lower(11,8)*upper(8,12)-lower(11,9)*upper(9,12) &
  -lower(11,10)*upper(10,12)
!
upper(11,13) = jacobian(11,13)-lower(11,1)*upper(1,13)-lower(11,2)*upper(2,13)-lower(11,3)*upper(3,13)-lower(11,4)*upper(4,13) &
  -lower(11,5)*upper(5,13)-lower(11,6)*upper(6,13)-lower(11,7)*upper(7,13)-lower(11,8)*upper(8,13)-lower(11,9)*upper(9,13) &
  -lower(11,10)*upper(10,13)
!
!
!
lower(12,1) = jacobian(12,1)/upper(1,1)
lower(12,2) = (jacobian(12,2)-lower(12,1)*upper(1,2))/upper(2,2)
lower(12,3) = (jacobian(12,3)-lower(12,1)*upper(1,3)-lower(12,2)*upper(2,3))/upper(3,3)
lower(12,4) = (jacobian(12,4)-lower(12,1)*upper(1,4)-lower(12,2)*upper(2,4)-lower(12,3)*upper(3,4))/upper(4,4)
lower(12,5) = (jacobian(12,5)-lower(12,1)*upper(1,5)-lower(12,2)*upper(2,5)-lower(12,3)*upper(3,5)-&
  lower(12,4)*upper(4,5))/upper(5,5)
!
lower(12,6) = (jacobian(12,6)-lower(12,1)*upper(1,6)-lower(12,2)*upper(2,6)-lower(12,3)*upper(3,6)-lower(12,4)*upper(4,6) &
  -lower(12,5)*upper(5,6))/upper(6,6)
!
lower(12,7) = (jacobian(12,7)-lower(12,1)*upper(1,7)-lower(12,2)*upper(2,7)-lower(12,3)*upper(3,7)-lower(12,4)*upper(4,7) &
  -lower(12,5)*upper(5,7)-lower(12,6)*upper(6,7))/upper(7,7)
!
lower(12,8) = (jacobian(12,8)-lower(12,1)*upper(1,8)-lower(12,2)*upper(2,8)-lower(12,3)*upper(3,8)-lower(12,4)*upper(4,8) &
  -lower(12,5)*upper(5,8)-lower(12,6)*upper(6,8)-lower(12,7)*upper(7,8))/upper(8,8)
!
lower(12,9) = (jacobian(12,9)-lower(12,1)*upper(1,9)-lower(12,2)*upper(2,9)-lower(12,3)*upper(3,9)-lower(12,4)*upper(4,9) &
  -lower(12,5)*upper(5,9)-lower(12,6)*upper(6,9)-lower(12,7)*upper(7,9)-lower(12,8)*upper(8,9))/upper(9,9)
!
lower(12,10) = (jacobian(12,10)-lower(12,1)*upper(1,10)-lower(12,2)*upper(2,10)-lower(12,3)*upper(3,10)-lower(12,4)*upper(4,10) &
  -lower(12,5)*upper(5,10)-lower(12,6)*upper(6,10)-lower(12,7)*upper(7,10)-lower(12,8)*upper(8,10)-&
  lower(12,9)*upper(9,10))/upper(10,10)
!
lower(12,11) = (jacobian(12,11)-lower(12,1)*upper(1,11)-lower(12,2)*upper(2,11)-lower(12,3)*upper(3,11)-lower(12,4)*upper(4,11) &
  -lower(12,5)*upper(5,11)-lower(12,6)*upper(6,11)-lower(12,7)*upper(7,11)-lower(12,8)*upper(8,11)-lower(12,9)*upper(9,11) &
  -lower(12,10)*upper(10,11))/upper(11,11)
!
upper(12,12) = jacobian(12,12)-lower(12,1)*upper(1,12)-lower(12,2)*upper(2,12)-lower(12,3)*upper(3,12)-lower(12,4)*upper(4,12) &
  -lower(12,5)*upper(5,12)-lower(12,6)*upper(6,12)-lower(12,7)*upper(7,12)-lower(12,8)*upper(8,12)-lower(12,9)*upper(9,12) &
  -lower(12,10)*upper(10,12)-lower(12,11)*upper(11,12)
!
upper(12,13) = jacobian(12,13)-lower(12,1)*upper(1,13)-lower(12,2)*upper(2,13)-lower(12,3)*upper(3,13)-lower(12,4)*upper(4,13) &
  -lower(12,5)*upper(5,13)-lower(12,6)*upper(6,13)-lower(12,7)*upper(7,13)-lower(12,8)*upper(8,13)-lower(12,9)*upper(9,13) &
  -lower(12,10)*upper(10,13)-lower(12,11)*upper(11,13)
!
!
!
lower(13,1) = jacobian(13,1)/upper(1,1)
lower(13,2) = (jacobian(13,2)-lower(13,1)*upper(1,2))/upper(2,2)
lower(13,3) = (jacobian(13,3)-lower(13,1)*upper(1,3)-lower(13,2)*upper(2,3))/upper(3,3)
lower(13,4) = (jacobian(13,4)-lower(13,1)*upper(1,4)-lower(13,2)*upper(2,4)-lower(13,3)*upper(3,4))/upper(4,4)
lower(13,5) = (jacobian(13,5)-lower(13,1)*upper(1,5)-lower(13,2)*upper(2,5)-lower(13,3)*upper(3,5)-&
  lower(13,4)*upper(4,5))/upper(5,5)
!
lower(13,6) = (jacobian(13,6)-lower(13,1)*upper(1,6)-lower(13,2)*upper(2,6)-lower(13,3)*upper(3,6)-lower(13,4)*upper(4,6) &
  -lower(13,5)*upper(5,6))/upper(6,6)
!
lower(13,7) = (jacobian(13,7)-lower(13,1)*upper(1,7)-lower(13,2)*upper(2,7)-lower(13,3)*upper(3,7)-lower(13,4)*upper(4,7) &
  -lower(13,5)*upper(5,7)-lower(13,6)*upper(6,7))/upper(7,7)
!
lower(13,8) = (jacobian(13,8)-lower(13,1)*upper(1,8)-lower(13,2)*upper(2,8)-lower(13,3)*upper(3,8)-lower(13,4)*upper(4,8) &
  -lower(13,5)*upper(5,8)-lower(13,6)*upper(6,8)-lower(13,7)*upper(7,8))/upper(8,8)
!
lower(13,9) = (jacobian(13,9)-lower(13,1)*upper(1,9)-lower(13,2)*upper(2,9)-lower(13,3)*upper(3,9)-lower(13,4)*upper(4,9) &
  -lower(13,5)*upper(5,9)-lower(13,6)*upper(6,9)-lower(13,7)*upper(7,9)-lower(13,8)*upper(8,9))/upper(9,9)
!
lower(13,10) = (jacobian(13,10)-lower(13,1)*upper(1,10)-lower(13,2)*upper(2,10)-lower(13,3)*upper(3,10)-lower(13,4)*upper(4,10) &
  -lower(13,5)*upper(5,10)-lower(13,6)*upper(6,10)-lower(13,7)*upper(7,10)-lower(13,8)*upper(8,10)-&
  lower(13,9)*upper(9,10))/upper(10,10)
!
lower(13,11) = (jacobian(13,11)-lower(13,1)*upper(1,11)-lower(13,2)*upper(2,11)-lower(13,3)*upper(3,11)-lower(13,4)*upper(4,11) &
  -lower(13,5)*upper(5,11)-lower(13,6)*upper(6,11)-lower(13,7)*upper(7,11)-lower(13,8)*upper(8,11)-lower(13,9)*upper(9,11) &
  -lower(13,10)*upper(10,11))/upper(11,11)
!
lower(13,12) = (jacobian(13,12)-lower(13,1)*upper(1,12)-lower(13,2)*upper(2,12)-lower(13,3)*upper(3,12)-lower(13,4)*upper(4,12) &
  -lower(13,5)*upper(5,12)-lower(13,6)*upper(6,12)-lower(13,7)*upper(7,12)-lower(13,8)*upper(8,12)-lower(13,9)*upper(9,12) &
  -lower(13,10)*upper(10,12)-lower(13,11)*upper(11,12))/upper(12,12)
!
upper(13,13) = jacobian(13,13)-lower(13,1)*upper(1,13)-lower(13,2)*upper(2,13)-lower(13,3)*upper(3,13)-lower(13,4)*upper(4,13) &
  -lower(13,5)*upper(5,13)-lower(13,6)*upper(6,13)-lower(13,7)*upper(7,13)-lower(13,8)*upper(8,13)-lower(13,9)*upper(9,13) &
  -lower(13,10)*upper(10,13)-lower(13,11)*upper(11,13)-lower(13,12)*upper(12,13)
!
!
!
!write(*,*) 'lower matrix'
!write(*,*) lower(1,1)
!write(*,*) lower(2,1),lower(2,2)
!write(*,*) lower(3,1),lower(3,2),lower(3,3)
!write(*,*) lower(4,1),lower(4,2),lower(4,3),lower(4,4)
!write(*,*) lower(5,1),lower(5,2),lower(5,3),lower(5,4),lower(5,5)
!write(*,*) lower(6,1),lower(6,2),lower(6,3),lower(6,4),lower(6,5),lower(6,6)
!write(*,*) lower(7,1),lower(7,2),lower(7,3),lower(7,4),lower(7,5),lower(7,6),lower(7,7)
!write(*,*) lower(8,1),lower(8,2),lower(8,3),lower(8,4),lower(8,5),lower(8,6),lower(8,7),lower(8,8)
!write(*,*) lower(9,1),lower(9,2),lower(9,3),lower(9,4),lower(9,5),lower(9,6),lower(9,7),lower(9,8),lower(9,9)
!write(*,*) lower(10,1),lower(10,2),lower(10,3),lower(10,4),lower(10,5),lower(10,6),lower(10,7),&
!  lower(10,8),lower(10,9),lower(10,10)
!write(*,*) lower(11,1),lower(11,2),lower(11,3),lower(11,4),lower(11,5),lower(11,6),lower(11,7),&
!  lower(11,8),lower(11,9),lower(11,10),lower(11,11)
!write(*,*) lower(12,1),lower(12,2),lower(12,3),lower(12,4),lower(12,5),lower(12,6),lower(12,7),&
!  lower(12,8),lower(12,9),lower(12,10),lower(12,11),lower(12,12)
!write(*,*) lower(13,1),lower(13,2),lower(13,3),lower(13,4),lower(13,5),lower(13,6),lower(13,7),&
!  lower(13,8),lower(13,9),lower(13,10),lower(13,11),lower(13,12),lower(13,13)
!write(*,*)
!write(*,*)
!write(*,*) 'upper matrix'
!write(*,*) upper(1,1),upper(1,2),upper(1,3),upper(1,4),upper(1,5),upper(1,6),upper(1,7),upper(1,8)&
!  ,upper(1,9),upper(1,10),upper(1,11),upper(1,12),upper(1,13)
!write(*,*) upper(2,2),upper(2,3),upper(2,4),upper(2,5),upper(2,6),upper(2,7),upper(2,8)&
!  ,upper(2,9),upper(2,10),upper(2,11),upper(2,12),upper(2,13)
!write(*,*) upper(3,3),upper(3,4),upper(3,5),upper(3,6),upper(3,7),upper(3,8)&
!  ,upper(3,9),upper(3,10),upper(3,11),upper(3,12),upper(3,13)
!write(*,*) upper(4,4),upper(4,5),upper(4,6),upper(4,7),upper(4,8)&
!  ,upper(4,9),upper(4,10),upper(4,11),upper(4,12),upper(4,13)
!write(*,*) upper(5,5),upper(5,6),upper(5,7),upper(5,8),upper(5,9),&
!  upper(5,10),upper(5,11),upper(5,12),upper(5,13)
!write(*,*) upper(6,6),upper(6,7),upper(6,8),upper(6,9),upper(6,10),&
!  upper(6,11),upper(6,12),upper(6,13)
!write(*,*) upper(7,7),upper(7,8),upper(7,9),upper(7,10),upper(7,11),&
!  upper(7,12),upper(7,13)
!write(*,*) upper(8,8),upper(8,9),upper(8,10),upper(8,11),upper(8,12),upper(8,13)
!write(*,*) upper(9,9),upper(9,10),upper(9,11),upper(9,12),upper(9,13)
!write(*,*) upper(10,10),upper(10,11),upper(10,12),upper(10,13)
!write(*,*) upper(11,11),upper(11,12),upper(11,13)
!write(*,*) upper(12,12),upper(12,13)
!write(*,*) upper(13,13)
!write(*,*)
!
!
!
!
!
!
z1 = (jacobian(1,1)-loc_rho)!/lower(1,1), rho

z2 = (jacobian(1,2)-loc_rho*loc_u-lower(2,1)*z1)!/lower(2,2), x-momentum

z3 = (jacobian(1,3)-loc_rho*loc_v-lower(3,1)*z1-lower(3,2)*z2)!/lower(3,3), y-momentum

z4 = (jacobian(1,4)-loc_rho*loc_w-lower(4,1)*z1-lower(4,2)*z2-lower(4,3)*z3)!/lower(4,4), z-momentum

z5 = (jacobian(1,5)-loc_pxx-lower(5,1)*z1-lower(5,2)*z2-lower(5,3)*z3-lower(5,4)*z4)!/lower(5,5), Pxx

z6 = (jacobian(1,6)-loc_pyy-lower(6,1)*z1-lower(6,2)*z2-lower(6,3)*z3-lower(6,4)*z4-lower(6,5)*z5)!/lower(6,6), Pyy

z7 = (jacobian(1,7)-loc_pzz-lower(7,1)*z1-lower(7,2)*z2-lower(7,3)*z3-lower(7,4)*z4-lower(7,5)*z5-lower(7,6)*z6)!/lower(7,7), Pzz

z8 = (jacobian(1,8)-loc_pxy-lower(8,1)*z1-lower(8,2)*z2-lower(8,3)*z3-lower(8,4)*z4-lower(8,5)*z5-lower(8,6)*z6-lower(8,7)*z7)!/lower(8,8), Pxy

z9 = (jacobian(1,9)-loc_pxz-lower(9,1)*z1-lower(9,2)*z2-lower(9,3)*z3-lower(9,4)*z4-lower(9,5)*z5-&
  lower(9,6)*z6-lower(9,7)*z8 -lower(9,8)*z8)!/lower(9,9), Pxz

z10 = (jacobian(1,10)-loc_pyz-lower(10,1)*z1-lower(10,2)*z2-lower(10,3)*z3-lower(10,4)*z4- &
  lower(10,5)*z5 -lower(10,6)*z6-lower(10,7)*z7 -lower(10,8)*z8 -lower(10,9)*z9)!/lower(10,10), Pyz

z11 = (jacobian(1,11)-loc_q_x-lower(11,1)*z1-lower(11,2)*z2-lower(11,3)*z3-lower(11,4)*z4- &
  lower(11,5)*z5-lower(11,6)*z6-lower(11,7)*z7 &
  -lower(11,8)*z8-lower(11,9)*z9-lower(11,10)*z10)!/lower(11,11), q_x

z12 = (jacobian(1,12)-loc_q_y-lower(12,1)*z1-lower(12,2)*z2-lower(12,3)*z3-lower(12,4)*z4- &
  lower(12,5)*z5-lower(12,6)*z6-lower(12,7)*z7 &
  -lower(12,8)*z8-lower(12,9)*z9-lower(12,10)*z10-lower(12,11)*z11)!/lower(12,12), q_y

z13 = (jacobian(1,13)-loc_q_z-lower(13,1)*z1-lower(13,2)*z2-lower(13,3)*z3-lower(13,4)*z4- &
  lower(13,5)*z5 -lower(13,6)*z6 -lower(13,7)*z7 &
  -lower(13,8)*z8-lower(13,9)*z9 -lower(13,10)*z10 -lower(13,11)*z11 -lower(13,12)*z12)!/lower(13,13), q_z
!
!write(*,*) 'z values',z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13
!
!
delta_lambda_z = z13/upper(13,13)
!
delta_lambda_y = (z12-upper(12,13)*delta_lambda_z)/upper(12,12)
!
delta_lambda_x = (z11-upper(11,13)*delta_lambda_z-upper(11,12)*delta_lambda_y)/upper(11,11)
!
delta_piyz = (z10-upper(10,13)*delta_lambda_z-upper(10,12)*delta_lambda_y - upper(10,11)*delta_lambda_x)/upper(10,10)
!
delta_pixz = (z9-upper(9,13)*delta_lambda_z-upper(9,12)*delta_lambda_y - upper(9,11)*delta_lambda_x-&
   upper(9,10)*delta_piyz)/upper(9,9)
!
delta_pixy = (z8-upper(8,13)*delta_lambda_z-upper(8,12)*delta_lambda_y - upper(8,11)*delta_lambda_x- upper(8,10)*delta_piyz &
        - upper(8,9)*delta_pixz)/upper(8,8)
!
delta_pizza = (z7-upper(7,13)*delta_lambda_z-upper(7,12)*delta_lambda_y - upper(7,11)*delta_lambda_x- upper(7,10)*delta_piyz &
  - upper(7,9)*delta_pixz - upper(7,8)*delta_pixy)/upper(7,7)
!
delta_piyy = (z6-upper(6,13)*delta_lambda_z - upper(6,12)*delta_lambda_y - upper(6,11)*delta_lambda_x - upper(6,10)*delta_piyz &
  - upper(6,9)*delta_pixz - upper(6,8)*delta_pixy - upper(6,7)*delta_pizza)/upper(6,6)
!
delta_pixx = (z5-upper(5,13)*delta_lambda_z-upper(5,12)*delta_lambda_y - upper(5,11)*delta_lambda_x- upper(5,10)*delta_piyz &
  - upper(5,9)*delta_pixz - upper(5,8)*delta_pixy - upper(5,7)*delta_pizza - upper(5,6)*delta_piyy)/upper(5,5)
!
delta_zeta_z = (z4-upper(4,13)*delta_lambda_z-upper(4,12)*delta_lambda_y - upper(4,11)*delta_lambda_x- upper(4,10)*delta_piyz &
  - upper(4,9)*delta_pixz - upper(4,8)*delta_pixy - upper(4,7)*delta_pizza - upper(4,6)*delta_piyy -&
  upper(4,5)*delta_pixx)/upper(4,4)
!
delta_zeta_y = (z3-upper(3,13)*delta_lambda_z-upper(3,12)*delta_lambda_y - upper(3,11)*delta_lambda_x- upper(3,10)*delta_piyz &
  - upper(3,9)*delta_pixz - upper(3,8)*delta_pixy - upper(3,7)*delta_pizza - upper(3,6)*delta_piyy- upper(3,5)*delta_pixx &
  - upper(3,4)*delta_zeta_z)/upper(3,3)
!
delta_zeta_x = (z2-upper(2,13)*delta_lambda_z-upper(2,12)*delta_lambda_y - upper(2,11)*delta_lambda_x- upper(2,10)*delta_piyz &
  - upper(2,9)*delta_pixz - upper(2,8)*delta_pixy - upper(2,7)*delta_pizza - upper(2,6)*delta_piyy - upper(2,5)*delta_pixx &
  - upper(2,4)*delta_zeta_z - upper(2,3)*delta_zeta_y)/upper(2,2)
!
delta_chi = (z1-upper(1,13)*delta_lambda_z-upper(1,12)*delta_lambda_y - upper(1,11)*delta_lambda_x- upper(1,10)*delta_piyz &
  - upper(1,9)*delta_pixz - upper(1,8)*delta_pixy - upper(1,7)*delta_pizza - upper(1,6)*delta_piyy - upper(1,5)*delta_pixx &
  - upper(1,4)*delta_zeta_z - upper(1,3)*delta_zeta_y - upper(1,2)*delta_zeta_x)/upper(1,1)
 !
max_diff=MAX(ABS(delta_chi),ABS(delta_zeta_x),ABS(delta_zeta_y),ABS(delta_zeta_z),ABS(delta_pixx),&
 ABS(delta_piyy),ABS(delta_pizza),ABS(delta_pixy),ABS(delta_pixz),&
 ABS(delta_piyz),ABS(delta_lambda_x),ABS(delta_lambda_y),ABS(delta_lambda_z))
!
!
!write(*,*) 'LM changes',delta_chi,delta_zeta_x,delta_zeta_y,delta_zeta_z,delta_pixx,delta_piyy,delta_pizza,delta_pixy,delta_pixz,&
!        delta_piyz,delta_lambda_x,delta_lambda_y,delta_lambda_z
!
!
loc_chi = loc_chi-delta_chi
loc_zeta_x = loc_zeta_x-delta_zeta_x
loc_zeta_y = loc_zeta_y-delta_zeta_y
loc_zeta_z = loc_zeta_z-delta_zeta_z
loc_pi_xx = loc_pi_xx-delta_pixx
loc_pi_yy = loc_pi_yy-delta_piyy
loc_pi_zz = loc_pi_zz-delta_pizza
loc_pi_xy = loc_pi_xy-delta_pixy
loc_pi_xz = loc_pi_xz-delta_pixz
loc_pi_yz = loc_pi_yz-delta_piyz
loc_lambda_x = loc_lambda_x-delta_lambda_x
loc_lambda_y = loc_lambda_y-delta_lambda_y
loc_lambda_z = loc_lambda_z-delta_lambda_z
!
!write(*,*) "new chi",loc_chi
!write(*,*) "new zeta x",loc_zeta_x
!write(*,*) "new zeta y",loc_zeta_y
!write(*,*) "new zeta z",loc_zeta_z
!write(*,*) "new pi xx",loc_pi_xx
!write(*,*) "new pi yy",loc_pi_yy
!write(*,*) "new pi zz",loc_pi_zz
!write(*,*) "new pi xy",loc_pi_xy
!write(*,*) "new pi xz",loc_pi_xz
!write(*,*) "new pi yz",loc_pi_yz
!write(*,*) "new lambda x",loc_lambda_x
!write(*,*) "new lambda y",loc_lambda_y
!write(*,*) "new lambda z",loc_lambda_z
!
END SUBROUTINE
!
! indices of the second term are inverted because I'm only solving for one matrix
! instead of both
!
!      lower(1,1) = SQRT(jacobian(1,1))
!      lower(2,1) = jacobian(2,1)/lower(1,1)
!      lower(3,1) = jacobian(3,1)/lower(1,1)
!      lower(4,1) = jacobian(4,1)/lower(1,1)
!      lower(5,1) = jacobian(5,1)/lower(1,1)
!      lower(6,1) = jacobian(6,1)/lower(1,1)
!      lower(7,1) = jacobian(7,1)/lower(1,1)
!      lower(8,1) = jacobian(8,1)/lower(1,1)
!      lower(9,1) = jacobian(9,1)/lower(1,1)
!      lower(10,1) = jacobian(10,1)/lower(1,1)
!      lower(11,1) = jacobian(11,1)/lower(1,1)
!      lower(12,1) = jacobian(12,1)/lower(1,1)
!      lower(13,1) = jacobian(13,1)/lower(1,1)
!!
!      lower(2,2) = SQRT(jacobian(2,2)-lower(2,1)**2)
!      lower(3,2) = (jacobian(2,3)-lower(2,1)*lower(3,1))/lower(2,2) !
!      lower(4,2) = (jacobian(2,4)-lower(2,1)*lower(4,1))/lower(2,2)
!      lower(5,2) = (jacobian(2,5)-lower(2,1)*lower(5,1))/lower(2,2)
!      lower(6,2) = (jacobian(2,6)-lower(2,1)*lower(6,1))/lower(2,2)
!      lower(7,2) = (jacobian(2,7)-lower(2,1)*lower(7,1))/lower(2,2)
!      lower(8,2) = (jacobian(2,8)-lower(2,1)*lower(8,1))/lower(2,2)
!      lower(9,2) = (jacobian(2,9)-lower(2,1)*lower(9,1))/lower(2,2)
!      lower(10,2) = (jacobian(2,10)-lower(2,1)*lower(10,1))/lower(2,2)
!      lower(11,2) = (jacobian(2,11)-lower(2,1)*lower(11,1))/lower(2,2)
!      lower(12,2) = (jacobian(2,12)-lower(2,1)*lower(12,1))/lower(2,2)
!      lower(13,2) = (jacobian(2,13)-lower(2,1)*lower(13,1))/lower(2,2)
!
!      lower(3,3) = SQRT(jacobian(3,3)-lower(3,1)**2-lower(3,2)**2)
!      lower(4,3) = (jacobian(3,4)-lower(3,1)*lower(4,1)-lower(3,2)*lower(4,2))/lower(3,3)
!      lower(5,3) = (jacobian(3,5)-lower(3,1)*lower(5,1)-lower(3,2)*lower(5,2))/lower(3,3)
!      lower(6,3) = (jacobian(3,6)-lower(3,1)*lower(6,1)-lower(3,2)*lower(6,2))/lower(3,3)
!      lower(7,3) = (jacobian(3,7)-lower(3,1)*lower(7,1)-lower(3,2)*lower(7,2))/lower(3,3)
!      lower(8,3) = (jacobian(3,8)-lower(3,1)*lower(8,1)-lower(3,2)*lower(8,2))/lower(3,3)
!      lower(9,3) = (jacobian(3,9)-lower(3,1)*lower(9,1)-lower(3,2)*lower(9,2))/lower(3,3)
!      lower(10,3) = (jacobian(3,10)-lower(3,1)*lower(10,1)-lower(3,2)*lower(10,2))/lower(3,3)
!      lower(11,3) = (jacobian(3,11)-lower(3,1)*lower(11,1)-lower(3,2)*lower(11,2))/lower(3,3)
!      lower(12,3) = (jacobian(3,12)-lower(3,1)*lower(12,1)-lower(3,2)*lower(12,2))/lower(3,3)
!      lower(13,3) = (jacobian(3,13)-lower(3,1)*lower(13,1)-lower(3,2)*lower(13,2))/lower(3,3)
!!
!      lower(4,4) = SQRT(jacobian(4,4)-lower(4,1)**2-lower(4,2)**2 -lower(4,3)**2)
!      lower(5,4) = (jacobian(4,5)-lower(4,1)*lower(5,1)-lower(4,2)*lower(5,2)-lower(4,3)*lower(5,3))/lower(4,4)
!      lower(6,4) = (jacobian(4,6)-lower(4,1)*lower(6,1)-lower(4,2)*lower(6,2)-lower(4,3)*lower(6,3))/lower(4,4)
!      lower(7,4) = (jacobian(4,7)-lower(4,1)*lower(7,1)-lower(4,2)*lower(7,2)-lower(4,3)*lower(7,3))/lower(4,4)
!      lower(8,4) = (jacobian(4,8)-lower(4,1)*lower(8,1)-lower(4,2)*lower(8,2)-lower(4,3)*lower(8,3))/lower(4,4)
!      lower(9,4) = (jacobian(4,9)-lower(4,1)*lower(9,1)-lower(4,2)*lower(9,2)-lower(4,3)*lower(9,3))/lower(4,4)
!      lower(10,4) = (jacobian(4,10)-lower(4,1)*lower(10,1)-lower(4,2)*lower(10,2)-lower(4,3)*lower(10,3))/lower(4,4)
!      lower(11,3) = (jacobian(4,11)-lower(4,1)*lower(11,1)-lower(4,2)*lower(11,2)-lower(4,3)*lower(11,3))/lower(4,4)
!      lower(12,3) = (jacobian(4,12)-lower(4,1)*lower(12,1)-lower(4,2)*lower(12,2)-lower(4,3)*lower(12,3))/lower(4,4)
!      lower(13,3) = (jacobian(4,13)-lower(4,1)*lower(13,1)-lower(4,2)*lower(13,2)-lower(4,3)*lower(13,3))/lower(4,4)
!!
!      lower(5,5) = SQRT(jacobian(5,5)-lower(5,1)**2-lower(5,2)**2 -lower(5,3)**2 - lower(5,4)**2)
!      lower(6,5) = (jacobian(5,6)-lower(5,1)*lower(6,1)-lower(5,2)*lower(6,2)-lower(5,3)*lower(6,3)-lower(5,4)*lower(6,4))/lower(5,5)
!      lower(7,5) = (jacobian(5,7)-lower(5,1)*lower(7,1)-lower(5,2)*lower(7,2)-lower(5,3)*lower(7,3)-lower(5,4)*lower(7,4))/lower(5,5)
!      lower(8,5) = (jacobian(5,8)-lower(5,1)*lower(8,1)-lower(5,2)*lower(8,2)-lower(5,3)*lower(8,3)-lower(5,4)*lower(8,4))/lower(5,5)
!      lower(9,5) = (jacobian(5,9)-lower(5,1)*lower(9,1)-lower(5,2)*lower(9,2)-lower(5,3)*lower(9,3)-lower(5,4)*lower(9,4))/lower(5,5)
!      lower(10,5) = (jacobian(5,10)-lower(5,1)*lower(10,1)-lower(5,2)*lower(10,2)-lower(5,3)*lower(10,3)-lower(5,4)*lower(10,4))/lower(5,5)
!      lower(11,5) = (jacobian(5,11)-lower(5,1)*lower(11,1)-lower(5,2)*lower(11,2)-lower(5,3)*lower(11,3)-lower(5,4)*lower(11,4))/lower(5,5)
!      lower(12,5) = (jacobian(5,12)-lower(5,1)*lower(12,1)-lower(5,2)*lower(12,2)-lower(5,3)*lower(12,3)-lower(5,4)*lower(12,4))/lower(5,5)
!      lower(13,5) = (jacobian(5,13)-lower(5,1)*lower(13,1)-lower(5,2)*lower(13,2)-lower(5,3)*lower(13,3)-lower(5,4)*lower(13,4))/lower(5,5)
!!
!      lower(6,6) = SQRT(jacobian(6,6)-lower(6,1)**2-lower(6,2)**2 -lower(6,3)**2 - lower(6,4)**2- lower(6,5)**2)
!      lower(7,6) = (jacobian(6,7)-lower(6,1)*lower(7,1)-lower(6,2)*lower(7,2)-lower(6,3)*lower(7,3)-lower(6,4)*lower(7,4)-lower(6,5)*lower(7,5))/lower(6,6)
!      lower(8,6) = (jacobian(6,8)-lower(6,1)*lower(8,1)-lower(6,2)*lower(8,2)-lower(6,3)*lower(8,3)-lower(6,4)*lower(8,4)-lower(6,5)*lower(8,5))/lower(6,6)
!      lower(9,6) = (jacobian(6,9)-lower(6,1)*lower(9,1)-lower(6,2)*lower(9,2)-lower(6,3)*lower(9,3)-lower(6,4)*lower(9,4)-lower(6,5)*lower(9,5))/lower(6,6)
!      lower(10,6) = (jacobian(6,10)-lower(6,1)*lower(10,1)-lower(6,2)*lower(10,2)-lower(6,3)*lower(10,3)-lower(6,4)*lower(10,4)-lower(6,5)*lower(10,5))/lower(6,6)
!      lower(11,6) = (jacobian(6,11)-lower(6,1)*lower(11,1)-lower(6,2)*lower(11,2)-lower(6,3)*lower(11,3)-lower(6,4)*lower(11,4)-lower(6,5)*lower(11,5))/lower(6,6)
!      lower(12,6) = (jacobian(6,12)-lower(6,1)*lower(12,1)-lower(6,2)*lower(12,2)-lower(6,3)*lower(12,3)-lower(6,4)*lower(12,4)-lower(6,5)*lower(12,5))/lower(6,6)
!      lower(13,6) = (jacobian(6,13)-lower(6,1)*lower(13,1)-lower(6,2)*lower(13,2)-lower(6,3)*lower(13,3)-lower(6,4)*lower(13,4)-lower(6,5)*lower(13,5))/lower(6,6)
!!
!!
!!
!      lower(7,7) = SQRT(jacobian(7,7) - lower(7,1)**2 - lower(7,2)**2 - lower(7,3)**2 - lower(7,4)**2 - lower(7,5)**2 - lower(7,6)**2)
!!
!      lower(8,7) = (jacobian(7,8)-lower(7,1)*lower(8,1)-lower(7,2)*lower(8,2)-lower(7,3)*lower(8,3)-lower(7,4)*lower(8,4)&
!        -lower(7,5)*lower(8,5)-lower(7,6)*lower(8,6))/lower(7,7)
!!
!      lower(9,7) = (jacobian(7,9)-lower(7,1)*lower(9,1)-lower(7,2)*lower(9,2)-lower(7,3)*lower(9,3)-lower(7,4)*lower(9,4)&
!        -lower(7,5)*lower(9,5)-lower(7,6)*lower(9,6))/lower(7,7)
!!
!      lower(10,7) = (jacobian(7,10)-lower(7,1)*lower(10,1)-lower(7,2)*lower(10,2)-lower(7,3)*lower(10,3)-lower(7,4)*lower(10,4)&
!        -lower(7,5)*lower(10,5)-lower(7,6)*lower(10,6))/lower(7,7)
!!
!      lower(11,7) = (jacobian(7,11)-lower(7,1)*lower(11,1)-lower(7,2)*lower(11,2)-lower(7,3)*lower(11,3)-lower(7,4)*lower(11,4)&
!        -lower(7,5)*lower(11,5)-lower(7,6)*lower(11,6))/lower(7,7)
!!
!      lower(12,7) = (jacobian(7,12)-lower(7,1)*lower(12,1)-lower(7,2)*lower(12,2)-lower(7,3)*lower(12,3)-lower(7,4)*lower(12,4)&
!        -lower(7,5)*lower(12,5)-lower(7,6)*lower(12,6))/lower(7,7)
!!
!      lower(13,7) = (jacobian(7,13)-lower(7,1)*lower(13,1)-lower(7,2)*lower(13,2)-lower(7,3)*lower(13,3)-lower(7,4)*lower(13,4)&
!        -lower(7,5)*lower(13,5)-lower(7,6)*lower(13,6))/lower(7,7)
!!
!!
!!
!      lower(8,8) = SQRT(jacobian(8,8)-lower(8,1)**2-lower(8,2)**2 -lower(8,3)**2 - lower(8,4)**2- lower(8,5)**2 - lower(8,6)**2 - lower(8,7)**2)
!!
!      lower(9,8) = (jacobian(8,9) - lower(8,1)*lower(9,1) - lower(8,2)*lower(9,2)-lower(8,3)*lower(9,3)-lower(8,4)*lower(9,4)&
!        -lower(8,5)*lower(9,5) - lower(8,6)*lower(9,6) - lower(8,7)*lower(9,7))/lower(8,8)
!!
!      lower(10,8) = (jacobian(8,10)-lower(8,1)*lower(10,1)-lower(8,2)*lower(10,2)-lower(8,3)*lower(10,3)-lower(8,4)*lower(10,4)&
!        -lower(8,5)*lower(10,5) - lower(8,6)*lower(10,6) - lower(8,7)*lower(10,7))/lower(8,8)
!!
!      lower(11,8) = (jacobian(8,11)-lower(8,1)*lower(11,1)-lower(8,2)*lower(11,2)-lower(8,3)*lower(11,3)-lower(8,4)*lower(11,4)&
!        -lower(8,5)*lower(11,5) - lower(8,6)*lower(11,6) - lower(8,7)*lower(11,7))/lower(8,8)
!!
!      lower(12,8) = (jacobian(8,12)-lower(8,1)*lower(12,1)-lower(8,2)*lower(12,2)-lower(8,3)*lower(12,3)-lower(8,4)*lower(12,4)&
!        -lower(8,5)*lower(12,5) - lower(8,6)*lower(12,6) - lower(8,7)*lower(12,7))/lower(8,8)
!!
!      lower(13,8) = (jacobian(8,13)-lower(8,1)*lower(13,1)-lower(8,2)*lower(13,2)-lower(8,3)*lower(13,3)-lower(8,4)*lower(13,4)&
!        -lower(8,5)*lower(13,5) - lower(8,6)*lower(13,6) - lower(8,7)*lower(13,7))/lower(8,8)
!!
!!
!!
!      lower(9,9) = SQRT(jacobian(9,9) - lower(9,1)**2 - lower(9,2)**2 - lower(9,3)**2 - lower(9,4)**2- lower(9,5)**2 - lower(9,6)**2 &
!        - lower(9,7)**2 - lower(9,8)**2)
!!
!      lower(10,9) = (jacobian(9,10)-lower(9,1)*lower(10,1)-lower(9,2)*lower(10,2)-lower(9,3)*lower(10,3)-lower(9,4)*lower(10,4)&
!        -lower(9,5)*lower(10,5) - lower(9,6)*lower(10,6) - lower(9,7)*lower(10,7)- lower(9,8)*lower(10,8))/lower(9,9)
!!
!      lower(11,9) = (jacobian(9,11)-lower(9,1)*lower(11,1)-lower(9,2)*lower(11,2)-lower(9,3)*lower(11,3)-lower(9,4)*lower(11,4)&
!        -lower(9,5)*lower(11,5) - lower(9,6)*lower(11,6) - lower(9,7)*lower(11,7)- lower(9,8)*lower(11,8))/lower(9,9)
!!
!      lower(12,9) = (jacobian(9,12)-lower(9,1)*lower(12,1)-lower(9,2)*lower(12,2)- lower(9,3)*lower(12,3)-lower(9,4)*lower(12,4)&
!        -lower(9,5)*lower(12,5) - lower(9,6)*lower(12,6) - lower(9,7)*lower(12,7)- lower(9,8)*lower(12,8))/lower(9,9)
!!
!      lower(13,9) = (jacobian(9,13)-lower(9,1)*lower(13,1)-lower(9,2)*lower(13,2)- lower(9,3)*lower(13,3)-lower(9,4)*lower(13,4)&
!        -lower(9,5)*lower(13,5) - lower(9,6)*lower(13,6) - lower(9,7)*lower(13,7)- lower(9,8)*lower(13,8))/lower(9,9)
!!
!!
!!
!      lower(10,10) = SQRT(jacobian(10,10) - lower(10,1)**2 - lower(10,2)**2 - lower(10,3)**2 - lower(10,4)**2 - lower(10,5)**2 - lower(10,6)**2 &
!        - lower(10,7)**2 - lower(10,8)**2 - lower(10,9)**2)
!!
!      lower(11,10) = (jacobian(10,11)-lower(10,1)*lower(11,1)-lower(10,2)*lower(11,2)-lower(10,3)*lower(11,3)-lower(10,4)*lower(11,4)-lower(10,5)*lower(11,5) &
!        - lower(10,6)*lower(11,6) - lower(10,7)*lower(11,7) - lower(10,8)*lower(11,8)- lower(10,9)*lower(11,9))/lower(10,10)
!!
!      lower(12,10) =(jacobian(10,12)-lower(10,1)*lower(12,1)-lower(10,2)*lower(12,2)-lower(10,3)*lower(12,3)-lower(10,4)*lower(12,4)-lower(10,5)*lower(12,5) &
!        - lower(10,6)*lower(12,6) - lower(10,7)*lower(12,7) - lower(10,8)*lower(12,8)- lower(10,9)*lower(12,9))/lower(10,10)
!!
!      lower(13,10) =(jacobian(10,13)-lower(10,1)*lower(13,1)-lower(10,2)*lower(13,2)-lower(10,3)*lower(13,3)-lower(10,4)*lower(13,4)-lower(10,5)*lower(13,5) &
!        - lower(10,6)*lower(13,6) - lower(10,7)*lower(13,7) - lower(10,8)*lower(13,8)- lower(10,9)*lower(13,9))/lower(10,10)
!!
!!
!!
!      lower(11,11) = SQRT(jacobian(11,11)-lower(11,1)**2-lower(11,2)**2 -lower(11,3)**2 - lower(11,4)**2- lower(11,5)**2 - lower(11,6)**2 - lower(11,7)**2 &
!        - lower(11,8)**2 - lower(11,9)**2 - lower(11,10)**2)
!!
!      lower(12,11) = (jacobian(11,12)-lower(11,1)*lower(12,1)-lower(11,2)*lower(12,2)-lower(11,3)*lower(12,3)-lower(11,4)*lower(12,4)-lower(11,5)*lower(12,5) &
!        - lower(11,6)*lower(12,6) - lower(11,7)*lower(12,7) - lower(11,8)*lower(12,8)- lower(11,9)*lower(12,9)- lower(11,10)*lower(12,10))/lower(11,11)
!!
!      lower(13,11) = (jacobian(11,13)-lower(11,1)*lower(13,1)-lower(11,2)*lower(13,2)-lower(11,3)*lower(13,3)-lower(11,4)*lower(13,4)-lower(11,5)*lower(13,5) &
!        - lower(11,6)*lower(13,6) - lower(11,7)*lower(13,7) - lower(11,8)*lower(13,8)- lower(11,9)*lower(13,9)- lower(11,10)*lower(13,10))/lower(11,11)
!!
!!
!!
!      lower(12,12) = SQRT(jacobian(12,12)-lower(12,1)**2-lower(12,2)**2 -lower(12,3)**2 - lower(12,4)**2- lower(12,5)**2 - lower(12,6)**2 - lower(12,7)**2 &
!         - lower(12,8)**2 - lower(12,9)**2 - lower(12,10)**2 - lower(12,11)**2)
!!
!      lower(13,12) =(jacobian(12,13)-lower(12,1)*lower(13,1)-lower(12,2)*lower(13,2)-lower(12,3)*lower(13,3)-lower(12,4)*lower(13,4)-lower(12,5)*lower(13,5) &
!        - lower(12,6)*lower(13,6) - lower(12,7)*lower(13,7) - lower(12,8)*lower(13,8)- lower(12,9)*lower(13,9)- lower(12,10)*lower(13,10)&
!        - lower(12,11)*lower(13,11))/lower(12,12)
!!
!!
!!
!      lower(13,13) = SQRT(jacobian(13,13) - lower(13,1)**2 - lower(13,2)**2 - lower(13,3)**2 - lower(13,4)**2- lower(13,5)**2 - &
!        lower(13,6)**2 - lower(13,7)**2 - lower(13,8)**2 - lower(13,9)**2 - lower(13,10)**2 - lower(13,11)**2 - lower(13,12)**2)





!      IF (dimensions == 2) THEN


!      ELSE IF (dimensions == 3) THEN
!
! Calculate the majority of the L and U matrices. Use the Jacobian
! values for the top row of the U matrix.
!
!        l_21 = jacobian(2,1)/jacobian(1,1)
!        u_22 = jacobian(2,2)-l_21*jacobian(1,2)
!        u_23 = jacobian(2,3)-l_21*jacobian(1,3)
!        u_24 = jacobian(2,4)-l_21*jacobian(1,4)
!        u_25 = jacobian(2,5)-l_21*jacobian(1,5)
!
!        l_31 = jacobian(3,1)/jacobian(1,1)
!        l_32 = (jacobian(3,2) - l_31*jacobian(1,3))/u_22
!        u_33 = jacobian(3,3) - l_31*jacobian(1,3) - l_32*u_23
!        u_34 = jacobian(3,4) - l_31*jacobina(1,4) - l_32*u_24
!        u_35 = jacobian(3,5) - l_31*jacobina(1,5) - l_32*u_25
!
!        l_41 = jacobian(4,1)/jacobian(1,1)
!        l_42 = (jacobian(4,2)-l_41*jacobian(1,2))/u_22
!        l_43 = (jacobian(4,3)-l_41*jacobian(1,3) - l_42*u_23)/u_33
!        u_44 = jacobian(4,4)-l_41*jacobian(1,4) - l_42*u_24 - l_43*u_34
!        u_45 = jacobian(4,5)-l_41*jacobian(1,5) - l_42*u_25 - l_43*u_35
!
!        l_51 = jacobina(5,1)/jacobian(1,1)
!        l_52 = (jacobian(5,2) - l_51*jacobian(1,2))/u_22
!        l_53 = (jacobian(5,3) - l_51*jacobian(1,3) - l_52*u_23)/u_33
!        l_54 = (jacobian(5,4) - l_51*jacobian(1,4) - l_52*u_24 - l_53*u_34)/u_44
!        u_55 = jacobian(5,5) - l_51*jacobian(1,5) - l_52*u_25 - l_53*u_35 - l_54*u_45
!      END IF
