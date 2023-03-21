MODULE constants
!
! This stores constants of use later. These include defined constants
! from the user like tau and gama (gamma, but that's used in Fortran)
!
!
!
!
USE precise
IMPLICIT NONE
SAVE
REAL(KIND = dp) :: gas_const_R,prandtl,gama,reynolds,mach,speed_of_sound
real(kind = dp) :: cp_unit_air,cv_unit_air
!REAL(KIND = dp) :: tau_lower_limit = 0.51,tau_upper_limit = 20.0
REAL(KIND = dp) :: kb
REAL(KIND = dp),PARAMETER :: molec_diam_air = 3.579D-10,pi = 3.1415926535D0,&
  molec_mass_air = 28.97,avocados_number = 6.022D23,default_char_leng = 1.0D0
real(kind = dp),parameter :: k_boltz = 1.38064852D-23,kappa_convert=418.1D0!,prandtl = 0.71D0
real(kind = dp) :: molec_weight
REAL(KIND = dp) :: t_ref,rho_ref,p_ref,kappa_ref,mu_ref,non_dim_R,u_ref !dimensional reference values
REAL(KIND = dp) :: const_vol,const_press
real(kind=dp) :: rest_chi,rest_zeta_x,rest_zeta_y,rest_zeta_z,rest_pixx,rest_pixy,&
  rest_piyy,rest_pizz,rest_pixz,rest_piyz,rest_lambda_x,rest_lambda_y,rest_lambda_z
real(kind=dp) :: rest_press
!REAL

INTEGER :: dimensions,primary_flow_dir,cells_per_length

LOGICAL :: viscous,turbulent,bomb_file,body_forces,external_flow
!      LOGICAL :: energy


END MODULE
