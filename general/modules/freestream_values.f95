MODULE freestream_values
!
! This saves all the freestream values so they can be called in when
! needed.  Probably much better than trying to pass everything
! everywhere.
!
!
!
!
USE precise
IMPLICIT NONE
SAVE
!
! Constants
!
REAL(KIND=dp) :: temperature,pressure,v_inf,rho_inf,iso_therm_sos,mu_inf,kappa_inf
REAL(KIND=dp) :: stag_temperature,stag_pressure,stag_density,viscosity,nu
real(kind=dp) :: characteristic_time
real(kind=dp),allocatable :: lattice_reynolds(:)
CHARACTER(len=80) :: gas, turb_type
!
! LMs
!
real(kind=dp) :: fs_u,fs_v,fs_w,fs_vmag
real(kind=dp) :: fs_alpha,fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
  fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z
real(kind=dp) :: fs_pxx,fs_pyy,fs_pzz,fs_pxy,fs_pxz,fs_pyz,fs_qx,fs_qy,fs_qz
! random necessities
real(kind=dp) :: fs_beta_1,fs_beta_2,fs_mu,fs_kappa
!
! distribution functions
!
real(kind=dp),allocatable :: fs_fieq(:),fs_gieq(:)
!
!
!
real(kind=dp) :: su_chi,su_zeta_x,su_zeta_y,su_zeta_z,su_pixx,su_piyy,su_pizz,&
  su_pixy,su_pixz,su_piyz,su_lambda_x,su_lambda_y,su_lambda_z


END MODULE
