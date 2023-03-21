module inlet_outlet_elbm_data
!
! This stores values for all inlets and outlets in the simulations.
! fieq and gieq are stored as no information is known exterior to the domain,
! so equilibrium values at inlet/outlet conditions are used.  This will primarily
! be used in the boundary condition subroutines.
!
!
!
!
use precise
use nml_inlet_outlet
implicit none
save

private
public :: inl_chi,inl_zeta_x,inl_zeta_y,inl_zeta_z,&
  inl_pi_xx,inl_pi_yy,inl_pi_zz,inl_pi_xy,inl_pi_xz,inl_pi_yz,&
  inl_lambda_x,inl_lambda_y,inl_lambda_z
public :: outl_chi,outl_zeta_x,outl_zeta_y,outl_zeta_z,&
  outl_pi_xx,outl_pi_yy,outl_pi_zz,outl_pi_xy,outl_pi_xz,outl_pi_yz,&
  outl_lambda_x,outl_lambda_y,outl_lambda_z
public :: inl_fieq,inl_gieq,inl_beta1,outl_beta1,inl_beta2,outl_beta2
public :: inl_rho,inl_u,inl_v,inl_w,inl_temp,inl_press
public :: outl_rho,outl_u,outl_v,outl_w,outl_temp,outl_press
public :: outl_fieq,outl_gieq
public :: inl_alpha,outl_alpha
public :: inl_mu,inl_kappa,outl_mu,outl_kappa
!
!
!
!
!
real(kind=dp),allocatable :: inl_chi(:),inl_zeta_x(:),inl_zeta_y(:),inl_zeta_z(:),&
  inl_pi_xx(:),inl_pi_yy(:),inl_pi_zz(:),inl_pi_xy(:),inl_pi_xz(:),inl_pi_yz(:),&
  inl_lambda_x(:),inl_lambda_y(:),inl_lambda_z(:)
real(kind=dp),allocatable :: outl_chi(:),outl_zeta_x(:),outl_zeta_y(:),outl_zeta_z(:),&
  outl_pi_xx(:),outl_pi_yy(:),outl_pi_zz(:),outl_pi_xy(:),outl_pi_xz(:),outl_pi_yz(:),&
  outl_lambda_x(:),outl_lambda_y(:),outl_lambda_z(:)

real(kind=dp),allocatable :: inl_fieq(:,:),inl_gieq(:,:)
real(kind=dp),allocatable :: outl_fieq(:,:),outl_gieq(:,:)

real(kind=dp),allocatable :: inl_beta1(:),inl_beta2(:)
real(kind=dp),allocatable :: outl_beta1(:),outl_beta2(:)

real(kind=dp),allocatable :: inl_alpha(:),outl_alpha(:)
real(kind=dp),allocatable :: inl_mu(:),inl_kappa(:),outl_mu(:),outl_kappa(:)

real(kind=dp),allocatable :: inl_rho(:),inl_u(:),inl_v(:),inl_w(:),&
  inl_temp(:),inl_press(:)
real(kind=dp),allocatable :: outl_rho(:),outl_u(:),outl_v(:),outl_w(:),&
  outl_temp(:),outl_press(:)


end module
