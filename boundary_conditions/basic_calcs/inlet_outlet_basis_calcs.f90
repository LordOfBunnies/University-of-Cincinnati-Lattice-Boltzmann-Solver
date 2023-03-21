subroutine inlet_outlet_basis_calcs
!
! Find all the ELBM values necessary to probably solve the boundary conditions.
! This does not include interpolative bounceback, that is handled elsewhere.
!
! Called by: background_calcs
! Calls: quick_sutherlands,quick_therm_cond,calculate_lagrange_mults,
!    calculate_alpha,inlet_outlet_base_values
! External calls:
!
use precise
use linkwise
use constants
use grid_data
use nml_inlet_outlet
use freestream_values
use inlet_outlet_elbm_data
implicit none
integer :: i,j,a,b,c
real(kind=dp),allocatable :: loc_fieq(:),loc_gieq(:)
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_press
real(kind=dp) :: loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz
logical :: is_it_inlet

a = 0
b = 0
c = 0
!
! Inlets
!
!write(*,*) 'fuzzy winker doodle'

! Allocations first
!
allocate(inl_chi(num_inlets))
allocate(inl_zeta_x(num_inlets))
allocate(inl_zeta_y(num_inlets))
allocate(inl_pi_xx(num_inlets))
allocate(inl_pi_yy(num_inlets))
allocate(inl_pi_xy(num_inlets))
allocate(inl_lambda_x(num_inlets))
allocate(inl_lambda_y(num_inlets))
! 3D
if (dimensions == 3) then
  allocate(inl_zeta_z(num_inlets))
  allocate(inl_pi_zz(num_inlets))
  allocate(inl_pi_xz(num_inlets))
  allocate(inl_pi_yz(num_inlets))
  allocate(inl_lambda_z(num_inlets))

  inl_zeta_z = 0.0D0
  inl_pi_zz = -0.3D0
  inl_pi_xz = 0.0D0
  inl_pi_yz = 0.0D0
  inl_lambda_z = 0.0D0
end if
! Everything else, may not be needed
allocate(inl_fieq(dir,num_inlets))
allocate(inl_gieq(dir,num_inlets))
!
allocate(inl_alpha(num_inlets))
allocate(inl_mu(num_inlets))
allocate(inl_kappa(num_inlets))
allocate(inl_beta1(num_inlets))
allocate(inl_beta2(num_inlets))
!
allocate(inl_rho(num_inlets))
allocate(inl_u(num_inlets))
allocate(inl_v(num_inlets))
allocate(inl_w(num_inlets))
allocate(inl_temp(num_inlets))
allocate(inl_press(num_inlets))
!
allocate(loc_fieq(dir))
allocate(loc_gieq(dir))
!write(*,*) 'boot scoot'
inl_chi = -2.5D0
inl_zeta_x = 0.0D0
inl_zeta_y = 0.0D0
inl_pi_xx = -0.3D0
inl_pi_yy = -0.3D0
inl_pi_xy = 0.0D0
inl_lambda_x = 0.0D0
inl_lambda_y = 0.0D0
inl_alpha = 2.0D0

!
! Find simulations variables for inlets
!
do i = 1,num_inlets
!
!  call inlet_outlet_base_values(i,loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_press)
!


  if (inlet_type(1,i) == 1 .or. inlet_type(1,i) == 2 .or. inlet_type(1,i) == 4 .or. &
      inlet_type(1,i) ==7) then
    is_it_inlet = .true.
    ! XXXX This needs to be switched back on once actual iterations are needed
    call inlet_outlet_base_values(i,is_it_inlet,loc_rho,&
      loc_u,loc_v,loc_w,loc_temp,loc_press)
!
    loc_pxx = loc_rho*(loc_temp+loc_u**2)
    loc_pyy = loc_rho*(loc_temp+loc_v**2)
    loc_pzz = loc_rho*(loc_temp+loc_w**2)
    loc_pxy = loc_rho*loc_u*loc_v
    loc_pxz = loc_rho*loc_u*loc_w
    loc_pyz = loc_rho*loc_v*loc_w

    write(*,*) 'inlet values',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_press,non_dim_R
!

    call calculate_lagrangian_mults(loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,&
      loc_pxy,loc_pxz,loc_pyz,&
      inl_chi(i),inl_zeta_x(i),inl_zeta_y(i),inl_zeta_z(i),inl_pi_xx(i),inl_pi_yy(i),&
      inl_pi_zz(i),inl_pi_xy(i),inl_pi_xz(i),inl_pi_yz(i),&
      inl_lambda_x(i),inl_lambda_y(i),inl_lambda_z(i),inl_fieq(1:dir,i),a,b,c)

    inl_gieq = inl_fieq*(2*const_vol-dimensions)*loc_temp
    !call calculate_alpha(loc_fieq,loc_fieq,loc_fieq,inl_alpha(i))

    write(11,*) 'Inlet ',i,' Lagrange multiplier values',loc_rho,&
      loc_u,loc_v,loc_w,loc_temp,loc_press,inl_chi(i),inl_zeta_x(i),inl_zeta_y(i),&
      inl_zeta_z(i),inl_pi_xx(i),inl_pi_yy(i),inl_pi_zz(i),inl_pi_xy(i),inl_pi_xz(i),inl_pi_yz(i),&
      inl_lambda_x(i),inl_lambda_y(i),inl_lambda_z(i)
    write(11,*) 'Inlet',i,' alpha value',inl_alpha(i)
    write(11,*) 'Inlet',i,' fieq values',inl_fieq(1:dir,i)
    write(11,*) 'Inlet',i,' gieq values',inl_gieq(1:dir,i)
    write(*,*) 'Inlet ',i,' Lagrange multiplier values',loc_rho,&
      loc_u,loc_v,loc_w,loc_temp,loc_press,inl_chi(i),inl_zeta_x(i),inl_zeta_y(i),&
      inl_zeta_z(i),inl_pi_xx(i),inl_pi_yy(i),inl_pi_zz(i),inl_pi_xy(i),inl_pi_xz(i),inl_pi_yz(i),&
      inl_lambda_x(i),inl_lambda_y(i),inl_lambda_z(i)
!    write(*,*) 'Inlet',i,' alpha value',inl_alpha(i)
    write(*,*) 'Inlet',i,' fieq values',inl_fieq(1:dir,i)
    write(*,*) 'Inlet',i,' gieq values',inl_gieq(1:dir,i)
!  else if (inlet_type(1,i) == 2) then
!    call inlet_outlet_base_values(i,is_it_inlet,loc_rho,&
!      loc_u,loc_v,loc_w,loc_temp,loc_press)
!
  else if (inlet_type(1,i) == 3) then

    inl_chi(i) = fs_chi
    inl_zeta_x(i) = fs_zeta_x
    inl_zeta_y(i) = fs_zeta_y
    inl_pi_xx(i) = fs_pixx
    inl_pi_yy(i) = fs_piyy
    inl_pi_xy(i) = fs_pixy
    inl_lambda_x(i) = fs_lambda_x
    inl_lambda_y(i) = fs_lambda_y

    inl_mu(i) = fs_mu
    inl_kappa(i) = fs_kappa
    inl_alpha(i) = fs_alpha

    if (dimensions == 3) then
      inl_pi_zz(i) = fs_pizz
      inl_pi_xz(i) = fs_pixz
      inl_pi_yz(i) = fs_piyz
      inl_zeta_z(i) = fs_zeta_z
      inl_lambda_z(i) = fs_lambda_z
    end if

    inl_fieq(1:dir,i) = fs_fieq(1:dir)
    inl_gieq(1:dir,i) = fs_gieq(1:dir)
  end if
  !write(*,*) 'whoopsie'
!
!
end do

!
! Outlets
!
! Allocations first
!
!
allocate(outl_chi(num_outlets))
allocate(outl_zeta_x(num_outlets))
allocate(outl_zeta_y(num_outlets))
allocate(outl_pi_xx(num_outlets))
allocate(outl_pi_yy(num_outlets))
allocate(outl_pi_xy(num_outlets))
allocate(outl_lambda_x(num_outlets))
allocate(outl_lambda_y(num_outlets))
! 3D
if (dimensions == 3) then
  allocate(outl_zeta_z(num_outlets))
  allocate(outl_pi_zz(num_outlets))
  allocate(outl_pi_xz(num_outlets))
  allocate(outl_pi_yz(num_outlets))
  allocate(outl_lambda_z(num_outlets))

  outl_zeta_z = 0.0D0
  outl_pi_zz = -0.3D0
  outl_pi_xz = 0.0D0
  outl_pi_yz = 0.0D0
  outl_lambda_z = 0.0D0
end if
! Everything else, may not be needed
allocate(outl_fieq(dir,num_outlets))
allocate(outl_gieq(dir,num_outlets))
!
allocate(outl_alpha(num_outlets))
allocate(outl_mu(num_outlets))
allocate(outl_kappa(num_outlets))
allocate(outl_beta1(num_outlets))
allocate(outl_beta2(num_outlets))
!
allocate(outl_rho(num_outlets))
allocate(outl_u(num_outlets))
allocate(outl_v(num_outlets))
allocate(outl_w(num_outlets))
allocate(outl_temp(num_outlets))
allocate(outl_press(num_outlets))
!
outl_chi = -2.5D0
outl_zeta_x = 0.0D0
outl_zeta_y = 0.0D0
outl_pi_xx = -0.3D0
outl_pi_yy = -0.3D0
outl_pi_xy =  0.0D0
outl_lambda_x = 0.0D0
outl_lambda_y = 0.0D0
outl_alpha = 2.0D0
!
! Find all required conditions for the outlets
!
do i = 1,num_outlets

  if (outlet_type(1,i) == 1 .or. outlet_type(1,i) == 2 .or.  outlet_type(1,i) == 4 .or. &
      outlet_type(1,i) == 7) then
    is_it_inlet = .false.

!    write(*,*) 'wibbles!', i,outlet_flow_val(1,i),outlet_flow_val(2,i),outlet_flow_val(3,i),&
!      outlet_flow_val(4,i),outlet_flow_val(5,i)

    ! XXXX This needs to be switched back on once actual iterations are needed
    call inlet_outlet_base_values(i,is_it_inlet,loc_rho,&
      loc_u,loc_v,loc_w,loc_temp,loc_press)
!
    loc_pxx = loc_rho*(loc_temp+loc_u**2)
    loc_pyy = loc_rho*(loc_temp+loc_v**2)
    loc_pzz = loc_rho*(loc_temp+loc_w**2)
    loc_pxy = loc_rho*loc_u*loc_v
    loc_pxz = loc_rho*loc_u*loc_w
    loc_pyz = loc_rho*loc_v*loc_w

!    write(*,*) 'outlet values',loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_press,loc_pxx,&
!      loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,non_dim_R

    if (shifted) then
    outl_chi(i) = -2.5D0
    outl_zeta_x(i) = loc_u
    outl_zeta_y(i) = loc_v
    outl_zeta_z(i) = loc_w
    outl_pi_xx(i) = -0.3D0
    outl_pi_yy(i) = -0.3D0
    outl_pi_zz(i) = -0.3D0
    outl_pi_xy(i) = 0.0D0
    outl_pi_xz(i) = 0.0D0
    outl_pi_yz(i) = 0.0D0
    outl_lambda_x(i) = loc_u/25
    outl_lambda_y(i) = loc_v/25
    outl_lambda_z(i) = loc_w/25


    end if
!
    call calculate_lagrangian_mults(loc_rho,loc_u,loc_v,loc_w,loc_temp,loc_pxx,loc_pyy,loc_pzz,&
      loc_pxy,loc_pxz,loc_pyz,outl_chi(i),outl_zeta_x(i),outl_zeta_y(i),&
      outl_zeta_z(i),outl_pi_xx(i),outl_pi_yy(i),outl_pi_zz(i),outl_pi_xy(i),outl_pi_xz(i),outl_pi_yz(i),&
      outl_lambda_x(i),outl_lambda_y(i),outl_lambda_z(i),outl_fieq(1:dir,i),a,b,c)
!    call calculate_alpha(loc_fieq,loc_fieq,loc_fieq,outl_alpha(i))
    outl_gieq = outl_fieq*(2*const_vol-dimensions)*loc_temp

    write(11,*) 'Outlet ',i,' Lagrange multiplier values',loc_rho,&
      loc_u,loc_v,loc_w,loc_temp,loc_press,outl_chi(i),outl_zeta_x(i),outl_zeta_y(i),&
      outl_zeta_z(i),outl_pi_xx(i),outl_pi_yy(i),outl_pi_zz(i),outl_pi_xy(i),outl_pi_xz(i),outl_pi_yz(i),&
      outl_lambda_x(i),outl_lambda_y(i),outl_lambda_z(i)
    write(11,*) 'Outlet',i,' alpha value',outl_alpha(i)
    write(11,*) 'Outlet',i,' fieq values',outl_fieq(1:dir,i)
    write(11,*) 'Outlet',i,' gieq values',outl_gieq(1:dir,i)
    write(*,*) 'Outlet ',i,' Lagrange multiplier values',loc_rho,&
      loc_u,loc_v,loc_w,loc_temp,loc_press,outl_chi(i),outl_zeta_x(i),outl_zeta_y(i),&
      outl_zeta_z(i),outl_pi_xx(i),outl_pi_yy(i),outl_pi_zz(i),outl_pi_xy(i),outl_pi_xz(i),outl_pi_yz(i),&
      outl_lambda_x(i),outl_lambda_y(i),outl_lambda_z(i)
    write(*,*) 'Outlet',i,' alpha value',outl_alpha(i)
    write(*,*) 'Outlet',i,' fieq values',outl_fieq(1:dir,i)
    write(*,*) 'Outlet',i,' gieq values',outl_gieq(1:dir,i)


!  else if (outlet_type(1,i) == 2) then
!    call inlet_outlet_base_values(i,is_it_inlet,loc_rho,&
!      loc_u,loc_v,loc_w,loc_temp,loc_press)

  else if (outlet_type(1,i) == 3) then
    outl_chi(i) = fs_chi
    outl_zeta_x(i) = fs_zeta_x
    outl_zeta_y(i) = fs_zeta_y
    outl_pi_xx(i) = fs_pixx
    outl_pi_yy(i) = fs_piyy
    outl_pi_xy(i) = fs_pixy
    outl_lambda_x(i) = fs_lambda_x
    outl_lambda_y(i) = fs_lambda_y

    outl_mu(i) = fs_mu
    outl_kappa(i) = fs_kappa
    outl_alpha(i) = fs_alpha

    if (dimensions == 3) then
      outl_pi_zz(i) = fs_pizz
      outl_pi_xz(i) = fs_pixz
      outl_pi_yz(i) = fs_piyz
      outl_zeta_z(i) = fs_zeta_z
      outl_lambda_z(i) = fs_lambda_z
    end if

    outl_fieq(1:dir,i) = fs_fieq(1:dir)
    outl_gieq(1:dir,i) = fs_gieq(1:dir)
  end if


end do

!write(*,*) 'scruff says...'

end subroutine
