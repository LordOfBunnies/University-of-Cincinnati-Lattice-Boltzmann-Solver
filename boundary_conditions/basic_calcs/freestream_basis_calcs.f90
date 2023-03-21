subroutine freestream_basis_calcs
!
! Create basic freestream values for the more complicated ELBM procedures
!
!
!
! Called by: backgrounds_calcs
! Calls: quick_sutherlands,quick_therm_cond, calculate_lagrange_mults,calculate_alpha
!
use precise
use constants
use linkwise
use freestream_values
use grid_data, only: dir
use amr_processes, only: self
use quick_calcs
implicit none
integer :: i,a,b,c
real(kind=dp) :: loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pyz,loc_pxz
real(kind=dp) :: loc_fieq(dir),su_fieq(dir)
!real(kind=dp) :: loc_rho,loc_temp
real(kind=dp) :: loc_u,loc_v,loc_w
real(kind=dp) :: su_rho,su_u,su_v,su_w,su_temp,su_pxx,su_pyy,su_pzz,su_pxy,&
  su_pxz,su_pyz,su_qx,su_qy,su_qz,su_vmag

allocate(fs_fieq(dir))
allocate(fs_gieq(dir))

a = 0
b = 0
c = 0

fs_vmag = fs_u**2 + fs_v**2 + fs_w**2

if (.not. shifted) then

  fs_chi=-2.4D0
  fs_zeta_x=0.0D0
  fs_zeta_y=0.0D0
  fs_zeta_z=0.0D0
  fs_pixx=-0.40D0
  fs_piyy=-0.40D0
  fs_pizz=-0.40D0
  fs_pixy=0.0D0
  fs_pixz=0.0D0
  fs_piyz=0.0D0
  fs_lambda_x=0.0D0
  fs_lambda_y=0.0D0
  fs_lambda_z=0.0D0
  !initial alpha
  fs_alpha = 2.0D0

  fs_pxx = rho_inf*(temperature+fs_u**2)
  fs_pyy = rho_inf*(temperature+fs_v**2)
  fs_pzz = rho_inf*(temperature+fs_w**2)
  fs_pxy = rho_inf*fs_u*fs_v
  fs_pxz = rho_inf*fs_u*fs_w
  fs_pyz = rho_inf*fs_v*fs_w

  fs_qx = (rho_inf*fs_vmag + 5*rho_inf*temperature)*fs_u
  fs_qy = (rho_inf*fs_vmag + 5*rho_inf*temperature)*fs_v
  fs_qz = (rho_inf*fs_vmag + 5*rho_inf*temperature)*fs_w


  write(*,*) 'freestream stuff',rho_inf,fs_u,fs_v,fs_w,temperature
  call calculate_lagrangian_mults(rho_inf,fs_u,fs_v,fs_w,temperature,&
    fs_pxx,fs_pyy,fs_pzz,fs_pxy,fs_pxz,fs_pyz,&
    fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
    fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z,fs_fieq,a,b,c)

    fs_gieq = fs_fieq*(2*const_vol-dimensions)*temperature
!
! For supersonic flows, you may need to iterate a couple times to get freestream conditions
!
else
!
!
!
  fs_chi=-3.0D0
  fs_zeta_x=0.0D0
  fs_zeta_y=0.0D0
  fs_zeta_z=0.0D0
  fs_pixx=-0.40D0
  fs_piyy=-0.40D0
  fs_pizz=-0.40D0
  fs_pixy=0.0D0
  fs_pixz=0.0D0
  fs_piyz=0.0D0
  fs_lambda_x=0.0D0
  fs_lambda_y=0.0D0
  fs_lambda_z=0.0D0
  !initial alpha
  fs_alpha = 2.0D0

  select case(primary_flow_dir)
    case(1)
      fs_zeta_x = mach
      fs_lambda_x = -mach/50.0D0
    case(2)
      fs_zeta_y = mach
      fs_lambda_y = -mach/50.0D0
    case(3)
      fs_zeta_z = mach
      fs_lambda_z = -mach/50.0D0
    case(4)
      fs_zeta_x = -mach
      fs_lambda_x = mach/50.0D0
    case(5)
      fs_zeta_y = -mach
      fs_lambda_y = mach/50.0D0
    case(6)
      fs_zeta_z = -mach
      fs_lambda_z = mach/50.0D0
  end select

!  do i = 1,4

!    loc_u = real(i,kind=dp)/4.0D0*fs_u
!    loc_v = real(i,kind=dp)/4.0D0*fs_v
!    loc_w = real(i,kind=dp)/4.0D0*fs_w


    fs_pxx = rho_inf*(temperature+fs_u**2)
    fs_pyy = rho_inf*(temperature+fs_v**2)
    fs_pzz = rho_inf*(temperature+fs_w**2)
    fs_pxy = rho_inf*fs_u*fs_v
    fs_pxz = rho_inf*fs_u*fs_w
    fs_pyz = rho_inf*fs_v*fs_w

    fs_qx = (rho_inf*fs_vmag + 5*rho_inf*temperature)*fs_u
    fs_qy = (rho_inf*fs_vmag + 5*rho_inf*temperature)*fs_v
    fs_qz = (rho_inf*fs_vmag + 5*rho_inf*temperature)*fs_w


    write(*,*) 'freestream stuff',rho_inf,fs_u,fs_v,fs_w,temperature
    call calculate_lagrangian_mults(rho_inf,fs_u,fs_v,fs_w,temperature,&
      fs_pxx,fs_pyy,fs_pzz,fs_pxy,fs_pxz,fs_pyz,&
      fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
      fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z,fs_fieq,a,b,c)

      fs_gieq = fs_fieq*(2*const_vol-dimensions)*temperature

    if (self == 0) then
      write(*,*) 'Freestream Langrange multiplier values',fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
        fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z
      write(*,*) 'Freestream alpha value',fs_alpha
      write(*,*) 'freestream fieq values',fs_fieq
      write(*,*) 'freestream gieq values',fs_gieq
    end if

    a = 0
    b = 0
    c = 0

    if (mach > 1.0D0) then
      su_chi = -2.2D0
      su_zeta_X = 0.4D0
      su_zeta_y = 0.0D0
      su_zeta_z = 0.0D0
      su_pixx = -0.45D0
      su_piyy = -0.45D0
      su_pizz = -0.45D0
      su_pixy = 0.0D0
      su_pixz = 0.0D0
      su_piyz = 0.0D0
      su_lambda_x = 0.015D0
      su_lambda_y = 0.0D0
      su_lambda_z = 0.0D0

      su_rho = 1.2D0
      su_u = 0.80D0
      su_v = 0.0D0
      su_w = 0.0D0
      su_temp = 0.75D0
    else

      su_chi = -2.0D0
      su_zeta_X = 0.0D0
      su_zeta_y = 0.0D0
      su_zeta_z = 0.0D0
      su_pixx = -0.4D0
      su_piyy = -0.4D0
      su_pizz = -0.4D0
      su_pixy = 0.0D0
      su_pixz = 0.0D0
      su_piyz = 0.0D0
      su_lambda_x = 0.002D0
      su_lambda_y = 0.0D0
      su_lambda_z = 0.0D0

      su_rho = 1.2D0
      su_u = 0.4D0
      su_v = 0.0D0
      su_w = 0.0D0
      su_temp = 0.714D0

    end if

    su_vmag = su_u**2 + su_v**2 + su_w**2

    su_pxx = su_rho*(su_temp + su_u**2)
    su_pyy = su_rho*(su_temp + su_v**2)
    su_pzz = su_rho*(su_temp + su_w**2)
    su_pxy = su_rho*su_u*su_v
    su_pxz = su_rho*su_u*su_w
    su_pyz = su_rho*su_v*su_w

    su_qx = (su_rho*su_vmag + 5*su_rho*su_temp)*su_u
    su_qy = (su_rho*su_vmag + 5*su_rho*su_temp)*su_v
    su_qz = (su_rho*su_vmag + 5*su_rho*su_temp)*su_w

    write(*,*) 'startup stuff',su_rho,su_u,su_v,su_w,su_temp
    call calculate_lagrangian_mults(su_rho,su_u,su_v,su_w,su_temp,&
      su_pxx,su_pyy,su_pzz,su_pxy,su_pxz,su_pyz,&
      su_chi,su_zeta_x,su_zeta_y,su_zeta_z,su_pixx,su_piyy,&
      su_pizz,su_pixy,su_pixz,su_piyz,su_lambda_x,su_lambda_y,su_lambda_z,su_fieq,a,b,c)

    if (self == 0) then
      write(*,*) 'Startup Langrange multiplier values',su_chi,su_zeta_x,su_zeta_y,su_zeta_z,&
        su_pixx,su_piyy,su_pizz,su_pixy,su_pixz,su_piyz,su_lambda_x,su_lambda_y,su_lambda_z
    end if
!  end do

end if

!
! For testing purposes
!
!fs_u = 0.72D0
!fs_v = -0.18D0
!fs_w = 0.2D0
!loc_rho = .9
!loc_temp = .69
!loc_pxx = loc_rho*(loc_temp+fs_u**2)
!loc_pyy = loc_rho*(loc_temp+fs_v**2)
!loc_pzz = loc_rho*(loc_temp+fs_w**2)
!loc_pxy = loc_rho*fs_u*fs_v
!loc_pxz = loc_rho*fs_u*fs_w
!loc_pyz = loc_rho*fs_v*fs_w
!call calculate_lagrangian_mults(loc_rho,fs_u,fs_v,fs_w,loc_temp,&
!  loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,&
!  fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
!  fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z,loc_fieq,a,b,c)
!fs_fieq = loc_fieq
!call calculate_alpha(fs_fieq,fs_fieq,fs_fieq,fs_alpha)
!fs_gieq = fs_fieq*(2*const_vol-dimensions)*temperature

!if (self == 0) then
!  write(*,*) fs_fieq
!  write(*,*) fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
!  fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z
!end if

write(11,*) 'Freestream Langrange multiplier values',fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
  fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z
write(11,*) 'Freestream alpha value',fs_alpha
write(11,*) 'freestream fieq values',fs_fieq
write(11,*) 'freestream gieq values',fs_gieq
!write(*,*) 'Freestream values ',self,rho_inf,fs_u,fs_v,fs_w,temperature
!write(*,*) 'Freestream Langrange multiplier values',self,fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,fs_piyy,&
!  fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z
!write(*,*) 'Freestream alpha value',self,fs_alpha
!write(*,*) 'freestream fieq values',self,fs_fieq
!write(*,*) 'freestream gieq values',self,fs_gieq

CALL quick_sutherlands(fs_mu,temperature,.false.)
CALL quick_therm_cond(fs_kappa,rho_inf,temperature)
fs_beta_1 = fs_alpha/((2.0D0*fs_mu)/(rho_inf*temperature)+1.0D0)
fs_beta_2 = 2.0D0/((2.0D0*fs_kappa)/(rho_inf*const_press*temperature)+1.0D0)
!write(*,*) 'coefficient values',fs_alpha,fs_beta_1,fs_beta_2
!
!
!
end subroutine
