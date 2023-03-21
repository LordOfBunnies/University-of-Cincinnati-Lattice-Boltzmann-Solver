SUBROUTINE interlevel_interpolation(lvl,fine)
!
! Performs interpolation between levels on timesteps where there are
! more than one level to be run.
!
! Called by: loop
! Calls:
! External calls: amrex_fillpatch,amrex_fill_boundary
!
use precise
use constants
use linkwise
use grid_data
use freestream_values
use quick_calcs
use derivatives
use amrex_base_module
use amrex_amr_module
use amr_info_holder
use mpi
use fill_stuff, only: edge_bounds,lo_bc,hi_bc,lo_bc_fi,hi_bc_fi,edge_bounds_distro
use amr_processes
implicit none
integer :: fine,i,j,lvl,a,b,c
!integer :: box_lo(4),box_hi(4)
integer :: upwind,downwind,center,adj_count,n_paired_nodes
integer :: second_upwind,second_downwind,second_central
logical :: x_up_clear,x_down_clear,y_up_clear,y_down_clear,z_up_clear,z_down_clear
logical :: x_up_close,x_down_close,y_up_close,y_down_close,z_up_close,z_down_close
logical :: x_null,y_null,z_null,valid_interp_nodes(27)
logical,allocatable :: needs_interp(:,:,:),node_pairs(:)
real(kind=dp) :: dTdx,dTdy,dTdz,loc_fieq,k_fineq,k_gineq,loc_fineq,loc_gineq
real(kind=dp) :: loc_gieq,temp_der_sum,crs_fieq,crs_gieq
real(kind=dp) :: beta_1_c,beta_1_f,beta_2_f,beta_2_c,mu_c,mu_f,kappa_c,kappa_f
real(kind=dp) :: beta_1_rat,beta_2_rat,loc_time,loc_epsilon
real(kind=dp) :: interp_fix(dir),interp_fiy(dir),interp_fiz(dir)
real(kind=dp) :: interp_gix(dir),interp_giy(dir),interp_giz(dir)
real(kind=dp) :: qabc(27)
!real(kind=dp) :: interp_chi_x,interp_chi_y,interp_chi_z,interp_zetax_x,&
!  interp_zetax_y,interp_zetax_z,interp_zetay_x,interp_zetay_y,interp_zetay_z,&
!  interp_zetaz_x,interp_zetaz_y,interp_zetaz_z
!real(kind=dp) :: interp_pixx_x,interp_pixx_y,interp_pixx_z,interp_piyy_x,&
!  interp_piyy_y,interp_piyy_z,interp_pizz_x,interp_pizz_y,interp_pizz_z,&
!  interp_pixy_x,interp_pixy_y,interp_pixy_z,interp_pixz_x,interp_pixz_y,&
!  interp_pixz_z,interp_piyz_x,interp_piyz_y,interp_piyz_z
!real(kind=dp) :: interp_lambdax_x,interp_lambdax_y,interp_lambdax_z,interp_lambday_x,&
!  interp_lambday_y,interp_lambday_z,interp_lambdaz_x,interp_lambdaz_y,interp_lambdaz_z
real(kind=dp) :: rho_sum,u_sum,v_sum,w_sum,energy_sum, pxx,pyy,pzz,pxy,pxz,pyz,fieq_skip(dir),&
  v_mag,loc_fi_diff(dir)
real(kind=dp) :: dir_sum_fi,dir_sum_gi,dummy_time
real(kind=dp),allocatable :: old_fi(:)

!real(kind=dp) :: q_xxx,q_yyy,q_zzz,q_xxy,q_xxz,q_xyy,q_xzz,q_xyz,q_yyz,q_yzz
!real(kind=dp) :: p_xx,p_yy,p_zz,p_xy,p_xz,p_yz

real(kind=dp),contiguous,pointer :: crs_temp(:,:,:,:),crs_rho(:,:,:,:),&
  crs_chi(:,:,:,:),crs_zeta_x(:,:,:,:),&
  crs_zeta_y(:,:,:,:),crs_zeta_z(:,:,:,:),crs_pi_xx(:,:,:,:),crs_pi_yy(:,:,:,:),crs_pi_zz(:,:,:,:),&
  crs_pi_yz(:,:,:,:),crs_pi_xz(:,:,:,:),crs_pi_xy(:,:,:,:),crs_lambda_x(:,:,:,:),crs_lambda_y(:,:,:,:),&
  crs_lambda_z(:,:,:,:),crs_gi(:,:,:,:),crs_fi(:,:,:,:)
real(kind=dp),contiguous,pointer :: crs_u(:,:,:,:),crs_v(:,:,:,:),crs_w(:,:,:,:),crs_alpha(:,:,:,:)
!
type(amrex_multifab) :: tmp_temp,tmp_chi,tmp_zeta_x,tmp_zeta_y,tmp_zeta_z,tmp_pi_xx,&
  tmp_pi_yy,tmp_pi_zz,tmp_pi_xy,tmp_pi_xz,tmp_pi_yz,tmp_lambda_z,tmp_lambda_y,tmp_lambda_x,&
  tmp_rho
type(amrex_multifab) :: tmp_u,tmp_v,tmp_w,tmp_alpha
type(amrex_multifab) :: tmp_fi,tmp_gi
type(amrex_mfiter) :: mfi
type(amrex_box) :: ducks

if (self == 0) then
  write(*,*) 'interlevel interpolation',lvl,self
end if
!p_xx = 0.0D0
!p_yy = 0.0D0
!p_zz = 0.0D0
!p_xy = 0.0D0
!p_xz = 0.0D0
!p_yz = 0.0D0
!
!q_xxx = 0.0D0
!q_yyy = 0.0D0
!q_zzz = 0.0D0
!q_xxy = 0.0D0
!q_xxz = 0.0D0
!q_xyy = 0.0D0
!q_xzz = 0.0D0
!q_xyz = 0.0D0
!q_yyz = 0.0D0
!q_yzz = 0.0D0
!
!
!do i = 1,dir
!  p_xx = p_xx - rcx(i)**2*loc_fidiff(i)
!  p_yy = p_yy - rcy(i)**2*loc_fidiff(i)
!  p_zz = p_zz - rcz(i)**2*loc_fidiff(i)
!  p_xy = p_xy - rcx(i)*rcy(i)*loc_fidiff(i)
!  p_xz = p_xz - rcx(i)*rcz(i)*loc_fidiff(i)
!  p_yz = p_yz - rcy(i)*rcz(i)*loc_fidiff(i)
!
!  q_xxx = q_xxx - rcx(i)**3 * loc_fidiff(i)
!  q_yyy = q_yyy - rcy(i)**3 * loc_fidiff(i)
!  q_zzz = q_zzz - rcz(i)**3 * loc_fidiff(i)
!  q_xxy = q_xxy - rcx(i)**2 * rcy(i) * loc_fidiff(i)
!  q_xxz = q_xxz - rcx(i)**2 * rcz(i) * loc_fidiff(i)
!  q_xyy = q_xyy - rcx(i) * rcy(i)**2 * loc_fidiff(i)
!  q_xzz = q_xzz - rcx(i) * rcz(i)**2 * loc_fidiff(i)
!  q_xyz = q_xyz - rcx(i) * rcy(i) * rcz(i) * loc_fidiff(i)
!  q_yyz = q_yyz - rcy(i)**2 * rcz(i) * loc_fidiff(i)
!  q_yzz = q_yzz - rcz(i)**2 * rcy(i) * loc_fidiff(i)
!
!end do

if (mod(timestep(lvl),2) == 1) then
  loc_time = time(lvl)+dt(lvl)
else
  return
!  switcheroo = 1
end if

allocate(old_fi(dir))

dummy_time = time(lvl-1) + 0.00001D0

if (lvl /= 0) then
!write(*,*) 'time check',loc_time,time_prev(lvl),time,time_prev
!
! create the mfiter so we can go through all the blocks
!

!
! Find the local coarse conditions interpolated onto the finer mesh
!
  call amrex_multifab_build(tmp_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_rho,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_u,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_v,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_w,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_alpha,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
!
  call amrex_multifab_build(tmp_chi,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_zeta_x,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_zeta_y,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_zeta_z,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_pi_xx,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_pi_yy,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_pi_zz,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_pi_xy,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_pi_xz,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_pi_yz,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_lambda_x,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_lambda_y,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_lambda_z,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
!
  call amrex_multifab_build(tmp_fi,mffi(lvl)%ba,mffi(lvl)%dm,dir,nghosts_mid,node_based_3d)
  call amrex_multifab_build(tmp_gi,mfgi(lvl)%ba,mfgi(lvl)%dm,dir,nghosts_mid,node_based_3d)
!
!
!
! write(*,*) 'mrow',lvl,lo_bc,hi_bc,amrex_ref_ratio(lvl-1),loc_time
  call amrex_fillcoarsepatch(tmp_temp,time_prev(lvl-1),mftemp(lvl-1),dummy_time,mftemp(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,dummy_time,1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  write(*,*) 'woof',lo_bc,hi_bc
!
!
!
  call amrex_fillcoarsepatch(tmp_rho,time_prev(lvl-1),mfrho(lvl-1),time(lvl-1),mfrho(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_u,time_prev(lvl-1),mfu_vel(lvl-1),time(lvl-1),mfu_vel(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_v,time_prev(lvl-1),mfv_vel(lvl-1),time(lvl-1),mfv_vel(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_w,time_prev(lvl-1),mfw_vel(lvl-1),time(lvl-1),mfw_vel(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_alpha,time_prev(lvl-1),mfalpha(lvl-1),time(lvl-1),mfalpha(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_chi,time_prev(lvl-1),mfchi(lvl-1),time(lvl-1),mfchi(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_zeta_x,time_prev(lvl-1),mfzeta_x(lvl-1),time(lvl-1),mfzeta_x(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_zeta_y,time_prev(lvl-1),mfzeta_y(lvl-1),time(lvl-1),mfzeta_y(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_zeta_z,time_prev(lvl-1),mfzeta_z(lvl-1),time(lvl-1),mfzeta_z(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_pi_xx,time_prev(lvl-1),mfpi_xx(lvl-1),time(lvl-1),mfpi_xx(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_pi_yy,time_prev(lvl-1),mfpi_yy(lvl-1),time(lvl-1),mfpi_yy(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_pi_zz,time_prev(lvl-1),mfpi_zz(lvl-1),time(lvl-1),mfpi_zz(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_pi_xy,time_prev(lvl-1),mfpi_xy(lvl-1),time(lvl-1),mfpi_xy(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_pi_xz,time_prev(lvl-1),mfpi_xz(lvl-1),time(lvl-1),mfpi_xz(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_pi_yz,time_prev(lvl-1),mfpi_yz(lvl-1),time(lvl-1),mfpi_yz(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_lambda_x,time_prev(lvl-1),mflambda_x(lvl-1),time(lvl-1),mflambda_x(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_lambda_y,time_prev(lvl-1),mflambda_y(lvl-1),time(lvl-1),mflambda_y(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(tmp_lambda_z,time_prev(lvl-1),mflambda_z(lvl-1),time(lvl-1),mflambda_z(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!!
!!
!  write(*,*) 'meow',self
!!
!  call amrex_fillcoarsepatch(tmp_fi,time_prev(lvl-1),mffi(lvl-1),loc_time,mffi(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,loc_time,1,1,dir,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
!  call amrex_fillcoarsepatch(tmp_gi,time_prev(lvl-1),mfgi(lvl-1),loc_time,mfgi(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,loc_time,1,1,dir,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)

   call amrex_fillcoarsepatch(tmp_fi,time_prev(lvl-1),mffi(lvl-1),time(lvl-1),mffi(lvl-1),&
    amrex_geom(lvl-1),edge_bounds_distro,amrex_geom(lvl),edge_bounds_distro,time(lvl-1),1,1,dir,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)

  call amrex_fillcoarsepatch(tmp_gi,time_prev(lvl-1),mfgi(lvl-1),time(lvl-1),mfgi(lvl-1),&
    amrex_geom(lvl-1),edge_bounds_distro,amrex_geom(lvl),edge_bounds_distro,time(lvl-1),1,1,dir,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
! Fill the ghost nodes that overlap other boxes on the same level
!
  call mffi(lvl)% fill_boundary(amrex_geom(lvl),1,dir)
  call mfgi(lvl)% fill_boundary(amrex_geom(lvl),1,dir)

  call mfrho(lvl)% fill_boundary(amrex_geom(lvl))
  call mfu_vel(lvl)% fill_boundary(amrex_geom(lvl))
  call mfv_vel(lvl)% fill_boundary(amrex_geom(lvl))
  call mftemp(lvl)% fill_boundary(amrex_geom(lvl))
  call mfalpha(lvl)% fill_boundary(amrex_geom(lvl))

  call mfchi(lvl)% fill_boundary(amrex_geom(lvl))
  call mfzeta_x(lvl)% fill_boundary(amrex_geom(lvl))
  call mfzeta_y(lvl)% fill_boundary(amrex_geom(lvl))
  call mfpi_xx(lvl)% fill_boundary(amrex_geom(lvl))
  call mfpi_xy(lvl)% fill_boundary(amrex_geom(lvl))
  call mfpi_yy(lvl)% fill_boundary(amrex_geom(lvl))
  call mflambda_x(lvl)% fill_boundary(amrex_geom(lvl))
  call mflambda_y(lvl)% fill_boundary(amrex_geom(lvl))

  if (dimensions == 3) then
    call mfw_vel(lvl)% fill_boundary(amrex_geom(lvl))
    call mfzeta_z(lvl)% fill_boundary(amrex_geom(lvl))
    call mfpi_xz(lvl)% fill_boundary(amrex_geom(lvl))
    call mfpi_yz(lvl)% fill_boundary(amrex_geom(lvl))
    call mfpi_zz(lvl)% fill_boundary(amrex_geom(lvl))
    call mflambda_z(lvl)% fill_boundary(amrex_geom(lvl))
  end if
!  write(*,*) 'bark',self
!
!
! Use fill_patch on the finer to get the flow information down to it
!
!
! Go through all the boxes on the level to interpolate between the levels
!
!
  call amrex_mfiter_build(mfi,mffi(lvl),tiling=.false.)

  do while (mfi%next())
!    write(*,*) 'inside the first box sans problems'
    ducks = mfi%validbox()

    crs_temp => tmp_temp%dataptr(mfi)
    crs_rho => tmp_rho%dataptr(mfi)
    crs_chi => tmp_chi%dataptr(mfi)
    crs_zeta_X => tmp_zeta_x%dataptr(mfi)
    crs_zeta_y => tmp_zeta_y%dataptr(mfi)
    crs_zeta_z => tmp_zeta_z%dataptr(mfi)
    crs_pi_xx => tmp_pi_xx%dataptr(mfi)
    crs_pi_xy => tmp_pi_xy%dataptr(mfi)
    crs_pi_xz => tmp_pi_xz%dataptr(mfi)
    crs_pi_yy => tmp_pi_yy%dataptr(mfi)
    crs_pi_zz => tmp_pi_zz%dataptr(mfi)
    crs_pi_yz => tmp_pi_yz%dataptr(mfi)
    crs_lambda_x => tmp_lambda_x%dataptr(mfi)
    crs_lambda_y => tmp_lambda_y%dataptr(mfi)
    crs_lambda_z => tmp_lambda_z%dataptr(mfi)

    crs_fi => tmp_fi%dataptr(mfi)
    crs_gi => tmp_gi%dataptr(mfi)

    crs_u => tmp_u%dataptr(mfi)
    crs_v => tmp_v%dataptr(mfi)
    crs_w => tmp_w%dataptr(mfi)
    crs_alpha => tmp_alpha%dataptr(mfi)

    rho => mfrho(lvl)%dataptr(mfi)
    temp => mftemp(lvl)%dataptr(mfi)
    fi => mffi(lvl)%dataptr(mfi)
    gi => mfgi(lvl)%dataptr(mfi)
    state => mfstate(lvl)%dataptr(mfi)
    alpha => mfalpha(lvl)%dataptr(mfi)

    u_vel => mfu_vel(lvl)%dataptr(mfi)
    v_vel => mfv_vel(lvl)%dataptr(mfi)
    w_vel => mfw_vel(lvl)%dataptr(mfi)

    chi => mfchi(lvl)%dataptr(mfi)
    zeta_X => mfzeta_x(lvl)%dataptr(mfi)
    zeta_y => mfzeta_y(lvl)%dataptr(mfi)
    zeta_z => mfzeta_z(lvl)%dataptr(mfi)
    pi_xx => mfpi_xx(lvl)%dataptr(mfi)
    pi_xy => mfpi_xy(lvl)%dataptr(mfi)
    pi_xz => mfpi_xz(lvl)%dataptr(mfi)
    pi_yy => mfpi_yy(lvl)%dataptr(mfi)
    pi_zz => mfpi_zz(lvl)%dataptr(mfi)
    pi_yz => mfpi_yz(lvl)%dataptr(mfi)
    lambda_x => mflambda_x(lvl)%dataptr(mfi)
    lambda_y => mflambda_y(lvl)%dataptr(mfi)
    lambda_z => mflambda_z(lvl)%dataptr(mfi)

    if (allocated(needs_interp)) deallocate(needs_interp)
    allocate(needs_interp(ducks%lo(1)-nghosts_mid:ducks%hi(1)+nghosts_mid,&
      ducks%lo(2)-nghosts_mid:ducks%hi(2)+nghosts_mid,&
      ducks%lo(3)-nghosts_mid:ducks%hi(3)+nghosts_mid))
    needs_interp = .false.
!
! Interpolates to all ghost nodes, first the reconstruction nodes are done, those fine nodes that overlap
!   with coarse nodes.  Then nodes on edges, faces, and then at the center of cubes.  These values need to be
!
!
    do c = ducks%lo(3)-nghosts_mid,ducks%hi(3)+nghosts_mid
      do b = ducks%lo(2)-nghosts_mid,ducks%hi(2)+nghosts_mid
        do a = ducks%lo(1)-nghosts_mid,ducks%hi(1)+nghosts_mid
          if (state(a,b,c,1) >=10000000 .and. state(a,b,c,1) <= 1100000000) then
            old_fi(1:dir) = fi(a,b,c,1:dir)
            valid_interp_nodes =.false.

! Copy the interpolated coarse info into the ghost nodes, fi and gi are replaced as part of the
!   interpolation process
!            rho(a,b,c,1) = crs_rho(a,b,c,1)
!            u_vel(a,b,c,1) = crs_u(a,b,c,1)
!            v_vel(a,b,c,1) = crs_v(a,b,c,1)
!            w_vel(a,b,c,1) = crs_w(a,b,c,1)
!            temp(a,b,c,1) = crs_temp(a,b,c,1)
!            alpha(a,b,c,1) = crs_alpha(a,b,c,1)
!
            chi(a,b,c,1) = crs_chi(a,b,c,1)
            zeta_x(a,b,c,1) = crs_zeta_x(a,b,c,1)
            zeta_y(a,b,c,1) = crs_zeta_y(a,b,c,1)
            zeta_z(a,b,c,1) = crs_zeta_z(a,b,c,1)
            pi_xx(a,b,c,1) = crs_pi_xx(a,b,c,1)
            pi_yy(a,b,c,1) = crs_pi_yy(a,b,c,1)
            pi_zz(a,b,c,1) = crs_pi_zz(a,b,c,1)
            pi_xy(a,b,c,1) = crs_pi_xy(a,b,c,1)
            pi_xz(a,b,c,1) = crs_pi_xz(a,b,c,1)
            pi_yz(a,b,c,1) = crs_pi_yz(a,b,c,1)
            lambda_x(a,b,c,1) = crs_lambda_x(a,b,c,1)
            lambda_y(a,b,c,1) = crs_lambda_y(a,b,c,1)
            lambda_z(a,b,c,1) = crs_lambda_z(a,b,c,1)


!            if (mod(a,2) == 0 .and. mod(b,2) == 0 .and. mod(c,2) == 0) then
!
!            second_upwind = mod(state(a,b,c,1),1000000)/100000
!            second_downwind = mod(state(a,b,c,1),100000)/10000
!            second_central = mod(state(a,b,c,1),10000)/1000
!            upwind = mod(state(a,b,c,1),1000)/100
!            downwind = mod(state(a,b,c,1),100)/10
!            center = mod(state(a,b,c,1),10)

            call quick_sutherlands(mu_c,crs_temp(a,b,c,1),.false.)
            call quick_sutherlands(mu_f,temp(a,b,c,1),.true.)
            call quick_therm_cond(kappa_c,crs_rho(a,b,c,1),crs_temp(a,b,c,1))
            call quick_therm_cond(kappa_f,rho(a,b,c,1),temp(a,b,c,1))

              rho_sum = 0.0D0
              u_sum = 0.0D0
              v_sum = 0.0D0
              w_sum = 0.0D0
              energy_sum = 0.0D0
              do i = 1,dir
                fi(a,b,c,i) = crs_fi(a,b,c,i)
                gi(a,b,c,i) = crs_gi(a,b,c,i)



                if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0 .or. isnan(fi(a,b,c,i)) .or. &
                    isnan(gi(a,b,c,i))) then
                  if (.not. shifted) then
                  fi(a,b,c,i) = fieq_comp_39(crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
                    crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
                    crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
                    crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1),i)
                  gi(a,b,c,i) = fi(a,b,c,i)*((2*const_vol-dimensions)*crs_temp(a,b,c,1))
                  else
                  fi(a,b,c,i) = fieq_comp_39_shifted(crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
                    crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
                    crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
                    crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1),i)
                  gi(a,b,c,i) = fi(a,b,c,i)*((2*const_vol-dimensions)*crs_temp(a,b,c,1))
                  end if
                end if

                rho_sum = rho_sum + fi(a,b,c,i)
                u_sum = u_sum + rcx(i)*fi(a,b,c,i)
                v_sum = v_sum + rcy(i)*fi(a,b,c,i)
                w_sum = w_sum + rcz(i)*fi(a,b,c,i)
                energy_sum = energy_sum + cmag(i)*fi(a,b,c,i)+gi(a,b,c,i)


              end do

!            if (a == 1388 .and. b == 1013 .and. c == 376) then
!              write(*,*) 'wrinkle feet',fi(a,b,c,1:dir)
!              write(*,*) 'droopy jowls',gi(a,b,c,1:dir)
!            end if

              rho(a,b,c,1) = rho_sum
              u_vel(a,b,c,1) = u_sum/rho_sum
              v_vel(a,b,c,1) = v_sum/rho_sum
              w_vel(a,b,c,1) = w_sum/rho_sum
              v_mag = u_vel(a,b,c,1)**2 + v_vel(a,b,c,1)**2 + w_vel(a,b,c,1)**2
              temp(a,b,c,1) = (energy_sum - rho(a,b,c,1)*v_mag)/(2*const_vol*rho(a,b,c,1))

              pxx = rho(a,b,c,1)*(temp(a,b,c,1) + u_vel(a,b,c,1)**2)
              pyy = rho(a,b,c,1)*(temp(a,b,c,1) + v_vel(a,b,c,1)**2)
              pzz = rho(a,b,c,1)*(temp(a,b,c,1) + w_vel(a,b,c,1)**2)
              pxy = rho(a,b,c,1)*u_vel(a,b,c,1)*v_vel(a,b,c,1)
              pxz = rho(a,b,c,1)*u_vel(a,b,c,1)*w_vel(a,b,c,1)
              pyz = rho(a,b,c,1)*v_vel(a,b,c,1)*w_vel(a,b,c,1)
!
              CALL calculate_lagrangian_mults(rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
                w_vel(a,b,c,1),temp(a,b,c,1),pxx,pyy,pzz,pxy,pxz,pyz,&
                chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),fieq_skip,a,b,c)
!
!              do i = 1,dir
!                fieq_skip(i) = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                  pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!                  pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
!              end do

!              do i = 1,dir
!                if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!                  fi(a,b,c,i) = fieq_skip(i)
!                  gi(a,b,c,i) = (2.0D0*const_vol - dimensions)*temp(a,b,c,1)*fieq_skip(i)
!                end if
!              end do
!            if (a == 290 .and. b == 105 .and. c == 127) then
!              write(*,*) 'inter fis 2',fi(a,b,c,1:dir)
!              write(*,*) 'inter gis 2',gi(a,b,c,1:dir)
!              write(*,*) 'inter coarse fieqs 2',fieq_skip(1:dir)
!              !write(*,*) 'inter coarse gieqs 2',crs_gi(a,b,c,1:dir)
!            end if
              loc_fi_diff = fieq_skip(1:dir)-fi(a,b,c,1:dir)

              if (turbulent) then
                CALL calculate_alpha(fi(a,b,c,1:dir),loc_fi_diff(1:dir),alpha(a,b,c,1),temp(a,b,c,1),&
                  state(a,b,c,1),lvl,a,b,c)
              else
                alpha(a,b,c,1) = 2.0D0
              end if

          end if
        end do
      end do
    end do

  end do

!
!
!
!
!
!
!


!    do c = ducks%lo(3)-nghosts_mid,ducks%hi(3)+nghosts_mid
!      do b = ducks%lo(2)-nghosts_mid,ducks%hi(2)+nghosts_mid
!        do a = ducks%lo(1)-nghosts_mid,ducks%hi(1)+nghosts_mid
!          if (state(a,b,c,1) >=10000000 .and. state(a,b,c,1) <= 1100000000) then
!
!            if (a == 290 .and. b == 105 .and. c == 127) then
!              write(*,*) 'inter fis 1',fi(a,b,c,1:dir)
!              write(*,*) 'inter gis 1',gi(a,b,c,1:dir)
!              write(*,*) 'inter coarse fis 1',crs_fi(a,b,c,1:dir)
!              write(*,*) 'inter coarse gis 1',crs_gi(a,b,c,1:dir)
!            end if
!!
!! Last ditch effort to get valid data into interpolation region
!!
!              if (any(isnan(fi(a,b,c,1:dir))) .or. any(isnan(gi(a,b,c,1:dir)))) then
!!                write(*,*) 'spackle',a,b,c,ducks%lo,ducks%hi,lvl
!!                write(*,*) 'shazbot',fi(a,b,c,1:dir)
!!                write(*,*) 'shingle',gi(a,b,c,1:dir)
!                if (.not. shifted) then
!                do i = 1,dir
!                  if (isnan(fi(a,b,c,i)) .and. .not. isnan(gi(a,b,c,i))) then
!                    fi(a,b,c,i) = gi(a,b,c,i)/((2*const_vol-dimensions)*crs_temp(a,b,c,1))
!                  else if (isnan(gi(a,b,c,i)) .and. .not. isnan(fi(a,b,c,i))) then
!                    gi(a,b,c,i) = fi(a,b,c,i)*((2*const_vol-dimensions)*crs_temp(a,b,c,1))
!                  else
!! Safety check
!                    if (state(a,b,c,i) >= 0) then
!                      if (.not. isnan(fi(a+cx(i),b+cy(i),c+cz(i),i))) then
!                        fi(a,b,c,i) = fi(a+cx(i),b+cy(i),c+cz(i),i)
!                      else if (.not. isnan(fi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i))) then
!                        fi(a,b,c,i) = fi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i)
!                      else
!                        fi(a,b,c,i) = 0.08D0
!                      end if
!! Opposite direction safety check
!                    else if (state(a,b,c,opp(i)) >= 0) then
!                      if (.not. isnan(gi(a+cx(i),b+cy(i),c+cz(i),i))) then
!                        gi(a,b,c,i) = gi(a+cx(i),b+cy(i),c+cz(i),i)
!                      else if (.not. isnan(gi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i))) then
!                        gi(a,b,c,i) = gi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i)
!                      else
!                        gi(a,b,c,i) = 0.112D0
!                      end if
!! Failsafe
!                    else
!!                      fi(a,b,c,i) = 0.08D0
!!                      gi(a,b,c,i) = 0.112D0
!                      fi(a,b,c,i) = fs_fieq(i)
!                      gi(a,b,c,i) = fs_gieq(i)
!                    end if
!                  end if
!                end do
!!
!! Shifted lattice sanity checks
!!
!                else
!                  do i = 1,dir
!                  if (isnan(fi(a,b,c,i)) .and. .not. isnan(gi(a,b,c,i))) then
!                    fi(a,b,c,i) = gi(a,b,c,i)/((2*const_vol-dimensions)*crs_temp(a,b,c,1))
!                  else if (isnan(gi(a,b,c,i)) .and. .not. isnan(fi(a,b,c,i))) then
!                    gi(a,b,c,i) = fi(a,b,c,i)*((2*const_vol-dimensions)*crs_temp(a,b,c,1))
!                  else
!! Safety check
!                    if (state(a,b,c,i) >= 0) then
!                      if (.not. isnan(fi(a+cx(i),b+cy(i),c+cz(i),i))) then
!                        fi(a,b,c,i) = fi(a+cx(i),b+cy(i),c+cz(i),i)
!                      else if (.not. isnan(fi(a+cx(erp(i)),b+cy(erp(i)),c+cz(erp(i)),i))) then
!                        fi(a,b,c,i) = fi(a+cx(erp(i)),b+cy(erp(i)),c+cz(erp(i)),i)
!                      else
!                        fi(a,b,c,i) = 0.08D0
!                      end if
!! Opposite direction safety check
!                    else if (state(a,b,c,erp(i)) >= 0) then
!                      if (.not. isnan(gi(a+cx(i),b+cy(i),c+cz(i),i))) then
!                        gi(a,b,c,i) = gi(a+cx(i),b+cy(i),c+cz(i),i)
!                      else if (.not. isnan(gi(a+cx(erp(i)),b+cy(erp(i)),c+cz(erp(i)),i))) then
!                        gi(a,b,c,i) = gi(a+cx(erp(i)),b+cy(erp(i)),c+cz(erp(i)),i)
!                      else
!                        gi(a,b,c,i) = 0.112D0
!                      end if
!! Failsafe
!                    else
!!                      fi(a,b,c,i) = 0.08D0
!!                      gi(a,b,c,i) = 0.112D0
!                      fi(a,b,c,i) = fs_fieq(i)
!                      gi(a,b,c,i) = fs_gieq(i)
!                    end if
!                  end if
!                end do
!                end if
!              end if
!
!
!!
!!            end if
!          end if
!        end do
!      end do
!    end do
!
!  end do
!
! Destroy the temporary things
!
  call amrex_mfiter_destroy(mfi)

  call amrex_multifab_destroy(tmp_rho)
  call amrex_multifab_destroy(tmp_u)
  call amrex_multifab_destroy(tmp_v)
  call amrex_multifab_destroy(tmp_w)
  call amrex_multifab_destroy(tmp_temp)
  call amrex_multifab_destroy(tmp_alpha)

  call amrex_multifab_destroy(tmp_chi)
  call amrex_multifab_destroy(tmp_zeta_x)
  call amrex_multifab_destroy(tmp_zeta_y)
  call amrex_multifab_destroy(tmp_zeta_z)
  call amrex_multifab_destroy(tmp_pi_xx)
  call amrex_multifab_destroy(tmp_pi_yy)
  call amrex_multifab_destroy(tmp_pi_zz)
  call amrex_multifab_destroy(tmp_pi_xy)
  call amrex_multifab_destroy(tmp_pi_xz)
  call amrex_multifab_destroy(tmp_pi_yz)
  call amrex_multifab_destroy(tmp_lambda_x)
  call amrex_multifab_destroy(tmp_lambda_y)
  call amrex_multifab_destroy(tmp_lambda_z)

  call amrex_multifab_destroy(tmp_fi)
  call amrex_multifab_destroy(tmp_gi)



end if
deallocate(old_fi)

END SUBROUTINE

!                  write(*,*) 'how the fluff?',fi(a,b,c,1:dir),a,b,c
!                  write(*,*) 'though I fixed this',gi(a,b,c,1:dir),a,b,c
!                  write(*,*) 'coarse fis',crs_fi(a,b,c,1:dir)
!                  write(*,*) 'coarse gis',crs_gi(a,b,c,1:dir)


!            write(*,*) 'wuf wuf'
!           if (a == 51 .and. b == 37 .and. c == 2 .and. lvl == 1) then
!
!             write(*,*) 'ARGH!!!',state(a,b,c,1),state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!                 state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),a,b,c,ducks%lo,ducks%hi,self
!             write(*,*) 'temperature stuff',temp(a,b,c,1),temp(a+1,b,c,1),temp(a,b+1,c,1),&
!               temp(a,b,c+1,1),temp(a-1,b,c,1),temp(a,b-1,c,1),temp(a,b,c-1,1)
!             write(*,*) 'coarse stuff',crs_temp(a,b,c,1),crs_temp(a+1,b,c,1),crs_temp(a,b+1,c,1),&
!               crs_temp(a,b,c+1,1),crs_temp(a-1,b,c,1),crs_temp(a,b-1,c,1),crs_temp(a,b,c-1,1)
!           end if
!
!
!            write(*,*) 'interpolated fi and gi',fi(a,b,c,i),loc_fieq,gi(a,b,c,i),loc_gieq,loc_gineq,beta_1_rat,beta_2_rat,i,lvl
!            write(*,*) 'coarse values',crs_temp(a,b,c,1),crs_fi(a,b,c,i),crs_fieq,crs_gi(a,b,c,i),crs_gieq
!            write(*,*) 'normal LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!                    zeta_y(a,b,c,1),zeta_z(a,b,c,1),pi_xx(a,b,c,1),pi_yy(a,b,c,1),&
!                    pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz(a,b,c,1),&
!                    lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!            write(*,*) 'coarse LMs',crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
!                    crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
!                    crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
!                    crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1)


!
!if (MODULO(timestep(lvl),2) == 0) then
!!  call amrex_average_down()
!
!
!
!else
!
!
!
!end if
!
!
!
!  call amrex_mfiter_build(mfi,mffi(lvl),tiling=.false.)
!
!  do while (mfi%next())
!
!    ducks = mfi%validbox()
!
!    crs_temp => tmp_temp%dataptr(mfi)
!    crs_rho => tmp_rho%dataptr(mfi)
!    crs_chi => tmp_chi%dataptr(mfi)
!    crs_zeta_X => tmp_zeta_x%dataptr(mfi)
!    crs_zeta_y => tmp_zeta_y%dataptr(mfi)
!    crs_zeta_z => tmp_zeta_z%dataptr(mfi)
!    crs_pi_xx => tmp_pi_xx%dataptr(mfi)
!    crs_pi_xy => tmp_pi_xy%dataptr(mfi)
!    crs_pi_xz => tmp_pi_xz%dataptr(mfi)
!    crs_pi_yy => tmp_pi_yy%dataptr(mfi)
!    crs_pi_zz => tmp_pi_zz%dataptr(mfi)
!    crs_pi_yz => tmp_pi_yz%dataptr(mfi)
!    crs_lambda_x => tmp_lambda_x%dataptr(mfi)
!    crs_lambda_y => tmp_lambda_y%dataptr(mfi)
!    crs_lambda_z => tmp_lambda_z%dataptr(mfi)
!
!    crs_fi => tmp_fi%dataptr(mfi)
!    crs_gi => tmp_gi%dataptr(mfi)
!
!    crs_u => tmp_u%dataptr(mfi)
!    crs_v => tmp_v%dataptr(mfi)
!    crs_w => tmp_w%dataptr(mfi)
!
!    rho => mfrho(lvl)%dataptr(mfi)
!    temp => mftemp(lvl)%dataptr(mfi)
!    fi => mffi(lvl)%dataptr(mfi)
!    gi => mfgi(lvl)%dataptr(mfi)
!    state => mfstate(lvl)%dataptr(mfi)
!
!    u_vel => mfu_vel(lvl)%dataptr(mfi)
!    v_vel => mfv_vel(lvl)%dataptr(mfi)
!    w_vel => mfw_vel(lvl)%dataptr(mfi)
!
!    chi => mfchi(lvl)%dataptr(mfi)
!    zeta_X => mfzeta_x(lvl)%dataptr(mfi)
!    zeta_y => mfzeta_y(lvl)%dataptr(mfi)
!    zeta_z => mfzeta_z(lvl)%dataptr(mfi)
!    pi_xx => mfpi_xx(lvl)%dataptr(mfi)
!    pi_xy => mfpi_xy(lvl)%dataptr(mfi)
!    pi_xz => mfpi_xz(lvl)%dataptr(mfi)
!    pi_yy => mfpi_yy(lvl)%dataptr(mfi)
!    pi_zz => mfpi_zz(lvl)%dataptr(mfi)
!    pi_yz => mfpi_yz(lvl)%dataptr(mfi)
!    lambda_x => mflambda_x(lvl)%dataptr(mfi)
!    lambda_y => mflambda_y(lvl)%dataptr(mfi)
!    lambda_z => mflambda_z(lvl)%dataptr(mfi)
!!
!! Interpolates to all ghost nodes
!!
!    do c = ducks%lo(3)-nghosts,ducks%hi(3)+nghosts
!      do b = ducks%lo(2)-nghosts,ducks%hi(2)+nghosts
!        do a = ducks%lo(1)-nghosts,ducks%hi(1)+nghosts
!          if (state(a,b,c,1) >=1000 .and. state(a,b,c,1) <= 50999) then
!
!
!
!          end if
!
!        end do
!      end do
!    end do
!  end do

!
!
!
!  call amrex_fillpatch(mfrho(lvl),time_prev(lvl-1),mfrho(lvl-1),loc_time,mfrho(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfrho(lvl),loc_time,mfrho(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfu_vel(lvl),time_prev(lvl-1),mfu_vel(lvl-1),loc_time,mfu_vel(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfu_vel(lvl),loc_time,mfu_vel(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfv_vel(lvl),time_prev(lvl-1),mfv_vel(lvl-1),loc_time,mfv_vel(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfv_vel(lvl),loc_time,mfv_vel(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfw_vel(lvl),time_prev(lvl-1),mfw_vel(lvl-1),loc_time,mfw_vel(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfw_vel(lvl),loc_time,mfw_vel(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mftemp(lvl),time_prev(lvl-1),mftemp(lvl-1),loc_time,mftemp(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mftemp(lvl),loc_time,mftemp(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfchi(lvl),time_prev(lvl-1),mfchi(lvl-1),loc_time,mfchi(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfchi(lvl),loc_time,mfchi(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfzeta_x(lvl),time_prev(lvl-1),mfzeta_x(lvl-1),loc_time,mfzeta_x(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfzeta_x(lvl),loc_time,mfzeta_x(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfzeta_y(lvl),time_prev(lvl-1),mfzeta_y(lvl-1),loc_time,mfzeta_y(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfzeta_y(lvl),loc_time,mfzeta_y(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfzeta_z(lvl),time_prev(lvl-1),mfzeta_z(lvl-1),loc_time,mfzeta_z(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfzeta_z(lvl),loc_time,mfzeta_z(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfpi_xx(lvl),time_prev(lvl-1),mfpi_xx(lvl-1),loc_time,mfpi_xx(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfpi_xx(lvl),loc_time,mfpi_xx(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfpi_yy(lvl),time_prev(lvl-1),mfpi_yy(lvl-1),loc_time,mfpi_yy(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfpi_yy(lvl),loc_time,mfpi_yy(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfpi_zz(lvl),time_prev(lvl-1),mfpi_zz(lvl-1),loc_time,mfpi_zz(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfpi_zz(lvl),loc_time,mfpi_zz(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfpi_xy(lvl),time_prev(lvl-1),mfpi_xy(lvl-1),loc_time,mfpi_xy(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfpi_xy(lvl),loc_time,mfpi_xy(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfpi_xz(lvl),time_prev(lvl-1),mfpi_xz(lvl-1),loc_time,mfpi_xz(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfpi_xz(lvl),loc_time,mfpi_xz(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mfpi_yz(lvl),time_prev(lvl-1),mfpi_yz(lvl-1),loc_time,mfpi_yz(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mfpi_yz(lvl),loc_time,mfpi_yz(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mflambda_x(lvl),time_prev(lvl-1),mflambda_x(lvl-1),loc_time,mflambda_x(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mflambda_x(lvl),loc_time,mflambda_x(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mflambda_y(lvl),time_prev(lvl-1),mflambda_y(lvl-1),loc_time,mflambda_y(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mflambda_y(lvl),loc_time,mflambda_y(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  call amrex_fillpatch(mflambda_z(lvl),time_prev(lvl-1),mflambda_z(lvl-1),loc_time,mflambda_z(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,time_prev(lvl),mflambda_z(lvl),loc_time,mflambda_z(lvl),amrex_geom(lvl),edge_bounds,&
!    loc_time,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
!            if (allocated(node_pairs)) deallocate(node_pairs)
!            allocate(node_pairs(27))
!
!            if (needs_interp(a,b,c)) then
!              node_pairs = .false.
!              n_paired_nodes = 0
!              do i = 2,27
!
!!                if (a == 58 .and. b == 67 .and. c == 5) then
!!                  write(*,*) 'adjacent nodes',state(a+cx_27(i),b+cy_27(i),c+cz_27(i),1),&
!!                    state(a+cx_27(opp_27(i)),b+cy_27(opp_27(i)),c+cz_27(opp_27(i)),1),&
!!                    needs_interp(a+cx_27(i),b+cy_27(i),c+cz_27(i)),&
!!                    needs_interp(a+cx_27(opp_27(i)),b+cy_27(opp_27(i)),c+cz_27(opp_27(i)))
!!                end if
!                if (state(a+cx_27(i),b+cy_27(i),c+cz_27(i),1) >= 0 .and. &
!                  state(a+cx_27(opp_27(i)),b+cy_27(opp_27(i)),c+cz_27(opp_27(i)),1) >=0 .and. &
!                  .not. needs_interp(a+cx_27(i),b+cy_27(i),c+cz_27(i)) .and. &
!                  .not. needs_interp(a+cx_27(opp_27(i)),b+cy_27(opp_27(i)),c+cz_27(opp_27(i)))) then
!!
!                  n_paired_nodes = n_paired_nodes + 1
!                  node_pairs(i) = .true.
!                end if
!              end do
!!
!!
!!
!              if (n_paired_nodes > 0) then
!
!                do i = 1,dir
!                  dir_sum_fi = 0.0D0
!                  dir_sum_gi = 0.0D0
!                  do j = 2,27
!                    if (node_pairs(j)) then
!                      dir_sum_fi = dir_sum_fi + fi(a+cx_27(j),b+cy_27(j),c+cz_27(j),i)
!                      dir_sum_gi = dir_sum_gi + gi(a+cx_27(j),b+cy_27(j),c+cz_27(j),i)
!                    end if
!                  end do
!
!                  fi(a,b,c,i) = dir_sum_fi/n_paired_nodes
!                  gi(a,b,c,i) = dir_sum_gi/n_paired_nodes
!
!                end do
!!
!!
!!
!              else
!                do i = 1,dir
!                  dir_sum_fi = 0.0D0
!                  dir_sum_gi = 0.0D0
!                  adj_count = 0
!                  do j = 2,7
!                    if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= 0 .and. &
!                      .not. needs_interp(a+cx_27(j),b+cy_27(j),c+cz_27(j))) then
!                        dir_sum_fi = dir_sum_fi + fi(a+cx_27(j),b+cy_27(j),c+cz_27(j),i)
!                        dir_sum_gi = dir_sum_gi + gi(a+cx_27(j),b+cy_27(j),c+cz_27(j),i)
!
!                        adj_count = adj_count+1
!                    end if
!                  end do
!
!                  if (dir_sum_fi < very_small_number .or. dir_sum_gi < very_small_number) then
!                    fi(a,b,c,i) = (fi(a,b,c,i) + crs_fi(a,b,c,i))/2.0D0
!                    gi(a,b,c,i) = (gi(a,b,c,i) + crs_gi(a,b,c,i))/2.0D0
!                  else
!                    fi(a,b,c,i) = dir_sum_fi/adj_count
!                    gi(a,b,c,i) = dir_sum_gi/adj_count
!                  end if
!
!                end do
!              end if
!            end if

!             if (a == 60 .and. b == 66 .and. c == 64) then
!                write(*,*) 'interlevel fis',fi(a,b,c,1:dir),a,b,c
!                write(*,*) 'interlevel gis',gi(a,b,c,1:dir),a,b,c
!                write(*,*) 'interlevel states',state(a,b,c,1:dir),a,b,c
!!                write(*,*) 'paired interp',n_paired_nodes,node_pairs
!!                write(*,*) 'interlevel lagrange',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!!                  pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!!                  pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!!                write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!!                  w_vel(a,b,c,1),temp(a,b,c,1)
!
!              end if

!              if (a == 57 .and. b == 67 .and. c == 5) then
!                write(*,*) 'interlevel fis',fi(a,b,c,1:dir),a,b,c
!                write(*,*) 'interlevel gis',gi(a,b,c,1:dir),a,b,c
!                write(*,*) 'interlevel states',state(a,b,c,1:dir),a,b,c
!                write(*,*) 'paired interp',n_paired_nodes,node_pairs
!!                write(*,*) 'interlevel lagrange',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!!                  pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!!                  pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!!                write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!!                  w_vel(a,b,c,1),temp(a,b,c,1)
!
!              end if

!            if (mod(a,2) /= 0 .or. mod(b,2) /= 0 .or. mod(c,2) /= 0) then













!            beta_1_c = 1.0D0/((2.0D0*mu_c)/(crs_rho(a,b,c,1)*crs_temp(a,b,c,1))+1)
!            beta_1_f = 1.0D0/((2.0D0*mu_f)/(rho(a,b,c,1)*temp(a,b,c,1))+1)
!            beta_2_c = 1.0D0/((2.0D0*kappa_c)/(crs_rho(a,b,c,1)*(const_vol+1)*crs_temp(a,b,c,1))+1)
!            beta_2_f = 1.0D0/((2.0D0*kappa_f)/(rho(a,b,c,1)*(const_vol+1)*temp(a,b,c,1))+1)
!            beta_1_c = beta_1_calculator(mu_c,crs_rho(a,b,c,1),crs_temp(a,b,c,1),lvl-1,timestep(amr_max_lvl))
!            beta_1_f = beta_1_calculator(mu_f,rho(a,b,c,1),temp(a,b,c,1),lvl,timestep(amr_max_lvl))
!            beta_2_c = beta_2_calculator(kappa_c,crs_rho(a,b,c,1),crs_temp(a,b,c,1),lvl-1,timestep(amr_max_lvl))
!            beta_2_f = beta_2_calculator(kappa_f,rho(a,b,c,1),temp(a,b,c,1),lvl,timestep(amr_max_lvl))
!
!            beta_1_rat = beta_1_c/(2.0D0*beta_1_f)
!            beta_2_rat = beta_2_c/(2.0D0*beta_2_f)
!
!            k_fineq = elbm_epsilon(lvl)*crs_rho(a,b,c,1)/(12.0D0*crs_temp(a,b,c,1)**2*beta_2_c) * &
!              (beta_1_rat-beta_2_rat)
!            k_gineq = elbm_epsilon(lvl)*crs_rho(a,b,c,1)*(2*const_vol-dimensions)/(2.0D0*beta_2_c) * &
!              (beta_1_rat + beta_2_rat)
!!
!            select case(second_central)
!!
!! central difference derivatives good for all directions
!!
!                case(7)
!                  dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                    crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                    .true.,.false.,.false.)
!                  dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                    crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                    .true.,.false.,.false.)
!                  dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                    crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                    .true.,.false.,.false.)
!!
!! x-direction central diff der bad
!!
!                case(6)
!                  dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                     crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1) /),&
!                    .true.,.false.,.false.)
!                  dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                    crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1) /),&
!                    .true.,.false.,.false.)
!! x-direction derivatives
!                  if (second_upwind == 1 .or. second_upwind == 3 .or. &
!                      second_upwind == 5 .or. second_upwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 1 .or. second_downwind == 3 .or. &
!                           second_downwind == 5 .or. second_downwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 1 .or. center == 3 .or. center == 5 .or. center == 7) then
!                    dTdx = (crs_temp(a+2,b,c,1)-crs_temp(a-2,b,c,1))/2.0D0
!                  else if (upwind == 1 .or. upwind == 3 .or. upwind == 5 .or. upwind == 7) then
!                    dTdx = crs_temp(a,b,c,1)-crs_temp(a-2,b,c,1)
!                  else if (downwind == 1 .or. downwind == 3 .or. downwind == 5 .or. downwind == 7) then
!                    dTdx = crs_temp(a+2,b,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdx = 0.0D0
!                  end if
!!
!! y-direction central diff der bad
!!
!                case(5)
!                  dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                     crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1) /),&
!                    .true.,.false.,.false.)
!                  dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                     crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1) /),&
!                    .true.,.false.,.false.)
!! y-direction derivatives
!                  if (center == 2 .or. center == 3 .or. center == 6 .or. center == 7) then
!                    dTdy = (crs_temp(a,b+2,c,1)-crs_temp(a,b-2,c,1))/2.0D0
!                  else if (second_upwind == 2 .or. second_upwind == 3 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                       crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1) /),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 2 .or. second_downwind == 3 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                      crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                      .false.,.false.,.true.)
!                  else if (upwind == 2 .or. upwind == 3 .or. upwind == 6 .or. upwind == 7) then
!                    dTdy = crs_temp(a,b,c,1)-crs_temp(a,b-2,c,1)
!                  else if (downwind == 2 .or. downwind == 3 .or. downwind == 6 .or. downwind == 7) then
!                    dTdy = crs_temp(a,b+2,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdy = 0.0D0
!                  end if
!!
!! x- and y-direction central diff ders bad
!!
!                case(4)
!                  dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                     crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                    .true.,.false.,.false.)
!! x-direction derivatives
!                  if (second_upwind == 1 .or. second_upwind == 3 .or. &
!                      second_upwind == 5 .or. second_upwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 1 .or. second_downwind == 3 .or. &
!                           second_downwind == 5 .or. second_downwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 1 .or. center == 3 .or. center == 5 .or. center == 7) then
!                    dTdx = (crs_temp(a+2,b,c,1)-crs_temp(a-2,b,c,1))/2.0D0
!                  else if (upwind == 1 .or. upwind == 3 .or. upwind == 5 .or. upwind == 7) then
!                    dTdx = crs_temp(a,b,c,1)-crs_temp(a-2,b,c,1)
!                  else if (downwind == 1 .or. downwind == 3 .or. downwind == 5 .or. downwind == 7) then
!                    dTdx = crs_temp(a+2,b,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdx = 0.0D0
!                  end if
!! y-direction derivatives
!                  if (second_upwind == 2 .or. second_upwind == 3 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                      crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 2 .or. second_downwind == 3 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                      crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 2 .or. center == 3 .or. center == 6 .or. center == 7) then
!                    dTdy = (crs_temp(a,b+2,c,1)-crs_temp(a,b-2,c,1))/2.0D0
!                  else if (upwind == 2 .or. upwind == 3 .or. upwind == 6 .or. upwind == 7) then
!                    dTdy = crs_temp(a,b,c,1)-crs_temp(a,b-2,c,1)
!                  else if (downwind == 2 .or. downwind == 3 .or. downwind == 6 .or. downwind == 7) then
!                    dTdy = crs_temp(a,b+2,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdy = 0.0D0
!                  end if
!!
!! z-direction central diff der bad
!!
!                case(3)
!                  dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                    crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                    .true.,.false.,.false.)
!                  dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                    crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                    .true.,.false.,.false.)
!! z-direction derivatives
!                  if (second_upwind == 4 .or. second_upwind == 5 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 4 .or. second_downwind == 5 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 4 .or. center == 5 .or. center == 6 .or. center == 7) then
!                    dTdz = (crs_temp(a,b,c+2,1)-crs_temp(a,b,c-2,1))/2.0D0
!                  else if (upwind == 4 .or. upwind == 5 .or. upwind == 6 .or. upwind == 7) then
!                    dTdz = crs_temp(a,b,c,1)-crs_temp(a,b,c-2,1)
!                  else if (downwind == 4 .or. downwind == 5 .or. downwind == 6 .or. downwind == 7) then
!                    dTdz = crs_temp(a,b,c+2,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdz = 0.0D0
!                  end if
!!
!! x- and z-direction central diff der bad
!!
!                case(2)
!                  dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                    crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                    .true.,.false.,.false.)
!! z-direction derivatives
!                  if (second_upwind == 4 .or. second_upwind == 5 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 4 .or. second_downwind == 5 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 4 .or. center == 5 .or. center == 6 .or. center == 7) then
!                    dTdz = (crs_temp(a,b,c+2,1)-crs_temp(a,b,c-2,1))/2.0D0
!                  else if (upwind == 4 .or. upwind == 5 .or. upwind == 6 .or. upwind == 7) then
!                    dTdz = crs_temp(a,b,c,1)-crs_temp(a,b,c-2,1)
!                  else if (downwind == 4 .or. downwind == 5 .or. downwind == 6 .or. downwind == 7) then
!                    dTdz = crs_temp(a,b,c+2,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdz = 0.0D0
!                  end if
!! x-direction derivatives
!                  if (second_upwind == 1 .or. second_upwind == 3 .or. &
!                      second_upwind == 5 .or. second_upwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 1 .or. second_downwind == 3 .or. &
!                           second_downwind == 5 .or. second_downwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                    .false.,.false.,.true.)
!                  else if (center == 1 .or. center == 3 .or. center == 5 .or. center == 7) then
!                    dTdx = (crs_temp(a+2,b,c,1)-crs_temp(a-2,b,c,1))/2.0D0
!                  else if (upwind == 1 .or. upwind == 3 .or. upwind == 5 .or. upwind == 7) then
!                    dTdx = crs_temp(a,b,c,1)-crs_temp(a-2,b,c,1)
!                  else if (downwind == 1 .or. downwind == 3 .or. downwind == 5 .or. downwind == 7) then
!                    dTdx = crs_temp(a+2,b,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdx = 0.0D0
!                  end if
!!
!! y- and z-direction central diff der bad
!!
!                case(1)
!                  dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                    crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                    .true.,.false.,.false.)
!! y-direction derivatives
!                  if (second_upwind == 2 .or. second_upwind == 3 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                      crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 2 .or. second_downwind == 3 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                      crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 2 .or. center == 3 .or. center == 6 .or. center == 7) then
!                    dTdy = (crs_temp(a,b+2,c,1)-crs_temp(a,b-2,c,1))/2.0D0
!                  else if (upwind == 2 .or. upwind == 3 .or. upwind == 6 .or. upwind == 7) then
!                    dTdy = crs_temp(a,b,c,1)-crs_temp(a,b-2,c,1)
!                  else if (downwind == 2 .or. downwind == 3 .or. downwind == 6 .or. downwind == 7) then
!                    dTdy = crs_temp(a,b+2,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdy = 0.0D0
!                  end if
!! z-direction derivatives
!                  if (second_upwind == 4 .or. second_upwind == 5 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 4 .or. second_downwind == 5 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 4 .or. center == 5 .or. center == 6 .or. center == 7) then
!                    dTdz = (crs_temp(a,b,c+2,1)-crs_temp(a,b,c-2,1))/2.0D0
!                  else if (upwind == 4 .or. upwind == 5 .or. upwind == 6 .or. upwind == 7) then
!                    dTdz = crs_temp(a,b,c,1)-crs_temp(a,b,c-2,1)
!                  else if (downwind == 4 .or. downwind == 5 .or. downwind == 6 .or. downwind == 7) then
!                    dTdz = crs_temp(a,b,c+2,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdz = 0.0D0
!                  end if
!!
!! all direction central diff ders bad
!!
!                case(0)
!! x_direction derivatives
!                  if (second_upwind == 1 .or. second_upwind == 3 .or. &
!                      second_upwind == 5 .or. second_upwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 1 .or. second_downwind == 3 .or. &
!                           second_downwind == 5 .or. second_downwind == 7) then
!                    dTdx = second_order_known((/crs_temp(a,b,c,1),crs_temp(a+2,b,c,1),crs_temp(a+4,b,c,1),&
!                      crs_temp(a-2,b,c,1),crs_temp(a-4,b,c,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 1 .or. center == 3 .or. center == 5 .or. center == 7) then
!                    dTdx = (crs_temp(a+2,b,c,1)-crs_temp(a-2,b,c,1))/2.0D0
!                  else if (upwind == 1 .or. upwind == 3 .or. upwind == 5 .or. upwind == 7) then
!                    dTdx = crs_temp(a,b,c,1)-crs_temp(a-2,b,c,1)
!                  else if (downwind == 1 .or. downwind == 3 .or. downwind == 5 .or. downwind == 7) then
!                    dTdx = crs_temp(a+2,b,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdx = 0.0D0
!                  end if
!! y-direction derivatives
!                  if (second_upwind == 2 .or. second_upwind == 3 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                      crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 2 .or. second_downwind == 3 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdy = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b+2,c,1),crs_temp(a,b+4,c,1),&
!                    crs_temp(a,b-2,c,1),crs_temp(a,b-4,c,1)/),&
!                    .false.,.false.,.true.)
!                  else if (center == 2 .or. center == 3 .or. center == 6 .or. center == 7) then
!                    dTdy = (crs_temp(a,b+2,c,1)-crs_temp(a,b-2,c,1))/2.0D0
!                  else if (upwind == 2 .or. upwind == 3 .or. upwind == 6 .or. upwind == 7) then
!                    dTdy = crs_temp(a,b,c,1)-crs_temp(a,b-2,c,1)
!                  else if (downwind == 2 .or. downwind == 3 .or. downwind == 6 .or. downwind == 7) then
!                    dTdy = crs_temp(a,b+2,c,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdy = 0.0D0
!                  end if
!! z-direction derivatives
!                  if (second_upwind == 4 .or. second_upwind == 5 .or. second_upwind == 6 .or. &
!                    second_upwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.true.,.false.)
!                  else if (second_downwind == 4 .or. second_downwind == 5 .or. second_downwind == 6 .or. &
!                    second_downwind == 7) then
!                    dTdz = second_order_known((/crs_temp(a,b,c,1),crs_temp(a,b,c+2,1),crs_temp(a,b,c+4,1),&
!                      crs_temp(a,b,c-2,1),crs_temp(a,b,c-4,1)/),&
!                      .false.,.false.,.true.)
!                  else if (center == 4 .or. center == 5 .or. center == 6 .or. center == 7) then
!                    dTdz = (crs_temp(a,b,c+2,1)-crs_temp(a,b,c-2,1))/2.0D0
!                  else if (upwind == 4 .or. upwind == 5 .or. upwind == 6 .or. upwind == 7) then
!                    dTdz = crs_temp(a,b,c,1)-crs_temp(a,b,c-2,1)
!                  else if (downwind == 4 .or. downwind == 5 .or. downwind == 6 .or. downwind == 7) then
!                    dTdz = crs_temp(a,b,c+2,1)-crs_temp(a,b,c,1)
!                  else
!                    dTdz = 0.0D0
!                  end if
!! The shit's hit the fan
!                case default
!                  dTdx = 0.0D0
!                  dTdy = 0.0D0
!                  dTdz = 0.0D0
!              end select
!!
!!
!! Go through all directions for interlevel interpolations
!            !write(*,*) 'temperature derivatives',dTdx,dTdy,dTdz
!!            if (lvl == 3 .and. a == 223 .and. b == 128 .and. c == 124) then
!!              write(*,*) 'ghost node info',state(a,b,c,1),state(a,b,c,2),fi(a,b,c,1),fi(a,b,c,2),fi(a,b,c,5),&
!!                fi(a,b,c,16),fi(a,b,c,19),fi(a,b,c,34),fi(a,b,c,37),u_vel(a,b,c,1),rho(a,b,c,1),lvl,a,b,c,self
!!              write(*,*) 'fluid node info',state(a+2,b,c,1),state(a+2,b,c,2),fi(a+2,b,c,1),fi(a+2,b,c,2),fi(a+2,b,c,5),&
!!                fi(a+2,b,c,16),fi(a+2,b,c,19),fi(a+2,b,c,34),fi(a+2,b,c,37),u_vel(a+2,b,c,1),rho(a+2,b,c,1),lvl,a+2,b,c,self
!!              write(*,*) 'coarse values',crs_u(a,b,c,1),crs_u(a+2,b,c,1),crs_v(a,b,c,1),crs_v(a+2,b,c,1),&
!!                crs_w(a,b,c,1),crs_w(a+2,b,c,1),crs_temp(a,b,c,1),crs_temp(a+2,b,c,1)
!!              write(*,*) 'local values',u_vel(a,b,c,1),u_vel(a+2,b,c,1),v_vel(a,b,c,1),v_vel(a+2,b,c,1),&
!!                w_vel(a,b,c,1),w_vel(a+2,b,c,1),temp(a,b,c,1),temp(a+2,b,c,1)
!!            end if
!
!
!
!              select case (dir)
!                case(19)
!
!                do i = 1,dir
!                  crs_fieq = fieq_comp_19(crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
!                    crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
!                    crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
!                    crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1),i)
!
!                  crs_gieq = (2*const_vol-dimensions)*crs_temp(a,b,c,1)*crs_fieq
!
!                  loc_fieq = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!                    zeta_y(a,b,c,1),zeta_z(a,b,c,1),pi_xx(a,b,c,1),pi_yy(a,b,c,1),&
!                    pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz(a,b,c,1),&
!                    lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
!
!                  loc_gieq = (2*const_vol-dimensions)*temp(a,b,c,1)*loc_fieq
!
!                  select case (i)
!                    case(1)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(2)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdx - &
!                        6.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*dTdx
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(3)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdy - &
!                        6.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*dTdy
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(4)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdz - &
!                        6.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*dTdz
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(5)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdx - &
!                        6.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        k_gineq*dTdx
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(6)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdy - &
!                        6.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        k_gineq*dTdy
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(7)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdz - &
!                        6.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        k_gineq*dTdz
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(8)
!                      temp_der_sum = dTdx+dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        2.0D0*temp_der_sum - &
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(9)
!                      temp_der_sum = dTdy-dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(10)
!                      temp_der_sum = -dTdx+dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(11)
!                      temp_der_sum = dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(12)
!                      temp_der_sum = dTdx+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(13)
!                      temp_der_sum = dTdx-dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(14)
!                      temp_der_sum = dTdx+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(15)
!                      temp_der_sum = -dTdx+dTdz
!                      loc_fineq =  beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(16)
!                      temp_der_sum = dTdx-dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(17)
!                      temp_der_sum = dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(18)
!                      temp_der_sum = dTdx+dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq * &
!                        ((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!                    case(19)
!                      temp_der_sum = -dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((4.0D0-12.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        2.0D0*temp_der_sum-&
!                        3.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = loc_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = loc_gieq+ loc_gineq
!
!                  end select
!                end do
!                case(39)
!!                  if (shifted) then
!!
!!                    end do
!!                  else
!!
!!                    end do
!!                  end if
!!
!!
!!
!!
!!
!!                  if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0 .or. fi(a,b,c,i) > 1.0D0 &
!!                    .or. gi(a,b,c,i) > 1.0D0) then
!!                    write(*,*) 'fi or gi problem in interlevel interp',fi(a,b,c,i),gi(a,b,c,i),&
!!                      crs_fi(a,b,c,i),crs_gi(a,b,c,i),a,b,c,i,&
!!                      lvl,ducks%lo,ducks%hi
!!                  end if
!!!
!!                  if (crs_fieq < 0.0D0 .or. crs_gieq < 0.0D0 .or. crs_fieq > 1.0D0 &
!!                    .or. crs_gieq > 1.0D0) then
!!                    write(*,*) 'coarse fi or gi problem in interlevel interp',crs_fieq,crs_gieq,a,b,c,i,lvl,&
!!                        crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
!!                        crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
!!                        crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
!!                        crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1)
!!                  end if
!!!
!!                  if (loc_fieq < 0.0D0 .or. loc_gieq < 0.0D0 .or. loc_fieq > 1.0D0 &
!!                    .or. loc_gieq > 1.0D0) then
!!                    write(*,*) 'local fieq or gieq problem in interlevel interp',loc_fieq,loc_gieq,a,b,c,i,lvl,&
!!                        rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),pi_xx(a,b,c,1),pi_yy(a,b,c,1),&
!!                        pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz(a,b,c,1),&
!!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!!                  end if
!!
!!
!! Decide whether its the small or large lattice
!!
!!
!                  if (.not. shifted) then
!
!                  do i = 1,dir
!                      crs_fieq = fieq_comp_39(crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
!                        crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
!                        crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
!                        crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1),i)
!
!                      crs_gieq = (2*const_vol-dimensions)*crs_temp(a,b,c,1)*crs_fieq
!
!                      loc_fieq = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),pi_xx(a,b,c,1),pi_yy(a,b,c,1),&
!                        pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz(a,b,c,1),&
!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
!
!                      loc_gieq = (2*const_vol-dimensions)*temp(a,b,c,1)*loc_fieq
!
!                    select case (i)
!!
!                    case(1)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!
!!
!                    case(2)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq)+&
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdx-&
!                        6.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*dTdx
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(3)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdy-6.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*dTdy
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(4)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdz-6.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*dTdz
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(5)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdx-6.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - k_gineq*dTdx
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(6)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdy-6.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - k_gineq*dTdy
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(7)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((1.0D0-3.0D0*crs_temp(a,b,c,1))*dTdz-6.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - k_gineq*dTdz
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(8)
!                      temp_der_sum = dTdx+dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(9)
!                      temp_der_sum = -dTdx+dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(10)
!                      temp_der_sum = -dTdx-dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(11)
!                      temp_der_sum = dTdx-dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(12)
!                      temp_der_sum = dTdx+dTdy-dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum+ &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(13)
!                      temp_der_sum = -dTdx+dTdy-dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(14)
!                      temp_der_sum = dTdx+dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)- &
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(15)
!                      temp_der_sum = dTdx+dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((5.0D0-15.0D0*crs_temp(a,b,c,1))*temp_der_sum + &
!                        4.0D0*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(16)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((8.0D0-6.0D0*crs_temp(a,b,c,1))*dTdx-12.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        2.0D0*k_gineq*dTdx
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(17)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((8.0D0-6.0D0*crs_temp(a,b,c,1))*dTdy-12.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        2.0D0*k_gineq*dTdy
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(18)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((8.0D0-6.0D0*crs_temp(a,b,c,1))*dTdz-12.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        2.0D0*k_gineq*dTdz
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(19)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((8.0D0-6.0D0*crs_temp(a,b,c,1))*dTdx-12.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        2.0D0*k_gineq*dTdx
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(20)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((8.0D0-6.0D0*crs_temp(a,b,c,1))*dTdy-12.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        2.0D0*k_gineq*dTdy
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(21)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((8.0D0-6.0D0*crs_temp(a,b,c,1))*dTdz-12.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)- &
!                        2.0D0*k_gineq*dTdz
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(22)
!                      temp_der_sum = dTdx+dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(23)
!                      temp_der_sum = dTdx+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(24)
!                      temp_der_sum = dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(25)
!                      temp_der_sum = -dTdx+dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(26)
!                      temp_der_sum = dTdx+dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(27)
!                      temp_der_sum = dTdx-dTdy
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(28)
!                      temp_der_sum = dTdx-dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(29)
!                      temp_der_sum = dTdx+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq * ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(30)
!                      temp_der_sum = -dTdx+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq * ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(31)
!                      temp_der_sum = dTdy-dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(32)
!                      temp_der_sum = dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*&
!                        ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(33)
!                      temp_der_sum = -dTdy+dTdz
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((32.0D0-24.0D0*crs_temp(a,b,c,1))*temp_der_sum+&
!                        16.0D0*temp_der_sum-&
!                        6.0D0*crs_temp(a,b,c,1)*temp_der_sum)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        2.0D0*k_gineq*temp_der_sum
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(34)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((27.0D0-9.0D0*crs_temp(a,b,c,1))*dTdx-18.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        3.0D0*k_gineq*dTdx
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(35)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((27.0D0-9.0D0*crs_temp(a,b,c,1))*dTdy-18.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        3.0D0*k_gineq*dTdy
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(36)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*&
!                        ((27.0D0-9.0D0*crs_temp(a,b,c,1))*dTdz-18.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        3.0D0*k_gineq*dTdz
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(37)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((27.0D0-9.0D0*crs_temp(a,b,c,1))*dTdx-18.0D0*crs_temp(a,b,c,1)*dTdx)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        3.0D0*k_gineq*dTdx
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(38)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((27.0D0-9.0D0*crs_temp(a,b,c,1))*dTdy-18.0D0*crs_temp(a,b,c,1)*dTdy)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        3.0D0*k_gineq*dTdy
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(39)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*((27.0D0-9.0D0*crs_temp(a,b,c,1))*dTdz-18.0D0*crs_temp(a,b,c,1)*dTdz)
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)-&
!                        3.0D0*k_gineq*dTdz
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!
!                    end select
!
!                    end do
!                  else
!
!                    qabc = fill_qone_interlevel(dTdx,dTdy,dTdz)
!
!                    do i = 1,dir
!                      crs_fieq = fieq_comp_39_shifted(crs_rho(a,b,c,1),crs_chi(a,b,c,1),crs_zeta_x(a,b,c,1),&
!                        crs_zeta_y(a,b,c,1),crs_zeta_z(a,b,c,1),crs_pi_xx(a,b,c,1),crs_pi_yy(a,b,c,1),&
!                        crs_pi_zz(a,b,c,1),crs_pi_xy(a,b,c,1),crs_pi_xz(a,b,c,1),crs_pi_yz(a,b,c,1),&
!                        crs_lambda_x(a,b,c,1),crs_lambda_y(a,b,c,1),crs_lambda_z(a,b,c,1),i)
!
!                      crs_gieq = (2*const_vol-dimensions)*crs_temp(a,b,c,1)*crs_fieq
!
!                      loc_fieq = fieq_comp_39_shifted(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),pi_xx(a,b,c,1),pi_yy(a,b,c,1),&
!                        pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz(a,b,c,1),&
!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
!
!                      loc_gieq = (2*const_vol-dimensions)*temp(a,b,c,1)*loc_fieq
!
!                    select case (i)
!!
!                    case(1)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq)+&
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(2)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq)+&
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(3)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(4)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(5)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(6)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(7)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) - &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) - &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(8)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(9)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(10)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(11)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(12)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(13)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(14)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(15)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(16)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(17)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(18)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(19)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+loc_gineq
!!
!                    case(20)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(21)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+ &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(22)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(23)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(24)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(25)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(26)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(27)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(28)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(29)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq * (sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(30)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq *(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(31)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(32)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(33)
!
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(34)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(35)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(36)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+&
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(37)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(38)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq)+ &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!!
!                    case(39)
!                      loc_fineq = beta_1_rat*(crs_fi(a,b,c,i)-crs_fieq) + &
!                        k_fineq*(sum(yi(:,i)*qabc) - 3*crs_temp(a,b,c,1)*sum(kroen_val(:,i)*qabc))
!                      fi(a,b,c,i) = crs_fieq+ loc_fineq
!
!                      loc_gineq = beta_1_rat*(crs_gi(a,b,c,i)-crs_gieq) + &
!                        k_gineq*(cx(i)*dTdx+cy(i)*dTdy+cz(i)*dTdz)
!                      gi(a,b,c,i) = crs_gieq+ loc_gineq
!
!                    end select
!
!                    end do
!
!                  end if
!
!              end select
!!
!!
!
!            end if
!          end if
!
!        end do
!      end do
!    end do
!    do c = ducks%lo(3)-nghosts_mid,ducks%hi(3)+nghosts_mid
!      do b = ducks%lo(2)-nghosts_mid,ducks%hi(2)+nghosts_mid
!        do a = ducks%lo(1)-nghosts_mid,ducks%hi(1)+nghosts_mid
!          if (state(a,b,c,1) >=10000000 .and. state(a,b,c,1) <= 1100000000) then !checks if a true ghost node
!          if (mod(a,2) == 0 .and. mod(b,2) == 0 .and. mod(c,2) /= 0) then
!!
!! Edge in z
!!
!            if (state(a,b,c,7) >= 0 .and. state(a,b,c,39) >= 0) then
!              z_up_clear = .true.
!              z_up_close = .true.
!            else if (state(a,b,c,7) < 0) then
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                !write(*,*) 'snarf +z'
!                z_up_close = .false.
!                z_up_clear = .false.
!
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!
!!              if (a == 244 .and. b == 130 .and. c == 145) then
!!                write(*,*) 'checking fis for ghost nodex',fi(a,b,c,1:dir),a,b,c,self
!!                write(*,*) 'checking gis for ghost node',gi(a,b,c,1:dir),a,b,c,self
!!                write(*,*) 'weird states for ghost nodes',state(a,b,c,1:dir),a,b,c,self
!!              end if
!              cycle
!            else
!              z_up_clear = .false.
!              z_up_close = .true.
!            end if
!
!            if (state(a,b,c,4) >= 0 .and. state(a,b,c,36) >= 0) then
!              z_down_clear = .true.
!              z_down_close = .true.
!            else if (state(a,b,c,4) < 0) then
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                !write(*,*) 'snarf -z'
!
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!              z_down_close = .false.
!              z_down_clear = .false.
!
!              cycle
!            else
!              z_down_clear = .false.
!              z_down_close = .true.
!            end if
!!
!            if (z_up_clear .and. z_down_clear ) then
!              do i = 1,dir
!                fi(a,b,c,i) = (9.0D0*(fi(a,b,c+1,i)+fi(a,b,c-1,i))-&
!                                 fi(a,b,c-3,i)-fi(a,b,c+3,i))/16.0D0
!                gi(a,b,c,i) = (9.0D0*(gi(a,b,c+1,i)+gi(a,b,c-1,i))-&
!                                 gi(a,b,c+3,i)-gi(a,b,c-3,i))/16.0D0
!              end do
!
!
!            else if (z_up_clear .and. z_down_close .and. .not. z_down_clear) then
!              do i = 1,dir
!                fi(a,b,c,i) = 3.0D0/8.0D0*fi(a,b,c+1,i)+3.0D0/4.0D0*fi(a,b,c-1,i)-fi(a,b,c-3,i)/8.0D0
!                gi(a,b,c,i) = 3.0D0/8.0D0*gi(a,b,c+1,i)+3.0D0/4.0D0*gi(a,b,c-1,i)-gi(a,b,c-3,i)/8.0D0
!              end do
!!            else if (z_down_clear .and. z_up_close .and. .not. z_up_clear) then
!              do i = 1,dir
!                fi(a,b,c,i) = 3.0D0/8.0D0*fi(a,b,c-1,i)+3.0D0/4.0D0*fi(a,b,c+1,i)-fi(a,b,c+3,i)/8.0D0
!                gi(a,b,c,i) = 3.0D0/8.0D0*gi(a,b,c-1,i)+3.0D0/4.0D0*gi(a,b,c+1,i)-gi(a,b,c+3,i)/8.0D0
!              end do
!!              chi(a,b,c,1) =      3.0D0/8.0D0*chi(a,b,c-1,1)+3.0D0/4.0D0*chi(a,b,c+1,1)-chi(a,b,c+3,1)/8.0D0
!!              zeta_x(a,b,c,1) =   3.0D0/8.0D0*zeta_x(a,b,c-1,1)+3.0D0/4.0D0*zeta_x(a,b,c+1,1)-zeta_x(a,b,c+3,1)/8.0D0
!!              zeta_y(a,b,c,1) =   3.0D0/8.0D0*zeta_y(a,b,c-1,1)+3.0D0/4.0D0*zeta_y(a,b,c+1,1)-zeta_y(a,b,c+3,1)/8.0D0
!!              zeta_z(a,b,c,1) =   3.0D0/8.0D0*zeta_z(a,b,c-1,1)+3.0D0/4.0D0*zeta_z(a,b,c+1,1)-zeta_z(a,b,c+3,1)/8.0D0
!!              pi_xx(a,b,c,1) =    3.0D0/8.0D0*pi_xx(a,b,c-1,1)+3.0D0/4.0D0*pi_xx(a,b,c+1,1)-pi_xx(a,b,c+3,1)/8.0D0
!!              pi_yy(a,b,c,1) =    3.0D0/8.0D0*pi_yy(a,b,c-1,1)+3.0D0/4.0D0*pi_yy(a,b,c+1,1)-pi_yy(a,b,c+3,1)/8.0D0
!!              pi_zz(a,b,c,1) =    3.0D0/8.0D0*pi_zz(a,b,c-1,1)+3.0D0/4.0D0*pi_zz(a,b,c+1,1)-pi_zz(a,b,c+3,1)/8.0D0
!!              pi_xy(a,b,c,1) =    3.0D0/8.0D0*pi_xy(a,b,c-1,1)+3.0D0/4.0D0*pi_xy(a,b,c+1,1)-pi_xy(a,b,c+3,1)/8.0D0
!!              pi_xz(a,b,c,1) =    3.0D0/8.0D0*pi_xz(a,b,c-1,1)+3.0D0/4.0D0*pi_xz(a,b,c+1,1)-pi_xz(a,b,c+3,1)/8.0D0
!!              pi_yz(a,b,c,1) =    3.0D0/8.0D0*pi_yz(a,b,c-1,1)+3.0D0/4.0D0*pi_yz(a,b,c+1,1)-pi_yz(a,b,c+3,1)/8.0D0
!!              lambda_x(a,b,c,1) = 3.0D0/8.0D0*lambda_x(a,b,c-1,1)+3.0D0/4.0D0*lambda_x(a,b,c+1,1)-lambda_x(a,b,c+3,1)/8.0D0
!!              lambda_y(a,b,c,1) = 3.0D0/8.0D0*lambda_y(a,b,c-1,1)+3.0D0/4.0D0*lambda_y(a,b,c+1,1)-lambda_y(a,b,c+3,1)/8.0D0
!!              lambda_z(a,b,c,1) = 3.0D0/8.0D0*lambda_z(a,b,c-1,1)+3.0D0/4.0D0*lambda_z(a,b,c+1,1)-lambda_z(a,b,c+3,1)/8.0D0
!            else if (z_up_close .and. z_down_close) then
!              do i=1,dir
!                fi(a,b,c,i) = (fi(a,b,c+1,i)+fi(a,b,c-1,i))/2.0D0
!                gi(a,b,c,i) = (gi(a,b,c+1,i)+gi(a,b,c-1,i))/2.0D0
!              end do
!!              chi(a,b,c,1) = (chi(a,b,c+1,1)+chi(a,b,c-1,1))/2.0D0
!!              zeta_x(a,b,c,1) = (zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))/2.0D0
!!              zeta_y(a,b,c,1) = (zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))/2.0D0
!!              zeta_z(a,b,c,1) = (zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))/2.0D0
!!              pi_xx(a,b,c,1) =(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))/2.0D0
!!              pi_yy(a,b,c,1) =(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))/2.0D0
!!              pi_zz(a,b,c,1) =(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))/2.0D0
!!              pi_xy(a,b,c,1) =(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))/2.0D0
!!              pi_xz(a,b,c,1) =(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))/2.0D0
!!              pi_yz(a,b,c,1) =(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))/2.0D0
!!              lambda_x(a,b,c,1) = (lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))/2.0D0
!!              lambda_y(a,b,c,1) = (lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))/2.0D0
!!              lambda_z(a,b,c,1) = (lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))/2.0D0
!            else
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!            end if
!!            if (a == 244 .and. b == 130 .and. c == 145) then
!!              write(*,*) 'checking fis for ghost nodez',fi(a,b,c,1:dir),a,b,c
!!              write(*,*) 'checking gis for ghost node',gi(a,b,c,1:dir),a,b,c
!!            end if
!
!!
!! Edge in y
!!
!          else if (mod(a,2) == 0 .and. mod(c,2) == 0 .and. mod(b,2) /= 0) then
!            if (state(a,b,c,6) >= 0 .and. state(a,b,c,38) >= 0) then
!              y_up_clear = .true.
!              y_up_close = .true.
!            else if (state(a,b,c,6) < 0) then
!              y_up_clear = .false.
!              y_up_close = .false.
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                !write(*,*) 'snarf +y'
!
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!              cycle
!            else
!              y_up_close = .true.
!              y_up_clear = .false.
!            end if
!            if (state(a,b,c,3) >= 0 .and. state(a,b,c,35) >= 0) then
!              y_down_clear = .true.
!              y_down_close = .true.
!            else if (state(a,b,c,3) < 0) then
!              y_down_close = .false.
!              y_down_clear = .false.
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                !write(*,*) 'snarf -y'
!
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!              cycle
!            else
!              y_down_close = .true.
!              y_down_clear = .false.
!            end if
!!
!            if (y_up_clear .and. y_down_clear ) then
!              do i = 1,dir
!                fi(a,b,c,i) = (9.0D0*(fi(a,b+1,c,i)+fi(a,b-1,c,i))-&
!                                 fi(a,b-3,c,i)-fi(a,b+3,c,i))/16.0D0
!                gi(a,b,c,i) = (9.0D0*(gi(a,b+1,c,i)+gi(a,b-1,c,i))-&
!                                 gi(a,b+3,c,i)-gi(a,b-3,c,i))/16.0D0
!              end do
!            else if (y_up_clear .and. y_down_close .and. .not. y_down_clear) then
!              do i = 1,dir
!                fi(a,b,c,i) = 3.0D0/8.0D0*fi(a,b+1,c,i)+3.0D0/4.0D0*fi(a,b-1,c,i)-fi(a,b-3,c,i)/8.0D0
!                gi(a,b,c,i) = 3.0D0/8.0D0*gi(a,b+1,c,i)+3.0D0/4.0D0*gi(a,b-1,c,i)-gi(a,b-3,c,i)/8.0D0
!              end do
!
!            else if (y_down_clear .and. y_up_close .and. .not. y_up_clear) then
!              do i = 1,dir
!                fi(a,b,c,i) = 3.0D0/8.0D0*fi(a,b-1,c,i)+3.0D0/4.0D0*fi(a,b+1,c,i)-fi(a,b+3,c,i)/8.0D0
!                gi(a,b,c,i) = 3.0D0/8.0D0*gi(a,b-1,c,i)+3.0D0/4.0D0*gi(a,b+1,c,i)-gi(a,b+3,c,i)/8.0D0
!              end do
!            else if (y_down_close .and. y_up_close) then
!              do i=1,dir
!                fi(a,b,c,i) = (fi(a,b+1,c,i)+fi(a,b-1,c,i))/2.0D0
!                gi(a,b,c,i) = (gi(a,b+1,c,i)+gi(a,b-1,c,i))/2.0D0
!              end do

!            else
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!            end if
!!
!! Edge in x
!!
!          else if (mod(c,2) == 0 .and. mod(b,2) == 0 .and. mod(a,2) /= 0) then
!            if (state(a,b,c,5) >= 0 .and. state(a,b,c,37) >= 0) then
!              x_up_clear = .true.
!              x_up_close = .true.
!            else if (state(a,b,c,5) < 0) then
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!               !write(*,*) 'snarf +x'
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!              cycle
!              x_up_close = .false.
!              x_up_clear = .false.
!            else
!              x_up_close = .true.
!              x_up_clear = .false.
!            end if
!
!            if (state(a,b,c,2) >= 0 .and. state(a,b,c,34) >= 0) then
!              x_down_clear = .true.
!              x_down_close = .true.
!            else if (state(a,b,c,2) < 0) then
!              x_down_close = .false.
!              x_down_clear = .false.
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                !write(*,*) 'snarf -x'
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!              cycle
!            else
!              x_down_close = .true.
!              x_down_clear = .false.
!            end if
!!
!            if (x_up_clear .and. x_down_clear ) then
!              do i = 1,dir
!                fi(a,b,c,i) = (9.0D0*(fi(a+1,b,c,i)+fi(a-1,b,c,i))-&
!                                 fi(a-3,b,c,i)-fi(a+3,b,c,i))/16.0D0
!                gi(a,b,c,i) = (9.0D0*(gi(a+1,b,c,i)+gi(a-1,b,c,i))-&
!                                 gi(a+3,b,c,i)-gi(a-3,b,c,i))/16.0D0
!              end do
!!              chi(a,b,c,1) = (9.0D0*(chi(a+1,b,c,1)+chi(a-1,b,c,1))-&
!!                chi(a+3,b,c,1)-chi(a-3,b,c,1))/16.0D0
!!              zeta_x(a,b,c,1) = (9.0D0*(zeta_x(a+1,b,c,1)+zeta_x(a-1,b,c,1))-&
!!                zeta_x(a+3,b,c,1)-zeta_x(a-3,b,c,1))/16.0D0
!!              zeta_y(a,b,c,1) = (9.0D0*(zeta_y(a+1,b,c,1)+zeta_y(a-1,b,c,1))-&
!!                zeta_y(a+3,b,c,1)-zeta_y(a-3,b,c,1))/16.0D0
!!              zeta_z(a,b,c,1) = (9.0D0*(zeta_z(a+1,b,c,1)+zeta_z(a-1,b,c,1))-&
!!                zeta_z(a+3,b,c,1)-zeta_z(a-3,b,c,1))/16.0D0
!!              pi_xx(a,b,c,1) = (9.0D0*(pi_xx(a+1,b,c,1)+pi_xx(a-1,b,c,1))-&
!!                pi_xx(a+3,b,c,1)-pi_xx(a-3,b,c,1))/16.0D0
!!              pi_yy(a,b,c,1) = (9.0D0*(pi_yy(a+1,b,c,1)+pi_yy(a-1,b,c,1))-&
!!                pi_yy(a-3,b,c,1)-pi_yy(a+3,b,c,1))/16.0D0
!!              pi_zz(a,b,c,1) = (9.0D0*(pi_zz(a+1,b,c,1)+pi_zz(a-1,b,c,1))-&
!!                pi_zz(a+3,b,c,1)-pi_zz(a-3,b,c,1))/16.0D0
!!              pi_xy(a,b,c,1) = (9.0D0*(pi_xy(a+1,b,c,1)+pi_xy(a-1,b,c,1))-&
!!                pi_xy(a+3,b,c,1)-pi_xy(a-3,b,c,1))/16.0D0
!!              pi_xz(a,b,c,1) = (9.0D0*(pi_xz(a+1,b,c,1)+pi_xz(a-1,b,c,1))-&
!!                pi_xz(a+3,b,c,1)-pi_xz(a-3,b,c,1))/16.0D0
!!              pi_yz(a,b,c,1) = (9.0D0*(pi_yz(a+1,b,c,1)+pi_yz(a-1,b,c,1))-&
!!                pi_yz(a+3,b,c,1)-pi_yz(a-3,b,c,1))/16.0D0
!!              lambda_x(a,b,c,1) = (9.0D0*(lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))-&
!!                lambda_x(a+3,b,c,1)-lambda_x(a-3,b,c,1))/16.0D0
!!              lambda_y(a,b,c,1) = (9.0D0*(lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))-&
!!                lambda_y(a+3,b,c,1)-lambda_y(a-3,b,c,1))/16.0D0
!!              lambda_z(a,b,c,1) = (9.0D0*(lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))-&
!!                lambda_z(a+3,b,c,1)-lambda_z(a-3,b,c,1))/16.0D0
!            else if (x_up_clear .and. x_down_close .and. .not. x_down_clear) then
!              do i = 1,dir
!                fi(a,b,c,i) = 3.0D0/8.0D0*fi(a+1,b,c,i)+3.0D0/4.0D0*fi(a-1,b,c,i)-fi(a-3,b,c,i)/8.0D0
!                gi(a,b,c,i) = 3.0D0/8.0D0*gi(a+1,b,c,i)+3.0D0/4.0D0*gi(a-1,b,c,i)-gi(a-3,b,c,i)/8.0D0
!              end do
!!              chi(a,b,c,1)      = 3.0D0/8.0D0*chi(a+1,b,c,1)     +3.0D0/4.0D0*chi(a-1,b,c,1)     -chi(a-3,b,c,1)/8.0D0
!!              zeta_x(a,b,c,1)   = 3.0D0/8.0D0*zeta_x(a+1,b,c,1)  +3.0D0/4.0D0*zeta_x(a-1,b,c,1)  -zeta_x(a-3,b,c,1)/8.0D0
!!              zeta_y(a,b,c,1)   = 3.0D0/8.0D0*zeta_y(a+1,b,c,1)  +3.0D0/4.0D0*zeta_y(a-1,b,c,1)  -zeta_y(a-3,b,c,1)/8.0D0
!!              zeta_z(a,b,c,1)   = 3.0D0/8.0D0*zeta_z(a+1,b,c,1)  +3.0D0/4.0D0*zeta_z(a-1,b,c,1)  -zeta_z(a-3,b,c,1)/8.0D0
!!              pi_xx(a,b,c,1)    = 3.0D0/8.0D0*pi_xx(a+1,b,c,1)   +3.0D0/4.0D0*pi_xx(a-1,b,c,1)   -pi_xx(a-3,b,c,1)/8.0D0
!!              pi_yy(a,b,c,1)    = 3.0D0/8.0D0*pi_yy(a+1,b,c,1)   +3.0D0/4.0D0*pi_yy(a-1,b,c,1)   -pi_yy(a-3,b,c,1)/8.0D0
!!              pi_zz(a,b,c,1)    = 3.0D0/8.0D0*pi_zz(a+1,b,c,1)   +3.0D0/4.0D0*pi_zz(a-1,b,c,1)   -pi_zz(a-3,b,c,1)/8.0D0
!!              pi_xy(a,b,c,1)    = 3.0D0/8.0D0*pi_xy(a+1,b,c,1)   +3.0D0/4.0D0*pi_xy(a-1,b,c,1)   -pi_xy(a-3,b,c,1)/8.0D0
!!              pi_xz(a,b,c,1)    = 3.0D0/8.0D0*pi_xz(a+1,b,c,1)   +3.0D0/4.0D0*pi_xz(a-1,b,c,1)   -pi_xz(a-3,b,c,1)/8.0D0
!!              pi_yz(a,b,c,1)    = 3.0D0/8.0D0*pi_yz(a+1,b,c,1)   +3.0D0/4.0D0*pi_yz(a-1,b,c,1)   -pi_yz(a-3,b,c,1)/8.0D0
!!              lambda_x(a,b,c,1) = 3.0D0/8.0D0*lambda_x(a+1,b,c,1)+3.0D0/4.0D0*lambda_x(a-1,b,c,1)-lambda_x(a-3,b,c,1)/8.0D0
!!              lambda_y(a,b,c,1) = 3.0D0/8.0D0*lambda_y(a+1,b,c,1)+3.0D0/4.0D0*lambda_y(a-1,b,c,1)-lambda_y(a-3,b,c,1)/8.0D0
!!              lambda_z(a,b,c,1) = 3.0D0/8.0D0*lambda_z(a+1,b,c,1)+3.0D0/4.0D0*lambda_z(a-1,b,c,1)-lambda_z(a-3,b,c,1)/8.0D0
!            else if (x_down_clear .and. x_up_close .and. .not. x_up_clear) then
!              do i = 1,dir
!                fi(a,b,c,i) = 3.0D0/8.0D0*fi(a-1,b,c,i)+3.0D0/4.0D0*fi(a+1,b,c,i)-fi(a+3,b,c,i)/8.0D0
!                gi(a,b,c,i) = 3.0D0/8.0D0*gi(a-1,b,c,i)+3.0D0/4.0D0*gi(a+1,b,c,i)-gi(a+3,b,c,i)/8.0D0
!              end do
!!              chi(a,b,c,1)      = 3.0D0/8.0D0*chi(a-1,b,c,1)     +3.0D0/4.0D0*chi(a+1,b,c,1)     -chi(a+3,b,c,1)/8.0D0
!!              zeta_x(a,b,c,1)   = 3.0D0/8.0D0*zeta_x(a-1,b,c,1)  +3.0D0/4.0D0*zeta_x(a+1,b,c,1)  -zeta_x(a+3,b,c,1)/8.0D0
!!              zeta_y(a,b,c,1)   = 3.0D0/8.0D0*zeta_y(a-1,b,c,1)  +3.0D0/4.0D0*zeta_y(a+1,b,c,1)  -zeta_y(a+3,b,c,1)/8.0D0
!!              zeta_z(a,b,c,1)   = 3.0D0/8.0D0*zeta_z(a-1,b,c,1)  +3.0D0/4.0D0*zeta_z(a+1,b,c,1)  -zeta_z(a+3,b,c,1)/8.0D0
!!              pi_xx(a,b,c,1)    = 3.0D0/8.0D0*pi_xx(a-1,b,c,1)   +3.0D0/4.0D0*pi_xx(a+1,b,c,1)   -pi_xx(a+3,b,c,1)/8.0D0
!!              pi_yy(a,b,c,1)    = 3.0D0/8.0D0*pi_yy(a-1,b,c,1)   +3.0D0/4.0D0*pi_yy(a+1,b,c,1)   -pi_yy(a+3,b,c,1)/8.0D0
!!              pi_zz(a,b,c,1)    = 3.0D0/8.0D0*pi_zz(a-1,b,c,1)   +3.0D0/4.0D0*pi_zz(a+1,b,c,1)   -pi_zz(a+3,b,c,1)/8.0D0
!!              pi_xy(a,b,c,1)    = 3.0D0/8.0D0*pi_xy(a-1,b,c,1)   +3.0D0/4.0D0*pi_xy(a+1,b,c,1)   -pi_xy(a+3,b,c,1)/8.0D0
!!              pi_xz(a,b,c,1)    = 3.0D0/8.0D0*pi_xz(a-1,b,c,1)   +3.0D0/4.0D0*pi_xz(a+1,b,c,1)   -pi_xz(a+3,b,c,1)/8.0D0
!!              pi_yz(a,b,c,1)    = 3.0D0/8.0D0*pi_yz(a-1,b,c,1)   +3.0D0/4.0D0*pi_yz(a+1,b,c,1)   -pi_yz(a+3,b,c,1)/8.0D0
!!              lambda_x(a,b,c,1) = 3.0D0/8.0D0*lambda_x(a-1,b,c,1)+3.0D0/4.0D0*lambda_x(a+1,b,c,1)-lambda_x(a+3,b,c,1)/8.0D0
!!              lambda_y(a,b,c,1) = 3.0D0/8.0D0*lambda_y(a-1,b,c,1)+3.0D0/4.0D0*lambda_y(a+1,b,c,1)-lambda_y(a+3,b,c,1)/8.0D0
!!              lambda_z(a,b,c,1) = 3.0D0/8.0D0*lambda_z(a-1,b,c,1)+3.0D0/4.0D0*lambda_z(a+1,b,c,1)-lambda_z(a+3,b,c,1)/8.0D0
!            else if (x_down_close .and. x_up_close) then
!              do i=1,dir
!                fi(a,b,c,i) = (fi(a+1,b,c,i)+fi(a-1,b,c,i))/2.0D0
!                gi(a,b,c,i) = (gi(a+1,b,c,i)+gi(a-1,b,c,i))/2.0D0
!              end do
!!              chi(a,b,c,1)    =        (chi(a+1,b,c,1)+     chi(a-1,b,c,1))/2.0D0
!!              zeta_x(a,b,c,1) =     (zeta_x(a+1,b,c,1)+  zeta_x(a-1,b,c,1))/2.0D0
!!              zeta_y(a,b,c,1) =     (zeta_y(a+1,b,c,1)+  zeta_y(a-1,b,c,1))/2.0D0
!!              zeta_z(a,b,c,1) =     (zeta_z(a+1,b,c,1)+  zeta_z(a-1,b,c,1))/2.0D0
!!              pi_xx(a,b,c,1) =       (pi_xx(a+1,b,c,1)+   pi_xx(a-1,b,c,1))/2.0D0
!!              pi_yy(a,b,c,1) =       (pi_yy(a+1,b,c,1)+   pi_yy(a-1,b,c,1))/2.0D0
!!              pi_zz(a,b,c,1) =       (pi_zz(a+1,b,c,1)+   pi_zz(a-1,b,c,1))/2.0D0
!!              pi_xy(a,b,c,1) =       (pi_xy(a+1,b,c,1)+   pi_xy(a-1,b,c,1))/2.0D0
!!              pi_xz(a,b,c,1) =       (pi_xz(a+1,b,c,1)+   pi_xz(a-1,b,c,1))/2.0D0
!!              pi_yz(a,b,c,1) =       (pi_yz(a+1,b,c,1)+   pi_yz(a-1,b,c,1))/2.0D0
!!              lambda_x(a,b,c,1) = (lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))/2.0D0
!!              lambda_y(a,b,c,1) = (lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))/2.0D0
!!              lambda_z(a,b,c,1) = (lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))/2.0D0
!            else
!              needs_interp(a,b,c) = .true.
!              do i =1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!            end if
!            !write(*,*) 'values out on edges',fi(a,b,c,1),gi(a,b,c,1),a,b,c
!          end if
!
!          do i = 1,dir
!            if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!              fi(a,b,c,i) = crs_fi(a,b,c,i)
!              gi(a,b,c,i) = crs_gi(a,b,c,i)
!            end if
!          end do
!
!          end if ! checks if node is a true ghost node
!
!
!        end do
!      end do
!    end do
!    !write(*,*) 'Ruff ruff!',self
!!
!!
!! Face interpolations
!!
!!
!    do c = ducks%lo(3)-nghosts_mid,ducks%hi(3)+nghosts_mid
!      do b = ducks%lo(2)-nghosts_mid,ducks%hi(2)+nghosts_mid
!        do a = ducks%lo(1)-nghosts_mid,ducks%hi(1)+nghosts_mid
!          if (state(a,b,c,1) >=10000000 .and. state(a,b,c,1) <= 1100000000) then !checks if a true ghost node
!            if (mod(a,2) == 0 .and. mod(b,2) /= 0 .and. mod(c,2) /= 0) then
!            x_null = .false.
!            y_null = .false.
!            z_null = .false.
!!
!! Face in yz
!!
!            if (state(a,b,c,7) >= 0 .and. state(a,b,c,39) >= 0) then
!              z_up_clear = .true.
!              z_up_close = .true.
!            else if (state(a,b,c,7) >= 0) then
!              z_up_close = .true.
!              z_up_clear = .false.
!            else
!              z_up_clear = .false.
!              z_up_close = .false.
!            end if
!
!            if (state(a,b,c,4) >= 0 .and. state(a,b,c,36) >= 0) then
!              z_down_clear = .true.
!              z_down_close = .true.
!            else if (state(a,b,c,4) >= 0) then
!              z_down_close = .true.
!              z_down_clear = .false.
!            else
!              z_down_clear = .false.
!              z_down_close = .false.
!            end if
!
!            if (state(a,b,c,6) >= 0 .and. state(a,b,c,38) >= 0) then
!              y_up_clear = .true.
!              y_up_close = .true.
!            else if (state(a,b,c,6) < 0) then
!              y_up_close = .false.
!              y_up_clear = .false.
!            else
!              y_up_close = .true.
!              y_up_clear = .false.
!            end if
!
!            if (state(a,b,c,3) >= 0 .and. state(a,b,c,35) >= 0) then
!              y_down_clear = .true.
!              y_down_close = .true.
!            else if (state(a,b,c,3) < 0) then
!              y_down_close = .false.
!              y_down_clear = .false.
!            else
!              y_down_close = .true.
!              y_down_clear = .false.
!            end if
!
!            if (.not. y_down_close .or. .not. y_up_close) then
!              y_null = .true.
!            end if
!
!            if (.not. z_down_close .or. .not. z_up_close) then
!              z_null = .true.
!            end if
!!
!            if (z_up_clear .and. z_down_clear ) then
!              do i = 1,dir
!                interp_fiz(i) = (9.0D0*(fi(a,b,c+1,i)+fi(a,b,c-1,i))-&
!                                 fi(a,b,c-3,i)-fi(a,b,c+3,i))/16.0D0
!
!                interp_giz(i) = (9.0D0*(gi(a,b,c+1,i)+gi(a,b,c-1,i))-&
!                                 gi(a,b,c+3,i)-gi(a,b,c-3,i))/16.0D0
!              end do
!!              interp_chi_z = (9.0D0*(chi(a,b,c+1,1)+chi(a,b,c-1,1))-&
!!                chi(a,b,c+3,1)-chi(a,b,c-3,1))/16.0D0
!!              interp_zetax_z = (9.0D0*(zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))-&
!!                zeta_x(a,b,c+3,1)-zeta_x(a,b,c-3,1))/16.0D0
!!              interp_zetay_z = (9.0D0*(zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))-&
!!                zeta_y(a,b,c+3,1)-zeta_y(a,b,c-3,1))/16.0D0
!!              interp_zetaz_z = (9.0D0*(zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))-&
!!                zeta_z(a,b,c+3,1)-zeta_z(a,b,c-3,1))/16.0D0
!!              interp_pixx_z = (9.0D0*(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))-&
!!                pi_xx(a,b,c+3,1)-pi_xx(a,b,c-3,1))/16.0D0
!!              interp_piyy_z = (9.0D0*(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))-&
!!                pi_yy(a,b,c+3,1)-pi_yy(a,b,c-3,1))/16.0D0
!!              interp_pizz_z = (9.0D0*(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))-&
!!                pi_zz(a,b,c+3,1)-pi_zz(a,b,c-3,1))/16.0D0
!!              interp_pixy_z = (9.0D0*(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))-&
!!                pi_xy(a,b,c+3,1)-pi_xy(a,b,c-3,1))/16.0D0
!!              interp_pixz_z = (9.0D0*(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))-&
!!                pi_xz(a,b,c+3,1)-pi_xz(a,b,c-3,1))/16.0D0
!!              interp_piyz_z = (9.0D0*(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))-&
!!                pi_yz(a,b,c+3,1)-pi_yz(a,b,c-3,1))/16.0D0
!!              interp_lambdax_z = (9.0D0*(lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))-&
!!                lambda_x(a,b,c+3,1)-lambda_x(a,b,c-3,1))/16.0D0
!!              interp_lambday_z = (9.0D0*(lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))-&
!!                lambda_y(a,b,c+3,1)-lambda_y(a,b,c-3,1))/16.0D0
!!              interp_lambdaz_z = (9.0D0*(lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))-&
!!                lambda_z(a,b,c+3,1)-lambda_z(a,b,c-3,1))/16.0D0
!            else if (z_up_clear .and. z_down_close .and. .not. z_down_clear) then
!              do i = 1,dir
!                interp_fiz(i) = 3.0D0/8.0D0*fi(a,b,c+1,i)+3.0D0/4.0D0*fi(a,b,c-1,i)-fi(a,b,c-3,i)/8.0D0
!                interp_giz(i) = 3.0D0/8.0D0*gi(a,b,c+1,i)+3.0D0/4.0D0*gi(a,b,c-1,i)-gi(a,b,c-3,i)/8.0D0
!              end do
!!              interp_chi_z =     3.0D0/8.0D0*chi(a,b,c+1,1)+3.0D0/4.0D0*chi(a,b,c-1,1)-chi(a,b,c-3,1)/8.0D0
!!              interp_zetax_z =   3.0D0/8.0D0*zeta_x(a,b,c+1,1)+3.0D0/4.0D0*zeta_x(a,b,c-1,1)-zeta_x(a,b,c-3,1)/8.0D0
!!              interp_zetay_z =   3.0D0/8.0D0*zeta_y(a,b,c+1,1)+3.0D0/4.0D0*zeta_y(a,b,c-1,1)-zeta_y(a,b,c-3,1)/8.0D0
!!              interp_zetaz_z =   3.0D0/8.0D0*zeta_z(a,b,c+1,1)+3.0D0/4.0D0*zeta_z(a,b,c-1,1)-zeta_z(a,b,c-3,1)/8.0D0
!!              interp_pixx_z =    3.0D0/8.0D0*pi_xx(a,b,c+1,1)+3.0D0/4.0D0*pi_xx(a,b,c-1,1)-pi_xx(a,b,c-3,1)/8.0D0
!!              interp_piyy_z =    3.0D0/8.0D0*pi_yy(a,b,c+1,1)+3.0D0/4.0D0*pi_yy(a,b,c-1,1)-pi_yy(a,b,c-3,1)/8.0D0
!!              interp_pizz_z =    3.0D0/8.0D0*pi_zz(a,b,c+1,1)+3.0D0/4.0D0*pi_zz(a,b,c-1,1)-pi_zz(a,b,c-3,1)/8.0D0
!!              interp_pixy_z =    3.0D0/8.0D0*pi_xy(a,b,c+1,1)+3.0D0/4.0D0*pi_xy(a,b,c-1,1)-pi_xy(a,b,c-3,1)/8.0D0
!!              interp_pixz_z =    3.0D0/8.0D0*pi_xz(a,b,c+1,1)+3.0D0/4.0D0*pi_xz(a,b,c-1,1)-pi_xz(a,b,c-3,1)/8.0D0
!!              interp_piyz_z =    3.0D0/8.0D0*pi_yz(a,b,c+1,1)+3.0D0/4.0D0*pi_yz(a,b,c-1,1)-pi_yz(a,b,c-3,1)/8.0D0
!!              interp_lambdax_z = 3.0D0/8.0D0*lambda_x(a,b,c+1,1)+3.0D0/4.0D0*lambda_x(a,b,c-1,1)-lambda_x(a,b,c-3,1)/8.0D0
!!              interp_lambday_z = 3.0D0/8.0D0*lambda_y(a,b,c+1,1)+3.0D0/4.0D0*lambda_y(a,b,c-1,1)-lambda_y(a,b,c-3,1)/8.0D0
!!              interp_lambdaz_z = 3.0D0/8.0D0*lambda_z(a,b,c+1,1)+3.0D0/4.0D0*lambda_z(a,b,c-1,1)-lambda_z(a,b,c-3,1)/8.0D0
!            else if (z_down_clear .and. z_up_close .and. .not. z_up_clear) then
!              do i = 1,dir
!                interp_fiz(i) = 3.0D0/8.0D0*fi(a,b,c-1,i)+3.0D0/4.0D0*fi(a,b,c+1,i)-fi(a,b,c+3,i)/8.0D0
!                interp_giz(i) = 3.0D0/8.0D0*gi(a,b,c-1,i)+3.0D0/4.0D0*gi(a,b,c+1,i)-gi(a,b,c+3,i)/8.0D0
!              end do
!!              interp_chi_z =      3.0D0/8.0D0*chi(a,b,c-1,1)+3.0D0/4.0D0*chi(a,b,c+1,1)-chi(a,b,c+3,1)/8.0D0
!!              interp_zetax_z =   3.0D0/8.0D0*zeta_x(a,b,c-1,1)+3.0D0/4.0D0*zeta_x(a,b,c+1,1)-zeta_x(a,b,c+3,1)/8.0D0
!!              interp_zetay_z =   3.0D0/8.0D0*zeta_y(a,b,c-1,1)+3.0D0/4.0D0*zeta_y(a,b,c+1,1)-zeta_y(a,b,c+3,1)/8.0D0
!!              interp_zetaz_z =   3.0D0/8.0D0*zeta_z(a,b,c-1,1)+3.0D0/4.0D0*zeta_z(a,b,c+1,1)-zeta_z(a,b,c+3,1)/8.0D0
!!              interp_pixx_z =    3.0D0/8.0D0*pi_xx(a,b,c-1,1)+3.0D0/4.0D0*pi_xx(a,b,c+1,1)-pi_xx(a,b,c+3,1)/8.0D0
!!              interp_piyy_z =    3.0D0/8.0D0*pi_yy(a,b,c-1,1)+3.0D0/4.0D0*pi_yy(a,b,c+1,1)-pi_yy(a,b,c+3,1)/8.0D0
!!              interp_pizz_z =    3.0D0/8.0D0*pi_zz(a,b,c-1,1)+3.0D0/4.0D0*pi_zz(a,b,c+1,1)-pi_zz(a,b,c+3,1)/8.0D0
!!              interp_pixy_z =    3.0D0/8.0D0*pi_xy(a,b,c-1,1)+3.0D0/4.0D0*pi_xy(a,b,c+1,1)-pi_xy(a,b,c+3,1)/8.0D0
!!              interp_pixz_z =    3.0D0/8.0D0*pi_xz(a,b,c-1,1)+3.0D0/4.0D0*pi_xz(a,b,c+1,1)-pi_xz(a,b,c+3,1)/8.0D0
!!              interp_piyz_z =    3.0D0/8.0D0*pi_yz(a,b,c-1,1)+3.0D0/4.0D0*pi_yz(a,b,c+1,1)-pi_yz(a,b,c+3,1)/8.0D0
!!              interp_lambdax_z= 3.0D0/8.0D0*lambda_x(a,b,c-1,1)+3.0D0/4.0D0*lambda_x(a,b,c+1,1)-lambda_x(a,b,c+3,1)/8.0D0
!!              interp_lambday_z= 3.0D0/8.0D0*lambda_y(a,b,c-1,1)+3.0D0/4.0D0*lambda_y(a,b,c+1,1)-lambda_y(a,b,c+3,1)/8.0D0
!!              interp_lambdaz_z= 3.0D0/8.0D0*lambda_z(a,b,c-1,1)+3.0D0/4.0D0*lambda_z(a,b,c+1,1)-lambda_z(a,b,c+3,1)/8.0D0
!            else if (z_up_close .and. z_down_close) then
!              do i=1,dir
!                interp_fiz(i) = (fi(a,b,c+1,i)+fi(a,b,c-1,i))/2.0D0
!                interp_giz(i) = (gi(a,b,c+1,i)+gi(a,b,c-1,i))/2.0D0
!              end do
!!              interp_chi_z = (chi(a,b,c+1,1)+chi(a,b,c-1,1))/2.0D0
!!              interp_zetax_z = (zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))/2.0D0
!!              interp_zetay_z = (zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))/2.0D0
!!              interp_zetaz_z = (zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))/2.0D0
!!              interp_pixx_z =(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))/2.0D0
!!              interp_piyy_z =(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))/2.0D0
!!              interp_pizz_z =(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))/2.0D0
!!              interp_pixy_z =(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))/2.0D0
!!              interp_pixz_z =(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))/2.0D0
!!              interp_piyz_z =(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))/2.0D0
!!              interp_lambdax_z = (lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))/2.0D0
!!              interp_lambday_z = (lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))/2.0D0
!!              interp_lambdaz_z = (lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))/2.0D0
!            else
!              do i =1,dir
!                interp_fiz(i) = crs_fi(a,b,c,i)
!                interp_giz(i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!            if (y_up_clear .and. y_down_clear) then
!              !if (a == 58 .and. b == 67 .and. c == 5) write(*,*) 'giggle me elmo!'
!              do i = 1,dir
!                interp_fiy(i) = (9.0D0*(fi(a,b+1,c,i)+fi(a,b-1,c,i))-&
!                                 fi(a,b-3,c,i)-fi(a,b+3,c,i))/16.0D0
!                interp_giy(i) = (9.0D0*(gi(a,b+1,c,i)+gi(a,b-1,c,i))-&
!                                 gi(a,b+3,c,i)-gi(a,b-3,c,i))/16.0D0
!              end do
!!              interp_chi_y = (9.0D0*(chi(a,b+1,c,1)+chi(a,b-1,c,1))-&
!!                chi(a,b+3,c,1)-chi(a,b-3,c,1))/16.0D0
!!              interp_zetax_y = (9.0D0*(zeta_x(a,b+1,c,1)+zeta_x(a,b-1,c,1))-&
!!                zeta_x(a,b+3,c,1)-zeta_x(a,b-3,c,1))/16.0D0
!!              interp_zetay_y = (9.0D0*(zeta_y(a,b+1,c,1)+zeta_y(a,b-1,c,1))-&
!!                zeta_y(a,b+3,c,1)-zeta_y(a,b-3,c,1))/16.0D0
!!              interp_zetaz_y = (9.0D0*(zeta_z(a,b+1,c,1)+zeta_z(a,b-1,c,1))-&
!!                zeta_z(a,b+3,c,1)-zeta_z(a,b-3,c,1))/16.0D0
!!              interp_pixx_y = (9.0D0*(pi_xx(a,b+1,c,1)+pi_xx(a,b-1,c,1))-&
!!                pi_xx(a,b+3,c,1)-pi_xx(a,b-3,c,1))/16.0D0
!!              interp_piyy_y = (9.0D0*(pi_yy(a,b+1,c,1)+pi_yy(a,b-1,c,1))-&
!!                pi_yy(a,b-3,c,1)-pi_yy(a,b+3,c,1))/16.0D0
!!              interp_pizz_y = (9.0D0*(pi_zz(a,b+1,c,1)+pi_zz(a,b-1,c,1))-&
!!                pi_zz(a,b+3,c,1)-pi_zz(a,b-3,c,1))/16.0D0
!!              interp_pixy_y = (9.0D0*(pi_xy(a,b+1,c,1)+pi_xy(a,b-1,c,1))-&
!!                pi_xy(a,b+3,c,1)-pi_xy(a,b-3,c,1))/16.0D0
!!              interp_pixz_y = (9.0D0*(pi_xz(a,b+1,c,1)+pi_xz(a,b-1,c,1))-&
!!                pi_xz(a,b+3,c,1)-pi_xz(a,b-3,c,1))/16.0D0
!!              interp_piyz_y = (9.0D0*(pi_yz(a,b+1,c,1)+pi_yz(a,b-1,c,1))-&
!!                pi_yz(a,b+3,c,1)-pi_yz(a,b-3,c,1))/16.0D0
!!              interp_lambdax_y = (9.0D0*(lambda_x(a,b+1,c,1)+lambda_x(a,b-1,c,1))-&
!!                lambda_x(a,b+3,c,1)-lambda_x(a,b-3,c,1))/16.0D0
!!              interp_lambday_y = (9.0D0*(lambda_y(a,b+1,c,1)+lambda_y(a,b-1,c,1))-&
!!                lambda_y(a,b+3,c,1)-lambda_y(a,b-3,c,1))/16.0D0
!!              interp_lambdaz_y = (9.0D0*(lambda_z(a,b+1,c,1)+lambda_z(a,b-1,c,1))-&
!!                lambda_z(a,b+3,c,1)-lambda_z(a,b-3,c,1))/16.0D0
!            else if (y_up_clear .and. y_down_close .and. .not. y_down_clear) then
!              do i = 1,dir
!                interp_fiy(i) = 3.0D0/8.0D0*fi(a,b+1,c,i)+3.0D0/4.0D0*fi(a,b-1,c,i)-fi(a,b-3,c,i)/8.0D0
!                interp_giy(i) = 3.0D0/8.0D0*gi(a,b+1,c,i)+3.0D0/4.0D0*gi(a,b-1,c,i)-gi(a,b-3,c,i)/8.0D0
!              end do
!!              interp_chi_y = 3.0D0/8.0D0*chi(a,b+1,c,1)     +3.0D0/4.0D0*chi(a,b-1,c,1)     -chi(a,b-3,c,1)/8.0D0
!!              interp_zetax_y = 3.0D0/8.0D0*zeta_x(a,b+1,c,1)  +3.0D0/4.0D0*zeta_x(a,b-1,c,1)  -zeta_x(a,b-3,c,1)/8.0D0
!!              interp_zetay_y = 3.0D0/8.0D0*zeta_y(a,b+1,c,1)  +3.0D0/4.0D0*zeta_y(a,b-1,c,1)  -zeta_y(a,b-3,c,1)/8.0D0
!!              interp_zetaz_y = 3.0D0/8.0D0*zeta_z(a,b+1,c,1)  +3.0D0/4.0D0*zeta_z(a,b-1,c,1)  -zeta_z(a,b-3,c,1)/8.0D0
!!              interp_pixx_y = 3.0D0/8.0D0*pi_xx(a,b+1,c,1)   +3.0D0/4.0D0*pi_xx(a,b-1,c,1)   -pi_xx(a,b-3,c,1)/8.0D0
!!              interp_piyy_y = 3.0D0/8.0D0*pi_yy(a,b+1,c,1)   +3.0D0/4.0D0*pi_yy(a,b-1,c,1)   -pi_yy(a,b-3,c,1)/8.0D0
!!              interp_pizz_y = 3.0D0/8.0D0*pi_zz(a,b+1,c,1)   +3.0D0/4.0D0*pi_zz(a,b-1,c,1)   -pi_zz(a,b-3,c,1)/8.0D0
!!              interp_pixy_y = 3.0D0/8.0D0*pi_xy(a,b+1,c,1)   +3.0D0/4.0D0*pi_xy(a,b-1,c,1)   -pi_xy(a,b-3,c,1)/8.0D0
!!              interp_pixz_y = 3.0D0/8.0D0*pi_xz(a,b+1,c,1)   +3.0D0/4.0D0*pi_xz(a,b-1,c,1)   -pi_xz(a,b-3,c,1)/8.0D0
!!              interp_piyz_y = 3.0D0/8.0D0*pi_yz(a,b+1,c,1)   +3.0D0/4.0D0*pi_yz(a,b-1,c,1)   -pi_yz(a,b-3,c,1)/8.0D0
!!              interp_lambdax_y = 3.0D0/8.0D0*lambda_x(a,b+1,c,1)+3.0D0/4.0D0*lambda_x(a,b-1,c,1)-lambda_x(a,b-3,c,1)/8.0D0
!!              interp_lambday_y = 3.0D0/8.0D0*lambda_y(a,b+1,c,1)+3.0D0/4.0D0*lambda_y(a,b-1,c,1)-lambda_y(a,b-3,c,1)/8.0D0
!!              interp_lambdaz_y = 3.0D0/8.0D0*lambda_z(a,b+1,c,1)+3.0D0/4.0D0*lambda_z(a,b-1,c,1)-lambda_z(a,b-3,c,1)/8.0D0
!            else if (y_down_clear .and. y_up_close .and. .not. y_up_clear) then
!              do i = 1,dir
!                interp_fiy(i) = 3.0D0/8.0D0*fi(a,b-1,c,i)+3.0D0/4.0D0*fi(a,b+1,c,i)-fi(a,b+3,c,i)/8.0D0
!                interp_giy(i) = 3.0D0/8.0D0*gi(a,b-1,c,i)+3.0D0/4.0D0*gi(a,b+1,c,i)-gi(a,b+3,c,i)/8.0D0
!              end do
!!              interp_chi_y = 3.0D0/8.0D0*chi(a,b-1,c,1)     +3.0D0/4.0D0*chi(a,b+1,c,1)     -chi(a,b+3,c,1)/8.0D0
!!              interp_zetax_y = 3.0D0/8.0D0*zeta_x(a,b-1,c,1)  +3.0D0/4.0D0*zeta_x(a,b+1,c,1)  -zeta_x(a,b+3,c,1)/8.0D0
!!              interp_zetay_y = 3.0D0/8.0D0*zeta_y(a,b-1,c,1)  +3.0D0/4.0D0*zeta_y(a,b+1,c,1)  -zeta_y(a,b+3,c,1)/8.0D0
!!              interp_zetaz_y = 3.0D0/8.0D0*zeta_z(a,b-1,c,1)  +3.0D0/4.0D0*zeta_z(a,b+1,c,1)  -zeta_z(a,b+3,c,1)/8.0D0
!!              interp_pixx_y = 3.0D0/8.0D0*pi_xx(a,b-1,c,1)   +3.0D0/4.0D0*pi_xx(a,b+1,c,1)   -pi_xx(a,b+3,c,1)/8.0D0
!!              interp_piyy_y = 3.0D0/8.0D0*pi_yy(a,b-1,c,1)   +3.0D0/4.0D0*pi_yy(a,b+1,c,1)   -pi_yy(a,b+3,c,1)/8.0D0
!!              interp_pizz_y = 3.0D0/8.0D0*pi_zz(a,b-1,c,1)   +3.0D0/4.0D0*pi_zz(a,b+1,c,1)   -pi_zz(a,b+3,c,1)/8.0D0
!!              interp_pixy_y = 3.0D0/8.0D0*pi_xy(a,b-1,c,1)   +3.0D0/4.0D0*pi_xy(a,b+1,c,1)   -pi_xy(a,b+3,c,1)/8.0D0
!!              interp_pixz_y = 3.0D0/8.0D0*pi_xz(a,b-1,c,1)   +3.0D0/4.0D0*pi_xz(a,b+1,c,1)   -pi_xz(a,b+3,c,1)/8.0D0
!!              interp_piyz_y = 3.0D0/8.0D0*pi_yz(a,b-1,c,1)   +3.0D0/4.0D0*pi_yz(a,b+1,c,1)   -pi_yz(a,b+3,c,1)/8.0D0
!!              interp_lambdax_y = 3.0D0/8.0D0*lambda_x(a,b-1,c,1)+3.0D0/4.0D0*lambda_x(a,b+1,c,1)-lambda_x(a,b+3,c,1)/8.0D0
!!              interp_lambday_y = 3.0D0/8.0D0*lambda_y(a,b-1,c,1)+3.0D0/4.0D0*lambda_y(a,b+1,c,1)-lambda_y(a,b+3,c,1)/8.0D0
!!              interp_lambdaz_y = 3.0D0/8.0D0*lambda_z(a,b-1,c,1)+3.0D0/4.0D0*lambda_z(a,b+1,c,1)-lambda_z(a,b+3,c,1)/8.0D0
!            else if (y_up_close .and. y_down_close) then
!              do i=1,dir
!                interp_fiy(i) = (fi(a,b+1,c,i)+fi(a,b-1,c,i))/2.0D0
!                interp_giy(i) = (gi(a,b+1,c,i)+gi(a,b-1,c,i))/2.0D0
!              end do
!!              interp_chi_y    =       (chi(a,b+1,c,1)+     chi(a,b-1,c,1))/2.0D0
!!              interp_zetax_y =     (zeta_x(a,b+1,c,1)+  zeta_x(a,b-1,c,1))/2.0D0
!!              interp_zetay_y =     (zeta_y(a,b+1,c,1)+  zeta_y(a,b-1,c,1))/2.0D0
!!              interp_zetaz_y =     (zeta_z(a,b+1,c,1)+  zeta_z(a,b-1,c,1))/2.0D0
!!              interp_pixx_y =       (pi_xx(a,b+1,c,1)+   pi_xx(a,b-1,c,1))/2.0D0
!!              interp_piyy_y =       (pi_yy(a,b+1,c,1)+   pi_yy(a,b-1,c,1))/2.0D0
!!              interp_pizz_y =       (pi_zz(a,b+1,c,1)+   pi_zz(a,b-1,c,1))/2.0D0
!!              interp_pixy_y =       (pi_xy(a,b+1,c,1)+   pi_xy(a,b-1,c,1))/2.0D0
!!              interp_pixz_y =       (pi_xz(a,b+1,c,1)+   pi_xz(a,b-1,c,1))/2.0D0
!!              interp_piyz_y =       (pi_yz(a,b+1,c,1)+   pi_yz(a,b-1,c,1))/2.0D0
!!              interp_lambdax_y = (lambda_x(a,b+1,c,1)+lambda_x(a,b-1,c,1))/2.0D0
!!              interp_lambday_y = (lambda_y(a,b+1,c,1)+lambda_y(a,b-1,c,1))/2.0D0
!!              interp_lambdaz_y = (lambda_z(a,b+1,c,1)+lambda_z(a,b-1,c,1))/2.0D0
!            else
!              !needs_interp(a,b,c) = .true.
!              do i =1,dir
!                interp_fiy(i) = crs_fi(a,b,c,i)
!                interp_giy(i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!            if (.not. y_null .and. .not. z_null) then
!
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fiz(i)+interp_fiy(i))/2.0D0
!                gi(a,b,c,i) = (interp_giz(i)+interp_giy(i))/2.0D0
!              end do
!!              chi(a,b,c,1)    = (interp_chi_z + interp_chi_y)/2.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_z + interp_zetax_y)/2.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_z + interp_zetay_y)/2.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_z + interp_zetaz_y)/2.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_z + interp_pixx_y)/2.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_z + interp_piyy_y)/2.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_z + interp_pizz_y)/2.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_z + interp_pixy_y)/2.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_z + interp_pixz_y)/2.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_z + interp_piyz_y)/2.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_z + interp_lambdax_y)/2.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_z + interp_lambday_y)/2.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_z + interp_lambdaz_y)/2.0D0
!            else if (z_null .and. .not. y_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fiy(i)
!                gi(a,b,c,i) = interp_giy(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_y
!!              zeta_x(a,b,c,1) =  interp_zetax_y
!!              zeta_y(a,b,c,1) =  interp_zetay_y
!!              zeta_z(a,b,c,1) =  interp_zetaz_y
!!              pi_xx(a,b,c,1) =  interp_pixx_y
!!              pi_yy(a,b,c,1) =  interp_piyy_y
!!              pi_zz(a,b,c,1) =  interp_pizz_y
!!              pi_xy(a,b,c,1) =  interp_pixy_y
!!              pi_xz(a,b,c,1) =  interp_pixz_y
!!              pi_yz(a,b,c,1) =  interp_piyz_y
!!              lambda_x(a,b,c,1) =interp_lambdax_y
!!              lambda_y(a,b,c,1) =interp_lambday_y
!!              lambda_z(a,b,c,1) =interp_lambdaz_y
!            else if (y_null .and. .not. z_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fiz(i)
!                gi(a,b,c,i) = interp_giz(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_z
!!              zeta_x(a,b,c,1) =  interp_zetax_z
!!              zeta_y(a,b,c,1) =  interp_zetay_z
!!              zeta_z(a,b,c,1) =  interp_zetaz_z
!!              pi_xx(a,b,c,1) =  interp_pixx_z
!!              pi_yy(a,b,c,1) =  interp_piyy_z
!!              pi_zz(a,b,c,1) =  interp_pizz_z
!!              pi_xy(a,b,c,1) =  interp_pixy_z
!!              pi_xz(a,b,c,1) =  interp_pixz_z
!!              pi_yz(a,b,c,1) =  interp_piyz_z
!!              lambda_x(a,b,c,1) =interp_lambdax_z
!!              lambda_y(a,b,c,1) =interp_lambday_z
!!              lambda_z(a,b,c,1) =interp_lambdaz_z
!            else
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!!            if (a == 58 .and. b == 67 .and. c == 5) then
!!              write(*,*) 'fi interlevel inside',fi(a,b,c,1:dir)
!!              write(*,*) 'gi interlevel inside',gi(a,b,c,1:dir)
!!              write(*,*) 'coarse fi',crs_fi(a,b,c,1:dir)
!!              write(*,*) 'coarse gi',crs_gi(a,b,c,1:dir)
!!              write(*,*) 'interp fiy',interp_fiy
!!              write(*,*) 'interp giy',interp_giy
!!              write(*,*) 'interp fiz',interp_fiz
!!              write(*,*) 'interp giz',interp_giz
!!              write(*,*) 'y bools',y_up_clear,y_up_close,y_down_clear,y_down_close,&
!!                y_null
!!              write(*,*) 'z_bools',z_up_clear,z_up_close,z_down_clear,z_down_close,&
!!                z_null
!!              write(*,*) 'z fi values for 37',fi(a,b,c-1,37),fi(a,b,c-3,37),fi(a,b,c+1,37),&
!!                fi(a,b,c+3,37)
!!              write(*,*) 'y fi values for 37',fi(a,b-1,c,37),fi(a,b-3,c,37),fi(a,b+1,c,37),&
!!                fi(a,b+3,c,37)
!!            end if
!
!!
!! Face in xy
!!
!          else if (mod(a,2) /= 0 .and. mod(c,2) == 0 .and. mod(b,2) /= 0) then
!            x_null = .false.
!            y_null = .false.
!            z_null = .false.
!
!            if (state(a,b,c,6) >= 0 .and. state(a,b,c,38) >= 0) then
!              y_up_clear = .true.
!              y_up_close = .true.
!            else if (state(a,b,c,6) < 0) then
!              y_up_close = .false.
!              y_up_clear = .false.
!            else
!              y_up_close = .true.
!              y_up_clear = .false.
!            end if
!            if (state(a,b,c,3) >= 0 .and. state(a,b,c,35) >= 0) then
!              y_down_clear = .true.
!              y_down_close = .true.
!            else if (state(a,b,c,3) < 0) then
!              y_down_close = .false.
!              y_down_clear = .false.
!            else
!              y_down_close = .true.
!              y_down_clear = .false.
!            end if
!
!            if (state(a,b,c,5) >= 0 .and. state(a,b,c,37) >= 0) then
!              x_up_clear = .true.
!              x_up_close = .true.
!            else if (state(a,b,c,5) < 0) then
!              x_up_close = .false.
!              x_up_clear = .false.
!            else
!              x_up_close = .true.
!              x_up_clear = .false.
!            end if
!            if (state(a,b,c,2) >= 0 .and. state(a,b,c,34) >= 0) then
!              x_down_clear = .true.
!              x_down_close = .true.
!            else if (state(a,b,c,2) < 0) then
!              x_down_clear = .false.
!              x_down_close = .false.
!            else
!              x_down_close = .true.
!              x_down_clear = .false.
!            end if
!!
!!
!!
!            if (.not. y_down_close .or. .not. y_up_close) then
!              y_null = .true.
!            end if
!            if (.not. x_down_close .or. .not. x_up_close) then
!              x_null = .true.
!            end if
!!
!!
!!
!            if (y_up_clear .and. y_down_clear) then
!              do i = 1,dir
!                interp_fiy(i) = (9.0D0*(fi(a,b+1,c,i)+fi(a,b-1,c,i))-&
!                                 fi(a,b-3,c,i)-fi(a,b+3,c,i))/16.0D0
!                interp_giy(i) = (9.0D0*(gi(a,b+1,c,i)+gi(a,b-1,c,i))-&
!                                 gi(a,b+3,c,i)-gi(a,b-3,c,i))/16.0D0
!              end do
!!              interp_chi_y = (9.0D0*(chi(a,b+1,c,1)+chi(a,b-1,c,1))-&
!!                chi(a,b+3,c,1)-chi(a,b-3,c,1))/16.0D0
!!              interp_zetax_y = (9.0D0*(zeta_x(a,b+1,c,1)+zeta_x(a,b-1,c,1))-&
!!                zeta_x(a,b+3,c,1)-zeta_x(a,b-3,c,1))/16.0D0
!!              interp_zetay_y = (9.0D0*(zeta_y(a,b+1,c,1)+zeta_y(a,b-1,c,1))-&
!!                zeta_y(a,b+3,c,1)-zeta_y(a,b-3,c,1))/16.0D0
!!              interp_zetaz_y = (9.0D0*(zeta_z(a,b+1,c,1)+zeta_z(a,b-1,c,1))-&
!!                zeta_z(a,b+3,c,1)-zeta_z(a,b-3,c,1))/16.0D0
!!              interp_pixx_y = (9.0D0*(pi_xx(a,b+1,c,1)+pi_xx(a,b-1,c,1))-&
!!                pi_xx(a,b+3,c,1)-pi_xx(a,b-3,c,1))/16.0D0
!!              interp_piyy_y = (9.0D0*(pi_yy(a,b+1,c,1)+pi_yy(a,b-1,c,1))-&
!!                pi_yy(a,b-3,c,1)-pi_yy(a,b+3,c,1))/16.0D0
!!              interp_pizz_y = (9.0D0*(pi_zz(a,b+1,c,1)+pi_zz(a,b-1,c,1))-&
!!                pi_zz(a,b+3,c,1)-pi_zz(a,b-3,c,1))/16.0D0
!!              interp_pixy_y = (9.0D0*(pi_xy(a,b+1,c,1)+pi_xy(a,b-1,c,1))-&
!!                pi_xy(a,b+3,c,1)-pi_xy(a,b-3,c,1))/16.0D0
!!              interp_pixz_y = (9.0D0*(pi_xz(a,b+1,c,1)+pi_xz(a,b-1,c,1))-&
!!                pi_xz(a,b+3,c,1)-pi_xz(a,b-3,c,1))/16.0D0
!!              interp_piyz_y = (9.0D0*(pi_yz(a,b+1,c,1)+pi_yz(a,b-1,c,1))-&
!!                pi_yz(a,b+3,c,1)-pi_yz(a,b-3,c,1))/16.0D0
!!              interp_lambdax_y = (9.0D0*(lambda_x(a,b+1,c,1)+lambda_x(a,b-1,c,1))-&
!!                lambda_x(a,b+3,c,1)-lambda_x(a,b-3,c,1))/16.0D0
!!              interp_lambday_y = (9.0D0*(lambda_y(a,b+1,c,1)+lambda_y(a,b-1,c,1))-&
!!                lambda_y(a,b+3,c,1)-lambda_y(a,b-3,c,1))/16.0D0
!!              interp_lambdaz_y = (9.0D0*(lambda_z(a,b+1,c,1)+lambda_z(a,b-1,c,1))-&
!!                lambda_z(a,b+3,c,1)-lambda_z(a,b-3,c,1))/16.0D0
!            else if (y_up_clear .and. y_down_close .and. .not. y_down_clear) then
!              do i = 1,dir
!                interp_fiy(i) = 3.0D0/8.0D0*fi(a,b+1,c,i)+3.0D0/4.0D0*fi(a,b-1,c,i)-fi(a,b-3,c,i)/8.0D0
!                interp_giy(i) = 3.0D0/8.0D0*gi(a,b+1,c,i)+3.0D0/4.0D0*gi(a,b-1,c,i)-gi(a,b-3,c,i)/8.0D0
!              end do
!!              interp_chi_y = 3.0D0/8.0D0*chi(a,b+1,c,1)     +3.0D0/4.0D0*chi(a,b-1,c,1)     -chi(a,b-3,c,1)/8.0D0
!!              interp_zetax_y = 3.0D0/8.0D0*zeta_x(a,b+1,c,1)  +3.0D0/4.0D0*zeta_x(a,b-1,c,1)  -zeta_x(a,b-3,c,1)/8.0D0
!!              interp_zetay_y = 3.0D0/8.0D0*zeta_y(a,b+1,c,1)  +3.0D0/4.0D0*zeta_y(a,b-1,c,1)  -zeta_y(a,b-3,c,1)/8.0D0
!!              interp_zetaz_y = 3.0D0/8.0D0*zeta_z(a,b+1,c,1)  +3.0D0/4.0D0*zeta_z(a,b-1,c,1)  -zeta_z(a,b-3,c,1)/8.0D0
!!              interp_pixx_y = 3.0D0/8.0D0*pi_xx(a,b+1,c,1)   +3.0D0/4.0D0*pi_xx(a,b-1,c,1)   -pi_xx(a,b-3,c,1)/8.0D0
!!              interp_piyy_y = 3.0D0/8.0D0*pi_yy(a,b+1,c,1)   +3.0D0/4.0D0*pi_yy(a,b-1,c,1)   -pi_yy(a,b-3,c,1)/8.0D0
!!              interp_pizz_y = 3.0D0/8.0D0*pi_zz(a,b+1,c,1)   +3.0D0/4.0D0*pi_zz(a,b-1,c,1)   -pi_zz(a,b-3,c,1)/8.0D0
!!              interp_pixy_y = 3.0D0/8.0D0*pi_xy(a,b+1,c,1)   +3.0D0/4.0D0*pi_xy(a,b-1,c,1)   -pi_xy(a,b-3,c,1)/8.0D0
!!              interp_pixz_y = 3.0D0/8.0D0*pi_xz(a,b+1,c,1)   +3.0D0/4.0D0*pi_xz(a,b-1,c,1)   -pi_xz(a,b-3,c,1)/8.0D0
!!              interp_piyz_y = 3.0D0/8.0D0*pi_yz(a,b+1,c,1)   +3.0D0/4.0D0*pi_yz(a,b-1,c,1)   -pi_yz(a,b-3,c,1)/8.0D0
!!              interp_lambdax_y = 3.0D0/8.0D0*lambda_x(a,b+1,c,1)+3.0D0/4.0D0*lambda_x(a,b-1,c,1)-lambda_x(a,b-3,c,1)/8.0D0
!!              interp_lambday_y = 3.0D0/8.0D0*lambda_y(a,b+1,c,1)+3.0D0/4.0D0*lambda_y(a,b-1,c,1)-lambda_y(a,b-3,c,1)/8.0D0
!!              interp_lambdaz_y = 3.0D0/8.0D0*lambda_z(a,b+1,c,1)+3.0D0/4.0D0*lambda_z(a,b-1,c,1)-lambda_z(a,b-3,c,1)/8.0D0
!            else if (y_down_clear .and. y_up_close .and. .not. y_up_clear) then
!              do i = 1,dir
!                interp_fiy(i) = 3.0D0/8.0D0*fi(a,b-1,c,i)+3.0D0/4.0D0*fi(a,b+1,c,i)-fi(a,b+3,c,i)/8.0D0
!                interp_giy(i) = 3.0D0/8.0D0*gi(a,b-1,c,i)+3.0D0/4.0D0*gi(a,b+1,c,i)-gi(a,b+3,c,i)/8.0D0
!              end do
!!              interp_chi_y = 3.0D0/8.0D0*chi(a,b-1,c,1)     +3.0D0/4.0D0*chi(a,b+1,c,1)     -chi(a,b+3,c,1)/8.0D0
!!              interp_zetax_y = 3.0D0/8.0D0*zeta_x(a,b-1,c,1)  +3.0D0/4.0D0*zeta_x(a,b+1,c,1)  -zeta_x(a,b+3,c,1)/8.0D0
!!              interp_zetay_y = 3.0D0/8.0D0*zeta_y(a,b-1,c,1)  +3.0D0/4.0D0*zeta_y(a,b+1,c,1)  -zeta_y(a,b+3,c,1)/8.0D0
!!              interp_zetaz_y = 3.0D0/8.0D0*zeta_z(a,b-1,c,1)  +3.0D0/4.0D0*zeta_z(a,b+1,c,1)  -zeta_z(a,b+3,c,1)/8.0D0
!!              interp_pixx_y = 3.0D0/8.0D0*pi_xx(a,b-1,c,1)   +3.0D0/4.0D0*pi_xx(a,b+1,c,1)   -pi_xx(a,b+3,c,1)/8.0D0
!!              interp_piyy_y = 3.0D0/8.0D0*pi_yy(a,b-1,c,1)   +3.0D0/4.0D0*pi_yy(a,b+1,c,1)   -pi_yy(a,b+3,c,1)/8.0D0
!!              interp_pizz_y = 3.0D0/8.0D0*pi_zz(a,b-1,c,1)   +3.0D0/4.0D0*pi_zz(a,b+1,c,1)   -pi_zz(a,b+3,c,1)/8.0D0
!!              interp_pixy_y = 3.0D0/8.0D0*pi_xy(a,b-1,c,1)   +3.0D0/4.0D0*pi_xy(a,b+1,c,1)   -pi_xy(a,b+3,c,1)/8.0D0
!!              interp_pixz_y = 3.0D0/8.0D0*pi_xz(a,b-1,c,1)   +3.0D0/4.0D0*pi_xz(a,b+1,c,1)   -pi_xz(a,b+3,c,1)/8.0D0
!!              interp_piyz_y = 3.0D0/8.0D0*pi_yz(a,b-1,c,1)   +3.0D0/4.0D0*pi_yz(a,b+1,c,1)   -pi_yz(a,b+3,c,1)/8.0D0
!!              interp_lambdax_y = 3.0D0/8.0D0*lambda_x(a,b-1,c,1)+3.0D0/4.0D0*lambda_x(a,b+1,c,1)-lambda_x(a,b+3,c,1)/8.0D0
!!              interp_lambday_y = 3.0D0/8.0D0*lambda_y(a,b-1,c,1)+3.0D0/4.0D0*lambda_y(a,b+1,c,1)-lambda_y(a,b+3,c,1)/8.0D0
!!              interp_lambdaz_y = 3.0D0/8.0D0*lambda_z(a,b-1,c,1)+3.0D0/4.0D0*lambda_z(a,b+1,c,1)-lambda_z(a,b+3,c,1)/8.0D0
!            else if (y_down_close .and. y_up_close) then
!              do i=1,dir
!                interp_fiy(i) = (fi(a,b+1,c,i)+fi(a,b-1,c,i))/2.0D0
!                interp_giy(i) = (gi(a,b+1,c,i)+gi(a,b-1,c,i))/2.0D0
!              end do
!!              interp_chi_y    =       (chi(a,b+1,c,1)+     chi(a,b-1,c,1))/2.0D0
!!              interp_zetax_y =     (zeta_x(a,b+1,c,1)+  zeta_x(a,b-1,c,1))/2.0D0
!!              interp_zetay_y =     (zeta_y(a,b+1,c,1)+  zeta_y(a,b-1,c,1))/2.0D0
!!              interp_zetaz_y =     (zeta_z(a,b+1,c,1)+  zeta_z(a,b-1,c,1))/2.0D0
!!              interp_pixx_y =       (pi_xx(a,b+1,c,1)+   pi_xx(a,b-1,c,1))/2.0D0
!!              interp_piyy_y =       (pi_yy(a,b+1,c,1)+   pi_yy(a,b-1,c,1))/2.0D0
!!              interp_pizz_y =       (pi_zz(a,b+1,c,1)+   pi_zz(a,b-1,c,1))/2.0D0
!!              interp_pixy_y =       (pi_xy(a,b+1,c,1)+   pi_xy(a,b-1,c,1))/2.0D0
!!              interp_pixz_y =       (pi_xz(a,b+1,c,1)+   pi_xz(a,b-1,c,1))/2.0D0
!!              interp_piyz_y =       (pi_yz(a,b+1,c,1)+   pi_yz(a,b-1,c,1))/2.0D0
!!              interp_lambdax_y = (lambda_x(a,b+1,c,1)+lambda_x(a,b-1,c,1))/2.0D0
!!              interp_lambday_y = (lambda_y(a,b+1,c,1)+lambda_y(a,b-1,c,1))/2.0D0
!!              interp_lambdaz_y = (lambda_z(a,b+1,c,1)+lambda_z(a,b-1,c,1))/2.0D0
!            else
!              !needs_interp(a,b,c) = .true.
!              do i =1,dir
!                interp_fiy(i) = crs_fi(a,b,c,i)
!                interp_giy(i) = crs_gi(a,b,c,i)
!              end do
!
!            end if
!
!            if (x_up_clear .and. x_down_clear ) then
!              do i = 1,dir
!                interp_fix(i) = (9.0D0*(fi(a+1,b,c,i)+fi(a-1,b,c,i))-&
!                                 fi(a-3,b,c,i)-fi(a+3,b,c,i))/16.0D0
!
!                interp_gix(i) = (9.0D0*(gi(a+1,b,c,i)+gi(a-1,b,c,i))-&
!                                 gi(a+3,b,c,i)-gi(a-3,b,c,i))/16.0D0
!              end do
!!              interp_chi_x = (9.0D0*(chi(a+1,b,c,1)+chi(a-1,b,c,1))-&
!!                chi(a+3,b,c,1)-chi(a-3,b,c,1))/16.0D0
!!              interp_zetax_x = (9.0D0*(zeta_x(a+1,b,c,1)+zeta_x(a-1,b,c,1))-&
!!                zeta_x(a+3,b,c,1)-zeta_x(a-3,b,c,1))/16.0D0
!!              interp_zetay_x = (9.0D0*(zeta_y(a+1,b,c,1)+zeta_y(a-1,b,c,1))-&
!!                zeta_y(a+3,b,c,1)-zeta_y(a-3,b,c,1))/16.0D0
!!              interp_zetaz_x = (9.0D0*(zeta_z(a+1,b,c,1)+zeta_z(a-1,b,c,1))-&
!!                zeta_z(a+3,b,c,1)-zeta_z(a-3,b,c,1))/16.0D0
!!              interp_pixx_x = (9.0D0*(pi_xx(a+1,b,c,1)+pi_xx(a-1,b,c,1))-&
!!                pi_xx(a+3,b,c,1)-pi_xx(a-3,b,c,1))/16.0D0
!!              interp_piyy_x = (9.0D0*(pi_yy(a+1,b,c,1)+pi_yy(a-1,b,c,1))-&
!!                pi_yy(a-3,b,c,1)-pi_yy(a+3,b,c,1))/16.0D0
!!              interp_pizz_x = (9.0D0*(pi_zz(a+1,b,c,1)+pi_zz(a-1,b,c,1))-&
!!                pi_zz(a+3,b,c,1)-pi_zz(a-3,b,c,1))/16.0D0
!!              interp_pixy_x = (9.0D0*(pi_xy(a+1,b,c,1)+pi_xy(a-1,b,c,1))-&
!!                pi_xy(a+3,b,c,1)-pi_xy(a-3,b,c,1))/16.0D0
!!              interp_pixz_x = (9.0D0*(pi_xz(a+1,b,c,1)+pi_xz(a-1,b,c,1))-&
!!                pi_xz(a+3,b,c,1)-pi_xz(a-3,b,c,1))/16.0D0
!!              interp_piyz_x = (9.0D0*(pi_yz(a+1,b,c,1)+pi_yz(a-1,b,c,1))-&
!!                pi_yz(a+3,b,c,1)-pi_yz(a-3,b,c,1))/16.0D0
!!              interp_lambdax_x = (9.0D0*(lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))-&
!!                lambda_x(a+3,b,c,1)-lambda_x(a-3,b,c,1))/16.0D0
!!              interp_lambday_x = (9.0D0*(lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))-&
!!                lambda_y(a+3,b,c,1)-lambda_y(a-3,b,c,1))/16.0D0
!!              interp_lambdaz_x = (9.0D0*(lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))-&
!!                lambda_z(a+3,b,c,1)-lambda_z(a-3,b,c,1))/16.0D0
!            else if (x_up_clear .and. x_down_close .and. .not. x_down_clear) then
!              do i = 1,dir
!                interp_fix(i) = 3.0D0/8.0D0*fi(a+1,b,c,i)+3.0D0/4.0D0*fi(a-1,b,c,i)-fi(a-3,b,c,i)/8.0D0
!                interp_gix(i) = 3.0D0/8.0D0*gi(a+1,b,c,i)+3.0D0/4.0D0*gi(a-1,b,c,i)-gi(a-3,b,c,i)/8.0D0
!              end do
!!              interp_chi_x      = 3.0D0/8.0D0*chi(a+1,b,c,1)     +3.0D0/4.0D0*chi(a-1,b,c,1)     -chi(a-3,b,c,1)/8.0D0
!!              interp_zetax_x   = 3.0D0/8.0D0*zeta_x(a+1,b,c,1)  +3.0D0/4.0D0*zeta_x(a-1,b,c,1)  -zeta_x(a-3,b,c,1)/8.0D0
!!              interp_zetay_x   = 3.0D0/8.0D0*zeta_y(a+1,b,c,1)  +3.0D0/4.0D0*zeta_y(a-1,b,c,1)  -zeta_y(a-3,b,c,1)/8.0D0
!!              interp_zetaz_x   = 3.0D0/8.0D0*zeta_z(a+1,b,c,1)  +3.0D0/4.0D0*zeta_z(a-1,b,c,1)  -zeta_z(a-3,b,c,1)/8.0D0
!!              interp_pixx_x    = 3.0D0/8.0D0*pi_xx(a+1,b,c,1)   +3.0D0/4.0D0*pi_xx(a-1,b,c,1)   -pi_xx(a-3,b,c,1)/8.0D0
!!              interp_piyy_x    = 3.0D0/8.0D0*pi_yy(a+1,b,c,1)   +3.0D0/4.0D0*pi_yy(a-1,b,c,1)   -pi_yy(a-3,b,c,1)/8.0D0
!!              interp_pizz_x    = 3.0D0/8.0D0*pi_zz(a+1,b,c,1)   +3.0D0/4.0D0*pi_zz(a-1,b,c,1)   -pi_zz(a-3,b,c,1)/8.0D0
!!              interp_pixy_x    = 3.0D0/8.0D0*pi_xy(a+1,b,c,1)   +3.0D0/4.0D0*pi_xy(a-1,b,c,1)   -pi_xy(a-3,b,c,1)/8.0D0
!!              interp_pixz_x    = 3.0D0/8.0D0*pi_xz(a+1,b,c,1)   +3.0D0/4.0D0*pi_xz(a-1,b,c,1)   -pi_xz(a-3,b,c,1)/8.0D0
!!              interp_piyz_x    = 3.0D0/8.0D0*pi_yz(a+1,b,c,1)   +3.0D0/4.0D0*pi_yz(a-1,b,c,1)   -pi_yz(a-3,b,c,1)/8.0D0
!!              interp_lambdax_x = 3.0D0/8.0D0*lambda_x(a+1,b,c,1)+3.0D0/4.0D0*lambda_x(a-1,b,c,1)-lambda_x(a-3,b,c,1)/8.0D0
!!              interp_lambday_x = 3.0D0/8.0D0*lambda_y(a+1,b,c,1)+3.0D0/4.0D0*lambda_y(a-1,b,c,1)-lambda_y(a-3,b,c,1)/8.0D0
!!              interp_lambdaz_x = 3.0D0/8.0D0*lambda_z(a+1,b,c,1)+3.0D0/4.0D0*lambda_z(a-1,b,c,1)-lambda_z(a-3,b,c,1)/8.0D0
!            else if (x_down_clear .and. x_up_close .and. .not. x_up_clear) then
!              do i = 1,dir
!                interp_fix(i) = 3.0D0/8.0D0*fi(a-1,b,c,i)+3.0D0/4.0D0*fi(a+1,b,c,i)-fi(a+3,b,c,i)/8.0D0
!                interp_gix(i) = 3.0D0/8.0D0*gi(a-1,b,c,i)+3.0D0/4.0D0*gi(a+1,b,c,i)-gi(a+3,b,c,i)/8.0D0
!              end do
!!              interp_chi_x      = 3.0D0/8.0D0*chi(a-1,b,c,1)     +3.0D0/4.0D0*chi(a+1,b,c,1)     -chi(a+3,b,c,1)/8.0D0
!!              interp_zetax_x   = 3.0D0/8.0D0*zeta_x(a-1,b,c,1)  +3.0D0/4.0D0*zeta_x(a+1,b,c,1)  -zeta_x(a+3,b,c,1)/8.0D0
!!              interp_zetay_x   = 3.0D0/8.0D0*zeta_y(a-1,b,c,1)  +3.0D0/4.0D0*zeta_y(a+1,b,c,1)  -zeta_y(a+3,b,c,1)/8.0D0
!!              interp_zetaz_x   = 3.0D0/8.0D0*zeta_z(a-1,b,c,1)  +3.0D0/4.0D0*zeta_z(a+1,b,c,1)  -zeta_z(a+3,b,c,1)/8.0D0
!!              interp_pixx_x    = 3.0D0/8.0D0*pi_xx(a-1,b,c,1)   +3.0D0/4.0D0*pi_xx(a+1,b,c,1)   -pi_xx(a+3,b,c,1)/8.0D0
!!              interp_piyy_x    = 3.0D0/8.0D0*pi_yy(a-1,b,c,1)   +3.0D0/4.0D0*pi_yy(a+1,b,c,1)   -pi_yy(a+3,b,c,1)/8.0D0
!!              interp_pizz_x    = 3.0D0/8.0D0*pi_zz(a-1,b,c,1)   +3.0D0/4.0D0*pi_zz(a+1,b,c,1)   -pi_zz(a+3,b,c,1)/8.0D0
!!              interp_pixy_x    = 3.0D0/8.0D0*pi_xy(a-1,b,c,1)   +3.0D0/4.0D0*pi_xy(a+1,b,c,1)   -pi_xy(a+3,b,c,1)/8.0D0
!!              interp_pixz_x    = 3.0D0/8.0D0*pi_xz(a-1,b,c,1)   +3.0D0/4.0D0*pi_xz(a+1,b,c,1)   -pi_xz(a+3,b,c,1)/8.0D0
!!              interp_piyz_x    = 3.0D0/8.0D0*pi_yz(a-1,b,c,1)   +3.0D0/4.0D0*pi_yz(a+1,b,c,1)   -pi_yz(a+3,b,c,1)/8.0D0
!!              interp_lambdax_x = 3.0D0/8.0D0*lambda_x(a-1,b,c,1)+3.0D0/4.0D0*lambda_x(a+1,b,c,1)-lambda_x(a+3,b,c,1)/8.0D0
!!              interp_lambday_x = 3.0D0/8.0D0*lambda_y(a-1,b,c,1)+3.0D0/4.0D0*lambda_y(a+1,b,c,1)-lambda_y(a+3,b,c,1)/8.0D0
!!              interp_lambdaz_x = 3.0D0/8.0D0*lambda_z(a-1,b,c,1)+3.0D0/4.0D0*lambda_z(a+1,b,c,1)-lambda_z(a+3,b,c,1)/8.0D0
!            else if (x_up_close .and. x_down_close) then
!              do i=1,dir
!                interp_fix(i) = (fi(a+1,b,c,i)+fi(a-1,b,c,i))/2.0D0
!                interp_gix(i) = (gi(a+1,b,c,i)+gi(a-1,b,c,i))/2.0D0
!              end do
!!              interp_chi_x    =        (chi(a+1,b,c,1)+     chi(a-1,b,c,1))/2.0D0
!!              interp_zetax_x =     (zeta_x(a+1,b,c,1)+  zeta_x(a-1,b,c,1))/2.0D0
!!              interp_zetay_x =     (zeta_y(a+1,b,c,1)+  zeta_y(a-1,b,c,1))/2.0D0
!!              interp_zetaz_x =     (zeta_z(a+1,b,c,1)+  zeta_z(a-1,b,c,1))/2.0D0
!!              interp_pixx_x =       (pi_xx(a+1,b,c,1)+   pi_xx(a-1,b,c,1))/2.0D0
!!              interp_piyy_x =       (pi_yy(a+1,b,c,1)+   pi_yy(a-1,b,c,1))/2.0D0
!!              interp_pizz_x =       (pi_zz(a+1,b,c,1)+   pi_zz(a-1,b,c,1))/2.0D0
!!              interp_pixy_x =       (pi_xy(a+1,b,c,1)+   pi_xy(a-1,b,c,1))/2.0D0
!!              interp_pixz_x =       (pi_xz(a+1,b,c,1)+   pi_xz(a-1,b,c,1))/2.0D0
!!              interp_piyz_x =       (pi_yz(a+1,b,c,1)+   pi_yz(a-1,b,c,1))/2.0D0
!!              interp_lambdax_x = (lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))/2.0D0
!!              interp_lambday_x = (lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))/2.0D0
!!              interp_lambdaz_x = (lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))/2.0D0
!            else
!              !needs_interp(a,b,c) = .true.
!              do i =1,dir
!                interp_fix(i) = crs_fi(a,b,c,i)
!                interp_gix(i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!            if (.not. x_null .and. .not. y_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fix(i)+interp_fiy(i))/2.0D0
!                gi(a,b,c,i) = (interp_gix(i)+interp_giy(i))/2.0D0
!!                if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!!                  write(*,*) 'big badda boomxy',fi(a,b,c,i),gi(a,b,c,i),a,b,c,i
!!                end if
!              end do
!!              chi(a,b,c,1)    = (interp_chi_x + interp_chi_y)/2.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_x + interp_zetax_y)/2.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_x + interp_zetay_y)/2.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_x + interp_zetaz_y)/2.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_x + interp_pixx_y)/2.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_x + interp_piyy_y)/2.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_x + interp_pizz_y)/2.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_x + interp_pixy_y)/2.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_x + interp_pixz_y)/2.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_x + interp_piyz_y)/2.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_x + interp_lambdax_y)/2.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_x + interp_lambday_y)/2.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_x + interp_lambdaz_y)/2.0D0
!            else if (x_null .and. .not. y_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fiy(i)
!                gi(a,b,c,i) = interp_giy(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_y
!!              zeta_x(a,b,c,1) =  interp_zetax_y
!!              zeta_y(a,b,c,1) =  interp_zetay_y
!!              zeta_z(a,b,c,1) =  interp_zetaz_y
!!              pi_xx(a,b,c,1) =  interp_pixx_y
!!              pi_yy(a,b,c,1) =  interp_piyy_y
!!              pi_zz(a,b,c,1) =  interp_pizz_y
!!              pi_xy(a,b,c,1) =  interp_pixy_y
!!              pi_xz(a,b,c,1) =  interp_pixz_y
!!              pi_yz(a,b,c,1) =  interp_piyz_y
!!              lambda_x(a,b,c,1) =interp_lambdax_y
!!              lambda_y(a,b,c,1) =interp_lambday_y
!!              lambda_z(a,b,c,1) =interp_lambdaz_y
!            else if (y_null .and. .not. x_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fix(i)
!                gi(a,b,c,i) = interp_gix(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_x
!!              zeta_x(a,b,c,1) =  interp_zetax_x
!!              zeta_y(a,b,c,1) =  interp_zetay_x
!!              zeta_z(a,b,c,1) =  interp_zetaz_x
!!              pi_xx(a,b,c,1) =  interp_pixx_x
!!              pi_yy(a,b,c,1) =  interp_piyy_x
!!              pi_zz(a,b,c,1) =  interp_pizz_x
!!              pi_xy(a,b,c,1) =  interp_pixy_x
!!              pi_xz(a,b,c,1) =  interp_pixz_x
!!              pi_yz(a,b,c,1) =  interp_piyz_x
!!              lambda_x(a,b,c,1) =interp_lambdax_x
!!              lambda_y(a,b,c,1) =interp_lambday_x
!!              lambda_z(a,b,c,1) =interp_lambdaz_x
!            else
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!!            if (a == 479 .and. b == 239 .and. c == 72) then
!!              write(*,*) 'checking values for ghost node',fi(a,b,c,1:dir),a,b,c
!!              write(*,*) 'checking gis',gi(a,b,c,1:dir),a,b,c
!!              write(*,*) 'checking fix interp',interp_fix,a,b,c
!!              write(*,*) 'checking gix interp',interp_gix,a,b,c
!!              write(*,*) 'checking fiy interp',interp_fiy,a,b,c
!!              write(*,*) 'checking giy interp',interp_giy,a,b,c
!!              write(*,*) 'interp states',state(a,b,c,1:dir),a,b,c
!!              write(*,*) 'checking bools',x_up_clear,x_up_close,x_down_clear,x_down_close,&
!!                y_up_clear,y_up_close,y_down_clear,y_down_close,a,b,c
!!            end if
!!
!! Face in xz
!!
!          else if (mod(c,2) /= 0 .and. mod(b,2) == 0 .and. mod(a,2) /= 0) then
!            x_null = .false.
!            y_null = .false.
!            z_null = .false.
!
!            if (state(a,b,c,5) >= 0 .and. state(a,b,c,37) >= 0) then
!              x_up_clear = .true.
!              x_up_close = .true.
!            else if (state(a,b,c,5) >= 0) then
!              x_up_close = .true.
!              x_up_clear = .false.
!            else
!              x_up_close = .false.
!              x_up_clear = .false.
!            end if
!            if (state(a,b,c,2) >= 0 .and. state(a,b,c,34) >= 0) then
!              x_down_clear = .true.
!              x_down_close = .true.
!            else if (state(a,b,c,2) >= 0) then
!              x_down_close = .true.
!              x_down_clear = .false.
!            else
!              x_down_close = .false.
!              x_down_clear = .false.
!            end if
!
!            if (state(a,b,c,7) >= 0 .and. state(a,b,c,39) >= 0) then
!              z_up_clear = .true.
!              z_up_close = .true.
!            else if (state(a,b,c,7) >= 0) then
!              z_up_close = .true.
!              z_up_clear = .false.
!            else
!              z_up_clear = .false.
!              z_up_close = .false.
!            end if
!            if (state(a,b,c,4) >= 0 .and. state(a,b,c,36) >= 0) then
!              z_down_clear = .true.
!              z_down_close = .true.
!            else if (state(a,b,c,4) >= 0) then
!              z_down_close = .true.
!              z_down_clear = .false.
!            else
!              z_down_clear = .false.
!              z_down_close = .false.
!            end if
!!
!!
!!
!            if (.not. z_down_close .or. .not. z_up_close) then
!              z_null = .true.
!            end if
!            if (.not. x_down_close .or. .not. x_up_close) then
!              x_null = .true.
!            end if
!!
!!
!!
!            if (x_up_clear .and. x_down_clear ) then
!              do i = 1,dir
!                interp_fix(i) = (9.0D0*(fi(a+1,b,c,i)+fi(a-1,b,c,i))-&
!                                 fi(a-3,b,c,i)-fi(a+3,b,c,i))/16.0D0
!
!                interp_gix(i) = (9.0D0*(gi(a+1,b,c,i)+gi(a-1,b,c,i))-&
!                                 gi(a+3,b,c,i)-gi(a-3,b,c,i))/16.0D0
!              end do
!!              interp_chi_x = (9.0D0*(chi(a+1,b,c,1)+chi(a-1,b,c,1))-&
!!                chi(a+3,b,c,1)-chi(a-3,b,c,1))/16.0D0
!!              interp_zetax_x = (9.0D0*(zeta_x(a+1,b,c,1)+zeta_x(a-1,b,c,1))-&
!!                zeta_x(a+3,b,c,1)-zeta_x(a-3,b,c,1))/16.0D0
!!              interp_zetay_x = (9.0D0*(zeta_y(a+1,b,c,1)+zeta_y(a-1,b,c,1))-&
!!                zeta_y(a+3,b,c,1)-zeta_y(a-3,b,c,1))/16.0D0
!!              interp_zetaz_x = (9.0D0*(zeta_z(a+1,b,c,1)+zeta_z(a-1,b,c,1))-&
!!                zeta_z(a+3,b,c,1)-zeta_z(a-3,b,c,1))/16.0D0
!!              interp_pixx_x = (9.0D0*(pi_xx(a+1,b,c,1)+pi_xx(a-1,b,c,1))-&
!!                pi_xx(a+3,b,c,1)-pi_xx(a-3,b,c,1))/16.0D0
!!              interp_piyy_x = (9.0D0*(pi_yy(a+1,b,c,1)+pi_yy(a-1,b,c,1))-&
!!                pi_yy(a-3,b,c,1)-pi_yy(a+3,b,c,1))/16.0D0
!!              interp_pizz_x = (9.0D0*(pi_zz(a+1,b,c,1)+pi_zz(a-1,b,c,1))-&
!!                pi_zz(a+3,b,c,1)-pi_zz(a-3,b,c,1))/16.0D0
!!              interp_pixy_x = (9.0D0*(pi_xy(a+1,b,c,1)+pi_xy(a-1,b,c,1))-&
!!                pi_xy(a+3,b,c,1)-pi_xy(a-3,b,c,1))/16.0D0
!!              interp_pixz_x = (9.0D0*(pi_xz(a+1,b,c,1)+pi_xz(a-1,b,c,1))-&
!!                pi_xz(a+3,b,c,1)-pi_xz(a-3,b,c,1))/16.0D0
!!              interp_piyz_x = (9.0D0*(pi_yz(a+1,b,c,1)+pi_yz(a-1,b,c,1))-&
!!                pi_yz(a+3,b,c,1)-pi_yz(a-3,b,c,1))/16.0D0
!!              interp_lambdax_x = (9.0D0*(lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))-&
!!                lambda_x(a+3,b,c,1)-lambda_x(a-3,b,c,1))/16.0D0
!!              interp_lambday_x = (9.0D0*(lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))-&
!!                lambda_y(a+3,b,c,1)-lambda_y(a-3,b,c,1))/16.0D0
!!              interp_lambdaz_x = (9.0D0*(lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))-&
!!                lambda_z(a+3,b,c,1)-lambda_z(a-3,b,c,1))/16.0D0
!            else if (x_up_clear .and. x_down_close .and. .not. x_down_clear) then
!              do i = 1,dir
!                interp_fix(i) = 3.0D0/8.0D0*fi(a+1,b,c,i)+3.0D0/4.0D0*fi(a-1,b,c,i)-fi(a-3,b,c,i)/8.0D0
!                interp_gix(i) = 3.0D0/8.0D0*gi(a+1,b,c,i)+3.0D0/4.0D0*gi(a-1,b,c,i)-gi(a-3,b,c,i)/8.0D0
!              end do
!!              interp_chi_x      = 3.0D0/8.0D0*chi(a+1,b,c,1)     +3.0D0/4.0D0*chi(a-1,b,c,1)     -chi(a-3,b,c,1)/8.0D0
!!              interp_zetax_x   = 3.0D0/8.0D0*zeta_x(a+1,b,c,1)  +3.0D0/4.0D0*zeta_x(a-1,b,c,1)  -zeta_x(a-3,b,c,1)/8.0D0
!!              interp_zetay_x   = 3.0D0/8.0D0*zeta_y(a+1,b,c,1)  +3.0D0/4.0D0*zeta_y(a-1,b,c,1)  -zeta_y(a-3,b,c,1)/8.0D0
!!              interp_zetaz_x   = 3.0D0/8.0D0*zeta_z(a+1,b,c,1)  +3.0D0/4.0D0*zeta_z(a-1,b,c,1)  -zeta_z(a-3,b,c,1)/8.0D0
!!              interp_pixx_x    = 3.0D0/8.0D0*pi_xx(a+1,b,c,1)   +3.0D0/4.0D0*pi_xx(a-1,b,c,1)   -pi_xx(a-3,b,c,1)/8.0D0
!!              interp_piyy_x    = 3.0D0/8.0D0*pi_yy(a+1,b,c,1)   +3.0D0/4.0D0*pi_yy(a-1,b,c,1)   -pi_yy(a-3,b,c,1)/8.0D0
!!              interp_pizz_x    = 3.0D0/8.0D0*pi_zz(a+1,b,c,1)   +3.0D0/4.0D0*pi_zz(a-1,b,c,1)   -pi_zz(a-3,b,c,1)/8.0D0
!!              interp_pixy_x    = 3.0D0/8.0D0*pi_xy(a+1,b,c,1)   +3.0D0/4.0D0*pi_xy(a-1,b,c,1)   -pi_xy(a-3,b,c,1)/8.0D0
!!              interp_pixz_x    = 3.0D0/8.0D0*pi_xz(a+1,b,c,1)   +3.0D0/4.0D0*pi_xz(a-1,b,c,1)   -pi_xz(a-3,b,c,1)/8.0D0
!!              interp_piyz_x    = 3.0D0/8.0D0*pi_yz(a+1,b,c,1)   +3.0D0/4.0D0*pi_yz(a-1,b,c,1)   -pi_yz(a-3,b,c,1)/8.0D0
!!              interp_lambdax_x = 3.0D0/8.0D0*lambda_x(a+1,b,c,1)+3.0D0/4.0D0*lambda_x(a-1,b,c,1)-lambda_x(a-3,b,c,1)/8.0D0
!!              interp_lambday_x = 3.0D0/8.0D0*lambda_y(a+1,b,c,1)+3.0D0/4.0D0*lambda_y(a-1,b,c,1)-lambda_y(a-3,b,c,1)/8.0D0
!!              interp_lambdaz_x = 3.0D0/8.0D0*lambda_z(a+1,b,c,1)+3.0D0/4.0D0*lambda_z(a-1,b,c,1)-lambda_z(a-3,b,c,1)/8.0D0
!            else if (x_down_clear .and. x_up_close .and. .not. x_up_clear) then
!              do i = 1,dir
!                interp_fix(i) = 3.0D0/8.0D0*fi(a-1,b,c,i)+3.0D0/4.0D0*fi(a+1,b,c,i)-fi(a+3,b,c,i)/8.0D0
!                interp_gix(i) = 3.0D0/8.0D0*gi(a-1,b,c,i)+3.0D0/4.0D0*gi(a+1,b,c,i)-gi(a+3,b,c,i)/8.0D0
!              end do
!!              interp_chi_x      = 3.0D0/8.0D0*chi(a-1,b,c,1)     +3.0D0/4.0D0*chi(a+1,b,c,1)     -chi(a+3,b,c,1)/8.0D0
!!              interp_zetax_x   = 3.0D0/8.0D0*zeta_x(a-1,b,c,1)  +3.0D0/4.0D0*zeta_x(a+1,b,c,1)  -zeta_x(a+3,b,c,1)/8.0D0
!!              interp_zetay_x   = 3.0D0/8.0D0*zeta_y(a-1,b,c,1)  +3.0D0/4.0D0*zeta_y(a+1,b,c,1)  -zeta_y(a+3,b,c,1)/8.0D0
!!              interp_zetaz_x   = 3.0D0/8.0D0*zeta_z(a-1,b,c,1)  +3.0D0/4.0D0*zeta_z(a+1,b,c,1)  -zeta_z(a+3,b,c,1)/8.0D0
!!              interp_pixx_x    = 3.0D0/8.0D0*pi_xx(a-1,b,c,1)   +3.0D0/4.0D0*pi_xx(a+1,b,c,1)   -pi_xx(a+3,b,c,1)/8.0D0
!!              interp_piyy_x    = 3.0D0/8.0D0*pi_yy(a-1,b,c,1)   +3.0D0/4.0D0*pi_yy(a+1,b,c,1)   -pi_yy(a+3,b,c,1)/8.0D0
!!              interp_pizz_x    = 3.0D0/8.0D0*pi_zz(a-1,b,c,1)   +3.0D0/4.0D0*pi_zz(a+1,b,c,1)   -pi_zz(a+3,b,c,1)/8.0D0
!!              interp_pixy_x    = 3.0D0/8.0D0*pi_xy(a-1,b,c,1)   +3.0D0/4.0D0*pi_xy(a+1,b,c,1)   -pi_xy(a+3,b,c,1)/8.0D0
!!              interp_pixz_x    = 3.0D0/8.0D0*pi_xz(a-1,b,c,1)   +3.0D0/4.0D0*pi_xz(a+1,b,c,1)   -pi_xz(a+3,b,c,1)/8.0D0
!!              interp_piyz_x    = 3.0D0/8.0D0*pi_yz(a-1,b,c,1)   +3.0D0/4.0D0*pi_yz(a+1,b,c,1)   -pi_yz(a+3,b,c,1)/8.0D0
!!              interp_lambdax_x = 3.0D0/8.0D0*lambda_x(a-1,b,c,1)+3.0D0/4.0D0*lambda_x(a+1,b,c,1)-lambda_x(a+3,b,c,1)/8.0D0
!!              interp_lambday_x = 3.0D0/8.0D0*lambda_y(a-1,b,c,1)+3.0D0/4.0D0*lambda_y(a+1,b,c,1)-lambda_y(a+3,b,c,1)/8.0D0
!!              interp_lambdaz_x = 3.0D0/8.0D0*lambda_z(a-1,b,c,1)+3.0D0/4.0D0*lambda_z(a+1,b,c,1)-lambda_z(a+3,b,c,1)/8.0D0
!            else if (x_down_close .and. x_up_close) then
!              do i=1,dir
!                interp_fix(i) = (fi(a+1,b,c,i)+fi(a-1,b,c,i))/2.0D0
!                interp_gix(i) = (gi(a+1,b,c,i)+gi(a-1,b,c,i))/2.0D0
!              end do
!!              interp_chi_x    =        (chi(a+1,b,c,1)+     chi(a-1,b,c,1))/2.0D0
!!              interp_zetax_x =     (zeta_x(a+1,b,c,1)+  zeta_x(a-1,b,c,1))/2.0D0
!!              interp_zetay_x =     (zeta_y(a+1,b,c,1)+  zeta_y(a-1,b,c,1))/2.0D0
!!              interp_zetaz_x =     (zeta_z(a+1,b,c,1)+  zeta_z(a-1,b,c,1))/2.0D0
!!              interp_pixx_x =       (pi_xx(a+1,b,c,1)+   pi_xx(a-1,b,c,1))/2.0D0
!!              interp_piyy_x =       (pi_yy(a+1,b,c,1)+   pi_yy(a-1,b,c,1))/2.0D0
!!              interp_pizz_x =       (pi_zz(a+1,b,c,1)+   pi_zz(a-1,b,c,1))/2.0D0
!!              interp_pixy_x =       (pi_xy(a+1,b,c,1)+   pi_xy(a-1,b,c,1))/2.0D0
!!              interp_pixz_x =       (pi_xz(a+1,b,c,1)+   pi_xz(a-1,b,c,1))/2.0D0
!!              interp_piyz_x =       (pi_yz(a+1,b,c,1)+   pi_yz(a-1,b,c,1))/2.0D0
!!              interp_lambdax_x = (lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))/2.0D0
!!              interp_lambday_x = (lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))/2.0D0
!!              interp_lambdaz_x = (lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))/2.0D0
!            else
!              do i =1,dir
!                interp_fix(i) = crs_fi(a,b,c,i)
!                interp_gix(i) = crs_gi(a,b,c,i)
!              end do
!
!            end if
!
!            if (z_up_clear .and. z_down_clear ) then
!              do i = 1,dir
!                interp_fiz(i) = (9.0D0*(fi(a,b,c+1,i)+fi(a,b,c-1,i))-&
!                                 fi(a,b,c-3,i)-fi(a,b,c+3,i))/16.0D0
!
!                interp_giz(i) = (9.0D0*(gi(a,b,c+1,i)+gi(a,b,c-1,i))-&
!                                 gi(a,b,c+3,i)-gi(a,b,c-3,i))/16.0D0
!              end do
!!              interp_chi_z = (9.0D0*(chi(a,b,c+1,1)+chi(a,b,c-1,1))-&
!!                chi(a,b,c+3,1)-chi(a,b,c-3,1))/16.0D0
!!              interp_zetax_z = (9.0D0*(zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))-&
!!                zeta_x(a,b,c+3,1)-zeta_x(a,b,c-3,1))/16.0D0
!!              interp_zetay_z = (9.0D0*(zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))-&
!!                zeta_y(a,b,c+3,1)-zeta_y(a,b,c-3,1))/16.0D0
!!              interp_zetaz_z = (9.0D0*(zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))-&
!!                zeta_z(a,b,c+3,1)-zeta_z(a,b,c-3,1))/16.0D0
!!              interp_pixx_z = (9.0D0*(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))-&
!!                pi_xx(a,b,c+3,1)-pi_xx(a,b,c-3,1))/16.0D0
!!              interp_piyy_z = (9.0D0*(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))-&
!!                pi_yy(a,b,c+3,1)-pi_yy(a,b,c-3,1))/16.0D0
!!              interp_pizz_z = (9.0D0*(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))-&
!!                pi_zz(a,b,c+3,1)-pi_zz(a,b,c-3,1))/16.0D0
!!              interp_pixy_z = (9.0D0*(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))-&
!!                pi_xy(a,b,c+3,1)-pi_xy(a,b,c-3,1))/16.0D0
!!              interp_pixz_z = (9.0D0*(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))-&
!!                pi_xz(a,b,c+3,1)-pi_xz(a,b,c-3,1))/16.0D0
!!              interp_piyz_z = (9.0D0*(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))-&
!!                pi_yz(a,b,c+3,1)-pi_yz(a,b,c-3,1))/16.0D0
!!              interp_lambdax_z = (9.0D0*(lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))-&
!!                lambda_x(a,b,c+3,1)-lambda_x(a,b,c-3,1))/16.0D0
!!              interp_lambday_z = (9.0D0*(lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))-&
!!                lambda_y(a,b,c+3,1)-lambda_y(a,b,c-3,1))/16.0D0
!!              interp_lambdaz_z = (9.0D0*(lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))-&
!!                lambda_z(a,b,c+3,1)-lambda_z(a,b,c-3,1))/16.0D0
!            else if (z_up_clear .and. z_down_close .and. .not. z_down_clear) then
!              do i = 1,dir
!                interp_fiz(i) = 3.0D0/8.0D0*fi(a,b,c+1,i)+3.0D0/4.0D0*fi(a,b,c-1,i)-fi(a,b,c-3,i)/8.0D0
!                interp_giz(i) = 3.0D0/8.0D0*gi(a,b,c+1,i)+3.0D0/4.0D0*gi(a,b,c-1,i)-gi(a,b,c-3,i)/8.0D0
!              end do
!!              interp_chi_z =     3.0D0/8.0D0*chi(a,b,c+1,1)+3.0D0/4.0D0*chi(a,b,c-1,1)-chi(a,b,c-3,1)/8.0D0
!!              interp_zetax_z =   3.0D0/8.0D0*zeta_x(a,b,c+1,1)+3.0D0/4.0D0*zeta_x(a,b,c-1,1)-zeta_x(a,b,c-3,1)/8.0D0
!!              interp_zetay_z =   3.0D0/8.0D0*zeta_y(a,b,c+1,1)+3.0D0/4.0D0*zeta_y(a,b,c-1,1)-zeta_y(a,b,c-3,1)/8.0D0
!!              interp_zetaz_z =   3.0D0/8.0D0*zeta_z(a,b,c+1,1)+3.0D0/4.0D0*zeta_z(a,b,c-1,1)-zeta_z(a,b,c-3,1)/8.0D0
!!              interp_pixx_z =    3.0D0/8.0D0*pi_xx(a,b,c+1,1)+3.0D0/4.0D0*pi_xx(a,b,c-1,1)-pi_xx(a,b,c-3,1)/8.0D0
!!              interp_piyy_z =    3.0D0/8.0D0*pi_yy(a,b,c+1,1)+3.0D0/4.0D0*pi_yy(a,b,c-1,1)-pi_yy(a,b,c-3,1)/8.0D0
!!              interp_pizz_z =    3.0D0/8.0D0*pi_zz(a,b,c+1,1)+3.0D0/4.0D0*pi_zz(a,b,c-1,1)-pi_zz(a,b,c-3,1)/8.0D0
!!              interp_pixy_z =    3.0D0/8.0D0*pi_xy(a,b,c+1,1)+3.0D0/4.0D0*pi_xy(a,b,c-1,1)-pi_xy(a,b,c-3,1)/8.0D0
!!              interp_pixz_z =    3.0D0/8.0D0*pi_xz(a,b,c+1,1)+3.0D0/4.0D0*pi_xz(a,b,c-1,1)-pi_xz(a,b,c-3,1)/8.0D0
!!              interp_piyz_z =    3.0D0/8.0D0*pi_yz(a,b,c+1,1)+3.0D0/4.0D0*pi_yz(a,b,c-1,1)-pi_yz(a,b,c-3,1)/8.0D0
!!              interp_lambdax_z = 3.0D0/8.0D0*lambda_x(a,b,c+1,1)+3.0D0/4.0D0*lambda_x(a,b,c-1,1)-lambda_x(a,b,c-3,1)/8.0D0
!!              interp_lambday_z = 3.0D0/8.0D0*lambda_y(a,b,c+1,1)+3.0D0/4.0D0*lambda_y(a,b,c-1,1)-lambda_y(a,b,c-3,1)/8.0D0
!!              interp_lambdaz_z = 3.0D0/8.0D0*lambda_z(a,b,c+1,1)+3.0D0/4.0D0*lambda_z(a,b,c-1,1)-lambda_z(a,b,c-3,1)/8.0D0
!            else if (z_down_clear .and. z_up_close .and. .not. z_up_clear) then
!              do i = 1,dir
!                interp_fiz(i) = 3.0D0/8.0D0*fi(a,b,c-1,i)+3.0D0/4.0D0*fi(a,b,c+1,i)-fi(a,b,c+3,i)/8.0D0
!                interp_giz(i) = 3.0D0/8.0D0*gi(a,b,c-1,i)+3.0D0/4.0D0*gi(a,b,c+1,i)-gi(a,b,c+3,i)/8.0D0
!              end do
!!              interp_chi_z =      3.0D0/8.0D0*chi(a,b,c-1,1)+3.0D0/4.0D0*chi(a,b,c+1,1)-chi(a,b,c+3,1)/8.0D0
!!              interp_zetax_z =   3.0D0/8.0D0*zeta_x(a,b,c-1,1)+3.0D0/4.0D0*zeta_x(a,b,c+1,1)-zeta_x(a,b,c+3,1)/8.0D0
!!              interp_zetay_z =   3.0D0/8.0D0*zeta_y(a,b,c-1,1)+3.0D0/4.0D0*zeta_y(a,b,c+1,1)-zeta_y(a,b,c+3,1)/8.0D0
!!              interp_zetaz_z =   3.0D0/8.0D0*zeta_z(a,b,c-1,1)+3.0D0/4.0D0*zeta_z(a,b,c+1,1)-zeta_z(a,b,c+3,1)/8.0D0
!!              interp_pixx_z =    3.0D0/8.0D0*pi_xx(a,b,c-1,1)+3.0D0/4.0D0*pi_xx(a,b,c+1,1)-pi_xx(a,b,c+3,1)/8.0D0
!!              interp_piyy_z =    3.0D0/8.0D0*pi_yy(a,b,c-1,1)+3.0D0/4.0D0*pi_yy(a,b,c+1,1)-pi_yy(a,b,c+3,1)/8.0D0
!!              interp_pizz_z =    3.0D0/8.0D0*pi_zz(a,b,c-1,1)+3.0D0/4.0D0*pi_zz(a,b,c+1,1)-pi_zz(a,b,c+3,1)/8.0D0
!!              interp_pixy_z =    3.0D0/8.0D0*pi_xy(a,b,c-1,1)+3.0D0/4.0D0*pi_xy(a,b,c+1,1)-pi_xy(a,b,c+3,1)/8.0D0
!!              interp_pixz_z =    3.0D0/8.0D0*pi_xz(a,b,c-1,1)+3.0D0/4.0D0*pi_xz(a,b,c+1,1)-pi_xz(a,b,c+3,1)/8.0D0
!!              interp_piyz_z =    3.0D0/8.0D0*pi_yz(a,b,c-1,1)+3.0D0/4.0D0*pi_yz(a,b,c+1,1)-pi_yz(a,b,c+3,1)/8.0D0
!!              interp_lambdax_z= 3.0D0/8.0D0*lambda_x(a,b,c-1,1)+3.0D0/4.0D0*lambda_x(a,b,c+1,1)-lambda_x(a,b,c+3,1)/8.0D0
!!              interp_lambday_z= 3.0D0/8.0D0*lambda_y(a,b,c-1,1)+3.0D0/4.0D0*lambda_y(a,b,c+1,1)-lambda_y(a,b,c+3,1)/8.0D0
!!              interp_lambdaz_z= 3.0D0/8.0D0*lambda_z(a,b,c-1,1)+3.0D0/4.0D0*lambda_z(a,b,c+1,1)-lambda_z(a,b,c+3,1)/8.0D0
!            else if (z_up_close .and. z_down_close) then
!              do i=1,dir
!                interp_fiz(i) = (fi(a,b,c+1,i)+fi(a,b,c-1,i))/2.0D0
!                interp_giz(i) = (gi(a,b,c+1,i)+gi(a,b,c-1,i))/2.0D0
!              end do
!!              interp_chi_z = (chi(a,b,c+1,1)+chi(a,b,c-1,1))/2.0D0
!!              interp_zetax_z = (zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))/2.0D0
!!              interp_zetay_z = (zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))/2.0D0
!!              interp_zetaz_z = (zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))/2.0D0
!!              interp_pixx_z =(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))/2.0D0
!!              interp_piyy_z =(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))/2.0D0
!!              interp_pizz_z =(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))/2.0D0
!!              interp_pixy_z =(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))/2.0D0
!!              interp_pixz_z =(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))/2.0D0
!!              interp_piyz_z =(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))/2.0D0
!!              interp_lambdax_z = (lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))/2.0D0
!!              interp_lambday_z = (lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))/2.0D0
!!              interp_lambdaz_z = (lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))/2.0D0
!            else
!              do i =1,dir
!                interp_fiz(i) = crs_fi(a,b,c,i)
!                interp_giz(i) = crs_gi(a,b,c,i)
!              end do
!
!            end if
!
!            if (.not. x_null .and. .not. z_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fix(i)+interp_fiz(i))/2.0D0
!                gi(a,b,c,i) = (interp_gix(i)+interp_giz(i))/2.0D0
!!                if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!!                  write(*,*) 'big badda boomxz',fi(a,b,c,i),gi(a,b,c,i),a,b,c,i
!!                end if
!              end do
!!              chi(a,b,c,1)    = (interp_chi_z + interp_chi_x)/2.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_z + interp_zetax_x)/2.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_z + interp_zetay_x)/2.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_z + interp_zetaz_x)/2.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_z + interp_pixx_x)/2.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_z + interp_piyy_x)/2.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_z + interp_pizz_x)/2.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_z + interp_pixy_x)/2.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_z + interp_pixz_x)/2.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_z + interp_piyz_x)/2.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_z + interp_lambdax_x)/2.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_z + interp_lambday_x)/2.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_z + interp_lambdaz_x)/2.0D0
!            else if (z_null .and. .not. x_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fix(i)
!                gi(a,b,c,i) = interp_gix(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_x
!!              zeta_x(a,b,c,1) =  interp_zetax_x
!!              zeta_y(a,b,c,1) =  interp_zetay_x
!!              zeta_z(a,b,c,1) =  interp_zetaz_x
!!              pi_xx(a,b,c,1) =  interp_pixx_x
!!              pi_yy(a,b,c,1) =  interp_piyy_x
!!              pi_zz(a,b,c,1) =  interp_pizz_x
!!              pi_xy(a,b,c,1) =  interp_pixy_x
!!              pi_xz(a,b,c,1) =  interp_pixz_x
!!              pi_yz(a,b,c,1) =  interp_piyz_x
!!              lambda_x(a,b,c,1) =interp_lambdax_x
!!              lambda_y(a,b,c,1) =interp_lambday_x
!!              lambda_z(a,b,c,1) =interp_lambdaz_x
!            else if (x_null .and. .not. z_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fiz(i)
!                gi(a,b,c,i) = interp_giz(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_z
!!              zeta_x(a,b,c,1) =  interp_zetax_z
!!              zeta_y(a,b,c,1) =  interp_zetay_z
!!              zeta_z(a,b,c,1) =  interp_zetaz_z
!!              pi_xx(a,b,c,1) =  interp_pixx_z
!!              pi_yy(a,b,c,1) =  interp_piyy_z
!!              pi_zz(a,b,c,1) =  interp_pizz_z
!!              pi_xy(a,b,c,1) =  interp_pixy_z
!!              pi_xz(a,b,c,1) =  interp_pixz_z
!!              pi_yz(a,b,c,1) =  interp_piyz_z
!!              lambda_x(a,b,c,1) =interp_lambdax_z
!!              lambda_y(a,b,c,1) =interp_lambday_z
!!              lambda_z(a,b,c,1) =interp_lambdaz_z
!            else
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!
!            end if
!
!
!            !write(*,*) 'values out on faces',fi(a,b,c,1),gi(a,b,c,1),a,b,c,x_up_clear,x_down_clear,&
!            !  y_up_clear,y_down_clear,z_up_clear,z_down_clear
!          end if
!
!          do i = 1,dir
!            if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!              fi(a,b,c,i) = crs_fi(a,b,c,i)
!              gi(a,b,c,i) = crs_gi(a,b,c,i)
!            end if
!          end do
!
!          end if ! checks if node is a true ghost node
!
!        end do
!      end do
!    end do
!  !write(*,*) 'arf arf arf',self
!    do c = ducks%lo(3)-nghosts_mid,ducks%hi(3)+nghosts_mid
!      do b = ducks%lo(2)-nghosts_mid,ducks%hi(2)+nghosts_mid
!        do a = ducks%lo(1)-nghosts_mid,ducks%hi(1)+nghosts_mid
!          if (state(a,b,c,1) >=10000000 .and. state(a,b,c,1) <= 1100000000) then
!            if (mod(a,2) /= 0 .and. mod(b,2) /= 0 .and. mod(c,2) /= 0) then
!
!            x_null = .false.
!            y_null = .false.
!            z_null = .false.
!
!            if (state(a,b,c,5) >= 0 .and. state(a,b,c,37) >= 0) then
!              x_up_clear = .true.
!              x_up_close = .true.
!            else if (state(a,b,c,5) >= 0) then
!              x_up_close = .true.
!              x_up_clear = .false.
!            else
!              x_up_close = .false.
!              x_up_clear = .false.
!            end if
!            if (state(a,b,c,2) >= 0 .and. state(a,b,c,34) >= 0) then
!              x_down_clear = .true.
!              x_down_close = .true.
!            else if (state(a,b,c,2) >= 0) then
!              x_down_close = .true.
!              x_down_clear = .false.
!            else
!              x_down_close = .false.
!              x_down_clear = .false.
!            end if
!
!            if (state(a,b,c,6) >= 0 .and. state(a,b,c,38) >= 0) then
!              y_up_clear = .true.
!              y_up_close = .true.
!            else if (state(a,b,c,6) < 0) then
!              y_up_close = .false.
!              y_up_clear = .false.
!            else
!              y_up_close = .true.
!              y_up_clear = .false.
!            end if
!            if (state(a,b,c,3) >= 0 .and. state(a,b,c,35) >= 0) then
!              y_down_clear = .true.
!              y_down_close = .true.
!            else if (state(a,b,c,3) < 0) then
!              y_down_close = .false.
!              y_down_clear = .false.
!            else
!              y_down_close = .true.
!              y_down_clear = .false.
!            end if
!
!            if (state(a,b,c,7) >= 0 .and. state(a,b,c,39) >= 0) then
!              z_up_clear = .true.
!              z_up_close = .true.
!            else if (state(a,b,c,7) >= 0) then
!              z_up_close = .true.
!              z_up_clear = .false.
!            else
!              z_up_clear = .false.
!              z_up_close = .false.
!            end if
!            if (state(a,b,c,4) >= 0 .and. state(a,b,c,36) >= 0) then
!              z_down_clear = .true.
!              z_down_close = .true.
!            else if (state(a,b,c,4) >= 0) then
!              z_down_close = .true.
!              z_down_clear = .false.
!            else
!              z_down_clear = .false.
!              z_down_close = .false.
!            end if
!
!            if (.not. x_down_close .or. .not. x_up_close) then
!              x_null = .true.
!            end if
!            if (.not. y_down_close .or. .not. y_up_close) then
!              y_null = .true.
!            end if
!            if (.not. z_down_close .or. .not. z_up_close) then
!              z_null = .true.
!            end if
!!
!            if (x_up_clear .and. x_down_clear ) then
!              do i = 1,dir
!                interp_fix(i) = (9.0D0*(fi(a+1,b,c,i)+fi(a-1,b,c,i))-&
!                                 fi(a-3,b,c,i)-fi(a+3,b,c,i))/16.0D0
!
!                interp_gix(i) = (9.0D0*(gi(a+1,b,c,i)+gi(a-1,b,c,i))-&
!                                 gi(a+3,b,c,i)-gi(a-3,b,c,i))/16.0D0
!              end do
!!              interp_chi_x = (9.0D0*(chi(a+1,b,c,1)+chi(a-1,b,c,1))-&
!!                chi(a+3,b,c,1)-chi(a-3,b,c,1))/16.0D0
!!              interp_zetax_x = (9.0D0*(zeta_x(a+1,b,c,1)+zeta_x(a-1,b,c,1))-&
!!                zeta_x(a+3,b,c,1)-zeta_x(a-3,b,c,1))/16.0D0
!!              interp_zetay_x = (9.0D0*(zeta_y(a+1,b,c,1)+zeta_y(a-1,b,c,1))-&
!!                zeta_y(a+3,b,c,1)-zeta_y(a-3,b,c,1))/16.0D0
!!              interp_zetaz_x = (9.0D0*(zeta_z(a+1,b,c,1)+zeta_z(a-1,b,c,1))-&
!!                zeta_z(a+3,b,c,1)-zeta_z(a-3,b,c,1))/16.0D0
!!              interp_pixx_x = (9.0D0*(pi_xx(a+1,b,c,1)+pi_xx(a-1,b,c,1))-&
!!                pi_xx(a+3,b,c,1)-pi_xx(a-3,b,c,1))/16.0D0
!!              interp_piyy_x = (9.0D0*(pi_yy(a+1,b,c,1)+pi_yy(a-1,b,c,1))-&
!!                pi_yy(a-3,b,c,1)-pi_yy(a+3,b,c,1))/16.0D0
!!              interp_pizz_x = (9.0D0*(pi_zz(a+1,b,c,1)+pi_zz(a-1,b,c,1))-&
!!                pi_zz(a+3,b,c,1)-pi_zz(a-3,b,c,1))/16.0D0
!!              interp_pixy_x = (9.0D0*(pi_xy(a+1,b,c,1)+pi_xy(a-1,b,c,1))-&
!!                pi_xy(a+3,b,c,1)-pi_xy(a-3,b,c,1))/16.0D0
!!              interp_pixz_x = (9.0D0*(pi_xz(a+1,b,c,1)+pi_xz(a-1,b,c,1))-&
!!                pi_xz(a+3,b,c,1)-pi_xz(a-3,b,c,1))/16.0D0
!!              interp_piyz_x = (9.0D0*(pi_yz(a+1,b,c,1)+pi_yz(a-1,b,c,1))-&
!!                pi_yz(a+3,b,c,1)-pi_yz(a-3,b,c,1))/16.0D0
!!              interp_lambdax_x = (9.0D0*(lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))-&
!!                lambda_x(a+3,b,c,1)-lambda_x(a-3,b,c,1))/16.0D0
!!              interp_lambday_x = (9.0D0*(lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))-&
!!                lambda_y(a+3,b,c,1)-lambda_y(a-3,b,c,1))/16.0D0
!!              interp_lambdaz_x = (9.0D0*(lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))-&
!!                lambda_z(a+3,b,c,1)-lambda_z(a-3,b,c,1))/16.0D0
!            else if (x_up_clear .and. x_down_close .and. .not. x_down_clear) then
!              do i = 1,dir
!                interp_fix(i) = 3.0D0/8.0D0*fi(a+1,b,c,i)+3.0D0/4.0D0*fi(a-1,b,c,i)-fi(a-3,b,c,i)/8.0D0
!                interp_gix(i) = 3.0D0/8.0D0*gi(a+1,b,c,i)+3.0D0/4.0D0*gi(a-1,b,c,i)-gi(a-3,b,c,i)/8.0D0
!              end do
!!              interp_chi_x      = 3.0D0/8.0D0*chi(a+1,b,c,1)     +3.0D0/4.0D0*chi(a-1,b,c,1)     -chi(a-3,b,c,1)/8.0D0
!!              interp_zetax_x   = 3.0D0/8.0D0*zeta_x(a+1,b,c,1)  +3.0D0/4.0D0*zeta_x(a-1,b,c,1)  -zeta_x(a-3,b,c,1)/8.0D0
!!              interp_zetay_x   = 3.0D0/8.0D0*zeta_y(a+1,b,c,1)  +3.0D0/4.0D0*zeta_y(a-1,b,c,1)  -zeta_y(a-3,b,c,1)/8.0D0
!!              interp_zetaz_x   = 3.0D0/8.0D0*zeta_z(a+1,b,c,1)  +3.0D0/4.0D0*zeta_z(a-1,b,c,1)  -zeta_z(a-3,b,c,1)/8.0D0
!!              interp_pixx_x    = 3.0D0/8.0D0*pi_xx(a+1,b,c,1)   +3.0D0/4.0D0*pi_xx(a-1,b,c,1)   -pi_xx(a-3,b,c,1)/8.0D0
!!              interp_piyy_x    = 3.0D0/8.0D0*pi_yy(a+1,b,c,1)   +3.0D0/4.0D0*pi_yy(a-1,b,c,1)   -pi_yy(a-3,b,c,1)/8.0D0
!!              interp_pizz_x    = 3.0D0/8.0D0*pi_zz(a+1,b,c,1)   +3.0D0/4.0D0*pi_zz(a-1,b,c,1)   -pi_zz(a-3,b,c,1)/8.0D0
!!              interp_pixy_x    = 3.0D0/8.0D0*pi_xy(a+1,b,c,1)   +3.0D0/4.0D0*pi_xy(a-1,b,c,1)   -pi_xy(a-3,b,c,1)/8.0D0
!!              interp_pixz_x    = 3.0D0/8.0D0*pi_xz(a+1,b,c,1)   +3.0D0/4.0D0*pi_xz(a-1,b,c,1)   -pi_xz(a-3,b,c,1)/8.0D0
!!              interp_piyz_x    = 3.0D0/8.0D0*pi_yz(a+1,b,c,1)   +3.0D0/4.0D0*pi_yz(a-1,b,c,1)   -pi_yz(a-3,b,c,1)/8.0D0
!!              interp_lambdax_x = 3.0D0/8.0D0*lambda_x(a+1,b,c,1)+3.0D0/4.0D0*lambda_x(a-1,b,c,1)-lambda_x(a-3,b,c,1)/8.0D0
!!              interp_lambday_x = 3.0D0/8.0D0*lambda_y(a+1,b,c,1)+3.0D0/4.0D0*lambda_y(a-1,b,c,1)-lambda_y(a-3,b,c,1)/8.0D0
!!              interp_lambdaz_x = 3.0D0/8.0D0*lambda_z(a+1,b,c,1)+3.0D0/4.0D0*lambda_z(a-1,b,c,1)-lambda_z(a-3,b,c,1)/8.0D0
!            else if (x_down_clear .and. x_up_close .and. .not. x_up_clear) then
!              do i = 1,dir
!                interp_fix(i) = 3.0D0/8.0D0*fi(a-1,b,c,i)+3.0D0/4.0D0*fi(a+1,b,c,i)-fi(a+3,b,c,i)/8.0D0
!                interp_gix(i) = 3.0D0/8.0D0*gi(a-1,b,c,i)+3.0D0/4.0D0*gi(a+1,b,c,i)-gi(a+3,b,c,i)/8.0D0
!              end do
!!              interp_chi_x      = 3.0D0/8.0D0*chi(a-1,b,c,1)     +3.0D0/4.0D0*chi(a+1,b,c,1)     -chi(a+3,b,c,1)/8.0D0
!!              interp_zetax_x   = 3.0D0/8.0D0*zeta_x(a-1,b,c,1)  +3.0D0/4.0D0*zeta_x(a+1,b,c,1)  -zeta_x(a+3,b,c,1)/8.0D0
!!              interp_zetay_x   = 3.0D0/8.0D0*zeta_y(a-1,b,c,1)  +3.0D0/4.0D0*zeta_y(a+1,b,c,1)  -zeta_y(a+3,b,c,1)/8.0D0
!!              interp_zetaz_x   = 3.0D0/8.0D0*zeta_z(a-1,b,c,1)  +3.0D0/4.0D0*zeta_z(a+1,b,c,1)  -zeta_z(a+3,b,c,1)/8.0D0
!!              interp_pixx_x    = 3.0D0/8.0D0*pi_xx(a-1,b,c,1)   +3.0D0/4.0D0*pi_xx(a+1,b,c,1)   -pi_xx(a+3,b,c,1)/8.0D0
!!              interp_piyy_x    = 3.0D0/8.0D0*pi_yy(a-1,b,c,1)   +3.0D0/4.0D0*pi_yy(a+1,b,c,1)   -pi_yy(a+3,b,c,1)/8.0D0
!!              interp_pizz_x    = 3.0D0/8.0D0*pi_zz(a-1,b,c,1)   +3.0D0/4.0D0*pi_zz(a+1,b,c,1)   -pi_zz(a+3,b,c,1)/8.0D0
!!              interp_pixy_x    = 3.0D0/8.0D0*pi_xy(a-1,b,c,1)   +3.0D0/4.0D0*pi_xy(a+1,b,c,1)   -pi_xy(a+3,b,c,1)/8.0D0
!!              interp_pixz_x    = 3.0D0/8.0D0*pi_xz(a-1,b,c,1)   +3.0D0/4.0D0*pi_xz(a+1,b,c,1)   -pi_xz(a+3,b,c,1)/8.0D0
!!              interp_piyz_x    = 3.0D0/8.0D0*pi_yz(a-1,b,c,1)   +3.0D0/4.0D0*pi_yz(a+1,b,c,1)   -pi_yz(a+3,b,c,1)/8.0D0
!!              interp_lambdax_x = 3.0D0/8.0D0*lambda_x(a-1,b,c,1)+3.0D0/4.0D0*lambda_x(a+1,b,c,1)-lambda_x(a+3,b,c,1)/8.0D0
!!              interp_lambday_x = 3.0D0/8.0D0*lambda_y(a-1,b,c,1)+3.0D0/4.0D0*lambda_y(a+1,b,c,1)-lambda_y(a+3,b,c,1)/8.0D0
!!              interp_lambdaz_x = 3.0D0/8.0D0*lambda_z(a-1,b,c,1)+3.0D0/4.0D0*lambda_z(a+1,b,c,1)-lambda_z(a+3,b,c,1)/8.0D0
!            else if (x_up_close .and. x_down_close) then
!              do i=1,dir
!                interp_fix(i) = (fi(a+1,b,c,i)+fi(a-1,b,c,i))/2.0D0
!                interp_gix(i) = (gi(a+1,b,c,i)+gi(a-1,b,c,i))/2.0D0
!              end do
!!              interp_chi_x    =        (chi(a+1,b,c,1)+     chi(a-1,b,c,1))/2.0D0
!!              interp_zetax_x =     (zeta_x(a+1,b,c,1)+  zeta_x(a-1,b,c,1))/2.0D0
!!              interp_zetay_x =     (zeta_y(a+1,b,c,1)+  zeta_y(a-1,b,c,1))/2.0D0
!!              interp_zetaz_x =     (zeta_z(a+1,b,c,1)+  zeta_z(a-1,b,c,1))/2.0D0
!!              interp_pixx_x =       (pi_xx(a+1,b,c,1)+   pi_xx(a-1,b,c,1))/2.0D0
!!              interp_piyy_x =       (pi_yy(a+1,b,c,1)+   pi_yy(a-1,b,c,1))/2.0D0
!!              interp_pizz_x =       (pi_zz(a+1,b,c,1)+   pi_zz(a-1,b,c,1))/2.0D0
!!              interp_pixy_x =       (pi_xy(a+1,b,c,1)+   pi_xy(a-1,b,c,1))/2.0D0
!!              interp_pixz_x =       (pi_xz(a+1,b,c,1)+   pi_xz(a-1,b,c,1))/2.0D0
!!              interp_piyz_x =       (pi_yz(a+1,b,c,1)+   pi_yz(a-1,b,c,1))/2.0D0
!!              interp_lambdax_x = (lambda_x(a+1,b,c,1)+lambda_x(a-1,b,c,1))/2.0D0
!!              interp_lambday_x = (lambda_y(a+1,b,c,1)+lambda_y(a-1,b,c,1))/2.0D0
!!              interp_lambdaz_x = (lambda_z(a+1,b,c,1)+lambda_z(a-1,b,c,1))/2.0D0
!            else
!              do i =1,dir
!                interp_fix(i) = crs_fi(a,b,c,i)
!                interp_gix(i) = crs_gi(a,b,c,i)
!              end do
!            end if
!!
!            if (y_up_clear .and. y_down_clear) then
!              do i = 1,dir
!                interp_fiy(i) = (9.0D0*(fi(a,b+1,c,i)+fi(a,b-1,c,i))-&
!                                 fi(a,b-3,c,i)-fi(a,b+3,c,i))/16.0D0
!                interp_giy(i) = (9.0D0*(gi(a,b+1,c,i)+gi(a,b-1,c,i))-&
!                                 gi(a,b+3,c,i)-gi(a,b-3,c,i))/16.0D0
!              end do
!!              interp_chi_y = (9.0D0*(chi(a,b+1,c,1)+chi(a,b-1,c,1))-&
!!                chi(a,b+3,c,1)-chi(a,b-3,c,1))/16.0D0
!!              interp_zetax_y = (9.0D0*(zeta_x(a,b+1,c,1)+zeta_x(a,b-1,c,1))-&
!!                zeta_x(a,b+3,c,1)-zeta_x(a,b-3,c,1))/16.0D0
!!              interp_zetay_y = (9.0D0*(zeta_y(a,b+1,c,1)+zeta_y(a,b-1,c,1))-&
!!                zeta_y(a,b+3,c,1)-zeta_y(a,b-3,c,1))/16.0D0
!!              interp_zetaz_y = (9.0D0*(zeta_z(a,b+1,c,1)+zeta_z(a,b-1,c,1))-&
!!                zeta_z(a,b+3,c,1)-zeta_z(a,b-3,c,1))/16.0D0
!!              interp_pixx_y = (9.0D0*(pi_xx(a,b+1,c,1)+pi_xx(a,b-1,c,1))-&
!!                pi_xx(a,b+3,c,1)-pi_xx(a,b-3,c,1))/16.0D0
!!              interp_piyy_y = (9.0D0*(pi_yy(a,b+1,c,1)+pi_yy(a,b-1,c,1))-&
!!                pi_yy(a,b-3,c,1)-pi_yy(a,b+3,c,1))/16.0D0
!!              interp_pizz_y = (9.0D0*(pi_zz(a,b+1,c,1)+pi_zz(a,b-1,c,1))-&
!!                pi_zz(a,b+3,c,1)-pi_zz(a,b-3,c,1))/16.0D0
!!              interp_pixy_y = (9.0D0*(pi_xy(a,b+1,c,1)+pi_xy(a,b-1,c,1))-&
!!                pi_xy(a,b+3,c,1)-pi_xy(a,b-3,c,1))/16.0D0
!!              interp_pixz_y = (9.0D0*(pi_xz(a,b+1,c,1)+pi_xz(a,b-1,c,1))-&
!!                pi_xz(a,b+3,c,1)-pi_xz(a,b-3,c,1))/16.0D0
!!              interp_piyz_y = (9.0D0*(pi_yz(a,b+1,c,1)+pi_yz(a,b-1,c,1))-&
!!                pi_yz(a,b+3,c,1)-pi_yz(a,b-3,c,1))/16.0D0
!!              interp_lambdax_y = (9.0D0*(lambda_x(a,b+1,c,1)+lambda_x(a,b-1,c,1))-&
!!                lambda_x(a,b+3,c,1)-lambda_x(a,b-3,c,1))/16.0D0
!!              interp_lambday_y = (9.0D0*(lambda_y(a,b+1,c,1)+lambda_y(a,b-1,c,1))-&
!!                lambda_y(a,b+3,c,1)-lambda_y(a,b-3,c,1))/16.0D0
!!              interp_lambdaz_y = (9.0D0*(lambda_z(a,b+1,c,1)+lambda_z(a,b-1,c,1))-&
!!                lambda_z(a,b+3,c,1)-lambda_z(a,b-3,c,1))/16.0D0
!            else if (y_up_clear .and. y_down_close .and. .not. y_down_clear) then
!              do i = 1,dir
!                interp_fiy(i) = 3.0D0/8.0D0*fi(a,b+1,c,i)+3.0D0/4.0D0*fi(a,b-1,c,i)-fi(a,b-3,c,i)/8.0D0
!                interp_giy(i) = 3.0D0/8.0D0*gi(a,b+1,c,i)+3.0D0/4.0D0*gi(a,b-1,c,i)-gi(a,b-3,c,i)/8.0D0
!              end do
!!              interp_chi_y = 3.0D0/8.0D0*chi(a,b+1,c,1)     +3.0D0/4.0D0*chi(a,b-1,c,1)     -chi(a,b-3,c,1)/8.0D0
!!              interp_zetax_y = 3.0D0/8.0D0*zeta_x(a,b+1,c,1)  +3.0D0/4.0D0*zeta_x(a,b-1,c,1)  -zeta_x(a,b-3,c,1)/8.0D0
!!              interp_zetay_y = 3.0D0/8.0D0*zeta_y(a,b+1,c,1)  +3.0D0/4.0D0*zeta_y(a,b-1,c,1)  -zeta_y(a,b-3,c,1)/8.0D0
!!              interp_zetaz_y = 3.0D0/8.0D0*zeta_z(a,b+1,c,1)  +3.0D0/4.0D0*zeta_z(a,b-1,c,1)  -zeta_z(a,b-3,c,1)/8.0D0
!!              interp_pixx_y = 3.0D0/8.0D0*pi_xx(a,b+1,c,1)   +3.0D0/4.0D0*pi_xx(a,b-1,c,1)   -pi_xx(a,b-3,c,1)/8.0D0
!!              interp_piyy_y = 3.0D0/8.0D0*pi_yy(a,b+1,c,1)   +3.0D0/4.0D0*pi_yy(a,b-1,c,1)   -pi_yy(a,b-3,c,1)/8.0D0
!!              interp_pizz_y = 3.0D0/8.0D0*pi_zz(a,b+1,c,1)   +3.0D0/4.0D0*pi_zz(a,b-1,c,1)   -pi_zz(a,b-3,c,1)/8.0D0
!!              interp_pixy_y = 3.0D0/8.0D0*pi_xy(a,b+1,c,1)   +3.0D0/4.0D0*pi_xy(a,b-1,c,1)   -pi_xy(a,b-3,c,1)/8.0D0
!!              interp_pixz_y = 3.0D0/8.0D0*pi_xz(a,b+1,c,1)   +3.0D0/4.0D0*pi_xz(a,b-1,c,1)   -pi_xz(a,b-3,c,1)/8.0D0
!!              interp_piyz_y = 3.0D0/8.0D0*pi_yz(a,b+1,c,1)   +3.0D0/4.0D0*pi_yz(a,b-1,c,1)   -pi_yz(a,b-3,c,1)/8.0D0
!!              interp_lambdax_y = 3.0D0/8.0D0*lambda_x(a,b+1,c,1)+3.0D0/4.0D0*lambda_x(a,b-1,c,1)-lambda_x(a,b-3,c,1)/8.0D0
!!              interp_lambday_y = 3.0D0/8.0D0*lambda_y(a,b+1,c,1)+3.0D0/4.0D0*lambda_y(a,b-1,c,1)-lambda_y(a,b-3,c,1)/8.0D0
!!              interp_lambdaz_y = 3.0D0/8.0D0*lambda_z(a,b+1,c,1)+3.0D0/4.0D0*lambda_z(a,b-1,c,1)-lambda_z(a,b-3,c,1)/8.0D0
!            else if (y_down_clear .and. y_up_close .and. .not. y_up_clear) then
!              do i = 1,dir
!                interp_fiy(i) = 3.0D0/8.0D0*fi(a,b-1,c,i)+3.0D0/4.0D0*fi(a,b+1,c,i)-fi(a,b+3,c,i)/8.0D0
!                interp_giy(i) = 3.0D0/8.0D0*gi(a,b-1,c,i)+3.0D0/4.0D0*gi(a,b+1,c,i)-gi(a,b+3,c,i)/8.0D0
!              end do
!!              interp_chi_y = 3.0D0/8.0D0*chi(a,b-1,c,1)     +3.0D0/4.0D0*chi(a,b+1,c,1)     -chi(a,b+3,c,1)/8.0D0
!!              interp_zetax_y = 3.0D0/8.0D0*zeta_x(a,b-1,c,1)  +3.0D0/4.0D0*zeta_x(a,b+1,c,1)  -zeta_x(a,b+3,c,1)/8.0D0
!!              interp_zetay_y = 3.0D0/8.0D0*zeta_y(a,b-1,c,1)  +3.0D0/4.0D0*zeta_y(a,b+1,c,1)  -zeta_y(a,b+3,c,1)/8.0D0
!!              interp_zetaz_y = 3.0D0/8.0D0*zeta_z(a,b-1,c,1)  +3.0D0/4.0D0*zeta_z(a,b+1,c,1)  -zeta_z(a,b+3,c,1)/8.0D0
!!              interp_pixx_y = 3.0D0/8.0D0*pi_xx(a,b-1,c,1)   +3.0D0/4.0D0*pi_xx(a,b+1,c,1)   -pi_xx(a,b+3,c,1)/8.0D0
!!              interp_piyy_y = 3.0D0/8.0D0*pi_yy(a,b-1,c,1)   +3.0D0/4.0D0*pi_yy(a,b+1,c,1)   -pi_yy(a,b+3,c,1)/8.0D0
!!              interp_pizz_y = 3.0D0/8.0D0*pi_zz(a,b-1,c,1)   +3.0D0/4.0D0*pi_zz(a,b+1,c,1)   -pi_zz(a,b+3,c,1)/8.0D0
!!              interp_pixy_y = 3.0D0/8.0D0*pi_xy(a,b-1,c,1)   +3.0D0/4.0D0*pi_xy(a,b+1,c,1)   -pi_xy(a,b+3,c,1)/8.0D0
!!              interp_pixz_y = 3.0D0/8.0D0*pi_xz(a,b-1,c,1)   +3.0D0/4.0D0*pi_xz(a,b+1,c,1)   -pi_xz(a,b+3,c,1)/8.0D0
!!              interp_piyz_y = 3.0D0/8.0D0*pi_yz(a,b-1,c,1)   +3.0D0/4.0D0*pi_yz(a,b+1,c,1)   -pi_yz(a,b+3,c,1)/8.0D0
!!              interp_lambdax_y = 3.0D0/8.0D0*lambda_x(a,b-1,c,1)+3.0D0/4.0D0*lambda_x(a,b+1,c,1)-lambda_x(a,b+3,c,1)/8.0D0
!!              interp_lambday_y = 3.0D0/8.0D0*lambda_y(a,b-1,c,1)+3.0D0/4.0D0*lambda_y(a,b+1,c,1)-lambda_y(a,b+3,c,1)/8.0D0
!!              interp_lambdaz_y = 3.0D0/8.0D0*lambda_z(a,b-1,c,1)+3.0D0/4.0D0*lambda_z(a,b+1,c,1)-lambda_z(a,b+3,c,1)/8.0D0
!            else if (y_up_close .and. y_down_close) then
!              do i=1,dir
!                interp_fiy(i) = (fi(a,b+1,c,i)+fi(a,b-1,c,i))/2.0D0
!                interp_giy(i) = (gi(a,b+1,c,i)+gi(a,b-1,c,i))/2.0D0
!              end do
!!              interp_chi_y    =       (chi(a,b+1,c,1)+     chi(a,b-1,c,1))/2.0D0
!!              interp_zetax_y =     (zeta_x(a,b+1,c,1)+  zeta_x(a,b-1,c,1))/2.0D0
!!              interp_zetay_y =     (zeta_y(a,b+1,c,1)+  zeta_y(a,b-1,c,1))/2.0D0
!!              interp_zetaz_y =     (zeta_z(a,b+1,c,1)+  zeta_z(a,b-1,c,1))/2.0D0
!!              interp_pixx_y =       (pi_xx(a,b+1,c,1)+   pi_xx(a,b-1,c,1))/2.0D0
!!              interp_piyy_y =       (pi_yy(a,b+1,c,1)+   pi_yy(a,b-1,c,1))/2.0D0
!!              interp_pizz_y =       (pi_zz(a,b+1,c,1)+   pi_zz(a,b-1,c,1))/2.0D0
!!              interp_pixy_y =       (pi_xy(a,b+1,c,1)+   pi_xy(a,b-1,c,1))/2.0D0
!!              interp_pixz_y =       (pi_xz(a,b+1,c,1)+   pi_xz(a,b-1,c,1))/2.0D0
!!              interp_piyz_y =       (pi_yz(a,b+1,c,1)+   pi_yz(a,b-1,c,1))/2.0D0
!!              interp_lambdax_y = (lambda_x(a,b+1,c,1)+lambda_x(a,b-1,c,1))/2.0D0
!!              interp_lambday_y = (lambda_y(a,b+1,c,1)+lambda_y(a,b-1,c,1))/2.0D0
!!              interp_lambdaz_y = (lambda_z(a,b+1,c,1)+lambda_z(a,b-1,c,1))/2.0D0
!            else
!              do i =1,dir
!                interp_fiy(i) = crs_fi(a,b,c,i)
!                interp_giy(i) = crs_gi(a,b,c,i)
!              end do
!            end if
!!
!            if (z_up_clear .and. z_down_clear ) then
!              do i = 1,dir
!                interp_fiz(i) = (9.0D0*(fi(a,b,c+1,i)+fi(a,b,c-1,i))-&
!                                 fi(a,b,c-3,i)-fi(a,b,c+3,i))/16.0D0
!
!                interp_giz(i) = (9.0D0*(gi(a,b,c+1,i)+gi(a,b,c-1,i))-&
!                                 gi(a,b,c+3,i)-gi(a,b,c-3,i))/16.0D0
!              end do
!!              interp_chi_z = (9.0D0*(chi(a,b,c+1,1)+chi(a,b,c-1,1))-&
!!                chi(a,b,c+3,1)-chi(a,b,c-3,1))/16.0D0
!!              interp_zetax_z = (9.0D0*(zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))-&
!!                zeta_x(a,b,c+3,1)-zeta_x(a,b,c-3,1))/16.0D0
!!              interp_zetay_z = (9.0D0*(zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))-&
!!                zeta_y(a,b,c+3,1)-zeta_y(a,b,c-3,1))/16.0D0
!!              interp_zetaz_z = (9.0D0*(zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))-&
!!                zeta_z(a,b,c+3,1)-zeta_z(a,b,c-3,1))/16.0D0
!!              interp_pixx_z = (9.0D0*(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))-&
!!                pi_xx(a,b,c+3,1)-pi_xx(a,b,c-3,1))/16.0D0
!!              interp_piyy_z = (9.0D0*(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))-&
!!                pi_yy(a,b,c+3,1)-pi_yy(a,b,c-3,1))/16.0D0
!!              interp_pizz_z = (9.0D0*(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))-&
!!                pi_zz(a,b,c+3,1)-pi_zz(a,b,c-3,1))/16.0D0
!!              interp_pixy_z = (9.0D0*(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))-&
!!                pi_xy(a,b,c+3,1)-pi_xy(a,b,c-3,1))/16.0D0
!!              interp_pixz_z = (9.0D0*(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))-&
!!                pi_xz(a,b,c+3,1)-pi_xz(a,b,c-3,1))/16.0D0
!!              interp_piyz_z = (9.0D0*(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))-&
!!                pi_yz(a,b,c+3,1)-pi_yz(a,b,c-3,1))/16.0D0
!!              interp_lambdax_z = (9.0D0*(lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))-&
!!                lambda_x(a,b,c+3,1)-lambda_x(a,b,c-3,1))/16.0D0
!!              interp_lambday_z = (9.0D0*(lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))-&
!!                lambda_y(a,b,c+3,1)-lambda_y(a,b,c-3,1))/16.0D0
!!              interp_lambdaz_z = (9.0D0*(lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))-&
!!                lambda_z(a,b,c+3,1)-lambda_z(a,b,c-3,1))/16.0D0
!            else if (z_up_clear .and. z_down_close .and. .not. z_down_clear) then
!              do i = 1,dir
!                interp_fiz(i) = 3.0D0/8.0D0*fi(a,b,c+1,i)+3.0D0/4.0D0*fi(a,b,c-1,i)-fi(a,b,c-3,i)/8.0D0
!                interp_giz(i) = 3.0D0/8.0D0*gi(a,b,c+1,i)+3.0D0/4.0D0*gi(a,b,c-1,i)-gi(a,b,c-3,i)/8.0D0
!              end do
!!              interp_chi_z =     3.0D0/8.0D0*chi(a,b,c+1,1)+3.0D0/4.0D0*chi(a,b,c-1,1)-chi(a,b,c-3,1)/8.0D0
!!              interp_zetax_z =   3.0D0/8.0D0*zeta_x(a,b,c+1,1)+3.0D0/4.0D0*zeta_x(a,b,c-1,1)-zeta_x(a,b,c-3,1)/8.0D0
!!              interp_zetay_z =   3.0D0/8.0D0*zeta_y(a,b,c+1,1)+3.0D0/4.0D0*zeta_y(a,b,c-1,1)-zeta_y(a,b,c-3,1)/8.0D0
!!              interp_zetaz_z =   3.0D0/8.0D0*zeta_z(a,b,c+1,1)+3.0D0/4.0D0*zeta_z(a,b,c-1,1)-zeta_z(a,b,c-3,1)/8.0D0
!!              interp_pixx_z =    3.0D0/8.0D0*pi_xx(a,b,c+1,1)+3.0D0/4.0D0*pi_xx(a,b,c-1,1)-pi_xx(a,b,c-3,1)/8.0D0
!!              interp_piyy_z =    3.0D0/8.0D0*pi_yy(a,b,c+1,1)+3.0D0/4.0D0*pi_yy(a,b,c-1,1)-pi_yy(a,b,c-3,1)/8.0D0
!!              interp_pizz_z =    3.0D0/8.0D0*pi_zz(a,b,c+1,1)+3.0D0/4.0D0*pi_zz(a,b,c-1,1)-pi_zz(a,b,c-3,1)/8.0D0
!!              interp_pixy_z =    3.0D0/8.0D0*pi_xy(a,b,c+1,1)+3.0D0/4.0D0*pi_xy(a,b,c-1,1)-pi_xy(a,b,c-3,1)/8.0D0
!!              interp_pixz_z =    3.0D0/8.0D0*pi_xz(a,b,c+1,1)+3.0D0/4.0D0*pi_xz(a,b,c-1,1)-pi_xz(a,b,c-3,1)/8.0D0
!!              interp_piyz_z =    3.0D0/8.0D0*pi_yz(a,b,c+1,1)+3.0D0/4.0D0*pi_yz(a,b,c-1,1)-pi_yz(a,b,c-3,1)/8.0D0
!!              interp_lambdax_z = 3.0D0/8.0D0*lambda_x(a,b,c+1,1)+3.0D0/4.0D0*lambda_x(a,b,c-1,1)-lambda_x(a,b,c-3,1)/8.0D0
!!              interp_lambday_z = 3.0D0/8.0D0*lambda_y(a,b,c+1,1)+3.0D0/4.0D0*lambda_y(a,b,c-1,1)-lambda_y(a,b,c-3,1)/8.0D0
!!              interp_lambdaz_z = 3.0D0/8.0D0*lambda_z(a,b,c+1,1)+3.0D0/4.0D0*lambda_z(a,b,c-1,1)-lambda_z(a,b,c-3,1)/8.0D0
!            else if (z_down_clear .and. z_up_close .and. .not. z_up_clear) then
!              do i = 1,dir
!                interp_fiz(i) = 3.0D0/8.0D0*fi(a,b,c-1,i)+3.0D0/4.0D0*fi(a,b,c+1,i)-fi(a,b,c+3,i)/8.0D0
!                interp_giz(i) = 3.0D0/8.0D0*gi(a,b,c-1,i)+3.0D0/4.0D0*gi(a,b,c+1,i)-gi(a,b,c+3,i)/8.0D0
!              end do
!!              interp_chi_z =      3.0D0/8.0D0*chi(a,b,c-1,1)+3.0D0/4.0D0*chi(a,b,c+1,1)-chi(a,b,c+3,1)/8.0D0
!!              interp_zetax_z =   3.0D0/8.0D0*zeta_x(a,b,c-1,1)+3.0D0/4.0D0*zeta_x(a,b,c+1,1)-zeta_x(a,b,c+3,1)/8.0D0
!!              interp_zetay_z =   3.0D0/8.0D0*zeta_y(a,b,c-1,1)+3.0D0/4.0D0*zeta_y(a,b,c+1,1)-zeta_y(a,b,c+3,1)/8.0D0
!!              interp_zetaz_z =   3.0D0/8.0D0*zeta_z(a,b,c-1,1)+3.0D0/4.0D0*zeta_z(a,b,c+1,1)-zeta_z(a,b,c+3,1)/8.0D0
!!              interp_pixx_z =    3.0D0/8.0D0*pi_xx(a,b,c-1,1)+3.0D0/4.0D0*pi_xx(a,b,c+1,1)-pi_xx(a,b,c+3,1)/8.0D0
!!              interp_piyy_z =    3.0D0/8.0D0*pi_yy(a,b,c-1,1)+3.0D0/4.0D0*pi_yy(a,b,c+1,1)-pi_yy(a,b,c+3,1)/8.0D0
!!              interp_pizz_z =    3.0D0/8.0D0*pi_zz(a,b,c-1,1)+3.0D0/4.0D0*pi_zz(a,b,c+1,1)-pi_zz(a,b,c+3,1)/8.0D0
!!              interp_pixy_z =    3.0D0/8.0D0*pi_xy(a,b,c-1,1)+3.0D0/4.0D0*pi_xy(a,b,c+1,1)-pi_xy(a,b,c+3,1)/8.0D0
!!              interp_pixz_z =    3.0D0/8.0D0*pi_xz(a,b,c-1,1)+3.0D0/4.0D0*pi_xz(a,b,c+1,1)-pi_xz(a,b,c+3,1)/8.0D0
!!              interp_piyz_z =    3.0D0/8.0D0*pi_yz(a,b,c-1,1)+3.0D0/4.0D0*pi_yz(a,b,c+1,1)-pi_yz(a,b,c+3,1)/8.0D0
!!              interp_lambdax_z= 3.0D0/8.0D0*lambda_x(a,b,c-1,1)+3.0D0/4.0D0*lambda_x(a,b,c+1,1)-lambda_x(a,b,c+3,1)/8.0D0
!!              interp_lambday_z= 3.0D0/8.0D0*lambda_y(a,b,c-1,1)+3.0D0/4.0D0*lambda_y(a,b,c+1,1)-lambda_y(a,b,c+3,1)/8.0D0
!!              interp_lambdaz_z= 3.0D0/8.0D0*lambda_z(a,b,c-1,1)+3.0D0/4.0D0*lambda_z(a,b,c+1,1)-lambda_z(a,b,c+3,1)/8.0D0
!            else if (z_up_close .and. z_down_close) then
!              do i=1,dir
!                interp_fiz(i) = (fi(a,b,c+1,i)+fi(a,b,c-1,i))/2.0D0
!                interp_giz(i) = (gi(a,b,c+1,i)+gi(a,b,c-1,i))/2.0D0
!              end do
!!              interp_chi_z = (chi(a,b,c+1,1)+chi(a,b,c-1,1))/2.0D0
!!              interp_zetax_z = (zeta_x(a,b,c+1,1)+zeta_x(a,b,c-1,1))/2.0D0
!!              interp_zetay_z = (zeta_y(a,b,c+1,1)+zeta_y(a,b,c-1,1))/2.0D0
!!              interp_zetaz_z = (zeta_z(a,b,c+1,1)+zeta_z(a,b,c-1,1))/2.0D0
!!              interp_pixx_z =(pi_xx(a,b,c+1,1)+pi_xx(a,b,c-1,1))/2.0D0
!!              interp_piyy_z =(pi_yy(a,b,c+1,1)+pi_yy(a,b,c-1,1))/2.0D0
!!              interp_pizz_z =(pi_zz(a,b,c+1,1)+pi_zz(a,b,c-1,1))/2.0D0
!!              interp_pixy_z =(pi_xy(a,b,c+1,1)+pi_xy(a,b,c-1,1))/2.0D0
!!              interp_pixz_z =(pi_xz(a,b,c+1,1)+pi_xz(a,b,c-1,1))/2.0D0
!!              interp_piyz_z =(pi_yz(a,b,c+1,1)+pi_yz(a,b,c-1,1))/2.0D0
!!              interp_lambdax_z = (lambda_x(a,b,c+1,1)+lambda_x(a,b,c-1,1))/2.0D0
!!              interp_lambday_z = (lambda_y(a,b,c+1,1)+lambda_y(a,b,c-1,1))/2.0D0
!!              interp_lambdaz_z = (lambda_z(a,b,c+1,1)+lambda_z(a,b,c-1,1))/2.0D0
!            else
!              do i =1,dir
!                interp_fiz(i) = crs_fi(a,b,c,i)
!                interp_giz(i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!            if (.not. x_null .and. .not. y_null .and. .not. z_null) then
!
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fix(i) + interp_fiy(i) + interp_fiz(i))/3.0D0
!                gi(a,b,c,i) = (interp_gix(i) + interp_giy(i) + interp_giz(i))/3.0D0
!!                if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!!                  write(*,*) 'big badda boomxyz',fi(a,b,c,i),gi(a,b,c,i),a,b,c,i
!!                end if
!              end do
!!              chi(a,b,c,1)    = (interp_chi_z + interp_chi_y + interp_chi_x)/3.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_z + interp_zetax_y + interp_zetax_x)/3.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_z + interp_zetay_y + interp_zetay_x)/3.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_z + interp_zetaz_y + interp_zetaz_x)/3.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_z + interp_pixx_y + interp_pixx_x)/3.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_z + interp_piyy_y + interp_piyy_x)/3.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_z + interp_pizz_y + interp_pizz_x)/3.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_z + interp_pixy_y + interp_pixy_x)/3.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_z + interp_pixz_y + interp_pixz_x)/3.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_z + interp_piyz_y + interp_piyz_x)/3.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_z + interp_lambdax_y + interp_lambdax_x)/3.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_z + interp_lambday_y + interp_lambday_x)/3.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_z + interp_lambdaz_y + interp_lambdaz_x)/3.0D0
!            else if (x_null .and. .not. y_null .and. .not. z_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fiy(i)+interp_fiz(i))/2
!                gi(a,b,c,i) = (interp_giy(i)+interp_giz(i))/2
!              end do
!!              chi(a,b,c,1)    = (interp_chi_z + interp_chi_y)/2.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_z + interp_zetax_y)/2.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_z + interp_zetay_y)/2.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_z + interp_zetaz_y)/2.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_z + interp_pixx_y)/2.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_z + interp_piyy_y)/2.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_z + interp_pizz_y)/2.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_z + interp_pixy_y)/2.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_z + interp_pixz_y)/2.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_z + interp_piyz_y)/2.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_z + interp_lambdax_y)/2.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_z + interp_lambday_y)/2.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_z + interp_lambdaz_y)/2.0D0
!            else if (y_null .and. .not. x_null .and. .not. z_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fix(i)+interp_fiz(i))/2
!                gi(a,b,c,i) = (interp_gix(i)+interp_giz(i))/2
!              end do
!!              chi(a,b,c,1)    = (interp_chi_z + interp_chi_x)/2.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_z + interp_zetax_x)/2.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_z + interp_zetay_x)/2.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_z + interp_zetaz_x)/2.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_z + interp_pixx_x)/2.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_z + interp_piyy_x)/2.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_z + interp_pizz_x)/2.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_z + interp_pixy_x)/2.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_z + interp_pixz_x)/2.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_z + interp_piyz_x)/2.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_z + interp_lambdax_x)/2.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_z + interp_lambday_x)/2.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_z + interp_lambdaz_x)/2.0D0
!            else if (z_null .and. .not. y_null .and. .not. x_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = (interp_fiy(i)+interp_fix(i))/2
!                gi(a,b,c,i) = (interp_giy(i)+interp_gix(i))/2
!              end do
!!              chi(a,b,c,1)    = (interp_chi_y + interp_chi_x)/2.0D0
!!              zeta_x(a,b,c,1) = (interp_zetax_y + interp_zetax_x)/2.0D0
!!              zeta_y(a,b,c,1) = (interp_zetay_y + interp_zetay_x)/2.0D0
!!              zeta_z(a,b,c,1) = (interp_zetaz_y + interp_zetaz_x)/2.0D0
!!              pi_xx(a,b,c,1) =  (interp_pixx_y + interp_pixx_x)/2.0D0
!!              pi_yy(a,b,c,1) =  (interp_piyy_y + interp_piyy_x)/2.0D0
!!              pi_zz(a,b,c,1) =  (interp_pizz_y + interp_pizz_x)/2.0D0
!!              pi_xy(a,b,c,1) =  (interp_pixy_y + interp_pixy_x)/2.0D0
!!              pi_xz(a,b,c,1) =  (interp_pixz_y + interp_pixz_x)/2.0D0
!!              pi_yz(a,b,c,1) =  (interp_piyz_y + interp_piyz_x)/2.0D0
!!              lambda_x(a,b,c,1) = (interp_lambdax_y + interp_lambdax_x)/2.0D0
!!              lambda_y(a,b,c,1) = (interp_lambday_y + interp_lambday_x)/2.0D0
!!              lambda_z(a,b,c,1) = (interp_lambdaz_y + interp_lambdaz_x)/2.0D0
!            else if (x_null .and. y_null .and. .not. z_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fiz(i)
!                gi(a,b,c,i) = interp_giz(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_z
!!              zeta_x(a,b,c,1) =  interp_zetax_z
!!              zeta_y(a,b,c,1) =  interp_zetay_z
!!              zeta_z(a,b,c,1) =  interp_zetaz_z
!!              pi_xx(a,b,c,1) =  interp_pixx_z
!!              pi_yy(a,b,c,1) =  interp_piyy_z
!!              pi_zz(a,b,c,1) =  interp_pizz_z
!!              pi_xy(a,b,c,1) =  interp_pixy_z
!!              pi_xz(a,b,c,1) =  interp_pixz_z
!!              pi_yz(a,b,c,1) =  interp_piyz_z
!!              lambda_x(a,b,c,1) =interp_lambdax_z
!!              lambda_y(a,b,c,1) =interp_lambday_z
!!              lambda_z(a,b,c,1) =interp_lambdaz_z
!            else if (x_null .and. z_null .and. .not. y_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fiy(i)
!                gi(a,b,c,i) = interp_giy(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_y
!!              zeta_x(a,b,c,1) =  interp_zetax_y
!!              zeta_y(a,b,c,1) =  interp_zetay_y
!!              zeta_z(a,b,c,1) =  interp_zetaz_y
!!              pi_xx(a,b,c,1) =  interp_pixx_y
!!              pi_yy(a,b,c,1) =  interp_piyy_y
!!              pi_zz(a,b,c,1) =  interp_pizz_y
!!              pi_xy(a,b,c,1) =  interp_pixy_y
!!              pi_xz(a,b,c,1) =  interp_pixz_y
!!              pi_yz(a,b,c,1) =  interp_piyz_y
!!              lambda_x(a,b,c,1) =interp_lambdax_y
!!              lambda_y(a,b,c,1) =interp_lambday_y
!!              lambda_z(a,b,c,1) =interp_lambdaz_y
!            else if (y_null .and. z_null .and. .not. x_null) then
!              do i = 1,dir
!                fi(a,b,c,i) = interp_fix(i)
!                gi(a,b,c,i) = interp_gix(i)
!              end do
!!              chi(a,b,c,1)    =  interp_chi_x
!!              zeta_x(a,b,c,1) =  interp_zetax_x
!!              zeta_y(a,b,c,1) =  interp_zetay_x
!!              zeta_z(a,b,c,1) =  interp_zetaz_x
!!              pi_xx(a,b,c,1) =  interp_pixx_x
!!              pi_yy(a,b,c,1) =  interp_piyy_x
!!              pi_zz(a,b,c,1) =  interp_pizz_x
!!              pi_xy(a,b,c,1) =  interp_pixy_x
!!              pi_xz(a,b,c,1) =  interp_pixz_x
!!              pi_yz(a,b,c,1) =  interp_piyz_x
!!              lambda_x(a,b,c,1) =interp_lambdax_x
!!              lambda_y(a,b,c,1) =interp_lambday_x
!!              lambda_z(a,b,c,1) =interp_lambdaz_x
!            else
!              needs_interp(a,b,c) = .true.
!              do i = 1,dir
!                fi(a,b,c,i) = crs_fi(a,b,c,i)
!                gi(a,b,c,i) = crs_gi(a,b,c,i)
!              end do
!            end if
!
!            !write(*,*) 'values out at centers',fi(a,b,c,1),gi(a,b,c,1),a,b,c
!          end if
!
!          do i = 1,dir
!            if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
!              fi(a,b,c,i) = crs_fi(a,b,c,i)
!              gi(a,b,c,i) = crs_gi(a,b,c,i)
!            end if
!          end do
!
!          end if
!        end do
!      end do
!    end do
