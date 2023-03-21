MODULE amr_info_holder
!
! Hold the information for all major calculations used by amrex.
!
!
!
!
USE precise
use amrex_amr_module
use amrex_base_module
use output_data, only : box_zones
use nml_output, only : multigrid_save_info
IMPLICIT NONE
SAVE
private
!
! Contiguous pointers
public :: rho,fi,fout,u_vel,v_vel,w_vel,temp,alpha,gi,gout,press,state
public :: chi,zeta_x,zeta_y,zeta_z,lambda_x,lambda_y,lambda_z,pi_xx,&
  pi_yy,pi_zz,pi_xy,pi_xz,pi_yz
!
! Time information
public :: time,time_prev,dt,ext_dx,timestep,max_timesteps,elbm_epsilon,init_epsilon,startup_iterations
!
! Multifabs
public :: mfrho,mfu_vel,mfv_vel,mfw_vel,mffi,mffout,mfgi,mfgout,&
  mftemp,mfalpha,mfpress,mfstate
public :: mfchi,mfzeta_x,mfzeta_y,mfzeta_z,mflambda_x,&
  mflambda_y,mflambda_z,mfpi_xx,mfpi_yy,mfpi_zz,mfpi_xy,mfpi_xz,mfpi_yz
!
! Misc
public :: n_boxes_lvl,zone_storage,grid_info,fixed_lvl,ggm_save_info
public :: nghosts,nghosts_state,nghosts_mid
public :: x_comps,y_comps,z_comps,nodal_domain_limits
!
! Input information
public :: plot_interval,regrid_interval,ref_rat,block_factor,max_grid,proper_shell
public :: boi_lo,boi_hi,box_lvl,coarse_box,amr,num_boi,multigrid,startup,cumulative_zones
public :: regrid_times,regrid_timesteps,coarse_regrid_done
!
! Useful AMR stuff
!
real(kind=dp),allocatable :: boi_lo(:,:),boi_hi(:,:)
integer :: num_boi,cumulative_zones
integer,allocatable :: box_lvl(:),nodal_domain_limits(:,:,:),coarse_box(:)
logical :: amr,multigrid,startup,coarse_regrid_done
!
! Arrays iterated in, pointed to with mfiters and the multifabs
!
REAL(KIND=dp),contiguous,pointer :: rho(:,:,:,:),&
  fi(:,:,:,:),fout(:,:,:,:),u_vel(:,:,:,:),v_vel(:,:,:,:),w_vel(:,:,:,:),&
  temp(:,:,:,:),alpha(:,:,:,:),gi(:,:,:,:),gout(:,:,:,:),press(:,:,:,:)
!      real(kind=dp),contiguous,pointer :: big_qbar(:,:,:,:),little_qbar(:,:,:,:),&
!        pressure_tensor(:,:,:,:),shear_tensor(:,:,:,:)
real(kind=dp),contiguous,pointer :: chi(:,:,:,:),zeta_x(:,:,:,:),zeta_y(:,:,:,:),&
  zeta_z(:,:,:,:),lambda_x(:,:,:,:),lambda_y(:,:,:,:),lambda_z(:,:,:,:),&
  pi_xx(:,:,:,:),pi_yy(:,:,:,:),pi_zz(:,:,:,:),pi_xy(:,:,:,:),pi_xz(:,:,:,:),&
  pi_yz(:,:,:,:)
INTEGER,CONTIGUOUS,POINTER :: state(:,:,:,:)
!
! Time and dx information
!
real(kind=dp),allocatable :: time(:),time_prev(:),dt(:),ext_dx(:),elbm_epsilon(:),&
  init_epsilon(:),regrid_times(:),regrid_timesteps(:,:)
integer,allocatable :: timestep(:),max_timesteps(:)
!
! Create the multifabs required for amrex
! These are Allocated in init_things()
!
TYPE(amrex_multifab),ALLOCATABLE :: mfrho(:),mffi(:),mffout(:),&
  mfu_vel(:),mfv_vel(:),mfw_vel(:),mftemp(:),mfalpha(:),mfgi(:),&
  mfgout(:),mfpress(:)
type(amrex_imultifab),allocatable :: mfstate(:)
type(amrex_multifab),ALLOCATABLE :: mfchi(:),mfzeta_x(:),mfzeta_y(:),mfzeta_z(:),&
  mflambda_x(:),mflambda_y(:),mflambda_z(:),mfpi_xx(:),mfpi_yy(:),mfpi_zz(:),&
  mfpi_xy(:),mfpi_xz(:),mfpi_yz(:)
!
! Number of boxes at each level
!
integer, allocatable :: n_boxes_lvl(:)
!
! User defined type used to hold data for output
!
type(box_zones),allocatable :: zone_storage(:),grid_info(:,:)
type(box_zones),allocatable :: ggm_save_info(:,:)
!
!
!
!type(multigrid_save_info),allocatable :: multigrid_save(:),ggm_save_data(:)
!
!
!
integer :: nghosts,nghosts_state,nghosts_mid
integer :: x_comps,y_comps,z_comps,startup_iterations
integer :: plot_interval,regrid_interval,ref_rat,block_factor,max_grid,proper_shell,fixed_lvl

END module
