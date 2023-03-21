MODULE amr_processes
!
! Module containing all the virtual processes used for AMR by
! Amrex.
!
!
!
use iso_c_binding
use precise
use timing
use constants
use amrex_amr_module
use amr_info_holder
use ggm_stuff
use amrex_base_module
use amrex_fort_module, only: rt => amrex_real
use grid_data
use fill_stuff, only : edge_bounds,edge_bounds_distro,lo_bc,hi_bc,lo_bc_fi,hi_bc_fi
IMPLICIT NONE
SAVE
PRIVATE
PUBLIC :: init_amrex,from_coarse,from_scratch,clear_lvl,remake,&
  tag_nodes,fill_qone_interlevel
public :: count_big_zones,get_real_coords,fieq_comp_19,fieq_comp_39,total_big_zones
public :: overlapping_boxes,ml_overlapping_boxes,fieq_comp_39_shifted
public :: self,nprocs,commune !MPI information
public :: amr_max_lvl,finest_lvl !amr useful info
public :: magic_box,mpi_amrex_retrieve,magic_real_box,get_box_sizes,magic_box_exclusive
public :: node_based_3d,node_based_2d,direction_tester
public :: beta_1_calculator,beta_2_calculator,epsilon_calculator,zoom_zoom_calculator
public :: chi_estimate,zeta_x_estimate,zeta_y_estimate,zeta_z_estimate,pi_xx_estimate,&
  pi_yy_estimate,pi_zz_estimate,pi_xy_estimate,pi_xz_estimate,pi_yz_estimate,lambda_x_estimate,&
  lambda_y_estimate,lambda_z_estimate

integer :: self,nprocs,commune,amr_max_lvl,finest_lvl
!integer,parameter :: ref_rat = 2
logical,parameter,dimension(3) :: node_based_3d = (/.true.,.true.,.true./)
logical,parameter,dimension(2) :: node_based_2d = (/.true.,.true./)

CONTAINS

SUBROUTINE init_amrex
!
!
!
!
!
!
!use ggm_stuff
IMPLICIT NONE

INTEGER :: i,lev,finer
type(amrex_parmparse) :: pp

WRITE(11,*) 'Initializing virtual functions'
!write(*,*) 'Initializing virtual functions',dir
CALL amrex_init_virtual_functions(from_scratch,from_coarse,&
  remake,clear_lvl,tag_nodes)
finer = amrex_get_finest_level()
write(11,*) 'Virtual functions initialized.'
!write(*,*) 'Virtual functions initialized.',finer
!write(11,*) 'Preparing to read input.'
!CALL input_read()
!WRITE(11,*) 'Input successfully read.'
call mpi_amrex_retrieve

!write(*,*) 'Allocating multifabs up to level',amrex_max_level
allocate(mffi(0:amrex_max_level))
allocate(mffout(0:amrex_max_level))
allocate(mfgi(0:amrex_max_level))
allocate(mfgout(0:amrex_max_level))
!
! Flow variables
!
allocate(mfrho(0:amrex_max_level))
allocate(mfu_vel(0:amrex_max_level))
allocate(mfv_vel(0:amrex_max_level))
allocate(mftemp(0:amrex_max_level))
allocate(mfpress(0:amrex_max_level))
!
! Allocate the Lagrange multipliers
!
allocate(mfchi(0:amrex_max_level))
allocate(mfzeta_x(0:amrex_max_level))
allocate(mfzeta_y(0:amrex_max_level))
allocate(mfpi_xx(0:amrex_max_level))
allocate(mfpi_yy(0:amrex_max_level))
allocate(mfpi_xy(0:amrex_max_level))
allocate(mflambda_x(0:amrex_max_level))
allocate(mflambda_y(0:amrex_max_level))
!
! Allocating things which don't fit nicely into other categories
!
allocate(mfstate(0:amrex_max_level))
allocate(mfalpha(0:amrex_max_level))
!
! Other random things
!
allocate(time(0:amrex_max_level))
allocate(time_prev(0:amrex_max_level))
allocate(dt(0:amrex_max_level))
allocate(rdx(0:amrex_max_level))
allocate(ext_dx(0:amrex_max_level))
allocate(timestep(0:amrex_max_level))
allocate(max_timesteps(0:amrex_max_level))
allocate(elbm_epsilon(0:amrex_max_level))
allocate(init_epsilon(0:amrex_max_level))
allocate(nodal_domain_limits(3,2,0:amrex_max_level))
!allocate(dx(0:max_level))
!
! For 3D runs
!
if (dimensions == 3) then
  allocate(mfw_vel(0:amrex_max_level))
  allocate(mfzeta_z(0:amrex_max_level))
  allocate(mfpi_xz(0:amrex_max_level))
  allocate(mfpi_yz(0:amrex_max_level))
  allocate(mfpi_zz(0:amrex_max_level))
  allocate(mflambda_z(0:amrex_max_level))
end if

if (use_ggm) then
  allocate(mfgm(0:amrex_max_level))
  allocate(mfggm(0:amrex_max_level))
end if

!ALLOCATE(
!allocate(nzones_step(nsteps))

!WRITE(*,*) 'Variables successfully allocated.'
!WRITE(*,*) lbound(mfrho),ubound(mfrho)
!DO i = 0,n_levels-1
!  CALL amrex_multifab_build(mfrho(0),
!END DO

!CALL amrex_init_virtual_functions(from_scratch,from_coarse,&
!  remake,clear_lvl,tag_nodes)


! write a multifab build
      !CALL amrex_mfiter_build(mfi,mfrho(0),tiling=.false.)

rdx = amrex_geom%dx(1)

END SUBROUTINE
!
!
!
!
!
!
!
!
SUBROUTINE from_coarse(lvl,t,boxes,distro) bind(c)
!
!
!
!use fill_stuff, only : fill_patch,fill_coarse
!use amr_grid_writing, only: amr_grid_zones
IMPLICIT NONE
TYPE(amrex_boxarray) :: box
TYPE(amrex_distromap) :: dmap
TYPE(c_ptr), intent(in), value :: boxes,distro
INTEGER,intent(in),value :: lvl
!INTEGER,parameter :: s_comp=1,d_comp=1
integer :: ncomp,nghost,base_comp,finer,dummy,new_fine,ier
REAL(KIND=dp),intent(in),value :: t
real(kind=dp) :: finest_time,dummy_time

finer = amrex_get_finest_level()
finest_time = time(finer)
WRITE(*,*) 'Using from_coarse',finer,lvl,t,nghosts

box = boxes
dmap = distro

dummy_time = time(lvl-1) + 0.00001D0

!call amrex_print(dmap)

call clear_lvl(lvl)

call mpi_barrier(commune,ier)

time(lvl) = t
dummy = 0
!s_comp = 1
!d_comp = 1
!ncomp = 1
!nghost = 0

!ALLOCATE(mfrho_old(0:n_levels-1))
!mfrho_old(level) = mfrho(level)
!
! Build the flow variables
!
ncomp = 1
call amrex_multifab_build(mfrho(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfu_vel(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfv_vel(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mftemp(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpress(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
if (dimensions == 3) then
  call amrex_multifab_build(mfw_vel(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!
! Build the LMs
!
call amrex_multifab_build(mfchi(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_x(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_y(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xx(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xy(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_yy(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_x(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_y(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)

! For the 3D elements
if (dimensions == 3) then
  call amrex_multifab_build(mfzeta_z(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_xz(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_yz(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_zz(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mflambda_z(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!
! Other needs
!
call amrex_multifab_build(mfalpha(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)

if (use_ggm) then
  call amrex_multifab_build(mfgm(lvl),box,dmap,num_ggm_vars,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfggm(lvl),box,dmap,num_ggm_vars,nghosts_mid,node_based_3d)
end if
!
! Build the multi-directional variables
!
ncomp = dir
call amrex_imultifab_build(mfstate(lvl),box,dmap,ncomp,nghosts_state,node_based_3d)
!
!
!
ncomp = dir
call amrex_multifab_build(mffi(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mffout(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgi(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgout(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)



!
!
write(*,*) 'pre-fill',time_prev(lvl-1),time(lvl-1),t,time_prev(lvl),time(lvl)
!
!
  call amrex_fillcoarsepatch(mftemp(lvl),time_prev(lvl-1),mftemp(lvl-1),time(lvl-1),mftemp(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfrho(lvl),time_prev(lvl-1),mfrho(lvl-1),time(lvl-1),mfrho(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfu_vel(lvl),time_prev(lvl-1),mfu_vel(lvl-1),time(lvl-1),mfu_vel(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfv_vel(lvl),time_prev(lvl-1),mfv_vel(lvl-1),time(lvl-1),mfv_vel(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfalpha(lvl),time_prev(lvl-1),mfalpha(lvl-1),time(lvl-1),mfalpha(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfchi(lvl),time_prev(lvl-1),mfchi(lvl-1),time(lvl-1),mfchi(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfzeta_x(lvl),time_prev(lvl-1),mfzeta_x(lvl-1),time(lvl-1),mfzeta_x(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfzeta_y(lvl),time_prev(lvl-1),mfzeta_y(lvl-1),time(lvl-1),mfzeta_y(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfpi_xx(lvl),time_prev(lvl-1),mfpi_xx(lvl-1),time(lvl-1),mfpi_xx(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfpi_yy(lvl),time_prev(lvl-1),mfpi_yy(lvl-1),time(lvl-1),mfpi_yy(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfpi_xy(lvl),time_prev(lvl-1),mfpi_xy(lvl-1),time(lvl-1),mfpi_xy(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mflambda_x(lvl),time_prev(lvl-1),mflambda_x(lvl-1),time(lvl-1),mflambda_x(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mflambda_y(lvl),time_prev(lvl-1),mflambda_y(lvl-1),time(lvl-1),mflambda_y(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
!  call amrex_fillcoarsepatch(mfalpha(lvl),time_prev(lvl-1),mfalpha(lvl-1),t,mfalpha(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,t,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
if (dimensions == 3) then
!
  call amrex_fillcoarsepatch(mfw_vel(lvl),time_prev(lvl-1),mfw_vel(lvl-1),time(lvl-1),mfw_vel(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfzeta_z(lvl),time_prev(lvl-1),mfzeta_z(lvl-1),time(lvl-1),mfzeta_z(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfpi_zz(lvl),time_prev(lvl-1),mfpi_zz(lvl-1),time(lvl-1),mfpi_zz(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfpi_xz(lvl),time_prev(lvl-1),mfpi_xz(lvl-1),time(lvl-1),mfpi_xz(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mfpi_yz(lvl),time_prev(lvl-1),mfpi_yz(lvl-1),time(lvl-1),mfpi_yz(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
  call amrex_fillcoarsepatch(mflambda_z(lvl),time_prev(lvl-1),mflambda_z(lvl-1),time(lvl-1),mflambda_z(lvl-1),&
    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,time(lvl-1),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
end if
!!
!!
!!
  call amrex_fillcoarsepatch(mffi(lvl),time_prev(lvl-1),mffi(lvl-1),time(lvl-1),mffi(lvl-1),&
    amrex_geom(lvl-1),edge_bounds_distro,amrex_geom(lvl),edge_bounds_distro,time(lvl-1),1,1,dir,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
  call amrex_fillcoarsepatch(mfgi(lvl),time_prev(lvl-1),mfgi(lvl-1),time(lvl-1),mfgi(lvl-1),&
    amrex_geom(lvl-1),edge_bounds_distro,amrex_geom(lvl),edge_bounds_distro,time(lvl-1),1,1,dir,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
  call amrex_fillcoarsepatch(mffout(lvl),time_prev(lvl-1),mffout(lvl-1),time(lvl-1),mffout(lvl-1),&
    amrex_geom(lvl-1),edge_bounds_distro,amrex_geom(lvl),edge_bounds_distro,time(lvl-1),1,1,dir,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
  call amrex_fillcoarsepatch(mfgout(lvl),time_prev(lvl-1),mfgout(lvl-1),time(lvl-1),mfgout(lvl-1),&
    amrex_geom(lvl-1),edge_bounds_distro,amrex_geom(lvl),edge_bounds_distro,time(lvl-1),1,1,dir,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
!  call amrex_fillcoarsepatch(mffi(lvl),time_prev(lvl-1),mffi(lvl-1),t,mffi(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,t,1,1,dir,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!!
!  call amrex_fillcoarsepatch(mfgi(lvl),time_prev(lvl-1),mfgi(lvl-1),t,mfgi(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,t,1,1,dir,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!!
!  call amrex_fillcoarsepatch(mffout(lvl),time_prev(lvl-1),mffout(lvl-1),t,mffout(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,t,1,1,dir,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!!
!  call amrex_fillcoarsepatch(mfgout(lvl),time_prev(lvl-1),mfgout(lvl-1),t,mfgout(lvl-1),&
!    amrex_geom(lvl-1),edge_bounds,amrex_geom(lvl),edge_bounds,t,1,1,dir,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
END SUBROUTINE
!
!
!
!
!
!
!
!
!
SUBROUTINE from_scratch(level,t,boxes,distro) bind(c)
!
! Creates all the multifabs which were allocated in init_amrex
!
IMPLICIT NONE
TYPE(amrex_boxarray) :: box
TYPE(amrex_distromap) :: dmap
TYPE(c_ptr), intent(in), value :: boxes,distro
TYPE(amrex_mfiter) :: mfi
!      TYPE(amrex_box) :: b1
!      REAL(KIND=dp),CONTIGUOUS,POINTER ::
INTEGER,intent(in),value :: level
INTEGER :: ncomp,finer,ier!,nghost
REAL(amrex_real),intent(in),value :: t

finer = amrex_get_finest_level()

WRITE(*,*) 'Using from_scratch',level,nghosts
box = boxes
dmap = distro

call clear_lvl(level)

call mpi_barrier(commune,ier)

!write(*,*) 'returning in scratch from clear_lvl'
time(level) = t
!write(*,*) 'time updated'
ncomp = 1
!write(*,*) 'initializing flow variables'
!write(*,*) 'building multifabs on level ',level,' up to level ',finer,&
!  ' for max level ',amrex_max_level
call amrex_multifab_build(mfrho(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfu_vel(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfv_vel(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mftemp(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpress(level),box,dmap,ncomp,nghosts_mid,node_based_3d)

if (dimensions == 3) then
  call amrex_multifab_build(mfw_vel(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!
! Build the LMs
!
call amrex_multifab_build(mfchi(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_x(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_y(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xx(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xy(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_yy(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_x(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_y(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
!
! For the 3D elements
!
if (dimensions == 3) then
  call amrex_multifab_build(mfzeta_z(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_xz(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_yz(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_zz(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mflambda_z(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!
! Other needs
!
call amrex_multifab_build(mfalpha(level),box,dmap,ncomp,nghosts_mid,node_based_3d)

if (use_ggm) then
  call amrex_multifab_build(mfgm(level),box,dmap,num_ggm_vars,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfggm(level),box,dmap,num_ggm_vars,nghosts_mid,node_based_3d)
end if
!
! Build the multi-directional variables
!
ncomp = dir
call amrex_imultifab_build(mfstate(level),box,dmap,ncomp,nghosts_state,node_based_3d)
call amrex_multifab_build(mffi(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mffout(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgi(level),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgout(level),box,dmap,ncomp,nghosts_mid,node_based_3d)



END SUBROUTINE
!
!
!
!
!
!
!
!
!
!
!
SUBROUTINE clear_lvl(level) bind(c)
!
!
!
IMPLICIT NONE
INTEGER,intent(in),value :: level
integer :: finer
finer = amrex_get_finest_level()
WRITE(*,*) 'Using clear_lvl',level
!
! Flow variables
!
!write(*,*) 'clearing flow variables'
call amrex_multifab_destroy(mfrho(level))
call amrex_multifab_destroy(mfu_vel(level))
call amrex_multifab_destroy(mfv_vel(level))
call amrex_multifab_destroy(mftemp(level))
call amrex_multifab_destroy(mfpress(level))
!
! LMs
!
!write(*,*) 'clearing lagrange multipliers'
call amrex_multifab_destroy(mfchi(level))
call amrex_multifab_destroy(mfzeta_x(level))
call amrex_multifab_destroy(mfzeta_y(level))
call amrex_multifab_destroy(mfpi_xx(level))
call amrex_multifab_destroy(mfpi_yy(level))
call amrex_multifab_destroy(mfpi_xy(level))
call amrex_multifab_destroy(mflambda_x(level))
call amrex_multifab_destroy(mflambda_y(level))
!
! Distro funcs
!
!write(*,*) 'clearing distro funcs'
call amrex_multifab_destroy(mffi(level))
call amrex_multifab_destroy(mfgi(level))
call amrex_multifab_destroy(mffout(level))
call amrex_multifab_destroy(mfgout(level))
!
! Other
!
!write(*,*) 'clearing other stuff'
call amrex_multifab_destroy(mfalpha(level))
!write(*,*) 'clearing state'
call amrex_imultifab_destroy(mfstate(level))

! 3D stuff
!write(*,*) 'clearing 3D stuff'
if (dimensions == 3) then
  call amrex_multifab_destroy(mfw_vel(level))
  call amrex_multifab_destroy(mfzeta_z(level))
  call amrex_multifab_destroy(mfpi_xz(level))
  call amrex_multifab_destroy(mfpi_yz(level))
  call amrex_multifab_destroy(mfpi_zz(level))
  call amrex_multifab_destroy(mflambda_z(level))
end if

if (use_ggm) then
  call amrex_multifab_destroy(mfgm(level))
  call amrex_multifab_destroy(mfggm(level))
end if

END SUBROUTINE
!
!
!
!
!
!
!
!
SUBROUTINE tag_nodes(level,cp,t,settag,cleartag) bind(c)
!
! This is the virtual class for tagging nodes before refinement/coarsening
!
use ggm_stuff
IMPLICIT NONE
INTEGER,INTENT(IN),value :: level
INTEGER :: taglow(4),taghigh(4),finer,max_ggm_lvl,i,ier
REAL(amrex_real),intent(in),value :: t
logical :: initial
TYPE(c_ptr),intent(in),value :: cp
CHARACTER(kind = c_char),intent(in),value :: settag,cleartag
CHARACTER(kind = c_char),contiguous,pointer :: tagarr(:,:,:,:)

TYPE(amrex_tagboxarray) :: tag
TYPE(amrex_box) :: box_1
TYPE(amrex_mfiter) :: mfi

tag = cp
finer = amrex_get_finest_level()
WRITE(*,*) 'Using tag_nodes',level,t,finer,amrex_max_level
!
! For initial multigrid using the Coarse Box of Interest information.
!
if (t < very_small_number) then
  initial = .true.
  !write(*,*) 'steve'
!  if (level >= maxval(box_lvl)) return


  do i = 1,num_boi

    if (level >= coarse_box(i)) cycle

    CALL amrex_mfiter_build(mfi,mfrho(level),tiling=.false.)
    do while (mfi%next())
      box_1 = mfi%validbox()
!    if (level > 0) then
!      call amrex_print(box_1)
!    end if
      tagarr => tag%dataptr(mfi)
      taglow = lbound(tagarr)
      taghigh = ubound(tagarr)

      call initial_multigrid_creation(level,settag,cleartag,tagarr,&
        box_1%lo,box_1%hi,taglow,taghigh,initial)
    end do
  end do
  call amrex_mfiter_destroy(mfi)
!
!
!
else if (t <= coarse_time + dt(0) .and. t >= coarse_time - dt(0) &
  .and. coarse_time > very_small_number .and. .not. coarse_regrid_done) then
!  write(*,*) 'spargle margle!',level
initial = .false.
  do i = 1,num_boi
    if (level >= box_lvl(i)) cycle

    CALL amrex_mfiter_build(mfi,mfrho(level),tiling=.false.)
    do while (mfi%next())
      box_1 = mfi%validbox()
!    call amrex_print(box_1)

      tagarr => tag%dataptr(mfi)
      taglow = lbound(tagarr)
      taghigh = ubound(tagarr)

      call initial_multigrid_creation(level,settag,cleartag,tagarr,&
        box_1%lo,box_1%hi,taglow,taghigh,initial)
    end do

    call amrex_mfiter_destroy(mfi)

  end do

!  coarse_regrid_done = .true.
!
! For the AMR part using GGM
!
else

CALL amrex_mfiter_build(mfi,mfrho(level),tiling=.false.)
!write(*,*) 'giggles all around!!!'
  do while (mfi%next())
    box_1 = mfi%validbox()
!    call amrex_print(box_1)

    tagarr => tag%dataptr(mfi)
    taglow = lbound(tagarr)
    taghigh = ubound(tagarr)
!
! For pure gradient refinement
!
    if (gradient_refine_only) then
      call gradient_amr(level,settag,cleartag,tagarr,&
        box_1%lo,box_1%hi,taglow,taghigh,mfi)
    else

      call ggm_calculator(level,settag,cleartag,tagarr,&
        box_1%lo,box_1%hi,taglow,taghigh,mfi)

    end if

  end do

call mpi_barrier(commune,ier)
call amrex_mfiter_destroy(mfi)

end if
!
END SUBROUTINE
!
!
!
! The heart of the AMR part of this, remakes a level.
!
!
!
!
!
SUBROUTINE remake(lvl,t,boxes,distro) bind(c)
!
! Remakes a level
!
IMPLICIT NONE
INTEGER,intent(in),value :: lvl
INTEGER :: ncomp,dummy,ier
REAL(KIND=dp),intent(in),value :: t
logical :: tofrom

TYPE(amrex_boxarray) :: box
TYPE(amrex_distromap) :: dmap
TYPE(amrex_multifab) :: mfrho_temp,mfu_vel_temp,mfv_vel_temp,mfw_vel_temp,&
  mftemp_temp,mfpress_temp,mfalpha_temp
type(amrex_multifab) :: mfchi_temp,mfzeta_x_temp,mfzeta_y_temp,mfzeta_z_temp,&
  mfpi_xx_temp,mfpi_yy_temp,mfpi_zz_temp,mfpi_xy_temp,mfpi_xz_temp,mfpi_yz_temp,&
  mflambda_x_temp,mflambda_y_temp,mflambda_z_temp
type(amrex_multifab) :: mffi_temp,mfgi_temp,mffout_temp,mfgout_temp
type(amrex_imultifab) :: imfstate_temp
TYPE(c_ptr), intent(in), value :: boxes,distro

WRITE(*,*) 'Using remake',lvl
box = boxes
dmap = distro

ncomp = 1
dummy = 0

call amrex_multifab_build(mfrho_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfu_vel_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfv_vel_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mftemp_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
if (dimensions == 3) then
  call amrex_multifab_build(mfw_vel_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
end if

call amrex_multifab_build(mfalpha_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
!!
!! Build the LMs
!!
call amrex_multifab_build(mfchi_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_x_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_y_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xx_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xy_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_yy_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_x_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_y_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
!
!! For the 3D elements
if (dimensions == 3) then
  call amrex_multifab_build(mfzeta_z_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_xz_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_yz_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_zz_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mflambda_z_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!!
!! Other needs
!!
!call amrex_multifab_build(mfalpha_temp,box,dmap,ncomp,nghost,node_based_3d)
!!
!! Build the multi-directional variables
!!
ncomp = dir
call amrex_imultifab_build(imfstate_temp,box,dmap,ncomp,nghosts_state,node_based_3d)
call amrex_multifab_build(mffi_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mffout_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgi_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgout_temp,box,dmap,ncomp,nghosts_mid,node_based_3d)
!
! tofrom is a switch to copy info to (if true) the temporary state imultifab, or
!   from (if false).  It's ugly, but there's not much choice here.
!
!tofrom = .true.
!call state_copy(imfstate_temp,tofrom)
!
!
!
if (lvl == 0) then
  call amrex_fillpatch(mfrho_temp,time_prev(lvl),mfrho(lvl),time(lvl),mfrho(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfu_vel_temp,time_prev(lvl),mfu_vel(lvl),time(lvl),mfu_vel(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfv_vel_temp,time_prev(lvl),mfv_vel(lvl),time(lvl),mfv_vel(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfw_vel_temp,time_prev(lvl),mfw_vel(lvl),time(lvl),mfw_vel(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mftemp_temp,time_prev(lvl),mftemp(lvl),time(lvl),mftemp(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfalpha_temp,time_prev(lvl),mfalpha(lvl),time(lvl),mfalpha(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
!
!
!
  call amrex_fillpatch(mfchi_temp,time_prev(lvl),mfchi(lvl),time(lvl),mfchi(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfzeta_x_temp,time_prev(lvl),mfzeta_x(lvl),time(lvl),mfzeta_x(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfzeta_y_temp,time_prev(lvl),mfzeta_y(lvl),time(lvl),mfzeta_y(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfzeta_z_temp,time_prev(lvl),mfzeta_z(lvl),time(lvl),mfzeta_z(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfpi_xx_temp,time_prev(lvl),mfpi_xx(lvl),time(lvl),mfpi_xx(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfpi_yy_temp,time_prev(lvl),mfpi_yy(lvl),time(lvl),mfpi_yy(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfpi_zz_temp,time_prev(lvl),mfpi_zz(lvl),time(lvl),mfpi_zz(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfpi_xy_temp,time_prev(lvl),mfpi_xy(lvl),time(lvl),mfpi_xy(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfpi_xz_temp,time_prev(lvl),mfpi_xz(lvl),time(lvl),mfpi_xz(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mfpi_yz_temp,time_prev(lvl),mfpi_yz(lvl),time(lvl),mfpi_yz(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mflambda_x_temp,time_prev(lvl),mflambda_x(lvl),time(lvl),mflambda_x(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mflambda_y_temp,time_prev(lvl),mflambda_y(lvl),time(lvl),mflambda_y(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mflambda_z_temp,time_prev(lvl),mflambda_z(lvl),time(lvl),mflambda_z(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,1)
  call amrex_fillpatch(mffi_temp,time_prev(lvl),mffi(lvl),time(lvl),mffi(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,dir)
  call amrex_fillpatch(mfgi_temp,time_prev(lvl),mfgi(lvl),time(lvl),mfgi(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,dir)
  call amrex_fillpatch(mffout_temp,time_prev(lvl),mffout(lvl),time(lvl),mffout(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,dir)
  call amrex_fillpatch(mfgout_temp,time_prev(lvl),mfgout(lvl),time(lvl),mfgout(lvl),amrex_geom(lvl),&
    edge_bounds,time(lvl),1,1,dir)
!
else
!  write(*,*) 'stink'
  call amrex_fillpatch(mfrho_temp,time_prev(lvl-1),mfrho(lvl-1),time(lvl-1),mfrho(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfrho(lvl),time(lvl),mfrho(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  write(*,*) 'hack hack cough'
  call amrex_fillpatch(mfu_vel_temp,time_prev(lvl-1),mfu_vel(lvl-1),time(lvl-1),mfu_vel(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfu_vel(lvl),time(lvl),mfu_vel(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  write(*,*) 'ms. dodger?'
  call amrex_fillpatch(mfv_vel_temp,time_prev(lvl-1),mfv_vel(lvl-1),time(lvl-1),mfv_vel(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfv_vel(lvl),time(lvl),mfv_vel(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  write(*,*) 'sparklefingers'
  call amrex_fillpatch(mftemp_temp,time_prev(lvl-1),mftemp(lvl-1),time(lvl-1),mftemp(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mftemp(lvl),time(lvl),mftemp(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mfalpha_temp,time_prev(lvl-1),mfalpha(lvl-1),time(lvl-1),mfalpha(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfalpha(lvl),time(lvl),mfalpha(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!  write(*,*) 'sissypants'
!  call amrex_fillpatch(mfpress_temp,t,mftemp(lvl-1),t,mftemp(lvl-1),amrex_geom(lvl-1),&
!    edge_bounds,t,mftemp(lvl),t,mftemp(lvl),amrex_geom(lvl),edge_bounds,&
!    t,1,1,1,amrex_ref_ratio(lvl-1),&
!    amrex_interp_node_bilinear,lo_bc,hi_bc)
!    write(*,*) 'special day'
  call amrex_fillpatch(mfchi_temp,time_prev(lvl-1),mfchi(lvl-1),time(lvl-1),mfchi(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfchi(lvl),time(lvl),mfchi(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!    write(*,*) 'rumple hob knob'
  call amrex_fillpatch(mfzeta_x_temp,time_prev(lvl-1),mfzeta_x(lvl-1),time(lvl-1),mfzeta_x(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfzeta_x(lvl),time(lvl),mfzeta_x(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mfzeta_y_temp,time_prev(lvl-1),mfzeta_y(lvl-1),time(lvl-1),mfzeta_y(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfzeta_y(lvl),time(lvl),mfzeta_y(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mfpi_xx_temp,time_prev(lvl-1),mfpi_xx(lvl-1),time(lvl-1),mfpi_xx(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfpi_xx(lvl),time(lvl),mfpi_xx(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mfpi_yy_temp,time_prev(lvl-1),mfpi_yy(lvl-1),time(lvl-1),mfpi_yy(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfpi_yy(lvl),time(lvl),mfpi_yy(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mfpi_xy_temp,time_prev(lvl-1),mfpi_xy(lvl-1),time(lvl-1),mfpi_xy(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mfpi_xy(lvl),time(lvl),mfpi_xy(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mflambda_x_temp,time_prev(lvl-1),mflambda_x(lvl-1),time(lvl-1),mflambda_x(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mflambda_x(lvl),time(lvl),mflambda_x(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
  call amrex_fillpatch(mflambda_y_temp,time_prev(lvl-1),mflambda_y(lvl-1),time(lvl-1),mflambda_y(lvl-1),amrex_geom(lvl-1),&
    edge_bounds,time_prev(lvl),mflambda_y(lvl),time(lvl),mflambda_y(lvl),amrex_geom(lvl),edge_bounds,&
    time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
    amrex_interp_node_bilinear,lo_bc,hi_bc)
!
!
!
  if (dimensions == 3) then
    call amrex_fillpatch(mfw_vel_temp,time_prev(lvl-1),mfw_vel(lvl-1),time(lvl-1),mfw_vel(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mfw_vel(lvl),time(lvl),mfw_vel(lvl),amrex_geom(lvl),edge_bounds,&
      time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc,hi_bc)
!
    call amrex_fillpatch(mfpi_zz_temp,time_prev(lvl-1),mfpi_zz(lvl-1),time(lvl-1),mfpi_zz(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mfpi_zz(lvl),time(lvl),mfpi_zz(lvl),amrex_geom(lvl),edge_bounds,&
      time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc,hi_bc)
!
    call amrex_fillpatch(mfzeta_z_temp,time_prev(lvl-1),mfzeta_z(lvl-1),time(lvl-1),mfzeta_z(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,t,mfzeta_z(lvl),time(lvl),mfzeta_z(lvl),amrex_geom(lvl),edge_bounds,&
      time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc,hi_bc)
!
    call amrex_fillpatch(mfpi_xz_temp,time_prev(lvl-1),mfpi_xz(lvl-1),time(lvl-1),mfpi_xz(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mfpi_xz(lvl),time(lvl),mfpi_xz(lvl),amrex_geom(lvl),edge_bounds,&
      time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc,hi_bc)
!
    call amrex_fillpatch(mfpi_yz_temp,time_prev(lvl-1),mfpi_yz(lvl-1),time(lvl-1),mfpi_yz(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mfpi_yz(lvl),time(lvl),mfpi_yz(lvl),amrex_geom(lvl),edge_bounds,&
      time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc,hi_bc)
!
    call amrex_fillpatch(mflambda_z_temp,time_prev(lvl-1),mflambda_z(lvl-1),time(lvl-1),mflambda_z(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mflambda_z(lvl),time(lvl),mflambda_z(lvl),amrex_geom(lvl),edge_bounds,&
      time(lvl),1,1,1,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc,hi_bc)
  end if
!
    call amrex_fillpatch(mffi_temp,time_prev(lvl-1),mffi(lvl-1),time(lvl-1),mffi(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mffi(lvl),time(lvl),mffi(lvl),amrex_geom(lvl),edge_bounds_distro,&
      time(lvl),1,1,dir,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
    call amrex_fillpatch(mfgi_temp,time_prev(lvl-1),mfgi(lvl-1),time(lvl-1),mfgi(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mfgi(lvl),time(lvl),mfgi(lvl),amrex_geom(lvl),edge_bounds_distro,&
      time(lvl),1,1,dir,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
    call amrex_fillpatch(mffout_temp,time_prev(lvl-1),mffout(lvl-1),time(lvl-1),mffout(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mffout(lvl),time(lvl),mffout(lvl),amrex_geom(lvl),edge_bounds_distro,&
      time(lvl),1,1,dir,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)
!
    call amrex_fillpatch(mfgout_temp,time_prev(lvl-1),mfgout(lvl-1),time(lvl-1),mfgout(lvl-1),amrex_geom(lvl-1),&
      edge_bounds,time_prev(lvl),mfgout(lvl),time(lvl),mfgout(lvl),amrex_geom(lvl),edge_bounds_distro,&
      time(lvl),1,1,dir,amrex_ref_ratio(lvl-1),&
      amrex_interp_node_bilinear,lo_bc_fi,hi_bc_fi)

end if
!write(*,*) 'speak booboo',lvl
!
! Nuke what's there
!
call clear_lvl(lvl)
call mpi_barrier(commune,ier)
!
!
!
ncomp = 1
call amrex_multifab_build(mfrho(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfu_vel(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfv_vel(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mftemp(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpress(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)

if (dimensions == 3) then
  call amrex_multifab_build(mfw_vel(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!
! Build the LMs
!
call amrex_multifab_build(mfchi(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_x(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfzeta_y(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xx(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_xy(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfpi_yy(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_x(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mflambda_y(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
!
! For the 3D elements
!
if (dimensions == 3) then
  call amrex_multifab_build(mfzeta_z(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_xz(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_yz(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfpi_zz(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mflambda_z(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
end if
!
! Other needs
!
call amrex_multifab_build(mfalpha(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)

if (use_ggm) then
  call amrex_multifab_build(mfgm(lvl),box,dmap,num_ggm_vars,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfggm(lvl),box,dmap,num_ggm_vars,nghosts_mid,node_based_3d)
end if
!
! Build the multi-directional variables
!
ncomp = dir
call amrex_imultifab_build(mfstate(lvl),box,dmap,ncomp,nghosts_state,node_based_3d)
!
call amrex_multifab_build(mffi(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mffout(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgi(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
call amrex_multifab_build(mfgout(lvl),box,dmap,ncomp,nghosts_mid,node_based_3d)
!
!
!
!tofrom = .false.
!call state_copy(imfstate_temp,tofrom)
!
!
!
call mfrho(lvl)%copy(mfrho_temp,1,1,1,nghosts_mid)
call mfu_vel(lvl)%copy(mfu_vel_temp,1,1,1,nghosts_mid)
call mfv_vel(lvl)%copy(mfv_vel_temp,1,1,1,nghosts_mid)
call mftemp(lvl)%copy(mftemp_temp,1,1,1,nghosts_mid)
!call mfpress(lvl)%copy(mfpress_temp,1,1,1,nghosts_mid)
!
call mfchi(lvl)%copy(mfchi_temp,1,1,1,nghosts_mid)
call mfzeta_x(lvl)%copy(mfzeta_x_temp,1,1,1,nghosts_mid)
call mfzeta_y(lvl)%copy(mfzeta_y_temp,1,1,1,nghosts_mid)
call mfpi_xx(lvl)%copy(mfpi_xx_temp,1,1,1,nghosts_mid)
call mfpi_yy(lvl)%copy(mfpi_yy_temp,1,1,1,nghosts_mid)
call mfpi_xy(lvl)%copy(mfpi_xy_temp,1,1,1,nghosts_mid)
call mflambda_x(lvl)%copy(mflambda_x_temp,1,1,1,nghosts_mid)
call mflambda_y(lvl)%copy(mflambda_y_temp,1,1,1,nghosts_mid)
!
call mfalpha(lvl)%copy(mfalpha_temp,1,1,1,nghosts_mid)
!
call mffi(lvl)%copy(mffi_temp,1,1,dir,nghosts_mid)
call mfgi(lvl)%copy(mfgi_temp,1,1,dir,nghosts_mid)
call mffout(lvl)%copy(mffout_temp,1,1,dir,nghosts_mid)
call mfgout(lvl)%copy(mfgout_temp,1,1,dir,nghosts_mid)
!
if (dimensions == 3) then
  call mfw_vel(lvl)%copy(mfw_vel_temp,1,1,1,nghosts_mid)
  call mfzeta_z(lvl)%copy(mfzeta_z_temp,1,1,1,nghosts_mid)
  call mfpi_xz(lvl)%copy(mfpi_xz_temp,1,1,1,nghosts_mid)
  call mfpi_yz(lvl)%copy(mfpi_yz_temp,1,1,1,nghosts_mid)
  call mfpi_zz(lvl)%copy(mfpi_zz_temp,1,1,1,nghosts_mid)
  call mflambda_z(lvl)%copy(mflambda_z_temp,1,1,1,nghosts_mid)
end if
!
!
!
call amrex_multifab_destroy(mfrho_temp)
call amrex_multifab_destroy(mfu_vel_temp)
call amrex_multifab_destroy(mfv_vel_temp)
call amrex_multifab_destroy(mftemp_temp)
!call amrex_multifab_destroy(mfpress_temp)
!
! LMs
!
!write(*,*) 'clearing lagrange multipliers'
call amrex_multifab_destroy(mfchi_temp)
call amrex_multifab_destroy(mfzeta_x_temp)
call amrex_multifab_destroy(mfzeta_y_temp)
call amrex_multifab_destroy(mfpi_xx_temp)
call amrex_multifab_destroy(mfpi_yy_temp)
call amrex_multifab_destroy(mfpi_xy_temp)
call amrex_multifab_destroy(mflambda_x_temp)
call amrex_multifab_destroy(mflambda_y_temp)
!
! Distro funcs
!
!write(*,*) 'clearing distro funcs'
call amrex_multifab_destroy(mffi_temp)
call amrex_multifab_destroy(mfgi_temp)
call amrex_multifab_destroy(mffout_temp)
call amrex_multifab_destroy(mfgout_temp)
!
! Other
!
call amrex_multifab_destroy(mfalpha_temp)
call amrex_imultifab_destroy(imfstate_temp)
! 3D stuff
if (dimensions == 3) then
  call amrex_multifab_destroy(mfw_vel_temp)
  call amrex_multifab_destroy(mfzeta_z_temp)
  call amrex_multifab_destroy(mfpi_xz_temp)
  call amrex_multifab_destroy(mfpi_yz_temp)
  call amrex_multifab_destroy(mfpi_zz_temp)
  call amrex_multifab_destroy(mflambda_z_temp)
end if

!write(*,*) 'pop, six, squish, unh uh, cicero, Lipschitz',lvl,self
!call update_grid_info(lvl,dummy)
!!call amr_grid_zones(lvl,dummy)
!!
!call fluid_or_solid_from_coarse(lvl)
!!
!call nullify_lvl(lvl,dummy)
!call ghost_lvl(lvl,dummy)
!call coarse_fine_lvl(lvl,dummy)
END SUBROUTINE
!
!
!
!
!
!
!
!
!
!
!
function count_big_zones(lvl) result (n_big_zones)
  type(amrex_mfiter) :: mfi
  !type(amrex_box) = brx
  integer :: lvl,n_big_zones

  n_big_zones = 0

!  write(*,*) 'smuggle',lvl,self
  call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
  do while (mfi%next())
    n_big_zones = n_big_zones + 1
!    call amrex_print()
!    write(*,*) 'found one!', n_big_zones
  end do
!  write(*,*) 'number of zones on each',self,lvl,n_big_zones
  call amrex_mfiter_destroy(mfi)
end function
!
!
! Same as count_big_zones, but includes the MPI function to give every processor the total
!   number of zones at a given level.
!
!
function total_big_zones(lvl) result (lvl_zones)
  use mpi
!
  integer,intent(in) :: lvl
  integer :: lvl_zones,loc_zones,ier

  type(amrex_mfiter) :: mfi

  lvl_zones = 0
  loc_zones = 0

!  write(*,*) 'big zones level ',lvl

  call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
  do while (mfi%next())
    loc_zones = loc_zones +1
  end do
!  write(*,*) 'doodad',self,loc_zones,lvl
  call amrex_mfiter_destroy(mfi)

  CALL MPI_Allreduce(loc_zones,lvl_zones,1,MPI_INTEGER, MPI_SUM,commune,ier)
!  write(*,*) 'local zones',lvl,loc_zones,lvl_zones,self
end function
!
!
!
function get_box_sizes (n_boxes,lvl) result (loc_sizes)
  integer,intent(in) :: n_boxes,lvl
  integer :: loc_sizes(3,n_boxes),counter
  type(amrex_mfiter) :: mfi
  type(amrex_box) :: ox

  call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)

  counter = 1
  do while (mfi%next())
    ox = mfi%validbox()

    loc_sizes(1,counter) = ox%hi(1)-ox%lo(1)+1
    loc_sizes(2,counter) = ox%hi(2)-ox%lo(2)+1
    loc_sizes(3,counter) = ox%hi(3)-ox%lo(3)+1

    counter = counter + 1
  end do

  call amrex_mfiter_destroy(mfi)

end function
!
!
!
function get_real_coords (geo_lvl,loc) result (rxyz)
  class(amrex_geometry),intent(in) :: geo_lvl
  integer,intent(in) :: loc(3)
  integer :: i
  real(kind=dp) :: rxyz(3)
!  write(*,*) 'local dx',geo_lvl%dx(1)
  do i = 1,3
!    write(*,*) 'all grid info',amrex_problo(i),i,loc(i),geo_lvl%domain%lo(i),&
!      geo_lvl%dx(i),self
    rxyz(i) = amrex_problo(i) + (REAL(loc(i),kind=dp) - &
    geo_lvl%domain%lo(i))*geo_lvl%dx(i)
  end do
!  write(*,*) 'get real coords output',rxyz
end function
!
! Determine if a node is inside a box
!
function magic_box (low,high,loc_nd) result (inside)
  integer,intent(in) :: low(3),high(3),loc_nd(3)
!  type(amrex_box), intent(in) :: bucks
  logical :: inside
  inside = all( loc_nd(1:dimensions) >=low(1:dimensions) .and. &
    loc_nd(1:dimensions) <= high(1:dimensions))
end function
!
!
!
function magic_box_exclusive (low,high,loc_nd) result (inside)
  integer,intent(in) :: low(3),high(3),loc_nd(3)
!  type(amrex_box), intent(in) :: bucks
  logical :: inside
  inside = all( loc_nd(1:dimensions) > low(1:dimensions) .and. &
    loc_nd(1:dimensions) < high(1:dimensions))

end function
!
! Determine if a point is inside a box
!
function magic_real_box(lo_corner,hi_corner,loc_pt) result (inside)
  real(kind=dp),intent(in) :: lo_corner(3),hi_corner(3),loc_pt(3)
  logical :: inside

  inside = all(loc_pt(1:dimensions) >= lo_corner(1:dimensions) .and. &
    loc_pt(1:dimensions) <= hi_corner(1:dimensions))
end function
!
!
!
function direction_tester(a,b,c,lox_box,hit_box) result (dir_safe)
!
use linkwise, only: cx,cy,cz
!
integer,intent(in) :: lox_box(dimensions),hit_box(dimensions),a,b,c
logical :: dir_safe(dir)
integer :: i

dir_safe = .false.

do i = 1,dir
  if (a+cx(i) >= lox_box(1) - nghosts_mid .and. b+cy(i) >= lox_box(2) - nghosts_mid &
    .and. c+cz(i) >= lox_box(3)-nghosts_mid .and. a+cx(i) <= hit_box(1)+nghosts_mid &
    .and. b+cy(i) <= hit_box(2)+nghosts_mid .and. c+cz(i) <= hit_box(3)+nghosts_mid) then
!
    dir_safe(i) = .true.
!
  end if
end do

end function
!
!
!
subroutine overlapping_boxes(lo_corner_1,lo_corner_2,hi_corner_1,hi_corner_2,&
  over_lo,over_hi,over)
  integer,intent(in) :: lo_corner_1(3),lo_corner_2(3),hi_corner_1(3),hi_corner_2(3)
  integer,intent(out) :: over_lo(3),over_hi(3)
  logical,intent(out) :: over
  logical :: lo_over(3,2)
  logical :: hi_over(3,2)

  lo_over = .false.
  hi_over = .false.

  if (lo_corner_1(1) >= lo_corner_2(1) .and. lo_corner_1(1) <= hi_corner_2(1)) then
    lo_over(1,1) = .true.
  end if
  if (lo_corner_1(2) >= lo_corner_2(2) .and. lo_corner_1(2) <= hi_corner_2(2)) then
    lo_over(2,1) = .true.
  end if
  if (lo_corner_1(3) >= lo_corner_2(3) .and. lo_corner_1(3) <= hi_corner_2(3)) then
    lo_over(3,1) = .true.
  end if

  if (lo_corner_2(1) >= lo_corner_1(1) .and. lo_corner_2(1) <= hi_corner_1(1)) then
    lo_over(1,2) = .true.
  end if
  if (lo_corner_2(2) >= lo_corner_1(2) .and. lo_corner_2(2) <= hi_corner_1(2)) then
    lo_over(2,2) = .true.
  end if
  if (lo_corner_2(3) >= lo_corner_1(3) .and. lo_corner_2(3) <= hi_corner_1(3)) then
    lo_over(3,2) = .true.
  end if

  if (hi_corner_1(1) >= lo_corner_2(1) .and. hi_corner_1(1) <= hi_corner_2(1)) then
    hi_over(1,1) = .true.
  end if
  if (hi_corner_1(2) >= lo_corner_2(2) .and. hi_corner_1(2) <= hi_corner_2(2)) then
    hi_over(2,1) = .true.
  end if
  if (hi_corner_1(3) >= lo_corner_2(3) .and. hi_corner_1(3) <= hi_corner_2(3)) then
    hi_over(3,1) = .true.
  end if

  if (hi_corner_2(1) >= lo_corner_1(1) .and. hi_corner_2(1) <= hi_corner_1(1)) then
    hi_over(1,2) = .true.
  end if
  if (lo_corner_2(2) >= lo_corner_1(2) .and. hi_corner_2(2) <= hi_corner_1(2)) then
    hi_over(2,2) = .true.
  end if
  if (hi_corner_2(3) >= lo_corner_1(3) .and. hi_corner_2(3) <= hi_corner_1(3)) then
    hi_over(3,2) = .true.
  end if

  if (any(lo_over(1,1:2)) .and. any(lo_over(2,1:2)) .and. any(lo_over(3,1:2))) then
    if (lo_over(1,1) .and. .not. lo_over(1,2)) then
      over_lo(1) = lo_corner_1(1)
    else if (lo_over(1,2) .and. .not. lo_over(1,1)) then
      over_lo(1) = lo_corner_2(1)
    else if (lo_over(1,1) .and.  lo_over(1,2)) then
      over_lo(1) = lo_corner_1(1)
    end if

    if (lo_over(2,1) .and. .not. lo_over(2,2)) then
      over_lo(2) = lo_corner_1(2)
    else if (lo_over(2,2) .and. .not. lo_over(2,1)) then
      over_lo(2) = lo_corner_2(2)
    else if (lo_over(2,1) .and.  lo_over(2,2)) then
      over_lo(2) = lo_corner_1(2)
    end if

    if (lo_over(3,1) .and. .not. lo_over(3,2)) then
      over_lo(3) = lo_corner_1(3)
    else if (lo_over(1,2) .and. .not. lo_over(3,1)) then
      over_lo(3) = lo_corner_2(3)
    else if (lo_over(3,1) .and.  lo_over(3,2)) then
      over_lo(3) = lo_corner_1(3)
    end if

    if (hi_over(1,1) .and. .not. hi_over(1,2)) then
      over_hi(1) = hi_corner_2(1)
    else if (hi_over(1,2) .and. .not. hi_over(1,1)) then
      over_hi(1) = hi_corner_1(1)
    else if (hi_over(1,1) .and.  hi_over(1,2)) then
      over_hi(1) = hi_corner_1(1)
    end if

    if (hi_over(2,1) .and. .not. hi_over(2,2)) then
      over_hi(2) = hi_corner_2(2)
    else if (hi_over(2,2) .and. .not. hi_over(2,1)) then
      over_hi(2) = hi_corner_1(2)
    else if (hi_over(2,1) .and.  hi_over(2,2)) then
      over_hi(2) = hi_corner_1(2)
    end if

    if (hi_over(3,1) .and. .not. hi_over(3,2)) then
      over_hi(3) = hi_corner_2(3)
    else if (hi_over(3,2) .and. .not. hi_over(3,1)) then
      over_hi(3) = hi_corner_1(3)
    else if (hi_over(3,1) .and.  hi_over(3,2)) then
      over_hi(3) = hi_corner_1(3)
    end if

    over = .true.
  else
   over_lo = -1
   over_hi = -1

    over = .false.
  end if

end subroutine
!
!
!
subroutine ml_overlapping_boxes(lo_corner_in_1,lo_corner_in_2,hi_corner_in_1,hi_corner_in_2,&
  over_lo,over_hi,over,lvl_1,lvl_2)
  integer,intent(IN) :: lo_corner_in_1(3),lo_corner_in_2(3),hi_corner_in_1(3),hi_corner_in_2(3),&
    lvl_1,lvl_2
  integer,intent(out) :: over_lo(3),over_hi(3)
  integer :: lo_corner_1(3),lo_corner_2(3),hi_corner_2(3),hi_corner_1(3),lvl_diff
  logical,intent(out) :: over
  logical :: lo_over(3,2)
  logical :: hi_over(3,2)

  lo_over = .false.
  hi_over = .false.

  if (lvl_1 == lvl_2) then
    call overlapping_boxes(lo_corner_in_1,lo_corner_in_2,hi_corner_in_1,hi_corner_in_2,&
      over_lo,over_hi,over)
    return
  end if

  if (lvl_1 > lvl_2) then
    lvl_diff = lvl_1-lvl_2
    lo_corner_2 = lo_corner_in_2*(lvl_diff*2)
    hi_corner_2 = hi_corner_in_2*(lvl_diff*2)
    lo_corner_1 = lo_corner_in_1
    hi_corner_1 = hi_corner_in_1
!    write(*,*) 'new corners',lo_corner_2,hi_corner_2
!    write(*,*) 'old corners',lo_corner_1,hi_corner_1
  else if (lvl_2 > lvl_1) then
    lvl_diff = lvl_2-lvl_1
    lo_corner_1 = lo_corner_in_1*(lvl_diff*2)
    hi_corner_1 = hi_corner_in_1*(lvl_diff*2)
    lo_corner_2 = lo_corner_in_2
    hi_corner_2 = hi_corner_in_2
!    write(*,*) 'new corners',lo_corner_1,hi_corner_1
!    write(*,*) 'old corners',lo_corner_2,hi_corner_2
  end if
!
! Determine what is inside of where
!
  if (lo_corner_1(1) >= lo_corner_2(1) .and. lo_corner_1(1) <= hi_corner_2(1)) then
    lo_over(1,1) = .true.
  end if
  if (lo_corner_1(2) >= lo_corner_2(2) .and. lo_corner_1(2) <= hi_corner_2(2)) then
    lo_over(2,1) = .true.
  end if
  if (lo_corner_1(3) >= lo_corner_2(3) .and. lo_corner_1(3) <= hi_corner_2(3)) then
    lo_over(3,1) = .true.
  end if

  if (lo_corner_2(1) >= lo_corner_1(1) .and. lo_corner_2(1) <= hi_corner_1(1)) then
    lo_over(1,2) = .true.
  end if
  if (lo_corner_2(2) >= lo_corner_1(2) .and. lo_corner_2(2) <= hi_corner_1(2)) then
    lo_over(2,2) = .true.
  end if
  if (lo_corner_2(3) >= lo_corner_1(3) .and. lo_corner_2(3) <= hi_corner_1(3)) then
    lo_over(3,2) = .true.
  end if

  if (hi_corner_1(1) >= lo_corner_2(1) .and. hi_corner_1(1) <= hi_corner_2(1)) then
    hi_over(1,1) = .true.
  end if
  if (hi_corner_1(2) >= lo_corner_2(2) .and. hi_corner_1(2) <= hi_corner_2(2)) then
    hi_over(2,1) = .true.
  end if
  if (hi_corner_1(3) >= lo_corner_2(3) .and. hi_corner_1(3) <= hi_corner_2(3)) then
    hi_over(3,1) = .true.
  end if

  if (hi_corner_2(1) >= lo_corner_1(1) .and. hi_corner_2(1) <= hi_corner_1(1)) then
    hi_over(1,2) = .true.
  end if
  if (hi_corner_2(2) >= lo_corner_1(2) .and. hi_corner_2(2) <= hi_corner_1(2)) then
    hi_over(2,2) = .true.
  end if
  if (hi_corner_2(3) >= lo_corner_1(3) .and. hi_corner_2(3) <= hi_corner_1(3)) then
    hi_over(3,2) = .true.
  end if
!
! Ensure there is overlap, then assign corners for the overlapping region
!
  if (any(lo_over(1,1:2)) .and. any(lo_over(2,1:2)) .and. any(lo_over(3,1:2))) then
!    write(*,*) 'overlap exists'
    if (lo_over(1,1) .and. .not. lo_over(1,2)) then
      over_lo(1) = lo_corner_1(1)
    else if (lo_over(1,2) .and. .not. lo_over(1,1)) then
      over_lo(1) = lo_corner_2(1)
    else if (lo_over(1,1) .and.  lo_over(1,2)) then
      over_lo(1) = lo_corner_1(1)
    end if

    if (lo_over(2,1) .and. .not. lo_over(2,2)) then
      over_lo(2) = lo_corner_1(2)
    else if (lo_over(2,2) .and. .not. lo_over(2,1)) then
      over_lo(2) = lo_corner_2(2)
    else if (lo_over(2,1) .and.  lo_over(2,2)) then
      over_lo(2) = lo_corner_1(2)
    end if

    if (lo_over(3,1) .and. .not. lo_over(3,2)) then
      over_lo(3) = lo_corner_1(3)
    else if (lo_over(3,2) .and. .not. lo_over(3,1)) then
      over_lo(3) = lo_corner_2(3)
    else if (lo_over(3,1) .and.  lo_over(3,2)) then
      over_lo(3) = lo_corner_1(3)
    end if

    if (hi_over(1,1) .and. .not. hi_over(1,2)) then
      over_hi(1) = hi_corner_1(1)
    else if (hi_over(1,2) .and. .not. hi_over(1,1)) then
      over_hi(1) = hi_corner_2(1)
    else if (hi_over(1,1) .and.  hi_over(1,2)) then
      over_hi(1) = hi_corner_2(1)
    end if

    if (hi_over(2,1) .and. .not. hi_over(2,2)) then
      over_hi(2) = hi_corner_1(2)
    else if (hi_over(2,2) .and. .not. hi_over(2,1)) then
      over_hi(2) = hi_corner_2(2)
    else if (hi_over(2,1) .and.  hi_over(2,2)) then
      over_hi(2) = hi_corner_2(2)
    end if

    if (hi_over(3,1) .and. .not. hi_over(3,2)) then
      over_hi(3) = hi_corner_1(3)
    else if (hi_over(3,2) .and. .not. hi_over(3,1)) then
      over_hi(3) = hi_corner_2(3)
    else if (hi_over(3,1) .and.  hi_over(3,2)) then
      over_hi(3) = hi_corner_2(3)
    end if

    over = .true.
    over_lo = over_lo/(lvl_diff*2)
    over_hi = over_hi/(lvl_diff*2)
!    write(*,*) 'overlapping region',over_lo,over_hi
  else
   over_lo = -1
   over_hi = -1

    over = .false.
  end if

end subroutine
!
!
!
!
!
subroutine mpi_amrex_retrieve
  integer :: ier
!  write(*,*) 'doodles...'
  self = amrex_parallel_myproc()
  nprocs = amrex_parallel_nprocs()
  commune = amrex_parallel_communicator()
  amr_max_lvl = amrex_max_level
  finest_lvl = amrex_get_finest_level()
end subroutine
!
!
!
!
!
function fieq_comp_19 (loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z,loc_dir) result (loc_fieq)
  integer :: loc_dir
  real(kind=dp) :: loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z
  real(kind=dp) :: loc_fieq

select case (loc_dir)
  case(1)
    loc_fieq = loc_rho*EXP(loc_chi)
  case(2)
    loc_fieq = loc_rho*EXP(loc_chi + loc_zeta_x + loc_pi_xx + loc_lambda_x)
  case(3)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_y + loc_pi_yy + loc_lambda_y)
  case(4)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_z + loc_pi_zz + loc_lambda_z)
  case(5)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_x + loc_pi_xx - loc_lambda_x)
  case(6)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_y + loc_pi_yy - loc_lambda_y)
  case(7)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_z + loc_pi_zz - loc_lambda_z)
  case(8)
! 8 -  (+1,+1,0)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_x + loc_zeta_y + &
      loc_pi_xx + loc_pi_yy + loc_pi_xy + &
      2.0D0*loc_lambda_x + 2.0D0*loc_lambda_y)
  case(9)
! 9 -  (0,+1,-1)
    loc_fieq = loc_rho*EXP(loc_chi + loc_zeta_y - loc_zeta_z + &
      loc_pi_yy + loc_pi_zz - loc_pi_yz + &
      2.0D0*loc_lambda_y - 2.0D0*loc_lambda_z)
  case(10)
! 10 - (-1,+1,0)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_x + loc_zeta_y + &
      loc_pi_xx + loc_pi_yy - loc_pi_xy - &
      2.0D0*loc_lambda_x + 2.0D0*loc_lambda_y)
  case(11)
! 11 - (0,+1,+1)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_y + loc_zeta_z + &
      loc_pi_yy + loc_pi_zz + loc_pi_yz +&
      2.0D0*loc_lambda_y + 2.0D0*loc_lambda_z)
  case(12)
! 12 - (+1,0,+1)
    loc_fieq = loc_rho*EXP(loc_chi + loc_zeta_x + loc_zeta_z + &
      loc_pi_xx + loc_pi_zz + loc_pi_xz +&
      2.0D0*loc_lambda_x + 2.0D0*loc_lambda_z)
  case(13)
! 13 - (+1,0,-1)
    loc_fieq = loc_rho*EXP(loc_chi + loc_zeta_x - loc_zeta_z + &
      loc_pi_xx + loc_pi_zz - loc_pi_xz +&
      2.0D0*loc_lambda_x - 2.0D0*loc_lambda_z)
  case(14)
! 14 - (-1,0,-1)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_x - loc_zeta_z + &
      loc_pi_xx + loc_pi_zz + loc_pi_xz - &
      2.0D0*loc_lambda_x - 2.0D0*loc_lambda_z)
  case(15)
! 15 - (-1,0,+1)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_x + loc_zeta_z + &
      loc_pi_xx + loc_pi_zz - loc_pi_xz - &
      2.0D0*loc_lambda_x + 2.0D0*loc_lambda_z)
  case(16)
! 16 - (+1,-1,0)
    loc_fieq = loc_rho*EXP(loc_chi + loc_zeta_x - loc_zeta_y + &
      loc_pi_xx + loc_pi_yy - loc_pi_xy +&
      2.0D0*loc_lambda_x - 2.0D0*loc_lambda_y)
  case(17)
! 17 - (0,-1,-1)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_y - loc_zeta_z + &
      loc_pi_yy + loc_pi_zz + loc_pi_yz - &
      2.0D0*loc_lambda_y - 2.0D0*loc_lambda_z)
  case(18)
! 18 - (-1,-1,0)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_x - loc_zeta_y + &
      loc_pi_xx + loc_pi_yy + loc_pi_xy  - &
      2.0D0*loc_lambda_x - 2.0D0*loc_lambda_y)
  case(19)
! 19 - (0,-1,+1)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_y + loc_zeta_z + &
       loc_pi_yy + loc_pi_zz - loc_pi_yz - &
       2.0D0*loc_lambda_y + 2.0D0*loc_lambda_z)
end select
end function
!
!
!
!
!
function fieq_comp_39 (loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z,loc_dir) result (loc_fieq)
  integer :: loc_dir
  real(kind=dp) :: loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z
  real(kind=dp) :: loc_rho
  real(kind=dp) :: loc_fieq

  select case (loc_dir)
  case(1)
    loc_fieq = loc_rho*EXP(loc_chi)
  case(2)
    loc_fieq = loc_rho*EXP(loc_chi + loc_zeta_x + loc_pi_xx + loc_lambda_x)
  case(3)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_y + loc_pi_yy + loc_lambda_y)
  case(4)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_z + loc_pi_zz + loc_lambda_z)
  case(5)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_x + loc_pi_xx - loc_lambda_x)
  case(6)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_y + loc_pi_yy - loc_lambda_y)
  case(7)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_z + loc_pi_zz - loc_lambda_z)
  case(8)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_x + loc_zeta_y + loc_zeta_z + &
      loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy + loc_pi_xz + loc_pi_yz +&
      3.0D0*loc_lambda_x + 3.0D0*loc_lambda_y + 3.0D0*loc_lambda_z)
  case(9)
    loc_fieq = loc_rho*EXP(loc_chi - loc_zeta_x + loc_zeta_y + loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy - loc_pi_xz + loc_pi_yz +&
    3.0D0*(-loc_lambda_x + loc_lambda_y + loc_lambda_z))
  case(10)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_x - loc_zeta_y + loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy - loc_pi_xz - loc_pi_yz -&
    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y + 3.0D0*loc_lambda_z)
  case(11)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_x - loc_zeta_y + loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy + loc_pi_xz - loc_pi_yz +&
    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y + 3.0D0*loc_lambda_z)
  case(12)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_x + loc_zeta_y - loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy - loc_pi_xz - loc_pi_yz + &
    3.0D0*loc_lambda_x + 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
  case(13)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_x + loc_zeta_y - loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy + loc_pi_xz - loc_pi_yz - &
    3.0D0*loc_lambda_x + 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
  case(14)
    loc_fieq =loc_rho*EXP(loc_chi - loc_zeta_x - loc_zeta_y - loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz + loc_pi_xy + loc_pi_xz + loc_pi_yz -&
    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
  case(15)
    loc_fieq =loc_rho*EXP(loc_chi + loc_zeta_x - loc_zeta_y - loc_zeta_z + &
    loc_pi_xx + loc_pi_yy + loc_pi_zz - loc_pi_xy - loc_pi_xz + loc_pi_yz +&
    3.0D0*loc_lambda_x - 3.0D0*loc_lambda_y - 3.0D0*loc_lambda_z)
  case(16)
    loc_fieq = loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x + 4.0D0*loc_pi_xx + 8.0D0*loc_lambda_x)
  case(17)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_y + 4.0D0*loc_pi_yy + 8.0D0*loc_lambda_y)
  case(18)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_z + 4.0D0*loc_pi_zz + 8.0D0*loc_lambda_z)
  case(19)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x + 4.0D0*loc_pi_xx - 8.0D0*loc_lambda_x)
  case(20)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_y + 4.0D0*loc_pi_yy - 8.0D0*loc_lambda_y)
  case(21)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_z + 4.0D0*loc_pi_zz - 8.0D0*loc_lambda_z)
  case(22)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_y + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy + 4.0D0*loc_pi_xy + &
    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_y)
  case(23)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_z + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_xz +&
    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_z)
  case(24)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_y + 2.0D0*loc_zeta_z + &
    4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_yz +&
    16.0D0*loc_lambda_y + 16.0D0*loc_lambda_z)
  case(25)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_y + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy - 4.0D0*loc_pi_xy - &
    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_y)
  case(26)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_y + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy + 4.0D0*loc_pi_xy  - &
    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_y)
  case(27)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_y + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_yy - 4.0D0*loc_pi_xy +&
    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_y)
  case(28)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_z + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_xz +&
    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_z)
  case(29)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x - 2.0D0*loc_zeta_z + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_xz - &
    16.0D0*loc_lambda_x - 16.0D0*loc_lambda_z)
  case(30)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_x + 2.0D0*loc_zeta_z + &
    4.0D0*loc_pi_xx + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_xz - &
    16.0D0*loc_lambda_x + 16.0D0*loc_lambda_z)
  case(31)
    loc_fieq =loc_rho*EXP(loc_chi + 2.0D0*loc_zeta_y - 2.0D0*loc_zeta_z + &
     4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_yz + &
     16.0D0*loc_lambda_y - 16.0D0*loc_lambda_z)
  case(32)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_y - 2.0D0*loc_zeta_z + &
    4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz + 4.0D0*loc_pi_yz - &
    16.0D0*loc_lambda_y - 16.0D0*loc_lambda_z)
  case(33)
    loc_fieq =loc_rho*EXP(loc_chi - 2.0D0*loc_zeta_y + 2.0D0*loc_zeta_z + &
     4.0D0*loc_pi_yy + 4.0D0*loc_pi_zz - 4.0D0*loc_pi_yz - &
    16.0D0*loc_lambda_y + 16.0D0*loc_lambda_z)
  case(34)
    loc_fieq =loc_rho*EXP(loc_chi + 3.0D0*loc_zeta_x + 9.0D0*loc_pi_xx + 27.0D0*loc_lambda_x)
  case(35)
    loc_fieq =loc_rho*EXP(loc_chi + 3.0D0*loc_zeta_y + 9.0D0*loc_pi_yy + 27.0D0*loc_lambda_y)
  case(36)
    loc_fieq =loc_rho*EXP(loc_chi + 3.0D0*loc_zeta_z + 9.0D0*loc_pi_zz + 27.0D0*loc_lambda_z)
  case(37)
    loc_fieq =loc_rho*EXP(loc_chi - 3.0D0*loc_zeta_x + 9.0D0*loc_pi_xx - 27.0D0*loc_lambda_x)
  case(38)
    loc_fieq =loc_rho*EXP(loc_chi - 3.0D0*loc_zeta_y + 9.0D0*loc_pi_yy - 27.0D0*loc_lambda_y)
  case(39)
    loc_fieq =loc_rho*EXP(loc_chi - 3.0D0*loc_zeta_z + 9.0D0*loc_pi_zz - 27.0D0*loc_lambda_z)
  end select

end function
!
!
!
!
!
function fieq_comp_39_shifted (loc_rho,loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z,loc_dir) result (loc_fieq)

use linkwise

  integer :: loc_dir
  real(kind=dp) :: loc_chi,loc_zeta_x,loc_zeta_y,loc_zeta_z,loc_pi_xx,&
    loc_pi_yy,loc_pi_zz,loc_pi_xy,loc_pi_xz,loc_pi_yz,loc_lambda_x,loc_lambda_y,&
    loc_lambda_z
  real(kind=dp) :: loc_rho
  real(kind=dp) :: loc_fieq

  loc_fieq = loc_rho*exp(loc_chi + cx(loc_dir)*loc_zeta_x + cy(loc_dir)*loc_zeta_y + cz(loc_dir)*loc_zeta_z +&
    cx2(loc_dir)*loc_pi_xx + cy2(loc_dir)*loc_pi_yy + cz2(loc_dir)*loc_pi_zz + cxcy(loc_dir)*loc_pi_xy +&
    cxcz(loc_dir)*loc_pi_xz + cycz(loc_dir)*loc_pi_yz + cxc2(loc_dir)*loc_lambda_x + &
    cyc2(loc_dir)*loc_lambda_y + czc2(loc_dir)*loc_lambda_z)

end function
!
!
!
function beta_1_calculator (loc_mu,loc_rho,loc_temp,lvl,steptime) result (beta_1)
use timing
!
!
!
integer,intent(in) :: steptime,lvl
real(kind=dp) :: beta_1
real(kind=dp),intent(in) :: loc_rho,loc_temp,loc_mu

if (steptime < startup_timesteps) then
  beta_1 = init_beta_1 + beta_1_scaling*steptime
else
  beta_1 = 1.0D0/((2.0D0**(lvl+1)*loc_mu)/(loc_rho*loc_temp)+1.0D0)
end if
!beta_1 = 1.0D0/(2.0D0**(lvl+1)*loc_mu/(loc_rho*loc_temp)+1.0D0)
!beta_1 = 0.75D0
end function

function beta_2_calculator (loc_kappa,loc_rho,loc_temp,lvl,steptime) result (beta_2)
use timing
!
!
!
integer,intent(in) :: steptime,lvl
real(kind=dp) :: beta_2
real(kind=dp),intent(in) :: loc_rho,loc_temp,loc_kappa

if (steptime < startup_timesteps) then
  beta_2 = init_beta_2 + beta_2_scaling*steptime
else
  beta_2 = 1.0D0/((2.0D0*loc_kappa)/(loc_rho*const_press*loc_temp)+1.0D0)
end if
!beta_2 = 1.0D0/(2.0D0**(lvl+1)*loc_kappa/(loc_rho*const_press*loc_temp)+1.0D0)
!beta_2 = 0.75D0
end function
!
!
!
function epsilon_calculator(loc_u,loc_v,loc_w,loc_temp,steptime,a,b,c,lvl) result (loc_epsilon)
use timing
!
!
!
real(kind=dp),intent(in) :: loc_u,loc_v,loc_w,loc_temp
real(kind=dp) :: loc_epsilon,v_mag,unit_velocity
integer,intent(in) :: steptime,lvl,a,b,c

!v_mag = sqrt(loc_u**2 + loc_v**2 + loc_w**2)
!
!unit_velocity = zoom_zoom_calculator(v_mag,loc_temp)
!
!if (unit_velocity < 1.0D0) unit_velocity = 1.0D0

if (steptime < startup_timesteps) then
!  loc_epsilon = ((real(steptime,kind=dp)+1.0D0)/real(startup_timesteps,kind=dp))*&
!    (elbm_epsilon(lvl)/(1.0D0/unit_velocity))
  loc_epsilon = ((real(steptime,kind=dp)+1.0D0)/real(startup_timesteps,kind=dp))*&
    elbm_epsilon(lvl)
else
!  loc_epsilon = elbm_epsilon(lvl)/(1.0D0/unit_velocity)
  loc_epsilon = elbm_epsilon(lvl)
end if
!write(*,*) 'local epsilon',loc_u,loc_v,loc_w,loc_temp,steptime,a,b,c,lvl,startup_timesteps,loc_epsilon
end function
!
!
!
function zoom_zoom_calculator (v_mag,loc_temp) result (unit_velocity)
!
!
!
real(kind=dp),intent(in) :: v_mag,loc_temp
real(kind=dp) :: unit_velocity,unit_temp,loc_sos

unit_temp = loc_temp*t_ref
loc_sos = sqrt(gama*gas_const_r*unit_temp)

unit_velocity = v_mag*loc_sos

end function
!
!
!
function fill_qone_interlevel(dTdx,dTdy,dTdz) result (qabc)
real(kind=dp) :: qabc(27),dTdx,dTdy,dTdz

qabc(1) = 3*dTdx
qabc(2) = dTdy
qabc(3) = dTdz
qabc(4) = dTdy
qabc(5) = dTdx
qabc(6) = 0.0D0
qabc(7) = dTdz
qabc(8) = 0.0D0
qabc(9) = dTdx
qabc(10) = dTdy
qabc(11) = dTdx
qabc(12) = 0.0D0
qabc(13) = dTdx
qabc(14) = 3*dTdy
qabc(15) = dTdz
qabc(16) = 0.0D0
qabc(17) = dTdz
qabc(18) = dTdy
qabc(19) = dTdz
qabc(20) = 0.0D0
qabc(21) = dTdx
qabc(22) = 0.0D0
qabc(23) = dTdz
qabc(24) = dTdy
qabc(25) = dTdx
qabc(26) = dTdy
qabc(27) = 3*dTdz

end function
!
!
!
!
!
function zeta_x_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_zeta_x)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_zeta_x

new_zeta_x = 4.68D0 - 1.54D0*loc_rho + 0.285D0*loc_rho**2 - 2.17D0*loc_u + 2.19D0*loc_u**2 - &
  0.1D0*loc_v + loc_v**2 - 0.2D0*loc_w + loc_w**2 - 9.53D0*loc_temp + 6.66D0*loc_temp**2

end function
!
!
!
function chi_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_chi)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_chi

new_chi = -9.1D0 + 0.348D0*loc_rho - 0.052D0*loc_rho**2 + 1.684D0*loc_u - 1.847D0*loc_u**2 + &
  0.237D0*loc_v - 1.2D0*loc_v**2 + 0.237D0*loc_w - 1.2D0*loc_w**2 + &
  17.07D0*loc_temp - 10.73D0*loc_temp**2

end function
!
!
!
function zeta_y_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_zeta_y)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_zeta_y

if (abs(loc_v) <= sp_plus_epsilon) then
  new_zeta_y = 0.0D0
else
new_zeta_y = 0.69D0 + 0.525D0*loc_rho - 0.065D0*loc_rho**2 + 0.552D0*loc_u - 0.31D0*loc_u**2 + &
  1.45D0*loc_v - 0.199D0*loc_v**2 + 0.184D0*loc_w - 0.216D0*loc_w**2 - &
  2.38D0*loc_temp + 0.816D0*loc_temp**2
end if

end function
!
!
!
function zeta_z_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_zeta_z)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_zeta_z

if (abs(loc_w) <= sp_plus_epsilon) then
  new_zeta_z = 0.0D0
else
new_zeta_z = -0.73D0 - 0.405D0*loc_rho + 0.0504D0*loc_rho**2 - 0.264D0*loc_u + 0.078D0*loc_u**2 + &
  0.353D0*loc_v - 0.144D0*loc_v**2 + 1.29D0*loc_w + 0.262D0*loc_w**2 + &
  2.52D0*loc_temp - 1.098D0*loc_temp**2
end if

end function
!
!
!
function pi_xx_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_pi_xx)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_pi_xx

new_pi_xx = 1.73D0 - 0.705D0*loc_rho + 0.249D0*loc_rho**2 + 1.38D0*loc_u - 0.8898D0*loc_u**2 - &
  4.52D0*loc_temp + 2.586D0*loc_temp**2

end function
!
!
!
function pi_yy_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_pi_yy)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_pi_yy

new_pi_yy = 1.501D0 + 0.247D0*loc_rho - 0.0199D0*loc_rho**2 + 0.519D0*loc_u - 0.352D0*loc_u**2 - &
  0.00666D0*loc_v - 0.638D0*loc_v**2 - 0.0137D0*loc_w - 0.192D0*loc_w**2 - &
  4.76D0*loc_temp + 2.552D0*loc_temp**2

end function
!
!
!
function pi_zz_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_pi_zz)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_pi_zz

new_pi_zz = 0.8022D0 + 0.0269D0*loc_rho + 0.0553D0*loc_rho**2 + 0.728D0*loc_u - 0.406D0*loc_u**2 - &
  0.02342D0*loc_v - 0.485D0*loc_v**2 - 0.0823D0*loc_w - 0.3951D0*loc_w**2 - &
  3.334D0*loc_temp + 1.896D0*loc_temp**2

end function
!
!
!
function pi_xy_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_pi_xy)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_pi_xy

if (abs(loc_u) <= sp_plus_epsilon .or. abs(loc_v) <= sp_plus_epsilon) then
  new_pi_xy = 0.0D0
else
new_pi_xy = 0.8022D0 + 0.0269D0*loc_rho + 0.0553D0*loc_rho**2 + 0.728D0*loc_u - 0.406D0*loc_u**2 - &
  0.02342D0*loc_v - 0.485D0*loc_v**2 - 0.0823D0*loc_w - 0.3951D0*loc_w**2 - &
  3.334D0*loc_temp + 1.896D0*loc_temp**2
end if

end function
!
!
!
function pi_xz_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_pi_xz)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_pi_xz

new_pi_xz = 0.0D0

!0.8022D0 + 0.0269D0*loc_rho + 0.0553D0*loc_rho**2 + 0.728D0*loc_u - 0.406D0*loc_u**2 - &
!  0.02342D0*loc_v - 0.485D0*loc_v**2 - 0.0823D0*loc_w - 0.3951D0*loc_w**2 - &
!  3.334D0*loc_temp + 1.896D0*loc_temp**2

end function
!
!
!
function pi_yz_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_pi_yz)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_pi_yz

new_pi_yz = 0.0D0

!0.8022D0 + 0.0269D0*loc_rho + 0.0553D0*loc_rho**2 + 0.728D0*loc_u - 0.406D0*loc_u**2 - &
!  0.02342D0*loc_v - 0.485D0*loc_v**2 - 0.0823D0*loc_w - 0.3951D0*loc_w**2 - &
!  3.334D0*loc_temp + 1.896D0*loc_temp**2

end function
!
!
!
function lambda_x_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_lambda_x)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_lambda_x

new_lambda_x = -1.942D0 + 0.0716D0*loc_rho - 0.0368D0*loc_rho**2 - 0.0956D0*loc_u + 0.107D0*loc_u**2 - &
  0.0174D0*loc_v + 0.1096D0*loc_v**2 + 0.0269D0*loc_w - 0.03079D0*loc_w**2 + &
  4.091D0*loc_temp - 2.161D0*loc_temp**2

end function
!
!
!
function lambda_y_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_lambda_y)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_lambda_y

new_lambda_y = 0.0D0
!-1.942D0 + 0.0716D0*loc_rho - 0.0368D0*loc_rho**2 - 0.0956D0*loc_u + 0.107D0*loc_u**2 - &
!  0.0174D0*loc_v + 0.1096D0*loc_v**2 + 0.0269D0*loc_w - 0.03079D0*loc_w**2 + &
!  4.091D0*loc_temp - 2.161D0*loc_temp**2

end function
!
!
!
function lambda_z_estimate(loc_rho,loc_u,loc_v,loc_w,loc_temp) result (new_lambda_z)
!
real(kind=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_temp,new_lambda_z

new_lambda_z = 0.0D0
!-1.942D0 + 0.0716D0*loc_rho - 0.0368D0*loc_rho**2 - 0.0956D0*loc_u + 0.107D0*loc_u**2 - &
!  0.0174D0*loc_v + 0.1096D0*loc_v**2 + 0.0269D0*loc_w - 0.03079D0*loc_w**2 + &
!  4.091D0*loc_temp - 2.161D0*loc_temp**2

end function



end module








