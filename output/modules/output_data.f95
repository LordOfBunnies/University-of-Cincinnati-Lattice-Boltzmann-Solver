MODULE output_data
!
! This is to hold the data used to output the zones to the CGNS file.
!
!
!
!
! output_info data
! 1) information saved
! 2) save number
! 3) direction link_number
! 4) link_number
!
! output_grid data
! 1) x,y,z
! 2) link number
! 3) save number
!
! output_planes_info
! 1) information saved
! 2) save number
! 3) direction link_number
! 4) link_number
!
! output_planes_grid
! 1) x,y,z
! 2) link number
! 3) save number
!
USE precise
use nml_output
IMPLICIT NONE
SAVE
private
public :: x_verts,y_verts,z_verts,box_zones,mic_weights
public :: save_steps,ggm_steps,base_info,ggm_base_info

real(kind=dp),allocatable :: x_verts(:,:,:),y_verts(:,:,:),&
  z_verts(:,:,:),mic_weights(:,:)

type(timing_info_steps),allocatable :: save_steps(:),ggm_steps(:)

type(base_data),allocatable :: base_info(:),ggm_base_info(:)
!integer,allocatable :: n_outsteps(:),step_lvl_saves(:)

type,public :: box_zones
  integer,allocatable :: lower(:,:),higher(:,:),little_zones(:)
!  integer,allocatable :: ghost_lo(:,:),ghost_hi(:,:)
  real(kind=dp),allocatable :: lo_corner(:,:),hi_corner(:,:)
  integer :: big_zones
  integer,allocatable :: output_number(:,:),zone_start(:),zone_end(:)!,zone_lvl(:)

  integer,allocatable :: box_owner(:),nodes_in_box(:)

  type(multigrid_save_info),allocatable :: multigrid_save(:),ggm_save_data(:)

end type



END MODULE


!  real(kind=dp),allocatable,pointer :: win_start_u_vel(:),win_start_v_vel(:),win_start_w_vel(:),&
!    win_start_temp(:),win_start_rho(:),win_start_fi(:,:)
!  real(kind=dp),allocatable,pointer :: win_start_chi(:),win_start_zeta_x(:),win_start_zeta_y(:),&
!    win_start_zeta_z(:),win_start_pi_xx(:),win_start_pi_yy(:),win_start_pi_xy(:),win_start_pi_zz(:),&
!    win_start_pi_xz(:),win_start_pi_yz(:),win_start_lambda_x(:),win_start_lambda_y(:),&
!    win_start_lambda_z(:)
