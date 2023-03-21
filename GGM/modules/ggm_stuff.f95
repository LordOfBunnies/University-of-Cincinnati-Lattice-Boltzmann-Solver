MODULE ggm_stuff
!
! Stores data for doing GGM stuff
!
!
!
!
USE precise
use amrex_base_module
use amrex_amr_module
IMPLICIT NONE
SAVE
real(kind=dp),contiguous,pointer :: gm(:,:,:,:),ggm(:,:,:,:)
real(kind=dp),allocatable :: holder(:,:,:,:),ggm_limit(:,:)
integer :: ggm_order,num_ggm_vars,ggm_dir,num_regrids,regrid_basis
integer,allocatable :: ggm_variables(:),ggm_max_lvl(:),ggm_min_lvl(:)!,regrid_timesteps(:,:)
real(kind=dp),allocatable :: ggm_output_freq(:),ggm_tag_limit(:)!,regrid_times(:)
real(kind=dp) :: ggm_dt,ggm_sensitivity
logical :: upwind,downwind,central,ggm_output,use_ggm

!
! PGS variables
!
!integer,allocatable :: pgs_tags(:,:,:,:)
integer :: pgs_max_scan
logical :: use_pgs

!
! For pure gradient comparison
!
real(kind=dp) :: gradient_limit
real(kind=dp) :: ggm_time_total,regrid_time_total
integer :: gradient_selection
logical :: gradient_refine_only

type(amrex_multifab),allocatable :: mfgm(:),mfggm(:)



END MODULE

!REAL(KIND = dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!REAL(KIND = dp) :: drhodx,drhody,drhodz,dpdx,dpdy,dpdz
!REAL(KIND = dp) :: dqdx,dqdy,dqdz
!REAL(KIND = dp) :: split_limit,join_limit
!REAL(KIND = dp),ALLOCATABLE :: ggm(:),gm(:),q(:),q2(:)
!INTEGER :: ggm_order, num_ggm_vars, ggm_dir,ggm_output_freq
!INTEGER,ALLOCATABLE :: ggm_variables(:)
!LOGICAL :: use_ggm
