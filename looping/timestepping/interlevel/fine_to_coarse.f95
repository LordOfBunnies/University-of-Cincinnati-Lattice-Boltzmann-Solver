subroutine fine_to_coarse(lvl,fine)
!
! Average down to the lower level
!
!
! Called by: Calculate_flow_variables
! Calls:
!
use precise
use constants
use amr_info_holder
use grid_data
use amr_processes, only: self
use amrex_base_module
use amrex_amr_module
implicit none
integer :: lvl,fine
!
! Can't average down on lvl 0
!
if (self == 0) write(*,*) 'averaging down',lvl
!
!
!
if (lvl == 0) return
!
!
if (mod(timestep(lvl),2) == 1) return

  call amrex_average_down_nodes(mfrho(lvl),mfrho(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfu_vel(lvl),mfu_vel(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfv_vel(lvl),mfv_vel(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mftemp(lvl),mftemp(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfalpha(lvl),mfalpha(lvl-1),1,1,(/2,2,2/))

  call amrex_average_down_nodes(mfchi(lvl),mfchi(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfzeta_x(lvl),mfzeta_x(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfzeta_y(lvl),mfzeta_y(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfpi_xx(lvl),mfpi_xx(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfpi_yy(lvl),mfpi_yy(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mfpi_xy(lvl),mfpi_xy(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mflambda_x(lvl),mflambda_x(lvl-1),1,1,(/2,2,2/))
  call amrex_average_down_nodes(mflambda_y(lvl),mflambda_y(lvl-1),1,1,(/2,2,2/))

!
! 3-Dimensional junk
!
  if (dimensions == 3) then
    call amrex_average_down_nodes(mfpi_zz(lvl),mfpi_zz(lvl-1),1,1,(/2,2,2/))
    call amrex_average_down_nodes(mfpi_xz(lvl),mfpi_xz(lvl-1),1,1,(/2,2,2/))
    call amrex_average_down_nodes(mfpi_yz(lvl),mfpi_yz(lvl-1),1,1,(/2,2,2/))
    call amrex_average_down_nodes(mflambda_z(lvl),mflambda_z(lvl-1),1,1,(/2,2,2/))
    call amrex_average_down_nodes(mfzeta_z(lvl),mfzeta_z(lvl-1),1,1,(/2,2,2/))

    call amrex_average_down_nodes(mfw_vel(lvl),mfw_vel(lvl-1),1,1,(/2,2,2/))

  end if
!
!
!
  call amrex_average_down_nodes(mffi(lvl),mffi(lvl-1),1,dir,(/2,2,2/))
  call amrex_average_down_nodes(mfgi(lvl),mfgi(lvl-1),1,dir,(/2,2,2/))
  call amrex_average_down_nodes(mffout(lvl),mffout(lvl-1),1,dir,(/2,2,2/))
  call amrex_average_down_nodes(mfgout(lvl),mfgout(lvl-1),1,dir,(/2,2,2/))
end subroutine


!if (mod(timestep(lvl),2) == 0) then
!
!  write(*,*) 'averaging down',self
!  call amrex_average_down_nodes_w_ghosts(mfrho(lvl),mfrho(lvl-1),1,1,(/2,2,2/),nghosts)
!!  write(*,*) 'after rho',self
!  call amrex_average_down_nodes_w_ghosts(mfu_vel(lvl),mfu_vel(lvl-1),1,1,(/2,2,2/),nghosts)
!!  write(*,*) 'after u',self
!  call amrex_average_down_nodes_w_ghosts(mfv_vel(lvl),mfv_vel(lvl-1),1,1,(/2,2,2/),nghosts)
!!  write(*,*) 'after v',self
!  call amrex_average_down_nodes_w_ghosts(mftemp(lvl),mftemp(lvl-1),1,1,(/2,2,2/),nghosts)
!  !call amrex_average_down_nodal(mfalpha(lvl),mfalpha(lvl-1),1,1,amrex_ref_ratio(lvl))
!!  write(*,*) 'fiddlesticks',self
!!
!! LMs
!!
!  call amrex_average_down_nodes_w_ghosts(mfchi(lvl),mfchi(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mfzeta_x(lvl),mfzeta_x(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mfzeta_y(lvl),mfzeta_y(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mfpi_xx(lvl),mfpi_xx(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mfpi_yy(lvl),mfpi_yy(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mfpi_xy(lvl),mfpi_xy(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mflambda_x(lvl),mflambda_x(lvl-1),1,1,(/2,2,2/),nghosts)
!  call amrex_average_down_nodes_w_ghosts(mflambda_y(lvl),mflambda_y(lvl-1),1,1,(/2,2,2/),nghosts)
!!  write(*,*) 'LMs averaged'
!!
!! 3-Dimensional junk
!!
!  if (dimensions == 3) then
!    call amrex_average_down_nodes_w_ghosts(mfpi_zz(lvl),mfpi_zz(lvl-1),1,1,(/2,2,2/),nghosts)
!    call amrex_average_down_nodes_w_ghosts(mfpi_xz(lvl),mfpi_xz(lvl-1),1,1,(/2,2,2/),nghosts)
!    call amrex_average_down_nodes_w_ghosts(mfpi_yz(lvl),mfpi_yz(lvl-1),1,1,(/2,2,2/),nghosts)
!    call amrex_average_down_nodes_w_ghosts(mflambda_z(lvl),mflambda_z(lvl-1),1,1,(/2,2,2/),nghosts)
!    call amrex_average_down_nodes_w_ghosts(mfzeta_z(lvl),mfzeta_z(lvl-1),1,1,(/2,2,2/),nghosts)
!
!    call amrex_average_down_nodes_w_ghosts(mfw_vel(lvl),mfw_vel(lvl-1),1,1,(/2,2,2/),nghosts)
!
!  end if
!
!
!  write(*,*) 'averaging down',self
!  call amrex_average_down_nodal(mfrho(lvl),mfrho(lvl-1),(/2,2,2/),nghosts)
!!  write(*,*) 'after rho',self
!  call amrex_average_down_nodal(mfu_vel(lvl),mfu_vel(lvl-1),(/2,2,2/),nghosts)
!!  write(*,*) 'after u',self
!  call amrex_average_down_nodal(mfv_vel(lvl),mfv_vel(lvl-1),(/2,2,2/),nghosts)
!!  write(*,*) 'after v',self
!  call amrex_average_down_nodal(mftemp(lvl),mftemp(lvl-1),(/2,2,2/),nghosts)
!  !call amrex_average_down_nodal(mfalpha(lvl),mfalpha(lvl-1),1,1,amrex_ref_ratio(lvl))
!!  write(*,*) 'fiddlesticks',self
!!
!! LMs
!!
!  call amrex_average_down_nodal(mfchi(lvl),mfchi(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mfzeta_x(lvl),mfzeta_x(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mfzeta_y(lvl),mfzeta_y(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mfpi_xx(lvl),mfpi_xx(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mfpi_yy(lvl),mfpi_yy(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mfpi_xy(lvl),mfpi_xy(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mflambda_x(lvl),mflambda_x(lvl-1),(/2,2,2/),nghosts)
!  call amrex_average_down_nodal(mflambda_y(lvl),mflambda_y(lvl-1),(/2,2,2/),nghosts)
!!  write(*,*) 'LMs averaged'
!!
!! 3-Dimensional junk
!!
!  if (dimensions == 3) then
!    call amrex_average_down_nodal(mfpi_zz(lvl),mfpi_zz(lvl-1),(/2,2,2/),nghosts)
!    call amrex_average_down_nodal(mfpi_xz(lvl),mfpi_xz(lvl-1),(/2,2,2/),nghosts)
!    call amrex_average_down_nodal(mfpi_yz(lvl),mfpi_yz(lvl-1),(/2,2,2/),nghosts)
!    call amrex_average_down_nodal(mflambda_z(lvl),mflambda_z(lvl-1),(/2,2,2/),nghosts)
!    call amrex_average_down_nodal(mfzeta_z(lvl),mfzeta_z(lvl-1),(/2,2,2/),nghosts)
!
!    call amrex_average_down_nodal(mfw_vel(lvl),mfw_vel(lvl-1),(/2,2,2/),nghosts)
!
!  end if

