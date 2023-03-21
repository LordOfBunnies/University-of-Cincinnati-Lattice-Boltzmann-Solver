subroutine inviscid_wall(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc,&
      i,bc_num,loc_state,inviscid_dir,loc_cx,loc_cy,loc_cz)
!
! This is mostly a placeholder for later work as non-boundary inviscid walls are
!   more complicated (I think)
!
!
!
!
!
use precise
use constants
use grid_data
implicit none
INTEGER,INTENT(IN) :: i,loc_state(dir),bc_num
integer :: loc_cx,loc_cy,loc_cz,inviscid_dir
real(kind=dp) :: fout_bc,fout_opp_bc,gout_bc,gout_opp_bc






end subroutine
