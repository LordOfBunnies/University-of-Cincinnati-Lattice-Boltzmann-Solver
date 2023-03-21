module bc_info
!
! This has been created to contain some of the boundary condition information to
!   put things in an external interface to allow cool Fortran tricks
!
!
!
!
!
use precise
use constants
use grid_data
implicit none
save

private
public :: bc_switchboard



contains
!
!
!
subroutine bc_switchboard(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc,&
  loc_dir,bc_num,out_dir,loc_state,loc_cx,loc_cy,loc_cz)
!
! Puts the right call in the right place
!
!
! Called by: Streaming
! Calls: rho_v_bc,press_temp_bc,freestream_bc,v_temp_bc,null_bc,freestream_edge_bc,
!   bounceback,inviscid_boundary
! External calls:
!
implicit none
integer :: bc_num,loc_dir
integer,optional :: loc_state(:),out_dir,loc_cx,loc_cy,loc_cz
real(kind=dp),intent(out) :: fout_bc,fout_opp_bc,gout_bc,gout_opp_bc
logical :: inlet

select case(bc_num)
  case(-1)
    call rho_v_bc()
  case(-2)
    call press_temp_bc()
  case(-3)
    call freestream_bc(fout_bc,gout_bc,loc_dir)
  case(-4)
    call v_temp_bc()
  case(-5)
    call null_bc
  case(-8)
    call freestream_edge_bc(fout_bc,gout_bc,loc_dir)
  case(-10)
    call bounceback(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc)
  case(-11)
    call inviscid_edge(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc,&
      loc_dir,bc_num,loc_state,out_dir,loc_cx,loc_cy,loc_cz)
  case(-12)
    call inviscid_wall(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc,&
      loc_dir,bc_num,loc_state,out_dir,loc_cx,loc_cy,loc_cz)
end select

end subroutine


end module
