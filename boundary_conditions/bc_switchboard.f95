subroutine bc_switchboard(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc,&
  loc_dir,bc_num)
!
! Puts the right call in the right place
!
!
! Called by: Streaming
! Calls: rho_v_bc,press_temp_bc,freestream_bc,v_temp_bc,null_bc,freestream_edge_bc,
!   bounceback,inviscid_boundary
! External calls:
!
use precise
use constants

implicit none
integer :: bc_num,loc_dir
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
    call inviscid_edge(loc_dir,2,loc_dir)
  case(-12)
    call inviscid_wall(loc_dir,2,loc_dir)
end select

end subroutine
