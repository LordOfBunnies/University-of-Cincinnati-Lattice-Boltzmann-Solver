subroutine startup_transient_damping(loc_fout,loc_fi,loc_fi_opp,loc_fieq_opp,fout_adj,&
  loc_gout,loc_gi,loc_gi_opp,loc_gieq_opp,&
  gout_adj,loc_state,loc_state_opp,lvl)
!
! This helps damp out some of the major transients which occur during startup
!
! Equilibrium values
!
! loc_*out is the out population going into the wall
! loc_*i is fi/gi the population going into the wall
! loc_*i_opp is the missing population in the opposite direction of the wall
! loc_*ieq is the equilibrium population opposite the wall
! loc_*ieq_opp is the equilibrium population going into the wall
! loc_state is into the wall
! loc_state_opp is away from the wall
!
! *_adj terms are adjacent values to ensure that the modifications to *i are not
!   lost during streaming
!
! Called by: interpolative_bounceback_bc
! Calls:
!
!
use precise
use amr_info_holder, only: startup_iterations,elbm_epsilon,init_epsilon


implicit none
real(kind=dp) :: loc_fout,loc_fi,loc_fi_opp,loc_fieq_opp,&
  loc_gout,loc_gi,loc_gi_opp,loc_gieq_opp
real(kind=dp) :: fout_adj,gout_adj
real(kind=dp) :: adjustment,damping_ratio
integer,intent(in) :: loc_state,loc_state_opp,lvl

damping_ratio = elbm_epsilon(lvl)!1.0D0/real(startup_iterations,kind=dp)
!
! Away from the wall
!
!if (loc_fout > loc_fi_opp) then
! damping value positive, for directions like 34 where the bounced back value
!   is small and needs to be increased slowly
  loc_fout = loc_fi_opp + damping_ratio*(loc_fout-loc_fi_opp)
!else
! damping value negative, for directions like 37 where the bounced back value
!   is large, but need to be decreased slowly
!  loc_fout = loc_fi_opp + damping_ratio*(loc_fi_opp-loc_fout)
!end if
!
! Into the wall
!
!if (loc_fieq_opp > loc_fi) then
! For directions into the wall, this will tend to pull them towards the
!   equilibrium value, whether it is larger or smaller than the relevant
!   fi
  loc_fi = loc_fi + damping_ratio*(loc_fieq_opp-loc_fi)
!else
!  loc_fi = loc_fi + damping_ratio*(loc_fi-loc_fieq_opp)
!end if
!
! Away from the wall
!
!if (loc_gout > loc_gi_opp) then

!  loc_gout = loc_gi_opp + damping_ratio*(loc_gi_opp-loc_gout)
!else
  loc_gout = loc_gi_opp + damping_ratio*(loc_gout-loc_gi_opp)
!end if
!
! Into the wall
!
!if (loc_gieq_opp > loc_gi) then
  loc_gi = loc_gi + damping_ratio*(loc_gieq_opp-loc_gi)
!else
!  loc_gi = loc_gi + damping_ratio*(loc_gi-loc_gieq_opp)
!end if

if (loc_state_opp >=0) then
  fout_adj = loc_fi
  gout_adj = loc_gi
end if

end subroutine
