subroutine gradient_amr(lvl,settag,cleartag,tagarr,tux_lo,tux_hi,&
  tag_lo,tag_hi,mfi)
!
! Second order gradient derivative to determime if the mesh is well
! resolved at a given node.
!
! Calls: populate_q,calculate_grad_mag,calculate_gradgrad_mag
! Called by: tag_nodes in amr_processes
!
USE precise
use output_data, only: ggm_steps
USE ggm_stuff
use constants
use amrex_amr_module
use amrex_base_module
use amr_info_holder
use amr_processes
use mpi
IMPLICIT NONE
INTEGER :: lvl,g,a,b,c,i,fine,ier
integer,intent(in) :: tux_lo(3),tux_hi(3),tag_lo(4),tag_hi(4)
character(kind=c_char),intent(in) :: settag,cleartag
REAL(KIND = dp) :: elapsed_ggm_time
logical :: for_output

character(kind=c_char),intent(inout) :: tagarr(tag_lo(1):tag_hi(1),&
  tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))

type(amrex_box) :: tux
type(amrex_mfiter) :: mfi
!
! Analyze what needs to be done
!
!write(*,*) 'steve?',lvl
fine = amrex_get_finest_level()
!for_output = .false.

if(allocated(holder)) deallocate(holder)
allocate(holder(tux_lo(1)-nghosts_mid:tux_hi(1)+nghosts_mid,&
  tux_lo(2)-nghosts_mid:tux_hi(2)+nghosts_mid,&
  tux_lo(3)-nghosts_mid:tux_hi(3)+nghosts_mid,num_ggm_vars))
!
  DO g = 1,num_ggm_vars

      call populate_q(lvl,tux_lo,tux_hi,g,mfi)

      call calculate_grad_tag(lvl,settag,cleartag,tux_lo,tux_hi,tag_lo,tag_hi,&
        tagarr,g,mfi)

    end do

deallocate(holder)

END SUBROUTINE

