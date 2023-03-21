subroutine ggm_output_shortcut(g,fine)
!
! To get the output of the GGM routines so you can see what's going on with the calculation.
!
!
!
! Called by: lb_output
! Calls: populate_q,calculate_grad_mag,calculate_gradgrad_mag
! External calls: amrex_multifab_build,amrex_mfiter_build
!
use precise
use constants
use ggm_stuff
use amrex_amr_module
use amrex_base_module
use amr_info_holder
use amr_processes
use mpi

implicit none
integer :: lvl,g,max_lvl,ier
INTEGER :: taglow(4),taghigh(4),fine

logical :: for_output
!
!CHARACTER(kind = c_char),intent(in),value :: settag,cleartag
!CHARACTER(kind = c_char),contiguous,pointer :: tagarr(:,:,:,:)

!TYPE(amrex_tagboxarray) :: tag
TYPE(amrex_box) :: larks
TYPE(amrex_mfiter) :: mfi

!do g = 1,num_ggm_vars

if (fine > ggm_max_lvl(g)) then
  max_lvl = ggm_max_lvl(g)
else
  max_lvl = fine
end if
!
!
!
!write(*,*) 'cuddly kitty!',g,lvl,max_lvl,ggm_min_lvl(g)
!do lvl = ggm_min_lvl(g),max_lvl
do lvl = 0,max_lvl
!
!
!
!  call amrex_multifab_build(mfgm(lvl),mffi(lvl)%ba,mffi(lvl)%dm,num_ggm_vars,&
!    nghosts,node_based_3d)
!  call amrex_multifab_build(mfggm(lvl),mffi(lvl)%ba,mffi(lvl)%dm,num_ggm_vars,&
!    0,node_based_3d)
!
!
!
  CALL amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
!
  do while (mfi%next())
    larks = mfi%validbox()
!    write(*,*) 'fuzzy bunnies!',self,lvl
    if(allocated(holder)) deallocate(holder)
    allocate(holder(larks%lo(1)-nghosts_mid:larks%hi(1)+nghosts_mid,&
      larks%lo(2)-nghosts_mid:larks%hi(2)+nghosts_mid,&
      larks%lo(3)-nghosts_mid:larks%hi(3)+nghosts_mid,g:g))

    call populate_q(lvl,larks%lo,larks%hi,g,mfi)
!    write(*,*) 'fluffy bears!',self,lvl
    !call mpi_barrier(commune,ier)

    call calculate_grad_mag(lvl,larks%lo,larks%hi,g,mfi)
!    write(*,*) 'soft deer!',self,lvl
    !call mpi_barrier(commune,ier)

    call calculate_ggm_for_output(lvl,larks%lo,larks%hi,g,mfi)
!    write(*,*) 'hi, my name is', self
    deallocate(holder)

  end do


!  write(*,*) 'ingle?',self
  call amrex_mfiter_destroy(mfi)

end do
call mpi_barrier(commune,ier)

end subroutine

!write(*,*) 'bonk?',self
!call ggm_data_write(g,fine)

!do lvl = ggm_min_lvl(g),max_lvl
!  call amrex_multifab_destroy(mfgm(lvl))
!  call amrex_multifab_destroy(mfggm(lvl))
!end do
