subroutine sanity_checks
!
!
!
!
!
!
!
use amr_info_holder, only: ggm_save_info,grid_info
use amr_processes, only: amr_max_lvl
use ggm_stuff
use precise
use nml_output
implicit none
integer :: max_lvl,i

do i = 1,num_saves

  if (save_shapes(2,i) > amr_max_lvl) then
    save_shapes(2,i) = amr_max_lvl
  end if

end do

if (use_ggm) then
  if (ggm_output) then
    do i = 1,num_ggm_vars

      if (ggm_max_lvl(i) > amr_max_lvl) then
        ggm_max_lvl(i) = amr_max_lvl
      end if

    end do
  end if
end if


end subroutine
