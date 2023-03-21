subroutine state_copy(lvl,imfstate_temp,tofrom)
!
! This copies old state information into the temporary multifab
!
!
!
! Called by: remake_level (amr_processes)
! Calls:
! External call: amrex_mfiter_build
!
use constants
use grid_data, only: dir
use amr_info_holder, only: mfstate,state,nghosts_state,nghosts
use amrex_amr_module
use amrex_base_module
implicit none
integer :: a,b,c,i
integer,intent(in) :: lvl
integer,contiguous,pointer :: state_temp(:,:,:,:)
logical :: tofrom

type(amrex_imultifab) :: imfstate_temp
type(amrex_mfiter) :: mfi
type(amrex_box) :: tax

!
! If true, copies information to the temporary multifab
!


if (tofrom) then
  CALL amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
  do while (mfi%next())
    tax = mfi%validbox()

    state => mfstate(lvl)%dataptr(mfi)
    state_temp => imfstate_temp%dataptr(mfi)
    do i = 1,dir
      do c = tax%lo(3)-nghosts_state,tax%hi(3)+nghosts_state
        do b = tax%lo(2)-nghosts_state,tax%hi(2)+nghosts_state
          do a = tax%lo(1)-nghosts_state,tax%hi(1)+nghosts_state

            state_temp(a,b,c,i) = state(a,b,c,i)
          end do
        end do
      end do
    end do
  end do
else
  CALL amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
  do while (mfi%next())
    tax = mfi%validbox()

    state => mfstate(lvl)%dataptr(mfi)
    state_temp => imfstate_temp%dataptr(mfi)
    do i = 1,dir
      do c = tax%lo(3)-nghosts_state,tax%hi(3)+nghosts_state
        do b = tax%lo(2)-nghosts_state,tax%hi(2)+nghosts_state
          do a = tax%lo(1)-nghosts_state,tax%hi(1)+nghosts_state

            state(a,b,c,i) = state_temp(a,b,c,i)
          end do
        end do
      end do
    end do
  end do
end if


end subroutine
