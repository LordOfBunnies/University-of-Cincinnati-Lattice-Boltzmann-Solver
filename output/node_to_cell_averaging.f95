subroutine node_to_cell_averaging(mfcell_rho,max_lvl)
!
!
!
!
!
!
!
use precise
use amr_info_holder
use amr_processes
use amrex_amr_module
use amrex_base_module


implicit none
integer :: lvl,a,b,c,max_lvl
integer :: valid_counter
real(kind=dp) :: cell_value
real(kind=dp),contiguous,pointer :: rho_cell(:,:,:,:)

type(amrex_multifab) :: mfcell_rho(0:max_lvl)
type(amrex_mfiter) :: mfi
type(amrex_box) :: dorx


do lvl = 0,max_lvl

  call amrex_mfiter_build(mfi,mfcell_rho(lvl),tiling=.false.)
  write(*,*) 'Thors day celebration!',self
  do while(mfi%next())
    dorx = mfi%validbox()

    state => mfstate(lvl)%dataptr(mfi)

    rho => mfrho(lvl)%dataptr(mfi)
    rho_cell => mfcell_rho(lvl)%dataptr(mfi)

    do c = dorx%lo(3),dorx%hi(3)
      do b = dorx%lo(2),dorx%hi(2)
        do a = dorx%lo(1),dorx%hi(1)
          valid_counter = 0
          cell_value = 0.0D0

          if (state(a,b,c,1) < -1000) then
            rho_cell(a,b,c,1) = 1.0D0
            cycle
          end if

          if (state(a,b,c,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a,b,c,1)
          end if
          if (state(a+1,b,c,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a+1,b,c,1)
          end if
          if (state(a,b+1,c,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a,b+1,c,1)
          end if
          if (state(a,b,c+1,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a,b,c+1,1)
          end if
          if (state(a+1,b+1,c,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a+1,b+1,c,1)
          end if
          if (state(a+1,b,c+1,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a+1,b,c+1,1)
          end if
          if (state(a,b+1,c+1,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a,b+1,c+1,1)
          end if
          if (state(a+1,b+1,c+1,1) >= 0) then
            valid_counter = valid_counter + 1
            cell_value = cell_value + rho(a+1,b+1,c+1,1)
          end if

          if (valid_counter == 0) then
            rho_cell(a,b,c,1) = 1.0D0
          else
            rho_cell(a,b,c,1) = cell_value/real(valid_counter,kind=dp)
          end if

!          write(*,*) 'nd/cc averaging',rho_cell(a,b,c,1),self,lvl

        end do
      end do
    end do



  end do

  call amrex_mfiter_destroy(mfi)

end do

end subroutine
