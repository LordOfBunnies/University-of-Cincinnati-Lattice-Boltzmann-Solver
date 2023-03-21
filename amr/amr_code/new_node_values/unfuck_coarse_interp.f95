subroutine unfuck_coarse_interp(lvl)
!
! Fuck you Amrex
!
!
!
!
!
!
use amrex_base_module
use mpi
use amrex_amr_module
use precise
USE grid_data
use quick_calcs
use amr_info_holder
use amr_processes
USE linkwise
use constants

implicit none
integer :: a,b,c,i,lvl

type(amrex_mfiter) :: mfi
type(amrex_box) :: box


call amrex_mfiter_build(mfi,mffout(lvl),tiling=.false.)
do while(mfi%next())

  fi => mffi(lvl)%dataptr(mfi)
  fout => mffout(lvl)%dataptr(mfi)
  gi => mfgi(lvl)%dataptr(mfi)
  gout => mfgout(lvl)%dataptr(mfi)

  box = mfi%validbox()

  DO c = box%lo(3)-nghosts_mid,box%hi(3)+nghosts_mid
    DO b = box%lo(2)-nghosts_mid,box%hi(2)+nghosts_mid
      DO a = box%lo(1)-nghosts_mid,box%hi(1)+nghosts_mid
        if (a == 374 .and. b == 256 .and. c == 256) then
          write(*,*) 'merde merde merde',fi(a,b,c,1:dir)
          write(*,*) 'skiite',fout(a,b,c,1:dir)
          write(*,*) 'pendejo',gi(a,b,c,1:dir)
          write(*,*) 'chinga-te',gout(a,b,c,1:dir)
        end if

        do i = 1,dir

        if (fi(a,b,c,i) < 0.0D0) then
          fi(a,b,c,i) = abs(fi(a,b,c,i))
          if (fi(a,b,c,i) < 1.0D-6) then
            fi(a,b,c,i) = 10.0D0*fi(a,b,c,i)
          end if

        end if

        if (fout(a,b,c,i) < 0.0D0) then
         write(*,*) 'bingo fuel',a,b,c,i
          fout(a,b,c,i) = abs(fout(a,b,c,i))
          if (fout(a,b,c,i) < 1.0D-6) then
            fout(a,b,c,i) = 10.0D0*fout(a,b,c,i)
          end if
        end if

        if (gi(a,b,c,i) < 0.0D0) then
          gi(a,b,c,i) = abs(gi(a,b,c,i))
          if (gi(a,b,c,i) < 1.0D-6) then
            gi(a,b,c,i) = 10.0D0*gi(a,b,c,i)
          end if
        end if

        if (gout(a,b,c,i) < 0.0D0) then
          gout(a,b,c,i) = abs(gout(a,b,c,i))
          if (gout(a,b,c,i) < 1.0D-6) then
            gout(a,b,c,i) = 10.0D0*gout(a,b,c,i)
          end if
        end if

        end do
      END DO
    END DO
  END DO

end do

call amrex_mfiter_destroy(mfi)

end subroutine



