subroutine damp_it(level_step,fine,lvl)
!
! Frankly my dear, I don't give a damp
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

INTEGER :: a,b,c,i,reach,fine,lvl
LOGICAL :: level_step(0:amr_max_lvl)

type(amrex_mfiter) :: mfi
type(amrex_box) :: fux


if (.not. shifted) then
if (lvl == 0) then
  reach = nghosts/2
else
  if (mod(timestep(lvl),2) == 0) then
    reach = nghosts
  else
    reach = nghosts
  end if
end if
else
if (lvl == 0) then
  reach = nghosts/2
else
  if (mod(timestep(lvl),2) == 0) then
    reach = nghosts
  else
    reach = nghosts
  end if
end if
end if


CALL amrex_mfiter_build(mfi,mffout(lvl),tiling=.false.)
do while(mfi%next())
  fi => mffi(lvl)%dataptr(mfi)
  gi => mfgi(lvl)%dataptr(mfi)
  fout => mffout(lvl)%dataptr(mfi)
  gout => mfgout(lvl)%dataptr(mfi)
  state => mfstate(lvl)%dataptr(mfi)

  fux = mfi%validbox()
!
!
!
  do i = 1,dir
    DO c = fux%lo(3)-reach,fux%hi(3)+reach
      DO b = fux%lo(2)-reach,fux%hi(2)+reach
        DO a = fux%lo(1)-reach,fux%hi(1)+reach
!
!
!
          if (fout(a,b,c,i) < 0.0D0) then

            if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) >= 0 .and. &
                .not. shifted) then

              if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + fout(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
              else if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  fout(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
              else if (fout(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                  fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                fout(a,b,c,i) = (fout(a-cx(i),b-cy(i),c-cz(i),i) + abs(fout(a,b,c,i)))/2
              else
                fout(a,b,c,i) = abs(fout(a,b,c,i))
              end if
!
!
!
!            else if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500 &
!               .and. shifted) then
!
!              if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + fout(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
!              else if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  fout(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
!                fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
!              else if (fout(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
!                  fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                fout(a,b,c,i) = (fout(a-cx(i),b-cy(i),c-cz(i),i) + abs(fout(a,b,c,i)))/2
!              else
!                fout(a,b,c,i) = abs(fout(a,b,c,i))
!              end if
!
!
!
            else if (.not. shifted) then
              if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) < 0) then
                if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                  fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
                else
                  fout(a,b,c,i) = abs(fout(a,b,c,i))
                end if
              else if (state(a,b,c,i) < 0 .and. state(a,b,c,opp(i)) >= 0) then
                if (fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                  fout(a,b,c,i) = (fout(a-cx(i),b-cy(i),c-cz(i),i) + abs(fout(a,b,c,i)))/2
                else
                  fout(a,b,c,i) = abs(fout(a,b,c,i))
                end if
              else
                fout(a,b,c,i) = abs(fout(a,b,c,i))
              end if
!
!
!
            else if (shifted) then
              if (a-cx(i) < fux%lo(1)-nghosts_mid .or. a-cx(i) < fux%lo(2)-nghosts_mid .or. &
                  c-cz(i) < fux%lo(3)-nghosts_mid .or. a+cx(i) > fux%hi(1)-nghosts_mid .or. &
                  b+cy(i) > fux%hi(2)-nghosts_mid .or. c+cz(i) > fux%hi(3)-nghosts_mid) then

                if (state(a,b,c,i) >= -500) then
                  fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
                else
                  fout(a,b,c,i) = abs(fout(a,b,c,i))
                end if
!
              else
                if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then

                  if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + fout(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
                  else if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      fout(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                    fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
                  else if (fout(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                      fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    fout(a,b,c,i) = (fout(a-cx(i),b-cy(i),c-cz(i),i) + abs(fout(a,b,c,i)))/2
                  else
                    fout(a,b,c,i) = abs(fout(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) < -500) then
                  if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                    fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
                  else
                    fout(a,b,c,i) = abs(fout(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) < -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then
                  if (fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    fout(a,b,c,i) = (fout(a-cx(i),b-cy(i),c-cz(i),i) + abs(fout(a,b,c,i)))/2
                  else
                    fout(a,b,c,i) = abs(fout(a,b,c,i))
                  end if
                end if
!                else
!                  fout(a,b,c,i) = abs(fout(a,b,c,i))
!                end if

              end if

            else
              fout(a,b,c,i) = abs(fout(a,b,c,i))
            end if

          end if
!
!
!
!
!
          if (fi(a,b,c,i) < 0.0D0) then

            if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) >= 0 .and. &
                .not. shifted) then

              if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + fi(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
              else if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  fi(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + abs(fi(a,b,c,i)))/2
              else if (fi(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                  fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                fi(a,b,c,i) = (fi(a-cx(i),b-cy(i),c-cz(i),i) + abs(fi(a,b,c,i)))/2
              else
                fi(a,b,c,i) = abs(fi(a,b,c,i))
              end if

!
!
            else if (.not. shifted) then
              if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) < 0) then
                if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                  fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + abs(fi(a,b,c,i)))/2
                else
                  fi(a,b,c,i) = abs(fi(a,b,c,i))
                end if
              else if (state(a,b,c,i) < 0 .and. state(a,b,c,opp(i)) >= 0) then
                if (fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                  fi(a,b,c,i) = (fi(a-cx(i),b-cy(i),c-cz(i),i) + abs(fi(a,b,c,i)))/2
                else
                  fi(a,b,c,i) = abs(fi(a,b,c,i))
                end if
              else
                fi(a,b,c,i) = abs(fi(a,b,c,i))
              end if
!
!
!
            else if (shifted) then
              if (a-cx(i) < fux%lo(1)-nghosts_mid .or. a-cx(i) < fux%lo(2)-nghosts_mid .or. &
                  c-cz(i) < fux%lo(3)-nghosts_mid .or. a+cx(i) > fux%hi(1)-nghosts_mid .or. &
                  b+cy(i) > fux%hi(2)-nghosts_mid .or. c+cz(i) > fux%hi(3)-nghosts_mid) then

                if (state(a,b,c,i) >= -500) then
                  fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + abs(fi(a,b,c,i)))/2
                else
                  fi(a,b,c,i) = abs(fi(a,b,c,i))
                end if
!
              else
                if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then

                  if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + fi(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
                  else if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      fi(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                    fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + abs(fi(a,b,c,i)))/2
                  else if (fi(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                      fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    fi(a,b,c,i) = (fi(a-cx(i),b-cy(i),c-cz(i),i) + abs(fi(a,b,c,i)))/2
                  else
                    fi(a,b,c,i) = abs(fi(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) < -500) then
                  if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                    fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + abs(fi(a,b,c,i)))/2
                  else
                    fi(a,b,c,i) = abs(fi(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) < -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then
                  if (fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    fi(a,b,c,i) = (fi(a-cx(i),b-cy(i),c-cz(i),i) + abs(fi(a,b,c,i)))/2
                  else
                    fi(a,b,c,i) = abs(fi(a,b,c,i))
                  end if
                end if
!              else
!                fi(a,b,c,i) = abs(fi(a,b,c,i))
!              end if

              end if
            else
              fi(a,b,c,i) = abs(fi(a,b,c,i))
            end if

          end if
!
!
!
!
!
          if (gout(a,b,c,i) < 0.0D0) then

            if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) >= 0 .and. &
                .not. shifted) then

              if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + gout(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
              else if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  gout(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + abs(gout(a,b,c,i)))/2
              else if (gout(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                  gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                gout(a,b,c,i) = (gout(a-cx(i),b-cy(i),c-cz(i),i) + abs(gout(a,b,c,i)))/2
              else
                gout(a,b,c,i) = abs(gout(a,b,c,i))
              end if
!
!
!
            else if (.not. shifted) then
              if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) < 0) then
                if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                  gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + abs(gout(a,b,c,i)))/2
                else
                  gout(a,b,c,i) = abs(gout(a,b,c,i))
                end if
              else if (state(a,b,c,i) < 0 .and. state(a,b,c,opp(i)) >= 0) then
                if (gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                  gout(a,b,c,i) = (gout(a-cx(i),b-cy(i),c-cz(i),i) + abs(gout(a,b,c,i)))/2
                else
                  gout(a,b,c,i) = abs(gout(a,b,c,i))
                end if
              else
                gout(a,b,c,i) = abs(gout(a,b,c,i))
              end if
!
!
!
            else if (shifted) then
              if (a-cx(i) < fux%lo(1)-nghosts_mid .or. a-cx(i) < fux%lo(2)-nghosts_mid .or. &
                  c-cz(i) < fux%lo(3)-nghosts_mid .or. a+cx(i) > fux%hi(1)-nghosts_mid .or. &
                  b+cy(i) > fux%hi(2)-nghosts_mid .or. c+cz(i) > fux%hi(3)-nghosts_mid) then

                if (state(a,b,c,i) >= -500) then
                  gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + abs(gout(a,b,c,i)))/2
                else
                  gout(a,b,c,i) = abs(gout(a,b,c,i))
                end if
!
              else
                if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then

                  if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + gout(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
                  else if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      gout(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                    gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + abs(gout(a,b,c,i)))/2
                  else if (gout(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                      gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    gout(a,b,c,i) = (gout(a-cx(i),b-cy(i),c-cz(i),i) + abs(gout(a,b,c,i)))/2
                  else
                    gout(a,b,c,i) = abs(gout(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) < -500) then
                  if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                    gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + abs(gout(a,b,c,i)))/2
                  else
                    gout(a,b,c,i) = abs(gout(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) < -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then
                  if (gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    gout(a,b,c,i) = (gout(a-cx(i),b-cy(i),c-cz(i),i) + abs(gout(a,b,c,i)))/2
                  else
                    gout(a,b,c,i) = abs(gout(a,b,c,i))
                  end if
                end if
!                else
!                  gout(a,b,c,i) = abs(gout(a,b,c,i))
!                end if

              end if
            else
              gout(a,b,c,i) = abs(gout(a,b,c,i))
            end if

          end if
!
!
!
!
!
!
!
!
          if (gi(a,b,c,i) < 0.0D0) then

            if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) >= 0 .and. &
                .not. shifted) then

              if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + gi(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
              else if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                  gi(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + abs(gi(a,b,c,i)))/2
              else if (gi(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                  gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                gi(a,b,c,i) = (gi(a-cx(i),b-cy(i),c-cz(i),i) + abs(gi(a,b,c,i)))/2
              else
                gi(a,b,c,i) = abs(gi(a,b,c,i))
              end if
!
!
!
            else if (.not. shifted) then
              if (state(a,b,c,i) >= 0 .and. state(a,b,c,opp(i)) < 0) then
                if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                  gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + abs(gi(a,b,c,i)))/2
                else
                  gi(a,b,c,i) = abs(gi(a,b,c,i))
                end if
              else if (state(a,b,c,i) < 0 .and. state(a,b,c,opp(i)) >= 0) then
                if (gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                  gi(a,b,c,i) = (gi(a-cx(i),b-cy(i),c-cz(i),i) + abs(gi(a,b,c,i)))/2
                else
                  gi(a,b,c,i) = abs(gi(a,b,c,i))
                end if
              else
                gi(a,b,c,i) = abs(gi(a,b,c,i))
              end if
!
!
!
            else if (shifted) then
              if (a-cx(i) < fux%lo(1)-nghosts_mid .or. a-cx(i) < fux%lo(2)-nghosts_mid .or. &
                  c-cz(i) < fux%lo(3)-nghosts_mid .or. a+cx(i) > fux%hi(1)-nghosts_mid .or. &
                  b+cy(i) > fux%hi(2)-nghosts_mid .or. c+cz(i) > fux%hi(3)-nghosts_mid) then

                if (state(a,b,c,i) >= -500) then
                  gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + abs(gi(a,b,c,i)))/2
                else
                  gi(a,b,c,i) = abs(gi(a,b,c,i))
                end if
!
              else
                if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then

                  if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + gi(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
                  else if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
                      gi(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
                    gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + abs(gi(a,b,c,i)))/2
                  else if (gi(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
                      gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    gi(a,b,c,i) = (gi(a-cx(i),b-cy(i),c-cz(i),i) + abs(gi(a,b,c,i)))/2
                  else
                    gi(a,b,c,i) = abs(gi(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) < -500) then
                  if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
                    gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + abs(gi(a,b,c,i)))/2
                  else
                    gi(a,b,c,i) = abs(gi(a,b,c,i))
                  end if
!
                else if (state(a,b,c,i) < -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= -500) then
                  if (gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
                    gi(a,b,c,i) = (gi(a-cx(i),b-cy(i),c-cz(i),i) + abs(gi(a,b,c,i)))/2
                  else
                    gi(a,b,c,i) = abs(gi(a,b,c,i))
                  end if
                end if
!                else
!                  gi(a,b,c,i) = abs(gi(a,b,c,i))
!                end if

              end if
            else
              gi(a,b,c,i) = abs(gi(a,b,c,i))
            end if

          end if

        end do
      end do
    end do

  end do

end do
call amrex_mfiter_destroy(mfi)

end subroutine




!
!              if (state(a,b,c,i) >= -500 .and. state(a-cx(i),b-cy(i),c-cz(i),1) < 0) then
!                if (fout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0) then
!                  fout(a,b,c,i) = (fout(a+cx(i),b+cy(i),c+cz(i),i) + abs(fout(a,b,c,i)))/2
!                else
!                  fout(a,b,c,i) = abs(fout(a,b,c,i))
!                end if
!              else if (state(a,b,c,i) < 0 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= 0) then
!                if (fout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                  fout(a,b,c,i) = (fout(a-cx(i),b-cy(i),c-cz(i),i) + abs(fout(a,b,c,i)))/2
!                else
!                  fout(a,b,c,i) = abs(fout(a,b,c,i))
!                end if
!              else
!
!              end if
!
!
!
!            else if (state(a,b,c,i) >= 0 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= 0 &
!               .and. shifted) then
!
!              if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + gout(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
!              else if (gout(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  gout(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
!                gout(a,b,c,i) = (gout(a+cx(i),b+cy(i),c+cz(i),i) + abs(gout(a,b,c,i)))/2
!              else if (gout(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
!                  gout(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                gout(a,b,c,i) = (gout(a-cx(i),b-cy(i),c-cz(i),i) + abs(gout(a,b,c,i)))/2
!              else
!                gout(a,b,c,i) = abs(gout(a,b,c,i))
!              end if

!
!
!            else if (state(a,b,c,i) >= 0 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= 0 &
!               .and. shifted) then
!
!              if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + gi(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
!              else if (gi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  gi(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
!                gi(a,b,c,i) = (gi(a+cx(i),b+cy(i),c+cz(i),i) + abs(gi(a,b,c,i)))/2
!              else if (gi(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
!                  gi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                gi(a,b,c,i) = (gi(a-cx(i),b-cy(i),c-cz(i),i) + abs(gi(a,b,c,i)))/2
!              else
!                gi(a,b,c,i) = abs(gi(a,b,c,i))
!              end if

!
!
!
!            else if (state(a,b,c,i) >= 0 .and. state(a-cx(i),b-cy(i),c-cz(i),1) >= 0 &
!               .and. shifted) then
!
!              if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + fi(a-cx(i),b-cy(i),c-cz(i),i))/2.0D0
!              else if (fi(a+cx(i),b+cy(i),c+cz(i),i) > 0.0D0 .and. &
!                  fi(a-cx(i),b-cy(i),c-cz(i),i) < 0.0D0) then
!                fi(a,b,c,i) = (fi(a+cx(i),b+cy(i),c+cz(i),i) + abs(fi(a,b,c,i)))/2
!              else if (fi(a+cx(i),b+cy(i),c+cz(i),i) < 0.0D0 .and. &
!                  fi(a-cx(i),b-cy(i),c-cz(i),i) > 0.0D0) then
!                fi(a,b,c,i) = (fi(a-cx(i),b-cy(i),c-cz(i),i) + abs(fi(a,b,c,i)))/2
!              else
!                fi(a,b,c,i) = abs(fi(a,b,c,i))
!              end if
!
