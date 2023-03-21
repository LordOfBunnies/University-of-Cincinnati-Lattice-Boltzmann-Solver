subroutine find_overlapped_coarse(lvl)
!
!
!
!
!
!
!
!
!
use precise
use constants
use grid_data, only: dir
use mpi
use linkwise
use amrex_amr_module
use amrex_base_module
use amr_info_holder, only: state,mfstate,grid_info,nghosts,nghosts_state,nghosts_mid
use amr_processes
implicit none
integer,intent(in) :: lvl
integer :: a,b,c,i
logical :: dir_safe(dir)

type(amrex_box) :: rox
type(amrex_mfiter) :: mfi

call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
!convert_counter = 0
do while (mfi%next())
  state => mfstate(lvl)%dataptr(mfi)
  rox = mfi%validbox()

    !See if the coarser box contains the finer box
    !inside_lo = magic_box(rox_lo,rox_hi,grid_info(lvl+1)%lower(1,gr))
    !inside_hi = magic_box(rox_lo,rox_hi,grid_info(lvl+1)%lhigher(1,gr))

    do c = rox%lo(3)-nghosts_state,rox%hi(3)+nghosts_state
      do b = rox%lo(2)-nghosts_state,rox%hi(2)+nghosts_state
        do a = rox%lo(1)-nghosts_state,rox%hi(1)+nghosts_state

!          if (a < amrex_geom(lvl)%domain%lo(1) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
!              b < amrex_geom(lvl)%domain%lo(2) .or. b > amrex_geom(lvl)%domain%hi(2)+1 .or. &
!              c < amrex_geom(lvl)%domain%lo(3) .or. c > amrex_geom(lvl)%domain%hi(3)+1) then
!            cycle
!          end if
          if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
              c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
              b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then

            if (state(a,b,c,1) >= -299 .and. state(a,b,c,1) <-200) then
              cycle
            else
              state(a,b,c,1) = -666
              cycle
            end if
          end if
          !write(*,*) 'vert',state(a,b,c,1),lvl,a,b,c

          !write(*,*) 'coordinates',a,b,c,lvl
          if (state(a,b,c,1) == -666) then
            dir_safe = direction_tester(a,b,c,rox%lo,rox%hi)
             !write(*,*) 'vert',dir_safe,lvl,a,b,c
            do i = 1,dir
              if (dir_safe(i)) then

                if (state(a+cx(i),b+cy(i),c+cz(i),1) >= 0) then
                  state(a,b,c,1) = -667
                  !write(*,*) 'find it!',a,b,c
!                  if (a == 108 .and. b == 72 .and. c == 35) then
!                    write(*,*) 'Master flash',state(a,b,c,1:39),rox%lo,rox%hi,i,a,b,c,lvl,self
!                  end if
                  exit
                end if
              end if
            end do
          end if
        end do
      end do
    end do
!write(*,*) 'cranky baby',self
!      write(*,*) 'nodes converted to semi-null',convert_counter,lvl
!
!
!
  do c = rox%lo(3)-nghosts_state,rox%hi(3)+nghosts_state
    do b = rox%lo(2)-nghosts_state,rox%hi(2)+nghosts_state
      do a = rox%lo(1)-nghosts_state,rox%hi(1)+nghosts_state
! Out of bounds automatically disqualified, don't overwrite outflow nodes
        if (a < amrex_geom(lvl)%domain%lo(1) .or. b < amrex_geom(lvl)%domain%lo(2) .or. &
              c < amrex_geom(lvl)%domain%lo(3) .or. a > amrex_geom(lvl)%domain%hi(1)+1 .or. &
              b > amrex_geom(lvl)%domain%hi(2)+1 .or. c > amrex_geom(lvl)%domain%hi(3)+1 ) then

          if (state(a,b,c,1) >= -299 .and. state(a,b,c,1) <-200) then
            cycle
          else
            state(a,b,c,1) = -666
            cycle
          end if
        end if

        if (state(a,b,c,1) == -667) then
          state(a,b,c,1) = (lvl+51)*1000
!          write(*,*) 'bingo',a,b,c
!              convert_counter = convert_counter+1
        end if
      end do
    end do
  end do
!  write(*,*) 'egg sammich',self
end do

call amrex_mfiter_destroy(mfi)

end subroutine

!          !if (a+cx(2) > rox%lo(1)-nghosts) then
!            if (state(a,b,c,2) >= 0 .or. state(a,b,c,16) >= 0 .or. &
!                state(a,b,c,34) >= 0) then
!              if (state(a,b,c,2) < 1000 .or. state(a,b,c,16) < 1000 .or. &
!                state(a,b,c,34) < 1000) then
!!
!                state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                cycle
!              end if
!            end if
!          !end if
!
!          !if (b+cy(3) > rox%lo(2)-nghosts) then
!            if (state(a,b,c,3) >= 0 .or. state(a,b,c,17) >= 0 .or. &
!                state(a,b,c,35) >= 0) then
!              if (state(a,b,c,3) < 1000 .or. state(a,b,c,17) < 1000 .or. &
!                state(a,b,c,35) < 1000) then
!!
!                state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                cycle
!              end if
!            end if
!          !end if
!
!          !if (c+cz(4) > rox%lo(3)-nghosts) then
!            if (state(a,b,c,4) >= 0 .or. state(a,b,c,18) >= 0 .or. &
!                state(a,b,c,36) >= 0) then
!              if (state(a,b,c,4) < 1000 .or. state(a,b,c,18) < 1000 .or. &
!                state(a,b,c,36) < 1000) then
!!
!                state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                cycle
!              end if
!            end if
!          !end if
!
!          !if (a+cx(5) < rox%hi(1)+nghosts) then
!            if (state(a,b,c,5) >= 0 .or. state(a,b,c,19) >= 0 .or. &
!                state(a,b,c,37) >= 0) then
!              if (state(a,b,c,5) < 1000 .or. state(a,b,c,19) < 1000 .or. &
!                state(a,b,c,37) < 1000) then
!!
!                state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                cycle
!              end if
!            end if
!          !end if
!
!          !if (b+cy(6) < rox%hi(2)+nghosts) then
!            if (state(a,b,c,6) >= 0 .or. state(a,b,c,20) >= 0 .or. &
!                state(a,b,c,38) >= 0) then
!              if (state(a,b,c,6) < 1000 .or. state(a,b,c,20) < 1000 .or. &
!                state(a,b,c,38) < 1000) then
!
!                state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                cycle
!              end if
!            end if
!          !end if
!
!          !if (c+cz(7) < rox%hi(3)+nghosts) then
!            if (state(a,b,c,7) >= 0 .or. state(a,b,c,21) >= 0 .or. &
!                state(a,b,c,39) >= 0) then
!              if (state(a,b,c,7) < 1000 .or. state(a,b,c,21) < 1000 .or. &
!                state(a,b,c,39) < 1000) then
!
!                state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                cycle
!              end if
!            end if
          !end if
!            if (a == 113 .and. b == 49 .and. c == 49 .and. lvl == 2) then
!              write(*,*) 'node status, nullify',state(a,b,c,1),a,b,c,lvl
!            end if
