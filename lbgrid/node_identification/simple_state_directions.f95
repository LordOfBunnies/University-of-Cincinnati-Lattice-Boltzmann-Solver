subroutine simple_state_directions(lvl)
!
! Finalize many of the requisite directions after find coarse/ghost/etc. nodes
!
!
!
! Called by: grid_gen
! Calls:
!
use precise
use amr_info_holder, only: state,mfstate,nghosts,nghosts_state,nghosts_mid
use amr_processes, only: self,commune
use amrex_amr_module
use grid_data
use linkwise
use mpi
implicit none
integer :: lvl,a,b,c,i,ier

type(amrex_mfiter) :: mfi
type(amrex_box) :: ducks

!write(*,*) 'simple state directions'

call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
do while (mfi%next())

  state => mfstate(lvl)%dataptr(mfi)
  ducks = mfi%validbox()
!
!
!
  do i = 2,dir

    do c = ducks%lo(3)-nghosts_state,ducks%hi(3)+nghosts_state
      do b = ducks%lo(2)-nghosts_state,ducks%hi(2)+nghosts_state
        do a = ducks%lo(1)-nghosts_state,ducks%hi(1)+nghosts_state
!          if (a == 1 .and. b == 32 .and. c == 32) then
!            write(*,*) 'ARGH!!!',state(a,b,c,1:dir),i,ducks%lo,ducks%hi,a,b,c,self
!          end if
!          if (a == 268 .and. b == 22 .and. c == 176) then
!            write(*,*) 'ARGH!!!',state(a,b,c,1:dir),i,ducks%lo,ducks%hi,a,b,c,self
!          end if

          if (a+cx(i) < ducks%lo(1) - nghosts_mid .or. a+cx(i) > ducks%hi(1) + nghosts_mid .or. &
              b+cy(i) < ducks%lo(2) - nghosts_mid .or. b+cy(i) > ducks%hi(2) + nghosts_mid .or. &
              c+cz(i) < ducks%lo(3) - nghosts_mid .or. c+cz(i) > ducks%hi(3) + nghosts_mid) then
            state(a,b,c,i) = -666
          else if (state(a,b,c,i) == -1001 .and. state(a,b,c,1) >= 0 .and. shifted) then
            state(a,b,c,i) = -666
          else if (state(a,b,c,i) < -1000 .and. state(a,b,c,1) >= 0) then
            cycle
          else if (state(a,b,c,i) < 0 .and. state(a,b,c,1) >= 0) then
            cycle
          else if (state(a,b,c,1) >= 0 .and. state(a+cx(i),b+cy(i),c+cz(i),1) >=0) then
            state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
          else if (state(a,b,c,1) == -666 .and. state(a+cx(i),b+cy(i),c+cz(i),1) >=0) then
            state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
            !write(*,*) 'Argo fuck yourself'

          else if (state(a,b,c,1) == -666) then
            state(a,b,c,i) = -666
!            if (a == 40 .and. b == 31 .and. c == 16) write(*,*) 'puffball'

          end if



        end do
      end do
    end do
  end do
end do
call amrex_mfiter_destroy(mfi)

call mpi_barrier(commune,ier)

end subroutine

!          if (a == 58 .and. b == 67 .and. c == 5 .and. i == 39) then
!            write(*,*) 'ARGH!!!',state(a,b,c,1:39),state(a+cx(i),b+cy(i),c+cz(i),1),&
!              ducks%lo,ducks%hi,a,b,c,lvl,self!,ducks%lo,ducks%hi,self
!          end if
!          if (a == 108 .and. b == 72 .and. c == 32 .and. i == 39) then
!            write(*,*) 'jingle hop',state(a,b,c,1:39),state(a+cx(i),b+cy(i),c+cz(i),1),&
!              ducks%lo,ducks%hi,a,b,c,lvl,self!,ducks%lo,ducks%hi,self
!            !write(*,*) 'cx test',cx,a,b,c
!          end if
!      if (a == 226 .and. b == 142 .and. c == 131) then
!        write(*,*) 'idiot check',state(a,b,c,1:39),ducks%lo,ducks%hi,i,a,b,c,lvl,self
!      end if

!XXXX search state(a+cx(i),b+cy(i),c+cz(i),i)
