subroutine nullify_nodes(fine,num)
!
! The goal is to make the nodes above finer grid into a null node.
! i.e. a node with a state value of -666.
!
!
! Called by: fluid_or_solid
! Calls:
! External calls: mpi_reduce,mpi_allgatherv,mpi_bcast
!
use precise
use constants
use linkwise
use grid_data
!use amrex_base_module
use amrex_amr_module
use amr_processes
use amr_info_holder,only: state,mfstate,grid_info,nghosts_state,nghosts,nghosts_mid
use mpi
implicit none
integer :: a,b,c,lvl,gr,ier,i,num
integer :: n_boxes,counter,fine,rec_buff(0:nprocs-1),disp_buff(0:nprocs-1)
integer :: over_lo(3),over_hi(3),rec_array(0:nprocs-1),disp_array(0:nprocs-1)
integer :: convert_counter
logical :: inside_lo,inside_hi,overlap,dir_safe(dir)


type(amrex_mfiter) :: mfi
type(amrex_box) :: rox

rec_buff = 1
do i = 0,nprocs-1
  disp_buff(i) = i
end do

!write(*,*) 'nullifying nodes'

!if (allocated(grid_info)) deallocate(grid_info)
!allocate(grid_info(0:amrex_max_level))

!
!
! XXXX fix this later, take the boxes and make everything within the box null or
! solid
!
!
do lvl = 0,fine-1
  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

!  convert_counter = 0

  do while (mfi%next())
    state => mfstate(lvl)%dataptr(mfi)
    rox = mfi%validbox()

    do gr = 1,grid_info(lvl+1,num)%big_zones
      !See if the coarser box contains the finer box
      overlap = .false.

      call ml_overlapping_boxes(rox%lo-nghosts_state,grid_info(lvl+1,num)%lower(1:3,gr),rox%hi+nghosts_state,&
        grid_info(lvl+1,num)%higher(1:3,gr),over_lo,over_hi,overlap,lvl,lvl+1)

      if (overlap) then
!
! This is for the overlapped region ONLY.  Later parts must do everything because of several
!   boxes can overlap
!
        do c = over_lo(3),over_hi(3)
          do b = over_lo(2),over_hi(2)
            do a = over_lo(1),over_hi(1)

              if (state(a,b,c,1) >= 0) then
                state(a,b,c,1) = -666
!                convert_counter = convert_counter+1
              end if

            end do
          end do
        end do

      else
        cycle
      end if !overlapping regions conditional

    end do



  end do
!  write(*,*) 'Number of nodes converted to null',convert_counter,lvl


  call amrex_mfiter_destroy(mfi)

end do

do lvl = 0,fine-1
  call simple_state_directions(lvl)
end do
!
!
! This is to catch coarse nodes at the edge of finer boxes.  It detects a nodes adjacent to
! valid fluid and turns it into a temporary semi-null value.  It then converts all semi-null
! values to proper fluid values.
!
!
do lvl = 0,fine-1
  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
  convert_counter = 0
  do while (mfi%next())
    state => mfstate(lvl)%dataptr(mfi)
    rox = mfi%validbox()

      !See if the coarser box contains the finer box
      !inside_lo = magic_box(rox_lo,rox_hi,grid_info(lvl+1)%lower(1,gr))
      !inside_hi = magic_box(rox_lo,rox_hi,grid_info(lvl+1)%lhigher(1,gr))

      do c = rox%lo(3)-nghosts_state,rox%hi(3)+nghosts_state
        do b = rox%lo(2)-nghosts_state,rox%hi(2)+nghosts_state
          do a = rox%lo(1)-nghosts_state,rox%hi(1)+nghosts_state

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

            if (state(a,b,c,1) == -666) then
              dir_safe = direction_tester(a,b,c,rox%lo,rox%hi)

              do i = 1,dir
                if (dir_safe(i)) then
                  if (state(a+cx(i),b+cy(i),c+cz(i),1) >= 0) then
                    state(a,b,c,1) = -667
                    exit
                  end if
                end if
              end do

            end if
          end do
        end do
      end do

!      write(*,*) 'nodes converted to semi-null',convert_counter,lvl
!
! Valid nodes on the coarse level, at the edge of the valid fine nodes.  These carry information
!   from the fine nodes to the coarse nodes.
!
      convert_counter = 0
      do c = rox%lo(3)-nghosts_state,rox%hi(3)+nghosts_state
        do b = rox%lo(2)-nghosts_state,rox%hi(2)+nghosts_state
          do a = rox%lo(1)-nghosts_state,rox%hi(1)+nghosts_state
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
!              convert_counter = convert_counter+1
            end if

!            if (a == 113 .and. b == 49 .and. c == 49 .and. lvl == 2) then
!              write(*,*) 'node status, nullify',state(a,b,c,1),a,b,c,lvl
!            end if

          end do
        end do
      end do


    end do



  call amrex_mfiter_destroy(mfi)

end do


end subroutine


!    do c = rox%lo(3),rox%hi(3)
!      do b = rox%lo(2),rox%hi(2)
!        do a = rox%lo(1),rox%hi(1)
!          inside = .false.
!          do gr = grid_info(lvl,num)%zone_start(self),grid_info(lvl,num)%zone_end(self)
!            inside = magic_box(grid_info(lvl,num)%lower(1:3,gr),&
!              grid_info(lvl,num)%higher(1:3,gr), &
!              (/2*a,2*b,2*c/))
!            if (inside .and. state(a,b,c,1)>=-1001) then
!              state(a,b,c,1) = -666
!            end if
!
!          end do
!        end do
!      end do
!    end do

!      write(*,*) 'nodes converted back to fluid',convert_counter,lvl

!     convert_counter = 0
!      do c = rox%lo(3)-nghosts,rox%hi(3)+nghosts
!        do b = rox%lo(2)-nghosts,rox%hi(2)+nghosts
!          do a = rox%lo(1)-nghosts,rox%hi(1)+nghosts
!            do i = 1,dir
!              if (a+cx(i) < rox%lo(1)-nghosts .or. b+cy(i) < rox%lo(2)-nghosts .or. &
!                  c+cz(i) < rox%lo(3)-nghosts .or. a+cx(i) > rox%hi(1)+nghosts .or. &
!                  b+cy(i) > rox%hi(2)+nghosts .or. c+cz(i) > rox%hi(3)+nghosts) then
!                if () then
!
!                end if
!              else
!                state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
!                convert_counter = convert_counter+1
!              end if
!
!          end do
!        end do
!      end do



!if (state(a,b,c,2) >= 0 .or. state(a,b,c,16) >= 0 .or. &
!                  state(a,b,c,34) >= 0) then
!                if (state(a,b,c,2) < 1000 .or. state(a,b,c,16) < 1000 .or. &
!                  state(a,b,c,34) < 1000) then
!!
!                  state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                  cycle
!                end if
!              end if
!            !end if
!
!            !if (b+cy(3) > rox%lo(2)-nghosts) then
!              if (state(a,b,c,3) >= 0 .or. state(a,b,c,17) >= 0 .or. &
!                  state(a,b,c,35) >= 0) then
!                if (state(a,b,c,3) < 1000 .or. state(a,b,c,17) < 1000 .or. &
!                  state(a,b,c,35) < 1000) then
!
!                  state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                  cycle
!                end if
!              end if
!            !end if
!
!            !if (c+cz(4) > rox%lo(3)-nghosts) then
!              if (state(a,b,c,4) >= 0 .or. state(a,b,c,18) >= 0 .or. &
!                  state(a,b,c,36) >= 0) then
!                if (state(a,b,c,4) < 1000 .or. state(a,b,c,18) < 1000 .or. &
!                  state(a,b,c,36) < 1000) then
!
!                  state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                  cycle
!                end if
!              end if
!            !end if
!
!            !if (a+cx(5) < rox%hi(1)+nghosts) then
!              if (state(a,b,c,5) >= 0 .or. state(a,b,c,19) >= 0 .or. &
!                  state(a,b,c,37) >= 0) then
!                if (state(a,b,c,5) < 1000 .or. state(a,b,c,19) < 1000 .or. &
!                  state(a,b,c,37) < 1000) then
!
!                  state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                  cycle
!                end if
!              end if
!            !end if
!
!            !if (b+cy(6) < rox%hi(2)+nghosts) then
!              if (state(a,b,c,6) >= 0 .or. state(a,b,c,20) >= 0 .or. &
!                  state(a,b,c,38) >= 0) then
!                 if (state(a,b,c,6) < 1000 .or. state(a,b,c,20) < 1000 .or. &
!                  state(a,b,c,38) < 1000) then
!
!                  state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                  cycle
!                end if
!              end if
!            !end if
!
!            !if (c+cz(7) < rox%hi(3)+nghosts) then
!              if (state(a,b,c,7) >= 0 .or. state(a,b,c,21) >= 0 .or. &
!                  state(a,b,c,39) >= 0) then
!                 if (state(a,b,c,7) < 1000 .or. state(a,b,c,21) < 1000 .or. &
!                  state(a,b,c,39) < 1000) then
!
!                  state(a,b,c,1) = -667
!!                convert_counter = convert_counter+1
!                  cycle
!                end if
!              end if


!
!              if (a+cx(i) < rox%lo(1)-nghosts .or. a+cx(i) > rox%hi(1)+nghosts &
!                .or. b+cy(i) <  rox%lo(2)-nghosts.or. b+cy(i) > rox%hi(2)+nghosts &
!                .or. c+cz(i) <  rox%lo(3)-nghosts.or. c+cz(i) > rox%hi(3)+nghosts) then
!
!                if (a+cx(i)+1 < rox%lo(1)-nghosts .or. a+cx(i)-1 > rox%hi(1)+nghosts &
!                  .or. b+cy(i)+1 <  rox%lo(2)-nghosts.or. b+cy(i)-1 > rox%hi(2)+nghosts &
!                  .or. c+cz(i)+1 <  rox%lo(3)-nghosts.or. c+cz(i)-1 > rox%hi(3)+nghosts) then
!
!                  if (a+cx(i)+2 < rox%lo(1)-nghosts .or. a+cx(i)-\2 > rox%hi(1)+nghosts &
!                    .or. b+cy(i)+2 <  rox%lo(2)-nghosts.or. b+cy(i)-2 > rox%hi(2)+nghosts &
!                    .or. c+cz(i)+2 <  rox%lo(3)-nghosts.or. c+cz(i)-2 > rox%hi(3)+nghosts) then
!
!                  else if () then
!
!                  else if () then
!
!
!                  end if
!
!                  else if () then
!
!                  else if () then
!                end if
!
!                else if () then
!              end if


!
            !end if





!do lvl = 0,fine
!  n_boxes = count_big_zones(lvl)
!  CALL MPI_Allreduce(n_boxes,grid_info(lvl,num)%big_zones,1,MPI_INTEGER, &
!  MPI_SUM,commune,ier)
!  if (allocated(grid_info(lvl,num)%lower)) deallocate(grid_info(lvl,num)%lower)
!  if (allocated(grid_info(lvl,num)%higher)) deallocate(grid_info(lvl,num)%higher)
!
!  allocate(grid_info(lvl,num)%lower(3,grid_info(lvl,num)%big_zones))
!  allocate(grid_info(lvl,num)%higher(3,grid_info(lvl,num)%big_zones))
!
!  if (allocated(grid_info(lvl,num)%little_zones)) deallocate(grid_info(lvl,num)%little_zones)
!  allocate(grid_info(lvl,num)%little_zones(0:nprocs-1))
!!
!  grid_info(lvl,num)%little_zones(self) = n_boxes
!!  write(*,*) 'off with their heads!',self,lvl,n_boxes,grid_info(lvl,num)%little_zones(self)
!!  write(*,*) 'local zones',self,'lvl',lvl,'n little zone',loc_zones,grid_info(lvl,num)%little_zones
!!  if (allocated(rec_array)) deallocate(rec_array)
!!  allocate(rec_array(0:nprocs-1))
!!  if (allocated(disp_array)) deallocate(disp_array)
!!  allocate(disp_array(0:nprocs-1))
!!  rec_array = 1
!!
!!  do i = 0,nprocs-1
!!    disp_array(i) = i
!!  end do
!!  write(*,*) 'receving array',self,rec_array
!!  write(*,*) 'displacement array',self,disp_array
!  call mpi_allgatherv(grid_info(lvl,num)%little_zones(self),1,&
!    MPI_INT,grid_info(lvl,num)%little_zones,rec_buff(0:nprocs-1),disp_buff(0:nprocs-1),&
!    MPI_INT,commune,ier)
!  call mpi_barrier(commune,ier)
!  if (allocated(grid_info(lvl,num)%zone_start)) deallocate(grid_info(lvl,num)%zone_start)
!  allocate(grid_info(lvl,num)%zone_start(0:nprocs-1))
!
!  if (allocated(grid_info(lvl,num)%zone_end)) deallocate(grid_info(lvl,num)%zone_end)
!  allocate(grid_info(lvl,num)%zone_end(0:nprocs-1))
!
!!  write(*,*) 'bourgoise',self,lvl,grid_info(lvl,num)%little_zones
!!
!!
!!
!  do i = 0,nprocs-1
!    if (i==0) then
!      grid_info(lvl,num)%zone_start(i) = 1
!      grid_info(lvl,num)%zone_end(i) = grid_info(lvl,num)%little_zones(i)
!    else
!      grid_info(lvl,num)%zone_start(i) = grid_info(lvl,num)%zone_end(i-1)+1
!      grid_info(lvl,num)%zone_end(i) = grid_info(lvl,num)%zone_start(i)+&
!        grid_info(lvl,num)%little_zones(i)-1
!    end if
!  end do
!!  write(*,*) 'casse-toi',self,lvl,n_boxes,grid_info(lvl,num)%zone_start,grid_info(lvl,num)%zone_end
!counter = grid_info(lvl,num)%zone_start(self)
!  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
!  do while (mfi%next())
!    rox = mfi%validbox()
!    grid_info(lvl,num)%lower(1,counter) = rox%lo(1)
!    grid_info(lvl,num)%higher(1,counter) = rox%hi(1)
!    grid_info(lvl,num)%lower(2,counter) = rox%lo(2)
!    grid_info(lvl,num)%higher(2,counter) = rox%hi(2)
!    grid_info(lvl,num)%lower(3,counter) = rox%lo(3)
!    grid_info(lvl,num)%higher(3,counter) = rox%hi(3)
!    counter = counter +1
!  end do
!  call amrex_mfiter_destroy(mfi)
!
!do i = 0,nprocs-1
!  rec_array(i) = grid_info(lvl,num)%little_zones(i)*3
!end do
!do i=0,nprocs-1
!  if (i==0) then
!    disp_array(i) = 0
!  else
!    disp_array(i) = disp_array(i-1)+grid_info(lvl,num)%little_zones(i-1)*3
!  end if
!end do
!!
!call mpi_allgatherv(grid_info(lvl,num)%lower(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%little_zones(self),MPI_INT,&
!  grid_info(lvl,num)%lower,rec_array,disp_array,&
!  MPI_INT,commune,ier)
!!
!call mpi_allgatherv(grid_info(lvl,num)%higher(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%little_zones(self),MPI_INT,&
!  grid_info(lvl,num)%higher,rec_array,disp_array,&
!  MPI_INT,commune,ier)
!!
!
!
!end do
