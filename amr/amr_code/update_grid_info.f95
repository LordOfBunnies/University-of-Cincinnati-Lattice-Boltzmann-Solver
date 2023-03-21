subroutine update_grid_info(lvl,num)
!
! This updates the grid info user-defined variable which hold information on
!   all the boxes and grids.  This will be used to short circuit some issues
!   with regridding so every node after a regrid doesn't have to be tested
!   again.
!
! Called by: regrid_question_mark
! Calls:
! External calls:
!
use precise
use mpi
use amrex_amr_module
use amrex_base_module
use amr_info_holder, only: grid_info,mfstate
use amr_processes
implicit none
integer :: lvl,lvl_zones
integer :: rec_array(0:nprocs-1),disp_array(0:nprocs-1),loc_zones
integer :: i,counter,ier,num

type(amrex_box) :: rox
type(amrex_mfiter) :: mfi
!
!
!
loc_zones = count_big_zones(lvl)
lvl_zones = total_big_zones(lvl)

CALL MPI_Allreduce(loc_zones,grid_info(lvl,num)%big_zones,1,MPI_INTEGER, &
  MPI_SUM,commune,ier)
!write(*,*) 'merde',self,loc_zones,lvl_zones
!
!
! Change all of these to zone_storage since we have to refresh all that anyway.
!
!
if (allocated(grid_info(lvl,num)%lower)) deallocate(grid_info(lvl,num)%lower)
if (allocated(grid_info(lvl,num)%higher)) deallocate(grid_info(lvl,num)%higher)
if (allocated(grid_info(lvl,num)%zone_start)) deallocate(grid_info(lvl,num)%zone_start)
if (allocated(grid_info(lvl,num)%zone_end)) deallocate(grid_info(lvl,num)%zone_end)

allocate(grid_info(lvl,num)%lower(3,lvl_zones))
allocate(grid_info(lvl,num)%higher(3,lvl_zones))
allocate(grid_info(lvl,num)%zone_start(0:nprocs-1))
allocate(grid_info(lvl,num)%zone_end(0:nprocs-1))
!
!
!
if (allocated(grid_info(lvl,num)%little_zones)) deallocate(grid_info(lvl,num)%little_zones)
allocate(grid_info(lvl,num)%little_zones(0:nprocs-1))
!
grid_info(lvl,num)%little_zones(self) = loc_zones
!
!if (allocated(rec_array)) deallocate(rec_array)
!allocate(rec_array(0:nprocs-1))
!if (allocated(disp_array)) deallocate(disp_array)
!allocate(disp_array(0:nprocs-1))
!
rec_array = 1
do i = 0,nprocs-1
  disp_array(i) = i
end do
!
call mpi_allgatherv(grid_info(lvl,num)%little_zones(self),1,&
  MPI_INT,grid_info(lvl,num)%little_zones(0:nprocs-1),&
  rec_array,disp_array,MPI_INT,commune,ier)
!
!
!
  do i = 0,nprocs-1
    if (i==0) then
      grid_info(lvl,num)%zone_start(i) = 1
      grid_info(lvl,num)%zone_end(i) = grid_info(lvl,num)%little_zones(i)
    else
      grid_info(lvl,num)%zone_start(i) = grid_info(lvl,num)%zone_end(i-1)+1
      grid_info(lvl,num)%zone_end(i) = grid_info(lvl,num)%zone_start(i)+&
        grid_info(lvl,num)%little_zones(i)-1
    end if
  end do
!  write(*,*) 'casse-toi',self,lvl,grid_info(lvl,num)%zone_start,grid_info(lvl,num)%zone_end

counter = grid_info(lvl,num)%zone_start(self)
  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
  do while (mfi%next())
    rox = mfi%validbox()
    grid_info(lvl,num)%lower(1,counter) = rox%lo(1)
    grid_info(lvl,num)%higher(1,counter) = rox%hi(1)
    grid_info(lvl,num)%lower(2,counter) = rox%lo(2)
    grid_info(lvl,num)%higher(2,counter) = rox%hi(2)
    grid_info(lvl,num)%lower(3,counter) = rox%lo(3)
    grid_info(lvl,num)%higher(3,counter) = rox%hi(3)
    counter = counter +1
  end do
  call amrex_mfiter_destroy(mfi)

do i = 0,nprocs-1
  rec_array(i) = grid_info(lvl,num)%little_zones(i)*3
end do
do i=0,nprocs-1
  if (i==0) then
    disp_array(i) = 0
  else
    disp_array(i) = disp_array(i-1)+grid_info(lvl,num)%little_zones(i-1)*3
  end if
end do
!
call mpi_allgatherv(grid_info(lvl,num)%lower(1:3,grid_info(lvl,num)%zone_start(self):&
  grid_info(lvl,num)%zone_end(self)),rec_array(self),MPI_INT,&
  grid_info(lvl,num)%lower(1:3,1:lvl_zones),rec_array(0:nprocs-1),disp_array(0:nprocs-1),&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(grid_info(lvl,num)%higher(1:3,grid_info(lvl,num)%zone_start(self):&
  grid_info(lvl,num)%zone_end(self)),rec_array(self),MPI_INT,&
  grid_info(lvl,num)%higher(1:3,1:lvl_zones),rec_array(0:nprocs-1),disp_array(0:nprocs-1),&
  MPI_INT,commune,ier)
!write(*,*) 'canard!',grid_info(lvl,num)%little_zones,loc_zones
!write(*,*) 'canard!',grid_info(lvl,num)%lower,grid_info(lvl,num)%higher, lvl,self
end subroutine
