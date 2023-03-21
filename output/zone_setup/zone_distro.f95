subroutine zone_distro(lvl)
!
! Distribute all necessary zonal information to every processor which will be participating
! in output.
!
! Called by: lb_output
! Calls:
! External calls: mpi_allgatherv, mpi_barrier
!
use precise
use mpi
use amr_info_holder, only : zone_storage
use amr_processes , only :self,commune,nprocs
implicit none
integer :: lvl,i,ier
integer :: n_elems,counter
!  integer :: stat(MPI_STATUS_SIZE)
integer,allocatable :: disp_array(:),rec_array(:)

!write(*,*) 'entering zone information distribution function.',self

counter = 0

n_elems = 3*zone_storage(lvl)%big_zones

if (allocated(disp_array)) deallocate(disp_array)
allocate(disp_array(0:nprocs-1))
if (allocated(rec_array)) deallocate(rec_array)
allocate(rec_array(0:nprocs-1))
do i = 0,nprocs-1
  rec_array(i) = zone_storage(lvl)%little_zones(i)*3
end do
do i=0,nprocs-1
  if (i==0) then
    disp_array(i) = 0
  else
    disp_array(i) = disp_array(i-1)+zone_storage(lvl)%little_zones(i-1)*3
  end if
end do
!  write(*,*) 'receiving array',self,lvl,rec_array
!  write(*,*) 'displacement array',self,lvl,disp_array
!
! The allgatherv distribute all necessary information for zone writing to every processor
!
call mpi_allgatherv(zone_storage(lvl)%lower(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_INT,&
  zone_storage(lvl)%lower,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(zone_storage(lvl)%higher(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_INT,&
  zone_storage(lvl)%higher,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(zone_storage(lvl)%lo_corner(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_DOUBLE_PRECISION,&
  zone_storage(lvl)%lo_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)
!
call mpi_allgatherv(zone_storage(lvl)%hi_corner(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_DOUBLE_PRECISION,&
  zone_storage(lvl)%hi_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)

call mpi_barrier(commune,ier)

deallocate(rec_array)
deallocate(disp_array)

!  write(*,*) 'lower test 3',self,'lower',zone_storage(lvl)%lower
!  write(*,*) 'higher test 3',self,'high',zone_storage(lvl)%higher
!  write(*,*) 'lo_corner check ',self,'lo_corner',zone_storage(lvl)%lo_corner
!  write(*,*) 'hi_corner check ',self,'hi_corner',zone_storage(lvl)%hi_corner
!  write(*,*) 'big_zones check ',self,'lvl',lvl,'big_zones',zone_storage(lvl)%big_zones
!  write(*,*) 'little_zones check ',self,'little_zones ',zone_storage(lvl)%little_zones
!  write(*,*) 'start and end',self,'lvl',lvl,zone_storage(lvl)%zone_start,&
!    zone_storage(lvl)%zone_end
end subroutine
