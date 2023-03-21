subroutine zone_box_setup(lvl)
!
! Shares and retrieves all the necessary data for CGNS output
!
! Called by: initial_output,output_setup,save_write
! Calls: zone_ranges,zone_distro
! External calls: mpi_allgatherv,mpi_barrier,mpi_allreduce
!
use precise
use mpi
use amr_processes, only : count_big_zones,self,commune,nprocs
use amr_info_holder, only : zone_storage
use nml_output, only : num_saves
implicit none
integer :: lvl,loc_zones,ier,large_zones,max_zones,i,req,j
integer :: tag
integer,allocatable :: rec_array(:),disp_array(:)

! Get the number of zones for each processor
loc_zones = count_big_zones(lvl)
!write(*,*) 'inside ZBS, n_loc_box',self,lvl,loc_zones
! Sum those zones and make sure every processor has a copy so all the other
! objects can be properly allocated

CALL MPI_Allreduce(loc_zones,zone_storage(lvl)%big_zones,1,MPI_INTEGER, &
  MPI_SUM,commune,ier)
!
! Distribute the number of zones owned by a processor at a given level.
! Then agglomerate and broadcast so all procs know everything about everyone
! else's zones.
!
if (allocated(zone_storage(lvl)%little_zones)) deallocate(zone_storage(lvl)%little_zones)
allocate(zone_storage(lvl)%little_zones(0:nprocs-1))
!
  zone_storage(lvl)%little_zones(self) = loc_zones
!  write(*,*) 'local zones',self,'lvl',lvl,'n little zone',loc_zones,zone_storage(lvl)%little_zones
  if (allocated(rec_array)) deallocate(rec_array)
  allocate(rec_array(0:nprocs-1))
  if (allocated(disp_array)) deallocate(disp_array)
  allocate(disp_array(0:nprocs-1))
  rec_array = 1
  do i = 0,nprocs-1
    disp_array(i) = i
  end do
!  write(*,*) 'receving array',self,rec_array
!  write(*,*) 'displacement array',self,disp_array
  call mpi_allgatherv(zone_storage(lvl)%little_zones(self),1,&
    MPI_INT,zone_storage(lvl)%little_zones,&
    rec_array,disp_array,MPI_INT,commune,ier)
!
!
!
  call mpi_barrier(commune,ier)
  deallocate(rec_array)
  deallocate(disp_array)
!  write(*,*) 'Do you hear the people sing?'
!  write(*,*) 'things and stuff',self,lvl,zone_storage(lvl)%little_zones

  if (allocated(zone_storage(lvl)%zone_start)) deallocate(zone_storage(lvl)%zone_start)
  allocate(zone_storage(lvl)%zone_start(0:nprocs-1))

  if (allocated(zone_storage(lvl)%zone_end)) deallocate(zone_storage(lvl)%zone_end)
  allocate(zone_storage(lvl)%zone_end(0:nprocs-1))
!
! Loop to assign start and end values to the zones each processor owns
!
  do i = 0,nprocs-1
    if (i==0) then
      zone_storage(lvl)%zone_start(i) = 1
      zone_storage(lvl)%zone_end(i) = zone_storage(lvl)%little_zones(i)
    else
      zone_storage(lvl)%zone_start(i) = zone_storage(lvl)%zone_end(i-1)+1
      zone_storage(lvl)%zone_end(i) = zone_storage(lvl)%zone_start(i)+&
        zone_storage(lvl)%little_zones(i)-1
    end if
  end do
!  write(*,*) 'zone start and end',self,lvl,zone_storage(lvl)%zone_start,zone_storage(lvl)%zone_end
!  write(*,*) 'Singing the songs of angry men!'
!
! allreduce sums the local zones for each processor to find the total number of zones on a level
!

  if (allocated(zone_storage(lvl)%lower)) deallocate(zone_storage(lvl)%lower)
  if (allocated(zone_storage(lvl)%higher)) deallocate(zone_storage(lvl)%higher)
  if (allocated(zone_storage(lvl)%lo_corner)) deallocate(zone_storage(lvl)%lo_corner)
  if (allocated(zone_storage(lvl)%hi_corner)) deallocate(zone_storage(lvl)%hi_corner)
  if (allocated(zone_storage(lvl)%output_number)) deallocate(zone_storage(lvl)%output_number)
!
!  max_zones = zone_storage(lvl)%big_zones/nprocs+1
!
!  write(*,*) "I'm",self,'number of big zones',zone_storage(lvl)%big_zones,&
!    'at level ',lvl
  allocate(zone_storage(lvl)%lower(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%higher(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%lo_corner(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%hi_corner(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%output_number(zone_storage(lvl)%big_zones,num_saves))
!
!
!
  call zone_ranges(lvl,zone_storage(lvl)%big_zones,zone_storage(lvl)%lower,&
    zone_storage(lvl)%higher,zone_storage(lvl)%lo_corner,zone_storage(lvl)%hi_corner,&
    zone_storage(lvl)%zone_start,zone_storage(lvl)%zone_end)
!
!
   call zone_distro(lvl)
!
!  write(*,*) self,zone_storage(lvl)%big_zones
!  write(*,*) self,zone_storage(lvl)%little_zones
!  write(*,*) self,zone_storage(lvl)%zone_start,zone_storage(lvl)%zone_end
!  write(*,*) self,zone_storage(lvl)%lower
!  write(*,*) self,zone_storage(lvl)%higher
!  write(*,*) self,zone_storage(lvl)%lo_corner
!  write(*,*) self,zone_storage(lvl)%hi_corner
!  write(*,*) self,zone_storage(lvl)%zones_owned

  !  zone_storage(lvl)%lower,&
  !  zone_storage(lvl)%higher)
end subroutine















