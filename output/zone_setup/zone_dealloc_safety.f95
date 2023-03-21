subroutine zone_dealloc_safety(fine,num)
!
! Deallocate the zone storage arrays after an output
!
!
! Called by: lb_output
! Calls:
!
use precise
use nml_output
use amr_info_holder, only : zone_storage
!use amrex_amr_module
implicit none
integer :: lvl,fine,num



do lvl = save_shapes(2,num),save_shapes(1,num),-1
! for initial output cases
  if (num == 0) then
    if (allocated(zone_storage(lvl)%zone_end)) deallocate(zone_storage(lvl)%zone_end)
    if (allocated(zone_storage(lvl)%zone_start)) deallocate(zone_storage(lvl)%zone_start)
    if (allocated(zone_storage(lvl)%little_zones)) deallocate(zone_storage(lvl)%little_zones)

    if (allocated(zone_storage(lvl)%lower)) deallocate(zone_storage(lvl)%lower)
    if (allocated(zone_storage(lvl)%higher)) deallocate(zone_storage(lvl)%higher)
    if (allocated(zone_storage(lvl)%lo_corner)) deallocate(zone_storage(lvl)%lo_corner)
    if (allocated(zone_storage(lvl)%hi_corner)) deallocate(zone_storage(lvl)%hi_corner)
    if (allocated(zone_storage(lvl)%output_number)) deallocate(zone_storage(lvl)%output_number)
! out of range, to prevent segfaults
  else if (lvl > save_shapes(2,num) .or. lvl < save_shapes(1,num)) then
    cycle
! clear the saved information, XXXX unneeded for static mesh, fix later
  else
    write(*,*) 'deallocating for safety',lvl

    if (allocated(zone_storage(lvl)%zone_end)) deallocate(zone_storage(lvl)%zone_end)
    if (allocated(zone_storage(lvl)%zone_start)) deallocate(zone_storage(lvl)%zone_start)
    if (allocated(zone_storage(lvl)%little_zones)) deallocate(zone_storage(lvl)%little_zones)

    if (allocated(zone_storage(lvl)%lower)) deallocate(zone_storage(lvl)%lower)
    if (allocated(zone_storage(lvl)%higher)) deallocate(zone_storage(lvl)%higher)
    if (allocated(zone_storage(lvl)%lo_corner)) deallocate(zone_storage(lvl)%lo_corner)
    if (allocated(zone_storage(lvl)%hi_corner)) deallocate(zone_storage(lvl)%hi_corner)
    if (allocated(zone_storage(lvl)%output_number)) deallocate(zone_storage(lvl)%output_number)
!    write(*,*) 'sub-parts deallocated',lvl
  end if
end do

if (allocated(zone_storage)) deallocate(zone_storage)
!write(*,*) 'main user defined variable deallocated'
end subroutine
