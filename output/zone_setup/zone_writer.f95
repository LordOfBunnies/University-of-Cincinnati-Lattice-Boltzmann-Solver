subroutine zone_writer(min_lvl,max_lvl,sv_num)
!
! Zone writing for the time accurate, deforming grid
!
!
! Called by: save_write
! Calls: zone_coord_setup
! External Calls: cg_zone_write_f
!
use precise
use amrex_base_module
use mpi
use nml_output
use cgns
use amr_info_holder , only: zone_storage
use amr_processes, only: self,nprocs,commune
implicit none
integer :: fine,i,sv_num,zone_counter,index_file,index_base,ier,index_coord
integer :: max_lvl,lvl,min_lvl
integer(cgsize_t) :: isize(3,3)
character(len=32) :: zone_number,zonename,format_nm_lvl_and_zone

index_file = sv_num+1
index_base = 1
!
! Set up the zones for output
!
format_nm_lvl_and_zone = '(I2.2,a1,I3.3)'
!
!
!
zone_counter = 0
do lvl = min_lvl,max_lvl

  call zone_box_setup(lvl)

  do i=1,zone_storage(lvl)%big_zones
    zone_counter = zone_counter+1
!    write(*,*) 'Sing with me',self,zone_counter,i
    isize(1,1) = zone_storage(lvl)%higher(1,i)-zone_storage(lvl)%lower(1,i)+1
    isize(2,1) = zone_storage(lvl)%higher(2,i)-zone_storage(lvl)%lower(2,i)+1
    isize(3,1) = zone_storage(lvl)%higher(3,i)-zone_storage(lvl)%lower(3,i)+1
    isize(1,2) = isize(1,1)-1
    isize(2,2) = isize(2,1)-1
    isize(3,2) = isize(3,1)-1
    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0
!    write(*,*) 'isize check',self,isize
!    write(*,*) 'pre-flight check',self,zone_counter,i,zone_storage(lvl)%big_zones
    write(zone_number,format_nm_lvl_and_zone) lvl,'_',zone_counter
    zonename = 'zone_lvl_'//TRIM(zone_number)
!    write(*,*) 'writing zone,',zonename,index_file,index_base,i,zone_counter

    zone_storage(lvl)%output_number(i,sv_num) = zone_counter

!    write(*,*) 'zone output number',self,zone_storage(lvl)%output_number(i,sv_num),zone_counter
!    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
!      Structured,zone_storage(lvl)%output_number(i,sv_num),ier)
    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
      Structured,zone_counter,ier)
!  write(*,*) 'after writing zone',index_file,index_base,zone_counter,zone_storage(lvl)%output_number(i,sv_num)
!    index_coord = index_zone
    call zone_coord_setup(index_file,index_base,zone_counter)
!  write(*,*) index_file,index_base,zone_counter
!    call mpi_barrier(commune,ier)
  end do

end do




end subroutine
