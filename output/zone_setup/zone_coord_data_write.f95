subroutine zone_coord_data_write(index_file,index_base,lvl_lo,lvl_hi,sv_num)
!
! Write the coordinates for the zone
!
! Called by: output_setup,
! Calls:
! External calls: cgp_coord_data_write_f,cgp_error_exit_f
!
use precise
use cgns
use amrex_amr_module
use amrex_base_module
use grid_data, only: rdx
use amr_info_holder, only: zone_storage
use amr_processes, only: self,nprocs,commune
use output_data, only: x_verts,y_verts,z_verts
implicit none
integer ::a,b,c,lvl,out_min,out_max,lvl_hi,lvl_lo,i,ier
integer :: index_file,index_base,index_zone,sv_num
integer :: coord_x,coord_y,coord_z
integer(cgsize_t) ::isize(3,3),loc_hi(3),loc_lo(3)

coord_x = 1
coord_y = 2
coord_z = 3

!write(*,*) 'inside zone coordinate data write subroutine',self,lvl_lo,lvl_hi
!write(*,*) 'amrex limits',amrex_problo,amrex_probhi

do lvl = lvl_lo,lvl_hi

  !write(*,*) 'amrex dx',amrex_geom(lvl)%dx(1),amrex_geom(lvl)%dx(2),amrex_geom(lvl)%dx(3)

  do i = zone_storage(lvl)%zone_start(self),zone_storage(lvl)%zone_end(self)

    isize(1,1) = zone_storage(lvl)%higher(1,i)-zone_storage(lvl)%lower(1,i)+1
    isize(2,1) = zone_storage(lvl)%higher(2,i)-zone_storage(lvl)%lower(2,i)+1
    isize(3,1) = zone_storage(lvl)%higher(3,i)-zone_storage(lvl)%lower(3,i)+1

    isize(1,2) = isize(1,1)-1
    isize(2,2) = isize(2,1)-1
    isize(3,2) = isize(3,1)-1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

!  write(*,*) 'isize array',isize
!    write(*,*) 'allocating vertices'
    allocate(x_verts(isize(1,1),isize(2,1),isize(3,1)))
    allocate(y_verts(isize(1,1),isize(2,1),isize(3,1)))
    allocate(z_verts(isize(1,1),isize(2,1),isize(3,1)))

!    WRITE(*,*) 'Hi Im processor ', self,' and my grid range is ',byx_lo(1),byx_lo(2),&
!      byx_lo(3),byx_hi(1),byx_hi(2),byx_hi(3),index_zone
    do c = 1,isize(3,1)
      do b = 1,isize(2,1)
        do a = 1,isize(1,1)
          x_verts(a,b,c) = zone_storage(lvl)%lo_corner(1,i) + (a-1)*rdx(lvl)
          y_verts(a,b,c) = zone_storage(lvl)%lo_corner(2,i) + (b-1)*rdx(lvl)
          z_verts(a,b,c) = zone_storage(lvl)%lo_corner(3,i) + (c-1)*rdx(lvl)
        end do
      end do
    end do

    loc_lo = (/1,1,1/)
    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)

    index_zone = zone_storage(lvl)%output_number(i,sv_num)
!    write(*,*) 'repeat index info',index_file,index_base,index_zone
!    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!
!    write(*,*) 'box corners',zone_storage(lvl)%lo_corner(1:3,i),x_verts(isize(1,1),isize(2,1),isize(3,1)),&
!      y_verts(isize(1,1),isize(2,1),isize(3,1)),z_verts(isize(1,1),isize(2,1),isize(3,1)),lvl
    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_x,&
      loc_lo,loc_hi,x_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'reading?'
!    write(*,*) 'x coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
    !CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,byx_lo,byx_hi,y_verts,ier)
    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,&
      loc_lo,loc_hi,y_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'y coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_z,&
      loc_lo,loc_hi, z_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'z coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    WRITE(*,*) 'Coordinate value write completed.'

    deallocate(x_verts)
    deallocate(y_verts)
    deallocate(z_verts)

  end do
end do

end subroutine
