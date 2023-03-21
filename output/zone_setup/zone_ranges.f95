subroutine zone_ranges(lvl,num_zones,lo_nodes,hi_nodes,min_corner,&
  max_corner,zone_starts,zone_ends)
!
! This gets the lowest corner and higher corner for all the zones to be called in
!
!
!
!
use precise
use mpi
use amr_processes, only : get_real_coords,self,nprocs
use amr_info_holder, only: mfrho
use amrex_base_module
use amrex_amr_module
!use mpi
implicit none
integer :: counter,lvl,num_zones,ind,zone_starts(0:nprocs-1),&
  zone_ends(0:nprocs-1)
integer :: lo_nodes(3,num_zones),hi_nodes(3,num_zones)
real(kind=dp) :: min_corner(3,num_zones),max_corner(3,num_zones)
real(kind=dp) :: poco_loco(3),poco_hoco(3)
!  logical :: is_full(num_zones)

type(amrex_mfiter) :: mfi
type(amrex_box) :: brx

!write(*,*) 'entering zone ranges sub',self,lvl

counter = zone_starts(self)
call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
do while (mfi%next())

  brx = mfi%validbox()

  lo_nodes(1,counter) = brx%lo(1)
  lo_nodes(2,counter) = brx%lo(2)
  lo_nodes(3,counter) = brx%lo(3)
  hi_nodes(1,counter) = brx%hi(1)
  hi_nodes(2,counter) = brx%hi(2)
  hi_nodes(3,counter) = brx%hi(3)
!    write(*,*) 'lo nodes',self,lvl,counter,lo_nodes(1,counter),lo_nodes(2,counter),lo_nodes(3,counter)
!    write(*,*) 'hi nodes',self,lvl,counter,hi_nodes(1,counter),hi_nodes(2,counter),hi_nodes(3,counter)
    poco_loco = get_real_coords(amrex_geom(lvl),brx%lo)
    poco_hoco = get_real_coords(amrex_geom(lvl),brx%hi)

!  min_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%lo)
!  max_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%hi)

    min_corner(1,counter) = poco_loco(1)
    min_corner(2,counter) = poco_loco(2)
    min_corner(3,counter) = poco_loco(3)
    max_corner(1,counter) = poco_hoco(1)
    max_corner(2,counter) = poco_hoco(2)
    max_corner(3,counter) = poco_hoco(3)
!  write(*,*) self,lvl,min_corner(1:3,counter),max_corner(1:3,counter)
  counter = counter +1

end do

call amrex_mfiter_destroy(mfi)
end subroutine zone_ranges
