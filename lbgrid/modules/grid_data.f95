MODULE grid_data
!
! Stores the data for the grid 
!
! Now will also store data for wall nodes. Right now they're being
! implemented to do interpolative bounceback, but eventually, they're
! going to be able to do surface microphones.
!
!
USE precise
!use amrex_amr_module
!use amrex_base_module

IMPLICIT NONE
SAVE
real(KIND=dp) :: x_geom_max,y_geom_max,z_geom_max,x_geom_min, &
          y_geom_min,z_geom_min
real(KIND=dp) :: x_grid_max,y_grid_max,z_grid_max,x_grid_min, &
          y_grid_min,z_grid_min
real(KIND=dp) :: longest_dist,grid_center(3)
integer :: num_x,num_y,num_z
integer :: dir
integer :: num_solid_dirs,diag_dx,num_card_solid_dirs,num_diag_solid_dirs

real(kind=dp),allocatable :: rdx(:)
!
! Allocatable arrays for big things where we don't know how much there is.
!
!INTEGER,ALLOCATABLE :: associated_tri(:,:,:,:)!,surface_lattice(:,:)

!REAL(KIND=dp),ALLOCATABLE :: wall_nodes(:,:),wall_dist(:)

!type(amrex_multifab) :: mf_asstri(:)
CONTAINS

function ray_parallel_adjuster(pt_1,pt_2,pt_3,normal_to_tri_x,&
    normal_to_tri_y,normal_to_tri_z,loc_dx) result (ray_end2)
!
! Changing from a subroutine to a function
!
real(kind=dp) :: pt_1,pt_2,pt_3
REAL(KIND = dp) :: ray_end2(3),normal_to_tri_x,&
  normal_to_tri_y,normal_to_tri_z,loc_dx

ray_end2(1) = pt_1 + 10.0D0*loc_dx*normal_to_tri_x
ray_end2(2) = pt_2 + 10.0D0*loc_dx*normal_to_tri_y
ray_end2(3) = pt_3 + 10.0D0*loc_dx*normal_to_tri_z



end function

function different_parallel_adjuster(normal_x,normal_y,normal_z,lvl) result (ray_end1)
!
use linkwise
!
real(kind=dp) :: normal_x,normal_y,normal_z,ray_end1(3)
integer :: lvl

ray_end1(1) = ray_end1(1) + rdx(lvl)*normal_x
ray_end1(2) = ray_end1(2) + rdx(lvl)*normal_y
ray_end1(3) = ray_end1(3) + rdx(lvl)*normal_z
end function

subroutine fluff_it_all(ray_end1,ray_end2,normal_x,normal_y,normal_z,lvl)
!
!
!
real(kind=dp) :: ray_end1(3),ray_end2(3)
real(kind=dp) :: normal_x,normal_y,normal_z
integer :: lvl


ray_end1(1) = ray_end1(1) + rdx(lvl)*normal_x
ray_end1(2) = ray_end1(2) + rdx(lvl)*normal_y
ray_end1(3) = ray_end1(3) + rdx(lvl)*normal_z

ray_end2(1) = ray_end2(1) - 100*rdx(lvl)*normal_z
ray_end2(2) = ray_end2(2) - 100*rdx(lvl)*normal_x
ray_end2(3) = ray_end2(3) - 100*rdx(lvl)*normal_y

end subroutine

function cgns_parallel_adjuster(lvl,lc) result (ray_end2)
!

real(kind=dp) :: ray_end2(3)
integer :: lvl,lc

ray_end2(1) = ray_end2(1) + lc*rdx(lvl)/10.0D0
ray_end2(2) = ray_end2(2) + lc*rdx(lvl)/10.0D0
ray_end2(3) = ray_end2(3) + lc*rdx(lvl)/10.0D0

end function

function cgns_node_push(lvl,loc_dir) result(ray_end1)
!
use linkwise
!
real(kind=dp) :: ray_end1(3)
integer :: lvl,loc_dir

ray_end1(1) = ray_end1(1) - rdx(lvl)*rcx(loc_dir)/10.0D0
ray_end1(2) = ray_end1(2) - rdx(lvl)*rcy(loc_dir)/10.0D0
ray_end1(3) = ray_end1(3) - rdx(lvl)*rcz(loc_dir)/10.0D0

end function

END MODULE
