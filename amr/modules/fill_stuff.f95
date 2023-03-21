module fill_stuff
!
!
!
!
use precise
use amrex_base_module
use amrex_amr_module
use iso_c_binding
use grid_data
implicit none
save
private
public :: edge_bounds,allocate_fi_bcs,edge_bounds_distro
public :: lo_bc,hi_bc,lo_bc_fi,hi_bc_fi

integer :: lo_bc(amrex_spacedim,1) = amrex_bc_int_dir
integer :: hi_bc(amrex_spacedim,1) = amrex_bc_int_dir
integer,allocatable :: lo_bc_fi(:,:),hi_bc_fi(:,:)
!integer :: ext_bc(amrex_spacedim,1) = amrex_bc_ext_dir
!integer ::

contains

subroutine edge_bounds(loc_mf,scomp,ncomp,time,loc_geom) bind(c)
use amrex_filcc_module,only : amrex_filcc
!use

real(kind=dp),contiguous,pointer :: loc_stuff(:,:,:,:)
integer :: loc_lo(4),loc_hi(4)

type(amrex_geometry) :: sub_geom
type(amrex_multifab) :: mother_fluffer
type(amrex_mfiter) :: mfi
type(c_ptr),value :: loc_mf,loc_geom
integer(c_int),value :: scomp,ncomp
real(amrex_real),value :: time

sub_geom = loc_geom
mother_fluffer = loc_mf
call amrex_mfiter_build(mfi,mother_fluffer,tiling=.false.)
do while(mfi%next())
  loc_stuff => mother_fluffer%dataptr(mfi)
  if (.not. sub_geom%domain%contains(loc_stuff)) then

    loc_lo = lbound(loc_stuff)
    loc_hi = ubound(loc_stuff)
    call amrex_filcc(loc_stuff,loc_lo,loc_hi,sub_geom%domain%lo,&
      sub_geom%domain%hi,sub_geom%dx,sub_geom%get_physical_location(loc_lo),&
      lo_bc,hi_bc)

  end if

end do

end subroutine
!
!
!
subroutine edge_bounds_distro(loc_mf,scomp,ncomp,time,loc_geom) bind(c)
use amrex_filcc_module,only : amrex_filcc
!use

real(kind=dp),contiguous,pointer :: loc_stuff(:,:,:,:)
integer :: loc_lo(4),loc_hi(4)

type(amrex_geometry) :: sub_geom
type(amrex_multifab) :: mother_fluffer
type(amrex_mfiter) :: mfi
type(c_ptr),value :: loc_mf,loc_geom
integer(c_int),value :: scomp,ncomp
real(amrex_real),value :: time

sub_geom = loc_geom
mother_fluffer = loc_mf
call amrex_mfiter_build(mfi,mother_fluffer,tiling=.false.)
do while(mfi%next())
  loc_stuff => mother_fluffer%dataptr(mfi)
  if (.not. sub_geom%domain%contains(loc_stuff)) then

    loc_lo = lbound(loc_stuff)
    loc_hi = ubound(loc_stuff)
    call amrex_filcc(loc_stuff,loc_lo,loc_hi,sub_geom%domain%lo,&
      sub_geom%domain%hi,sub_geom%dx,sub_geom%get_physical_location(loc_lo),&
      lo_bc_fi,hi_bc_fi)

  end if

end do

end subroutine
!
!
!
subroutine allocate_fi_bcs

if (allocated(lo_bc_fi)) deallocate(lo_bc_fi)
if (allocated(hi_bc_fi)) deallocate(hi_bc_fi)

allocate(lo_bc_fi(amrex_spacedim,dir))
allocate(hi_bc_fi(amrex_spacedim,dir))

lo_bc_fi = amrex_bc_int_dir
hi_bc_fi = amrex_bc_int_dir

end subroutine

end module
