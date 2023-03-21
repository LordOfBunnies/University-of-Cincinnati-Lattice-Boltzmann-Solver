subroutine hdf5_output_test(fine)
!
!
!
!
!
!
use mpi
use precise
use amrex_amr_module
use amrex_base_module
use amrex_plotfile_module
use output_data
use nml_output
use grid_data
use linkwise
use amr_info_holder
use amr_processes, only: self,commune,count_big_zones,total_big_zones
implicit none
integer :: lvl,max_lvl,fine,lvl_diff,num,n_lvls,ier,loc_zones
character(len=32) :: filename,ts_format,write_num!,var_name
real(kind=dp),contiguous,pointer :: rho_cells(:,:,:,:)

type(amrex_multifab),allocatable :: mfcell_rho(:)
type(amrex_string) :: var_name(1)
type(amrex_mfiter) :: mfi
type(amrex_box) :: corx

num = 1
ts_format = '(I9.9)'
write(write_num,ts_format) timestep(amrex_max_level)
filename = 'rho_test'//write_num
!filename = "derp_test"
call amrex_string_build(var_name(1),"rho")

if (fine > save_shapes(2,num)) then
  max_lvl = save_shapes(2,num)
else
  max_lvl = fine
end if

do lvl = 0,max_lvl
loc_zones = count_big_zones(lvl)
  write(*,*) 'num boxs per proc',loc_zones,self,lvl
end do

lvl_diff = max_lvl - save_shapes(1,num) + 1

if (allocated(mfcell_rho)) deallocate(mfcell_rho)
allocate(mfcell_rho(0:amrex_max_level))
!allocate(mfcell_rho(1))
!write(*,*) 'snarf',lvl_diff
!
! This needs to be cell based for the next part
!
!do lvl = save_shapes(1,num),max_lvl
do lvl = 0,max_lvl
  write(*,*) 'squash',max_lvl,self,lvl,time(fine),timestep
  call amrex_multifab_build(mfcell_rho(lvl),mfrho(lvl)%ba,mfrho(lvl)%dm,1,0,&
    (/.false.,.false.,.false./))
end do
!
!do lvl = 0,max_lvl
!  write(*,*) 'bagus',self,lvl
!  call amrex_average_down_nd_to_cc(mfcell_rho(lvl),1,mfrho(lvl),1,1,0)
!  write(*,*) 'bau',self,lvl
!end do

call node_to_cell_averaging(mfcell_rho(0:max_lvl),max_lvl)

!do lvl = 0,max_lvl
!
!  call amrex_mfiter_build(mfi,mfcell_rho(lvl),tiling=.false.)
!  do while(mfi%next())
!    corx=mfi%validbox()
!
!    rho_cells => mfcell_rho(lvl)%dataptr(mfi)
!    write(*,*) 'REDRUM!!!',rho_cells(corx%lo(1),corx%lo(2),corx%lo(3),1),&
!      rho_cells(corx%hi(1),corx%hi(2),corx%hi(3),1),corx%lo,corx%hi,self,lvl
!
!  end do
!
!  call amrex_mfiter_destroy(mfi)
!
!end do
!call amrex_multifab_build(mfcell_rho(1),mfrho(1)%ba,mfrho(1)%dm,1,0)
!call amrex_average_down_nd_to_cc(mfcell_rho(1),1,mfrho(1),1,1,0)
!write(*,*) 'bau',self,lvl,amrex_ref_ratio
!do lvl = save_shapes(1,num),max_lvl
n_lvls = amrex_get_numlevels()
!call amrex_write_HDF5plotfile(filename,n_lvls,mfrho,var_name,amrex_geom,time(fine),&
!  timestep,amrex_ref_ratio)
!write(*,*) 'jingle',self
call amrex_write_HDF5plotfile(filename,n_lvls,mfcell_rho,var_name,amrex_geom,time(fine),&
  timestep,amrex_ref_ratio)
!call amrex_write_plotfile(filename,n_lvls,mfcell_rho,var_name,amrex_geom,time(fine),&
!  timestep,amrex_ref_ratio)
call mpi_barrier(commune,ier)

!call amrex_write_plotfile(filename,n_lvls,mfcell_rho,var_name,amrex_geom,time(fine),&
!  timestep,amrex_ref_ratio)

!  call amrex_write_HDF5plotfile(filename,lvl_diff,mfcell_rho(save_shapes(1,num):max_lvl),&
!    var_name,amrex_geom(save_shapes(1,num):max_lvl),time(fine),timestep(fine),amrex_ref_ratio)
!write(*,*) 'cintanya'
!end do


do lvl = max_lvl,0,-1
  call amrex_multifab_destroy(mfcell_rho(lvl))
end do
deallocate(mfcell_rho)
end subroutine
