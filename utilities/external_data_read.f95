subroutine external_data_read
!
! Retrieves useful information that generally goes directly into
! Amrex. This puts it up at the module level to be used later.
!
! Called by: external_code_start
! Calls:
! External calls: amrex_parmparse_build,amrex_parmparse_destroy
!
use constants
use precise
use amrex_base_module
use amrex_amr_module
use amr_info_holder, only: plot_interval,regrid_interval,&
  ref_rat,block_factor,max_grid,proper_shell,nodal_domain_limits,fixed_lvl
implicit none
integer :: lvl
type(amrex_parmparse) :: pp

plot_interval = 400
regrid_interval = 40
ref_rat = 2
block_factor = 8
max_grid = 64
!is_period = .false.
proper_shell = 2
fixed_lvl = 0

!write(*,*) 'Reading data from command line inputs file.'

call amrex_parmparse_build(pp,"amr")
call pp%query("regrid_int",regrid_interval)
call pp%query("plot_int",plot_interval)
call pp%query("ref_ratio",ref_rat)
call pp%query("blocking_factor",block_factor)
call pp%query("max_grid_size",max_grid)
call pp%query("n_proper",proper_shell)
call pp%query("use_fixed_upto_level",fixed_lvl)
call amrex_parmparse_destroy(pp)

do lvl = 0,amrex_max_level
  nodal_domain_limits(1,1,lvl) = 0
  nodal_domain_limits(2,1,lvl) = 0
  nodal_domain_limits(3,1,lvl) = 0

  nodal_domain_limits(1,2,lvl) = amrex_geom(lvl)%domain%hi(1)+1
  nodal_domain_limits(2,2,lvl) = amrex_geom(lvl)%domain%hi(2)+1
  nodal_domain_limits(3,2,lvl) = amrex_geom(lvl)%domain%hi(3)+1

  !write(*,*) 'nodal domain',nodal_domain_limits(1:3,1:2,lvl)
end do



!call amrex_parmparse_build(pp,"geometry")
!call pp%query("is_periodic",is_period)
!call amrex_parmparse_build(pp)

write(*,*) 'Command line input file read.'

end subroutine
