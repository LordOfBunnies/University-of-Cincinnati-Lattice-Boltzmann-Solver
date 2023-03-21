subroutine link_info(loc_dx)
!
! Determines a couple random things which are needed later.
! Placed here because of initialization issues.
!
! Called by: grid_gen
! Calls:
!
use precise
use constants
use grid_data
use amr_processes, only: amr_max_lvl
use linkwise
use freestream_values

implicit none
integer :: lvl,i
real(kind=dp) :: loc_dx
!write(*,*) 'Finding link distances and lattice Reynolds numbers.'
!
!
!
allocate(link_dist(dir,0:amr_max_lvl))
allocate(lattice_reynolds(0:amr_max_lvl))
do i = 1,dir
  do lvl = 0,amr_max_lvl
    link_dist(i,lvl) = sqrt(cmag(i))*loc_dx/2.0D0**lvl
!    write(*,*) 'direction ',i,' distance',link_dist(i,lvl),' at ',lvl
  end do

end do
!
! Find the freestream based lattice Reynolds number
!
do lvl = 0,amr_max_lvl
  lattice_reynolds(lvl) = rho_ref*u_ref*rdx(lvl)/mu_ref

  !write(*,101) lattice_reynolds(lvl),lvl
end do

 101  format ('Lattice Reynolds number ',F13.4,' at level',I3)

end subroutine
