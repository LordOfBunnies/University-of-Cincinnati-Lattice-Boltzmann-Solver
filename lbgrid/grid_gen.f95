SUBROUTINE grid_gen(characteristic_length)
!
! Generates the grid for the use by the LB solver
!
! Called by: uclbs_main
! Calls: grid_creation,fluid_or_solid,out_of_bounds,state_directions
!
use grid_data
!use nml_bounding_box
use nml_geometry
use geom_data
use nml_output
use freestream_values
use linkwise
use constants
use precise
use startup
use amr_info_holder
use ggm_stuff
use amr_processes, only: commune
use amrex_base_module
use amrex_amr_module
IMPLICIT NONE
REAL(KIND=dp) :: characteristic_length,x_dist,y_dist,z_dist,&
          x_center,y_center,z_center
real(kind=dp) :: loc_dx
INTEGER :: a,b,istat,cur_lvl_refinement,fine,max_level,ier,lvl,i

CHARACTER(len=80) :: filename_out
filename_out = 'uclbs.out'
!
OPEN(FILE=filename_out,UNIT = 11,STATUS='OLD',FORM = 'FORMATTED',&
  ACTION='WRITE',POSITION='APPEND',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Problems opening output file, closing.'
  WRITE(*,*) istat
  CALL error_out
END IF

!write(*,*) 'Starting grid values'

fine = amrex_get_finest_level()
!max_level = amrex_max_level()

!      WRITE(*,*) 'Oink!'
!      LOGICAL :: external_flow
x_geom_max = -10000.0D0
x_geom_min =  10000.0D0
y_geom_max = -10000.0D0
y_geom_min =  10000.0D0
z_geom_max = -10000.0D0
z_geom_min =  10000.0D0

!IF (cgns_start .OR. mixed_start) THEN
!  bounding_method = 1000
!END IF
!
! Loop through the geometry to find the absolute maximum and minimum for
! every coordinate.  This can help define things later.
!
DO a=1,total_tris !triangle number
  DO b = 1,3       ! Vertex number
    IF (tri_vertices(1,b,a) >= x_geom_max) THEN
      x_geom_max = tri_vertices(1,b,a)
    END IF
    IF (tri_vertices(2,b,a) >= y_geom_max) THEN
      y_geom_max = tri_vertices(2,b,a)
    END IF
    IF (tri_vertices(3,b,a) >= z_geom_max) THEN
      z_geom_max = tri_vertices(3,b,a)
    END IF
    IF (tri_vertices(1,b,a) <= x_geom_min) THEN
      x_geom_min = tri_vertices(1,b,a)
    END IF
    IF (tri_vertices(2,b,a) <= y_geom_min) THEN
      y_geom_min = tri_vertices(2,b,a)
    END IF
    IF (tri_vertices(3,b,a) <= z_geom_min) THEN
      z_geom_min = tri_vertices(3,b,a)
    END IF

  END DO
!        WRITE(*,*) 'Current triangle number, ',a
!        WRITE(*,*) x_geom_max,x_geom_min,y_geom_max,y_geom_min,z_geom_max,&
!          z_geom_min
!        write(*,*) tri_vertices(1:3,1:3,a),a
END DO
!       WRITE(*,*) 'squib'

write(11,*) ''
write(11,301) x_geom_min
write(11,302) y_geom_min
write(11,303) z_geom_min
write(11,304) x_geom_max
write(11,305) y_geom_max
write(11,306) z_geom_max
write(11,*) ''
!      WRITE(*,*) 'squib'
x_dist = x_geom_max - x_geom_min
y_dist = y_geom_max - y_geom_min
z_dist = z_geom_max - z_geom_min
x_center = (x_geom_max+x_geom_min)/2.0D0
y_center = (y_geom_max+y_geom_min)/2.0D0
z_center = (z_geom_max+z_geom_min)/2.0D0
longest_dist = MAX(x_dist,y_dist,z_dist)
!
!
!
write(11,201) x_dist
write(11,202) y_dist
write(11,203) z_dist
write(11,204) x_center
write(11,205) y_center
write(11,206) z_center
write(11,207) longest_dist

x_grid_max = amrex_probhi(1)
x_grid_min = amrex_problo(1)
y_grid_max = amrex_probhi(2)
y_grid_min = amrex_problo(2)
z_grid_max = amrex_probhi(3)
z_grid_min = amrex_problo(3)

write(11,101) x_grid_min
write(11,102) y_grid_min
write(11,103) z_grid_min
write(11,104) x_grid_max
write(11,105) y_grid_max
write(11,106) z_grid_max
write(11,*) ''
!
! Begin calls to create the grid and find if the nodes are inside or outside
!
loc_dx = amrex_geom(0)%dx(1)


call link_info(loc_dx)

call timing_calcs

call sanity_checks

call fluid_or_solid(fine)

call mpi_barrier(commune,ier)

call outlfow_setup(fine)
!
call mpi_barrier(commune,ier)
!
! State values
! 0-max level = level of a valid fluid node
! 200-299 = outflow nodes, special case
! 1000-50999 = ghost nodes overlapping a box on the same level
! 51000-99999 = valid coarse nodes overlapped by a fine box
! 100000-9999999 = valid coarse nodes overlapped by the fine ghost nodes, state contains
!              derivative information
! 10M-1.1B = ghost nodes not overlapped by another box on the same level, state contains
!              derivative information
!
!
!
! Basic breakdown of state-directions
! 0-max level = fluid, valid flow node
! 200-299 = outflow boundary node
! 300+ = see above
! -1 - -99 = basic boundary conditions
! -666 = out of bounds link, skipped for almost all purposes
! -100 - -199 = inlet nodes
! -200 - -299 = outlet nodes
! -1000 = freestream nodes
! -1.1B - -1001 = solid, state value used for bounceback values
!
if (shifted) then
  call shifted_state_directions(fine)
else
  call state_directions(fine)
end if
!
call mpi_barrier(commune,ier)


if (allocated(grid_info)) deallocate(grid_info)
allocate(grid_info(0:amrex_max_level,0:num_saves))

if (use_ggm) then
  if (allocated(ggm_save_info)) deallocate(ggm_save_info)
  allocate(ggm_save_info(0:amrex_max_level,num_ggm_vars))
end if
!
do lvl = 0,fine
  do i = 0,num_saves
    call update_grid_info(lvl,i)
  end do
end do
call mpi_barrier(commune,ier)
!
call nullify_nodes(fine,0)
!
call mpi_barrier(commune,ier)
!

  do lvl = 0,fine
    call simple_state_directions(lvl)
  end do
!
call mpi_barrier(commune,ier)
!
call ghost_node_id(fine,0)
!
call mpi_barrier(commune,ier)
!
do lvl = 0,fine
  call simple_state_directions(lvl)
end do
write(*,*) 'fluffy bunnies'
call mpi_barrier(commune,ier)
!
!call coarse_fine_node_id(fine)
!write(*,*) 'after coarse/fine'

!call initial_output
!do lvl = 0,fine
!  call simple_state_directions(lvl)
!end do
!
!call inter_bb_setup
!
do lvl = 0,fine
  write(11,401) amrex_geom(lvl)%domain%lo,lvl
  write(11,402) amrex_geom(lvl)%domain%hi,lvl
end do
write(11,*) 'Grid created and node states determined moving on ',&
  'to linking the mesh together.'

close(11)


 301  FORMAT ('Minimum geometric x-location = ',F15.7)
 302  FORMAT ('Minimum geometric y-location = ',F15.7)
 303  FORMAT ('Minimum geometric z-location = ',F15.7)
 304  FORMAT ('Maximum geometric x-location = ',F15.7)
 305  FORMAT ('Maximum geometric y-location = ',F15.7)
 306  FORMAT ('Maximum geometric z-location = ',F15.7)
! 301  FORMAT ('Maximum geometric x-location = ',F11.7)

 101  FORMAT ('Min grid point x coordinate = ',F9.5)
 102  FORMAT ('Min grid point y coordinate = ',F9.5)
 103  FORMAT ('Min grid point z coordinate = ',F9.5)
 104  FORMAT ('Max grid point x coordinate = ',F9.5)
 105  FORMAT ('Max grid point y coordinate = ',F9.5)
 106  FORMAT ('Max grid point z coordinate = ',F9.5)

 201  FORMAT ('Geometric maximum x-distance = ',F9.5)
 202  FORMAT ('Geometric maximum y-distance = ',F9.5)
 203  FORMAT ('Geometric maximum z-distance = ',F9.5)
 204  FORMAT ('Geometric maximum x-center = ',F9.5)
 205  FORMAT ('Geometric maximum y-center = ',F9.5)
 206  FORMAT ('Geometric maximum z-center = ',F9.5)
 207  FORMAT ('Longest geometric distance = ',F9.5)
! 208  FORMAT ('Voxel size = ',F9.7)

 401 format ('Nodal grid minimum bounds ',3I6,' on level ',I2)
 402 format ('Nodal grid maximum bounds ',3I6,' on level ',I2)

end subroutine










!      CALL grid_creation
!
! This checks to see if points in the rectangular grid are out of bounds
! of the bounding box
!
!call out_of_bounds
!
! fluid_or_solid finds where something should be declared a fluid or solid node
!

!write(*,*) 'Preparing to call fluid_or_solid'
!cur_lvl_refinement = 0
!do i = 0,box_level
!  call inital_multigrid_creation
!  call fluid_or_solid
!
!end do
!call grid_creation

!
! dx will be the basic cell size everywhere to start
!
!      dx = characteristic_length/REAL(cells_per_length)
!      WRITE(11,*) ''
!      WRITE(11,208) dx
!      WRITE(11,*) ''

! 202  FORMAT(,F9.7)
!
! Start with rectangulard bounds mean max and min for every
! coordinate
!
!IF (bounding_method == 1) THEN
!
!  x_grid_max = x_center + (bounds(1)/2.0D0)*x_dist
!  y_grid_max = y_center + (bounds(2)/2.0D0)*y_dist
!  z_grid_max = z_center + (bounds(3)/2.0D0)*z_dist
!  x_grid_min = x_center - (bounds(1)/2.0D0)*x_dist
!  y_grid_min = y_center - (bounds(2)/2.0D0)*y_dist
!  z_grid_min = z_center - (bounds(3)/2.0D0)*z_dist
!
!ELSE IF (bounding_method == 2) THEN
!
!  x_grid_min = bounds(1)
!  y_grid_min = bounds(2)
!  z_grid_min = bounds(3)
!  x_grid_max = bounds(4)
!  y_grid_max = bounds(5)
!  z_grid_max = bounds(6)
!
!ELSE IF (bounding_method == 3) THEN
!
!  x_grid_max = bounds(1) + (bounds(4)/2.0D0)*x_dist
!  y_grid_max = bounds(2) + (bounds(5)/2.0D0)*y_dist
!  z_grid_max = bounds(3) + (bounds(6)/2.0D0)*z_dist
!  x_grid_min = bounds(1) - (bounds(4)/2.0D0)*x_dist
!  y_grid_min = bounds(2) - (bounds(5)/2.0D0)*y_dist
!  z_grid_min = bounds(3) - (bounds(6)/2.0D0)*z_dist
!!
!! Spherical bounds
!!
!ELSE IF (bounding_method == 11) THEN
!
!  x_grid_max = bounds(1) + bounds(4)*longest_dist
!  y_grid_max = bounds(2) + bounds(4)*longest_dist
!  z_grid_max = bounds(3) + bounds(4)*longest_dist
!  x_grid_min = bounds(1) - bounds(4)*longest_dist
!  y_grid_min = bounds(2) - bounds(4)*longest_dist
!  z_grid_min = bounds(3) - bounds(4)*longest_dist
!
!ELSE IF (bounding_method == 12) THEN
!
!  x_grid_max = bounds(1) + bounds(4)
!  y_grid_max = bounds(2) + bounds(4)
!  z_grid_max = bounds(3) + bounds(4)
!  x_grid_min = bounds(1) - bounds(4)
!  y_grid_min = bounds(2) - bounds(4)
!  z_grid_min = bounds(3) - bounds(4)
!  !
!!
!! Cylindrical bounds
!!
!ELSE IF (bounding_method == 21) THEN
!
!  IF (bounds(5) > 0.0) THEN
!    x_grid_min = bounds(1)
!    x_grid_max = bounds(1) + bounds(5)
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!
!  ELSE IF (bounds(5) < 0.0) THEN
!    x_grid_max = bounds(1)
!    x_grid_min = bounds(1) + bounds(5)
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!  ELSE
!    WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!    CALL error_out
!  END IF
!
!ELSE IF (bounding_method == 22) THEN
!  IF (bounds(5) > 0.0) THEN
!    x_grid_min = bounds(1) - bounds(4)
!    x_grid_max = bounds(1) + bounds(4)
!    y_grid_min = bounds(2)
!    y_grid_max = bounds(2) + bounds(5)
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!
!  ELSE IF (bounds(5) < 0.0) THEN
!    x_grid_max = bounds(1) + bounds(4)
!    x_grid_min = bounds(1) - bounds(4)
!    y_grid_max = bounds(2)
!    y_grid_min = bounds(2) + bounds(5)
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!  ELSE
!    WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!    CALL error_out
!  END IF
!
!
!ELSE IF (bounding_method == 23) THEN
!  IF (bounds(5) > 0.0) THEN
!    x_grid_min = bounds(1) - bounds(4)
!    x_grid_max = bounds(1) + bounds(4)
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3)
!    z_grid_max = bounds(3) + bounds(5)
!
!  ELSE IF (bounds(5) < 0.0) THEN
!    x_grid_max = bounds(1) + bounds(4)
!    x_grid_min = bounds(1) - bounds(4)
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3) + bounds(5)
!    z_grid_max = bounds(3)
!  ELSE
!    WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!    CALL error_out
!  END IF
!
!ELSE IF (bounding_method == 24) THEN
!
!  IF (bounds(5) > 0.0) THEN
!    x_grid_min = bounds(1)
!    x_grid_max = bounds(1) + bounds(5)*x_dist
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!
!  ELSE IF (bounds(5) < 0.0) THEN
!    x_grid_max = bounds(1)
!    x_grid_min = bounds(1) + bounds(5)*x_dist
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!  ELSE
!    WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!    CALL error_out
!  END IF
!
!ELSE IF (bounding_method == 25) THEN
!
!  IF (bounds(5) > 0.0) THEN
!    x_grid_min = bounds(1) - bounds(4)
!    x_grid_max = bounds(1) + bounds(4)
!    y_grid_min = bounds(2)
!    y_grid_max = bounds(2) + bounds(5)*y_dist
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!
!  ELSE IF (bounds(5) < 0.0) THEN
!    x_grid_max = bounds(1) + bounds(4)
!    x_grid_min = bounds(1) - bounds(4)
!    y_grid_max = bounds(2)
!    y_grid_min = bounds(2) + bounds(5)*y_dist
!    z_grid_min = bounds(3) - bounds(4)
!    z_grid_max = bounds(3) + bounds(4)
!  ELSE
!    WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!    CALL error_out
!  END IF
!
!ELSE IF (bounding_method == 26) THEN
!
!  IF (bounds(5) > 0.0) THEN
!    x_grid_min = bounds(1) - bounds(4)
!    x_grid_max = bounds(1) + bounds(4)
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3)
!    z_grid_max = bounds(3) + bounds(5)*z_dist
!
!  ELSE IF (bounds(5) < 0.0) THEN
!    x_grid_max = bounds(1) + bounds(4)
!    x_grid_min = bounds(1) - bounds(4)
!    y_grid_max = bounds(2) + bounds(4)
!    y_grid_min = bounds(2) - bounds(4)
!    z_grid_min = bounds(3) + bounds(5)*z_dist
!    z_grid_max = bounds(3)
!  ELSE
!    WRITE(*,*) "Distance of 0 chosen for cylindrical bounds."
!    CALL error_out
!  END IF
!
!ELSE IF (bounding_method == 1000) THEN
!  x_grid_max = x_geom_max
!  x_grid_min = x_geom_min
!  y_grid_max = y_geom_max
!  y_grid_min = y_geom_min
!  z_grid_max = z_geom_max
!  z_grid_min = z_geom_min
!
!!
!! Planar bounds
!!
!ELSE IF (bounding_method == -1) THEN
!! For x-y plane for 2D
!  x_grid_min = bounds(1)
!  x_grid_max = bounds(2)
!  y_grid_min = bounds(3)
!  y_grid_max = bounds(4)
!  z_grid_min = bounds(5)
!  z_grid_max = bounds(5)
!!
!ELSE IF (bounding_method == -2) THEN
!! For x-z plane for 2D
!  x_grid_min = bounds(1)
!  x_grid_max = bounds(2)
!  z_grid_min = bounds(3)
!  z_grid_max = bounds(4)
!  y_grid_min = bounds(5)
!  y_grid_max = bounds(5)
!!
!ELSE IF (bounding_method == -3) THEN
!! For y-z plane for 2D
!  y_grid_min = bounds(1)
!  y_grid_max = bounds(2)
!  z_grid_min = bounds(3)
!  z_grid_max = bounds(4)
!  x_grid_min = bounds(5)
!  x_grid_max = bounds(5)
!
!ELSE
!  WRITE(*,*) 'Failure applying boundaries, exiting.'
!  CALL error_out
!END IF
!
! Write out the max and min grid parameters to the output file
!
