SUBROUTINE loop
!
! This is where the magic happens, everything else it setup.
!
! Called by: uclbs_main
! Calls: valid_levels,calculate_flow_variables,collision,streaming,
!   interlevel_interpolation,steps_counter,lb_output
! External calls: mpi_barrier
!
USE linkwise
!USE loop_data
USE precise
USE startup
USE constants
USE freestream_values
USE grid_data
USE mpi
use timing
USE amrex_base_module
USE amrex_amr_module
use amr_processes, only: amr_max_lvl,commune,self,amr_max_lvl
use amr_info_holder, only: time,time_prev,timestep
IMPLICIT NONE
INTEGER :: i,istat, fine,coarsest_lvl_run,lvl,ier
LOGICAL :: level_step(0:amr_max_lvl)

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
!
! This is the main loop for the entire program.  Here, the distro
! functions are calculated, and moved with time.
!
! Major changes to the program including the use of the pressure and shear tensors
! Pressure tensor has 3  components:
! 1 - xx
! 2 - yy
! 3 - zz
!
! Shear tensor has 3 components:
! 1 - xy
! 2 - xz
! 3 - yz
!
! For LMs:
! chi is 1 component
! ksi is 3 components:
!   1 - x
!   2 - y
!   3 - z
! pi is 6 components:
!   1 - xx
!   2 - yy
!   3 - zz
!   4 - xy
!   5 - xz
!   6 - yz
! lambda is 3 components:
!   1 - x
!   2 - y
!   3 - z
!
!
! Grid values:
!
! 0-max level = level of a valid fluid node
! 200-299 = outflow nodes, special case
! 1000-50999 = ghost nodes overlapping a box on the same level
! 51000-99999 = valid coarse nodes overlapped by a fine box
! 100000-9999999 = valid coarse nodes overlapped by the fine ghost nodes, state contains
!              derivative information
! 10M-1.1B = ghost nodes not overlapped by another box on the same level, state contains
!              derivative information
!
! -1 - -99 = boundary nodes
! -100 - -199 = inlets
! -200 - -299 = outlets
! -666 = null nodes
! -1000 = freestream nodes
! -1001 = solid nodes
! -1001 - -1E9 = solid links
!
!
call mpi_barrier(commune,ier)

fine = amrex_get_finest_level()
DO WHILE (time(0) <= total_time)
!        WRITE(*,*) t
!
!
  fine = amrex_get_finest_level()
  CALL valid_levels(level_step,fine,coarsest_lvl_run)
!
!
!write(*,*) 'looping through level ',lvl
  call mpi_barrier(commune,ier)
  do lvl = fine,coarsest_lvl_run,-1
    CALL streaming(level_step,fine,lvl)
  end do
!
!
!
  call mpi_barrier(commune,ier)
  do lvl = fine,coarsest_lvl_run,-1
    CALL calculate_flow_variables(level_step,fine,lvl)
  end do
!
! Interpolate the ghost nodes which lie outside the boxes on a given level
!
  call mpi_barrier(commune,ier)
  do lvl = fine,coarsest_lvl_run,-1
    call interlevel_interpolation(lvl,fine)
  end do
!
! Fine level to coarse level interpolation
!
  call mpi_barrier(commune,ier)
  do lvl = fine,coarsest_lvl_run,-1
    call fine_to_coarse(lvl,fine)
  end do
!
! Collision, find the new fiout and giout
!
  call mpi_barrier(commune,ier)
  do lvl = fine,coarsest_lvl_run,-1
    CALL collision(level_step,fine,lvl)
  end do

  call mpi_barrier(commune,ier)
  do lvl = fine,coarsest_lvl_run,-1
    CALL damp_it(level_step,fine,lvl)
  end do
!  CALL interlevel_interpolation(level_step,fine)
  call mpi_barrier(commune,ier)

  CALL lb_output(fine)
  CALL steps_counter(level_step,fine)
  call regrid_question_mark(fine)

END DO

CLOSE(11)

END SUBROUTINE
!
!      omega = dt/tau
!      omega = 0.933
!      WRITE(*,*) KIND(fi)
!      WRITE(*,*) KIND(u_vel)
!      omega = 1/tau
!      omega = .50018
!      IF (.NOT. restart_start) THEN
!        t_start = 1
!      END IF
!      WRITE(*,*) 'Beginning timestepping.'
! This is to solve for the new values after having been moved in the previous timestep.
!
!
!        DO j = 1,link_number
!          rho_sum = 0.0D0
!          u_sum = 0.0D0
!          v_sum = 0.0D0
!          w_sum = 0.0D0
!          IF (link(1,j) > 0) THEN
!            DO i = 1,dir
!              rho_sum = rho_sum + fi(i,j)
!              IF (cx(i) > 0) THEN
!                u_sum = u_sum + fi(i,j)
!              ELSE IF (cx(i) < 0) THEN
!                u_sum = u_sum - fi(i,j)
!              END IF
!
!              IF (cy(i) > 0) THEN
!                v_sum = v_sum + fi(i,j)
!              ELSE IF (cy(i) < 0) THEN
!                v_sum = v_sum - fi(i,j)
!              END IF
!
!              IF (cz(i) > 0) THEN
!                w_sum = w_sum + fi(i,j)
!              ELSE IF (cz(i) < 0) THEN
!                w_sum = w_sum - fi(i,j)
!              END IF
!
!
!            END DO
!
!            rho(j) = rho_sum
!            u_vel(j) = u_sum/rho(j)
!            v_vel(j) = v_sum/rho(j)
!            w_vel(j) = w_sum/rho(j)
!          END IF

!        END DO
!
! This is the calculation step where the distro functions are computed
!
!
!        DO j = 1,link_number
!! If it's a valid fluid node, then we can calculate the
!          IF (link(1,j) > 0) THEN
!            DO i = 1,dir
!
!
!              cu = 3.0D0*(cx(i)*u_vel(j)+cy(i)*v_vel(j)+cz(i)*w_vel(j))
!              vel_mag = (u_vel(j)**2.0D0+v_vel(j)**2.0D0+w_vel(j)**2.0D0)
!              feq = wi(i)*rho(j)*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)
!
!              fout(i,j) = fi(i,j)+omega*(feq-fi(i,j))
!
!
!            END DO
!          END IF
!        END DO

!        WRITE(*,*) 'cu,velmag,feq,fout525600 ',cu,vel_mag,feq,fout(1,525600)
!
!
! This is the movement step where the  values are moved from fout to
! fi.  This will also include bounce-back, inlet, outlet, and freestream
! boundary application.
!
!

!        DO j = 1,link_number
!!          IF (.NOT. boundary_node(j)) THEN
!            DO i = 1,dir
!!
!! Normal fluid nodes, move values from fout to fi via their directions
!!
!              IF (link(i,j) > 0) THEN
!                fi(i,link(i,j)) = fout(i,j)
!!
!! Handle bounceback for solid edge nodes for non-CGNS cases
!!
!              ELSE IF (link(i,j) == -999) THEN
!                upstream_node = link(opp(i),j)
!                IF (upstream_node < 1) THEN
!                  fout(i,j) = fout(opp(i),j)
!                  fi(opp(i),j) = fout(i,j)
!                ELSE
!                  fi(opp(i),j) = fout(i,j)
!                END IF
!!
!! Handle interpolative bounceback
!!
!              ELSE IF (link(i,j) < -1000) THEN
!!                CALL bounceback(fout(i,j),fout(opp(i),j))
!!                WRITE(*,*) 'Spang'
!!                fi(opp(i),j) = fout(i,j)
!                CALL interpolative_bounceback(fout(i,j),j,i,link(i,j),fout(opp(i),j))
!                fi(opp(i),j) = fout(i,j)
!!
!! For inviscid and symmetry boundary conditions
!!
!              ELSE IF (link(i,j) == -11) THEN
!
!                CALL inviscid_boundary(i,inviscid_dir,j)
!                fi(inviscid_dir,j) = fout(i,j)
!
!! For those nodes interacting with inlets
!!
!              ELSE IF (link(i,j) <= -100 .AND. link(i,j) >= -199) THEN
!                CALL inlet_bc(fout(i,j),i,j,&
!                       link(i,j))
!                fi(opp(i),j) = fout(i,j)
!!                WRITE(*,*) 'Goofy'
!!                IF (fi(opp(i),j) <= 0.0) THEN
!!                  WRITE(*,*) 'Negative fi, inlet ',fi(opp(i),j),i
!!                END IF
!!
!! For those nodes interacting with outlets
!!
!              ELSE IF (link(i,j) <= -200 .AND. link(i,j) >= -299) THEN
!                CALL outlet_bc(fout(i,j),i,j,&
!                       link(i,j))
!!                WRITE(*,*) 'Troopy'
!                fi(opp(i),j) = fout(i,j)
!!                IF (fi(opp(i),j) <= 0.0) THEN
!!                  WRITE(*,*) 'Negative fi, outlet ',fi(opp(i),j),i
!!                END IF
!!                fi(opp(i),link(i,j,1)) = fi(opp(i),j)
!!
!! For those nodes interacting with freestream conditions
!!
!              ELSE IF (link(i,j) == -1000) THEN
!                CALL freestream_edge_bc(fout(i,j),i)
!                fi(opp(i),j) = fout(i,j)
!!                IF (fi(opp(i),j) <= 0.0) THEN
!!                  WRITE(*,*) 'Negative fi, edge ',fi(opp(i),j),i
!!                END IF
!!                fi(opp(i),link(i,j,1)) = fi(opp(i),j)
!!
!! For those nodes interacting with null nodes (should not exist, but will
!! catch more of this in debugging
!!
!
!              ELSE IF (link(i,j) == -666) THEN
!                CALL null_bc( )
!              END IF
!
!            END DO
!!          END IF
!        END DO

!        DO j= 1,link_number
!          IF (boundary_node(j)) THEN
!            DO i = 1,dir
!              IF (link(i,j,2) == 0) THEN
!                fi(i,link(i,j,1)) = fout(i,j)
!!
!!      END IF
!              ELSE IF (link(i,j,2) == 1) THEN
!!                CALL bounceback(fout(i,j),fout(opp(i),j))
!
!                fi(opp(i),j) = fout(i,j)
!!                IF (link(opp(i),j,2) == 1) THEN
!!                  fi(i,j) = fi(opp(i),j)
!!                END IF
!!                fi(i,link(i,j,1)) = fout(i,j)
!!                fi(opp(i),j) = fout(opp(i),j)
!!                fi(i,link(i,j,1)) =
!!
!! For inviscid and symmetry boundary conditions
!!
!              ELSE IF (link(i,j,2) == 11) THEN
!
!                CALL inviscid_boundary(i,inviscid_dir,j)
!                fi(inviscid_dir,j) = fout(i,j)
!
!! For those nodes interacting with inlets
!!
!              ELSE IF (link(i,j,2) >= 100 .AND. link(i,j,2) <= 199) THEN
!                CALL inlet_bc(fout(i,j),i,j,&
!                       link(i,j,2))
!!                IF (j == welp) THEN
!!!                  WRITE(*,*) fout(i,j),i,opp(i)
!!               END IF
!                fi(opp(i),j) = fout(i,j)
!                IF (fi(opp(i),j) <= 0.0) THEN
!                  WRITE(*,*) 'Negative fi, inlet ',fi(opp(i),j),i
!                END IF
!!                fi(opp(i),link(i,j,1)) = fi(opp(i),j)
!!
!! For those nodes interacting with outlets
!!
!              ELSE IF (link(i,j,2) >= 200 .AND. link(i,j,2) <= 299) THEN
!                CALL outlet_bc(fout(i,j),i,j,&
!                       link(i,j,2))
!
!                fi(opp(i),j) = fout(i,j)
!                IF (fi(opp(i),j) <= 0.0) THEN
!                  WRITE(*,*) 'Negative fi, outlet ',fi(opp(i),j),i
!                END IF
!!                fi(opp(i),link(i,j,1)) = fi(opp(i),j)
!!
!! For those nodes interacting with freestream conditions
!!
!              ELSE IF (link(i,j,2) == 1000) THEN
!                CALL freestream_edge_bc(fout(i,j),i)
!                fi(opp(i),j) = fout(i,j)
!                IF (fi(opp(i),j) <= 0.0) THEN
!                  WRITE(*,*) 'Negative fi, edge ',fi(opp(i),j),i
!                END IF
!!                fi(opp(i),link(i,j,1)) = fi(opp(i),j)
!!
!! For those nodes interacting with null nodes (should not exist, but will
!! catch more of this in debugging
!!
!
!              ELSE IF (link(i,j,2) == -1) THEN
!                CALL null_bc( )
!              END IF
!            END DO
!          END IF
!
!        END DO

! 200  FORMAT ('Timestep ',I8)
! 1000 FORMAT (F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',&
!        F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ',&
!        F17.15, ' ',F17.15, ' ',F17.15, ' ',F17.15, ' ')
! 1001 FORMAT (19I11)


!  write(*,*) 'local time',time
!  write(*,*) 'previous time',time_prev
!  write(*,*) 'timestep',timestep

!  do lvl = coarsest_lvl_run,fine
