SUBROUTINE mic_setup_2d
!
! This sets up microphones for 2d since the rules are slightly different
!
! Calls:
! Called by: mic_setup
!
USE grid_data
USE nml_output
USE output_data
USE precise
USE constants
!USE nml_bounding_box
use amr_processes, only: commune,magic_box,get_real_coords
use amr_info_holder, only: zone_storage,mfstate,state
use amrex_base_module
use amrex_amr_module
use mpi
IMPLICIT NONE
REAL(KIND =dp) :: total_dist,local_dist(4),loc_lo(3)
INTEGER :: i,j,invalid_mic_nodes,x_start,y_start,z_start,fine,mic_cur_lvl
logical :: inside

type(amrex_mfiter) :: mfi
type(amrex_box) :: borks

!      ALLOCATE(mic_links(4,num_mics))
ALLOCATE(mic_weights(4,num_mics))
ALLOCATE(mic_valid(num_mics))
allocate(mic_nodes(2,4,num_mics))
allocate(mic_owner(num_mics))

  mic_valid = .TRUE.


mic_setup_loop: DO i = 1,num_mics


  IF (mic_locations(1,i) > x_grid_max .OR. mic_locations(1,i) &
    <x_grid_min .OR. mic_locations(2,i) > y_grid_max .OR. &
    mic_locations(2,i) < y_grid_min .OR. mic_locations(3,i) >= &
    z_grid_max .OR. mic_locations(3,i) <= z_grid_min) THEN

    mic_valid(i) = .FALSE.
    WRITE(11,100) i
    CYCLE mic_setup_loop
  END IF

  if (fine >= mic_lvl(i)) then
    mic_cur_lvl = mic_lvl(i)
  else
    mic_cur_lvl = fine
  end if



  x_start = FLOOR((mic_locations(1,i)-x_grid_min)/rdx(mic_cur_lvl))
  y_start = FLOOR((mic_locations(2,i)-y_grid_min)/rdx(mic_cur_lvl))
  z_start = 0!FLOOR((mic_locations(3,i)-z_grid_min)/rdx(mic_cur_lvl))

  loc_lo = get_real_coords(amrex_geom(mic_cur_lvl),(/x_start,y_start,z_start/))
!        WRITE(*,*) x_start,y_start,z_start
!
! This find the local distance and local link numbers
!
  call amrex_mfiter_build(mfi,mfstate(mic_cur_lvl),tiling=.false.)
  state => mfstate(mic_cur_lvl)%dataptr(mfi)
  inside = .false.
  do while (mfi%next())
    borks = mfi%validbox()
    !inside = magic_box(borks%lo,borks%hi,mic_locations(1:3,i))
    inside = magic_box(borks%lo,borks%hi,(/x_start,y_start,z_start/))
    if (inside) then
      invalid_mic_nodes = 0
      total_dist = 0.0D0
      DO j = 1,4
!
! Go through all 8 possible nodes to be able to pull from
!
      SELECT CASE (j)
      CASE(1)
        IF (state(x_start,y_start,z_start,1) >= 0) THEN
          mic_nodes(1,j,i) = x_start
          mic_nodes(2,j,i) = y_start
!          mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start,y_start,z_start)
  !
! Find the local distance for weighting purposes
!
        local_dist(1) = SQRT((mic_locations(1,i)-loc_lo(1))**2+&
        (mic_locations(2,i)-loc_lo(2))**2+&
        (mic_locations(3,i)-loc_lo(3))**2)

!
        total_dist = total_dist + local_dist(1)
      ELSE
!        mic_links(j,i) = 0
        mic_weights(j,i) = 0.0D0
        invalid_mic_nodes = invalid_mic_nodes + 1
      END IF
!
    CASE(2)
      IF (state(x_start+1,y_start,z_start,1) >= 0) THEN
        mic_nodes(1,j,i) = x_start+1
        mic_nodes(2,j,i) = y_start
!        mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start)
!
! Find the local distance for weighting purposes
!
        local_dist(2) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
        (mic_locations(2,i)-loc_lo(2))**2)

        total_dist = total_dist+local_dist(2)
      ELSE
!        mic_links(j,i) = 0
        mic_weights(j,i) = 0.0D0
        invalid_mic_nodes = invalid_mic_nodes + 1
      END IF
    CASE(3)
      IF (state(x_start,y_start+1,z_start,1) >= 0) THEN
        mic_nodes(1,j,i) = x_start
        mic_nodes(2,j,i) = y_start+1
!        mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start)
!
! Find the local distance for weighting purposes
!
        local_dist(3) = SQRT((mic_locations(1,i)-loc_lo(1))**2+&
        (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl)+rdx(mic_cur_lvl))**2)

        total_dist = total_dist+local_dist(3)
      ELSE
!        mic_links(j,i) = 0
        mic_weights(j,i) = 0.0D0
        invalid_mic_nodes = invalid_mic_nodes + 1
      END IF

    CASE(4)
      IF (state(x_start+1,y_start+1,z_start,1) >= 0) THEN
        mic_nodes(1,j,i) = x_start+1
        mic_nodes(2,j,i) = y_start+1
        mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start,y_start,z_start+1)
!
! Find the local distance for weighting purposes
!
        local_dist(4) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
        (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2)

        total_dist = total_dist+local_dist(4)
      ELSE
!        mic_links(j,i) = 0
        mic_weights(j,i) = 0.0D0
        invalid_mic_nodes = invalid_mic_nodes + 1
      END IF
!    CASE(5)
!      IF (state(x_start+1,y_start+1,z_start,1) >= 0) THEN
!        mic_nodes(1,j,i) = x_start+1
!        mic_nodes(2,j,i) = y_start+1
!        mic_nodes(3,j,i) = z_start
!!        mic_links(j,i) = state(dir+1,x_start+1,y_start+1,z_start)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(5) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
!        (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
!        (mic_locations(3,i)-loc_lo(3))**2)
!
!        total_dist = total_dist+local_dist(5)
!      ELSE
!!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!    CASE(6)
!      IF (state(x_start+1,y_start,z_start+1,1) >= 0) THEN
!        mic_nodes(1,j,i) = x_start+1
!        mic_nodes(2,j,i) = y_start
!        mic_nodes(3,j,i) = z_start+1
!!        mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(6) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
!        (mic_locations(2,i)-loc_lo(2))**2+&
!        (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)
!
!        total_dist = total_dist+local_dist(6)
!      ELSE
!!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!    CASE(7)
!      IF (state(x_start,y_start+1,z_start+1,1) >= 0) THEN
!        mic_nodes(1,j,i) = x_start
!        mic_nodes(2,j,i) = y_start+1
!        mic_nodes(3,j,i) = z_start+1
!!        mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(7) = SQRT((mic_locations(1,i)-loc_lo(1))**2+&
!        (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
!        (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)
!
!        total_dist = total_dist+local_dist(7)
!      ELSE
!!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!
!      END IF
!    CASE(8)
!      IF (state(x_start+1,y_start+1,z_start+1,1) >= 0) THEN
!        mic_nodes(1,j,i) = x_start+1
!        mic_nodes(2,j,i) = y_start+1
!        mic_nodes(3,j,i) = z_start+1
!!        mic_links(j,i) = state(dir+1,x_start+1,y_start+1,&
!!          z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(8) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
!        (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
!        (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)
!
!        total_dist = total_dist+local_dist(8)
!      ELSE
!!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!
!      END IF
    END SELECT
    END DO

    IF (invalid_mic_nodes == 4) THEN
      mic_valid(i) = .FALSE.
      WRITE(11,101) i
      WRITE(11,100) i
      CYCLE mic_setup_loop
    END IF

    WRITE(11,*) 'Finished determining microphone distances.'

    DO j = 1,4
      if (mic_nodes(1,j,i)> 0) then
!    IF (mic_links(j,i) > 0) THEN
        mic_weights(j,i) = local_dist(j)/total_dist
!            WRITE(*,*) mic_weights(j,i),local_dist(j),total_dist
      ELSE
        mic_weights(j,i) = 0.0D0
      END IF

    END DO

    WRITE(11,*) 'Finished determining microphone weighting.'
    else
      CYCLE
    end if
  end do


end do mic_setup_loop

!      IF (bounding_method == -1) THEN
!        IF (mic_locations(1,i) > x_grid_max .OR. mic_locations(1,i) &
!          <x_grid_min .OR. mic_locations(2,i) > y_grid_max .OR. &
!          mic_locations(2,i) < y_grid_min) THEN
!
!          mic_valid(i) = .FALSE.
!          WRITE(11,100) i
!          CYCLE mic_setup_loop
!        END IF
!
!!
!! Find the iteration parameters
!!
!        x_start = FLOOR((mic_locations(1,i)-x_grid_min)/dx)
!        y_start = FLOOR((mic_locations(2,i)-y_grid_min)/dx)
!        z_start = 1
!
!        invalid_mic_nodes = 0
!        total_dist = 0.0D0
!
!        DO j = 1,4
!!
!! Go through all 8 possible nodes to be able to pull from
!!
!          SELECT CASE (j)
!          CASE(1)
!            IF (state(1,x_start,y_start,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!              y_start,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start,&
!              y_start,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(2)
!            IF (state(1,x_start+1,y_start,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!              y_start,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start+1,&
!              y_start,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(3)
!            IF (state(1,x_start+1,y_start+1,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start+1,y_start+1,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!              y_start+1,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start+1,&
!              y_start+1,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(4)
!            IF (state(1,x_start,y_start+1,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!              y_start+1,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start,&
!              y_start+1,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          END SELECT
!        END DO
!        IF (invalid_mic_nodes == 8) THEN
!          mic_valid(i) = .FALSE.
!          WRITE(11,101) i
!          WRITE(11,100) i
!          CYCLE mic_setup_loop
!        END IF
!
!        WRITE(11,*) 'Finished determining microphone distances.'
!
!        DO j = 1,8
!          IF (mic_links(j,i) > 0) THEN
!            mic_weights(j,i) = local_dist(j)/total_dist
!!            WRITE(*,*) mic_weights(j,i),local_dist(j),total_dist
!          ELSE
!            mic_weights(j,i) = 0.0D0
!          END IF
!
!        END DO
!
!        WRITE(11,*) 'Finished determining microphone weighting.'
!
!      ELSE IF (bounding_method == -2) THEN
!        IF (mic_locations(1,i) > x_grid_max .OR. mic_locations(1,i) &
!          <x_grid_min .OR. mic_locations(3,i) > &
!          z_grid_max .OR. mic_locations(3,i) < z_grid_min) THEN
!
!          mic_valid(i) = .FALSE.
!          WRITE(11,100) i
!          CYCLE mic_setup_loop
!        END IF
!!
!! Find the iteration parameters
!!
!        x_start = FLOOR((mic_locations(1,i)-x_grid_min)/dx)
!        y_start = 1
!        z_start = FLOOR((mic_locations(3,i)-z_grid_min)/dx)
!
!        invalid_mic_nodes = 0
!        total_dist = 0.0D0
!
!        DO j = 1,4
!!
!! Go through all 8 possible nodes to be able to pull from
!!
!          SELECT CASE (j)
!          CASE(1)
!            IF (state(1,x_start,y_start,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!              y_start,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!              y_start,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(2)
!            IF (state(1,x_start+1,y_start,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!              y_start,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start+1,&
!              y_start,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(3)
!            IF (state(1,x_start+1,y_start,z_start+1) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!              y_start,z_start+1))**2.0D0+(mic_locations(3,i)-grid(3,x_start+1,&
!              y_start,z_start+1))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(4)
!            IF (state(1,x_start,y_start,z_start+1) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!              y_start,z_start+1))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!              y_start,z_start+1))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          END SELECT
!        END DO
!        IF (invalid_mic_nodes == 8) THEN
!          mic_valid(i) = .FALSE.
!          WRITE(11,101) i
!          WRITE(11,100) i
!          CYCLE mic_setup_loop
!        END IF
!
!        WRITE(11,*) 'Finished determining microphone distances.'
!
!        DO j = 1,8
!          IF (mic_links(j,i) > 0) THEN
!            mic_weights(j,i) = local_dist(j)/total_dist
!!            WRITE(*,*) mic_weights(j,i),local_dist(j),total_dist
!          ELSE
!            mic_weights(j,i) = 0.0D0
!          END IF
!
!        END DO
!
!        WRITE(11,*) 'Finished determining microphone weighting.'
!
!
!      ELSE IF (bounding_method == -3) THEN
!        IF (mic_locations(2,i) > y_grid_max .OR. &
!          mic_locations(2,i) < y_grid_min .OR. mic_locations(3,i) > &
!          z_grid_max .OR. mic_locations(3,i) < z_grid_min) THEN
!
!          mic_valid(i) = .FALSE.
!          WRITE(11,100) i
!          CYCLE mic_setup_loop
!        END IF
!
!!
!! Find the iteration parameters
!!
!        x_start = 1
!        y_start = FLOOR((mic_locations(2,i)-y_grid_min)/dx)
!        z_start = FLOOR((mic_locations(3,i)-z_grid_min)/dx)
!
!        invalid_mic_nodes = 0
!        total_dist = 0.0D0
!
!        DO j = 1,4
!!
!! Go through all 8 possible nodes to be able to pull from
!!
!          SELECT CASE (j)
!          CASE(1)
!            IF (state(1,x_start,y_start,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(2,i)-grid(2,x_start,&
!              y_start,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!              y_start,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(2)
!            IF (state(1,x_start,y_start+1,z_start) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(2,i)-grid(2,x_start,&
!              y_start+1,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!              y_start+1,z_start))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(3)
!            IF (state(1,x_start,y_start+1,z_start+1) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(2,i)-grid(2,x_start,&
!              y_start+1,z_start+1))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!              y_start+1,z_start+1))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          CASE(4)
!            IF (state(1,x_start,y_start,z_start+1) == 0) THEN
!              mic_links(j,i) = state(dir+1,x_start,y_start,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!              local_dist(1) = SQRT((mic_locations(2,i)-grid(2,x_start,&
!              y_start,z_start+1))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!              y_start,z_start+1))**2.0D0)
!
!!
!              total_dist = total_dist + local_dist(1)
!            ELSE
!              mic_links(j,i) = 0
!              mic_weights(j,i) = 0.0D0
!              invalid_mic_nodes = invalid_mic_nodes + 1
!            END IF
!          END SELECT
!        END DO
!        IF (invalid_mic_nodes == 4) THEN
!          mic_valid(i) = .FALSE.
!          WRITE(11,101) i
!          WRITE(11,100) i
!          CYCLE mic_setup_loop
!        END IF
!
!        WRITE(11,*) 'Finished determining microphone distances.'
!
!        DO j = 1,8
!          IF (mic_links(j,i) > 0) THEN
!            mic_weights(j,i) = local_dist(j)/total_dist
!!            WRITE(*,*) mic_weights(j,i),local_dist(j),total_dist
!          ELSE
!            mic_weights(j,i) = 0.0D0
!          END IF
!
!        END DO
!
!        WRITE(11,*) 'Finished determining microphone weighting.'
!
!      END IF
!
!      END DO mic_setup_loop

 100  FORMAT ('Microphone ',I5,' is not in a valid location')
 101  FORMAT ('Microphone ',I5,' is surrounded by invalid nodes',&
        ' and is therefore invalid.')
      END SUBROUTINE mic_setup_2d
