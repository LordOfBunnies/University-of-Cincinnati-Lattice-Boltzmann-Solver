SUBROUTINE mic_setup
!
! This is to set up the mics and which nodes to pull from
!
! Called by: grid_gen
! Calls: mic_initialize
!
USE grid_data
USE nml_output
USE output_data
USE precise
USE constants
!USE nml_bounding_box
use amr_processes, only: nprocs,commune,magic_box,get_real_coords,self
use amr_info_holder, only: mfstate,state!,zone_storage
use amrex_base_module
use amrex_amr_module
use mpi
IMPLICIT NONE
INTEGER :: i,j,ier,invalid_mic_nodes,req_counter,&
  x_start,y_start,z_start,fine,mic_cur_lvl!,taste(mpi_status_size)
integer,allocatable :: req_array(:),stat_array(:,:)
REAL(KIND=dp) :: total_dist,local_dist(8),borks_lo(3),borks_hi(3),loc_lo(3)
logical :: inside

type(amrex_box) :: borks
type(amrex_mfiter) :: mfi


!
!
!
fine = amrex_get_finest_level()
!
! We will skip any non-fluid nodes and make sure to take the
!
! Mic node order
! (x,y,z)
! (x+1,y,z)
! (x,y+1,z)
! (x,y,z+1)
! (x+1,y+1,z)
! (x+1,y,z+1)
! (x,y+1,z+1)
! (x+1,y+1,z+1)
!
!      ALLOCATE(mic_links(8,num_mics))
!      ALLOCATE(mic_weights(8,num_mics))
!      ALLOCATE(mic_valid(num_mics))
!      mic_valid = .TRUE.

WRITE(11,*) ''
WRITE(11,*) 'Beginning mic setup'
!WRITE(*,*) 'Beginning mic setup'
WRITE(11,*) 'Start by determining distances from specified locations',&
  ' to 8 nearest nodes.'
IF (dimensions == 2) THEN
  CALL mic_setup_2d
  RETURN
ELSE
!  ALLOCATE(mic_links(8,num_mics))
  ALLOCATE(mic_weights(8,num_mics))
  ALLOCATE(mic_valid(num_mics))
  allocate(mic_nodes(3,8,num_mics))
  allocate(mic_owner(num_mics))

!  send_count = num_mics*(nprocs-1)
  allocate(stat_array(MPI_STATUS_SIZE,nprocs-1))
  allocate(req_array(nprocs-1))
  mic_owner = -1

  mic_valid = .TRUE.
END IF
! for the request array for the eventual wait_all
req_counter = 0
!
mic_setup_loop: DO i = 1,num_mics

!  write(*,*) mic_locations
  IF (mic_locations(1,i) > x_grid_max .OR. mic_locations(1,i) &
    <x_grid_min .OR. mic_locations(2,i) > y_grid_max .OR. &
    mic_locations(2,i) < y_grid_min .OR. mic_locations(3,i) > &
    z_grid_max .OR. mic_locations(3,i) < z_grid_min) THEN

    mic_valid(i) = .FALSE.
    WRITE(11,100) i
    CYCLE mic_setup_loop
  END IF
!
  if (fine >= mic_lvl(i)) then
    mic_cur_lvl = mic_lvl(i)
  else
    mic_cur_lvl = fine
  end if

  x_start = FLOOR((mic_locations(1,i)-amrex_problo(1))/rdx(mic_cur_lvl))
  y_start = FLOOR((mic_locations(2,i)-amrex_problo(2))/rdx(mic_cur_lvl))
  z_start = FLOOR((mic_locations(3,i)-amrex_problo(3))/rdx(mic_cur_lvl))

  loc_lo = get_real_coords(amrex_geom(mic_cur_lvl),(/x_start,y_start,z_start/))
!  WRITE(*,*) 'mic_locations',mic_locations(1,i),mic_locations(2,i)&
!    ,mic_locations(3,i),x_start,y_start,z_start,amrex_problo(1),amrex_problo(2),&
!    amrex_problo(3),mic_cur_lvl,rdx(mic_cur_lvl)
!
! This find the local distance and local link numbers
!
!  write(*,*) 'woof?'
  call amrex_mfiter_build(mfi,mfstate(mic_cur_lvl),tiling=.false.)
  state => mfstate(mic_cur_lvl)%dataptr(mfi)
  inside = .false.
  do while (mfi%next())
    borks = mfi%validbox()
    inside = magic_box(borks%lo,borks%hi,(/x_start,y_start,z_start/))
    if (inside) then
      invalid_mic_nodes = 0
      total_dist = 0.0D0
      mic_owner(i) = self
!      do j = 0,nprocs-1
!        if (j == self) then
!          cycle
!        end if
!        req_counter = req_counter+1
!        tag = 1000*i+j
!        req_array(req_counter) = tag
!        write(*,*) 'boo'
!        call mpi_isend(mic_owner(i),1,MPI_INT,j,tag,commune,req_array(req_counter),ier)
!      end do
    !end if
!      write(*,*) 'node inside box',borks%lo,borks%hi,x_start,y_start,z_start
      DO j = 1,8
!
! Go through all 8 possible nodes to be able to pull from
!
        SELECT CASE (j)
        CASE(1)
!          write(*,*) 'case 1'
          IF (state(x_start,y_start,z_start,1) >= 0) THEN
!            write(*,*) 'inside state info for case 1'
            mic_nodes(1,j,i) = x_start
            mic_nodes(2,j,i) = y_start
            mic_nodes(3,j,i) = z_start
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
            mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start)
!
! Find the local distance for weighting purposes
!
            local_dist(2) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
            (mic_locations(2,i)-loc_lo(2))**2+&
            (mic_locations(3,i)-loc_lo(3))**2)

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
            mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start)
!
! Find the local distance for weighting purposes
!
            local_dist(3) = SQRT((mic_locations(1,i)-loc_lo(1))**2+&
            (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
            (mic_locations(3,i)-loc_lo(3))**2)

            total_dist = total_dist+local_dist(3)
          ELSE
!        mic_links(j,i) = 0
            mic_weights(j,i) = 0.0D0
            invalid_mic_nodes = invalid_mic_nodes + 1
          END IF

        CASE(4)
          IF (state(x_start,y_start,z_start+1,1) >= 0) THEN
            mic_nodes(1,j,i) = x_start
            mic_nodes(2,j,i) = y_start
            mic_nodes(3,j,i) = z_start+1
!        mic_links(j,i) = state(dir+1,x_start,y_start,z_start+1)
!
! Find the local distance for weighting purposes
!
            local_dist(4) = SQRT((mic_locations(1,i)-loc_lo(1))**2+&
            (mic_locations(2,i)-loc_lo(2))**2+&
            (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)

            total_dist = total_dist+local_dist(4)
          ELSE
!        mic_links(j,i) = 0
            mic_weights(j,i) = 0.0D0
            invalid_mic_nodes = invalid_mic_nodes + 1
          END IF
        CASE(5)
          IF (state(x_start+1,y_start+1,z_start,1) >= 0) THEN
            mic_nodes(1,j,i) = x_start+1
            mic_nodes(2,j,i) = y_start+1
            mic_nodes(3,j,i) = z_start
!        mic_links(j,i) = state(dir+1,x_start+1,y_start+1,z_start)
!
! Find the local distance for weighting purposes
!
            local_dist(5) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
            (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
            (mic_locations(3,i)-loc_lo(3))**2)

            total_dist = total_dist+local_dist(5)
          ELSE
!        mic_links(j,i) = 0
            mic_weights(j,i) = 0.0D0
            invalid_mic_nodes = invalid_mic_nodes + 1
          END IF
        CASE(6)
          IF (state(x_start+1,y_start,z_start+1,1) >= 0) THEN
            mic_nodes(1,j,i) = x_start+1
            mic_nodes(2,j,i) = y_start
            mic_nodes(3,j,i) = z_start+1
!        mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start+1)
!
! Find the local distance for weighting purposes
!
            local_dist(6) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
            (mic_locations(2,i)-loc_lo(2))**2+&
            (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)

            total_dist = total_dist+local_dist(6)
          ELSE
!        mic_links(j,i) = 0
            mic_weights(j,i) = 0.0D0
            invalid_mic_nodes = invalid_mic_nodes + 1
          END IF
        CASE(7)
          IF (state(x_start,y_start+1,z_start+1,1) >= 0) THEN
            mic_nodes(1,j,i) = x_start
            mic_nodes(2,j,i) = y_start+1
            mic_nodes(3,j,i) = z_start+1
!        mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start+1)
!
! Find the local distance for weighting purposes
!
            local_dist(7) = SQRT((mic_locations(1,i)-loc_lo(1))**2+&
            (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
            (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)

            total_dist = total_dist+local_dist(7)
          ELSE
!        mic_links(j,i) = 0
            mic_weights(j,i) = 0.0D0
            invalid_mic_nodes = invalid_mic_nodes + 1

          END IF
        CASE(8)
          IF (state(x_start+1,y_start+1,z_start+1,1) >= 0) THEN
            mic_nodes(1,j,i) = x_start+1
            mic_nodes(2,j,i) = y_start+1
            mic_nodes(3,j,i) = z_start+1
!        mic_links(j,i) = state(dir+1,x_start+1,y_start+1,&
!          z_start+1)
!
! Find the local distance for weighting purposes
!
            local_dist(8) = SQRT((mic_locations(1,i)-loc_lo(1)+rdx(mic_cur_lvl))**2+&
            (mic_locations(2,i)-loc_lo(2)+rdx(mic_cur_lvl))**2+&
            (mic_locations(3,i)-loc_lo(3)+rdx(mic_cur_lvl))**2)

            total_dist = total_dist+local_dist(8)
          ELSE
!        mic_links(j,i) = 0
            mic_weights(j,i) = 0.0D0
            invalid_mic_nodes = invalid_mic_nodes + 1

          END IF
        END SELECT
      END DO

      IF (invalid_mic_nodes == 8) THEN
        mic_valid(i) = .FALSE.
        WRITE(11,101) i
        WRITE(11,100) i
        CYCLE mic_setup_loop
      END IF

      WRITE(11,*) 'Finished determining microphone distances.'
!      WRITE(*,*) 'Finished determining microphone distances.'
      DO j = 1,8
        if (mic_nodes(1,j,i)> 0) then
!    IF (mic_links(j,i) > 0) THEN
          mic_weights(j,i) = local_dist(j)/total_dist
!          WRITE(*,*) self,mic_weights(j,i),local_dist(j),total_dist
        ELSE
          mic_weights(j,i) = 0.0D0
        END IF

      END DO
!
! MPI wait for all the sends to finish
!
!      call mpi_waitall(nprocs-1,req_array,stat_array,ier)

      WRITE(11,*) 'Finished determining microphone weighting.'
      cycle mic_setup_loop
    else

      CYCLE
    end if
  end do

!tag = 1000*i+self
!!write(*,*) 'receiving things!',self
!call mpi_recv(mic_owner(i),1,MPI_INT,MPI_ANY_SOURCE,tag,&
!  commune,taste,ier)
END DO mic_setup_loop

!write(*,*) 'hitting the barrier',self
call mpi_barrier(commune,ier)

call mic_mpi_blast
!write(*,*) 'mics owned by',self,mic_owner
CALL mic_initialize
!      CALL SYSTEM('cd .',ier)
!      IF (ier /=0) THEN
!        WRITE(*,*) 'Something went wrong  microphone directory.'
!        WRITE(11,*) 'Something went wrong  microphone directory.'
!        CALL error_out
!      END IF
!      WRITE(*,*) 'Microphone initialization done.'
WRITE(11,*) 'Microphone initialization done.'
!WRITE(*,*) 'Microphone initialization done.',self
 100  FORMAT ('Microphone ',I5,' is not in a valid location')
 101  FORMAT ('Microphone ',I5,' is surrounded by invalid nodes',&
        ' and is therefore invalid.')

END SUBROUTINE


!  invalid_mic_nodes = 0
!  total_dist = 0.0D0
!  DO j = 1,8
!!
!! Go through all 8 possible nodes to be able to pull from
!!
!    SELECT CASE (j)
!    CASE(1)
!      IF (state(x_start,y_start,z_start,1) >= 0) THEN
!!        mic_links(j,i) = state(dir+1,x_start,y_start,z_start)
!  !
!! Find the local distance for weighting purposes
!!
!        local_dist(1) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!        y_start,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start,&
!        y_start,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!        y_start,z_start))**2.0D0)
!
!!
!        total_dist = total_dist + local_dist(1)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!!
!    CASE(2)
!      IF (state(x_start+1,y_start,z_start,1) >= 0) THEN
!        mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(2) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!        y_start,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start+1,&
!        y_start,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start+1,&
!        y_start,z_start))**2.0D0)
!
!        total_dist = total_dist+local_dist(2)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!    CASE(3)
!      IF (state(1,x_start,y_start+1,z_start) == 0) THEN
!        mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(3) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!        y_start+1,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start,&
!        y_start+1,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!        y_start+1,z_start))**2.0D0)
!
!        total_dist = total_dist+local_dist(3)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!
!    CASE(4)
!      IF (state(1,x_start,y_start,z_start+1) == 0) THEN
!        mic_links(j,i) = state(dir+1,x_start,y_start,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(4) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!        y_start,z_start+1))**2.0D0+(mic_locations(2,i)-grid(2,x_start,&
!        y_start,z_start+1))**2.0D0+(mic_locations(3,i)-grid(3,x_start,&
!        y_start,z_start+1))**2.0D0)
!
!        total_dist = total_dist+local_dist(4)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!    CASE(5)
!      IF (state(1,x_start+1,y_start+1,z_start) == 0) THEN
!        mic_links(j,i) = state(dir+1,x_start+1,y_start+1,z_start)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(5) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!        y_start+1,z_start))**2.0D0+(mic_locations(2,i)-grid(2,x_start+1,&
!        y_start+1,z_start))**2.0D0+(mic_locations(3,i)-grid(3,x_start+1,&
!        y_start+1,z_start))**2.0D0)
!
!        total_dist = total_dist+local_dist(5)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!    CASE(6)
!      IF (state(1,x_start+1,y_start,z_start+1) == 0) THEN
!        mic_links(j,i) = state(dir+1,x_start+1,y_start,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(6) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!        y_start,z_start+1))**2.0D0+(mic_locations(2,i)-grid(2,&
!        x_start+1,y_start,z_start+1))**2.0D0+(mic_locations(3,i)-&
!        grid(3,x_start+1,y_start,z_start+1))**2.0D0)
!
!        total_dist = total_dist+local_dist(6)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!      END IF
!    CASE(7)
!      IF (state(1,x_start,y_start+1,z_start+1) == 0) THEN
!        mic_links(j,i) = state(dir+1,x_start,y_start+1,z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(7) = SQRT((mic_locations(1,i)-grid(1,x_start,&
!        y_start+1,z_start+1))**2.0D0+(mic_locations(2,i)-grid(2,&
!        x_start,y_start+1,z_start+1))**2.0D0+(mic_locations(3,i)-&
!        grid(3,x_start,y_start+1,z_start+1))**2.0D0)
!
!        total_dist = total_dist+local_dist(7)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!
!      END IF
!    CASE(8)
!      IF (state(1,x_start+1,y_start+1,z_start+1) == 0) THEN
!        mic_links(j,i) = state(dir+1,x_start+1,y_start+1,&
!          z_start+1)
!!
!! Find the local distance for weighting purposes
!!
!        local_dist(8) = SQRT((mic_locations(1,i)-grid(1,x_start+1,&
!        y_start+1,z_start+1))**2.0D0+(mic_locations(2,i)-grid(2,&
!        x_start+1,y_start+1,z_start+1))**2.0D0+(mic_locations(3,i)-&
!        grid(3,x_start+1,y_start+1,z_start+1))**2.0D0)
!
!        total_dist = total_dist+local_dist(8)
!      ELSE
!        mic_links(j,i) = 0
!        mic_weights(j,i) = 0.0D0
!        invalid_mic_nodes = invalid_mic_nodes + 1
!
!      END IF
!    END SELECT
!  END DO
!   END IF
