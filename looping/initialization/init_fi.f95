      SUBROUTINE init_fi
!
! This will initialize all the major arrays for the first iteration
! of the looping procedure.  It will also take into account all the
!
! Called by: uclbs_main
! Calls:
!
USE linkwise
USE nml_init
USE grid_data
USE freestream_values
USE precise
USE constants
!      USE output_formatting
IMPLICIT NONE
INTEGER :: i,j,highest_priority,overlap_counter,istat
REAL(KIND = dp) :: cu
REAL(KIND = dp) :: vel_mag
LOGICAL :: inside

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
! ALLOCATE all the important things before we start working with them.
!
WRITE(11,*) ''
WRITE(11,*) 'Begin initializing grid.'
!      WRITE(11,*)
WRITE(11,*) ''


call multilevel_initialization


CLOSE(11)
END SUBROUTINE


!      ALLOCATE(u_vel(link_number))
!      ALLOCATE(v_vel(link_number))
!      ALLOCATE(w_vel(link_number))
!      ALLOCATE(rho(link_number))
!      ALLOCATE(fi(dir,link_number))
!      ALLOCATE(fout(dir,link_number))

      !      rho = rho_inf
!
!      IF (primary_flow_dir == 1) THEN
!        u_vel = v_inf
!        v_vel = 0.0D0
!        w_vel = 0.0D0
!      ELSE IF (primary_flow_dir == 2) THEN
!        u_vel = 0.0D0
!        v_vel = v_inf
!        w_vel = 0.0D0
!      ELSE IF (primary_flow_dir == 3) THEN
!        u_vel = 0.0D0
!        v_vel = 0.0D0
!        w_vel = v_inf
!      ELSE IF (primary_flow_dir == 4) THEN
!        u_vel = -v_inf
!        v_vel = 0.0D0
!        w_vel = 0.0D0
!      ELSE IF (primary_flow_dir == 5) THEN
!        u_vel = 0.0D0
!        v_vel = -v_inf
!        w_vel = 0.0D0
!      ELSE IF (primary_flow_dir == 6) THEN
!        u_vel = 0.0D0
!        v_vel = 0.0D0
!        w_vel = -v_inf
!      END IF
!!
!!
!!
!      CALL amr_initialization
!!
!!
!!
!      DO i = 1,link_number
!        highest_priority = 0
!        overlap_counter = 0
!        DO j = 1,num_init_regions
!
!          CALL inside_init(i,j,highest_priority,inside,overlap_counter)
!
!        END DO
!      END DO
!
!      WRITE(11,*) 'Determined nodes inside of initialization regions.'
!!
!! Assign
!!
!!      DO j = 1,link_number
!!        DO i = 1,dir
!!
!!          cu = 3.0D0*(cx(i)*u_vel(j)+cy(i)*v_vel(j)+cz(i)*w_vel(j))
!!          vel_mag = (u_vel(j)**2.0D0+v_vel(j)**2.0D0+w_vel(j)**2.0D0)
!!          fi(i,j) = wi(i)*rho(j)*(1.0D0 + cu + 0.5D0*cu**2.0D0 - 1.5D0*vel_mag)
!!!          feq(i,j) = wi(j)*rho(j)*(1 + cu + .5*cu**2 - 1.5*vel_mag)
!!        END DO
!!      END DO
!
!      CALL local_initialization
!!      WRITE(*,*) 'Cu, vel_mag', cu, vel_mag
!!      WRITE(*,*) 'Init to fi to ',fi(1,89306),fi(2,89306),fi(3,89306),&
!!        fi(4,89306),fi(5,89306),fi(6,89306),fi(7,89306),fi(8,89306),fi(9,89306)
!
!      WRITE(11,*) ''
!      WRITE(11,*) 'Finished initializing grid, moving on to timestepping.'
!      WRITE(11,*) ''
!      CALL output_setup
