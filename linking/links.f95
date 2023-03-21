SUBROUTINE links
!
! This will create the link lists that will be used once timestepping
! starts.  First every active node will be assigned a number, then the
! number will be associated to it's neighbors via it's positions on the
! link list. It will also say whether this link is fluid, solid, or
! boundary so things down the road will be able to correctly apply bounding conditions.
!
! Called by: uclbs_main
! Calls: mic_setup,boundary_links,output_setup,periodic_linking,error_out
!
use linkwise
use grid_data
use nml_inlet_outlet
use freestream_values
use constants
use precise
!use amr_info_holder
use amrex_base_module
use amrex_amr_module
use amr_processes,only:self
use cgns
IMPLICIT NONE
!      REAL ::
!      INTEGER,INTENT(IN) :: dimensions
INTEGER :: a,b,c,i,num_valid_nodes,link_error_checker,ier,index_file
INTEGER :: istat,lvl!,cox_lo(3),cox_hi(3),loc_dir!,weald
!      INTEGER :: valid_links(num_x,num_y,num_z)

CHARACTER(len=80) :: filename_out,filename

!type(amrex_mfiter) :: mfi
!type(amrex_box) :: cox

filename_out = 'uclbs.out'
!
OPEN(FILE=filename_out,UNIT = 11,STATUS='OLD',FORM = 'FORMATTED',&
  ACTION='WRITE',POSITION='APPEND',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Problems opening output file, closing.'
  WRITE(*,*) istat
  CALL error_out
END IF

WRITE(11,*) 'Starting linking.'
!
! Set up microphones
!
CALL mic_setup
WRITE(11,*) 'Microphones successfully set up.'
!WRITE(*,*) 'Microphones successfully set up.',self
!
! Find all the links at the boundaries and next to walls
!
!      CALL boundary_links
!
! Set up all the output conditions for later export via CGNS
!
CALL output_setup
!
CLOSE(11)

! XXXX temporary removal of subroutine calls to be reinstated after testing.
!CALL state_and_grid_output

!CALL ggm_output_setup

 200  FORMAT ('Total number of active nodes = ',I8)
 201  FORMAT ('Confirmed total of active nodes  =', I8)



! 1001 FORMAT (19I11)
END SUBROUTINE
!fine = amrex_get_finest_level()
!do lvl = fine,0,-1
!  CALL amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
!  do while(mfi%next())
!
!  cox = valid_box()
!  cox_lo = cox%lo
!  cox_hi = cox%hi
!  do loc_dir = 1,dir
!
!    z_loop: DO c = cox_lo(3),cox_hi(3)
!      y_loop: DO b = cox_lo(2),cox_hi(2)
!        x_loop: DO a = cox_lo(1),cox_hi(1)
!
! Make sure no nodes with no valid connections are included in final linking
! This works on the null nodes outside the bounding region
!
!            IF (state(1,a,b,c) == -1) THEN
!              null_dir_counter = 1
!              DO i = 2,dir
!                IF (state(i,a,b,c) == -1 .OR. state(i,a,b,c) &
!                  >= 1) THEN
!                  null_dir_counter = null_dir_counter + 1!
!
!                END IF

!              END DO
!              IF (null_dir_counter /= dir) THEN
!                state(dir+1,a,b,c) = link_number
!                link_number = link_number + 1
!              ELSE
!                state(dir+1,a,b,c) = 0
!              END IF
!!
!! This works on nodes declared solid
!!
!            ELSE IF (state(1,a,b,c) == 1) THEN
!              null_dir_counter = 1
!              DO i = 2,dir
!                IF (state(i,a,b,c) == -1 .OR. state(i,a,b,c) &
!                  >= 1) THEN
!                  null_dir_counter = null_dir_counter + 1
!
!                END IF
!
!              END DO
!              IF (null_dir_counter /= dir) THEN
!                state(dir+1,a,b,c) = link_number
!                link_number = link_number + 1
!              ELSE
!                state(dir+1,a,b,c) = 0
!              END IF
!!
!! This assigns the link number to the extra spot
!!
!      IF (state(1,a,b,c) == 0) THEN
!        state(a,b,c,dir+1) = link_number
!        link_number = link_number + 1
!      END IF
!          state(a,b,c,loc_dir) = state(a+cx(loc_dir),b+cy(loc_dir),&
!            c+cz(loc_dir),1)
!
!
!        END DO x_loop
!      END DO y_loop
!    END DO z_loop
!  end do
!end do

!z_loop_2: DO c = 1,num_z
!  y_loop_2: DO b = 1,num_y
!    x_loop_2: DO a = 1,num_x
!!
!! Weed out the null nodes as they were shown above then assign the
!!
!      IF (state(1,a,b,c) /= 0) THEN
!!              link_error_checker = link_error_checker + 1
!!              WRITE(*,*) link_error_checker
!        CYCLE x_loop_2
!      ELSE
!!              WRITE(*,*) link_error_checker,a,b,c,state(dir+1,a,b,c),state(1,a,b,c)
!        IF (link_error_checker == state(dir+1,a,b,c)) THEN
!!
!! First group of links assignments is for the node itself if it ever
!! becomes important for it to self-identify
!!
!!                WRITE(*,*) 'Sparkle',link_error_checker
!           link(1,link_error_checker) = link_error_checker
!!                link(1,link_error_checker) = state(1,a,b,c)
!          dir_loop: DO i = 2,dir
!!
!! Second declaration is taking data already determined in state and put
!! it into the link list
!!
!! links(direction, link #, 1) = link number in that direction
!! links(direction, link #, 2) = state of the link (fluid, solid, etc.)
!!
!! link states - 0 = fluid
!! 1 = solid
!! -1 = out of bounds
!! 1000 = freestream values bounds
!! 10X = inlet BC
!! 20X = outlet BC
!!
!! Wrap this in an if statement dealing with boundary conditions on those sides
!! It should be able to circumvent this making all edges null bcs
!!
!
!              IF (a+cx(i) < 1 .OR. a+cx(i) > num_x .OR. b+cy(i)&
!                <1 .OR. b+cy(i) >num_y .OR. c+cz(i) < 1 .OR. &
!                c+cz(i) > num_z .OR. state(i,a,b,c) /= 0) THEN
!
!
!                link(i,link_error_checker) = state(i,a,b,c)
!
!              ELSE
!                link(i,link_error_checker) = state(dir+1,a+&
!                  cx(i),b+cy(i),c+cz(i))
!              END IF
!!
!!
!!
!              IF (state(i,a,b,c) >=200 .AND. &
!                    state(i,a,b,c) <=299) THEN
!                IF (outlet_type(1,state(i,a,b,c)-200) == 5) THEN
!                  CALL periodic_linking(i,a,b,c,link_error_checker,&
!                    outlet_type(2,state(i,a,b,c)-200))
!                END IF
!
!              END IF
!!
!!
!          END DO dir_loop
!!
!! Put the geometry into an appropriate link list as well
!!
!          link_grid(1,link_error_checker) = grid(1,a,b,c)
!          link_grid(2,link_error_checker) = grid(2,a,b,c)
!          link_grid(3,link_error_checker) = grid(3,a,b,c)
!
!        ELSE
!          WRITE(*,*) 'Something has gone wrong in link assignment.'
!          CALL error_out
!        END IF
!!              WRITE(*,*) link_error_checker
!        link_error_checker = link_error_checker + 1
!      END IF
!
!    END DO x_loop_2
!  END DO y_loop_2
!END DO z_loop_2
!
!link_error_checker = link_error_checker - 1
!
!WRITE(11,201) link_error_checker


!  filename = 'test.cgns'
!
!  CALL cgp_open_f(TRIM(filename),CG_MODE_WRITE,index_file,ier)
!  IF (ier /= CG_OK) CALL cgp_error_exit_f
!  WRITE(11,*) 'Open CGNS file done'
!!  write(*,*) 'test file opened',self
!
!  CALL cgp_close_f(index_file,ier)
!  IF (ier /= CG_OK) CALL cgp_error_exit_f
!!  write(*,*) 'test file finished',self





