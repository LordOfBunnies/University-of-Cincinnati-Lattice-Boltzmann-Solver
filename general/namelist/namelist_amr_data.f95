SUBROUTINE namelist_amr_data
!
! This reads in relevant information for the AMR data
!
! Called by: namelist_global
! Calls: error_out, namelist_box_of_interest
!
USE precise
use amr_info_holder
use ggm_stuff
IMPLICIT NONE
integer :: num_boxes_of_interest,istat,select_gradient
real(kind=dp) :: grad_limit
!      INTEGER :: refinement_interval,maximum_size,&
!        blocking_factor,istat,box_level,ref_ratio
!      !INTEGER,ALLOCATABLE :: refinement_ratio(:)
!!      REAL(kind=dp) :: box_x_minimum,box_y_minimum,box_z_minimum,&
!!        box_x_maximum,box_y_maximum,box_z_maximum
!      LOGICAL :: x_periodic,y_periodic,z_periodic,box_of_interest
!Declare the namelist
NAMELIST /AMR_DATA/ num_boxes_of_interest,grad_limit,select_gradient

! Initialize the variables
num_boxes_of_interest = 1
grad_limit = 0.1D0
select_gradient = 7
!
! Read information from the NAMELIST file
!
READ(2,AMR_DATA,IOSTAT=istat)
IF (istat .NE. 0) THEN
  WRITE(*,*) 'Error reading namelist at AMR data read.'
  WRITE(*,*) istat
  CALL error_out
END IF
num_boi = num_boxes_of_interest
gradient_limit = grad_limit
gradient_selection = select_gradient
!
!write(*,*) 'Number of boxes of interest',num_boi
allocate(boi_lo(3,num_boi))
allocate(boi_hi(3,num_boi))
allocate(box_lvl(num_boi))
allocate(coarse_box(num_boi))

!coarse_refine = .false.
! Write information to the output file.
!
WRITE(11,*) 'AMR data for the run.'
write(11,101) num_boi
!
! Apply the read information into the amr variables
!
!      CALL amr_transfer_data_to_amrex(maximum_level, refinement_interval,&
!        maximum_size,blocking_factor, refinement_ratio,box_level,&
!        box_x_minimum,box_y_minimum,box_z_minimum,box_x_maximum,&
!        box_y_maximum,box_z_maximum,x_periodic,y_periodic,z_periodic)
!
IF (num_boxes_of_interest > 0) then
  CALL namelist_box_of_interest
end if
!WRITE(11,300) maximum_level
!WRITE(11,301) refinement_interval
!WRITE(11,302) maximum_size
!WRITE(11,303) blocking_factor
!WRITE(11,304) refinement_ratio
!
!WRITE(11,601) x_periodic
!WRITE(11,602) y_periodic
!WRITE(11,603) z_periodic

if (gradient_refine_only) then
  num_ggm_vars = 1
  if (allocated(ggm_variables)) deallocate(ggm_variables)
  allocate(ggm_variables(num_ggm_vars))

  ggm_variables(1) = gradient_selection
end if

 300  FORMAT ('Maximum level for the AMR = ',I3)
 301  FORMAT ('Mesh refined every ',I4,' steps')
 302  FORMAT ('Maximum size of grids for a processor = ', I4)
 303  FORMAT ('Minimum size of grid when refining = ',I4)
 304  FORMAT ('Ratio of nodes between levels = ',I4)

 404  FORMAT ('Error reading namelist at AMR data input read type ',I4)

 601  FORMAT ('Domain is periodic in the x-direction? ',L2)
 602  FORMAT ('Domain is periodic in the y-direction? ',L2)
 603  FORMAT ('Domain is periodic in the z-direction? ',L2)
 101  format ('Number of boxes of interest ',I3)

end subroutine
