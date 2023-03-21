SUBROUTINE namelist_box_of_interest
!
! Finds what level and the lower and upper corners of the box of interest
! is. The BoI is a region where the grid cannot coarsen, only refine.
!
! Called by: namelist_global
! Calls: error_out
!
USE precise
use amr_info_holder
use timing
IMPLICIT NONE
INTEGER :: istat,box_level,i,coarse_level
REAL(kind=dp) :: box_x_minimum,box_y_minimum,box_z_minimum,&
  box_x_maximum,box_y_maximum,box_z_maximum

NAMELIST /BOXES/ box_level,coarse_level,box_x_minimum,&
  box_y_minimum,box_z_minimum,box_x_maximum,box_y_maximum,&
  box_z_maximum

!write(*,*) 'Ooooooo spookey ghost!!!!'

coarse_regrid_done = .true.

box_level = 3
coarse_level = 1
box_x_minimum = -1.0
box_y_minimum = -1.0
box_z_minimum = -1.0
box_x_maximum = 1.0
box_y_maximum = 1.0
box_z_maximum = 1.0
do i = 1,num_boi
  READ(2,BOXES,IOSTAT=istat)
  IF (istat .NE. 0) THEN
    WRITE(*,*) 'Error reading namelist at box of interest data read.'
    WRITE(*,*) istat
    CALL error_out
  END IF

  boi_lo(1,i) = box_x_minimum
  boi_lo(2,i) = box_y_minimum
  boi_lo(3,i) = box_z_minimum
  boi_hi(1,i) = box_x_maximum
  boi_hi(2,i) = box_y_maximum
  boi_hi(3,i) = box_z_maximum
  coarse_box(i) = coarse_level
  box_lvl(i) = box_level

  write(11,*) 'Box of interest data, box number ',i
  write(11,*) 'Grid will not be coarsened within this box.'
  write(11,500) box_level
  write(11,507) coarse_level
  write(11,501) box_x_minimum
  write(11,502) box_y_minimum
  write(11,503) box_z_minimum
  write(11,504) box_x_maximum
  write(11,505) box_y_maximum
  write(11,506) box_z_maximum
  write(11,*) ''

  if (box_lvl(i) > coarse_box(i) .and. coarse_time > very_small_number .and. &
      coarse_regrid_done) then
    coarse_regrid_done = .false.
  end if

end do


 500  format ('Box of interest, minimum grid level = ',I2)
 501  format ('Box of interest, minimum x = ',F11.5)
 502  format ('Box of interest, minimum y = ',F11.5)
 503  format ('Box of interest, minimum z = ',F11.5)
 504  format ('Box of interest, maximum x = ',F11.5)
 505  format ('Box of interest, maximum y = ',F11.5)
 506  format ('Box of interest, maximum z = ',F11.5,/)
 507  format ('Coarse time box level = ',I2)
end subroutine
