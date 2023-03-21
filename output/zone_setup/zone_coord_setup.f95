subroutine zone_coord_setup(index_file,index_base,index_zone)
!
! Write the coordinates to the output file
!
! Called by: zone_writer
! Calls:
! External calls: cgp_coord_write_f
!
use cgns
implicit none
integer :: coord_x,coord_y,coord_z,ier
integer :: index_file,index_base,index_zone

coord_x = 1
coord_y = 2
coord_z = 3
!index_coord = 1
!  write(*,*) 'Beginning coordinate information write',index_file,index_base,index_zone,index_coord
  CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
   'CoordinateX',coord_x,ier)
!    write(*,*) self,ier
!    call cg_error_print_f
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!write(*,*) 'X-coordinate base !written',index_file,index_base,index_zone,coord_x

  CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
   'CoordinateY',coord_y,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!write(*,*) 'Y-coordinate base written',index_file,index_base,index_zone,coord_y

CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
 'CoordinateZ',coord_z,ier)
!write(*,*) 'Z-coordinate base written',index_file,index_base,index_zone,coord_z
!bz,coord_x,coord_y,coord_z
IF (ier /= CG_OK) CALL cgp_error_exit_f

end subroutine
