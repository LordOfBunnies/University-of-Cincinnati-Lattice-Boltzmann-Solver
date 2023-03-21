subroutine calcs_before_init
!
! Need the number of directions and the Mach number before initializing Amrex
!
!
! Called by: uclbs_main
! Calls: directions_and_weights
!
use precise
use constants
USE freestream_values
use fill_stuff, only: allocate_fi_bcs
implicit none
integer :: istat
character(len=80) :: filename_out
filename_out = 'uclbs.out'
!
open(FILE=filename_out,UNIT = 11,STATUS='OLD',FORM = 'FORMATTED',&
  ACTION='WRITE',POSITION='APPEND',ACCESS = 'SEQUENTIAL',IOSTAT=istat)
if (istat /= 0) THEN
  write(*,*) 'Problems opening output file, closing.'
  write(*,*) istat
  call error_out
end if
!
! Find the Mach number so the number of directions can be found before Amrex init
!
gas_const_R = 287.1
t_ref = temperature
speed_of_sound = SQRT(gama*gas_const_R*t_ref)
u_ref = speed_of_sound
!
! Basic LBM simulations are isothermal, this changes the speed of sound so
! this is something to return to later. For now, the velocities have been low
! because of this fact. XXXX
!
!write(*,*) v_inf
mach = v_inf/u_ref
!mach = 0.5D0
call directions_and_weights

!mach = v_inf/u_ref

call allocate_fi_bcs
!
!
!
write(11,*) 'Initial calculations done.'
write(11,105) t_ref
write(11,106) u_ref
write(11,103) speed_of_sound
write(11,104) mach
write(11,*) ''

!write(*,*) 'Initial calculations done.'


close(11)

 103  FORMAT ('Speed of sound = ',F11.5)
 104  FORMAT ('Mach number = ',F11.5)
 105  format ('Reference temperature = ', F8.4)
 106 format ('Reference velocity = ',F8.4)
end subroutine
