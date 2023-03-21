subroutine background_calcs(characteristic_length)
!
! This is just a sorting ground for background calculations
!
!
! Called by: uclbs_main
! Calls: basic_fluid_calcs, directions_and_weights, normalize_values
!
use precise
implicit none
real(KIND=dp) :: characteristic_length
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
! Do basic flow calculations
!
call basic_fluid_calcs
!
! Set up all the possible flow directions and their weightings
!
!call directions_and_weights
!
!
!
call viscosity_and_timing(characteristic_length)
!
! Normalize all the flow values before going further
!
call normalize_values(characteristic_length)
!
!
call freestream_basis_calcs

call inlet_outlet_basis_calcs

!
write(11,*) 'Background calculations done.'
write(11,*) ''
!write(*,*) 'Background calculations done.'
close(11)
!


      END SUBROUTINE
