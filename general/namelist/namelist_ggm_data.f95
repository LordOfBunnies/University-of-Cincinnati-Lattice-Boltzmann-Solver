SUBROUTINE namelist_ggm_data
!
! Reads the GGM data from the NAMELIST file
!
! Called by: namelist_global
! Calls: error_out
! XXXX
USE precise
USE ggm_stuff
IMPLICIT NONE
INTEGER :: order_of_accuracy,num_ggm_variables,ggm_direction
INTEGER :: ggm_variable,istat, i,min_ggm_lvl,max_ggm_lvl,pgs_reach
real(kind=dp) :: freq_ggm_output,ggm_interval,ggm_sense
logical :: output_ggm,upwind_calculation,downwind_calculation,central_diff,&
  pgs_on


NAMELIST /GGM_INFO/ order_of_accuracy,num_ggm_variables,&
  ggm_direction,upwind_calculation,downwind_calculation,central_diff,&
  output_ggm,ggm_interval,ggm_sense,pgs_reach,pgs_on

NAMELIST /GGM_DATA_STUFF/ ggm_variable

NAMELIST /GGM_OUTPUT_INFO/ freq_ggm_output,min_ggm_lvl,max_ggm_lvl

order_of_accuracy = 2
num_ggm_variables = 0
ggm_direction = 1
pgs_reach = 8
ggm_interval = 0.1D-4
ggm_sense = 0.03D0
upwind_calculation = .FALSE.
downwind_calculation = .FALSE.
central_diff = .TRUE.
output_ggm = .true.
pgs_on = .false.
!!
!! Read the baseline ggm information in
!!
READ(2,NML=GGM_INFO,IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Error reading GGM information, exiting.'
  WRITE(11,*) 'Error reading GGM information, exiting.'
  WRITE(*,*) istat
  CALL error_out
END IF
!!
!! Error check to make sure that
!!
IF (num_ggm_variables <= 0 .OR. num_ggm_variables >= 10) THEN
  use_ggm = .FALSE.
  RETURN
END IF

ggm_order = order_of_accuracy
num_ggm_vars = num_ggm_variables
ggm_dir = ggm_direction
pgs_max_scan = pgs_reach
upwind = upwind_calculation
downwind = downwind_calculation
central = central_diff
ggm_output = output_ggm
ggm_dt = ggm_interval
ggm_sensitivity = ggm_sense
use_pgs = pgs_on

if (.not. allocated(ggm_variables)) then
  ALLOCATE(ggm_variables(num_ggm_variables))
end if
if (ggm_output) then
  allocate(ggm_output_freq(num_ggm_variables))
  allocate(ggm_min_lvl(num_ggm_variables))
  allocate(ggm_max_lvl(num_ggm_variables))
end if
!!
!! Read the needed ggm variables in
!!
!! GGM variables as follows
!! 1 = density
!! 2 = u-velocity
!! 3 = v-velocity
!! 4 = w-velocity
!! 5 = Temperature
!! 6 = Pressure
!! 7 = Velocity magnitude
!! 8 = Dynamic pressure
!! 9 = Mach
!!
!! 10 = x-Vorticity
!! 11 = y-Vorticity
!! 12 = z-Vorticity
!! 13 = Vorticity magnitude
!!
!! 14 = xy-Shear Stress
!! 15 = xz-Shear Stress
!! 16 = yz-Shear Stress
!! 17 = Shear Stress magnitude
!!
!! 20 = alpha
!!
!! 29 = Q-criterion
!!
!!
DO i=1,num_ggm_vars
  ggm_variable = 0
  READ(2,NML=GGM_DATA_STUFF,IOSTAT=istat)
  IF (istat /= 0) THEN
    WRITE(*,*) 'Error reading GGM variable data, exiting.'
    WRITE(11,*) 'Error reading GGM variable data, exiting.'
    WRITE(*,*) istat
    CALL error_out
  END IF

  if (ggm_output) then

    freq_ggm_output = 0.1D-4
    min_ggm_lvl = 0
    max_ggm_lvl = 1
    read(2,NML=GGM_OUTPUT_INFO,IOSTAT=istat)
    if (istat /= 0) then
      WRITE(*,*) 'Error reading GGM output information, exiting.'
      WRITE(11,*) 'Error reading GGM output information, exiting.'
      WRITE(*,*) istat
      CALL error_out
    end if

    ggm_output_freq(i) = freq_ggm_output
    ggm_max_lvl(i) = max_ggm_lvl
    ggm_min_lvl(i) = min_ggm_lvl

  end if

  ggm_variables(i) = ggm_variable

END DO
!!
!! Error checker to make sure all the variables read in were valid
!!
DO i = 1,num_ggm_vars

  IF (ggm_variables(i) <= 0 .OR. ggm_variables(i) >= 30) THEN

    use_ggm = .FALSE.
    RETURN

  END IF

END DO

if (gradient_refine_only) then
  num_ggm_vars = 1
  if (allocated(ggm_variables)) deallocate(ggm_variables)
  allocate(ggm_variables(num_ggm_vars))

  ggm_variables(1) = gradient_selection
end if

END SUBROUTINE
