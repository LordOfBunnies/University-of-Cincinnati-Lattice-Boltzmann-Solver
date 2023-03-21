MODULE timing
!
! This holds the values for the times necessary to run the program
!
!
!
!
USE precise
IMPLICIT NONE
SAVE
REAL(KIND=dp) :: total_time,coarse_time,comp_time,cpu_start_time,&
  cpu_stop_time,ggm_start_time,ggm_stop_time
real(kind=dp) :: init_beta_1,init_beta_2,beta_1_scaling,beta_2_scaling
integer :: startup_timesteps



END MODULE
