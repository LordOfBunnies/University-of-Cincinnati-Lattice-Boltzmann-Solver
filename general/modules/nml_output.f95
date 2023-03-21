MODULE nml_output
!
! Saves information about the output options for the program
!
! mic_locations defines where mics are (x/y/z,#)
! mic_times define when the mics are active (start/stop,#)
! mic_type defines what the mics record
! mic_frequency defines how often they record (every X timesteps)
! mic_valid defines whether the microphone is actually active in the volume
!
! restarting is defined in the NAMELIST file and is used if the run is a restart
! write_restart_at_end tells the program to write a restart file at the end 
!   of the run or not
!
! save_time define when the program records info (start/stop,#)
! save_regions defines the boundaries of the regions to save
! restart_saves defines what times to output restart saves at
! save_shapes defines the highest and lowest levels as the output range
! save_frequency defines how often they record (every X timesteps)
!
USE precise
IMPLICIT NONE
SAVE
INTEGER :: num_mics,num_restart_saves,num_saves
LOGICAL :: restarting,write_restart_at_end

REAL(KIND=dp), ALLOCATABLE :: mic_locations(:,:),mic_times(:,:),save_frequency(:)
!read in
integer,allocatable :: mic_frequency(:),mic_lvl(:)
! found later
INTEGER, ALLOCATABLE :: mic_nodes(:,:,:),mic_owner(:)

REAL(KIND=dp), ALLOCATABLE :: save_times(:,:),restart_saves(:)!,save_regions(:,:) !XXXX may be needed later
INTEGER, ALLOCATABLE :: save_shapes(:,:)
LOGICAL, ALLOCATABLE :: save_data(:,:),mic_type(:,:),mic_valid(:),lm_saves(:,:)

!CONTAINS

!
! multigrid_save_info is to save required information for CGNS adpated mesh output.  It is declared allocatable
!   so that an instance can exist for every timestep
!
! cumulative_zones holds the information on how many total zones are written out.  This doesn't
!      increment if an adaptation has not occured since last output
! multigrid_counter holds the zone number for output purposes
! zones_per_step holds the number of zones that are output at each timestep
! amr_solution_names holds the FlowSolutionPointers for output in ZoneIterativeData.  When a
!       zone is active, a solution name is placed here.  When not used, Null is placed here instead.
! zone_point_names holds vectors of active zone names in BaseIterativeData.
! regrid_between_outputs says whether a regrid occured since the last output, if false, the same zones
!       can be used and new zone outputs do not need to be used.
! zone_used holds whether a zone was active during an output phase, and thus whether amr_solution_names
!       needs to be filled for that instance.
!
type,public :: base_data
  integer :: max_zones,adapted_zone_start_number
  integer,allocatable :: total_zones_at_ts(:),zones_per_lvl(:,:)
  character(len=32),allocatable ::  zone_pointer_names(:,:),amr_solution_names(:,:)
end type

type,public :: multigrid_save_info
  integer,allocatable :: multigrid_counter(:),zone_number(:),zone_lvls(:)
!  character(len=32),allocatable ::
  logical :: regrid_between_outputs
end type
!
! saved_timesteps are timesteps on the max level for a given where output is done
! parallel_bounds are bounds for data output so array writes don't overlap
!   (1,) is the minimum bound,(2,) is the maximum bound
! saved_times are floating times for output based on start time and frequency
! n_outsteps are the total number of outputs that should be done for each requested output
!
type,public :: timing_info_steps
  integer,allocatable :: saved_timestep(:)
  integer :: n_outsteps,output_counter
  real(kind=dp),allocatable :: saved_times(:)
  character(len=32),allocatable :: sol_names(:)
  logical,allocatable :: write_done(:)
end type


END MODULE
