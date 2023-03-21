subroutine mic_mpi_blast
!
!
!
!
!
!
!
use precise
use nml_output
use output_data
use mpi
use amr_processes, only: nprocs,self,commune
implicit none
integer :: i,j,buff_count,ier
integer :: taste(mpi_status_size)
integer,allocatable :: req_array(:)



!do i = 1,num_mics
!
!!  buff_count =
!  call mpi_bcast(mic_nodes(1:3,1:8,i),24,MPI_INT,mic_owner(i),commune,ier)
!  call mpi_bcast(mic_weights(1:8,i),8,MPI_DOUBLE_PRECISION,mic_owner(i),commune,ier)
!
!end do

end subroutine
