module ggm_amr_writing
!
!
!
!
!
!
!
use precise
use constants
use amr_info_holder
use amr_processes
use grid_data
use cgns
use mpi
use output_data
use amrex_base_module
use amrex_amr_module
implicit none
save
private !:: amr_big_zones
public :: ggm_amr_grid_zones, ggm_amr_grid_distro, ggm_amr_grid_ranges
public :: ggm_amr_grid_zone_writer, ggm_amr_grid_coord_setup, ggm_amr_save_writer
public :: ggm_amr_grid_output_nums,ggm_amr_grid_coord_data_write




contains
!
!
!
!
!
!
! This is zone_box_setup from the zone_writer module
subroutine ggm_amr_grid_zones(lvl,num)
!
!
!
implicit none
integer :: lvl,loc_zones,ier,large_zones,max_zones,i,req,j,num
integer :: tag,counter_thingy
integer,allocatable :: rec_array(:),disp_array(:)

! Get the number of zones for each processor
loc_zones = count_big_zones(lvl)
!! Sum those zones and make sure every processor has a copy so all the other
!! objects can be properly allocated
!
!if (loc_zones == 0) return
!
CALL MPI_Allreduce(loc_zones,ggm_save_info(lvl,num)%big_zones,1,MPI_INTEGER, &
  MPI_SUM,commune,ier)
!!
!!
!!
!ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zones_per_level(lvl) = &
!     ggm_save_info(lvl,num)%big_zones


!
! Distribute the number of zones owned by a processor at a given level.
! Then agglomerate and broadcast so all procs know everything about everyone
! else's zones.
!
if (allocated(ggm_save_info(lvl,num)%little_zones)) deallocate(ggm_save_info(lvl,num)%little_zones)
allocate(ggm_save_info(lvl,num)%little_zones(0:nprocs-1))
!
  ggm_save_info(lvl,num)%little_zones(self) = loc_zones
!  write(*,*) 'local zones',self,'lvl',lvl,'n little zone',loc_zones,ggm_save_info(lvl,num)%little_zones
  if (allocated(rec_array)) deallocate(rec_array)
  allocate(rec_array(0:nprocs-1))
  if (allocated(disp_array)) deallocate(disp_array)
  allocate(disp_array(0:nprocs-1))
  rec_array = 1
  do i = 0,nprocs-1
    disp_array(i) = i
  end do
!  write(*,*) 'receving array',self,rec_array
!  write(*,*) 'displacement array',self,disp_array
  call mpi_allgatherv(ggm_save_info(lvl,num)%little_zones(self),1,&
    MPI_INT,ggm_save_info(lvl,num)%little_zones,&
    rec_array,disp_array,MPI_INT,commune,ier)
!
!
!
  call mpi_barrier(commune,ier)
  deallocate(rec_array)
  deallocate(disp_array)
!  write(*,*) 'Do you hear the people sing?'
!  write(*,*) 'things and stuff',self,lvl,ggm_save_info(lvl,num)%little_zones


!
! Loop to assign start and end values to the zones each processor owns
! XXXX account for chance of no zones at a given level
!  if (ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%regrid_between_outputs) then

    if (allocated(ggm_save_info(lvl,num)%zone_start)) deallocate(ggm_save_info(lvl,num)%zone_start)
    allocate(ggm_save_info(lvl,num)%zone_start(0:nprocs-1))

    if (allocated(ggm_save_info(lvl,num)%zone_end)) deallocate(ggm_save_info(lvl,num)%zone_end)
    allocate(ggm_save_info(lvl,num)%zone_end(0:nprocs-1))

    do i = 0,nprocs-1
      if (i==0) then
        ggm_save_info(lvl,num)%zone_start(i) = ggm_base_info(num)%adapted_zone_start_number
        ggm_save_info(lvl,num)%zone_end(i) = ggm_base_info(num)%adapted_zone_start_number + &
          ggm_save_info(lvl,num)%little_zones(i) - 1
      else
        ggm_save_info(lvl,num)%zone_start(i) = ggm_save_info(lvl,num)%zone_end(i-1)+1
        ggm_save_info(lvl,num)%zone_end(i) = ggm_save_info(lvl,num)%zone_start(i)+&
          ggm_save_info(lvl,num)%little_zones(i)-1
      end if
    end do

    counter_thingy = ggm_save_info(lvl,num)%zone_start(0)
!    do i = 0,nprocs-1
      do j =  1,ggm_save_info(lvl,num)%big_zones
      !ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(j) =
        ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_number(j) = counter_thingy
        counter_thingy = counter_thingy + 1
      end do
!    end do

!    ggm_save_info(lvl,num)%zone_start(i) = base_info(num)%adapted_zone_start_number = &
!         ggm_save_info(lvl,num)%zone_start(nprocs-1) + ggm_save_info(lvl,num)%little_zones(nprocs-1)-1

    if (allocated(ggm_save_info(lvl,num)%lower)) deallocate(ggm_save_info(lvl,num)%lower)
    if (allocated(ggm_save_info(lvl,num)%higher)) deallocate(ggm_save_info(lvl,num)%higher)
    if (allocated(ggm_save_info(lvl,num)%lo_corner)) deallocate(ggm_save_info(lvl,num)%lo_corner)
    if (allocated(ggm_save_info(lvl,num)%hi_corner)) deallocate(ggm_save_info(lvl,num)%hi_corner)
    if (allocated(ggm_save_info(lvl,num)%output_number)) deallocate(ggm_save_info(lvl,num)%output_number)

    allocate(ggm_save_info(lvl,num)%lower(3,ggm_save_info(lvl,num)%zone_start(0):&
      ggm_save_info(lvl,num)%zone_end(nprocs-1)))
    allocate(ggm_save_info(lvl,num)%higher(3,ggm_save_info(lvl,num)%zone_start(0):&
      ggm_save_info(lvl,num)%zone_end(nprocs-1)))
    allocate(ggm_save_info(lvl,num)%lo_corner(3,ggm_save_info(lvl,num)%zone_start(0):&
      ggm_save_info(lvl,num)%zone_end(nprocs-1)))
    allocate(ggm_save_info(lvl,num)%hi_corner(3,ggm_save_info(lvl,num)%zone_start(0):&
      ggm_save_info(lvl,num)%zone_end(nprocs-1)))
    allocate(ggm_save_info(lvl,num)%output_number(ggm_steps(num)%n_outsteps,1))
!
    call ggm_amr_grid_ranges(lvl,num)
!
    call ggm_amr_grid_distro(lvl,num)

!  end if

end subroutine
!
!
!
!
!
!
subroutine ggm_amr_grid_distro(lvl,num)
!
!
!
!
implicit none
integer :: lvl,i,ier,num
integer :: n_elems,counter
!  integer :: stat(MPI_STATUS_SIZE)
integer,allocatable :: disp_array(:),rec_array(:)

!write(*,*) 'entering zone information distribution function.',self

counter = 0

n_elems = 3*ggm_save_info(lvl,num)%big_zones

if (allocated(disp_array)) deallocate(disp_array)
allocate(disp_array(0:nprocs-1))
if (allocated(rec_array)) deallocate(rec_array)
allocate(rec_array(0:nprocs-1))
do i = 0,nprocs-1
  rec_array(i) = ggm_save_info(lvl,num)%little_zones(i)*3
end do
do i=0,nprocs-1
  if (i==0) then
    disp_array(i) = 0
  else
    disp_array(i) = disp_array(i-1)+ggm_save_info(lvl,num)%little_zones(i-1)*3
  end if
end do
!  write(*,*) 'receiving array',self,lvl,rec_array
!  write(*,*) 'displacement array',self,lvl,disp_array
!
! The allgatherv distribute all necessary information for zone writing to every processor
!
!call mpi_allgatherv(ggm_save_info(lvl,num)%lower(1:3,ggm_save_info(lvl,num)%zone_start(self):&
!  ggm_save_info(lvl,num)%zone_end(self)),ggm_save_info(lvl,num)%little_zones(self),MPI_INT,&
!  ggm_save_info(lvl,num)%lower,rec_array,disp_array,&
!  MPI_INT,commune,ier)
!!
!call mpi_allgatherv(ggm_save_info(lvl,num)%higher(1:3,ggm_save_info(lvl,num)%zone_start(self):&
!  ggm_save_info(lvl,num)%zone_end(self)),ggm_save_info(lvl,num)%little_zones(self),MPI_INT,&
!  ggm_save_info(lvl,num)%higher,rec_array,disp_array,&
!  MPI_INT,commune,ier)
!!
!call mpi_allgatherv(ggm_save_info(lvl,num)%lo_corner(1:3,ggm_save_info(lvl,num)%zone_start(self):&
!  ggm_save_info(lvl,num)%zone_end(self)),ggm_save_info(lvl,num)%little_zones(self),MPI_DOUBLE_PRECISION,&
!  ggm_save_info(lvl,num)%lo_corner,rec_array,disp_array,&
!  MPI_DOUBLE_PRECISION,commune,ier)
!!
!call mpi_allgatherv(ggm_save_info(lvl,num)%hi_corner(1:3,ggm_save_info(lvl,num)%zone_start(self):&
!  ggm_save_info(lvl,num)%zone_end(self)),ggm_save_info(lvl,num)%little_zones(self),MPI_DOUBLE_PRECISION,&
!  ggm_save_info(lvl,num)%hi_corner,rec_array,disp_array,&
!  MPI_DOUBLE_PRECISION,commune,ier)
!
!
! XXXX account for chance of no zones at a given level
!
call mpi_allgatherv(ggm_save_info(lvl,num)%lower(1:3,ggm_save_info(lvl,num)%zone_start(self):&
  ggm_save_info(lvl,num)%zone_end(self)),rec_array(self),MPI_INT,&
  ggm_save_info(lvl,num)%lower,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(ggm_save_info(lvl,num)%higher(1:3,ggm_save_info(lvl,num)%zone_start(self):&
  ggm_save_info(lvl,num)%zone_end(self)),rec_array(self),MPI_INT,&
  ggm_save_info(lvl,num)%higher,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(ggm_save_info(lvl,num)%lo_corner(1:3,ggm_save_info(lvl,num)%zone_start(self):&
  ggm_save_info(lvl,num)%zone_end(self)),rec_array(self),MPI_DOUBLE_PRECISION,&
  ggm_save_info(lvl,num)%lo_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)
!
call mpi_allgatherv(ggm_save_info(lvl,num)%hi_corner(1:3,ggm_save_info(lvl,num)%zone_start(self):&
  ggm_save_info(lvl,num)%zone_end(self)),rec_array(self),MPI_DOUBLE_PRECISION,&
  ggm_save_info(lvl,num)%hi_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)

call mpi_barrier(commune,ier)

deallocate(rec_array)
deallocate(disp_array)


end subroutine
!
!
!
!
!
subroutine ggm_amr_grid_ranges(lvl,num)


integer :: counter,lvl,num!_zones,num,zone_starts(0:nprocs-1),zone_ends(0:nprocs-1)
!,ind,zone_starts(0:nprocs-1),zone_ends(0:nprocs-1)
!integer :: lo_nodes(3,num_zones),hi_nodes(3,num_zones)
!real(kind=dp) :: min_corner(3,num_zones),max_corner(3,num_zones)
real(kind=dp) :: poco_loco(3),poco_hoco(3)
!  logical :: is_full(num_zones)

type(amrex_mfiter) :: mfi
type(amrex_box) :: brx

!write(*,*) 'entering zone ranges sub',self,lvl

counter = ggm_save_info(lvl,num)%zone_start(self)
call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
do while (mfi%next())

  brx = mfi%validbox()

!  lo_nodes(1,counter) = brx%lo(1)
!  lo_nodes(2,counter) = brx%lo(2)
!  lo_nodes(3,counter) = brx%lo(3)
!  hi_nodes(1,counter) = brx%hi(1)
!  hi_nodes(2,counter) = brx%hi(2)
!  hi_nodes(3,counter) = brx%hi(3)
  ggm_save_info(lvl,num)%lower(1:3,counter) = brx%lo
  ggm_save_info(lvl,num)%higher(1:3,counter) = brx%hi
!    write(*,*) 'lo nodes',self,lvl,counter,lo_nodes(1,counter),lo_nodes(2,counter),lo_nodes(3,counter)
!    write(*,*) 'hi nodes',self,lvl,counter,hi_nodes(1,counter),hi_nodes(2,counter),hi_nodes(3,counter)
  poco_loco = get_real_coords(amrex_geom(lvl),brx%lo)
  poco_hoco = get_real_coords(amrex_geom(lvl),brx%hi)

!  min_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%lo)
!  max_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%hi)
  ggm_save_info(lvl,num)%lo_corner(1:3,counter) = poco_loco
  ggm_save_info(lvl,num)%hi_corner(1:3,counter) = poco_hoco
!  min_corner(1,counter) = poco_loco(1)
!  min_corner(2,counter) = poco_loco(2)
!  min_corner(3,counter) = poco_loco(3)
!  max_corner(1,counter) = poco_hoco(1)
!  max_corner(2,counter) = poco_hoco(2)
!  max_corner(3,counter) = poco_hoco(3)
!  write(*,*) self,lvl,min_corner(1:3,counter),max_corner(1:3,counter)
  counter = counter +1

end do

call amrex_mfiter_destroy(mfi)

end subroutine
!
!
!
!
!
subroutine ggm_amr_grid_zone_writer(lvl,num,index_file)
!
!
!
integer :: fine,i,n,num,zone_counter,index_file,index_base,ier
integer :: lvl,index_zone
integer(cgsize_t) :: isize(3,3)
!integer,intent(in),optional :: ggm_offset
character(len=32) :: zone_number,zonename,format_nm_lvl_and_zone

!if (present(ggm_offset)) then
!  index_file = sv_num+ggm_offset
!else
!index_file = num+1
!end if
index_base = 1
!
! Set up the zones for output
!
format_nm_lvl_and_zone = '(I2.2,a1,I3.3)'
!
!
!
zone_counter = ggm_base_info(num)%adapted_zone_start_number
!do lvl = min_lvl,max_lvl
!
!  call zone_box_setup(lvl)
! XXXX account for chance of no zones at a given level
  do n = 0,nprocs-1
    do i=ggm_save_info(lvl,num)%zone_start(n),ggm_save_info(lvl,num)%zone_end(n)
!    zone_counter = zone_counter+1
      ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_number(i) = i
      ggm_save_info(lvl,num)%multigrid_save(ggm_steps(num)%output_counter)%zone_lvls(i) = lvl
!    write(*,*) 'Sing with me',self,zone_counter,i
      isize(1,1) = ggm_save_info(lvl,num)%higher(1,i)-ggm_save_info(lvl,num)%lower(1,i)+1
      isize(2,1) = ggm_save_info(lvl,num)%higher(2,i)-ggm_save_info(lvl,num)%lower(2,i)+1
      isize(3,1) = ggm_save_info(lvl,num)%higher(3,i)-ggm_save_info(lvl,num)%lower(3,i)+1
      isize(1,2) = isize(1,1)-1
      isize(2,2) = isize(2,1)-1
      isize(3,2) = isize(3,1)-1
      isize(1,3) = 0
      isize(2,3) = 0
      isize(3,3) = 0
!    write(*,*) 'isize check',self,isize
!    write(*,*) 'pre-flight check',self,zone_counter,i,ggm_save_info(lvl,num)%big_zones
      write(zone_number,format_nm_lvl_and_zone) lvl,'_',i
      zonename = 'zone_lvl_'//TRIM(zone_number)
!    write(*,*) 'writing zone,',zonename,index_file,index_base,i,zone_counter

      ggm_save_info(lvl,num)%output_number(i,1) = i
    index_zone = i
!    write(*,*) 'zone output number',self,ggm_save_info(lvl,num)%output_number(i,sv_num),zone_counter
!    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
!      Structured,ggm_save_info(lvl,num)%output_number(i,sv_num),ier)
      CALL cg_zone_write_f(index_file,index_base,zonename,isize,Structured,index_zone,ier)
!  write(*,*) 'after writing zone',index_file,index_base,i,ggm_save_info(lvl,num)%output_number(i,sv_num)
!    write(*,*) 'after zone write, ggm ',index_file,index_base,i,self
!   index_coord = index_zone
      call ggm_amr_grid_coord_setup(index_file,index_base,i)
!  write(*,*) index_file,index_base,zone_counter
!    call mpi_barrier(commune,ier)
    end do
  end do

!ggm_base_info(num)%adapted_zone_start_number = ggm_base_info(num)%adapted_zone_start_number + &
!  ggm_save_info(lvl,num)%big_zones

!end do

end subroutine
!
!
!
!
!
subroutine ggm_amr_grid_coord_setup(index_file,index_base,index_zone)
!
!
!
implicit none
integer :: coord_x,coord_y,coord_z,ier!,lvl,num
integer :: index_file,index_base,index_zone
!integer,intent(in),optional :: ggm_offset

coord_x = 1
coord_y = 2
coord_z = 3
!
CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
  'CoordinateX',coord_x,ier)
!
IF (ier /= CG_OK) CALL cgp_error_exit_f
!
CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
  'CoordinateY',coord_y,ier)
IF (ier /= CG_OK) CALL cgp_error_exit_f
!
CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
  'CoordinateZ',coord_z,ier)
!
IF (ier /= CG_OK) CALL cgp_error_exit_f

end subroutine
!
!
!
!
!
subroutine ggm_amr_grid_coord_data_write(lvl,num,index_file)
!
!
!
implicit none
integer ::a,b,c,lvl,out_min,out_max,lvl_hi,lvl_lo,i,ier,num
integer :: index_file,index_base,index_zone,sv_num
integer :: coord_x,coord_y,coord_z
integer(cgsize_t) ::isize(3,3),loc_hi(3),loc_lo(3)
!integer,intent(in),optional :: ggm_offset

coord_x = 1
coord_y = 2
coord_z = 3
!index_file = num + 1
index_base = 1
!if (present(ggm_offset)) then
!  if (index_file < 1001) then
!    index_file = 1000 + index_file - 1
!  end if
!end if
!write(*,*) 'inside zone coordinate data write subroutine',self,lvl_lo,lvl_hi
!write(*,*) 'amrex limits',amrex_problo,amrex_probhi

!do lvl = lvl_lo,lvl_hi

  !write(*,*) 'amrex dx',amrex_geom(lvl)%dx(1),amrex_geom(lvl)%dx(2),amrex_geom(lvl)%dx(3)

  do i = ggm_save_info(lvl,num)%zone_start(self),ggm_save_info(lvl,num)%zone_end(self)

    isize(1,1) = ggm_save_info(lvl,num)%higher(1,i)-ggm_save_info(lvl,num)%lower(1,i)+1
    isize(2,1) = ggm_save_info(lvl,num)%higher(2,i)-ggm_save_info(lvl,num)%lower(2,i)+1
    isize(3,1) = ggm_save_info(lvl,num)%higher(3,i)-ggm_save_info(lvl,num)%lower(3,i)+1

    isize(1,2) = isize(1,1)-1
    isize(2,2) = isize(2,1)-1
    isize(3,2) = isize(3,1)-1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

!  write(*,*) 'isize array',isize
!    write(*,*) 'allocating vertices'
    allocate(x_verts(isize(1,1),isize(2,1),isize(3,1)))
    allocate(y_verts(isize(1,1),isize(2,1),isize(3,1)))
    allocate(z_verts(isize(1,1),isize(2,1),isize(3,1)))

!    WRITE(*,*) 'Hi Im processor ', self,' and my grid range is ',byx_lo(1),byx_lo(2),&
!      byx_lo(3),byx_hi(1),byx_hi(2),byx_hi(3),index_zone
    do c = 1,isize(3,1)
      do b = 1,isize(2,1)
        do a = 1,isize(1,1)
          x_verts(a,b,c) = ggm_save_info(lvl,num)%lo_corner(1,i) + (a-1)*rdx(lvl)
          y_verts(a,b,c) = ggm_save_info(lvl,num)%lo_corner(2,i) + (b-1)*rdx(lvl)
          z_verts(a,b,c) = ggm_save_info(lvl,num)%lo_corner(3,i) + (c-1)*rdx(lvl)
        end do
      end do
    end do

    loc_lo = (/1,1,1/)
    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)

!    index_zone = ggm_save_info(lvl,num)%output_number(i,1)
!    write(*,*) 'repeat index info',index_file,index_base,index_zone
!    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!
!    write(*,*) 'box corners',ggm_save_info(lvl,num)%lo_corner(1:3,i),x_verts(isize(1,1),isize(2,1),isize(3,1)),&
!      y_verts(isize(1,1),isize(2,1),isize(3,1)),z_verts(isize(1,1),isize(2,1),isize(3,1)),lvl
    CALL cgp_coord_write_data_f(index_file,index_base,i,coord_x,loc_lo,loc_hi,x_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'reading?'
!    write(*,*) 'x coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
    !CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,byx_lo,byx_hi,y_verts,ier)
    CALL cgp_coord_write_data_f(index_file,index_base,i,coord_y,loc_lo,loc_hi,y_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'y coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
    CALL cgp_coord_write_data_f(index_file,index_base,i,coord_z,loc_lo,loc_hi, z_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'z coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    WRITE(*,*) 'Coordinate value write completed.'

    deallocate(x_verts)
    deallocate(y_verts)
    deallocate(z_verts)

  end do
end subroutine
!
!
!
!
!

!
!
!
!
!
subroutine ggm_amr_grid_output_nums(num,min_lvl,max_lvl)
!
!
!
implicit none
integer,intent(in) :: num,min_lvl,max_lvl
integer :: lvl,num_lower_lvls,zone,i,ier
integer,allocatable :: rec_array(:),disp_array(:)
!
!
allocate(disp_array(0:nprocs-1))
allocate(rec_array(0:nprocs-1))
!

!
!
!
num_lower_lvls = 0

do lvl = min_lvl,max_lvl
  do zone = ggm_save_info(lvl,num)%zone_start(self),ggm_save_info(lvl,num)%zone_end(self)
    ggm_save_info(lvl,num)%output_number(zone,1) = zone + num_lower_lvls
!    write(*,*) 'output numbers',ggm_save_info(lvl,num)%zone_start(self),&
!      ggm_save_info(lvl,num)%zone_end(self),self,ggm_save_info(lvl,num)%output_number(zone,num)

  end do

  num_lower_lvls = num_lower_lvls + ggm_save_info(lvl,num)%big_zones
!
! Set up the receiving arrays (how much info for each processor) and the displacement array
!   (where that info goes)
!
  do i=0,nprocs-1

    rec_array(i) = ggm_save_info(lvl,num)%little_zones(i)
    if (i==0) then
      disp_array(i) = 0
    else
      disp_array(i) = disp_array(i-1)+ggm_save_info(lvl,num)%little_zones(i-1)
    end if

  end do

  !write(*,*) 'local output numbers',ggm_save_info(lvl,num)%output_number,self

  call mpi_allgatherv(ggm_save_info(lvl,num)%output_number(ggm_save_info(lvl,num)%zone_start(self):&
    ggm_save_info(lvl,num)%zone_end(self),1),ggm_save_info(lvl,num)%little_zones(self),MPI_INTEGER,&
    ggm_save_info(lvl,num)%output_number,rec_array,disp_array,&
    MPI_INTEGER,commune,ier)
!  write(*,*) 'local output numbers',ggm_save_info(lvl,num)%output_number,&
!    ggm_save_info(lvl,num)%little_zones,self
!  write(*,*) 'mpi junk',rec_array,disp_array,self
!
!
!
end do

deallocate(disp_array)
deallocate(rec_array)



end subroutine
!
!
!
!
!
subroutine ggm_amr_save_writer(lvl,num,fine,index_file,index_base)
!
!
!
!
!
!
use precise
use amrex_amr_module
use amrex_base_module
use zone_writing
use output_data
use nml_output
use grid_data
use linkwise
use ggm_stuff
use timing
use constants
use freestream_values
use cgns
use mpi
use amr_info_holder
use amr_processes, only: self,nprocs!,commune
implicit none
INTEGER :: index_base,index_zone,index_file,iphysdim,ier,index_field,dummy
integer :: fine,max_lvl,lvl,num,zee_count,total_zones,zer,n,index_sol
integer,allocatable :: level_zones(:)
integer*8 :: isize(3,3),loc_lo(4),loc_hi(4)
CHARACTER(len=32) :: zonename,solut_name,sol_num,ggm_field_value

type(amrex_box) :: bux
type(amrex_mfiter) :: mfi

if (fine > ggm_max_lvl(num)) then
  max_lvl = ggm_max_lvl(num)
else
  max_lvl = fine
end if

!do index_zone = 1,ggm_base_info(num)%zones_per_lvl(lvl,ggm_steps(num)%output_counter)
do n = 0,nprocs-1
  do index_zone = ggm_save_info(lvl,num)%zone_start(n),ggm_save_info(lvl,num)%zone_end(n)
    call cg_sol_write_f(index_file,index_base,index_zone,&
        ggm_steps(num)%sol_names(ggm_steps(num)%output_counter),Vertex,&
        index_sol,ier)
!    call cg_error_print_f
    write(*,*) 'after ggm solution',self,lvl,index_file,index_base,index_zone,&
        ggm_steps(num)%sol_names(save_steps(num)%output_counter),&
        ggm_steps(num)%output_counter
!
    IF (ier /= CG_OK) CALL cgp_error_exit_f
  end do
end do
  !    call cg_error_print_f
  !    write(*,*) 'after solution',self,index_file,index_base,&
  !      zone_storage(lvl)%output_number(index_zone,1),index_flow
  select case (ggm_variables(num))
    case (1)
      ggm_field_value = 'GGM Density'
    case (2)
      ggm_field_value = 'GGM u-Velocity'
    case (3)
      ggm_field_value = 'GGM v-Velocity'
    case (4)
      ggm_field_value = 'GGM w-Velocity'
    case (5)
      ggm_field_value = 'GGM Temperature'
    case (6)
      ggm_field_value = 'GGM Pressure'
    case (7)
      ggm_field_value = 'GGM Velocity Magnitude'
    case (8)
      ggm_field_value = 'GGM Dynamic Pressure'
    case (9)
      ggm_field_value = 'GGM Mach Number'
    case (10)
      ggm_field_value = 'GGM xy-Vorticity'
    case (11)
      ggm_field_value = 'GGM xz-Vorticity'
    case (12)
      ggm_field_value = 'GGM yz-Vorticity'
    case (13)
      ggm_field_value = 'GGM Vorticity Magnitude'
    case (14)
      ggm_field_value = 'GGM xy-Shear Stress'
    case (15)
      ggm_field_value = 'GGM xz-Shear Stress'
    case (16)
      ggm_field_value = 'GGM yz-Shear Stress'
    case (17)
      ggm_field_value = 'GGM Shear Stress Magnitude'
    case (20)
      ggm_field_value = 'GGM Alpha'
    case (29)
      ggm_field_value = 'GGM Q-criterion'

  end select
!
!
!
!
!do index_zone = 1,base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)

  index_field = 1
do n = 0,nprocs-1
  do index_zone = ggm_save_info(lvl,num)%zone_start(n),ggm_save_info(lvl,num)%zone_end(n)
    call cgp_field_write_f(index_file,index_base,index_zone,&
      index_sol,RealDouble,ggm_field_value,index_field,ier)

  end do
end do
!
!
!
!
!
index_zone = ggm_save_info(lvl,num)%zone_start(self)
!
call amrex_mfiter_build(mfi,mfggm(lvl),tiling=.false.)
do while (mfi%next())
  bux = mfi%validbox()
!    if (grid_info(lvl,num)%full(b)) then
!        write(*,*) self,grid_info(lvl,num)%full,c
!        write(*,*) self,grid_info(lvl,num)%zones_owned
!        write(*,*) 'bounds of rho',lbound(rho),ubound(rho)
!        write(*,*) 'bounds of the box',bux%lo,bux%hi
!        write(*,*) 'rho at low',rho(bux%lo(1),bux%lo(2),bux%lo(3),1)
!        write(*,*) 'rho at hi',rho(bux%hi(1),bux%hi(2),bux%hi(3),1)
  loc_lo = 1
  loc_hi(1) = bux%hi(1)-bux%lo(1)+1
  loc_hi(2) = bux%hi(2)-bux%lo(2)+1
  loc_hi(3) = bux%hi(3)-bux%lo(3)+1
  loc_hi(4) = 1

  ggm => mfggm(lvl)%dataptr(mfi)

  call cgp_field_write_data_f(index_file,index_base,index_zone,&
    index_sol,index_field,loc_lo,loc_hi,ggm(bux%lo(1):bux%hi(1),&
    bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
  if (ier /= CG_OK) call cgp_error_exit_f

  index_zone = index_zone+1
end do

  ggm_steps(num)%write_done(ggm_steps(num)%output_counter) = .true.
  ggm_steps(num)%output_counter = ggm_steps(num)%output_counter + 1

  write(11,*) 'GGM save successfully written.'
  write(*,*) 'GGM save successfully written.'

end subroutine
!
!
!
!
!

!subroutine zone_regrid_checks(num)
!use precise
!use constants
!use output_data, only: save_steps
!use timing
!use ggm_stuff, only: num_regrids
!use amr_info_holder, only: fixed_lvl
!implicit none
!integer :: i,j
!integer,intent(in) :: num
!
!do i = 1,num_regrids
!  do j = 2,save_steps(num)
!!
!! Check if a regrid occurs between timesteps
!!
!    if (regrid_time(i) >= save_steps(num)%saved_times(j-1) .and. &
!        regrid_time(i) <= save_steps(num)%saved_times(j)) then
!
!      do lvl = 0,amrex_max_level
!        if (lvl <= fixed_level) then
!          save_steps(num)%multigrid_save(i)%regrid_between_outputs(lvl) = .false.
!        else
!          save_steps(num)%multigrid_save(i)%regrid_between_outputs(lvl) = .true.
!        end if
!      end do
!      exit
!
!    else
!      save_steps(num)%multigrid_save(i)%regrid_between_outputs(0:amrex_max_level) = .false.
!    end if
!
!  end do
!end do
!
!end subroutine
!!
!!
!!
!!
!!
!subroutine zone_regrid_checks_ggm(num)
!use precise
!use constants
!use output_data, only: ggm_save_info
!use timing
!use ggm_stuff, only: num_regrids
!use amr_info_holder, only: fixed_lvl
!implicit none
!integer :: i,j,lvl
!integer,intent(in) :: num
!
!do i = 1,num_regrids
!  do j = 2,save_steps(num)
!!
!! Check if a regrid occurs between timesteps
!!
!    if (regrid_time(i) >= ggm_save_info(num)%saved_times(j-1) .and. &
!        regrid_time(i) <= ggm_save_info(num)%saved_times(j)) then
!
!      do lvl = 0,amrex_max_level
!        if (lvl <= fixed_level) then
!          ggm_save_info(num)%multigrid_save(i)%regrid_between_outputs(lvl) = .false.
!        else
!          ggm_save_info(num)%multigrid_save(i)%regrid_between_outputs(lvl) = .true.
!        end if
!      end do
!      exit
!
!    else
!      ggm_save_info(num)%multigrid_save(i)%regrid_between_outputs(0:amrex_max_level) = .false.
!    end if
!
!  end do
!end do
!
!end subroutine
!
end module



!  write(*,*) 'zone start and end',self,lvl,ggm_save_info(lvl,num)%zone_start,ggm_save_info(lvl,num)%zone_end
!  write(*,*) 'Singing the songs of angry men!'
!
! allreduce sums the local zones for each processor to find the total number of zones on a level
!


!  if (allocated(ggm_save_info(lvl,num)%output_number)) deallocate(ggm_save_info(lvl,num)%output_number)
!
!  max_zones = ggm_save_info(lvl,num)%big_zones/nprocs+1
!
!  write(*,*) "I'm",self,'number of big zones',ggm_save_info(lvl,num)%big_zones,&
!    'at level ',lvl

!  allocate(ggm_save_info(lvl,num)%output_number(ggm_save_info(lvl,num)%big_zones,num_saves))
!
!
!
!subroutine amr_grid_coord_data_write(index_file,index_base,lvl_lo,lvl_hi,sv_num,ggm_offset)
!!
!!
!!
!!
!implicit none
!integer ::a,b,c,lvl,out_min,out_max,lvl_hi,lvl_lo,i,ier
!integer :: index_file,index_base,index_zone,sv_num
!integer :: coord_x,coord_y,coord_z
!integer(cgsize_t) ::isize(3,3),loc_hi(3),loc_lo(3)
!integer,intent(in),optional :: ggm_offset
!
!coord_x = 1
!coord_y = 2
!coord_z = 3
!if (present(ggm_offset)) then
!  if (index_file < 1001) then
!    index_file = 1000 + index_file - 1
!  end if
!end if
!!write(*,*) 'inside zone coordinate data write subroutine',self,lvl_lo,lvl_hi
!!write(*,*) 'amrex limits',amrex_problo,amrex_probhi
!
!do lvl = lvl_lo,lvl_hi
!
!  !write(*,*) 'amrex dx',amrex_geom(lvl)%dx(1),amrex_geom(lvl)%dx(2),amrex_geom(lvl)%dx(3)
!
!  do i = ggm_save_info(lvl,num)%zone_start(self),ggm_save_info(lvl,num)%zone_end(self)
!
!    isize(1,1) = ggm_save_info(lvl,num)%higher(1,i)-ggm_save_info(lvl,num)%lower(1,i)+1
!    isize(2,1) = ggm_save_info(lvl,num)%higher(2,i)-ggm_save_info(lvl,num)%lower(2,i)+1
!    isize(3,1) = ggm_save_info(lvl,num)%higher(3,i)-ggm_save_info(lvl,num)%lower(3,i)+1
!
!    isize(1,2) = isize(1,1)-1
!    isize(2,2) = isize(2,1)-1
!    isize(3,2) = isize(3,1)-1
!
!    isize(1,3) = 0
!    isize(2,3) = 0
!    isize(3,3) = 0
!
!!  write(*,*) 'isize array',isize
!!    write(*,*) 'allocating vertices'
!    allocate(x_verts(isize(1,1),isize(2,1),isize(3,1)))
!    allocate(y_verts(isize(1,1),isize(2,1),isize(3,1)))
!    allocate(z_verts(isize(1,1),isize(2,1),isize(3,1)))
!
!!    WRITE(*,*) 'Hi Im processor ', self,' and my grid range is ',byx_lo(1),byx_lo(2),&
!!      byx_lo(3),byx_hi(1),byx_hi(2),byx_hi(3),index_zone
!    do c = 1,isize(3,1)
!      do b = 1,isize(2,1)
!        do a = 1,isize(1,1)
!          x_verts(a,b,c) = ggm_save_info(lvl,num)%lo_corner(1,i) + (a-1)*rdx(lvl)
!          y_verts(a,b,c) = ggm_save_info(lvl,num)%lo_corner(2,i) + (b-1)*rdx(lvl)
!          z_verts(a,b,c) = ggm_save_info(lvl,num)%lo_corner(3,i) + (c-1)*rdx(lvl)
!        end do
!      end do
!    end do
!
!    loc_lo = (/1,1,1/)
!    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)
!
!    index_zone = ggm_save_info(lvl,num)%output_number(i,sv_num)
!!    write(*,*) 'repeat index info',index_file,index_base,index_zone
!!    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!!
!!    write(*,*) 'box corners',ggm_save_info(lvl,num)%lo_corner(1:3,i),x_verts(isize(1,1),isize(2,1),isize(3,1)),&
!!      y_verts(isize(1,1),isize(2,1),isize(3,1)),z_verts(isize(1,1),isize(2,1),isize(3,1)),lvl
!    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_x,&
!      loc_lo,loc_hi,x_verts,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!    write(*,*) 'reading?'
!!    write(*,*) 'x coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!!    index_coord = index_coord+1
!    !CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,byx_lo,byx_hi,y_verts,ier)
!    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,&
!      loc_lo,loc_hi,y_verts,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!    write(*,*) 'y coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!!    index_coord = index_coord+1
!    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_z,&
!      loc_lo,loc_hi, z_verts,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!    write(*,*) 'z coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!!    WRITE(*,*) 'Coordinate value write completed.'
!
!    deallocate(x_verts)
!    deallocate(y_verts)
!    deallocate(z_verts)
!
!  end do
!end do
!
!
!end subroutine


!do lvl = save_shapes(1,num),max_lvl
!  do index_zone = 1,base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)
!!
!!    call cg_error_print_f
!    call cg_sol_write_f(index_file,index_base,ggm_save_info(lvl,num)%&
!      multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!      save_steps(num)%sol_names(save_steps(num)%output_counter),Vertex,&
!      save_steps(num)%output_counter,ier)
!
!!    call cg_error_print_f
!!    write(*,*) 'after solution',self,index_file,index_base,&
!!      ggm_save_info(lvl,num)%output_number(index_zone,1),index_flow
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!  end do
!!end do
!zee_count = 0
!!do lvl = save_shapes(1,num),max_lvl
!
!  do index_zone = 1,base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)
!
!    index_field = 1
!!    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!!      ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(1,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'Density',index_field,ier)
!      index_field = index_field+1
!    end if
!!    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!!      ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(2,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'U-velocity',index_field,ier)
!      index_field = index_field+1
!    end if
!!    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!!      ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(3,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'V-velocity',index_field,ier)
!      index_field = index_field+1
!    end if
!!    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!!      ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(4,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'W-velocity',index_field,ier)
!      index_field = index_field+1
!    end if
!!    write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!!      ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(5,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'Pressure',index_field,ier)
!      index_field = index_field+1
!    end if
!    !write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!    !  ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(6,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'Temperature',index_field,ier)
!      index_field = index_field+1
!    end if
!    !write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
!    !  ggm_save_info(lvl,num)%output_number(index_zone,num),self
!    if (save_data(7,num)) then
!      call cgp_field_write_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,RealDouble,'Alpha',index_field,ier)
!      index_field = index_field+1
!    end if
!
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!  end do
!!end do
!
!!do lvl = save_shapes(1,num),max_lvl
!index_zone = ggm_save_info(lvl,num)%zone_start(self)
!  call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
!  do while (mfi%next())
!    bux = mfi%validbox()
!!    if (ggm_save_info(lvl,num)%full(b)) then
!!        write(*,*) self,ggm_save_info(lvl,num)%full,c
!!        write(*,*) self,ggm_save_info(lvl,num)%zones_owned
!!        write(*,*) 'bounds of rho',lbound(rho),ubound(rho)
!!        write(*,*) 'bounds of the box',bux%lo,bux%hi
!!        write(*,*) 'rho at low',rho(bux%lo(1),bux%lo(2),bux%lo(3),1)
!!        write(*,*) 'rho at hi',rho(bux%hi(1),bux%hi(2),bux%hi(3),1)
!    loc_lo = 1
!    loc_hi(1) = bux%hi(1)-bux%lo(1)+1
!    loc_hi(2) = bux%hi(2)-bux%lo(2)+1
!    loc_hi(3) = bux%hi(3)-bux%lo(3)+1
!    loc_hi(4) = 1
!!
!!
!!
!    index_field = 1
!!        call cgp_field_write_data_f(index_file,index_base,&
!!          ggm_save_info(lvl,num)%zones_owned(c),&
!!          index_flow,index_field,loc_lo,loc_hi,rho(bux%lo(1):bux%hi(1),&
!!          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!    if (save_data(1,num)) then
!      rho => mfrho(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, rho',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,rho(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!    if (save_data(2,num)) then
!      u_vel => mfu_vel(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, u',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,u_vel(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!    if (save_data(3,num)) then
!      v_vel => mfv_vel(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, v',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,v_vel(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!    if (save_data(4,num)) then
!      w_vel => mfw_vel(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, w',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,w_vel(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!    if (save_data(5,num)) then
!      press => mfpress(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, pressure',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,press(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!    if (save_data(6,num)) then
!      temp => mftemp(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, temp',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,temp(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!
!    if (save_data(7,num)) then
!      alpha => mfalpha(lvl)%dataptr(mfi)
!      !write(*,*) 'field data write, temp',index_file,index_base,save_steps(num)%output_counter,index_field,&
!      !ggm_save_info(lvl,num)%output_number(zee_count,num),self
!      call cgp_field_write_data_f(index_file,index_base,&
!        ggm_save_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(index_zone),&
!        save_steps(num)%output_counter,index_field,loc_lo,loc_hi,alpha(bux%lo(1):bux%hi(1),&
!        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
!      index_field = index_field+1
!    end if
!!      write(*,*) 'field data',self,index_file,index_base,b,&
!!        index_flow,index_field,loc_lo,loc_hi
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    index_zone = index_zone+1
!  end do
