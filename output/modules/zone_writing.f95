module zone_writing
!
! All the zone writing subroutines so I can have the external interface and use
!   fun modern Fortran features
!
!
!
!
!


implicit none
save
private
public :: zone_writer,zone_coord_data_write,zone_ranges,zone_box_setup,&
  zone_dealloc_safety,zone_distro,zone_coord_setup,zone_output_num_setup
!public :: zone_ketchup,zone_backfill



!
!
!
!
!
contains
!
!
!
subroutine zone_writer(min_lvl,max_lvl,sv_num,ggm_offset)
!
! Zone writing for the time accurate, deforming grid
!
!
! Called by: save_write
! Calls: zone_coord_setup
! External Calls: cg_zone_write_f
!
use precise
use amrex_base_module
use mpi
use nml_output
use cgns
use amr_info_holder , only: zone_storage
use amr_processes, only: self,nprocs,commune
implicit none
integer :: fine,i,sv_num,zone_counter,index_file,index_base,ier,index_coord
integer :: max_lvl,lvl,min_lvl
integer(cgsize_t) :: isize(3,3)
integer,intent(in),optional :: ggm_offset
character(len=32) :: zone_number,zonename,format_nm_lvl_and_zone

if (present(ggm_offset)) then
  index_file = sv_num+ggm_offset
else
  index_file = sv_num+1
end if
index_base = 1
!
! Set up the zones for output
!
format_nm_lvl_and_zone = '(I2.2,a1,I3.3)'
!
!
!
zone_counter = 0
do lvl = min_lvl,max_lvl

  call zone_box_setup(lvl)

  do i=1,zone_storage(lvl)%big_zones
    zone_counter = zone_counter+1
!    write(*,*) 'Sing with me',self,zone_counter,i
    isize(1,1) = zone_storage(lvl)%higher(1,i)-zone_storage(lvl)%lower(1,i)+1
    isize(2,1) = zone_storage(lvl)%higher(2,i)-zone_storage(lvl)%lower(2,i)+1
    isize(3,1) = zone_storage(lvl)%higher(3,i)-zone_storage(lvl)%lower(3,i)+1
    isize(1,2) = isize(1,1)-1
    isize(2,2) = isize(2,1)-1
    isize(3,2) = isize(3,1)-1
    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0
!    write(*,*) 'isize check',self,isize
!    write(*,*) 'pre-flight check',self,zone_counter,i,zone_storage(lvl)%big_zones
    write(zone_number,format_nm_lvl_and_zone) lvl,'_',zone_counter
    zonename = 'zone_lvl_'//TRIM(zone_number)
!    write(*,*) 'writing zone,',zonename,index_file,index_base,i,zone_counter

    zone_storage(lvl)%output_number(i,sv_num) = zone_counter

!    write(*,*) 'zone output number',self,zone_storage(lvl)%output_number(i,sv_num),zone_counter
!    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
!      Structured,zone_storage(lvl)%output_number(i,sv_num),ier)
    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
      Structured,zone_counter,ier)
!  write(*,*) 'after writing zone',index_file,index_base,zone_counter,zone_storage(lvl)%output_number(i,sv_num)
!   index_coord = index_zone
    call zone_coord_setup(index_file,index_base,zone_counter)
!  write(*,*) index_file,index_base,zone_counter
!    call mpi_barrier(commune,ier)
  end do

end do

end subroutine zone_writer
!
!
!
!
!
subroutine zone_ranges(lvl,num_zones,lo_nodes,hi_nodes,min_corner,&
  max_corner,zone_starts,zone_ends)
!
! This gets the lowest corner and higher corner for all the zones to be called in
!
!
!
!
use precise
use mpi
use amr_processes, only : get_real_coords,self,nprocs
use amr_info_holder, only: mfrho
use amrex_base_module
use amrex_amr_module
!use mpi
implicit none
integer :: counter,lvl,num_zones,ind,zone_starts(0:nprocs-1),&
  zone_ends(0:nprocs-1)
integer :: lo_nodes(3,num_zones),hi_nodes(3,num_zones)
real(kind=dp) :: min_corner(3,num_zones),max_corner(3,num_zones)
real(kind=dp) :: poco_loco(3),poco_hoco(3)
!  logical :: is_full(num_zones)

type(amrex_mfiter) :: mfi
type(amrex_box) :: brx

!write(*,*) 'entering zone ranges sub',self,lvl

counter = zone_starts(self)
call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
do while (mfi%next())

  brx = mfi%validbox()

  lo_nodes(1,counter) = brx%lo(1)
  lo_nodes(2,counter) = brx%lo(2)
  lo_nodes(3,counter) = brx%lo(3)
  hi_nodes(1,counter) = brx%hi(1)
  hi_nodes(2,counter) = brx%hi(2)
  hi_nodes(3,counter) = brx%hi(3)
!    write(*,*) 'lo nodes',self,lvl,counter,lo_nodes(1,counter),lo_nodes(2,counter),lo_nodes(3,counter)
!    write(*,*) 'hi nodes',self,lvl,counter,hi_nodes(1,counter),hi_nodes(2,counter),hi_nodes(3,counter)
    poco_loco = get_real_coords(amrex_geom(lvl),brx%lo)
    poco_hoco = get_real_coords(amrex_geom(lvl),brx%hi)

!  min_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%lo)
!  max_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%hi)

    min_corner(1,counter) = poco_loco(1)
    min_corner(2,counter) = poco_loco(2)
    min_corner(3,counter) = poco_loco(3)
    max_corner(1,counter) = poco_hoco(1)
    max_corner(2,counter) = poco_hoco(2)
    max_corner(3,counter) = poco_hoco(3)
!  write(*,*) self,lvl,min_corner(1:3,counter),max_corner(1:3,counter)
  counter = counter +1

end do

call amrex_mfiter_destroy(mfi)
end subroutine zone_ranges
!
!
!
!
!
subroutine zone_distro(lvl)
!
! Distribute all necessary zonal information to every processor which will be participating
! in output.
!
! Called by: lb_output
! Calls:
! External calls: mpi_allgatherv, mpi_barrier
!
use precise
use mpi
use amr_info_holder, only : zone_storage
use amr_processes , only :self,commune,nprocs
implicit none
integer :: lvl,i,ier
integer :: n_elems,counter
!  integer :: stat(MPI_STATUS_SIZE)
integer,allocatable :: disp_array(:),rec_array(:)

!write(*,*) 'entering zone information distribution function.',self

counter = 0

n_elems = 3*zone_storage(lvl)%big_zones

if (allocated(disp_array)) deallocate(disp_array)
allocate(disp_array(0:nprocs-1))
if (allocated(rec_array)) deallocate(rec_array)
allocate(rec_array(0:nprocs-1))
do i = 0,nprocs-1
  rec_array(i) = zone_storage(lvl)%little_zones(i)*3
end do
do i=0,nprocs-1
  if (i==0) then
    disp_array(i) = 0
  else
    disp_array(i) = disp_array(i-1)+zone_storage(lvl)%little_zones(i-1)*3
  end if
end do
!  write(*,*) 'receiving array',self,lvl,rec_array
!  write(*,*) 'displacement array',self,lvl,disp_array
!
! The allgatherv distribute all necessary information for zone writing to every processor
!
call mpi_allgatherv(zone_storage(lvl)%lower(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_INT,&
  zone_storage(lvl)%lower,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(zone_storage(lvl)%higher(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_INT,&
  zone_storage(lvl)%higher,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(zone_storage(lvl)%lo_corner(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_DOUBLE_PRECISION,&
  zone_storage(lvl)%lo_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)
!
call mpi_allgatherv(zone_storage(lvl)%hi_corner(1:3,zone_storage(lvl)%zone_start(self):&
  zone_storage(lvl)%zone_end(self)),zone_storage(lvl)%little_zones(self),MPI_DOUBLE_PRECISION,&
  zone_storage(lvl)%hi_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)

call mpi_barrier(commune,ier)

deallocate(rec_array)
deallocate(disp_array)

!  write(*,*) 'lower test 3',self,'lower',zone_storage(lvl)%lower
!  write(*,*) 'higher test 3',self,'high',zone_storage(lvl)%higher
!  write(*,*) 'lo_corner check ',self,'lo_corner',zone_storage(lvl)%lo_corner
!  write(*,*) 'hi_corner check ',self,'hi_corner',zone_storage(lvl)%hi_corner
!  write(*,*) 'big_zones check ',self,'lvl',lvl,'big_zones',zone_storage(lvl)%big_zones
!  write(*,*) 'little_zones check ',self,'little_zones ',zone_storage(lvl)%little_zones
!  write(*,*) 'start and end',self,'lvl',lvl,zone_storage(lvl)%zone_start,&
!    zone_storage(lvl)%zone_end
end subroutine zone_distro
!
!
!
!
!
subroutine zone_dealloc_safety(fine,num,lvl_min,lvl_max)
!
! Deallocate the zone storage arrays after an output
!
!
! Called by: lb_output
! Calls:
!
use precise
use nml_output
use amr_info_holder, only : zone_storage
!use amrex_amr_module
implicit none
integer :: lvl
integer,intent(in) :: lvl_min,lvl_max,fine,num



do lvl = lvl_max,lvl_min,-1
! for initial output cases
  if (num == 0) then
    if (allocated(zone_storage(lvl)%zone_end)) deallocate(zone_storage(lvl)%zone_end)
    if (allocated(zone_storage(lvl)%zone_start)) deallocate(zone_storage(lvl)%zone_start)
    if (allocated(zone_storage(lvl)%little_zones)) deallocate(zone_storage(lvl)%little_zones)

    if (allocated(zone_storage(lvl)%lower)) deallocate(zone_storage(lvl)%lower)
    if (allocated(zone_storage(lvl)%higher)) deallocate(zone_storage(lvl)%higher)
    if (allocated(zone_storage(lvl)%lo_corner)) deallocate(zone_storage(lvl)%lo_corner)
    if (allocated(zone_storage(lvl)%hi_corner)) deallocate(zone_storage(lvl)%hi_corner)
    if (allocated(zone_storage(lvl)%output_number)) deallocate(zone_storage(lvl)%output_number)
! out of range, to prevent segfaults
  else if (lvl > lvl_max .or. lvl < lvl_min) then
    cycle
! clear the saved information, XXXX unneeded for static mesh, fix later
  else
    !write(*,*) 'deallocating for safety',lvl

    if (allocated(zone_storage(lvl)%zone_end)) deallocate(zone_storage(lvl)%zone_end)
    if (allocated(zone_storage(lvl)%zone_start)) deallocate(zone_storage(lvl)%zone_start)
    if (allocated(zone_storage(lvl)%little_zones)) deallocate(zone_storage(lvl)%little_zones)

    if (allocated(zone_storage(lvl)%lower)) deallocate(zone_storage(lvl)%lower)
    if (allocated(zone_storage(lvl)%higher)) deallocate(zone_storage(lvl)%higher)
    if (allocated(zone_storage(lvl)%lo_corner)) deallocate(zone_storage(lvl)%lo_corner)
    if (allocated(zone_storage(lvl)%hi_corner)) deallocate(zone_storage(lvl)%hi_corner)
    if (allocated(zone_storage(lvl)%output_number)) deallocate(zone_storage(lvl)%output_number)
!    write(*,*) 'sub-parts deallocated',lvl
  end if
end do

if (allocated(zone_storage)) deallocate(zone_storage)
!write(*,*) 'main user defined variable deallocated'
end subroutine zone_dealloc_safety
!
!
!
!
!
subroutine zone_coord_setup(index_file,index_base,index_zone,ggm_offset)
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
integer,intent(in),optional :: ggm_offset

coord_x = 1
coord_y = 2
coord_z = 3
if (present(ggm_offset)) then
  if (index_file < 1001) then
    index_file = 1000 + index_file - 1
  end if
end if
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

end subroutine zone_coord_setup
!
!
!
!
!
subroutine zone_coord_data_write(index_file,index_base,lvl_lo,lvl_hi,sv_num,ggm_offset)
!
! Write the coordinates for the zone
!
! Called by: output_setup,
! Calls:
! External calls: cgp_coord_data_write_f,cgp_error_exit_f
!
use precise
use cgns
use amrex_amr_module
use amrex_base_module
use grid_data, only: rdx
use amr_info_holder, only: zone_storage
use amr_processes, only: self,nprocs,commune
use output_data, only: x_verts,y_verts,z_verts
implicit none
integer ::a,b,c,lvl,out_min,out_max,lvl_hi,lvl_lo,i,ier
integer :: index_file,index_base,index_zone,sv_num
integer :: coord_x,coord_y,coord_z
integer(cgsize_t) ::isize(3,3),loc_hi(3),loc_lo(3)
integer,intent(in),optional :: ggm_offset

coord_x = 1
coord_y = 2
coord_z = 3
if (present(ggm_offset)) then
  if (index_file < 1001) then
    index_file = 1000 + index_file - 1
  end if
end if
!write(*,*) 'inside zone coordinate data write subroutine',self,lvl_lo,lvl_hi
!write(*,*) 'amrex limits',amrex_problo,amrex_probhi

do lvl = lvl_lo,lvl_hi

  !write(*,*) 'amrex dx',amrex_geom(lvl)%dx(1),amrex_geom(lvl)%dx(2),amrex_geom(lvl)%dx(3)

  do i = zone_storage(lvl)%zone_start(self),zone_storage(lvl)%zone_end(self)

    isize(1,1) = zone_storage(lvl)%higher(1,i)-zone_storage(lvl)%lower(1,i)+1
    isize(2,1) = zone_storage(lvl)%higher(2,i)-zone_storage(lvl)%lower(2,i)+1
    isize(3,1) = zone_storage(lvl)%higher(3,i)-zone_storage(lvl)%lower(3,i)+1

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
          x_verts(a,b,c) = zone_storage(lvl)%lo_corner(1,i) + (a-1)*rdx(lvl)
          y_verts(a,b,c) = zone_storage(lvl)%lo_corner(2,i) + (b-1)*rdx(lvl)
          z_verts(a,b,c) = zone_storage(lvl)%lo_corner(3,i) + (c-1)*rdx(lvl)
        end do
      end do
    end do

    loc_lo = (/1,1,1/)
    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)

    index_zone = zone_storage(lvl)%output_number(i,sv_num)
    write(*,*) 'repeat index info',index_file,index_base,index_zone
    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!
!    write(*,*) 'box corners',zone_storage(lvl)%lo_corner(1:3,i),x_verts(isize(1,1),isize(2,1),isize(3,1)),&
!      y_verts(isize(1,1),isize(2,1),isize(3,1)),z_verts(isize(1,1),isize(2,1),isize(3,1)),lvl
    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_x,&
      loc_lo,loc_hi,x_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'reading?'
    write(*,*) 'x coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
    !CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,byx_lo,byx_hi,y_verts,ier)
    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,&
      loc_lo,loc_hi,y_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
    write(*,*) 'y coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_z,&
      loc_lo,loc_hi, z_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
    write(*,*) 'z coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
    WRITE(*,*) 'Coordinate value write completed.'

    deallocate(x_verts)
    deallocate(y_verts)
    deallocate(z_verts)

  end do
end do

end subroutine zone_coord_data_write
!
!
!
!
!
subroutine zone_box_setup(lvl)
!
! Shares and retrieves all the necessary data for CGNS output
!
! Called by: initial_output,output_setup,save_write
! Calls: zone_ranges,zone_distro
! External calls: mpi_allgatherv,mpi_barrier,mpi_allreduce
!
use precise
use mpi
use amr_processes, only : count_big_zones,self,commune,nprocs
use amr_info_holder, only : zone_storage
use nml_output, only : num_saves
implicit none
integer :: lvl,loc_zones,ier,large_zones,max_zones,i,req,j
integer :: tag
integer,allocatable :: rec_array(:),disp_array(:)

! Get the number of zones for each processor
loc_zones = count_big_zones(lvl)
!write(*,*) 'inside ZBS, n_loc_box',self,lvl,loc_zones
! Sum those zones and make sure every processor has a copy so all the other
! objects can be properly allocated

CALL MPI_Allreduce(loc_zones,zone_storage(lvl)%big_zones,1,MPI_INTEGER, &
  MPI_SUM,commune,ier)
!write(*,*) 'whiffle'
!
! Distribute the number of zones owned by a processor at a given level.
! Then agglomerate and broadcast so all procs know everything about everyone
! else's zones.
!
if (allocated(zone_storage(lvl)%little_zones)) deallocate(zone_storage(lvl)%little_zones)
allocate(zone_storage(lvl)%little_zones(0:nprocs-1))
!
  zone_storage(lvl)%little_zones(self) = loc_zones
!  write(*,*) 'local zones',self,'lvl',lvl,'n little zone',loc_zones,zone_storage(lvl)%little_zones
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
  call mpi_allgatherv(zone_storage(lvl)%little_zones(self),1,&
    MPI_INT,zone_storage(lvl)%little_zones,&
    rec_array,disp_array,MPI_INT,commune,ier)
!
!
!
  call mpi_barrier(commune,ier)
  deallocate(rec_array)
  deallocate(disp_array)
!  write(*,*) 'Do you hear the people sing?'
!  write(*,*) 'things and stuff',self,lvl,zone_storage(lvl)%little_zones

  if (allocated(zone_storage(lvl)%zone_start)) deallocate(zone_storage(lvl)%zone_start)
  allocate(zone_storage(lvl)%zone_start(0:nprocs-1))

  if (allocated(zone_storage(lvl)%zone_end)) deallocate(zone_storage(lvl)%zone_end)
  allocate(zone_storage(lvl)%zone_end(0:nprocs-1))
!
! Loop to assign start and end values to the zones each processor owns
!
  do i = 0,nprocs-1
    if (i==0) then
      zone_storage(lvl)%zone_start(i) = 1
      zone_storage(lvl)%zone_end(i) = zone_storage(lvl)%little_zones(i)
    else
      zone_storage(lvl)%zone_start(i) = zone_storage(lvl)%zone_end(i-1)+1
      zone_storage(lvl)%zone_end(i) = zone_storage(lvl)%zone_start(i)+&
        zone_storage(lvl)%little_zones(i)-1
    end if
  end do
!  write(*,*) 'zone start and end',self,lvl,zone_storage(lvl)%zone_start,zone_storage(lvl)%zone_end
!  write(*,*) 'Singing the songs of angry men!'
!
! allreduce sums the local zones for each processor to find the total number of zones on a level
!

  if (allocated(zone_storage(lvl)%lower)) deallocate(zone_storage(lvl)%lower)
  if (allocated(zone_storage(lvl)%higher)) deallocate(zone_storage(lvl)%higher)
  if (allocated(zone_storage(lvl)%lo_corner)) deallocate(zone_storage(lvl)%lo_corner)
  if (allocated(zone_storage(lvl)%hi_corner)) deallocate(zone_storage(lvl)%hi_corner)
  if (allocated(zone_storage(lvl)%output_number)) deallocate(zone_storage(lvl)%output_number)
!
!  max_zones = zone_storage(lvl)%big_zones/nprocs+1
!
!  write(*,*) "I'm",self,'number of big zones',zone_storage(lvl)%big_zones,&
!    'at level ',lvl
  allocate(zone_storage(lvl)%lower(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%higher(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%lo_corner(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%hi_corner(3,zone_storage(lvl)%big_zones))
  allocate(zone_storage(lvl)%output_number(zone_storage(lvl)%big_zones,num_saves))
!
!
!
  call zone_ranges(lvl,zone_storage(lvl)%big_zones,zone_storage(lvl)%lower,&
    zone_storage(lvl)%higher,zone_storage(lvl)%lo_corner,zone_storage(lvl)%hi_corner,&
    zone_storage(lvl)%zone_start,zone_storage(lvl)%zone_end)
!
!
   call zone_distro(lvl)
!
!  write(*,*) self,zone_storage(lvl)%big_zones
!  write(*,*) self,zone_storage(lvl)%little_zones
!  write(*,*) self,zone_storage(lvl)%zone_start,zone_storage(lvl)%zone_end
!  write(*,*) self,zone_storage(lvl)%lower
!  write(*,*) self,zone_storage(lvl)%higher
!  write(*,*) self,zone_storage(lvl)%lo_corner
!  write(*,*) self,zone_storage(lvl)%hi_corner
!  write(*,*) self,zone_storage(lvl)%zones_owned

  !  zone_storage(lvl)%lower,&
  !  zone_storage(lvl)%higher)
end subroutine zone_box_setup
!
!
!
!
!
subroutine zone_output_num_setup(num,min_lvl,max_lvl)
!
!
!
!
!
use precise
use constants
use amr_info_holder, only: zone_storage
use amr_processes, only: self,nprocs,commune
use mpi
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
  do zone = zone_storage(lvl)%zone_start(self),zone_storage(lvl)%zone_end(self)
    zone_storage(lvl)%output_number(zone,num) = zone + num_lower_lvls
!    write(*,*) 'output numbers',zone_storage(lvl)%zone_start(self),&
!      zone_storage(lvl)%zone_end(self),self,zone_storage(lvl)%output_number(zone,num)

  end do

  num_lower_lvls = num_lower_lvls + zone_storage(lvl)%big_zones
!
! Set up the receiving arrays (how much info for each processor) and the displacement array
!   (where that info goes)
!
  do i=0,nprocs-1

    rec_array(i) = zone_storage(lvl)%little_zones(i)
    if (i==0) then
      disp_array(i) = 0
    else
      disp_array(i) = disp_array(i-1)+zone_storage(lvl)%little_zones(i-1)
    end if

  end do

  !write(*,*) 'local output numbers',zone_storage(lvl)%output_number,self

  call mpi_allgatherv(zone_storage(lvl)%output_number(zone_storage(lvl)%zone_start(self):&
    zone_storage(lvl)%zone_end(self),num),zone_storage(lvl)%little_zones(self),MPI_INTEGER,&
    zone_storage(lvl)%output_number,rec_array,disp_array,&
    MPI_INTEGER,commune,ier)
!  write(*,*) 'local output numbers',zone_storage(lvl)%output_number,&
!    zone_storage(lvl)%little_zones,self
!  write(*,*) 'mpi junk',rec_array,disp_array,self
!
!
!
end do

deallocate(disp_array)
deallocate(rec_array)

end subroutine zone_output_num_setup
!
!
!
!
!
!subroutine zone_ketchup()
!!
!!
!!
!use precise
!use constants
!use amr_info_holder, only: zone_storage!,multigrid_save,ggm_save_data
!implicit none
!integer :: i
!
!end subroutine
!
!subroutine zone_backfill()
!!
!!
!!
!use precise
!use constants
!use amr_info_holder, only: zone_storage!,multigrid_save,ggm_save_data
!implicit none
!integer :: i
!
!
!
!end subroutine



end module

!
!
!  allocate(disp_array(0:nprocs-1))
!  allocate(rec_array(0:nprocs-1))

!  do i = 0,nprocs-1
!    rec_array(i) = zone_storage(lvl)%little_zones(i)
!  end do
!
