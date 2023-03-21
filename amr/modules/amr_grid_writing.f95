module amr_grid_writing
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
public :: amr_grid_zones,amr_grid_distro,amr_grid_ranges,amr_grid_coord_data_write
public :: amr_grid_zone_writer,amr_grid_coord_setup,amr_save_writer




contains
!
!
!
!
!
!
! This is zone_box_setup from the zone_writer module
subroutine amr_grid_zones(lvl,num)
!
!
!
implicit none
integer :: lvl,loc_zones,ier,i,j,num
integer :: counter_thingy
integer,allocatable :: rec_array(:),disp_array(:)

! Get the number of zones for each processor
loc_zones = count_big_zones(lvl)
!! Sum those zones and make sure every processor has a copy so all the other
!! objects can be properly allocated
!
!if (loc_zones == 0) return
!
CALL MPI_Allreduce(loc_zones,grid_info(lvl,num)%big_zones,1,MPI_INTEGER, &
  MPI_SUM,commune,ier)
!!
!!
!!
!grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zones_per_level(lvl) = &
!     grid_info(lvl,num)%big_zones


!
! Distribute the number of zones owned by a processor at a given level.
! Then agglomerate and broadcast so all procs know everything about everyone
! else's zones.
!
if (allocated(grid_info(lvl,num)%little_zones)) deallocate(grid_info(lvl,num)%little_zones)
allocate(grid_info(lvl,num)%little_zones(0:nprocs-1))
!
  grid_info(lvl,num)%little_zones(self) = loc_zones
!  write(*,*) 'local zones',self,'lvl',lvl,'n little zone',loc_zones,grid_info(lvl,num)%little_zones
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
  call mpi_allgatherv(grid_info(lvl,num)%little_zones(self),1,&
    MPI_INT,grid_info(lvl,num)%little_zones,&
    rec_array,disp_array,MPI_INT,commune,ier)
!
!
!
  call mpi_barrier(commune,ier)
  deallocate(rec_array)
  deallocate(disp_array)
!  write(*,*) 'Do you hear the people sing?'
!  write(*,*) 'things and stuff',self,lvl,grid_info(lvl,num)%little_zones


!
! Loop to assign start and end values to the zones each processor owns
!
!  if (grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%regrid_between_outputs) then

    if (allocated(grid_info(lvl,num)%zone_start)) deallocate(grid_info(lvl,num)%zone_start)
    allocate(grid_info(lvl,num)%zone_start(0:nprocs-1))

    if (allocated(grid_info(lvl,num)%zone_end)) deallocate(grid_info(lvl,num)%zone_end)
    allocate(grid_info(lvl,num)%zone_end(0:nprocs-1))

    do i = 0,nprocs-1
      if (i==0) then
        grid_info(lvl,num)%zone_start(i) = base_info(num)%adapted_zone_start_number
        grid_info(lvl,num)%zone_end(i) = base_info(num)%adapted_zone_start_number + &
          grid_info(lvl,num)%little_zones(i) - 1

      else
        grid_info(lvl,num)%zone_start(i) = grid_info(lvl,num)%zone_end(i-1)+1
        grid_info(lvl,num)%zone_end(i) = grid_info(lvl,num)%zone_start(i)+&
          grid_info(lvl,num)%little_zones(i)-1
      end if
!      write(*,*) 'zone info ',self,i,lvl,base_info(num)%adapted_zone_start_number,&
!          grid_info(lvl,num)%little_zones(i),grid_info(lvl,num)%zone_start(i),&
!          grid_info(lvl,num)%zone_end(i)

!      if (i == self) then
!        write(*,*) 'le sigh',self,lvl,loc_zones,grid_info(lvl,num)%little_zones
!        write(*,*) 'zone starts ',self,grid_info(lvl,num)%zone_start
!        write(*,*) 'zone ends ',self,grid_info(lvl,num)%zone_end
!      end if
    end do
    write(*,*) 'little zones ',self,grid_info(lvl,num)%big_zones,grid_info(lvl,num)%little_zones
    write(*,*) 'zone starts ',self,grid_info(lvl,num)%zone_start
    write(*,*) 'zone ends ',self,grid_info(lvl,num)%zone_end
!
! XXXX
!
    counter_thingy = grid_info(lvl,num)%zone_start(0)
!    do i = 0,nprocs-1
!      do j =  grid_info(lvl,num)%zone_start(i),grid_info(lvl,num)%zone_end(i)
      do j = 1,grid_info(lvl,num)%big_zones
        grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(j) = counter_thingy
        write(*,*) 'beeble sneet',self,lvl,j,counter_thingy,grid_info(lvl,num)%&
          multigrid_save(save_steps(num)%output_counter)%zone_number(j)
        counter_thingy = counter_thingy + 1
      end do
!    end do

!    grid_info(lvl,num)%zone_start(i) = base_info(num)%adapted_zone_start_number = &
!         grid_info(lvl,num)%zone_start(nprocs-1) + grid_info(lvl,num)%little_zones(nprocs-1)-1

    if (allocated(grid_info(lvl,num)%lower)) deallocate(grid_info(lvl,num)%lower)
    if (allocated(grid_info(lvl,num)%higher)) deallocate(grid_info(lvl,num)%higher)
    if (allocated(grid_info(lvl,num)%lo_corner)) deallocate(grid_info(lvl,num)%lo_corner)
    if (allocated(grid_info(lvl,num)%hi_corner)) deallocate(grid_info(lvl,num)%hi_corner)
    if (allocated(grid_info(lvl,num)%output_number)) deallocate(grid_info(lvl,num)%output_number)

    allocate(grid_info(lvl,num)%lower(3,grid_info(lvl,num)%zone_start(0):&
      grid_info(lvl,num)%zone_end(nprocs-1)))
    allocate(grid_info(lvl,num)%higher(3,grid_info(lvl,num)%zone_start(0):&
      grid_info(lvl,num)%zone_end(nprocs-1)))
    allocate(grid_info(lvl,num)%lo_corner(3,grid_info(lvl,num)%zone_start(0):&
      grid_info(lvl,num)%zone_end(nprocs-1)))
    allocate(grid_info(lvl,num)%hi_corner(3,grid_info(lvl,num)%zone_start(0):&
      grid_info(lvl,num)%zone_end(nprocs-1)))
    allocate(grid_info(lvl,num)%output_number(save_steps(num)%n_outsteps,1))
!
    call amr_grid_ranges(lvl,num)
!
    call amr_grid_distro(lvl,num)
!
end subroutine
!
!
!
!
!
!
subroutine amr_grid_distro(lvl,num)
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

n_elems = 3*grid_info(lvl,num)%big_zones

if (allocated(disp_array)) deallocate(disp_array)
allocate(disp_array(0:nprocs-1))
if (allocated(rec_array)) deallocate(rec_array)
allocate(rec_array(0:nprocs-1))
do i = 0,nprocs-1
  rec_array(i) = grid_info(lvl,num)%little_zones(i)*3
end do
do i=0,nprocs-1
  if (i==0) then
    disp_array(i) = 0
  else
    disp_array(i) = disp_array(i-1)+grid_info(lvl,num)%little_zones(i-1)*3
  end if
end do
!  write(*,*) 'receiving array',self,lvl,num,rec_array
!  write(*,*) 'displacement array',self,lvl,num,disp_array
!
! The allgatherv distribute all necessary information for zone writing to every processor
!
!call mpi_allgatherv(grid_info(lvl,num)%lower(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%little_zones(self),MPI_INT,&
!  grid_info(lvl,num)%lower,rec_array,disp_array,&
!  MPI_INT,commune,ier)
!!
!call mpi_allgatherv(grid_info(lvl,num)%higher(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%little_zones(self),MPI_INT,&
!  grid_info(lvl,num)%higher,rec_array,disp_array,&
!  MPI_INT,commune,ier)
!!
!call mpi_allgatherv(grid_info(lvl,num)%lo_corner(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%little_zones(self),MPI_DOUBLE_PRECISION,&
!  grid_info(lvl,num)%lo_corner,rec_array,disp_array,&
!  MPI_DOUBLE_PRECISION,commune,ier)
!!
!call mpi_allgatherv(grid_info(lvl,num)%hi_corner(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%little_zones(self),MPI_DOUBLE_PRECISION,&
!  grid_info(lvl,num)%hi_corner,rec_array,disp_array,&
!  MPI_DOUBLE_PRECISION,commune,ier)
!
!write(*,*) 'doodle hickey',grid_info(lvl,num)%lower(1:3,grid_info(lvl,num)%zone_start(self):&
!  grid_info(lvl,num)%zone_end(self)),grid_info(lvl,num)%lower,self,lvl,grid_info(lvl,num)%zone_start(self),&
!  grid_info(lvl,num)%zone_end(self)
call mpi_allgatherv(grid_info(lvl,num)%lower(1:3,grid_info(lvl,num)%zone_start(self):&
  grid_info(lvl,num)%zone_end(self)),rec_array(self),MPI_INT,&
  grid_info(lvl,num)%lower,rec_array,disp_array,&
  MPI_INT,commune,ier)
  write(*,*) 'grumble'
!
call mpi_allgatherv(grid_info(lvl,num)%higher(1:3,grid_info(lvl,num)%zone_start(self):&
  grid_info(lvl,num)%zone_end(self)),rec_array(self),MPI_INT,&
  grid_info(lvl,num)%higher,rec_array,disp_array,&
  MPI_INT,commune,ier)
!
call mpi_allgatherv(grid_info(lvl,num)%lo_corner(1:3,grid_info(lvl,num)%zone_start(self):&
  grid_info(lvl,num)%zone_end(self)),rec_array(self),MPI_DOUBLE_PRECISION,&
  grid_info(lvl,num)%lo_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)
!
call mpi_allgatherv(grid_info(lvl,num)%hi_corner(1:3,grid_info(lvl,num)%zone_start(self):&
  grid_info(lvl,num)%zone_end(self)),rec_array(self),MPI_DOUBLE_PRECISION,&
  grid_info(lvl,num)%hi_corner,rec_array,disp_array,&
  MPI_DOUBLE_PRECISION,commune,ier)
!
!
!
call mpi_barrier(commune,ier)

deallocate(rec_array)
deallocate(disp_array)

end subroutine
!
!
!
!
!
subroutine amr_grid_ranges(lvl,num)


integer :: counter,lvl,num!_zones,num,zone_starts(0:nprocs-1),zone_ends(0:nprocs-1)
!,ind,zone_starts(0:nprocs-1),zone_ends(0:nprocs-1)
!integer :: lo_nodes(3,num_zones),hi_nodes(3,num_zones)
!real(kind=dp) :: min_corner(3,num_zones),max_corner(3,num_zones)
real(kind=dp) :: poco_loco(3),poco_hoco(3)
!  logical :: is_full(num_zones)

type(amrex_mfiter) :: mfi
type(amrex_box) :: brx

!write(*,*) 'entering zone ranges sub',self,lvl

counter = grid_info(lvl,num)%zone_start(self)
call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
do while (mfi%next())

  brx = mfi%validbox()

!  lo_nodes(1,counter) = brx%lo(1)
!  lo_nodes(2,counter) = brx%lo(2)
!  lo_nodes(3,counter) = brx%lo(3)
!  hi_nodes(1,counter) = brx%hi(1)
!  hi_nodes(2,counter) = brx%hi(2)
!  hi_nodes(3,counter) = brx%hi(3)

  grid_info(lvl,num)%lower(1:3,counter) = brx%lo
  grid_info(lvl,num)%higher(1:3,counter) = brx%hi
!    write(*,*) 'lo nodes',self,lvl,counter,lo_nodes(1,counter),lo_nodes(2,counter),lo_nodes(3,counter)
!    write(*,*) 'hi nodes',self,lvl,counter,hi_nodes(1,counter),hi_nodes(2,counter),hi_nodes(3,counter)
  poco_loco = get_real_coords(amrex_geom(lvl),brx%lo)
  poco_hoco = get_real_coords(amrex_geom(lvl),brx%hi)

!  min_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%lo)
!  max_corner(1:3,counter) = get_real_coords(amrex_geom(lvl),brx%hi)

!  min_corner(1,counter) = poco_loco(1)
!  min_corner(2,counter) = poco_loco(2)
!  min_corner(3,counter) = poco_loco(3)
!  max_corner(1,counter) = poco_hoco(1)
!  max_corner(2,counter) = poco_hoco(2)
!  max_corner(3,counter) = poco_hoco(3)
  grid_info(lvl,num)%lo_corner(1:3,counter) = poco_loco
  grid_info(lvl,num)%hi_corner(1:3,counter) = poco_hoco
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
subroutine amr_grid_zone_writer(lvl,num,index_file)
!
!
!
integer :: i,num,zone_counter,index_file,index_base,ier
integer :: lvl,n,index_zone
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
zone_counter = base_info(num)%adapted_zone_start_number

!write(*,*) 'pre-launch info',zone_counter,lvl,num,self
!do lvl = min_lvl,max_lvl
!
!  call zone_box_setup(lvl)
  do n = 0,nprocs-1
    do i=grid_info(lvl,num)%zone_start(n),grid_info(lvl,num)%zone_end(n)
!      zone_counter = zone_counter+1
      grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(i) = i
      grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_lvls(i) = lvl
!      write(*,*) 'zone stuff',self,n,i,lvl
!    write(*,*) 'Sing with me',self,zone_counter,i
      isize(1,1) = grid_info(lvl,num)%higher(1,i)-grid_info(lvl,num)%lower(1,i)+1
      isize(2,1) = grid_info(lvl,num)%higher(2,i)-grid_info(lvl,num)%lower(2,i)+1
      isize(3,1) = grid_info(lvl,num)%higher(3,i)-grid_info(lvl,num)%lower(3,i)+1
      isize(1,2) = isize(1,1)-1
      isize(2,2) = isize(2,1)-1
      isize(3,2) = isize(3,1)-1
      isize(1,3) = 0
      isize(2,3) = 0
      isize(3,3) = 0
!    write(*,*) 'isize check',self,isize
!    write(*,*) 'pre-flight check',self,n,i,grid_info(lvl,num)%lower(1:3,i),&
!      grid_info(lvl,num)%higher(1:3,i)
    write(zone_number,format_nm_lvl_and_zone) lvl,'_',i
    zonename = 'zone_lvl_'//TRIM(zone_number)


    grid_info(lvl,num)%output_number(i,1) = i
!    write(*,*) 'boo ',index_file,index_base,i,zonename,isize
!    call cg_goto_f(index_file,index_base,'Base_t',0,'Zone_t',i,'end')
!    write(*,*) 'jangle'
!    call cg_error_print_f()
!    write(*,*) 'zone output number',self,grid_info(lvl,num)%output_number(i,sv_num),zone_counter
!    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
!      Structured,grid_info(lvl,num)%output_number(i,sv_num),ier)
    !CALL cg_zone_write_f(index_file,index_base,zonename,isize,Structured,i,ier)
    index_zone = i
        write(*,*) 'writing zone, ',zonename,index_file,index_base,i,lvl,self
    CALL cg_zone_write_f(index_file,index_base,zonename,isize,Structured,index_zone,ier)
!    call cg_error_print_f()
!    write(*,*) 'writing zone, ',zonename,index_file,index_base,i,lvl,self
!  write(*,*) 'after writing zone',index_file,index_base,zone_counter,grid_info(lvl,num)%output_number(i,sv_num)
!   index_coord = index_zone
    call amr_grid_coord_setup(index_file,index_base,i)
!  write(*,*) 'clear for landing',index_file,index_base,n,i
!    call mpi_barrier(commune,ier)

!  write(*,*) 'writing zone 2,',self,zonename,index_file,index_base,i

      call mpi_barrier(commune,ier)
    end do
  end do

!do lvl = save_shapes(1,num),max_lvl
!
!  call amr_grid_coord_data_write(index_file,index_base,lvl,num)
!
!end do
!base_info(num)%adapted_zone_start_number = base_info(num)%adapted_zone_start_number + &
!  grid_info(lvl,num)%big_zones

!end do

end subroutine
!
!
!
!
!
!
subroutine amr_grid_coord_setup(index_file,index_base,index_zone)

use cgns
implicit none
integer :: coord_x,coord_y,coord_z,ier!,lvl,num
integer :: index_file,index_base,index_zone
!integer,intent(in),optional :: ggm_offset

coord_x = 1
coord_y = 2
coord_z = 3

!write(*,*) 'ruff ruff?',index_zone
!
  CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
   'CoordinateX',coord_x,ier)
!
!  write(*,*) 'x-coord written',index_zone,self
  call mpi_barrier(commune,ier)
!
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!
  CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
   'CoordinateY',coord_y,ier)
  IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!  write(*,*) 'y-coord written',index_zone,self
  call mpi_barrier(commune,ier)
!
CALL cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,&
  'CoordinateZ',coord_z,ier)

!write(*,*) 'z-coord written',index_zone,self
!
IF (ier /= CG_OK) CALL cgp_error_exit_f

end subroutine
!
!
!
!
!
subroutine amr_grid_coord_data_write(lvl,num,index_file)
!
!
!
implicit none
integer ::a,b,c,lvl,i,ier,num
integer :: index_file,index_base
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

  do i = grid_info(lvl,num)%zone_start(self),grid_info(lvl,num)%zone_end(self)

    isize(1,1) = grid_info(lvl,num)%higher(1,i)-grid_info(lvl,num)%lower(1,i)+1
    isize(2,1) = grid_info(lvl,num)%higher(2,i)-grid_info(lvl,num)%lower(2,i)+1
    isize(3,1) = grid_info(lvl,num)%higher(3,i)-grid_info(lvl,num)%lower(3,i)+1

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
          x_verts(a,b,c) = grid_info(lvl,num)%lo_corner(1,i) + (a-1)*rdx(lvl)
          y_verts(a,b,c) = grid_info(lvl,num)%lo_corner(2,i) + (b-1)*rdx(lvl)
          z_verts(a,b,c) = grid_info(lvl,num)%lo_corner(3,i) + (c-1)*rdx(lvl)
        end do
      end do
    end do

    loc_lo = (/1,1,1/)
    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)

!    index_zone = i
    write(*,*) 'repeat index info',index_file,index_base,i
!    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!
!    write(*,*) 'box corners',grid_info(lvl,num)%lo_corner(1:3,i),x_verts(isize(1,1),isize(2,1),isize(3,1)),&
!      y_verts(isize(1,1),isize(2,1),isize(3,1)),z_verts(isize(1,1),isize(2,1),isize(3,1)),lvl
    CALL cgp_coord_write_data_f(index_file,index_base,i,&
      coord_x,loc_lo,loc_hi,x_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'reading?'
!    write(*,*) 'x coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
    !CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,byx_lo,byx_hi,y_verts,ier)
    CALL cgp_coord_write_data_f(index_file,index_base,i,coord_y,&
      loc_lo,loc_hi,y_verts,ier)
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'y coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
    CALL cgp_coord_write_data_f(index_file,index_base,i,coord_z,&
      loc_lo,loc_hi, z_verts,ier)
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
subroutine amr_grid_output_nums(num,min_lvl,max_lvl)
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
  do zone = grid_info(lvl,num)%zone_start(self),grid_info(lvl,num)%zone_end(self)
    grid_info(lvl,num)%output_number(zone,1) = zone + num_lower_lvls
!    write(*,*) 'output numbers',grid_info(lvl,num)%zone_start(self),&
!      grid_info(lvl,num)%zone_end(self),self,grid_info(lvl,num)%output_number(zone,num)

  end do

  num_lower_lvls = num_lower_lvls + grid_info(lvl,num)%big_zones
!
! Set up the receiving arrays (how much info for each processor) and the displacement array
!   (where that info goes)
!
  do i=0,nprocs-1

    rec_array(i) = grid_info(lvl,num)%little_zones(i)
    if (i==0) then
      disp_array(i) = 0
    else
      disp_array(i) = disp_array(i-1)+grid_info(lvl,num)%little_zones(i-1)
    end if

  end do

  !write(*,*) 'local output numbers',grid_info(lvl,num)%output_number,self

  call mpi_allgatherv(grid_info(lvl,num)%output_number(grid_info(lvl,num)%zone_start(self):&
    grid_info(lvl,num)%zone_end(self),1),grid_info(lvl,num)%little_zones(self),MPI_INTEGER,&
    grid_info(lvl,num)%output_number,rec_array,disp_array,&
    MPI_INTEGER,commune,ier)
!  write(*,*) 'local output numbers',grid_info(lvl,num)%output_number,&
!    grid_info(lvl,num)%little_zones,self
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
subroutine amr_save_writer(lvl,num,fine,index_file,index_base)
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
use timing
use constants
use freestream_values
use cgns
use mpi
use amr_info_holder
use amr_processes, only: self,nprocs,commune
implicit none
INTEGER :: index_base,index_zone,index_file,index_field,ier,index_sol
integer :: n,fine,max_lvl,lvl,num,zee_count,total_zones,zer,counter
!integer,allocatable :: level_zones(:)
integer*8 :: loc_lo(4),loc_hi(4)
CHARACTER(len=32) :: zonename,char_save_number,solut_name,sol_num

type(amrex_box) :: bux
type(amrex_mfiter) :: mfi

if (fine > save_shapes(2,num)) then
  max_lvl = save_shapes(2,num)
else
  max_lvl = fine
end if

!write(*,*) 'writing save data'

!do lvl = save_shapes(1,num),max_lvl
! XXXX
!  do index_zone = 1,base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)
index_sol = save_steps(num)%output_counter
!
!
counter = 1

  do n = 0,nprocs-1
    do index_zone = grid_info(lvl,num)%zone_start(n),grid_info(lvl,num)%zone_end(n)
!    call cg_error_print_f
    write(*,*) 'before solution',self,lvl,index_file,index_base,index_zone,&
        save_steps(num)%sol_names(save_steps(num)%output_counter),&
        save_steps(num)%output_counter,index_sol,grid_info(lvl,num)%zone_start,&
        grid_info(lvl,num)%zone_end

      call cg_sol_write_f(index_file,index_base,index_zone,&
        save_steps(num)%sol_names(save_steps(num)%output_counter),Vertex,&
        index_sol,ier)
      IF (ier /= CG_OK) CALL cgp_error_exit_f
!
      grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_number(counter) = index_zone
      grid_info(lvl,num)%multigrid_save(save_steps(num)%output_counter)%zone_lvls(counter) = lvl
!

!    call cg_error_print_f
    write(*,*) 'after solution',self,lvl,index_file,index_base,index_zone,&
        save_steps(num)%sol_names(save_steps(num)%output_counter),&
        save_steps(num)%output_counter,index_sol

      counter = counter + 1
    end do
  end do
!end do
zee_count = 0
!do lvl = save_shapes(1,num),max_lvl

!  do index_zone = 1,base_info(num)%zones_per_lvl(lvl,save_steps(num)%output_counter)
  do n = 0,nprocs-1
    do index_zone = grid_info(lvl,num)%zone_start(n),grid_info(lvl,num)%zone_end(n)

    index_field = 1
!    write(*,*) 'output check, pre-write',index_file,index_base,index_zone,index_field,index_sol,lvl,self
    if (save_data(1,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'Density',index_field,ier)
      index_field = index_field+1
    end if
    write(*,*) 'output check, post-write',index_file,index_base,index_zone,index_field,index_sol,lvl,self
    if (save_data(2,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'U-velocity',index_field,ier)
      index_field = index_field+1
    end if

    if (save_data(3,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'V-velocity',index_field,ier)
      index_field = index_field+1
    end if

    if (save_data(4,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'W-velocity',index_field,ier)
      index_field = index_field+1
    end if

    if (save_data(5,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'Pressure',index_field,ier)
      index_field = index_field+1
    end if

    if (save_data(6,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'Temperature',index_field,ier)
      index_field = index_field+1
    end if
    !write(*,*) 'output check',index_file,index_base,save_steps(num)%output_counter,index_field,&
    !  grid_info(lvl,num)%output_number(index_zone,num),self
    if (save_data(7,num)) then
      call cgp_field_write_f(index_file,index_base,&
        index_zone,index_sol,RealDouble,'Alpha',index_field,ier)
      index_field = index_field+1
    end if
!    call cg_error_print_f()

    IF (ier /= CG_OK) CALL cgp_error_exit_f

!    write(*,*) 'deely bopper',index_zone,self,lvl

    call mpi_barrier(commune,ier)

    end do
  end do
!end do

!do lvl = save_shapes(1,num),max_lvl
index_zone = grid_info(lvl,num)%zone_start(self)
  call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
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
!
!
!
    index_field = 1
!        call cgp_field_write_data_f(index_file,index_base,&
!          grid_info(lvl,num)%zones_owned(c),&
!          index_flow,index_field,loc_lo,loc_hi,rho(bux%lo(1):bux%hi(1),&
!          bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
    if (save_data(1,num)) then
      rho => mfrho(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, rho',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,&
        index_zone,index_sol,index_field,loc_lo,loc_hi,rho(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if
    if (save_data(2,num)) then
      u_vel => mfu_vel(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, u',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,index_zone,&
        index_sol,index_field,loc_lo,loc_hi,u_vel(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if
    if (save_data(3,num)) then
      v_vel => mfv_vel(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, v',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,index_zone,&
        index_sol,index_field,loc_lo,loc_hi,v_vel(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if
    if (save_data(4,num)) then
      w_vel => mfw_vel(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, w',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,index_zone,&
        index_sol,index_field,loc_lo,loc_hi,w_vel(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if
    if (save_data(5,num)) then
      press => mfpress(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, pressure',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,index_zone,&
        index_sol,index_field,loc_lo,loc_hi,press(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if
    if (save_data(6,num)) then
      temp => mftemp(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, temp',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,index_zone,&
        index_sol,index_field,loc_lo,loc_hi,temp(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if

    if (save_data(7,num)) then
      alpha => mfalpha(lvl)%dataptr(mfi)
      !write(*,*) 'field data write, temp',index_file,index_base,save_steps(num)%output_counter,index_field,&
      !grid_info(lvl,num)%output_number(zee_count,num),self
      call cgp_field_write_data_f(index_file,index_base,index_zone,&
        index_sol,index_field,loc_lo,loc_hi,alpha(bux%lo(1):bux%hi(1),&
        bux%lo(2):bux%hi(2),bux%lo(3):bux%hi(3),1),ier)
      index_field = index_field+1
    end if
!      write(*,*) 'field data',self,index_file,index_base,b,&
!        index_flow,index_field,loc_lo,loc_hi
    IF (ier /= CG_OK) CALL cgp_error_exit_f
!    write(*,*) 'squizzle',index_zone,self,lvl
    index_zone = index_zone+1
  end do
!end do

  call amrex_mfiter_destroy(mfi)
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
!use output_data, only: ggm_steps
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
!    if (regrid_time(i) >= ggm_steps(num)%saved_times(j-1) .and. &
!        regrid_time(i) <= ggm_steps(num)%saved_times(j)) then
!
!      do lvl = 0,amrex_max_level
!        if (lvl <= fixed_level) then
!          ggm_steps(num)%multigrid_save(i)%regrid_between_outputs(lvl) = .false.
!        else
!          ggm_steps(num)%multigrid_save(i)%regrid_between_outputs(lvl) = .true.
!        end if
!      end do
!      exit
!
!    else
!      ggm_steps(num)%multigrid_save(i)%regrid_between_outputs(0:amrex_max_level) = .false.
!    end if
!
!  end do
!end do
!
!end subroutine
!
end module



!  write(*,*) 'zone start and end',self,lvl,grid_info(lvl,num)%zone_start,grid_info(lvl,num)%zone_end
!  write(*,*) 'Singing the songs of angry men!'
!
! allreduce sums the local zones for each processor to find the total number of zones on a level
!


!  if (allocated(grid_info(lvl,num)%output_number)) deallocate(grid_info(lvl,num)%output_number)
!
!  max_zones = grid_info(lvl,num)%big_zones/nprocs+1
!
!  write(*,*) "I'm",self,'number of big zones',grid_info(lvl,num)%big_zones,&
!    'at level ',lvl

!  allocate(grid_info(lvl,num)%output_number(grid_info(lvl,num)%big_zones,num_saves))
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
!  do i = grid_info(lvl,num)%zone_start(self),grid_info(lvl,num)%zone_end(self)
!
!    isize(1,1) = grid_info(lvl,num)%higher(1,i)-grid_info(lvl,num)%lower(1,i)+1
!    isize(2,1) = grid_info(lvl,num)%higher(2,i)-grid_info(lvl,num)%lower(2,i)+1
!    isize(3,1) = grid_info(lvl,num)%higher(3,i)-grid_info(lvl,num)%lower(3,i)+1
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
!          x_verts(a,b,c) = grid_info(lvl,num)%lo_corner(1,i) + (a-1)*rdx(lvl)
!          y_verts(a,b,c) = grid_info(lvl,num)%lo_corner(2,i) + (b-1)*rdx(lvl)
!          z_verts(a,b,c) = grid_info(lvl,num)%lo_corner(3,i) + (c-1)*rdx(lvl)
!        end do
!      end do
!    end do
!
!    loc_lo = (/1,1,1/)
!    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)
!
!    index_zone = grid_info(lvl,num)%output_number(i,sv_num)
!!    write(*,*) 'repeat index info',index_file,index_base,index_zone
!!    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!!
!!    write(*,*) 'box corners',grid_info(lvl,num)%lo_corner(1:3,i),x_verts(isize(1,1),isize(2,1),isize(3,1)),&
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
