subroutine initial_output
!
! Output the inital grid to ensure it is satisfactory. This also constitutes
! a major goal in testing to make sure Amrex, MPI, and CGNS all play well
! together.
!
!
! Called by:
! Calls:
! External calls:
!
use precise
use constants
use amrex_amr_module
use amrex_base_module
use amr_info_holder
use zone_writing
use amr_processes, only : self,nprocs,commune
use grid_data, only : rdx
use cgns
use mpi
implicit none
integer :: a,b,c,lvl,fine,num_verts,coord_x,coord_y,coord_z,i
real(kind=dp) :: init_loc(3)
integer*8 :: isize(3,3)
integer,allocatable :: node_number(:,:,:,:)
integer(cgsize_t) :: byx_lo(3),byx_hi(3),loc_lo(4),loc_hi(4),idata(2)
integer :: index_base,index_zone,index_coord,index_file,icelldim,iphysdim,&
  ier,index_flow,index_field,index_section,zone_counter,zee_count
CHARACTER(len=32) :: solut_name,filename,lvl_number,format_nm_lvl_and_zone
character(len=32) :: zonename,basename,zone_number

type(amrex_mfiter) :: mfi
type(amrex_box) :: byx

format_nm_lvl_and_zone = '(I2.2,a1,I3.3)'
write(*,*) 'inside initial output subroutine'

!WRITE(char_save_number,format_junk) a
!WRITE(solut_save_number,format_junk_ts) t

filename = 'initial_grid_output.cgns'
basename = 'drop the'
solut_name = 'State Information'
fine = amrex_get_finest_level()
allocate(zone_storage(0:amrex_max_level))

index_base = 1
index_file = 1
index_coord = 1
index_zone = 1
index_section = 1
index_flow = 1
index_field = 1

IF (dimensions == 3) THEN
  icelldim = 3
  iphysdim = 3
ELSE if (dimensions ==2 ) then
  icelldim = 2
  iphysdim = 3
END IF

do lvl = 0,fine
  call zone_box_setup(lvl)
end do
write(*,*) 'returning from box setup',self
call cgp_mpi_comm_f(commune,ier)
call cgp_pio_mode_f(CGP_INDEPENDENT,ier)
CALL cgp_open_f(filename,CG_MODE_WRITE,index_file,ier)
IF (ier /= CG_OK) CALL cgp_error_exit_f

CALL cg_base_write_f(index_file,basename,icelldim,iphysdim,&
  index_base,ier)
IF (ier /= CG_OK) CALL cgp_error_exit_f
!WRITE(11,*) 'Base write done'
!WRITE(*,*) 'Base write done'
zone_counter = 0

do lvl = 0,fine
  do i = 1,zone_storage(lvl)%big_zones
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
    write(*,*) 'isize check',self,isize
    write(*,*) 'pre-flight check',self,zone_counter,i,zone_storage(lvl)%big_zones,lvl
    write(zone_number,format_nm_lvl_and_zone) lvl,'_',zone_counter
    zonename = 'zone_lvl_'//TRIM(zone_number)

    zone_storage(lvl)%output_number(i,1) = zone_counter

!    write(*,*) 'zone output number',self,zone_storage(lvl)%output_number(i,loc_),zone_counter
    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
      Structured,zone_counter,ier)
!  write(*,*) 'after writing zone',index_file,index_base,zone_counter,index_coord
!    index_coord = index_zone
    call zone_coord_setup(index_file,index_base,zone_counter)
  write(*,*) 'indices',index_file,index_base,index_zone,index_coord
    call mpi_barrier(commune,ier)
  end do
end do

coord_x = 1
coord_y = 2
coord_z = 3

call zone_coord_data_write(index_file,index_base,0,fine,1)

! XXXX
do lvl = 0,fine
  do index_zone = 1,zone_storage(lvl)%big_zones
  write(*,*) self,index_file,index_base,index_zone,index_flow,solut_name
    zee_count = zee_count+1

    call cg_sol_write_f(index_file,index_base,zone_storage(lvl)%output_number(index_zone,1),&
      solut_name,Vertex,index_flow,ier)

    write(*,*) 'after solution',self,index_file,index_base,&
      zone_storage(lvl)%output_number(index_zone,1),index_flow
    IF (ier /= CG_OK) CALL cgp_error_exit_f
  end do
end do
!
! Write the field name on every process
!
zee_count = 0
do lvl = 0,fine
  do index_zone = 1,zone_storage(lvl)%big_zones

    call cgp_field_write_f(index_file,index_base,&
      zone_storage(lvl)%output_number(index_zone,1),index_flow,&
      INTEGER,'State and level',index_field,ier)

    IF (ier /= CG_OK) CALL cgp_error_exit_f
  end do
end do
!
! Output the data for every zone by processor
!
do lvl = 0,fine
zee_count = zone_storage(lvl)%zone_start(self)
  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
  do while (mfi%next())
!    c = c+1
    state => mfstate(lvl)%dataptr(mfi)

    byx = mfi%validbox()
!    write(*,*) 'box bounds',byx%lo,byx%hi,lvl

!    if (zone_storage(lvl)%full(b)) then
!        write(*,*) self,zone_storage(lvl)%full,c
!        write(*,*) self,zone_storage(lvl)%zones_owned
!        write(*,*) 'bounds of state',lbound(state),ubound(state)
!        write(*,*) 'bounds of the box',byx%lo,byx%hi
!        write(*,*) 'rho at low',rho(byx%lo(1),byx%lo(2),byx%lo(3),1)
!        write(*,*) 'rho at hi',rho(byx%hi(1),byx%hi(2),byx%hi(3),1)
    loc_lo = 1
    loc_hi(1) = byx%hi(1)-byx%lo(1)+1
    loc_hi(2) = byx%hi(2)-byx%lo(2)+1
    loc_hi(3) = byx%hi(3)-byx%lo(3)+1
    loc_hi(4) = 1
!        call cgp_field_write_data_f(index_file,index_base,&
!          zone_storage(lvl)%zones_owned(c),&
!          index_flow,index_field,loc_lo,loc_hi,rho(byx%lo(1):byx%hi(1),&
!          byx%lo(2):byx%hi(2),byx%lo(3):byx%hi(3),1),ier)
    call cgp_field_write_data_f(index_file,index_base,&
      zone_storage(lvl)%output_number(zee_count,1),&
      index_flow,index_field,loc_lo,loc_hi,state(byx%lo(1):byx%hi(1),&
      byx%lo(2):byx%hi(2),byx%lo(3):byx%hi(3),1),ier)
!      write(*,*) 'field data',self,index_file,index_base,b,&
!        index_flow,index_field,loc_lo,loc_hi
    IF (ier /= CG_OK) CALL cgp_error_exit_f
    zee_count = zee_count+1
  end do
  call amrex_mfiter_destroy(mfi)
end do


call cgp_close_f(index_file,ier)

call zone_dealloc_safety(fine,0,0,amrex_max_level)

write(*,*) 'initial output write completed',self,index_file

end subroutine
!do lvl = 0,fine
!  do i = zone_storage(lvl)%zone_start(self),zone_storage(lvl)%zone_end(self)
!
!    isize(1,1) = zone_storage(lvl)%higher(1,i)-zone_storage(lvl)%lower(1,i)+1
!    isize(2,1) = zone_storage(lvl)%higher(2,i)-zone_storage(lvl)%lower(2,i)+1
!    isize(3,1) = zone_storage(lvl)%higher(3,i)-zone_storage(lvl)%lower(3,i)+1
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
!
!!    index_zone = i
!
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
!          x_verts(a,b,c) = zone_storage(lvl)%lo_corner(1,i) + (a-1)*rdx(lvl)
!          y_verts(a,b,c) = zone_storage(lvl)%lo_corner(2,i) + (b-1)*rdx(lvl)
!          z_verts(a,b,c) = zone_storage(lvl)%lo_corner(3,i) + (c-1)*rdx(lvl)
!        end do
!      end do
!    end do
!
!    loc_lo = (/1,1,1/)
!    loc_hi = (/isize(1,1),isize(2,1),isize(3,1)/)
!!    write(*,*) 'repeat index info',index_file,index_base,index_zone
!!    write(*,*) 'preparing coordinate write',(/1,1,1/),(/isize(1,1),isize(2,1),isize(3,1)/)
!    index_zone = zone_storage(lvl)%output_number(i,1)
!!
!    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_x,&
!      loc_lo,loc_hi,x_verts,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!    write(*,*) 'reading?'
!!    write(*,*) 'x coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
!    !CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,byx_lo,byx_hi,y_verts,ier)
!    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_y,&
!      loc_lo,loc_hi,&
!      y_verts,ier)
!    IF (ier /= CG_OK) CALL cgp_error_exit_f
!!    write(*,*) 'y coordinate data written.',index_file,index_base,index_zone,coord_x,coord_y,coord_z
!    index_coord = index_coord+1
!    CALL cgp_coord_write_data_f(index_file,index_base,index_zone,coord_z,&
!      loc_lo,loc_hi,&
!      z_verts,ier)
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
! Write the solution name on every process
!

!do lvl = 0,fine
!
!  write(lvl_number,format_junk_lvl) lvl
!
!  call amrex_mfiter_built(mfi,mfrho(lvl),tiling=.false.)
!  do while (mfi%next())
!
!    byx = mfi%tilebox()
!    byx_lo = byx%lo
!    byx_hi = byx%hi
!
!    num_verts = (byx_hi(1)-byx_lo(1)+1)*(byx_hi(2)-byx_lo(2)+1)*(byx_hi(3)-&
!      byx_lo(3)+1)
!    isize(1,1) = byx_hi(1)-byx_lo(1)
!    isize(2,1) = byx_hi(2)-byx_lo(2)
!    isize(3,1) = byx_hi(3)-byx_lo(3)
!
!    isize(1,2) = isize(1,1)+1
!    isize(2,2) = isize(2,1)+1
!    isize(3,2) = isize(3,1)+1
!
!    isize(1,3) = 0
!    isize(2,3) = 0
!    isize(3,3) = 0
!
!    CALL cg_zone_write_f(index_file,index_base,zonename,isize,&
!      Structured,lvl+1,ier)
!    allocate(x_verts(byx_lo(1):byx_hi(1),byx_lo(2):byx_hi(2),byx_lo(3):byx_hi(3)))
!    allocate(y_verts(byx_lo(1):byx_hi(1),byx_lo(2):byx_hi(2),byx_lo(3):byx_hi(3)))
!    allocate(z_verts(byx_lo(1):byx_hi(1),byx_lo(2):byx_hi(2),byx_lo(3):byx_hi(3)))
!
!    init_loc = get_real_coords(amrex_geom,(/byx_lo(1),byx_lo(2),byx_lo(3)/))
!
!    do c = byx_lo(3),byx_hi(3)
!      do b = byx_lo(2),byx_hi(2)
!        do a = byx_lo(1),byx_hi(1)
!
!          x_verts(a,b,c) = init_loc(1)+(a-byx_lo(1))*rdx(lvl)
!          y_verts(a,b,c) = init_loc(2)+(b-byx_lo(2))*rdx(lvl)
!          z_verts(a,b,c) = init_loc(3)+(c-byx_lo(3))*rdx(lvl)
!
!        end do
!      end do
!    end do
!
!    CALL cg_coord_write_f(index_file,index_base,lvl+1,RealDouble,&
!     'CoordianteX',x_verts,n_coords,ier)
!    CALL cg_coord_write_f(index_file,index_base,lvl+1,RealDouble,&
!     'CoordianteY',y_verts,n_coords,ier)
!    CALL cg_coord_write_f(index_file,index_base,lvl+1,RealDouble,&
!     'CoordianteZ',z_verts,n_coords,ier)
!
!
!    deallocate(x_verts)
!    deallocate(y_verts)
!    deallocate(z_verts)
!
!
!
!  end do
!
!  n_coords = n_coords + 1
!
!end do

!CALL cgp_open_f(filename,CG_MODE_WRITE,index_file,ier)
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!!
!CALL cg_base_write_f(index_file,basename,icelldim,iphysdim,&
!  index_base,ier)
!IF (ier /= CG_OK) CALL cgp_error_exit_f

!CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
!  ier)
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!CALL cg_biter_write_f(index_file,index_base,'TimeIterValues',&
!  n_timesteps,ier)
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!write(*,*) 'Writing time iterations values to finalize file.'
!
!!CALL cg_simulation_type_write_f(index_file,index_base,TimeAccurate,&
!!  ier)
!!IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!CALL cg_goto_f(index_file,index_base,ier,'BaseIterativeData_t',&
!  1,'end')
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!boing = n_timesteps
!min_boing = 1
!write(*,*) 'preparing to write time values'
!
!CALL cgp_array_write_f('TimeValues',RealDouble,min_boing,boing,&
!  index_array,ier)
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!call cgp_array_write_data_f(index_array,min_boing,boing,&
!  time_vals,ier)
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!call cg_goto_f(index_file,index_base,'Zone_t',1,"end")
!IF (ier /= CG_OK) CALL cgp_error_exit_f
!
!write(*,*) 'time values written',self,time_vals,n_timesteps
!CALL cg_ziter_write_f(index_file,index_base,index_zone,&
! 'ZoneIterativeData',ier)
!IF (ier /= CG_OK) CALL cg_error_exit_f
!
! Go to the zone data
!
!CALL cg_goto_f(index_file,index_base,ier,'Zone_t',index_zone,&
!  'ZoneIterativeData_t',1,'end')
!IF (ier /= CG_OK) CALL cg_error_exit_f
!
! Tell it the solution names and their sizes
!
!CALL cg_array_write_f('FlowSolutionPointers',Character,2,idata,&
!  solut_name,ier)
!IF (ier /= CG_OK) CALL cg_error_exit_f
!
!do lvl = 0,fine
!  write(*,*) 'checking data before data write'
!  write(*,*) 'lower check',self,'lower',zone_storage(lvl)%lower
!  write(*,*) 'higher check',self,'high',zone_storage(lvl)%higher
!  write(*,*) 'lo_corner check',self,'lo_corner',zone_storage(lvl)%lo_corner
!  write(*,*) 'hi_corner check',self,'hi_corner',zone_storage(lvl)%hi_corner
!  write(*,*) 'big_zones check',self,'big_zones',zone_storage(lvl)%big_zones
!  write(*,*) 'output number',self,'output nums',zone_storage(lvl)%output_number
!end do
