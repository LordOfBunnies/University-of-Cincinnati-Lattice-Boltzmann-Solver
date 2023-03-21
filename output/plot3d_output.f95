subroutine plot3d_output(fine)
!
!
!
!
!
!
!
!
!
use precise
use constants
use mpi
use ggm_stuff
use amrex_amr_module
use amrex_base_module
use output_data
use nml_output
use grid_data
use linkwise
use freestream_values
use amr_info_holder
use amr_processes, only: self,nprocs,commune,count_big_zones,total_big_zones,get_real_coords,&
  get_box_sizes
implicit none
integer :: a,b,c,n,g,i,lvl,max_lvl,loc_zones,ier,total_blocks,max_boxes
integer :: fine,counter,int_dummy,n_comps
integer,allocatable :: lvl_boxes(:,:),rec_array(:),disp_array(:),num_blocks(:),&
  box_dims(:,:),block_starts(:,:),block_ends(:,:),x_dims(:),y_dims(:),z_dims(:),blanks(:,:,:)
real(kind=dp) :: start_loc(3),real_dummy
!real(kind=dp),contiguous,pointer :: rho_cells(:,:,:,:)
real,allocatable :: x_points(:,:,:),y_points(:,:,:),z_points(:,:,:)
character(len=32) :: grid_file,data_file,ts_format,grid_number,data_number

type(amrex_mfiter) :: mfi
type(amrex_box) :: mux
!type(amrex_multifab),allocatable :: mfcell_rho(:)
!
!
!
!write(*,*) 'Plot3D output phase'
!
!
!
real_dummy = 0.0D0
int_dummy = 0
n_comps = 6

if (fine > save_shapes(2,1)) then
  max_lvl = save_shapes(2,1)
else
  max_lvl = fine
end if
!
!if (allocated(mfcell_rho)) deallocate(mfcell_rho)
!allocate(mfcell_rho(0:amrex_max_level))
!
!
!
!do lvl = 0,max_lvl
!  !write(*,*) 'squash',max_lvl,self,lvl,time(fine),timestep
!  call amrex_multifab_build(mfcell_rho(lvl),mfrho(lvl)%ba,mfrho(lvl)%dm,1,0,&
!    (/.false.,.false.,.false./))
!end do

!call node_to_cell_averaging(mfcell_rho(0:max_lvl),max_lvl)

ts_format = '(I9.9)'

write(grid_number,ts_format) timestep(amrex_max_level)
write(data_number,ts_format) timestep(amrex_max_level)

grid_file = 'uclbs_grid_'//trim(grid_number)//'.x'
data_file = 'uclbs_data_'//trim(data_number)//'.q'
!write(*,*) grid_file, data_file
open(file=grid_file,unit=17,status='REPLACE',form='UNFORMATTED',action='WRITE',&
  ACCESS = 'APPEND',iostat=ier)
open(file=data_file,unit=18,status='REPLACE',form='UNFORMATTED',action='WRITE',&
  ACCESS = 'APPEND',iostat=ier)
!open(file='test.x',unit=19,status='REPLACE',form='UNFORMATTED',action='WRITE',&
!  ACCESS = 'SEQUENTIAL',iostat=ier)
!
!
!
if (allocated(lvl_boxes)) deallocate(lvl_boxes)
allocate(lvl_boxes(0:nprocs-1,0:amrex_max_level))
if (allocated(num_blocks)) deallocate(num_blocks)
allocate(num_blocks(0:max_lvl))
if (allocated(block_starts)) deallocate(block_starts)
allocate(block_starts(0:nprocs-1,0:amrex_max_level))
if (allocated(block_ends)) deallocate(block_ends)
allocate(block_ends(0:nprocs-1,0:amrex_max_level))
!
!
!
if (allocated(rec_array)) deallocate(rec_array)
allocate(rec_array(0:nprocs-1))
if (allocated(disp_array)) deallocate(disp_array)
allocate(disp_array(0:nprocs-1))
rec_array = 1
do i = 0,nprocs-1
  disp_array(i) = i
end do


total_blocks = 0
do lvl = 0,max_lvl
  loc_zones = count_big_zones(lvl)
!  write(*,*) 'local blocks',loc_zones,self,lvl
  CALL MPI_Allreduce(loc_zones,num_blocks(lvl),1,MPI_INTEGER,MPI_SUM,commune,ier)
  call mpi_allgatherv(loc_zones,1,MPI_INT,lvl_boxes(0:nprocs-1,lvl),&
    rec_array,disp_array,MPI_INT,commune,ier)

  total_blocks = total_blocks + num_blocks(lvl)

  do i = 0,nprocs-1
    if (i==0) then
      block_starts(i,lvl) = 1
      block_ends(i,lvl) = lvl_boxes(i,lvl)
    else
      block_starts(i,lvl) = block_ends(i-1,lvl) + 1
      block_ends(i,lvl) = block_starts(i,lvl)+lvl_boxes(i,lvl)-1
    end if
  end do
end do
!
!
!
if (self == 0) then
  write(17) total_blocks
  write(18) total_blocks
!  write(*,*) 'total blocks',total_blocks
end if

!
!
!
if (self == 0) then
    if (allocated(x_dims)) deallocate(x_dims)
    allocate(x_dims(total_blocks))
    if (allocated(y_dims)) deallocate(y_dims)
    allocate(y_dims(total_blocks))
    if (allocated(z_dims)) deallocate(z_dims)
    allocate(z_dims(total_blocks))
end if
counter = 1
do lvl = 0,max_lvl
  if (allocated(box_dims)) deallocate(box_dims)
  allocate(box_dims(3,num_blocks(lvl)))

  box_dims(1:3,block_starts(self,lvl):block_ends(self,lvl)) = get_box_sizes(lvl_boxes(self,lvl),lvl)

  do i = 0,nprocs-1
    rec_array(i) = lvl_boxes(i,lvl)*3
  end do
  do i=0,nprocs-1
    if (i==0) then
      disp_array(i) = 0
    else
      disp_array(i) = disp_array(i-1)+lvl_boxes(i-1,lvl)*3
    end if
  end do

!  write(*,*) 'test data',rec_array,disp_array,self,lvl
!  write(*,*) 'start/stop',block_starts(self,lvl),block_ends(self,lvl),self,lvl
!  write(*,*) 'block data',box_dims(1:3,block_starts(self,lvl):block_ends(self,lvl)),&
!    lvl_boxes(self,lvl)
!  WRITE(*,*) 'lvl boxes',lvl_boxes(0:nprocs-1,lvl),self

  call mpi_allgatherv(box_dims(1:3,block_starts(self,lvl):block_ends(self,lvl)),&
    rec_array(self),MPI_INT,box_dims(1:3,1:num_blocks(lvl)),rec_array(0:nprocs-1),&
    disp_array(0:nprocs-1),MPI_INT,commune,ier)
!  write(*,*) 'erp erp',ier,self
!
!
!
  if (self == 0) then

    do i = 1,num_blocks(lvl)
      x_dims(counter) = box_dims(1,i)
      y_dims(counter) = box_dims(2,i)
      z_dims(counter) = box_dims(3,i)
      counter = counter + 1
!      write(*,*) 'arcc test box',x_dims(i),y_dims(i),z_dims(i)
    end do
  end if

  call mpi_barrier(commune,ier)
!  do c = box_dims
  deallocate(box_dims)
end do
!
! Write out the dimensions
!
if (self == 0) then
    write(17) (x_dims(i),y_dims(i),z_dims(i),i=1,total_blocks)
    if (use_ggm .and. ggm_output) then
      write(18) (x_dims(i),y_dims(i),z_dims(i),i=1,total_blocks),7,0
    else
      write(18) (x_dims(i),y_dims(i),z_dims(i),i=1,total_blocks),6,0
    end if
!    write(*,*) (x_dims(i),y_dims(i),z_dims(i),i=1,total_blocks)
!    if (lvl == 0) then
!      write(19) x_dims(1),y_dims(1),z_dims(1)
!      write(*,*) x_dims(1),y_dims(1),z_dims(1)
!    end if
end if

call mpi_barrier(commune,ier)

close(17)
close(18)

if (use_ggm .and. ggm_output) then
  n_comps = n_comps + num_ggm_vars

  do g = 1,num_ggm_vars
!    write(*,*) 'bubble?',self
    call ggm_output_shortcut(g,fine)
  end do
end if
!
! Write data to the output files, this is unfortunately done a very stupid way
!
do lvl = 0,max_lvl
  do n = 0,nprocs-1
!
!
!
    call amrex_mfiter_build(mfi,mfrho(lvl),tiling=.false.)
    do while(mfi%next())
      mux = mfi%validbox()
      state => mfstate(lvl)%dataptr(mfi)
      rho => mfrho(lvl)%dataptr(mfi)
      u_vel => mfu_vel(lvl)%dataptr(mfi)
      v_vel => mfv_vel(lvl)%dataptr(mfi)
      w_vel => mfw_vel(lvl)%dataptr(mfi)
      temp => mftemp(lvl)%dataptr(mfi)
      alpha => mfalpha(lvl)%dataptr(mfi)
      if (use_ggm) ggm => mfggm(lvl)%dataptr(mfi)
      start_loc = get_real_coords(amrex_geom(lvl),mux%lo)
!
!
! Single processor open the file, so it goes to the end of the file correctly
!
!
      if (self == n) then

        open(file=grid_file,unit=17,status='OLD',form='UNFORMATTED',action='WRITE',&
          ACCESS = 'APPEND',iostat=ier)
        open(file=data_file,unit=18,status='OLD',form='UNFORMATTED',action='WRITE',&
          ACCESS = 'APPEND',iostat=ier)


        !write(18) mach,real_dummy,reynolds,time(fine)
        write(18) mach,real_dummy,reynolds,time(fine),gama,real_dummy,t_ref,int_dummy,&
          stag_temperature,const_vol,const_press,non_dim_R,non_dim_R,mach,time(fine),&
          timestep(fine)
!        write(*,*) 'Thor the noisy dog',self,lvl,mux%lo,mux%hi,start_loc,amrex_geom(lvl)%dx(1)

        if (allocated(x_points)) deallocate(x_points)
        allocate(x_points(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3)))
        if (allocated(y_points)) deallocate(y_points)
        allocate(y_points(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3)))
        if (allocated(z_points)) deallocate(z_points)
        allocate(z_points(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3)))
        if (allocated(blanks)) deallocate(blanks)
        allocate(blanks(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3)))

        do c = mux%lo(3),mux%hi(3)
          do b = mux%lo(2),mux%hi(2)
            do a = mux%lo(1),mux%hi(1)
              x_points(a,b,c) = start_loc(1)+amrex_geom(lvl)%dx(1)*(a-mux%lo(1))
              y_points(a,b,c) = start_loc(2)+amrex_geom(lvl)%dx(2)*(b-mux%lo(2))
              z_points(a,b,c) = start_loc(3)+amrex_geom(lvl)%dx(3)*(c-mux%lo(3))
              if (state(a,b,c,1) >= 0) then
                blanks(a,b,c) = 1
              else if (state(a,b,c,1) == -666) then
                blanks(a,b,c) = -1
              else
                blanks(a,b,c) = 0
              end if
!              write(*,*) 'grid coords',x_points(a,b,c),y_points(a,b,c),z_points(a,b,c),&
!                a,b,c,mux%lo,mux%hi
            end do
          end do
        end do

!        write(17) (((x_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!          (((y_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!          (((z_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
        write(17) (((x_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((y_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((z_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((blanks(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
if (use_ggm) then
        write(18) (((rho(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((u_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((v_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((w_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((temp(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((alpha(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((ggm(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
else
        write(18) (((rho(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((u_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((v_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((w_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((temp(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
          (((alpha(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
end if

!          ((((ggm(a,b,c,g),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!            g=1,num_ggm_vars)

        deallocate(x_points)
        deallocate(y_points)
        deallocate(z_points)
        deallocate(blanks)

        close(17)
        close(18)

      end if
    end do
  call mpi_barrier(commune,ier)
  call amrex_mfiter_destroy(mfi)

  end do

end do

!if (use_ggm .and. ggm_output) then
!  do lvl = 0,max_lvl
!    call amrex_multifab_destroy(mfgm(lvl))
!    call amrex_multifab_destroy(mfggm(lvl))
!  end do
!end if

!close(17)
!close(18)
!close(19)

deallocate(lvl_boxes)
deallocate(block_starts)
deallocate(block_ends)
!do lvl = max_lvl,0,-1
!  call amrex_multifab_destroy(mfcell_rho(lvl))
!end do
!
!deallocate(mfcell_rho)
end subroutine

!  do n = 1,max_boxes
!    mux = mfi%dataptr()
!
!    rho_cell => mfcell_rho(lvl)%dataptr(mfi)
!
!    do i = 0,nprocs-1
!      if (self == i) then
!        if (counter <= lvl_boxes(self,lvl)) then
!
!          start_loc = get_real_coords(amrex_geom(lvl),mux%lo)
!
!          do c = mux%lo(3),mux%hi(3)
!            do b = mux%lo(2),mux%hi(2)
!              do a = mux%lo(1),mux%hi(1)
!
!                write(17,*) start_loc(1)+amrex_geom(lvl)%dx(1)*(a-1),&
!                  start_loc(2)+amrex_geom(lvl)%dx(2)*(b-1),&
!                  start_loc(3)+amrex_geom(lvl)%dx(3)*(c-1)
!              end do
!            end do
!          end do
!
!
!        end if
!      end if
!
!
!    end do
!    call mpi_barrier(commune,ier)
!    counter = counter + 1
!  end do
!        write(18) (((rho(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
!        write(18) (((u_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
!        write(18) (((v_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
!        write(18) (((w_vel(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
!        write(18) (((temp(a,b,c,1),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))

!        do c = mux%lo(3),mux%hi(3)
!          do b = mux%lo(2),mux%hi(2)
!            do a = mux%lo(1),mux%hi(1)
!        write(17) x_points(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3)),&
!          y_points(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3)),&
!          z_points(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3))
!    write(17) box_dims(1,1:num_blocks(lvl)),box_dims(2,1:num_blocks(lvl)),&
!      box_dims(3,1:num_blocks(lvl))
!    write(18) box_dims(1,1:num_blocks(lvl)),box_dims(2,1:num_blocks(lvl)),&
!      box_dims(3,1:num_blocks(lvl))
!
!    deallocate(x_dims)
!    deallocate(y_dims)
!    deallocate(z_dims)
                !write(17) start_loc(1)+amrex_geom(lvl)%dx(1)*(a-1),&
                !    start_loc(2)+amrex_geom(lvl)%dx(2)*(b-1),&
                !    start_loc(3)+amrex_geom(lvl)%dx(3)*(c-1)
                !write(17) start_loc(x)+amrex_geom(lvl)%dx(x)*(a-1)
!                select case (x)
!                  case(1)
!                    write(17) start_loc(a)+amrex_geom(lvl)%dx(a)*(a-1)
!                  case(2)
!                    write(17) start_loc(b)+amrex_geom(lvl)%dx(b)*(b-1)
!                  case(3)
!                    write(17) start_loc(c)+amrex_geom(lvl)%dx(c)*(c-1)
!                end select
                !write(18) rho_cells(a,b,c,1)


                !  max_boxes = maxval(0:nprocs-1,lvl)
!    if (allocated(x_dims)) deallocate(x_dims)
!    allocate(x_dims(num_blocks(lvl)))
!    if (allocated(y_dims)) deallocate(y_dims)
!    allocate(y_dims(num_blocks(lvl)))
!    if (allocated(z_dims)) deallocate(z_dims)
!    allocate(z_dims(num_blocks(lvl)))
!  do i = 0,nprocs-1
!    if (self == i) then
!      do n = block_starts(i,lvl),block_ends(i,lvl)
!!          write(17) box_dims(1,n),box_dims(2,n),box_dims(3,n)
!!          write(18) box_dims(1,n),box_dims(2,n),box_dims(3,n)
!      end do
!    end if
!      CALL MPI_BARRIER(COMMUNE,IER)
!  end do

!        if (lvl == 0 .and. self == 0) then
!          write(19) (((x_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!            (((y_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!            (((z_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
!!          write(*,*) (((x_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!!            (((y_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3)),&
!!            (((z_points(a,b,c),a=mux%lo(1),mux%hi(1)),b=mux%lo(2),mux%hi(2)),c=mux%lo(3),mux%hi(3))
!        end if

        !write(18) rho_cells(mux%lo(1):mux%hi(1),mux%lo(2):mux%hi(2),mux%lo(3):mux%hi(3),1)
!            end do
!          end do
!        end do
