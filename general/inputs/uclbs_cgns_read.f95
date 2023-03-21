SUBROUTINE uclbs_cgns_read
!
! This reads in a CGNS file to create a new run of UCLBS
!
! Calls: error_out
! Called by: uclbs_main
! Outside calls: cg_open_f,cg_nbases_f,cg_nzones_f,cg_zone_read_f,cg_coord_read_f,
!                  cg_n_section_f, cg_elements_read_f,cg_close_f
!
USE precise
USE cgns
USE constants
USE geom_data
USE nml_inputs
USE nml_inlet_outlet
use amr_processes, only: self,nprocs,commune
use mpi
IMPLICIT NONE
integer :: i,j,l,g,k,p
INTEGER :: index_base,index_zone,index_file,icelldim,&
  iphysdim,ier,index_section,num_sections,nbounds
integer,allocatable :: num_bases(:),num_zones(:),pass_cells(:),cell_base(:),pass_nodes(:),&
  num_cells(:),count_start(:,:)

integer(cgsize_t),allocatable :: read_start(:),read_stop(:),read_group(:)
integer(cgsize_t),allocatable :: cell_start(:),cell_stop(:),cell_group(:),&
  section_type(:),istart(:),iend(:)
INTEGER(cgsize_t),ALLOCATABLE :: elems(:,:),sizes(:,:,:)!,isize(:,:,:)
INTEGER(cgsize_t) :: n_tris,counter,cell_min,&
  vertex_total,cell_total,isize(1,3)  !,idata(2)

INTEGER :: iparent_flag,itype,total_zones,total_bases,total_sections,first_tri,&
  total_section_cells,old_total
integer :: zone_counter,section_counter,tri_start,tri_diff,tri_old,cell_diff
INTEGER :: inlet_num_error_checker(98),outlet_num_error_checker(98),&
  inl_checkers,outl_checkers,inl_value,outl_value,nbase,nzone
integer :: base_disp(0:nprocs-1),base_rec(0:nprocs-1),node_disp(0:nprocs-1),&
  node_rec(0:nprocs-1)

REAL(KIND = dp) :: a1,a2,a3,b1,b2,b3,units
integer,allocatable :: rec_array(:),disp_array(:)

CHARACTER(len=64) :: filename,zonename,section_name
!      LOGICAL :: error_checker!,inlet_num_error_checker(98),outlet_num_error_checker(98)


icelldim = 2
iphysdim = 3
vertex_total = 0
cell_total = 0
cell_min = 1
!istart = 1
nbounds = 0
index_file = 1
index_base = 1
index_zone = 1
index_section = 1
n_tris = 0
cgns_tris = 0
inl_checkers = 0
outl_checkers = 0
total_zones = 0
total_bases = 0
total_sections = 0
zone_counter = 0
section_counter = 0


allocate(rec_array(0:nprocs-1))
allocate(disp_array(0:nprocs-1))
base_rec = 1
do i = 0,nprocs-1
  base_disp(i) = i
end do

inlet_num_error_checker = 0
outlet_num_error_checker = 0

allocate(num_bases(num_cgns_files))
allocate(num_zones(num_cgns_files))
allocate(sizes(1,3,num_cgns_files))
allocate(pass_nodes(0:nprocs-1))

DO g = 1,num_cgns_files

filename = cgns_input_files(g)

!
! Open the pre-existing file
!
  CALL cgp_open_f(filename,CG_MODE_READ,index_file,ier)
!      WRITE(*,*) ier
  IF (ier .NE. CG_OK) CALL cgp_error_exit_f

!
! Get the number of bases and zones
!
  CALL cg_nbases_f(index_file,nbase,ier)
  num_bases(g) = nbase
  total_bases = total_bases + num_bases(g)

!      WRITE(*,*) ier
  IF (ier .NE. CG_OK) CALL cgp_error_exit_f
  do i = 1,num_bases(g)
    CALL cg_nzones_f(index_file,index_base,nzone,ier)
    num_zones(g) = nzone
!      WRITE(*,*) ier
    IF (ier .NE. CG_OK) CALL cgp_error_exit_f
    total_zones = total_zones + num_zones(g)

    do j = 1,num_zones(g)
      CALL cg_nsections_f(index_file,index_base,index_zone,num_sections,ier)
      IF (ier .NE. CG_OK) CALL cgp_error_exit_f
      total_sections = total_sections + num_sections
 ! write(*,*) 'chips',g,i,j
    end do
  end do

  call cgp_close_f(index_file,ier)

end do

allocate(read_start(total_zones))
allocate(read_stop(total_zones))
allocate(read_group(total_zones))
allocate(cell_start(total_sections))
allocate(cell_stop(total_sections))
allocate(cell_group(total_sections))
allocate(section_type(total_sections))
allocate(cell_base(total_sections))
allocate(num_cells(total_sections))
allocate(istart(total_sections))
allocate(iend(total_sections))
!allocate(count_start(0:nprocs-1,total_sections))

index_file = 1
index_base = 1
index_zone = 1
index_section = 1

file_loop: DO g = 1,num_cgns_files
  filename = cgns_input_files(g)
  CALL cgp_open_f(filename,CG_MODE_READ,index_file,ier)
!
!  WRITE(*,*) 'Number of bases ',num_bases
!  WRITE(*,*) 'Number of zones ',num_zones
!  write(*,*) 'Number of sections',num_sections
  DO i = 1,num_bases(g)
    DO j = 1,num_zones(g)

      CALL cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
      zone_counter = zone_counter + 1

      first_tri = isize(1,2) + 1
!      write(*,*) 'isize',isize
!      vertex_total = vertex_total + isize(1,1)
!      cell_total = cell_total + isize(1,2)
!
! Assign the parallel reading locations for vertices
!
      read_group(zone_counter) = isize(1,1)/nprocs

      !do p = 0,nprocs-1
        if (self == 0) then
          read_start(zone_counter) = 1
          read_stop(zone_counter) = read_group(zone_counter)
          !node_disp(p) = 0
        else if (self == nprocs-1) then
          read_start(zone_counter) = self*read_group(zone_counter)+1
          read_stop(zone_counter) = isize(1,1)
          !node_disp(p) = disp_array(p-1) + read_stop(zone_counter) - &
          !  read_start(zone_counter)+1
        else
          read_start(zone_counter) = self*read_group(zone_counter)+1
          read_stop(zone_counter) = self*read_group(zone_counter) +&
            read_group(zone_counter)
          !node_disp(p) = disp_array(p-1) + read_stop(zone_counter) - &
          !  read_start(zone_counter)+1
        end if
        !node_rec(p) = read_stop(zone_counter) -read_start(zone_counter)+1
      !end do

      do p = 0,nprocs-1
        node_disp(p) = p
      end do
      node_rec = 1

      pass_nodes(self) = read_stop(zone_counter) - read_start(zone_counter) + 1

      call mpi_allgatherv(pass_nodes(self),1,MPI_INT,pass_nodes(0:nprocs-1),&
          node_rec(0:nprocs-1),node_disp(0:nprocs-1),MPI_INT,commune,ier)

!      write(*,*) 'boogie',node_disp,node_rec,pass_nodes,read_start,read_stop,self

!          CALL cg_nsection_f(indexf,index_base,index_zone,&
!            num_sections,ier)
!          DO k = 1,num_sections
!          END DO
        section_loop: do k = 1,num_sections

          section_counter = section_counter+1

          CALL cg_section_read_f(index_file,index_base,index_zone,k,section_name,&
          itype,istart(section_counter),iend(section_counter),nbounds,iparent_flag,ier)
!
          section_type(section_counter) = itype

          num_cells(section_counter) = iend(section_counter)-istart(section_counter)+1
!
! Assign the parallel reading locations for cells
!
          if (iend(section_counter)-istart(section_counter)+1 >= nprocs) then
            !write(*,*) 'wibble'
            cell_group(section_counter) = (iend(section_counter)-istart(section_counter)+1)/nprocs
            section_type(section_counter) = itype
            if (self == 0) then
              cell_start(section_counter)= istart(section_counter)
              cell_stop(section_counter) = istart(section_counter)+cell_group(section_counter)-1
            else if (self == nprocs-1) then
              cell_start(section_counter) = istart(section_counter) + self*cell_group(section_counter)
              cell_stop(section_counter) = iend(section_counter)
            else
              cell_start(section_counter) = istart(section_counter) + self*cell_group(section_counter)
              cell_stop(section_counter) = istart(section_counter) + self*cell_group(section_counter)&
                + cell_group(section_counter) - 1
            end if

          else
            !write(*,*) 'wobble'
            if (self < iend(section_counter)-istart(section_counter)) then
              cell_start(section_counter) = self+istart(section_counter)
              cell_stop(section_counter) = cell_start(section_counter)
            else
              cell_start(section_counter) = iend(section_counter)
              cell_stop(section_counter) = iend(section_counter)
            end if
          end if

!          do p = 0,nprocs-1
!            disp_array(p) = 0
!          end do
!          rec_array = 1
!
!
!
!          call mpi_allgatherv(iend(section_counter)-istart(section_counter)+1,1,MPI_INT,count_start(0:nprocs-1,section_counter),&
!          rec_array(0:nprocs-1),disp_array(0:nprocs-1),MPI_INT,commune,ier)

!        section_counter =section_counter +1
          IF (itype /= TRI_3) CYCLE section_loop
!        IF (TRIM(section_name) == 'TetElements' .OR. &
!          TRIM(section_name)=='PyramidElements' .OR. &
!          TRIM(section_name)=='HexElements' .OR. &
!          TRIM(section_name)=='QuadElements') CYCLE section_loop

!            diff = iend(section_counter) - istart(section_counter) + 1
!            WRITE(*,*) 20*diff
          IF (ALLOCATED(elems)) DEALLOCATE(elems)
          ALLOCATE(elems(3,cell_start(section_counter):cell_stop(section_counter)))
          elems = 0

        WRITE(11,*) 'Section ',k,' called ',section_name,' read in.'
        WRITE(11,*) 'Start ',istart(section_counter),' End ',iend(section_counter),' Num bounds ',nbounds

          CALL cgp_elements_read_data_f(index_file,index_base,index_zone,section_counter,&
            cell_start(section_counter),cell_stop(section_counter),&
            elems(1:3,cell_start(section_counter):cell_stop(section_counter)),ier)

!            WRITE(*,*) 'Elements of section ',k,' read in.'

          DO l = cell_start(k),cell_stop(k)

!              WRITE(*,*) elems(1,l),elems(2,l),elems(3,l),iparent_flag,&
!                iparent_data
            IF (elems(1,l) > 0 .and. self <= iend(section_counter)-istart(section_counter)) THEN
!              WRITE(*,*) elems(1,l),elems(2,l),elems(3,l),l,section_name,self
!                WRITE(*,*) elems(1,l)
              n_tris = n_tris + 1
            ELSE
              !IF (ALLOCATED(elems)) DEALLOCATE(elems)
              CYCLE section_loop
            END IF
          END DO
        !IF (ALLOCATED(elems)) DEALLOCATE(elems)

      end do section_loop !section/element loop
    END DO !zone loop
  END DO !base loop
  call cgp_close_f(index_file,ier)
end do file_loop!file loop

!
! Loop through all the bases and zones
!
! p!file loop

!        IF (ALLOCATED(elems)) DEALLOCATE(elems)
!  read_group = vertex_total/nprocs
!  cell_group = cell_total/nprocs
if (self == 0) cell_base=cell_start
call mpi_barrier(commune,ier)
call mpi_bcast(cell_base,total_sections,MPI_INT,0,commune,ier)
call mpi_barrier(commune,ier)
zone_counter = 0
section_counter = 0

!write(*,*) 'whiffle ball',cell_base

CALL MPI_Allreduce(n_tris,cgns_tris,1,MPI_INTEGER, &
  MPI_SUM,commune,ier)

total_tris = cgns_tris

index_base = 1
index_zone = 1
index_section = 1

WRITE(*,101) cgns_tris
ALLOCATE(tri_vertices(3,3,cgns_tris))
ALLOCATE(normal(3,cgns_tris))
ALLOCATE(tri_meaning(cgns_tris))
allocate(pass_cells(0:nprocs-1))
!      cgns_normals = 0.0
!      cgns_vertices = 0.0
!      tri_meaning = 0
counter = 0
tri_start = 1
zone_counter = 0
section_counter = 0
index_file = 1
index_base = 1
index_zone = 1
index_section = 1
total_section_cells = 0
old_total = 0

file_loop_2: DO g = 1,num_cgns_files
  filename = cgns_input_files(g)

!  write(*,*) 'cluck?'

  CALL cgp_open_f(filename,CG_MODE_READ,index_file,ier)

  IF (cgns_units(g) == 1) THEN
    WRITE(11,*) 'CGNS file units in meters, no adjustment needed.'
    WRITE(11,*) ''
    units = 1.0D0
  ELSE IF (cgns_units(g) == 2) THEN
    WRITE(11,*) 'Units in millimeters, need to divide vertex array ',&
      'by 1000.'
    units = 1000.0D0
  ELSE IF (cgns_units(g) == 3) THEN
    WRITE(11,*) 'Units in centimeters, need to divide vertex array ',&
      'by 100.'
    units = 100.0D0
  ELSE IF (cgns_units(g) == 4) THEN
    WRITE(11,*) 'Units in inches, need to divide vertex array ',&
      'by 39.37'
    units = 39.37D0
  ELSE IF (cgns_units(g) == 5) THEN
    WRITE(11,*) 'Units in feet, need to divide vertex array ',&
      'by 3.28'
    units = 3.28D0
  else
    units = 1.0D0
  END IF

!  write(*,*) 'peep?'

!  CALL cg_nbases_f(g,num_bases,ier)
!  WRITE(*,*) ier
!  IF (ier .NE. CG_OK) CALL cgp_error_exit_f
!
!  CALL cg_nzones_f(g,i,num_zones,ier)
!  WRITE(*,*) ier
!  IF (ier .NE. CG_OK) CALL cgp_error_exit_f

  DO i = 1,num_bases(g)
    DO j = 1,num_zones(g)

!          ALLOCATE(tri_vertices(3,3,cgns_tris))
!      write(*,*) 'meow?',g,i,j,isize
!      CALL cg_zone_read_f(g,i,j,zonename,isize,ier)
      zone_counter = zone_counter + 1
      ALLOCATE(nodes_x(isize(1,1)))
      ALLOCATE(nodes_y(isize(1,1)))
      ALLOCATE(nodes_z(isize(1,1)))
!      write(*,*) 'woof?'
!      CALL cgp_coord_read_data_f(index_file,index_base,index_zone, INT(1,cgsize_t),&
!        read_start(zone_counter),read_stop(zone_counter),&
!        nodes_x(read_start(zone_counter):read_stop(zone_counter)),ier)
!
!      CALL cgp_coord_read_data_f(index_file,index_base,index_zone,INT(2,cgsize_t),&
!        read_start(zone_counter),read_stop(zone_counter), &
!        nodes_y(read_start(zone_counter):read_stop(zone_counter)),ier)
!
!      CALL cgp_coord_read_data_f(index_file,index_base,index_zone,INT(3,cgsize_t), &
!        read_start(zone_counter),read_stop(zone_counter), &
!        nodes_z(read_start(zone_counter):read_stop(zone_counter)),ier)

      CALL cgp_coord_read_data_f(index_file,index_base,index_zone, INT(1,cgsize_t),&
        int(1,cgsize_t),isize(1,1),&
        nodes_x(int(1,cgsize_t):isize(1,1)),ier)

      CALL cgp_coord_read_data_f(index_file,index_base,index_zone,INT(2,cgsize_t),&
        int(1,cgsize_t),isize(1,1), &
        nodes_y(int(1,cgsize_t):isize(1,1)),ier)

      CALL cgp_coord_read_data_f(index_file,index_base,index_zone,INT(3,cgsize_t), &
        int(1,cgsize_t),isize(1,1), &
        nodes_z(int(1,cgsize_t):isize(1,1)),ier)
!
!
!
!      do p = 0,nprocs-1
!        if (p == 0) then
!          node_disp(p) = 0
!!        else if () then
!!
!        else
!          node_disp(p) = node_disp(p-1) + pass_nodes(p-1)
!        end if
!      end do

!      write(*,*) 'output tests',read_start(zone_counter),read_stop(zone_counter),&
!        pass_nodes,node_disp

!      call mpi_allgatherv(nodes_x(read_start(zone_counter):read_stop(zone_counter)),&
!        pass_nodes(self),MPI_DOUBLE_PRECISION,nodes_x,pass_nodes(0:nprocs-1),&
!        node_disp(0:nprocs-1),MPI_DOUBLE_PRECISION,commune,ier)
!      call mpi_allgatherv(nodes_y(read_start(zone_counter):read_stop(zone_counter)),&
!        pass_nodes(self),MPI_DOUBLE_PRECISION,nodes_y,pass_nodes(0:nprocs-1),&
!        node_disp(0:nprocs-1),MPI_DOUBLE_PRECISION,commune,ier)
!      call mpi_allgatherv(nodes_z(read_start(zone_counter):read_stop(zone_counter)),&
!        pass_nodes(self),MPI_DOUBLE_PRECISION,nodes_z,pass_nodes(0:nprocs-1),&
!        node_disp(0:nprocs-1),MPI_DOUBLE_PRECISION,commune,ier)
!      do p = 1,1293!isize(1,1)
!      write(*,*) 'nodes',nodes_x(p),nodes_y(p),nodes_z(p),p,self
!      end do
!
!
!write(*,*) 'moo?'
      section_loop_2: DO k = 1,num_sections
!      write(*,*) 'AWWWOOOOO',k
        section_counter = section_counter + 1

!      write(*,*) 'section stuff',section_counter,istart(section_counter),&
!        iend(section_counter),self
      CALL cg_section_read_f(index_file,index_base,index_zone,section_counter,section_name,&
          itype,istart(section_counter),iend(section_counter),nbounds,iparent_flag,ier)
!      write(*,*) 'bad bad',ier,section_counter,itype,section_name,self
!        CALL cg_section_read_f(indexf,i,j, &
!          k,section_name,itype,istart(section_counter),iend(section_counter),nbounds,&
!          iparent_flag,ier)

        IF (section_type(section_counter) /= TRI_3 .or. section_type(section_counter) /= 5) CYCLE section_loop_2
!            diff = iend(section_counter) - istart(section_counter) + 1
!            WRITE(*,*) 20*diff
        if (allocated(elems)) deallocate(elems)
        ALLOCATE(elems(3,istart(section_counter):iend(section_counter)))
!            ALLOCATE(normal(3,isize(1,2)))
            elems = 0
!        write(*,*) 'sneeple dook',section_counter, istart(section_counter),&
!          iend(section_counter),section_type(section_counter),self
        CALL cgp_elements_read_data_f(index_file,index_base,index_zone,section_counter,&
            istart(section_counter),iend(section_counter),&
            elems(1:3,istart(section_counter):iend(section_counter)),ier)
        IF (ier .NE. CG_OK) CALL cgp_error_exit_f
!        CALL cgp_elements_read_data_f(index_file,index_base,index_zone,section_counter,&
!            cell_start(section_counter),cell_stop(section_counter),&
!            elems(1:3,cell_start(section_counter):cell_stop(section_counter)),ier)
!        CALL cgp_elements_read_data_f(index_file,index_base,index_zone, &
!          section_counter,cell_start(section_counter),cell_stop(section_counter),&
!          elems(1:3,cell_start(section_counter):cell_stop(section_counter)),ier)

        counter = istart(section_counter) - cell_base(2)
!        counter = counter + (cell_start(section_counter)-cell_base(section_counter))
!        counter = cell_start(section_counter)

!        write(*,*) 'baaa',section_counter,cell_start(section_counter),cell_stop(section_counter),&
!          section_type(section_counter),cell_base(section_counter),section_counter,&
!          counter,old_total,ier,self

        DO l = istart(section_counter),iend(section_counter)
!          write(*,*) 'neigh',l,elems(1,l),elems(2,l),elems(3,l)

          IF (elems(1,l) >= 0) THEN
            counter = counter +1
            !WRITE(*,*) elems(1,l),elems(2,l),elems(3,l)
            tri_vertices(1,1,counter) = nodes_x(elems(1,l))/units
            tri_vertices(2,1,counter) = nodes_y(elems(1,l))/units
            tri_vertices(3,1,counter) = nodes_z(elems(1,l))/units

            tri_vertices(1,2,counter) = nodes_x(elems(2,l))/units
            tri_vertices(2,2,counter) = nodes_y(elems(2,l))/units
            tri_vertices(3,2,counter) = nodes_z(elems(2,l))/units

            tri_vertices(1,3,counter) = nodes_x(elems(3,l))/units
            tri_vertices(2,3,counter) = nodes_y(elems(3,l))/units
            tri_vertices(3,3,counter) = nodes_z(elems(3,l))/units

            !write(*,*) 'tri verts',tri_vertices(:,:,counter),counter,section_counter,l,self

            IF (section_name(1:9) == 'visc_wall') THEN
              tri_meaning(counter) = -1001
            ELSE IF (section_name(1:5) == 'inlet' .OR. &
              section_name(1:5) == 'Inlet' .OR. &
              section_name(1:5) == 'INLET') THEN

              READ(section_name(6:7),'(I2)') inl_value

              inlet_num_error_checker(inl_value) = &
                inl_value

              tri_meaning(counter) = -100-inl_value

!                  WRITE(*,*) section_name,tri_meaning(counter),inl_value

            ELSE IF (section_name(1:6) == 'outlet' .OR. &
              section_name(1:6) == 'Outlet' .OR. &
              section_name(1:6) == 'OUTLET') THEN

              READ(section_name(7:8),'(I2)') outl_value

              outlet_num_error_checker(outl_value) = &
                outl_value

              tri_meaning(counter) = -200-outl_value

!                  WRITE(*,*) section_name,tri_meaning(counter),outl_value

            ELSE IF (section_name(1:10) == 'freestream' .or. &
              section_name(1:10) == 'Freestream' .or. &
              section_name(1:10) == 'FREESTREAM') THEN
              tri_meaning(counter) = -1000
            ELSE IF (section_name == 'inviscid_wall') THEN
              tri_meaning(counter) = -11
            END IF
!
          ELSE
            !IF (ALLOCATED(elems)) DEALLOCATE(elems)
!                counter_start = l
            !exit
            !write(*,*) 'Marcus Aurelius'
            counter = counter + 1
            cycle
!            CYCLE section_loop_2
          END IF
        END DO
!write(*,*) 'quack?'
        if (allocated(elems)) DEALLOCATE(elems)
!        tri_old = tri_start
!        tri_diff = cell_stop(section_counter)-cell_start(section_counter)+1
!        CALL MPI_Allreduce(tri_diff,tri_start,1,MPI_INTEGER, &
!          MPI_SUM,commune,ier)
!
! Pass how many elements each processor will be passing, then solve for the diplacements
!
!        do p = 0,nprocs-1
!!          if (self == 0) then
!!            disp_array(p) = 0
!!          else
!            disp_array(p) = p
!!          end if
!
!        end do
!        rec_array = 1
!
!        if (iend(section_counter)-istart(section_counter)+1 >= nprocs) then
!          pass_cells(self) = cell_stop(section_counter)-cell_start(section_counter)+1
!        else
!          if (self <= iend(section_counter) - istart(section_counter)) then
!            pass_cells(self) = 1
!          else
!            pass_cells(self) = 0
!          end if
!        end if
!!
!        call mpi_allgatherv(pass_cells(self),1,MPI_INT,pass_cells(0:nprocs-1),&
!          rec_array(0:nprocs-1),disp_array(0:nprocs-1),MPI_INT,commune,ier)



!        rec_array(self) = cell_stop(section_counter)-cell_start(section_counter)+1
!        write(*,*) 'squeak',pass_cells,section_counter,istart(section_counter),&
!          iend(section_counter),self
!!        write(*,*) 'hoot?',rec_array(self),pass_cells(self),&
!!          cell_stop(section_counter),cell_start(section_counter),&
!!          tri_old,tri_start,counter,self
!!        call mpi_allgather(rec_array(self),1,MPI_INT,rec_array(0:nprocs-1),nprocs,&
!!          MPI_INT,commune,ier)
!        do p = 0,nprocs-1
!          if (sum(pass_cells) >=nprocs) then
!            if (p == 0) then
!              disp_array(p) = 0
!            else
!              disp_array(p) = disp_array(p-1) + pass_cells(p-1)
!            end if
!          else
!             if (p == 0) then
!              disp_array(p) = 0
!            else if (p < sum(pass_cells)) then
!              disp_array(p) = disp_array(p-1) + pass_cells(p-1)
!            else
!              disp_array(p) = disp_array(p-1)
!            end if
!          end if
!!          disp_array(p) = cell_base(section_counter)
!!          else
!!            disp_array(p) = rec_array(p) + disp_array(p-1)
!!          end if
!        end do
!        old_total = total_section_cells
!        total_section_cells = total_section_cells + sum(pass_cells)
        !call mpi_allgather(disp_array(self),1,MPI_INT,disp_array(0:nprocs-1),nprocs,&
        !  MPI_INT,commune,ier)
!        call mpi_allgather(rec_array(self),1,MPI_INT,rec_array(0:nprocs-1),nprocs,&
!          MPI_INT,commune,ier)

!
!write(*,*) 'ribbit?',disp_array,cell_diff,self
!write(*,*) 'croak?', pass_cells,self
!write(*,*) 'chirrup?',istart,iend,self
!write(*,*) 'honk?',cell_start,cell_stop,total_section_cells,old_total,self
!        cell_diff = cell_stop(section_counter)-cell_start(section_counter) + 1
!        call mpi_allgatherv(tri_meaning(old_total+disp_array(self)+1:old_total+disp_array(self)+pass_cells(self)), &
!          pass_cells(self),MPI_INT,tri_meaning(old_total+1:total_section_cells),&
!          pass_cells,disp_array,MPI_INT,commune,ier)
!
!        call mpi_barrier(commune,ier)
!
!write(*,*) 'hop?'
!        cell_diff = cell_diff*3
!        disp_array = disp_array*3
!        rec_array(self) = cell_stop(section_counter)-cell_start(section_counter)+1
!        rec_array = rec_array*3
!        pass_cells = pass_cells*3
!write(*,*) 'caw?'

!
!        write(*,*) 'chitter',old_total+disp_array(self)+1,old_total+disp_array(self)+pass_cells(self)
!        cell_diff = cell_diff*3
!        disp_array = disp_array*3
!        rec_array=rec_array*3
!        pass_cells = pass_cells*3
!        call MPI_allgatherv(tri_vertices(1:3,1:3,old_total+disp_array(self)+1:old_total+disp_array(self)+pass_cells(self)),&
!          9*pass_cells(self),MPI_DOUBLE_PRECISION,&
!          tri_vertices(1:3,1:3,old_total+1:total_section_cells),&
!          9*pass_cells,9*disp_array,MPI_DOUBLE_PRECISION,commune,ier)
!
!        call mpi_barrier(commune,ier)

        !do p = old_total+disp_array(self)+1, old_total+disp_array(self)+pass_cells(self)


!        call MPI_allgatherv(normal(1:3,old_total+disp_array(self)+1:old_total+disp_array(self)+pass_cells(self)),&
!          3*pass_cells(self),MPI_DOUBLE_PRECISION,&
!          normal(1:3,old_total+1:total_section_cells),&
!          3*pass_cells,3*disp_array,MPI_DOUBLE_PRECISION,commune,ier)

!        counter = tri_start
      END DO section_loop_2
!      do p = 1,1478
!        write(*,*) 'tri check',tri_vertices(1:3,1,p),tri_vertices(1:3,2,p),&
!          tri_vertices(1:3,3,p),tri_meaning(p),p,self
!
!      end do

      do p = 1,cgns_tris
        !do p = old_total+disp_array(0)+1,old_total+disp_array(nprocs-1)+pass_cells(nprocs-1)
!
! a vector is cgns_vertices 2 - cgns_vertices 1
! b vector is cgns_vertices 3 - cgns_vertices 1
!
!        Cross Produce to find the cgns_normals direction
!       |i  j  k |
!       |a1 a2 a3|
!       |b1 b2 b3|
!
! N = (a2*b3-a3*b2)i+(a3*b1-a1*b3)j+(a1*b2-b1*a2)k
!
            a1 = tri_vertices(1,2,p)-tri_vertices(1,1,p)
            a2 = tri_vertices(2,2,p)-tri_vertices(2,1,p)
            a3 = tri_vertices(3,2,p)-tri_vertices(3,1,p)

            b1 = tri_vertices(1,3,p)-tri_vertices(1,1,p)
            b2 = tri_vertices(2,3,p)-tri_vertices(2,1,p)
            b3 = tri_vertices(3,3,p)-tri_vertices(3,1,p)

            normal(1,p) = (a2*b3-a3*b2)*1000D0
            normal(2,p) = (a3*b1-a1*b3)*1000D0
            normal(3,p) = (a1*b2-b1*a2)*1000D0
!            write(*,*) 'normal test',a1,a2,a3,b1,b2,b3,normal(1,p),normal(2,p),normal(3,p),&
!              tri_vertices(1:3,1,p),tri_vertices(1:3,2,p),tri_vertices(1:3,3,p),tri_meaning(p),p,self
          if (abs(normal(1,p)) < plus_epsilon .and. abs(normal(2,p)) < plus_epsilon .and. &
               abs(normal(3,p)) < plus_epsilon ) then
            write(*,*) 'DANGER WILL ROBINSON, DANGER',self,p,tri_vertices(1:3,1,p),&
              tri_vertices(1:3,2,p),tri_vertices(1:3,3,p)
          end if
        end do

      DEALLOCATE(nodes_x)
      DEALLOCATE(nodes_y)
      DEALLOCATE(nodes_z)

    END DO
  END DO
call cgp_close_f(index_file,ier)
END DO file_loop_2

!do p = 1,1478
!write(*,*) 'normal test',normal(1,p),normal(2,p),normal(3,p),p,self
!end do


!do k = 1,total_sections
!  disp_array() =
!  rec_array() =
!  call mpi_allgatherv(  ,&
!    rec_array,disp_array,MPI_INT,commune,ier)
!  call MPI_allgatherv(tri_vertices(1:3,1:3, ),  tri_vertices,&
!    rec_array,disp_array,MPI_DOUBLE_PRECISION,commune,ier)
!  call MPI_allgatherv(,normal(1:3,),MPI_DOUBLE_PRECISION,normal,&
!    rec_array,disp_array,MPI_DOUBLE_PRECISION,commune,ier)
!end do
!write(*,*) 'google'
DO i = 1,98
  IF (inlet_num_error_checker(i) > num_inlets ) THEN
    WRITE(11,*) 'Number of inlets from cgns file exceeds number in',&
      ' uclbs.in. These must match or information will be lost.'
    CALL error_out
  END IF
  IF (outlet_num_error_checker(i) > num_outlets ) THEN
    WRITE(11,*) 'Number of outlets from cgns file exceeds number in',&
      ' uclbs.in. These must match or information will be lost.'
    CALL error_out
  END IF
  IF (inlet_num_error_checker(i) > 0) THEN
    inl_checkers = inl_checkers+1
  END IF
  IF (outlet_num_error_checker(i) > 0) THEN
    outl_checkers = outl_checkers + 1
  END IF
END DO
!write(*,*) 'dandy'
!IF (inl_checkers < num_inlets) THEN
!  WRITE(11,*) 'There are too few inlets specified in the cgns file or',&
!    ' too many in uclbs.in'
!  CALL error_out
!END IF
!IF (outl_checkers < num_outlets) THEN
!  WRITE(11,*) 'There are too few outlets specified in the cgns file or',&
!    ' too many in uclbs.in'
!  CALL error_out
!END IF
!write(*,*) 'sandalwood'

!      WRITE(*,*) tri_meaning

 101  FORMAT ('Total number of tris from CGNS files ',I8)

END SUBROUTINE



      !      WRITE(*,*) isize

!      CALL cg_goto_f(indexf,index_base,'end')
!
! Go through all the bases and zones to get
!
!  DO i = 1,num_bases
!    DO j = 1,num_zones
!
!      CALL cg_zone_read_f(g,i,j,zonename,isize,ier)
!      IF (ier .NE. CG_OK) CALL cgp_error_exit_f
!      WRITE(11,*) 'Number of cells in zone ',isize(1,1)
!
!      zone_counter = zone_counter + 1
!
!      ALLOCATE(nodes_x(isize(1,1)))
!      ALLOCATE(nodes_y(isize(1,1)))
!      ALLOCATE(nodes_z(isize(1,1)))
!!
!! Pull coordinates from the file based on their locations
!!
!      CALL cgp_coord_read_f(g,i,j, INT(1,cgsize_t),&
!        read_start(zone_counter),read_stop(zone_counter),&
!        nodes_x(read_start(zone_counter):read_stop(zone_counter)),ier)
!
!      CALL cgp_coord_read_f(g,i,j,INT(2,cgsize_t),&
!        read_start(zone_counter),read_stop(zone_counter), &
!        nodes_y(read_start(zone_counter):read_stop(zone_counter)),ier)
!
!      CALL cgp_coord_read_f(g,i,j,INT(3,cgsize_t), &
!        read_start(zone_counter),read_stop(zone_counter), &
!        nodes_z(read_start(zone_counter):read_stop(zone_counter)),ier)
!
!      IF (ier .NE. CG_OK) CALL cgp_error_exit_f
!
!      WRITE(11,*) 'Coordinates read in.'
!
!      CALL cg_nsections_f(g,i,j, num_sections,ier)
!      IF (ier .NE. CG_OK) CALL cgp_error_exit_f
!
!      WRITE(11,*) 'Number of sections ',num_sections
!
!!          ALLOCATE(elems(3,isize(1,2)))
!
!      section_loop: DO k = 1,num_sections
!
!!            ALLOCATE(elems(3,isize(1,1)))
!
!!        CALL cg_section_read_f(indexf,index_base,index_zone, &
!!          k,section_name,itype,istart(section_counter),iend(section_counter),nbounds,&
!!          iparent_flag,ier)
!!        IF (ier .NE. CG_OK) CALL cgp_error_exit_f
!        section_counter =section_counter +1
!        IF (itype /= TRI_3) CYCLE section_loop
!!        IF (TRIM(section_name) == 'TetElements' .OR. &
!!          TRIM(section_name)=='PyramidElements' .OR. &
!!          TRIM(section_name)=='HexElements' .OR. &
!!          TRIM(section_name)=='QuadElements') CYCLE section_loop
!
!!            diff = iend(section_counter) - istart(section_counter) + 1
!!            WRITE(*,*) 20*diff
!        ALLOCATE(elems(3,isize(1,2)))
!        elems = 0
!
!
!        WRITE(11,*) 'Section ',k,' called ',section_name,' read in.'
!        WRITE(11,*) 'Start ',istart(section_counter),' End ',iend(section_counter),' Num bounds ',&
!          nbounds
!
!        CALL cgp_elements_read_data_f(g,i,j,k,&
!          cell_start(section_counter),cell_stop(section_counter),&
!          elems(cell_start(section_counter):cell_stop(section_counter)),ier)
!
!!            WRITE(*,*) 'Elements of section ',k,' read in.'
!
!        DO l = cell_start,cell_stop
!
!!              WRITE(*,*) elems(1,l),elems(2,l),elems(3,l),iparent_flag,&
!!                iparent_data
!          IF (elems(1,l) > 0) THEN
!!                WRITE(*,*) elems(1,l)
!            n_tris = n_tris + 1
!          ELSE
!            IF (ALLOCATED(elems)) DEALLOCATE(elems)
!            CYCLE section_loop
!          END IF
!        END DO
!        IF (ALLOCATED(elems)) DEALLOCATE(elems)
!
!      END DO section_loop
!
!      deallocate(nodes_x)
!      deallocate(nodes_y)
!      deallocate(nodes_z)
!
!    END DO !zone loop
!  END DO !base loop
!  cgp_close_f(indexf,ier)
!END DO file_loop

! DO i = 1,num_bases(g)
!    DO j = 1,num_zones(g)
!
!      CALL cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
!      zone_counter = zone_counter + 1
!!      vertex_total = vertex_total + isize(1,1)
!!      cell_total = cell_total + isize(1,2)
!!
!! Assign the parallel reading locations for vertices
!!
!write(*,*) 'blub?',isize
!      read_group(zone_counter) = isize(1,1)/nprocs
!      if (self == 0) then
!        read_start(zone_counter) = 1
!        read_stop(zone_counter) = read_group(zone_counter)
!      else if (self == nprocs-1) then
!        read_start(zone_counter) = self*read_group(zone_counter)+1
!        read_stop(zone_counter) = isize(1,1)
!      else
!        read_start(zone_counter) = self*read_group(zone_counter)+1
!        read_stop(zone_counter) = self*read_group(zone_counter) +&
!          read_group(zone_counter)
!      end if
!write(*,*) 'spiffle',isize,zone_counter,zonename,index_file,index_base,index_zone
!!          CALL cg_nsection_f(indexf,index_base,index_zone,&
!!            num_sections,ier)
!!          DO k = 1,num_sections
!!          END DO
!        section_loop: do k = 1,num_sections
!          CALL cg_section_read_f(index_file,index_base,index_zone,k,section_name,&
!          itype,istart(section_counter),iend(section_counter),nbounds,iparent_flag,ier)
!!
!          section_counter = section_counter+1
!          write(*,*) 'salsa',index_file,index_base,index_zone,k,section_name,&
!          itype,istart(section_counter),iend(section_counter),nbounds,iparent_flag,ier,self
!!
!! Assign the parallel reading locations for cells
!!
!          cell_group(section_counter) = (iend(section_counter)-istart(section_counter))/nprocs
!          section_type(section_counter) = itype
!          if (self == 0) then
!            cell_start(section_counter)= istart(section_counter)
!            cell_stop(section_counter) = istart(section_counter)+cell_group(section_counter)
!          else if (self == nprocs-1) then
!            cell_start(section_counter) = istart(section_counter) + self*cell_group(section_counter) +1
!            cell_stop(section_counter) = iend(section_counter)
!          else
!            cell_start(section_counter) = istart(section_counter) + self*cell_group(section_counter) +1
!            cell_stop(section_counter) = istart(section_counter) + self*cell_group(section_counter)&
!              + cell_group(section_counter)
!          end if
!
!!        section_counter =section_counter +1
!          IF (itype /= TRI_3 .or. itype /= 5) cycle section_loop
!!            write(*,*) 'skipping junk'
!!          else!CYCLE section_loop
!          write(*,*) 'cell stuff',cell_start(section_counter),cell_stop(section_counter),self
!!        IF (TRIM(section_name) == 'TetElements' .OR. &
!!          TRIM(section_name)=='PyramidElements' .OR. &
!!          TRIM(section_name)=='HexElements' .OR. &
!!          TRIM(section_name)=='QuadElements') CYCLE section_loop
!
!!            diff = iend(section_counter) - istart(section_counter) + 1
!!            WRITE(*,*) 20*diff
!          IF (ALLOCATED(elems)) DEALLOCATE(elems)
!          ALLOCATE(elems(3,cell_start(section_counter):cell_stop(section_counter)))
!          elems = 0
!
!        WRITE(*,*) 'Section ',k,' called ',section_name,' read in.',self
!        WRITE(*,*) 'Start ',istart(section_counter),' End ',iend(section_counter),' Num bounds ',nbounds,self
!          index_section = k
!          CALL cgp_elements_read_data_f(index_file,index_base,index_zone,k,&
!            cell_start(section_counter),cell_stop(section_counter),&
!            elems(1:3,cell_start(section_counter):cell_stop(section_counter)),ier)
!!WRITE(*,*) elems(1,l),elems(2,l),elems(3,l),iparent_flag,&
!!                iparent_data
!            WRITE(*,*) 'Elements of section ',k,' read in.',cell_start(k),cell_stop(k),&
!              k,index_section,section_counter,ier,self
!            !write(*,*) 'fuck it',elems,self
!          DO l = cell_start(k),cell_stop(k)
!
!            WRITE(*,*) elems(1,l),elems(2,l),elems(3,l),l,section_name,self
!            !IF (elems(1,l) > 0) THEN
!            !    WRITE(*,*) elems(1,l),elems(2,l),elems(3,l)
!              n_tris = n_tris + 1
!!            ELSE
!!              IF (ALLOCATED(elems)) DEALLOCATE(elems)
!              !CYCLE section_loop
!            !END IF
!          END DO
!
!        cgns_tris = cgns_tris + n_tris
!!        end if
!
!!        call mpi_barrier(commune,ier)
!      end do section_loop !section/element loop
!    END DO !zone loop
!  END DO !base loop
!
!  sizes(1:1,1:3,g) = isize(1:1,1:3)
!
!  call cgp_close_f(index_file,ier)
!
!end do file_loo

!                WRITE(*,*) counter
!                WRITE(*,*) elems(1,l),elems(2,l),&
!                  elems(3,l)
!                WRITE(*,*) tri_vertices(1,1,counter),tri_vertices(2,1,counter),&
!                  tri_vertices(3,1,counter)
!                WRITE(*,*) tri_vertices(1,2,counter),tri_vertices(2,2,counter),&
!                  tri_vertices(3,2,counter)
!                WRITE(*,*) tri_vertices(1,3,counter),tri_vertices(2,3,counter),&
!                  tri_vertices(3,3,counter)




! a vector is cgns_vertices 2 - cgns_vertices 1
! b vector is cgns_vertices 3 - cgns_vertices 1
!
!        Cross Produce to find the cgns_normals direction
!       |i  j  k |
!       |a1 a2 a3|
!       |b1 b2 b3|
!
! N = (a2*b3-a3*b2)i+(a3*b1-a1*b3)j+(a1*b2-b1*a2)k
!
!
!
!            a1 = tri_vertices(1,2,counter)-tri_vertices(1,1,counter)
!            a2 = tri_vertices(2,2,counter)-tri_vertices(2,1,counter)
!            a3 = tri_vertices(3,2,counter)-tri_vertices(3,1,counter)
!
!            b1 = tri_vertices(1,3,counter)-tri_vertices(1,1,counter)
!            b2 = tri_vertices(2,3,counter)-tri_vertices(2,1,counter)
!            b3 = tri_vertices(3,3,counter)-tri_vertices(3,1,counter)
!
!            normal(1,counter) = (a2*b3-a3*b2)*1000D0
!            normal(2,counter) = (a3*b1-a1*b3)*1000D0
!            normal(3,counter) = (a1*b2-b1*a2)*1000D0
!                WRITE(*,*) normal(1,counter),normal(2,counter),&
!                  normal(3,counter)
!            WRITE(11,*) 'Normal direction ',normal(1,counter),normal(2,counter),&
!              normal(3,counter)


