module linkwise
!
! Stores all link and path information
!
!
!
!
USE precise
!use grid_data,only: rdx
IMPLICIT NONE
SAVE
private
public :: crossing_point_and_dist
public :: cx,cy,cz,rcx,rcy,rcz,link_dist,opp,cmag
public :: cxxx,cxxy,cxxz,cxyy,cxyz,cxzz,cyyy,czzz,cyyz,cyzz
public :: cxc2,cyc2,czc2,cx2,cy2,cz2,cxcy,cxcz,cycz
public :: yi,kroen_val,p_coeff,p_kroen,shifted,p_kr_val
public :: abs_latt_leng,use_shifted
public :: cx_27,cy_27,cz_27,opp_27,erp

integer,allocatable :: cx(:),cy(:),cz(:),opp(:),&
  periodic_bc_links(:,:)
REAL(KIND=dp), ALLOCATABLE :: cmag(:),rcx(:),rcy(:),rcz(:),link_dist(:,:),abs_latt_leng(:)
integer,allocatable :: cxxx(:),cxxy(:),cxxz(:),cxyy(:),cxzz(:),cxyz(:),cyyy(:),czzz(:),&
  cyyz(:),cyzz(:)
integer,allocatable :: cxc2(:),cyc2(:),czc2(:),cx2(:),cy2(:),cz2(:),cxcy(:),cxcz(:),cycz(:)
integer,allocatable :: yi(:,:),kroen_val(:,:),p_coeff(:,:),p_kroen(:),p_kr_val(:)
logical :: shifted,use_shifted

integer :: cx_27(1:27) = (/0, 1, 0, 0,-1, 0, 0, 1, 0,-1, 0, 1, 1,-1,-1, 1, 0,-1, 0,&
           1, 1, 1, 1,-1,-1,-1,-1/)
integer :: cy_27(1:27) = (/0, 0, 1, 0, 0,-1, 0, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,&
           1, 1,-1,-1, 1, 1,-1,-1/)
integer :: cz_27(1:27) = (/0, 0, 0, 1, 0, 0,-1, 0,-1, 0, 1, 1,-1,-1, 1, 0,-1, 0, 1,&
           1,-1, 1,-1, 1,-1, 1,-1/)
integer :: opp_27(1:27) = (/1, 5, 6, 7, 2, 3, 4, 18, 19, 16, 17, 14, 15, 12, 13,&
      10, 11, 8, 9,27,26,25,24,23,22,21,20/)

integer :: erp(1:39) = (/1,5, 6, 7, 2, 3, 4,14,15,12,13,10,11, 8, 9,19, 20,&
        21,16,17,18,26,29,32,27,22,25,30,23,28,33,24,31,37,38,39,&
        34,35,36/)

contains

integer function crossing_point_and_dist(ray_crossing,ray_end1,i,loc_dx,lvl) result (wall_dist)
!
! Called by: inter_bb_setup
!
  REAL(KIND = dp) :: ray_crossing(3),ray_end1(3),ray_intersect_length,loc_dx,dummy
  integer :: i,lvl!wall_dist,

  ray_intersect_length = SQRT((ray_end1(1)-ray_crossing(1))**2 &
         +(ray_end1(2)-ray_crossing(2))**2 + (ray_end1(3)-ray_crossing(3))**2)

  !dummy = ray_intersect_length/loc_dx
  wall_dist = -1001-FLOOR(ray_intersect_length/(link_dist(i,lvl))*1E9)!*1000000
end function

END MODULE
!  if (wall_dist == -1001) then
!    write(*,*) 'intersect length',ray_intersect_length,' link distance ',link_dist(i,lvl),&
!    ' and state number',wall_dist,'on level',lvl,'for direction',i
!  end if
!  write(*,*) 'intersect length',ray_intersect_length,' link distance ',link_dist(i,lvl),&
!    ' and state number',wall_dist,'on level',lvl,'for direction',i
!  write(*,*) 'details, ray end 1',ray_end1,' ray crossing ',ray_crossing,'direction',i,&
!    'local dx',link_dist(i,lvl)
!  wall_dist = -1001-floor(dummy)*100000
!  wall_dist = -1001-FLOOR(ray_intersect_length/(link_dist(i,lvl)*loc_dx)*1E9)!*1000000

!real(kind=dp),allocatable :: link_grid(:,:),link_dist(:)
!LOGICAL,ALLOCATABLE :: boundary_node(:)

!INTEGER :: link_number
!INTEGER,ALLOCATABLE :: link(:,:)
