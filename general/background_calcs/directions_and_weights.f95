subroutine directions_and_weights
!
! Define all the necessary directions and weightings for what we're
! doing.
!
! With ELBM, the weights are included in fieq as a result of the pressure tensor.
!
! Called by: background_calcs
! Calls: error_out
!uclbs_main
use linkwise
!use nml_bounding_box
use grid_data
!use loop_data
use constants
use precise
!use amr_processes, only :amr_max_lvl
use amr_info_holder, only: nghosts,nghosts_state,nghosts_mid
IMPLICIT NONE
!      INTEGER :: dimensions
real(kind=dp) :: center_2d,card_2d,diag_2d,center_3d,card_3d
real(kind=dp) :: center_fast_2d,card_close_fast_2d
real(kind=dp) :: diag_close_fast_2d,card_medium_fast_2d
real(kind=dp) :: card_far_fast_2d,diag_far_fast_2d,diag_close_fast_3d
real(kind=dp) :: center_fast_3d, card_close_fast_3d
real(kind=dp) :: card_medium_fast_3d,diag_3d
real(kind=dp) :: card_far_fast_3d, diag_far_fast_3d
INTEGER :: i,j,lvl

!write(*,*) 'direction assignment'

IF (dimensions == 2) THEN
  IF (mach <= 0.4D0) THEN
!
! Assign the directional weightings
!
    center_2d = 4.0D0/9.0D0
    card_2d = 1.0D0/9.0D0
    diag_2d = 1.0D0/36.0D0
    dir = 9
    nghosts = 1
!    cs = SQRT(1.0D0/3.0D0)
!
! Allocate and assign the values for the weightings
!
!    allocate(wi(dir))
    allocate(cx(dir))
    allocate(cy(dir))
    allocate(cz(dir))
    allocate(rcx(dir))
    allocate(rcy(dir))
    allocate(rcz(dir))
    allocate(cmag(dir))
    allocate(abs_latt_leng(dir))
    !allocate(link_dist(dir,0:amr_max_lvl))
    allocate(opp(dir))
!          wi = (/center_2d, card_2d, card_2d, card_2d, card_2d,&
!            diag_2d, diag_2d, diag_2d,diag_2d/)
!
! Determine the direction of the 2D plane
! -1 - is the x-y plane
! -2 - is the x-z plane
! -3 - is the y-z plane
!
! Values are 1-9
! 1 = (0,0)
! 2 = (+1,0)
! 3 = (0,+1)
! 4 = (-1,0)
! 5 = (0,-1)
! 6 = (+1,+1)
! 7 = (-1,+1)
! 8 = (-1,-1)
! 9 = (+1,-1)
!
!    IF (bounding_method == -1) THEN
! cz will be null
      cx = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
      cy = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)
      cz = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
      opp = (/0, 4, 5, 2, 3, 8, 9, 6, 7/)
!    ELSE IF (bounding_method == -2) THEN
!! cy will be null
!      cx = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
!      cy = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
!      cz = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)
!      opp = (/0, 4, 5, 2, 3, 8, 9, 6, 7/)
!    ELSE IF (bounding_method == -3) THEN
!! cx will be null
!      cx = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
!      cy = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
!      cz = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)
!      opp = (/0, 4, 5, 2, 3, 8, 9, 6, 7/)
!    END IF
    rcx = REAL(cx,kind=dp)
    rcy = REAL(cy,kind=dp)
    rcz = REAL(cz,kind=dp)
    DO i = 1,dir
      cmag(i) = rcx(i)**2+rcy(i)**2+rcz(i)**2
      abs_latt_leng(i) = SQRT(cmag(i))
!      do lvl = 0,amr_max_lvl
!        link_dist(i,lvl) = sqrt(cmag(i))/2.0D0**lvl
!      end do
    END DO

!          WRITE(11,*) ''
!          WRITE(11,201) dimensions
!          WRITE(11,202) dir
!          WRITE(11,103) cs
!          WRITE(11,*) ''
! 201  FORMAT ('Dimensions in this run = ',I8)
! 202  FORMAT ('Number of directions for this run = ',I8)
! 203  FORMAT ('Lattice Speed of Sound = ',F9.7)

  ELSE
    center_fast_2d = 91.0D0/324.0D0
    diag_close_fast_2d = 1.0D0/27.0D0
    card_close_fast_2d = 1.0D0/12.0D0
    card_medium_fast_2d = 7.0D0/360.0D0
    diag_far_fast_2d = 1.0D0/432.0D0
    card_far_fast_2d = 1.0D0/1620.0D0
!    cs = SQRT(2.8D0/3.0D0)
    dir = 21
    nghosts = 3

!    allocate(wi(dir))
    allocate(cx(dir))
    allocate(cy(dir))
    allocate(cz(dir))
    allocate(rcx(dir))
    allocate(rcy(dir))
    allocate(rcz(dir))
!    allocate(link_dist(dir,0:amr_max_lvl))
    allocate(cmag(dir))
    allocate(abs_latt_leng(dir))
    allocate(opp(dir))
!          wi = (/center_fast_2d,card_close_fast_2d,card_close_fast_2d,&
!            card_close_fast_2d,card_close_fast_2d,diag_close_fast_2d,&
!            diag_close_fast_2d,diag_close_fast_2d,diag_close_fast_2d,&
!            card_medium_fast_2d,card_medium_fast_2d,card_medium_fast_2d,&
!            card_medium_fast_2d,diag_far_fast_2d,diag_far_fast_2d,&
!            diag_far_fast_2d,diag_far_fast_2d,card_far_fast_2d,&
!            card_far_fast_2d,card_far_fast_2d,card_far_fast_2d/)
!
! 1 - (0,0)
! 2 - (+1,0)
! 3 - (0,+1)
! 4 - (-1,0)
! 5 - (0,-1)
! 6 - (+1,+1)
! 7 - (-1,+1)
! 8 - (-1,-1)
! 9 - (+1,-1)
! 10 - (+2,0)
! 11 - (0,+2)
! 12 - (-2,0)
! 13 - (0,-2)
! 14 - (+2,+2)
! 15 - (-2,+2)
! 16 - (-2,-2)
! 17 - (+2,-2)
! 18 - (+3,0)
! 19 - (0,+3)
! 20 - (-3,0)
! 21 - (0,-3)
!

!    IF (bounding_method == -1) THEN
! cz will be null
      cx = (/0, 1, 0, -1, 0, 1, -1, -1, 1, 2, 0, -2, 0, 2, -2,&
        -2, 2, 3, 0, -3, 0/)
      cy = (/0, 0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 2, 2, &
        -2, -2, 0, 3, 0, -3/)
      cz = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
        0, 0, 0, 0, 0/)
      opp = (/0, 4, 5, 2, 3, 8, 9, 6, 7, 12, 13, 10, 11, 16, 17,&
        14, 15, 20, 21, 18, 19/)
!    ELSE IF (bounding_method == -2) THEN
!! cy will be null
!      cx = (/0, 1, 0, -1, 0, 1, -1, -1, 1, 2, 0, -2, 0, 2, -2,&
!        -2, 2, 3, 0, -3, 0/)
!      cy = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!        0, 0, 0, 0, 0/)
!      cz = (/0, 0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 2, 2, &
!        -2, -2, 0, 3, 0, -3/)
!      opp = (/0, 4, 5, 2, 3, 8, 9, 6, 7, 12, 13, 10, 11, 16, 17,&
!        14, 15, 20, 21, 18, 19/)
!    ELSE IF (bounding_method == -3) THEN
!! cx will be null
!      cx = (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!        0, 0, 0, 0, 0/)
!      cy = (/0, 1, 0, -1, 0, 1, -1, -1, 1, 2, 0, -2, 0, 2, -2,&
!        -2, 2, 3, 0, -3, 0/)
!      cz = (/0, 0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 2, 2, &
!        -2, -2, 0, 3, 0, -3/)
!      opp = (/0, 4, 5, 2, 3, 8, 9, 6, 7, 12, 13, 10, 11, 16, 17,&
!        14, 15, 20, 21, 18, 19/)
!    ELSE
!      WRITE(*,*) 'How did this error get this far? Bye Bye.'
!      CALL error_out
!    END IF

    rcx = REAL(cx,kind=dp)
    rcy = REAL(cy,kind=dp)
    rcz = REAL(cz,kind=dp)

    DO i = 1,dir
      cmag(i) = rcx(i)**2+rcy(i)**2+rcz(i)**2
      abs_latt_leng(i) = SQRT(cmag(i))
!      do lvl = 0,amr_max_lvl
!        link_dist(i,lvl) = sqrt(cmag(i))/2.0D0**lvl
!      end do
    END DO

!          WRITE(11,*) ''
!          WRITE(11,200) 'Dimensions in this run = ',dimensions
!          WRITE(11,200) 'Number of directions for this run = ',dir
!          WRITE(11,100) 'Lattice Speed of Sound = ',cs
!          WRITE(11,*) ''

  END IF
ELSE IF (dimensions == 3) THEN
!
! Allocate and assign the values for the weightings
!
  IF (mach <= 0.0001D0) THEN !XXXX
!
! Assign the directional weightings
!
!    center_3d = 1.0D0/3.0D0
!    card_3d = 1.0D0/18.0D0
!    diag_3d = 1.0D0/36.0D0
!    cs = SQRT(1.0D0/3.0D0)
!    dir = 19 !for D3Q19
    dir = 27 !for D3Q27
    nghosts = 2
    nghosts_state = 6

!    allocate(wi(dir))
    allocate(cx(dir))
    allocate(cy(dir))
    allocate(cz(dir))
    allocate(rcx(dir))
    allocate(rcy(dir))
    allocate(rcz(dir))
!    allocate(link_dist(dir,0:amr_max_lvl))
    allocate(cmag(dir))
    allocate(abs_latt_leng(dir))
    allocate(opp(dir))
!          wi = (/center_3d,card_3d,card_3d,card_3d,card_3d,card_3d,&
!            card_3d,diag_3d,diag_3d,diag_3d,diag_3d,diag_3d,diag_3d,&
!            diag_3d,diag_3d,diag_3d,diag_3d,diag_3d,diag_3d/)
!
! 1 -  (0,0,0)
! 2 -  (+1,0,0)
! 3 -  (0,+1,0)
! 4 -  (0,0,+1)
! 5 -  (-1,0,0)
! 6 -  (0,-1,0)
! 7 -  (0,0,-1)
! 8 -  (+1,+1,0)
! 9 -  (0,+1,-1)
! 10 - (-1,+1,0)
! 11 - (0,+1,+1)
! 12 - (+1,0,+1)
! 13 - (+1,0,-1)
! 14 - (-1,0,-1)
! 15 - (-1,0,+1)
! 16 - (+1,-1,0)
! 17 - (0,-1,-1)
! 18 - (-1,-1,0)
! 19 - (0,-1,+1)
!                1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19
!    cx = (/0, 1, 0, 0,-1, 0, 0, 1, 0,-1, 0, 1, 1,-1,-1, 1, 0,-1, 0/)
!    cy = (/0, 0, 1, 0, 0,-1, 0, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1/)
!    cz = (/0, 0, 0, 1, 0, 0,-1, 0,-1, 0, 1, 1,-1,-1, 1, 0,-1, 0, 1/)
!    opp = (/0, 5, 6, 7, 2, 3, 4, 18, 19, 16, 17, 14, 15, 12, 13,&
!      10, 11, 8, 9/)
!    rcx = REAL(cx,kind=dp)
!    rcy = REAL(cy,kind=dp)
!    rcz = REAL(cz,kind=dp)
!
!    DO i = 1,dir
!      cmag(i) = rcx(i)**2+rcy(i)**2+rcz(i)**2
!      abs_latt_leng(i) = SQRT(cmag(i))
!!      do lvl = 0,amr_max_lvl
!!        link_dist(i,lvl) = sqrt(cmag(i))/2.0D0**lvl
!!      end do
!    END DO
! 1 -  (0,0,0)
! 2 -  (+1,0,0)
! 3 -  (0,+1,0)
! 4 -  (0,0,+1)
! 5 -  (-1,0,0)
! 6 -  (0,-1,0)
! 7 -  (0,0,-1)
! 8 -  (+1,+1,0)
! 9 -  (0,+1,-1)
! 10 - (-1,+1,0)
! 11 - (0,+1,+1)
! 12 - (+1,0,+1)
! 13 - (+1,0,-1)
! 14 - (-1,0,-1)
! 15 - (-1,0,+1)
! 16 - (+1,-1,0)
! 17 - (0,-1,-1)
! 18 - (-1,-1,0)
! 19 - (0,-1,+1)
! 20 - (+1,+1,+1)
! 21 - (+1,+1,-1)
! 22 - (+1,-1,+1)
! 23 - (+1,-1,-1)
! 24 - (-1,+1,+1)
! 25 - (-1,+1,-1)
! 26 - (-1,-1,+1)
! 27 - (-1,-1,-1)
!          1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
!         20,21,22,23,24,25,26,27
    cx = (/0, 1, 0, 0,-1, 0, 0, 1, 0,-1, 0, 1, 1,-1,-1, 1, 0,-1, 0,&
           1, 1, 1, 1,-1,-1,-1,-1/)
    cy = (/0, 0, 1, 0, 0,-1, 0, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,&
           1, 1,-1,-1, 1, 1,-1,-1/)
    cz = (/0, 0, 0, 1, 0, 0,-1, 0,-1, 0, 1, 1,-1,-1, 1, 0,-1, 0, 1,&
           1,-1, 1,-1, 1,-1, 1,-1/)
    opp = (/1, 5, 6, 7, 2, 3, 4, 18, 19, 16, 17, 14, 15, 12, 13,&
      10, 11, 8, 9,27,26,25,24,23,22,21,20/)
    rcx = REAL(cx,kind=dp)
    rcy = REAL(cy,kind=dp)
    rcz = REAL(cz,kind=dp)

    DO i = 1,dir
      cmag(i) = rcx(i)**2+rcy(i)**2+rcz(i)**2
      abs_latt_leng(i) = SQRT(cmag(i))
!      do lvl = 0,amr_max_lvl
!        link_dist(i,lvl) = sqrt(cmag(i))/2.0D0**lvl
!      end do
    END DO

!          WRITE(11,*) ''
!          WRITE(11,200) dimensions
!          WRITE(11,200) dir
!          WRITE(11,100) cs
!          WRITE(11,*) ''

  ELSE
    center_fast_3d = 1.0D0/12.0D0
    card_close_fast_3d = 1.0D0/12.0D0
    diag_close_fast_3d = 1.0D0/27.0D0
    card_medium_fast_3d = 2.0D0/135.0D0
    diag_far_fast_3d = 1.0D0/432.0D0
    card_far_fast_3d = 1.0D0/1620.0D0
!    cs = SQRT(2.8D0/3.0D0)
    dir = 39

    if (shifted) then
      nghosts = 8
      nghosts_mid = 10
      nghosts_state = 10
    else
      nghosts = 6
      nghosts_mid = 8
      nghosts_state = 8

    end if

!    allocate(wi(dir))
    allocate(cx(dir))
    allocate(cy(dir))
    allocate(cz(dir))
    allocate(rcx(dir))
    allocate(rcy(dir))
    allocate(rcz(dir))
!    allocate(link_dist(dir,0:amr_max_lvl))
    allocate(cmag(dir))
    allocate(abs_latt_leng(dir))
    allocate(opp(dir))

!          wi = (/center_fast_3d,card_close_fast_3d,card_close_fast_3d,&
!            card_close_fast_3d,card_close_fast_3d,card_close_fast_3d,&
!            card_close_fast_3d,diag_close_fast_3d,diag_close_fast_3d,&
!            diag_close_fast_3d,diag_close_fast_3d,diag_close_fast_3d,&
!            diag_close_fast_3d,diag_close_fast_3d,diag_close_fast_3d,&
!            card_medium_fast_3d,card_medium_fast_3d,card_medium_fast_3d,&
!            card_medium_fast_3d,card_medium_fast_3d,card_medium_fast_3d,&
!            diag_far_fast_3d,diag_far_fast_3d,diag_far_fast_3d,&
!            diag_far_fast_3d,diag_far_fast_3d,diag_far_fast_3d,&
!            diag_far_fast_3d,diag_far_fast_3d,diag_far_fast_3d,&
!            diag_far_fast_3d,diag_far_fast_3d,diag_far_fast_3d,&
!            card_far_fast_3d,card_far_fast_3d,card_far_fast_3d,&
!            card_far_fast_3d,card_far_fast_3d,card_far_fast_3d/)
!
!
! 1 - (0,0,0)
! 2 - (+1,0,0)
! 3 - (0,+1,0)
! 4 - (0,0,+1)
! 5 - (-1,0,0)
! 6 - (0,-1,0)
! 7 - (0,0,-1)
! 8 - (+1,+1,+1)
! 9 - (-1,+1,+1)
! 10 - (-1,-1,+1)
! 11 - (+1,-1,+1)
! 12 - (+1,+1,-1)
! 13 - (-1,+1,-1)
! 14 - (-1,-1,-1)
! 15 - (+1,-1,-1)
! 16 - (+2,0,0)
! 17 - (0,+2,0)
! 18 - (0,0,+2)
! 19 - (-2,0,0)
! 20 - (0,-2,0)
! 21 - (0,0,-2)
! 22 - (+2,+2,0)
! 23 - (+2,0,+2)
! 24 - (0,+2,+2)
! 25 - (-2,+2,0)
! 26 - (-2,-2,0)
! 27 - (+2,-2,0)
! 28 - (+2,0,-2)
! 29 - (-2,0,-2)
! 30 - (-2,0,+2)
! 31 - (0,+2,-2)
! 32 - (0,-2,-2)
! 33 - (0,-2,+2)
! 34 - (+3,0,0)
! 35 - (0,+3,0)
! 36 - (0,0,+3)
! 37 - (-3,0,0)
! 38 - (0,-3,0)
! 39 - (0,0,-3)
!
!                1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,
!                18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,
!                 35,36,37,38,39
    cx = (/0, 1, 0, 0,-1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 2, 0, 0,&
      -2, 0, 0, 2, 2, 0,-2,-2, 2, 2,-2,-2, 0, 0, 0, 3, 0, 0,-3,&
      0, 0/)
    cy = (/0, 0, 1, 0, 0,-1, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 2, 0,&
      0,-2, 0, 2, 0, 2, 2,-2,-2, 0, 0, 0, 2,-2,-2, 0, 3, 0, 0,&
      -3, 0/)
    cz = (/0, 0, 0, 1, 0, 0,-1, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 2,&
      0, 0,-2, 0, 2, 2, 0, 0, 0,-2,-2, 2,-2,-2, 2, 0, 0, 3, 0,&
      0,-3/)
!
!
!
      allocate(cxc2(dir))
      allocate(cyc2(dir))
      allocate(czc2(dir))
      allocate(cx2(dir))
      allocate(cy2(dir))
      allocate(cz2(dir))
      allocate(cxcy(dir))
      allocate(cxcz(dir))
      allocate(cycz(dir))

      allocate(cxxx(dir))
      allocate(cyyy(dir))
      allocate(czzz(dir))
      allocate(cxyz(dir))
      allocate(cxxy(dir))
      allocate(cxyy(dir))
      allocate(cxxz(dir))
      allocate(cxzz(dir))
      allocate(cyyz(dir))
      allocate(cyzz(dir))

      allocate(yi(27,dir))
      allocate(kroen_val(27,dir))
      allocate(p_coeff(27,dir))
      allocate(p_kroen(27))
      allocate(p_kr_val(27))
!
!
!
    if (mach >= 1.2D0 .or. use_shifted) then

      shifted = .true.

      if (primary_flow_dir == 1) then
        cx = cx + 1
      else if (primary_flow_dir == 2) then
        cy = cy + 1
      else if (primary_flow_dir == 3) then
        cz = cz + 1
      else if (primary_flow_dir == 4) then
        cx = cx - 1
      else if (primary_flow_dir == 5) then
        cy = cy - 1
      else if (primary_flow_dir == 6) then
        cz = cz - 1
      end if



      opp = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,&
        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,&
        32, 33, 34, 35, 36, 37, 38, 39/)



    else
      shifted = .false.
      opp = (/1,5, 6, 7, 2, 3, 4,14,15,12,13,10,11, 8, 9,19, 20,&
        21,16,17,18,26,29,32,27,22,25,30,23,28,33,24,31,37,38,39,&
        34,35,36/)
    end if
!
!
!

    rcx = REAL(cx,kind=dp)
    rcy = REAL(cy,kind=dp)
    rcz = REAL(cz,kind=dp)

    DO i = 1,dir
      cmag(i) = rcx(i)**2+rcy(i)**2+rcz(i)**2
      abs_latt_leng(i) = SQRT(cmag(i))

!      if (mach >= 1.2D0) then
        cx2(i) = cx(i)**2
        cy2(i) = cy(i)**2
        cz2(i) = cz(i)**2
        cxcy(i) = cx(i)*cy(i)
        cxcz(i) = cx(i)*cz(i)
        cycz(i) = cy(i)*cz(i)
        cxc2(i) = cx(i)*int(cmag(i))
        cyc2(i) = cy(i)*int(cmag(i))
        czc2(i) = cz(i)*int(cmag(i))

        cxxx(i) = cx(i)**3
        cyyy(i) = cy(i)**3
        czzz(i) = cz(i)**3
        cxyz(i) = cx(i) * cy(i) * cz(i)
        cxxy(i) = cx(i)**2 * cy(i)
        cxyy(i) = cx(i) * cy(i)**2
        cxxz(i) = cx(i)**2 * cz(i)
        cxzz(i) = cx(i) * cz(i)**2
        cyyz(i) = cy(i)**2 * cz(i)
        cyzz(i) = cz(i)**2 * cy(i)

!      end if
!      do lvl = 0,amr_max_lvl
!        link_dist(i,lvl) = sqrt(cmag(i))/2.0D0**lvl
!      end do
    END DO
!write(*,*) 'cx2',cx2
!write(*,*) 'cy2',cy2
!write(*,*) 'cz2',cz2
!write(*,*) 'cxcy',cxcy
!write(*,*) 'cxcz',cxcz
!write(*,*) 'cycz',cycz
!write(*,*) 'cxc2',cxc2
!write(*,*) 'cyc2',cyc2
!write(*,*) 'czc2',czc2
!
!write(*,*) 'cxxx',cxxx
!write(*,*) 'cyyy',cyyy
!write(*,*) 'czzz',czzz
!write(*,*) 'cxxy',cxxy
!write(*,*) 'cxxz',cxxz
!write(*,*) 'cxyy',cxyy
!write(*,*) 'cxzz',cxzz
!write(*,*) 'cxyz',cxyz
!write(*,*) 'cyyz',cyyz
!write(*,*) 'cyzz',cyzz


!
!1 = cxxx, delta = x
!2 = cxxy, delta = y
!3 = cxxz, delta = z
!4 = cxxy, delta = 0
!5 = cxyy, delta = 0
!6 = cxyz, delta = 0
!7 = cxxz, delta = 0
!8 = cxyz, delta = 0
!9 = cxzz, delta = 0
!10 = cxxy, delta = 0
!11 = cxyy, delta = 0
!12 = cxyz, delta = 0
!13 = cxyy, delta = x
!14 = cyyy, delta = y
!15 = cyyz, delta = z
!16 = cxyz, delta = 0
!17 = cyyz, delta = 0
!18 = cyzz, delta = 0
!19 = cxxz, delta = 0
!20 = cxyz, delta = 0
!21 = cxzz, delta = 0
!22 = cxyz, delta = 0
!23 = cyyz, delta = 0
!24 = cyzz, delta = 0
!25 = cxzz, delta = x
!26 = cyzz, delta = y
!27 = czzz, delta = z
!
!1 = cxx, pk = 1, pkv = cx
!2 = cxx, pk = 1, pkv = cy
!3 = cxx, pk = 1, pkv = cz
!4 = cxy, pk = 0, pkv = 0
!5 = cxy, pk = 0, pkv = 0
!6 = cxy, pk = 0, pkv = 0
!7 = cxz, pk = 0, pkv = 0
!8 = cxz, pk = 0, pkv = 0
!9 = cxz, pk = 0, pkv = 0
!10 =cxy, pk = 0, pkv = 0
!11 =cxy, pk = 0, pkv = 0
!12 =cxy, pk = 0, pkv = 0
!13 =cyy, pk = 1, pkv = cx
!14 =cyy, pk = 1, pkv = cy
!15 =cyy, pk = 1, pkv = cz
!16 =cyz, pk = 0, pkv = 0
!17 =cyz, pk = 0, pkv = 0
!18 =cyz, pk = 0, pkv = 0
!19 =cxz, pk = 0, pkv = 0
!20 =cxz, pk = 0, pkv = 0
!21 =cxz, pk = 0, pkv = 0
!22 =cyz, pk = 0, pkv = 0
!23 =cyz, pk = 0, pkv = 0
!24 =cyz, pk = 0, pkv = 0
!25 =czz, pk = 1, pkv = cx
!26 =czz, pk = 1, pkv = cy
!27 =czz, pk = 1, pkv = cz
!
!  write(*,*) 'boost'
!    if (mach >= 1.2D0) then
    do j = 1,dir
      do i = 1,27
        select case(i)
        case(1)
          yi(i,j) = cxxx(j)
          kroen_val(i,j) = cx(j)
          p_coeff(i,j) = cx(j)**2
          p_kroen(i) = 1
          p_kr_val = cx(j)
        case(2)
          yi(i,j) = cxxy(j)
          kroen_val(i,j) = cy(j)
          p_coeff(i,j) = cx(j)**2
          p_kroen(i) = 1
          p_kr_val = cy(j)
        case(3)
          yi(i,j) = cxxz(j)
          kroen_val(i,j) = cz(j)
          p_coeff(i,j) = cx(j)**2
          p_kroen(i) = 1
          p_kr_val = cz(j)
        case(4)
          yi(i,j) = cxxy(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cx(j)*cy(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(5)
          yi(i,j) = cxyy(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cx(j)*cy(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(6)
          yi(i,j) = cxyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(7)
          yi(i,j) = cxxz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cx(j)*cz(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(8)
          yi(i,j) = cxyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cx(j)*cz(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(9)
          yi(i,j) = cxzz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cx(j)*cz(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(10)
          yi(i,j) = cxxy(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(11)
          yi(i,j) = cxyy(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(12)
          yi(i,j) = cxyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(13)
          yi(i,j) = cxyy(j)
          kroen_val(i,j) = cx(j)
          p_coeff(i,j) = cy(j)**2
          p_kroen(i) = 1
          p_kr_val = cx(j)
        case(14)
          yi(i,j) = cyyy(j)
          kroen_val(i,j) = cy(j)
          p_coeff(i,j) = cy(j)**2
          p_kroen(i) = 1
          p_kr_val = cy(j)
        case(15)
          yi(i,j) = cyyz(j)
          kroen_val(i,j) = cz(j)
          p_coeff(i,j) = cy(j)**2
          p_kroen(i) = 1
          p_kr_val = cz(j)
        case(16)
          yi(i,j) = cxyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cz(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(17)
          yi(i,j) = cyyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cz(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(18)
          yi(i,j) = cyzz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cy(j)*cz(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(19)
          yi(i,j) = cxxz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cz(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(20)
          yi(i,j) = cxyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cz(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(21)
          yi(i,j) = cxzz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cz(j)*cx(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(22)
          yi(i,j) = cxyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cz(j)*cy(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(23)
          yi(i,j) = cyyz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cz(j)*cy(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(24)
          yi(i,j) = cyzz(j)
          kroen_val(i,j) = 0
          p_coeff(i,j) = cz(j)*cy(j)
          p_kroen(i) = 0
          p_kr_val = 0
        case(25)
          yi(i,j) = cxzz(j)
          kroen_val(i,j) = cx(j)
          p_coeff(i,j) = cz(j)**2
          p_kroen(i) = 1
          p_kr_val = cx(j)
        case(26)
          yi(i,j) = cyzz(j)
          kroen_val(i,j) = cy(j)
          p_coeff(i,j) = cz(j)**2
          p_kroen(i) = 1
          p_kr_val = cy(j)
        case(27)
          yi(i,j) = czzz(j)
          kroen_val(i,j) = cz(j)
          p_coeff(i,j) = cz(j)**2
          p_kroen(i) = 1
          p_kr_val = cz(j)
        end select
      end do
    end do
!    end if
!    write(*,*) 'paste'

  END IF
ELSE
  WRITE(*,*) 'How did this error get this far?'
  WRITE(11,404) dimensions,mach
  CALL error_out
END IF

write(11,*) ''
write(11,201) dimensions
write(11,202) dir
!write(11,203) cs
write(11,*) ''
write(11,501) cx
write(11,502) cy
write(11,503) cz
write(11,505) opp
!write(11,504) wi

 201  FORMAT ('Dimensions in this run = ',I8)
 202  FORMAT ('Number of directions for this run = ',I8)
 203  FORMAT ('Lattice Speed of Sound = ',F9.7)
 404  FORMAT ('Issues in assigning directions and weight, dimensions = ',&
        I4,' and Mach number = ',F5.3)
 501  FORMAT ('X-direction array cx,     ',39I3,' ')
 502  FORMAT ('Y-direction array cy,     ',39I3,' ')
 503  FORMAT ('Z-direction array cz,     ',39I3,' ')
 504  FORMAT ('Weighting array wi, ',39F6.4,' ')
 505  FORMAT ('Opposite direction array, ',39I3,' ')
end subroutine
