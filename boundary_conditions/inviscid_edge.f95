SUBROUTINE inviscid_edge(fout_bc,fout_opp_bc,gout_bc,gout_opp_bc,&
      i,bc_num,loc_state,inviscid_dir,loc_cx,loc_cy,loc_cz)
!
! This subroutine returns the dir required for inviscid and symmetry
! boundary conditions.
!
! Calls:
! Called by: loop
!
USE precise
USE constants
USE grid_data
USE linkwise
IMPLICIT NONE
INTEGER :: inviscid_dir
INTEGER,INTENT(IN) :: i,loc_state(dir),bc_num
integer :: loc_cx,loc_cy,loc_cz
real(kind=dp) :: fout_bc,fout_opp_bc,gout_bc,gout_opp_bc

!
! Need to determine which lattice we are using first
!
! Inviscid directions for the D2Q9 lattice
!
  loc_cx = 0
  loc_cy = 0
  loc_cz = 0
  SELECT CASE (i)
!
! Basic cardinal directions
!
  CASE(2)
    inviscid_dir = 5
    loc_cx = 0
    loc_cy = 0
    loc_cz = 0
  CASE(3)
    inviscid_dir = 6
    loc_cx = 0
    loc_cy = 0
    loc_cz = 0
  CASE(4)
    inviscid_dir = 7
    loc_cx = 0
    loc_cy = 0
    loc_cz = 0
  CASE(5)
    inviscid_dir = 2
    loc_cx = 0
    loc_cy = 0
    loc_cz = 0
  CASE(6)
    inviscid_dir = 3
    loc_cx = 0
    loc_cy = 0
    loc_cz = 0
  CASE(7)
    inviscid_dir = 4
    loc_cx = 0
    loc_cy = 0
    loc_cz = 0
!
! Close diagonal directions, requires three part if statement because it's in
! three directions
!
! +1,+1,+1
  CASE(8)
    IF (loc_state(2) == -11 .and. loc_state(3) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 14
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(2) == -11 .and. loc_state(3) == -11) then
      inviscid_dir = 10
      loc_cx = 0
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(2) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 13
      loc_cx = 0
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(3) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 15
      loc_cx = 1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(2) == -11) THEN
      inviscid_dir = 9
      loc_cx = 0
      loc_cy = 1
      loc_cz = 1
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 11
      loc_cx = 1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(4) == -11) then
      inviscid_dir = 12
      loc_cx = 1
      loc_cy = 1
      loc_cz = 0
    END IF
! -1,+1,+1
  CASE(9)
    IF (loc_state(5) == -11 .and. loc_state(3) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 15
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(5) == -11 .and. loc_state(3) == -11) then
      inviscid_dir = 12
      loc_cx = 0
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(5) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 11
      loc_cx = 0
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(3) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 14
      loc_cx = -1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(5) == -11) THEN
      inviscid_dir = 8
      loc_cx = 0
      loc_cy = 1
      loc_cz = 1
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 10
      loc_cx = -1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(4) == -11) then
      inviscid_dir = 13
      loc_cx = -1
      loc_cy = 1
      loc_cz = 0
    END IF
! -1,-1,+1
  CASE(10)
    IF (loc_state(5) == -11 .and. loc_state(6) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 12
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(5) == -11 .and. loc_state(6) == -11) then
      inviscid_dir = 8
      loc_cx = 0
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(5) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 15
      loc_cx = 0
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(6) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 13
      loc_cx = -1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(5) == -11) THEN
      inviscid_dir = 11
      loc_cx = 0
      loc_cy = -1
      loc_cz = 1
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 9
      loc_cx = -1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(4) == -11) then
      inviscid_dir = 14
      loc_cx = -1
      loc_cy = -1
      loc_cz = 0
    END IF
! +1,-1,+1
  CASE(11)
    IF (loc_state(2) == -11 .and. loc_state(6) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 13
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(2) == -11 .and. loc_state(6) == -11) then
      inviscid_dir = 9
      loc_cx = 0
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(2) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 14
      loc_cx = 0
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(6) == -11 .and. loc_state(4) == -11) then
      inviscid_dir = 12
      loc_cx = 1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(2) == -11) THEN
      inviscid_dir = 10
      loc_cx = 0
      loc_cy = -1
      loc_cz = 1
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 8
      loc_cx = 1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(4) == -11) then
      inviscid_dir = 15
      loc_cx = 1
      loc_cy = -1
      loc_cz = 0
    END IF
! +1,+1,-1
  CASE(12)
    IF (loc_state(2) == -11 .and. loc_state(3) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 10
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(2) == -11 .and. loc_state(3) == -11) then
      inviscid_dir = 14
      loc_cx = 0
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(2) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 9
      loc_cx = 0
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(3) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 11
      loc_cx = 1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(2) == -11) THEN
      inviscid_dir = 13
      loc_cx = 0
      loc_cy = 1
      loc_cz = -1
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 15
      loc_cx = 1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(7) == -11) then
      inviscid_dir = 8
      loc_cx = 1
      loc_cy = 1
      loc_cz = 0
    END IF
! -1,+1,-1
  CASE(13)
    IF (loc_state(5) == -11 .and. loc_state(3) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 11
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(5) == -11 .and. loc_state(3) == -11) then
      inviscid_dir = 15
      loc_cx = 0
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(5) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 8
      loc_cx = 0
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(3) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 10
      loc_cx = -1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(5) == -11) THEN
      inviscid_dir = 12
      loc_cx = 0
      loc_cy = 1
      loc_cz = -1
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 14
      loc_cx = -1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(7) == -11) then
      inviscid_dir = 9
      loc_cx = -1
      loc_cy = 1
      loc_cz = 0
    END IF
! -1,-1,-1
  CASE(14)
    IF (loc_state(5) == -11 .and. loc_state(6) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 8
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(5) == -11 .and. loc_state(6) == -11) then
      inviscid_dir = 12
      loc_cx = 0
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(5) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 11
      loc_cx = 0
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(6) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 9
      loc_cx = -1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(5) == -11) THEN
      inviscid_dir = 15
      loc_cx = 0
      loc_cy = -1
      loc_cz = -1
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 13
      loc_cx = -1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(7) == -11) then
      inviscid_dir = 10
      loc_cx = -1
      loc_cy = -1
      loc_cz = 0
    END IF
! +1,-1,-1
  CASE(15)
    IF (loc_state(2) == -11 .and. loc_state(6) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 9
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(2) == -11 .and. loc_state(6) == -11) then
      inviscid_dir = 13
      loc_cx = 0
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(2) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 10
      loc_cx = 0
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(6) == -11 .and. loc_state(7) == -11) then
      inviscid_dir = 8
      loc_cx = 1
      loc_cy = 0
      loc_cz = 0
    ELSE IF (loc_state(2) == -11) THEN
      inviscid_dir = 14
      loc_cx = 0
      loc_cy = -1
      loc_cz = -1
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 12
      loc_cx = 1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(7) == -11) then
      inviscid_dir = 11
      loc_cx = 1
      loc_cy = -1
      loc_cz = 0
    END IF
!
! Mid-range cardinal directions
!
  CASE(16)
    inviscid_dir = 19
    if (loc_state(2) == -11) then
      loc_cx = -1
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = 1
      loc_cy = 0
      loc_cz = 0
    end if

  CASE(17)
    inviscid_dir = 20
    if (loc_state(3) == -11) then
      loc_cx = 0
      loc_cy = -1
      loc_cz = 0
    else
      loc_cx = 0
      loc_cy = 1
      loc_cz = 0
    end if

  CASE(18)
    inviscid_dir = 21
    if (loc_state(4) == -11) then
    loc_cx = 0
    loc_cy = 0
    loc_cz = -1
    else
    loc_cx = 0
    loc_cy = 0
    loc_cz = 1
    end if

  CASE(19)
    inviscid_dir = 16
    if (loc_state(5) == -11) then
      loc_cx = 1
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = -1
      loc_cy = 0
      loc_cz = 0
    end if

  CASE(20)
    inviscid_dir = 17
    if (loc_state(6) == -11) then
      loc_cx = 0
      loc_cy = 1
      loc_cz = 0
    else
      loc_cx = 0
      loc_cy = -1
      loc_cz = 0
    end if

  CASE(21)
    inviscid_dir = 18
    if (loc_state(7) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 1
    else
      loc_cx = 0
      loc_cy = 0
      loc_cz = -1
    end if

!
! Mid-range diagonal directions
!
  CASE(22)
    IF (loc_state(2) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 26
      loc_cx = -1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(17) == -11 .and. loc_state(2) == -11) THEN
      inviscid_dir = 26
      loc_cx = -1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(16) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 26
      loc_cx = 1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(16) == -11 .and. loc_state(17) == -11) THEN
      inviscid_dir = 26
      loc_cx = 1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(2) == -11) then
      inviscid_dir = 25
      loc_cx = -1
      loc_cy = 2
      loc_cz = 0
    else if (loc_state(16) == -11) then
      inviscid_dir = 25
      loc_cx = 1
      loc_cy = 2
      loc_cz = 0
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 27
      loc_cx = 2
      loc_cy = -1
      loc_cz = 0
    ELSE IF (loc_state(17) == -11) THEN
      inviscid_dir = 27
      loc_cx = 2
      loc_cy = 1
      loc_cz = 0
    END IF
  CASE(23)
    IF (loc_state(2) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 29
      loc_cx = -1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(18) == -11 .and. loc_state(2) == -11) THEN
      inviscid_dir = 29
      loc_cx = -1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(16) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 29
      loc_cx = 1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(16) == -11 .and. loc_state(18) == -11) THEN
      inviscid_dir = 29
      loc_cx = 1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(2) == -11) then
      inviscid_dir = 30
      loc_cx = -1
      loc_cy = 0
      loc_cz = 2
    else if (loc_state(16) == -11) then
      inviscid_dir = 30
      loc_cx = 1
      loc_cy = 0
      loc_cz = 2
    ELSE IF (loc_state(4) == -11) THEN
      inviscid_dir = 28
      loc_cx = 2
      loc_cy = 0
      loc_cz = -1
    ELSE IF (loc_state(18) == -11) THEN
      inviscid_dir = 28
      loc_cx = 2
      loc_cy = 0
      loc_cz = 1
    END IF
!    IF (loc_state(16) == -11 .and. loc_state(18) == -11) THEN
!      inviscid_dir = 29
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(16) == -11) then
!      inviscid_dir = 30
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 2
!    ELSE IF (loc_state(18) == -11) THEN
!      inviscid_dir = 28
!      loc_cx = 2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(24)
    IF (loc_state(4) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = -1
      loc_cz = -1
    else if (loc_state(17) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = 1
      loc_cz = -1
    else if (loc_state(18) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = -1
      loc_cz = 1
    else if (loc_state(18) == -11 .and. loc_state(17) == -11) THEN
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = 1
      loc_cz = 1
    else if (loc_state(4) == -11) then
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = 2
      loc_cz = -1
    else if (loc_state(18) == -11) then
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = 2
      loc_cz = 1
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = -1
      loc_cz = 2
    ELSE IF (loc_state(17) == -11) THEN
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = 1
      loc_cz = 2
    END IF
!    IF (loc_state(17) == -11 .and. loc_state(18) == -11) THEN
!      inviscid_dir = 32
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(17) == -11) then
!      inviscid_dir = 33
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 2
!    ELSE IF (loc_state(18) == -11) THEN
!      inviscid_dir = 31
!      loc_cx = 0
!      loc_cy = 2
!      loc_cz = 0
!    END IF
  CASE(25)
    IF (loc_state(5) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 27
      loc_cx = 1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(17) == -11 .and. loc_state(5) == -11) THEN
      inviscid_dir = 27
      loc_cx = 1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(19) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 27
      loc_cx = -1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(19) == -11 .and. loc_state(17) == -11) THEN
      inviscid_dir = 27
      loc_cx = -1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(5) == -11) then
      inviscid_dir = 22
      loc_cx = 1
      loc_cy = 2
      loc_cz = 0
    else if (loc_state(19) == -11) then
      inviscid_dir = 22
      loc_cx = -1
      loc_cy = 2
      loc_cz = 0
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 26
      loc_cx = -2
      loc_cy = -1
      loc_cz = 0
    ELSE IF (loc_state(17) == -11) THEN
      inviscid_dir = 26
      loc_cx = -2
      loc_cy = 1
      loc_cz = 0
    END IF
!    IF (loc_state(19) == -11 .and. loc_state(17) == -11) THEN
!      inviscid_dir = 27
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(19) == -11) then
!      inviscid_dir = 22
!      loc_cx = 0
!      loc_cy = 2
!      loc_cz = 0
!    ELSE IF (loc_state(17) == -11) THEN
!      inviscid_dir = 26
!      loc_cx = -2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(26)
    IF (loc_state(5) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 22
      loc_cx = 1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(20) == -11 .and. loc_state(5) == -11) THEN
      inviscid_dir = 22
      loc_cx = 1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(19) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 22
      loc_cx = -1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(19) == -11 .and. loc_state(20) == -11) THEN
      inviscid_dir = 22
      loc_cx = -1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(5) == -11) then
      inviscid_dir = 27
      loc_cx = 1
      loc_cy = -2
      loc_cz = 0
    else if (loc_state(19) == -11) then
      inviscid_dir = 27
      loc_cx = -1
      loc_cy = -2
      loc_cz = 0
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 25
      loc_cx = -2
      loc_cy = 1
      loc_cz = 0
    ELSE IF (loc_state(20) == -11) THEN
      inviscid_dir = 25
      loc_cx = -2
      loc_cy = -1
      loc_cz = 0
    END IF
!    IF (loc_state(19) == -11 .and. loc_state(20) == -11) THEN
!      inviscid_dir = 22
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(19) == -11) then
!      inviscid_dir = 27
!      loc_cx = 0
!      loc_cy = -2
!      loc_cz = 0
!    ELSE IF (loc_state(20) == -11) THEN
!      inviscid_dir = 25
!      loc_cx = -2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(27)
    IF (loc_state(2) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 25
      loc_cx = -1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(20) == -11 .and. loc_state(2) == -11) THEN
      inviscid_dir = 25
      loc_cx = -1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(16) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 25
      loc_cx = 1
      loc_cy = 1
      loc_cz = 0
    else if (loc_state(16) == -11 .and. loc_state(20) == -11) THEN
      inviscid_dir = 25
      loc_cx = 1
      loc_cy = -1
      loc_cz = 0
    else if (loc_state(2) == -11) then
      inviscid_dir = 26
      loc_cx = -1
      loc_cy = -2
      loc_cz = 0
    else if (loc_state(16) == -11) then
      inviscid_dir = 26
      loc_cx = 1
      loc_cy = -2
      loc_cz = 0
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 22
      loc_cx = 2
      loc_cy = -1
      loc_cz = 0
    ELSE IF (loc_state(20) == -11) THEN
      inviscid_dir = 22
      loc_cx = 2
      loc_cy = 1
      loc_cz = 0
    END IF
!    IF (loc_state(16) == -11 .and. loc_state(20) == -11) THEN
!      inviscid_dir = 25
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(16) == -11) then
!      inviscid_dir = 26
!      loc_cx = 0
!      loc_cy = -2
!      loc_cz = 0
!    ELSE IF (loc_state(20) == -11) THEN
!      inviscid_dir = 22
!      loc_cx = 2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(28)
    IF (loc_state(2) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 30
      loc_cx = -1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(21) == -11 .and. loc_state(2) == -11) THEN
      inviscid_dir = 30
      loc_cx = -1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(16) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 30
      loc_cx = 1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(16) == -11 .and. loc_state(21) == -11) THEN
      inviscid_dir = 30
      loc_cx = 1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(2) == -11) then
      inviscid_dir = 29
      loc_cx = -1
      loc_cy = 0
      loc_cz = -2
    else if (loc_state(16) == -11) then
      inviscid_dir = 29
      loc_cx = 1
      loc_cy = 0
      loc_cz = -2
    ELSE IF (loc_state(7) == -11) THEN
      inviscid_dir = 23
      loc_cx = 2
      loc_cy = 0
      loc_cz = 1
    ELSE IF (loc_state(21) == -11) THEN
      inviscid_dir = 23
      loc_cx = 2
      loc_cy = 0
      loc_cz = -1
    END IF
!    IF (loc_state(16) == -11 .and. loc_state(21) == -11) THEN
!      inviscid_dir = 30
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(16) == -11) then
!      inviscid_dir = 29
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = -2
!    ELSE IF (loc_state(21) == -11) THEN
!      inviscid_dir = 23
!      loc_cx = 2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(29)
    IF (loc_state(5) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 23
      loc_cx = 1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(21) == -11 .and. loc_state(5) == -11) THEN
      inviscid_dir = 23
      loc_cx = 1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(19) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 23
      loc_cx = -1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(19) == -11 .and. loc_state(21) == -11) THEN
      inviscid_dir = 23
      loc_cx = -1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(5) == -11) then
      inviscid_dir = 28
      loc_cx = 1
      loc_cy = 0
      loc_cz = -2
    else if (loc_state(19) == -11) then
      inviscid_dir = 28
      loc_cx = -1
      loc_cy = 0
      loc_cz = -2
    ELSE IF (loc_state(7) == -11) THEN
      inviscid_dir = 30
      loc_cx = -2
      loc_cy = 0
      loc_cz = 1
    ELSE IF (loc_state(21) == -11) THEN
      inviscid_dir = 30
      loc_cx = -2
      loc_cy = 0
      loc_cz = -1
    END IF
!    IF (loc_state(19) == -11 .and. loc_state(21) == -11) THEN
!      inviscid_dir = 23
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(19) == -11) then
!      inviscid_dir = 28
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = -2
!    ELSE IF (loc_state(21) == -11) THEN
!      inviscid_dir = 30
!      loc_cx = -2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(30)
    IF (loc_state(5) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 28
      loc_cx = 1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(18) == -11 .and. loc_state(5) == -11) THEN
      inviscid_dir = 28
      loc_cx = 1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(19) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 28
      loc_cx = -1
      loc_cy = 0
      loc_cz = -1
    else if (loc_state(19) == -11 .and. loc_state(18) == -11) THEN
      inviscid_dir = 28
      loc_cx = -1
      loc_cy = 0
      loc_cz = 1
    else if (loc_state(5) == -11) then
      inviscid_dir = 23
      loc_cx = 1
      loc_cy = 0
      loc_cz = 2
    else if (loc_state(19) == -11) then
      inviscid_dir = 23
      loc_cx = -1
      loc_cy = 0
      loc_cz = 2
    ELSE IF (loc_state(4) == -11) THEN
      inviscid_dir = 29
      loc_cx = -2
      loc_cy = 0
      loc_cz = -1
    ELSE IF (loc_state(18) == -11) THEN
      inviscid_dir = 29
      loc_cx = -2
      loc_cy = 0
      loc_cz = 1
    END IF
!    IF (loc_state(19) == -11 .and. loc_state(18) == -11) THEN
!      inviscid_dir = 28
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(19) == -11) then
!      inviscid_dir = 23
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 2
!    ELSE IF (loc_state(18) == -11) THEN
!      inviscid_dir = 29
!      loc_cx = -2
!      loc_cy = 0
!      loc_cz = 0
!    END IF
  CASE(31)
    IF (loc_state(7) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = -1
      loc_cz = 1
    else if (loc_state(17) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = 1
      loc_cz = 1
    else if (loc_state(21) == -11 .and. loc_state(3) == -11) THEN
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = -1
      loc_cz = -1
    else if (loc_state(21) == -11 .and. loc_state(17) == -11) THEN
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = 1
      loc_cz = -1
    else if (loc_state(7) == -11) then
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = 2
      loc_cz = 1
    else if (loc_state(21) == -11) then
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = 2
      loc_cz = -1
    ELSE IF (loc_state(3) == -11) THEN
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = -1
      loc_cz = -2
    ELSE IF (loc_state(17) == -11) THEN
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = 1
      loc_cz = -2
    END IF
!    IF (loc_state(17) == -11 .and. loc_state(21) == -11) THEN
!      inviscid_dir = 33
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(17) == -11) then
!      inviscid_dir = 32
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = -2
!    ELSE IF (loc_state(21) == -11) THEN
!      inviscid_dir = 24
!      loc_cx = 0
!      loc_cy = 2
!      loc_cz = 0
!    END IF
  CASE(32)
    IF (loc_state(7) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = 1
      loc_cz = 1
    else if (loc_state(20) == -11 .and. loc_state(7) == -11) THEN
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = -1
      loc_cz = 1
    else if (loc_state(21) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = 1
      loc_cz = -1
    else if (loc_state(21) == -11 .and. loc_state(20) == -11) THEN
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = -1
      loc_cz = -1
    else if (loc_state(7) == -11) then
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = -2
      loc_cz = 1
    else if (loc_state(21) == -11) then
      inviscid_dir = 33
      loc_cx = 0
      loc_cy = -2
      loc_cz = -1
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = 1
      loc_cz = -2
    ELSE IF (loc_state(20) == -11) THEN
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = -1
      loc_cz = -2
    END IF
!    IF (loc_state(20) == -11 .and. loc_state(21) == -11) THEN
!      inviscid_dir = 24
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(20) == -11) then
!      inviscid_dir = 31
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = -2
!    ELSE IF (loc_state(21) == -11) THEN
!      inviscid_dir = 33
!      loc_cx = 0
!      loc_cy = -2
!      loc_cz = 0
!    END IF
  CASE(33)
    IF (loc_state(4) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = 1
      loc_cz = -1
    else if (loc_state(20) == -11 .and. loc_state(4) == -11) THEN
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = -1
      loc_cz = -1
    else if (loc_state(18) == -11 .and. loc_state(6) == -11) THEN
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = 1
      loc_cz = 1
    else if (loc_state(18) == -11 .and. loc_state(20) == -11) THEN
      inviscid_dir = 31
      loc_cx = 0
      loc_cy = -1
      loc_cz = 1
    else if (loc_state(4) == -11) then
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = -2
      loc_cz = -1
    else if (loc_state(18) == -11) then
      inviscid_dir = 32
      loc_cx = 0
      loc_cy = -2
      loc_cz = 1
    ELSE IF (loc_state(6) == -11) THEN
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = 1
      loc_cz = 2
    ELSE IF (loc_state(20) == -11) THEN
      inviscid_dir = 24
      loc_cx = 0
      loc_cy = -1
      loc_cz = 2
    END IF
!    IF (loc_state(20) == -11 .and. loc_state(18) == -11) THEN
!      inviscid_dir = 31
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 0
!    else if (loc_state(20) == -11) then
!      inviscid_dir = 24
!      loc_cx = 0
!      loc_cy = 0
!      loc_cz = 2
!    ELSE IF (loc_state(18) == -11) THEN
!      inviscid_dir = 32
!      loc_cx = 0
!      loc_cy = -2
!      loc_cz = 0
!    END IF
!
! Far cardinal directions
!
  CASE(34)
    inviscid_dir = 37
    if (loc_state(2) == -11 .and. loc_state(16) == -11) then
      loc_cx = 2
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(16) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = -2
      loc_cy = 0
      loc_cz = 0
    end if

  CASE(35)
    inviscid_dir = 38
    if (loc_state(3) == -11 .and. loc_state(17) == -11) then
      loc_cx = 0
      loc_cy = 2
      loc_cz = 0
    else if (loc_state(17) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = 0
      loc_cy = -2
      loc_cz = 0
    end if

  CASE(36)
    inviscid_dir = 39
    if (loc_state(4) == -11 .and. loc_state(18) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 2
    else if (loc_state(18) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = 0
      loc_cy = 0
      loc_cz = -2
    end if

  CASE(37)
    inviscid_dir = 34
    if (loc_state(5) == -11 .and. loc_state(19) == -11) then
      loc_cx = -2
      loc_cy = 0
      loc_cz = 0
    else if (loc_state(19) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = 2
      loc_cy = 0
      loc_cz = 0
    end if

  CASE(38)
    inviscid_dir = 35
    if (loc_state(6) == -11 .and. loc_state(20) == -11) then
      loc_cx = 0
      loc_cy = -2
      loc_cz = 0
    else if (loc_state(20) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = 0
      loc_cy = 2
      loc_cz = 0
    end if

  CASE(39)
    inviscid_dir = 36
    if (loc_state(7) == -11 .and. loc_state(21) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = -2
    else if (loc_state(21) == -11) then
      loc_cx = 0
      loc_cy = 0
      loc_cz = 0
    else
      loc_cx = 0
      loc_cy = 0
      loc_cz = 2
    end if
  END SELECT
!END IF


if (shifted) then
  if (primary_flow_dir == 1) then
    loc_cx = loc_cx + 1
  else if (primary_flow_dir == 2) then
    loc_cy = loc_cy + 1
  else if (primary_flow_dir == 3) then
    loc_cz = loc_cz + 1
  else if (primary_flow_dir == 4) then
    loc_cx = loc_cx - 1
  else if (primary_flow_dir == 5) then
    loc_cy = loc_cy - 1
  else if (primary_flow_dir == 6) then
    loc_cz = loc_cz - 1
  end if
end if

END SUBROUTINE

!IF (dir == 9) THEN
!  SELECT CASE (i)
!!
!! Cardinal directions for the lattice, shouldn't matter which plane
!!
!CASE(2)
!  inviscid_dir = 4
!CASE(3)
!  inviscid_dir = 5
!CASE(4)
!  inviscid_dir = 2
!CASE(5)
!  inviscid_dir = 3
!!
!! Diagonal directions, must determine which side the inviscid wall/
!! symmetry plane is on.
!!
!CASE(6)
!  IF (link(2,j) == -11) THEN
!    inviscid_dir = 7
!  ELSE IF (link(3,j) == -11) THEN
!    inviscid_dir = 9
!  END IF
!CASE(7)
!  IF (link(4,j) == -11) THEN
!    inviscid_dir = 6
!  ELSE IF (link(3,j) == -11) THEN
!    inviscid_dir = 8
!  END IF
!CASE(8)
!  IF (link(4,j) == -11) THEN
!    inviscid_dir = 9
!  ELSE IF (link(5,j) == -11) THEN
!    inviscid_dir = 7
!  END IF
!CASE(9)
!  IF (link(2,j) == -11) THEN
!    inviscid_dir = 8
!  ELSE IF (link(5,j) == -11) THEN
!    inviscid_dir = 6
!  END IF
!
!END SELECT
!!
!! Inviscid directions for the D2Q21 lattice
!!
!ELSE IF (dir == 21) THEN
!  SELECT CASE (i)
!  CASE(2)
!    inviscid_dir = 4
!  CASE(3)
!    inviscid_dir = 5
!  CASE(4)
!    inviscid_dir = 2
!  CASE(5)
!    inviscid_dir = 3
!!
!! Close diagonal directions
!!
!CASE (6)
!  IF (link(2,j) == -11) THEN
!    inviscid_dir = 7
!  ELSE IF (link(3,j) == -11) THEN
!    inviscid_dir = 9
!  END IF
!
!CASE(7)
!  IF (link(4,j) == -11) THEN
!    inviscid_dir = 6
!  ELSE IF (link(3,j) == -11) THEN
!    inviscid_dir = 8
!  END IF
!CASE(8)
!  IF (link(4,j) == -11) THEN
!    inviscid_dir = 9
!  ELSE IF (link(5,j) == -11) THEN
!    inviscid_dir = 7
!  END IF
!CASE(9)
!  IF (link(2,j) == -11) THEN
!    inviscid_dir = 8
!  ELSE IF (link(5,j) == -11) THEN
!    inviscid_dir = 6
!  END IF
!!
!! Mid-range cardinal direction
!!
!CASE(10)
!  inviscid_dir = 12
!CASE(11)
!  inviscid_dir = 13
!CASE(12)
!  inviscid_dir = 10
!CASE(13)
!  inviscid_dir = 11
!!
!! Mid-range diagonal directions
!!
!CASE(14)
!  IF (link(10,j) == -11) THEN
!    inviscid_dir = 15
!  ELSE IF (link(11,j) == -11) THEN
!    inviscid_dir = 17
!  END IF
!CASE(15)
!  IF (link(12,j) == -11) THEN
!    inviscid_dir = 14
!  ELSE IF (link(11,j) == -11) THEN
!    inviscid_dir = 16
!  END IF
!CASE(16)
!  IF (link(12,j) == -11) THEN
!    inviscid_dir = 17
!  ELSE IF (link(13,j) == -11) THEN
!    inviscid_dir = 15
!  END IF
!CASE(17)
!  IF (link(10,j) == -11) THEN
!    inviscid_dir = 16
!  ELSE IF (link(13,j) == -11) THEN
!    inviscid_dir = 14
!  END IF
!!
!! Far cardinal directions
!!
!CASE(18)
!  inviscid_dir = 20
!CASE(19)
!  inviscid_dir = 21
!CASE(20)
!  inviscid_dir = 18
!CASE(21)
!  inviscid_dir = 19
!END SELECT
!
!!
!! The inviscid directions for D3Q19 lattice
!!
!ELSE IF (dir == 19) THEN
!  SELECT CASE (i)
!!
!! Basic cardinal directions
!!
!CASE(2)
!  inviscid_dir = 5
!CASE(3)
!  inviscid_dir = 6
!CASE(4)
!  inviscid_dir = 7
!CASE(5)
!  inviscid_dir = 2
!CASE(6)
!  inviscid_dir = 3
!CASE(7)
!  inviscid_dir = 4
!!
!! Diagonal directions, requiring more information on which direction
!! the inviscid wall/symmetry plane is.
!!
!  CASE(8)
!    IF (link(2,j) == -11) THEN
!      inviscid_dir = 10
!    ELSE IF (link(3,j) == -11) THEN
!      inviscid_dir = 16
!    END IF
!  CASE(9)
!    IF (link(7,j) == -11) THEN
!      inviscid_dir = 11
!    ELSE IF (link(3,j) == -11) THEN
!      inviscid_dir = 17
!    END IF
!  CASE(10)
!    IF (link(5,j) == -11) THEN
!      inviscid_dir = 8
!    ELSE IF (link(3,j) == -11) THEN
!      inviscid_dir = 18
!    END IF
!  CASE(11)
!    IF (link(4,j) == -11) THEN
!      inviscid_dir = 9
!    ELSE IF (link(3,j) == -11) THEN
!      inviscid_dir = 19
!    END IF
!  CASE(12)
!    IF (link(2,j) == -11) THEN
!      inviscid_dir = 15
!    ELSE IF (link(4,j) == -11) THEN
!      inviscid_dir = 13
!    END IF
!  CASE(13)
!    IF (link(2,j) == -11) THEN
!      inviscid_dir = 14
!    ELSE IF (link(7,j) == -11) THEN
!      inviscid_dir = 12
!    END IF
!  CASE(14)
!    IF (link(5,j) == -11) THEN
!      inviscid_dir = 13
!    ELSE IF (link(7,j) == -11) THEN
!      inviscid_dir = 15
!    END IF
!  CASE(15)
!    IF (link(5,j) == -11) THEN
!      inviscid_dir = 12
!    ELSE IF (link(4,j) == -11) THEN
!      inviscid_dir = 14
!    END IF
!  CASE(16)
!    IF (link(2,j) == -11) THEN
!      inviscid_dir = 18
!    ELSE IF (link(6,j) == -11) THEN
!      inviscid_dir = 8
!    END IF
!  CASE(17)
!    IF (link(6,j) == -11) THEN
!      inviscid_dir = 9
!    ELSE IF (link(7,j) == -11) THEN
!      inviscid_dir = 19
!    END IF
!  CASE(18)
!    IF (link(5,j) == -11) THEN
!      inviscid_dir = 16
!    ELSE IF (link(6,j) == -11) THEN
!      inviscid_dir = 10
!    END IF
!  CASE(19)
!    IF (link(4,j) == -11) THEN
!      inviscid_dir = 17
!    ELSE IF (link(6,j) == -11) THEN
!      inviscid_dir = 11
!    END IF
!
!  END SELECT
!
!ELSE IF (dir == 39) THEN
