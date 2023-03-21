subroutine state_directions_lvl(lvl)
!
!
!
!
!
!
!
!
!
use grid_data
use linkwise
use nml_inlet_outlet
!use nml_bounding_box
use precise
use constants
use startup
use amrex_amr_module
use amrex_base_module
use amr_processes, only : get_real_coords,node_based_3d,self
use amr_info_holder, only: state,mfstate,nghosts,nghosts_mid
implicit none

REAL(kind=dp) :: loc_coords(3)!,min_coords(3),max_coords(3)
real(kind=dp),allocatable :: loc_dx(:)
!real(kind=dp),contiguous,pointer :: tri_catcher(:,:,:,:)
INTEGER :: a,b,c,i,j,pos_xvar_2d,pos_yvar_2d,pos_zvar_2d,&
  neg_xvar_2d,neg_yvar_2d,neg_zvar_2d,lvl,fine,bot_bound(3),top_bound(3)
!integer :: lvl
!LOGICAL :: skipper

type(amrex_mfiter) :: mfi
type(amrex_box) :: sax
!type(amrex_multifab) :: mf_tri_catcher

allocate(loc_dx(0:amrex_max_level))

num_solid_dirs = 0
num_card_solid_dirs = 0

IF (cgns_start .OR. mixed_start) THEN
  CALL state_directions_cgns
  RETURN
END IF

pos_xvar_2d = 2
neg_xvar_2d = 4
pos_zvar_2d = 3
neg_zvar_2d = 5
!IF (bounding_method == -1) THEN
!pos_yvar_2d = 3
!neg_yvar_2d = 5
!ELSE IF (bounding_method == -2) THEN
pos_yvar_2d = 0
neg_yvar_2d = 0
!ELSE IF (bounding_method == -3) THEN
!pos_yvar_2d = 2
!neg_yvar_2d = 4
!END IF

!do lvl = fine,0,-1
  !call amrex_multifab_build(mf_tri_catcher,mfstate(lvl)%ba,mfstate(lvl)%dm,1,nghosts,node_based_3d)

  CALL amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)
  do while(mfi%next())
  state => mfstate(lvl)%dataptr(mfi)
!  tri_catcher => mf_tri_catcher%dataptr(mfi)

  sax = mfi%validbox()
  bot_bound = sax%lo-nghosts_mid
  top_bound = sax%hi+nghosts_mid

!  write(*,*) 'grid bounds',amrex_geom(lvl)%domain%lo,amrex_geom(lvl)%domain%hi,lvl
!  write(*,*) 'state bounds ',bot_bound,top_bound,lvl

  dir_loop: DO i = 2,dir
  z_loop: DO c = sax%lo(3)-nghosts_mid,sax%hi(3)+nghosts_mid
    y_loop: DO b = sax%lo(2)-nghosts_mid,sax%hi(2)+nghosts_mid
      x_loop: DO a = sax%lo(1)-nghosts_mid,sax%hi(1)+nghosts_mid
!z_loop: DO c =1,amrex_geom(lvl)%domain%hi(3)
!  y_loop: DO b = 1,amrex_geom(lvl)%domain%hi(2)
!    x_loop: DO a = 1,amrex_geom(lvl)%domain%hi(1)
!        if (c == -1 .and. lvl == 3 .and. i == 39) then
!            write(*,*) 'state secrets',state(a+cx(i),b+cy(i),c+cz(i),1),i,state(a,b,c,1:dir),self
!          end if
        if (state(a,b,c,1) == -1001 .or. state(a,b,c,1) == -666) then
          state(a,b,c,i) = -666
          cycle
        end if
!          loc_coords = get_real_coords(amrex_geom(lvl),(/a,b,c/))
!
! In the -x direction
!
          IF (a+cx(i) < amrex_geom(lvl)%domain%lo(1)) THEN
!                  WRITE(*,*) 'Out of bounds in -x dir.'
            if (a+cx(i) < bot_bound(1)) then
              state(a,b,c,i) = -666
              cycle
            end if
! For inlets
            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if

            IF (inlet_side(1)) THEN
              DO j = 1,num_inlets
                IF (inlet_type(2,j) == 1) THEN
!                      WRITE(*,*) 'Inlet in positive x-dir at boundary.'
                  IF (inlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1

                  ELSE IF (inlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -100 - j
!                          WRITE(11,*) 'Assigning inlet number ',state(a,b,c,i),i,a,b,c
                  END IF
                END IF
              END DO
! For outlets
            ELSE IF (outlet_side(1)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 1) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
!                    if () then
!
!                    else
!                      state(a,b,c,i) = 200+j
!                    end if
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF
                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the -y direction
!
          ELSE IF (b+cy(i) < amrex_geom(lvl)%domain%lo(2)) THEN

            if (b+cy(i) < bot_bound(2)) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if

! For inlets
            IF (inlet_side(2)) THEN
              DO j = 1,num_inlets
                IF (inlet_type(2,j) == 2) THEN
                  IF (inlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1

                  ELSE IF (inlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -100 - j
                  END IF
                END IF
              END DO
! For outlets
            ELSE IF (outlet_side(2)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 2) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
                    !state(a,b,c,i) = 200+j
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF
                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the -z direction
!
          ELSE IF (c+cz(i) < amrex_geom(lvl)%domain%lo(3)) THEN

            if (c+cz(i) < bot_bound(3)) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if
! For inlets
            IF (inlet_side(3))THEN
            DO j = 1,num_inlets
              IF (inlet_type(2,j) == 3) THEN
                IF (inlet_type(1,j) == 10) THEN
                  state(a,b,c,i) = -10
!                        num_solid_dirs = num_solid_dirs + 1

                ELSE IF (inlet_type(1,j) == 11) THEN
                  state(a,b,c,i) = -11
                ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                ELSE
                  state(a,b,c,i) = -100 - j
                END IF
              END IF
            END DO
! For outlets
            ELSE IF (outlet_side(3)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 3) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
                    !state(a,b,c,i) = 200+j
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF
                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the +x direction
!
          ELSE IF (a+cx(i) > amrex_geom(lvl)%domain%hi(1)+1) THEN

            if (a+cx(i) >= top_bound(1)+1) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if
! For inlets
            IF (inlet_side(4)) THEN
              DO j = 1,num_inlets
                IF (inlet_type(2,j) == 4) THEN
                  IF (inlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1

                  ELSE IF (inlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -100 - j
                  END IF
                END IF
              END DO
! For outlets
            ELSE IF (outlet_side(4)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 4) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
                    !state(a,b,c,i) = 200+j
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF
                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the +y direction
!
          ELSE IF (b+cy(i) > amrex_geom(lvl)%domain%hi(2)+1) THEN

            if (b+cy(i) >= top_bound(2)+1) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if
! For inlets
            IF (inlet_side(5)) THEN
              DO j = 1,num_inlets
                IF (inlet_type(2,j) == 5) THEN
                  IF (inlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1

                  ELSE IF (inlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -100 - j
                  END IF
                END IF
              END DO
! For outlets
            ELSE IF (outlet_side(5)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 5) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
                    !state(a,b,c,i) = 200+j
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF
                END IF
              END DO
            ELSE

              state(a,b,c,i) = -1000
            END IF
!
! In the +z direction
!
          ELSE IF (c+cz(i) > amrex_geom(lvl)%domain%hi(3)+1) THEN

            if (c+cz(i) >= top_bound(3)+1) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000) then
              !write(*,*) 'boogie'
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if
! For inlets
            IF (inlet_side(6)) THEN
              DO j = 1,num_inlets
                IF (inlet_type(2,j) == 6) THEN
                  IF (inlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1

                  ELSE IF (inlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -100 - j
                  END IF
                END IF
              END DO
! For outlets
            ELSE IF (outlet_side(6)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 6) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
                     !num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
                    !state(a,b,c,i) = 200+j
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF
                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
              !write(*,*) 'woogie',state(a,b,c,i),state(a,b,c,1),state(a+cx(i),b+cy(i),c+cz(i),1),a,b,c,i
            END IF
!
! In the case of ghost nodes where the direction would leave all the nodes on that level
!
          else if (a+cx(i) < bot_bound(1) .or. a+cx(i) > top_bound(1)+1 &
            .or. b+cy(i) < bot_bound(2) .or. b+cy(i) > top_bound(2)+1 &
            .or. c+cz(i) < bot_bound(3) .or. c+cz(i) > top_bound(3)+1) then

            state(a,b,c,i) = -666
!
! If there's a surface
!
          ELSE
            IF (state(a+cx(i),b+cy(i),c+cz(i),1) == -1001) THEN
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
!
! In the case of null nodes
!
            ELSE IF (state(a+cx(i),b+cy(i),c+cz(i),1) == -666) THEN
              state(a,b,c,i) = -666
!
! Normal case
!
            ELSE
              state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
            END IF
          END IF

!          if (a == 108 .and. b == 72 .and. c == 35) then
!            write(*,*) 'checking states',state(a,b,c,1),state(a,b,c,i),&
!              state(a+cx(i),b+cy(i),c+cz(i),1),a,b,c,lvl,self
!          end if



        END DO x_loop
      END DO y_loop
    END DO z_loop
  END DO dir_loop


  end do !mfiter loop
  call amrex_mfiter_destroy(mfi)

!  call amrex_multifab_destroy(mf_tri_catcher)
!end do !level loop

end subroutine
