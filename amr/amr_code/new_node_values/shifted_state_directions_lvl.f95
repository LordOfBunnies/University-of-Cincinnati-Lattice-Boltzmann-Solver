subroutine shifted_state_directions_lvl(lvl)
!
! This is the define the state of the future links. This will create some
! false positive nodes and we will have to deal with them in another
!
!
! Called by: grid_gen
! Calls: error_out
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
IMPLICIT NONE
REAL(kind=dp) :: loc_coords(3)!,min_coords(3),max_coords(3)
real(kind=dp),allocatable :: loc_dx(:)
!real(kind=dp),contiguous,pointer :: tri_catcher(:,:,:,:)
INTEGER :: a,b,c,i,j,k,pos_xvar_2d,pos_yvar_2d,pos_zvar_2d,&
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
  CALL state_directions_cgns(fine)
  RETURN
END IF

!write(*,*) inlet_side,outlet_side,amrex_geom(0)%domain%lo(1),amrex_geom(0)%domain%hi(1)

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
!  write(*,*) 'state bounds ',bot_bound,top_bound,lvl,num_inlets,num_outlets

  dir_loop: DO i = 2,dir
  z_loop: DO c = sax%lo(3)-nghosts_mid,sax%hi(3)+nghosts_mid
    y_loop: DO b = sax%lo(2)-nghosts_mid,sax%hi(2)+nghosts_mid
      x_loop: DO a = sax%lo(1)-nghosts_mid,sax%hi(1)+nghosts_mid
!          if (a == -1 .and. b == 32 .and. c == 32) then
!            write(*,*) 'exterior node values',state(a,b,c,1),a,b,c
!          end if
!z_loop: DO c =1,amrex_geom(lvl)%domain%hi(3)
!  y_loop: DO b = 1,amrex_geom(lvl)%domain%hi(2)
!    x_loop: DO a = 1,amrex_geom(lvl)%domain%hi(1)
        if (state(a,b,c,1) == -1001 .or. state(a,b,c,1) == -666) then
          state(a,b,c,i) = -666
          cycle
        end if

!          if (a == 1 .and. b == 32 .and. c == 32) then
!            write(*,*) 'state secrets',state(a+cx(i),b+cy(i),c+cz(i),1),i,state(a,b,c,1:dir),self
!          end if
!        if (.not. shifted) then
!        if (i > 15 .and. i < 22) then
!          select case(i)
!            case(16)
!              if(state(a,b,c,2) < -1000) then
!                  state(a,b,c,16) = (state(a,b,c,2) + 1000)/2
!                cycle
!              end if
!            case(17)
!              if(state(a,b,c,3) < -1000) then
!                  state(a,b,c,17) = (state(a,b,c,3) + 1000)/2
!                cycle
!              end if
!            case(18)
!              if(state(a,b,c,4) < -1000) then
!                  state(a,b,c,18) = (state(a,b,c,4) + 1000)/2
!                cycle
!              end if
!            case(19)
!              if(state(a,b,c,5) < -1000) then
!                  state(a,b,c,19) = (state(a,b,c,5) + 1000)/2
!                cycle
!              end if
!            case(20)
!              if(state(a,b,c,6) < -1000) then
!                  state(a,b,c,20) = (state(a,b,c,6) + 1000)/2
!                  !write(*,*) 'jangly'
!                cycle
!              end if
!            case(21)
!              if(state(a,b,c,7) < -1000) then
!                  state(a,b,c,21) = (state(a,b,c,7) + 1000)/2
!                cycle
!              end if
!          end select
!
!        else if (i > 33) then
!          select case(i)
!            case(34)
!              if (state(a,b,c,2) < -1000) then
!                state(a,b,c,34) = (state(a,b,c,2) + 1000)/3
!                cycle
!              else if (state(a,b,c,16) < -1000) then
!                state(a,b,c,34) = (state(a,b,c,16) + 1000)*2/3
!                cycle
!              end if
!            case(35)
!              if (state(a,b,c,3) < -1000) then
!                state(a,b,c,35) = (state(a,b,c,3) + 1000)/3
!                cycle
!              else if (state(a,b,c,17) < -1000) then
!                state(a,b,c,35) = (state(a,b,c,35) + 1000)*2/3
!                cycle
!              end if
!            case(36)
!              if (state(a,b,c,4) < -1000) then
!                state(a,b,c,36) = (state(a,b,c,4) + 1000)/3
!                cycle
!              else if (state(a,b,c,18) < -1000) then
!                state(a,b,c,36) = (state(a,b,c,18) + 1000)*2/3
!                cycle
!              end if
!            case(37)
!              if (state(a,b,c,5) < -1000) then
!                state(a,b,c,37) = (state(a,b,c,5) + 1000)/3
!                cycle
!              else if (state(a,b,c,19) < -1000) then
!                state(a,b,c,37) = (state(a,b,c,19) + 1000)*2/3
!                cycle
!              end if
!            case(38)
!              if (state(a,b,c,6) < -1000) then
!                state(a,b,c,38) = (state(a,b,c,6) + 1000)/3
!                !if (a == 110 .and. b == 66 .and. c == 20) write(*,*) 'di mana'
!                !write(*,*) 'boof goof'
!                cycle
!              else if (state(a,b,c,20) < -1000) then
!                state(a,b,c,38) = (state(a,b,c,20) + 1000)*2/3
!                !if (a == 110 .and. b == 66 .and. c == 20) write(*,*) 'di sana'
!                cycle
!              end if
!            case(39)
!              if (state(a,b,c,7) < -1000) then
!                state(a,b,c,39) = (state(a,b,c,7) + 1000)/3
!                cycle
!              else if (state(a,b,c,21) < -1000) then
!                state(a,b,c,39) = (state(a,b,c,21) + 1000)*2/3
!                cycle
!              end if
!          end select
!
!        end if
!        end if
!        if (a == 26 .and. b == 32 .and. c == -3) then
!          write(*,*) 'out of bounds nodes',state(a,b,c,1:dir),a,b,c
!        end if
!          loc_coords = get_real_coords(amrex_geom(lvl),(/a,b,c/))
!
! In the -x direction
!
          IF (a-cx(i) < amrex_geom(lvl)%domain%lo(1)) THEN
!                  WRITE(*,*) 'Out of bounds in -x dir.',a,b,c
            if (a-cx(i) < bot_bound(1)) then
              state(a,b,c,i) = -666
!              write(*,*) 'spankers'
              cycle
            end if
! For inlets
            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000 .and. .not. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
!              write(*,*) 'tickle me elmo'
              cycle
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            end if

            IF (inlet_side(1)) THEN
              DO j = 1,num_inlets
                IF (inlet_type(2,j) == 1) THEN
                      !WRITE(*,*) 'Inlet in positive x-dir at boundary.'
                  IF (inlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1

                  ELSE IF (inlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (inlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -100 - j
!                    WRITE(*,*) 'Assigning inlet number ',state(a,b,c,i),i,a,b,c
                  END IF
!                  do k =2,dir
!                    if (cx(k) > 0 .and. primary_flow_dir == 1 .and. shifted) then
!                      if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                        state(a,b,c,k) = -300 - j
!                        !write(*,*) 'snoofle',a,b,c
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -300 - j
!                      end if
!                    end if
!                  end do
                END IF
              END DO

!              if (a == 1 .and. b == 16 .and. c == 20) then
!                write(*,*) 'bling bling',state(a,b,c,1:dir)
!              end if
! For outlets
            ELSE IF (outlet_side(1)) THEN
              DO j = 1,num_outlets
                IF (outlet_type(2,j) == 1) THEN
                  IF (outlet_type(1,j) == 10) THEN
                    state(a,b,c,i) = -10
!                          num_solid_dirs = num_solid_dirs + 1
                  else if (outlet_type(1,j) == 6) then
                    state(a,b,c,i) = state(a+cx(i),b+cy(i),c+cz(i),1)
                  ELSE IF (outlet_type(1,j) == 11) THEN
                    state(a,b,c,i) = -11
                  ELSE IF (outlet_type(1,j) == 8) THEN
                    state(a,b,c,i) = -1000
                  ELSE
                    state(a,b,c,i) = -200 - j
                  END IF

!                  do k =2,dir
!                    if (cx(k) > 0 .and. primary_flow_dir == 1 .and. shifted) then
!                      if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                        state(a,b,c,k) = -400 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -400 - j
!
!                      end if
!                    end if
!                  end do
                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
!              write(*,*) 'the fuck is going on.',a,b,c
            END IF
!
! In the -y direction
!
          ELSE IF (b-cy(i) < amrex_geom(lvl)%domain%lo(2)) THEN

            if (b-cy(i) < bot_bound(2)) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000 .and. .not. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
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
!                  do k =2,dir
!                    if (cy(k) > 0 .and. primary_flow_dir == 2 .and. shifted) then
!                      if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                        state(a,b,c,k) = -300 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -300 - j
!
!                      end if
!                    end if
!                  end do

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

!                  do k =2,dir
!                    if (cy(k) > 0 .and. primary_flow_dir == 2 .and. shifted) then
!                      if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                        state(a,b,c,k) = -400 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -400 - j
!                      end if
!                    end if
!                  end do

                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the -z direction
!
          ELSE IF (c-cz(i) < amrex_geom(lvl)%domain%lo(3)) THEN

            if (c-cz(i) < bot_bound(3)) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000 .and. .not. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
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

!                  do k =2,dir
!                    if (cz(k) > 0 .and. primary_flow_dir == 3 .and. shifted) then
!                      if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                        state(a,b,c,k) = -300 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -300 - j
!
!                      end if
!                    end if
!                  end do

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

!                  do k =2,dir
!                    if (cz(k) > 0 .and. primary_flow_dir == 3 .and. shifted) then
!                      if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                        state(a,b,c,k) = -400 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -400 - j
!
!                      end if
!                    end if
!                  end do

                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the +x direction
!
          ELSE IF (a-cx(i) > amrex_geom(lvl)%domain%hi(1)+1) THEN

            if (a-cx(i) >= top_bound(1)+1) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000 .and. .not. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
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

!                  do k =2,dir
!                    if (cx(k) < 0 .and. primary_flow_dir == 4 .and. shifted) then
!                      if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                        state(a,b,c,k) = -300 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -300 - j
!                      end if
!                    end if
!                  end do

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

!                  do k =2,dir
!                    if (cx(k) < 0 .and. primary_flow_dir == 4 .and. shifted) then
!                      if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                        state(a,b,c,k) = -400 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -400 - j
!                      end if
!                    end if
!                  end do

                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
            END IF
!
! In the +y direction
!
          ELSE IF (b-cy(i) > amrex_geom(lvl)%domain%hi(2)+1) THEN

            if (b-cy(i) >= top_bound(2)+1) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000 .and. .not. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
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

!                  do k =2,dir
!                    if (cy(k) < 0 .and. primary_flow_dir == 5 .and. shifted) then
!                      if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                        state(a,b,c,k) = -300 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -300 - j
!                      end if
!                    end if
!                  end do

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

!                  do k =2,dir
!                    if (cy(k) < 0 .and. primary_flow_dir == 5 .and. shifted) then
!                      if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                        state(a,b,c,k) = -400 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -400 - j
!                      end if
!                    end if
!                  end do

                END IF
              END DO
            ELSE

              state(a,b,c,i) = -1000
            END IF
!
! In the +z direction
!
          ELSE IF (c-cz(i) > amrex_geom(lvl)%domain%hi(3)+1) THEN

            if (c-cz(i) >= top_bound(3)+1) then
              state(a,b,c,i) = -666
              cycle
            end if

            if (state(a+cx(i),b+cy(i),c+cz(i),1) <-1000 .and. .not. shifted) then
              !write(*,*) 'boogie'
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
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

!                  do k =2,dir
!                    if (cz(k) < 0 .and. primary_flow_dir == 6 .and. shifted) then
!                      if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                        state(a,b,c,k) = -300 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -300 - j
!                      end if
!                    end if
!                  end do

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

!                  do k =2,dir
!                    if (cz(k) < 0 .and. primary_flow_dir == 6 .and. shifted) then
!                      if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                        state(a,b,c,k) = -400 - j
!                      else if (state(a,b,c,i) == -11 .or. state(a,b,c,i) == -1000 .or. &
!                          state(a,b,c,i) == -10) then
!                        state(a,b,c,k) = -400 - j
!                      end if
!                    end if
!                  end do

                END IF
              END DO
            ELSE
              state(a,b,c,i) = -1000
              !write(*,*) 'woogie',state(a,b,c,i),state(a,b,c,1),state(a+cx(i),b+cy(i),c+cz(i),1),a,b,c,i
            END IF
!
! In the case of ghost nodes where the direction would leave all the nodes on that level
!
          else if (a-cx(i) < bot_bound(1) .or. a-cx(i) > top_bound(1) &
            .or. b-cy(i) < bot_bound(2) .or. b-cy(i) > top_bound(2) &
            .or. c-cz(i) < bot_bound(3) .or. c-cz(i) > top_bound(3)) then

            state(a,b,c,i) = -666
!
! If there's a surface
!
          ELSE
            IF (state(a+cx(i),b+cy(i),c+cz(i),1) == -1001 .and. .not. shifted) THEN
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
            else if (state(a-cx(i),b-cy(i),c-cz(i),1) <-1000 .and. shifted) then
              call inter_bb_setup(state(a,b,c,i),state(a,b,c,1),i,a,b,c,lvl)
              cycle
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
!          if (a == 1 .and. b == 16 .and. c == 20) then
!            write(*,*) 'state secrets',state(a+cx(i),b+cy(i),c+cz(i),1),i,state(a,b,c,1:dir),self
!          end if
!          if (a == 40 .and. b == 30 .and. c == 16 .and. i == 39) then
!            write(*,*) 'How the fluff?',state(a,b,c,1:39),sax%lo,sax%hi,a,b,c,lvl,self
!          end if
!          if (a == 1 .and. b == 1 .and. c == -1 .and. lvl == 0) then

!          end if
        END DO x_loop
      END DO y_loop
    END DO z_loop
  END DO dir_loop
!
! For shifted lattices, make sure the directions are accounted for.
!
!  if (shifted) then
!  DO i = 2,dir
!    DO c = sax%lo(3)-nghosts_mid,sax%hi(3)+nghosts_mid
!      DO b = sax%lo(2)-nghosts_mid,sax%hi(2)+nghosts_mid
!        DO a = sax%lo(1)-nghosts_mid,sax%hi(1)+nghosts_mid
!
!!
!          if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199 .and. &
!                primary_flow_dir ==1 .and. cx(i) < 0) then
!
!            do k = 2,dir
!              if (cx(k) > 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!          else if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199 .and. &
!                primary_flow_dir ==2 .and. cy(i) < 0) then
!
!            do k = 2,dir
!              if (cy(k) > 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          else if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199 .and. &
!              primary_flow_dir == 3 .and. cz(i) < 0) then
!
!          do k = 2,dir
!            if (cz(k) > 0) then
!              state(a,b,c,k) = state(a,b,c,i) -200
!            end if
!          end do
!
!          else if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199 .and. &
!                primary_flow_dir ==4 .and. cx(i) > 0) then
!
!            do k = 2,dir
!              if (cx(k) < 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!            else if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199 .and. &
!                primary_flow_dir ==5 .and. cy(i) > 0) then
!
!            do k = 2,dir
!              if (cy(k) < 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!            else if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199 .and. &
!                primary_flow_dir == 6 .and. cz(i) > 0) then
!
!            do k = 2,dir
!              if (cz(k) < 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          end if
!!
!          if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299 .and. &
!                primary_flow_dir ==1 .and. cx(i) > 0) then
!
!            do k = 2,dir
!              if (cx(k) < 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!          else if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299 .and. &
!                primary_flow_dir ==2 .and. cy(i) > 0) then
!
!            do k = 2,dir
!              if (cy(k) < 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          else if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299 .and. &
!              primary_flow_dir == 3 .and. cz(i) > 0) then
!
!            do k = 2,dir
!              if (cz(k) < 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          else if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299 .and. &
!                primary_flow_dir ==4 .and. cx(i) < 0) then
!
!            do k = 2,dir
!              if (cx(k) > 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          else if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299 .and. &
!                primary_flow_dir ==5 .and. cy(i) < 0) then
!
!            do k = 2,dir
!              if (cy(k) > 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          else if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299 .and. &
!                primary_flow_dir == 6 .and. cz(i) < 0) then
!
!            do k = 2,dir
!              if (cz(k) > 0) then
!                state(a,b,c,k) = state(a,b,c,i) -200
!              end if
!            end do
!
!          end if
!
!
!        end do
!      end do
!    end do
!  end do
!
!  end if


  end do !mfiter loop
  call amrex_mfiter_destroy(mfi)

!  call amrex_multifab_destroy(mf_tri_catcher)
!end do !level loop

!WRITE(11,101) num_solid_dirs
!WRITE(11,102) num_card_solid_dirs
!WRITE(*,101) num_solid_dirs
!WRITE(*,102) num_card_solid_dirs

! 101  FORMAT('Total number of solid directions ',I8)
! 102  FORMAT('Number of cardinal solid directions ',I8)

END SUBROUTINE
!          if (a == 1 .and. b == 31 .and. c == 1) then
!            write(*,*) 'state secrets',i,state(a,b,c,i)
!          end if
!      whelp = 1
          !if (state(a,b,c,i) == 0 .and. lvl /= 0) write(*,*) 'RED ALERT!!!',a,b,c,i,lvl
                     !if (state(a,b,c,i) == 0 .and. lvl /= 0) write(*,*) 'RED ALERT!!!',a,b,c,i,lvl
!              if (state(a,b,c,i) == 0 .and. lvl /= 0) write(*,*) 'RED ALERT!!!',state(a,b,c,1),a,b,c,i,&
!                lvl,sax%lo-nghosts,sax%hi+nghosts
!            write(*,*) 'null direction at edge of level',a+cx(i),b+cy(i),c+cz(i),&
!              bot_bound,top_bound
!      WRITE(11,*) state(1,1,17,17),state(2,1,17,17),state(3,1,17,17),state(4,1,17,17),&
!        state(5,1,17,17),state(6,1,17,17),state(7,1,17,17),state(8,1,17,17),state(9,1,17,17),&
!        state(10,1,17,17),state(11,1,17,17),state(12,1,17,17),state(13,1,17,17),state(14,1,17,17),&
!        state(16,1,17,17),state(17,1,17,17),state(18,1,17,17),state(19,1,17,17),state(20,1,17,17)



!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1001
!              num_solid_dirs = num_solid_dirs + 1


!              IF (i >= 2 .AND. i <= 7) THEN
!                num_card_solid_dirs = num_card_solid_dirs + 1
!                    ELSE IF (i == 3) THEN
!                      num_card_solid_dirs = num_card_solid_dirs + 1
!                    ELSE IF (i == 4) THEN
!                      num_card_solid_dirs = num_card_solid_dirs + 1
!                    ELSE IF (i == 5) THEN
!                      num_card_solid_dirs = num_card_solid_dirs + 1
!                    ELSE IF (i == 6) THEN
!                      num_card_solid_dirs = num_card_solid_dirs + 1
!                    ELSE IF (i == 7) THEN
!                      num_card_solid_dirs = num_card_solid_dirs + 1
!              END IF





!
! Rectangular bounding boxes
!
!        IF (bounding_method >=-3 .AND. bounding_method <= 3) THEN
!          skipper = .FALSE.
!
!
!
! Spherical bounding boxes
!
!
!        ELSE IF (bounding_method == 11 .OR. &
!             bounding_method == 12) THEN
!          IF (a+cx(i) < 1 .OR. state(1,a+cx(i),b+cy(i),c+cz(i))&
!                   == -1001 .AND. state(a,b,c,1) /= -1) THEN
!            IF (cx(i) < 0) THEN
!! For inlets
!              DO j = 1,num_inlets
!                IF (inlet_type(2,j) == 1) THEN
!                  IF (inlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (inlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!! For outlets
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 1) THEN
!                  IF (outlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (outlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (outlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!            END IF
!!
!! In the positive y direction
!!
!          ELSE IF (b+cy(i) < 1 .OR. state(a+cx(i),b+cy(i),c+cz(i),i)&
!                     == -1 .AND. state(a,b,c,1) /= -1) THEN
!            IF (cy(i) < 0) THEN
!! For inlets
!              DO j = 1,num_inlets
!                IF (inlet_type(2,j) == 2) THEN
!                  IF (inlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (inlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!! For outlets
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 2) THEN
!                  IF (outlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (outlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (outlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!            END IF
!!
!! In the positive z-direction
!!
!          ELSE IF (c+cz(i) < 1 .OR. state(a+cx(i),b+cy(i),c+cz(i),i)&
!                   == -1 .AND. state(a,b,c,1) /= -1) THEN
!            IF (cz(i) < 0) THEN
!! For inlets
!              DO j = 1,num_inlets
!                IF (inlet_type(2,j) == 3) THEN
!                  IF (inlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (inlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!! For outlets
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 3) THEN
!                  IF (outlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (outlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (outlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!            END IF
!!
!! Check on the negative x direction stuff
!!
!          ELSE IF (a+cx(i) > amrex_geom(lvl)%domain%hi(1) .OR. state(a+cx(i),b+cy(i),c+cz(i),i)&
!                   == -1 .AND. state(a,b,c,1) /= -1) THEN
!            IF (cx(i) > 0) THEN
!! For inlets
!              DO j = 1,num_inlets
!                IF (inlet_type(2,j) == 4) THEN
!                  IF (inlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    !num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (inlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!! For outlets
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 4) THEN
!                  IF (outlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    !num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (outlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (outlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!            END IF
!!
!! Check on the negative y direction stuff
!!
!          ELSE IF (b+cy(i) > amrex_geom(lvl)%domain%hi(2) .OR. state(a+cx(i),b+cy(i),c+cz(i),1)&
!                   == -1 .AND. state(a,b,c,1) /= -1) THEN
!            IF (cy(i) > 0) THEN
!! For inlets
!              DO j = 1,num_inlets
!                IF (inlet_type(2,j) == 5) THEN
!                  IF (inlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    !num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (inlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!! For outlets
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 5) THEN
!                  IF (outlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    !num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (outlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (outlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!            END IF
!!
!! Check on the negative z direction stuff
!!
!          ELSE IF (c+cz(i) > amrex_geom(lvl)%domain%hi(3) .OR. state(a+cx(i),b+cy(i),c+cz(i),1)&
!                   == -1 .AND. state(a,b,c,1) /= -1) THEN
!            IF (cz(i) > 0) THEN
!! For inlets
!              DO j = 1,num_inlets
!                IF (inlet_type(2,j) == 6) THEN
!                  IF (inlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    !num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (inlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!! For outlets
!              DO j = 1,num_outlets
!                IF (outlet_type(2,j) == 6) THEN
!                  IF (outlet_type(1,j) == 10) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                    !num_solid_dirs = num_solid_dirs + 1
!                  ELSE IF (outlet_type(1,j) == 11) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                  ELSE IF (outlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                  ELSE
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                  END IF
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                END IF
!              END DO
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!            END IF
!!
!!
!! This says what the state of the node in every direction is. Very
!! important for the linking portion later
!!
!          ELSE
!            IF (state(a+cx(i),b+cy(i),c+cz(i),1) == -1001) THEN
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1001
!              !num_solid_dirs = num_solid_dirs + 1
!            ELSE IF (state(a+cx(i),b+cy(i),c+cz(i),1) == -666) THEN
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -666
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = state(a+cx(i),b+cy(i),c+cz(i),1)
!            END IF
!          END IF
!!
!!
!! Cylindrical bounding boxes
!!
!!
!        ELSE IF (bounding_method >=21 .OR. bounding_method <=26) THEN
!          IF (a+cx(i) < amrex_geom(lvl)%domain%lo(1)) THEN
!            DO j = 1,num_inlets
!              IF (inlet_type(2,j) == 1) THEN
!                IF (inlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (inlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (inlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!! For outlets
!            DO j = 1,num_outlets
!              IF (outlet_type(2,j) == 1) THEN
!                IF (outlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (outlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (outlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!
!          ELSE IF (b+cy(i) < amrex_geom(lvl)%domain%lo(2)) THEN
!            DO j = 1,num_inlets
!              IF (inlet_type(2,j) == 2) THEN
!                IF (inlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (inlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (inlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!! For outlets
!            DO j = 1,num_outlets
!              IF (outlet_type(2,j) == 2) THEN
!                IF (outlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (outlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (outlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!
!          ELSE IF (c+cz(i) < amrex_geom(lvl)%domain%lo(3)) THEN
!            DO j = 1,num_inlets
!              IF (inlet_type(2,j) == 3) THEN
!                IF (inlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (inlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (inlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!! For outlets
!            DO j = 1,num_outlets
!              IF (outlet_type(2,j) == 3) THEN
!                IF (outlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (outlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (outlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!
!          ELSE IF (a+cx(i) > amrex_geom(lvl)%domain%hi(1)) THEN
!            DO j = 1,num_inlets
!              IF (inlet_type(2,j) == 4) THEN
!                IF (inlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (inlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (inlet_type(1,j) == 8) THEN
!                    state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!! For outlets
!            DO j = 1,num_outlets
!              IF (outlet_type(2,j) == 5) THEN
!                IF (outlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (outlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (outlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!
!          ELSE IF (b+cy(i) > amrex_geom(lvl)%domain%hi(2)) THEN
!            DO j = 1,num_inlets
!              IF (inlet_type(2,j) == 5) THEN
!                IF (inlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (inlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (inlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!! For outlets
!            DO j = 1,num_outlets
!              IF (outlet_type(2,j) == 5) THEN
!                IF (outlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (outlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (outlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!
!          ELSE IF (c+cz(i) > amrex_geom(lvl)%domain%hi(3)) THEN
!            DO j = 1,num_inlets
!              IF (inlet_type(2,j) == 6) THEN
!                IF (inlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (inlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (inlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -100 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!! For outlets
!            DO j = 1,num_outlets
!              IF (outlet_type(2,j) == 6) THEN
!                IF (outlet_type(1,j) == 10) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -10
!                  !num_solid_dirs = num_solid_dirs + 1
!                ELSE IF (outlet_type(1,j) == 11) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -11
!                ELSE IF (outlet_type(1,j) == 8) THEN
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!                ELSE
!                  state(a+cx(i),b+cy(i),c+cz(i),1) = -200 - j
!                END IF
!              ELSE
!                state(a+cx(i),b+cy(i),c+cz(i),1) = -1000
!              END IF
!            END DO
!          ELSE
!            IF (state(1,a+cx(i),b+cy(i),c+cz(i)) == -1001) THEN
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -1001
!              !num_solid_dirs = num_solid_dirs + 1
!            ELSE IF (state(1,a+cx(i),b+cy(i),c+cz(i)) == -666) THEN
!              state(a+cx(i),b+cy(i),c+cz(i),1) = -666
!            ELSE
!              state(a+cx(i),b+cy(i),c+cz(i),1) = state(a+cx(i),b+cy(i),c+cz(i),1)
!              if (state(a+cx(i),b+cy(i),c+cz(i),1) == -1001 and state(a,b,c,1) >=0) then
!                call inter_bb_setup(state(a+cx(i),b+cy(i),c+cz(i),1),state(a,b,c,1),i,a,b,c,lvl)
!              end if
!            END IF
!          END IF
!        END IF
