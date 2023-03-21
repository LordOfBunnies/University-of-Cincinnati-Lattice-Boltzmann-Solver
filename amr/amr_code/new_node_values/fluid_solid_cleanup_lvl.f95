subroutine fluid_solid_cleanup_lvl(lvl)
!
! This is used to help clean up extraneous elements
!
! Called by: fluid_or_solid
! Calls:
!
USE precise
USE constants
USE grid_data
USE geom_data
USE linkwise
use amrex_amr_module
use amrex_base_module
use amr_info_holder, only: mfstate,state,nghosts,nghosts_mid
use mpi
use amr_processes, only: node_based_3d,node_based_2d,self
IMPLICIT NONE
INTEGER :: i,j,k,a,skip_counter,solid_counter,&
  weird_counter,fluid_counter,null_counter,lvl,fine
integer :: n_dir,state_bounds_lo(4),state_bounds_hi(4)
integer :: box_lbound(3),box_ubound(3)
integer :: total_solid_convert_counter,total_fluid_convert_counter
!integer,contiguous,pointer :: state(:,:,:,:)
!integer :: state(bot_bound(1):top_bound(1),bot_bound(2):top_bound(2),&
!  bot_bound(3):top_bound(3)),bot_bound(4),top_bound(4)
LOGICAL :: skips(7)

type(amrex_mfiter) :: mfi
type(amrex_box) :: pox

if (dimensions == 3) then
  n_dir = 6
else if (dimensions == 2) then
  n_dir = 4
end if
!
! Go through all of the nodes to check if any are wildly out of place
!

!do lvl = 0,fine
  call amrex_mfiter_build(mfi,mfstate(lvl),tiling=.false.)

!  write(*,*) 'Fluid/solid cleanup on ',lvl,' processor ',self

  do while (mfi%next())
    state => mfstate(lvl)%dataptr(mfi)
    pox = mfi%validbox()

    box_lbound = pox%lo-nghosts_mid
    box_ubound = pox%hi+nghosts_mid
    state_bounds_lo = lbound(state)
    state_bounds_hi = ubound(state)
!    write(*,*) 'state lower bounds',state_bounds_lo,'state upper bounds',state_bounds_hi,&
!      'at level ',lvl
!    write(*,*) 'box coords',pox%lo,pox%hi
!    write(*,*) 'box bounds',box_lbound,box_ubound,'at level',lvl,self


    DO k = pox%lo(3)-nghosts_mid,pox%hi(3)+nghosts_mid
      DO j = pox%lo(2)-nghosts_mid,pox%hi(2)+nghosts_mid
        DO i = pox%lo(1)-nghosts_mid,pox%hi(1)+nghosts_mid

      skips = .FALSE.
      skip_counter = 0
      solid_counter = 0
      fluid_counter = 0
      weird_counter = 0
      null_counter = 0

      DO a = 2,7
!
! Test the feasability of using all the possible directions
! Test the +x direction
!
        IF (i + cx_27(a) < state_bounds_lo(1) .or. i+cx_27(a) < box_lbound(1)) THEN
          skips(5) = .TRUE.
          skip_counter = skip_counter + 1
        ELSE IF (i + cx_27(a) > state_bounds_hi(1) .or. i+cx_27(a) > box_ubound(1)) THEN
          skips(2) = .TRUE.
          skip_counter = skip_counter + 1
        ELSE IF ( j + cy_27(a) <  state_bounds_lo(2) .or. j+cy_27(a) < box_lbound(2)) THEN
          skips(6) = .TRUE.
          skip_counter = skip_counter + 1
        ELSE IF ( j + cy_27(a) > state_bounds_hi(2) .or. j+cy_27(a) > box_ubound(2)) THEN
          skips(3) = .TRUE.
          skip_counter = skip_counter + 1
        ELSE IF ( k + cz_27(a) <  state_bounds_lo(3) .or. k+cz_27(a) < box_lbound(3)) THEN
          skips(7) = .TRUE.
          skip_counter = skip_counter + 1
        ELSE IF ( k + cz_27(a) > state_bounds_hi(3) .or. k+cz_27(a) > box_ubound(3)) THEN
          skips(4) = .TRUE.
          skip_counter = skip_counter + 1
        END IF
      END DO
!
! If the direction is usable, test the possibilities
!
        IF (.NOT. skips(2)) THEN
          IF (state(i+cx_27(2),j,k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i+cx_27(2),j,k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i+cx_27(2),j,k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF
        END IF

!
! Test the -x direction
!
        IF (.NOT. skips(5)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i+cx_27(5),j,k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i+cx_27(5),j,k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i+cx_27(5),j,k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF

        END IF
!
! Test the +y direction
!
        IF (.NOT. skips(3)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j+cy_27(3),k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j+cy_27(3),k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j+cy_27(3),k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF

        END IF
!
! Test the -y direction
!
        IF (.NOT. skips(6)) THEN
          !write(*,*) 'coords',i,j,k,lvl
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j+cy_27(6),k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j+cy_27(6),k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j+cy_27(6),k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF

        END IF
!
! Test the +z direction
!
        IF (.NOT. skips(4)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j,k+cz_27(4),1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j,k+cz_27(4),1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j,k+cz_27(4),1) == -1) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF
        END IF
!
! Test the -z direction
!
        IF (.NOT. skips(7)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j,k+cz_27(7),1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j,k+cz_27(7),1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j,k+cz_27(7),1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF
        END IF


            IF (state(i,j,k,1) == -1001 .AND. fluid_counter >= 5) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,101)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter >= 5) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,100)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter == 4 .AND.&
              skip_counter == 1) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,102)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter == 3 .AND.&
              skip_counter == 2) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,102)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter == 3 .AND.&
              skip_counter == 3) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,102)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 4 .AND.&
              skip_counter == 1) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,103)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 3 .AND.&
              skip_counter == 2) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,103)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 3 .AND.&
              skip_counter == 3) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,103)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 6) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,107)
            ELSE IF (state(i,j,k,1) == -666 .AND. solid_counter == 6) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,106)
            ELSE IF (state(i,j,k,1) == -1001 .AND. null_counter >= 5) THEN
              state(i,j,k,1) = -666
!              WRITE(11,104)
            ELSE IF (state(i,j,k,1) >= 0 .AND. null_counter >= 5) THEN
              state(i,j,k,1) = -666
!              WRITE(11,105)
            ELSE IF (state(i,j,k,1) == -1001 .AND. null_counter >= 4 .AND.&
              skip_counter == 1) THEN
              state(i,j,k,1) = -666
!              WRITE(11,104)
            ELSE IF (state(i,j,k,1) >= 0 .AND. null_counter >= 4 .AND. &
              skip_counter == 1) THEN
              state(i,j,k,1) = -666
!              WRITE(11,105)
            ELSE IF (state(i,j,k,1) == -1001 .AND. null_counter >= 3 .AND. &
              skip_counter == 2) THEN
              state(i,j,k,1) = -666
!              WRITE(11,104)
            ELSE IF (state(i,j,k,1) >= 0 .AND. null_counter >= 3 .AND. &
              skip_counter == 2) THEN
              state(i,j,k,1) = -666
!              WRITE(11,105)


            END IF
          END DO
        END DO
      END DO

      DO k = pox%hi(3)+nghosts_mid,pox%lo(3)-nghosts_mid,-1
        DO j = pox%hi(2)+nghosts_mid,pox%lo(2)-nghosts_mid,-1
          DO i = pox%hi(1)+nghosts_mid,pox%lo(1)-nghosts_mid,-1
            skips = .FALSE.
            skip_counter = 0
            solid_counter = 0
            fluid_counter = 0
            weird_counter = 0
            null_counter = 0

            DO a = 2,7
!
! Test the feasability of using all the possible directions
! Test the +x direction
!
              IF (i + cx_27(a) < state_bounds_lo(1) .or. i+cx_27(a) < box_lbound(1)) THEN
                skips(5) = .TRUE.
                skip_counter = skip_counter + 1
              ELSE IF (i + cx_27(a) > state_bounds_hi(1) .or. i+cx_27(a) > box_ubound(1)) THEN
                skips(2) = .TRUE.
                skip_counter = skip_counter + 1
              ELSE IF (j + cy_27(a) < state_bounds_lo(2) .or. j+cy_27(a) < box_lbound(2)) THEN
                skips(6) = .TRUE.
                skip_counter = skip_counter + 1
              ELSE IF (j + cy_27(a) > state_bounds_hi(2) .or. j+cy_27(a) > box_ubound(2)) THEN
                skips(3) = .TRUE.
                skip_counter = skip_counter + 1
              ELSE IF (k + cz_27(a) < state_bounds_lo(3) .or. k+cz_27(a) < box_lbound(3)) THEN
                skips(7) = .TRUE.
                skip_counter = skip_counter + 1
              ELSE IF (k + cz_27(a) > state_bounds_hi(3) .or. k+cz_27(a) > box_ubound(3)) THEN
                skips(4) = .TRUE.
                skip_counter = skip_counter + 1
              END IF
            END DO
!
! If the direction is usable, test the possibilities
!
        IF (.NOT. skips(2)) THEN
          IF (state(i+cx_27(2),j,k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i+cx_27(2),j,k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i+cx_27(2),j,k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF
        END IF

!
! Test the -x direction
!
        IF (.NOT. skips(5)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i+cx_27(5),j,k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i+cx_27(5),j,k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i+cx_27(5),j,k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF

        END IF
!
! Test the +y direction
!
        IF (.NOT. skips(3)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j+cy_27(3),k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j+cy_27(3),k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j+cy_27(3),k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF

        END IF
!
! Test the -y direction
!
        IF (.NOT. skips(6)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j+cy_27(6),k,1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j+cy_27(6),k,1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j+cy_27(6),k,1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF

        END IF
!
! Test the +z direction
!
        IF (.NOT. skips(4)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j,k+cz_27(4),1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j,k+cz_27(4),1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j,k+cz_27(4),1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF
        END IF
!
! Test the -z direction
!
        IF (.NOT. skips(7)) THEN
!
! If the direction is usable, test the possibilities
!
          IF (state(i,j,k+cz_27(7),1) >= 0) THEN
            fluid_counter = fluid_counter + 1
          ELSE IF (state(i,j,k+cz_27(7),1) == -1001) THEN
            solid_counter = solid_counter +1
          ELSE IF (state(i,j,k+cz_27(7),1) == -666) THEN
            null_counter = null_counter +1
          ELSE
            weird_counter = weird_counter + 1
          END IF
        END IF


            IF (state(i,j,k,1) == -1001 .AND. fluid_counter >= 5) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,101)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter >= 5) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,100)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter == 4 .AND.&
              skip_counter == 1) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,102)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter == 3 .AND.&
              skip_counter == 2) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,102)
            ELSE IF (state(i,j,k,1) >= 0 .AND. solid_counter == 3 .AND.&
              skip_counter == 3) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,102)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 4 .AND.&
              skip_counter == 1) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,103)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 3 .AND.&
              skip_counter == 2) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,103)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 3 .AND.&
              skip_counter == 3) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,103)
            ELSE IF (state(i,j,k,1) == -1001 .AND. fluid_counter == 6) THEN
              state(i,j,k,1) = lvl
!              WRITE(11,107)
            ELSE IF (state(i,j,k,1) == -666 .AND. solid_counter == 6) THEN
              state(i,j,k,1) = -1001
!              WRITE(11,106)
            ELSE IF (state(i,j,k,1) == -1001 .AND. null_counter >= 5) THEN
              state(i,j,k,1) = -666
!              WRITE(11,104)
            ELSE IF (state(i,j,k,1) >= 0 .AND. null_counter >= 5) THEN
              state(i,j,k,1) = -666
!              WRITE(11,105)
            ELSE IF (state(i,j,k,1) == -1001 .AND. null_counter >= 4 .AND.&
              skip_counter == 1) THEN
              state(i,j,k,1) = -666
!              WRITE(11,104)
            ELSE IF (state(i,j,k,1) >= 0 .AND. null_counter >= 4 .AND. &
              skip_counter == 1) THEN
              state(i,j,k,1) = -666
!              WRITE(11,105)
            ELSE IF (state(i,j,k,1) == -1001 .AND. null_counter >= 3 .AND. &
              skip_counter == 2) THEN
              state(i,j,k,1) = -666
!              WRITE(11,104)
            ELSE IF (state(i,j,k,1) >= 0 .AND. null_counter >= 3 .AND. &
              skip_counter == 2) THEN
              state(i,j,k,1) = -666


            END IF
          END DO
        END DO
      END DO

  end do

!  call amrex_multifab_destroy(temp_state)
  call amrex_mfiter_destroy(mfi)
!end do


!
! Go through all of the nodes to check if any are wildly out of place
!


 100  FORMAT ('Isolated node converted from fluid to solid.')
 101  FORMAT ('Isolated node converted from solid to fluid.')
 102  FORMAT ('Wall node converted from fluid to solid.')
 103  FORMAT ('Wall node converted from solid to fluid.')
 104  FORMAT ('Isolated node converted from solid to null.')
 105  FORMAT ('Isolated node converted from fluid to null.')
 106  FORMAT ('Isolated null node converted to solid.')
 107  FORMAT ('Isolated null node converted to fluid.')

END SUBROUTINE
