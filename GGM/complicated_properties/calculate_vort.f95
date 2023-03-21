SUBROUTINE calculate_vort(lvl,tux_lo,tux_hi,g,choice,mfi)
!
! Calculate the vorticity
!
! Calls:
! Called by: populate_q
!
USE ggm_stuff
use derivatives
USE precise
USE grid_data
USE constants
use amr_info_holder
use amrex_amr_module
use amrex_base_module
IMPLICIT NONE
INTEGER :: a,b,c
integer,intent(in) :: tux_lo(3),tux_hi(3),g,choice,lvl
integer :: substate(5),sub_u(5),sub_v(5),sub_w(5),reach
REAL(KIND=dp) :: dudy,dvdx,dwdx,dudz,dvdz,dwdy,shear_mag
!real(kind=dp) :: holder(tux_lo(1)-nghosts:tux_hi(1)+nghosts,&
!  tux_lo(2)-nghosts:tux_hi(2)+nghosts,tux_lo(3)-nghosts:tux_hi(3)+nghosts,&
!  num_ggm_vars)

type(amrex_mfiter) :: mfi

reach = nghosts_mid-2

!      DO b = 1,link_number
if (dimensions == 2) then

else
  select case (choice)
    case(1) !xy shear
      u_vel => mfu_vel(lvl)%dataptr(mfi)
      v_vel => mfv_vel(lvl)%dataptr(mfi)
      state => mfstate(lvl)%dataptr(mfi)

      do c = tux_lo(3)-reach,tux_hi(3)+reach
        do b = tux_lo(2)-reach,tux_hi(2)+reach
          do a = tux_lo(1)-reach,tux_hi(1)+reach

            dudy = derivative_windage_secondord((/u_vel(a,b,c,1),u_vel(a,b+1,c,1),u_vel(a,b+2,c,1),&
                     u_vel(a,b-1,c,1),u_vel(a,b-2,c,1)/),(/state(a,b,c,1),state(a,b,c,3),&
                     state(a,b,c,17),state(a,b,c,6),state(a,b,c,20)/),v_vel(a,b,c,1),upwind,downwind,central)
            dvdx = derivative_windage_secondord((/v_vel(a,b,c,1),v_vel(a+1,b,c,1),v_vel(a+2,b,c,1),&
                     v_vel(a-1,b,c,1),v_vel(a-2,b,c,1)/),(/state(a,b,c,1),state(a,b,c,2),&
                     state(a,b,c,16),state(a,b,c,5),state(a,b,c,19)/),u_vel(a,b,c,1),upwind,downwind,central)

            holder(a,b,c,g) = dudy-dvdx

          end do
        end do
      end do
    case(2) !xz shear
      u_vel => mfu_vel(lvl)%dataptr(mfi)
      w_vel => mfw_vel(lvl)%dataptr(mfi)
      state => mfstate(lvl)%dataptr(mfi)

      do c = tux_lo(3)-reach,tux_hi(3)+reach
        do b = tux_lo(2)-reach,tux_hi(2)+reach
          do a = tux_lo(1)-reach,tux_hi(1)+reach

            dudz = derivative_windage_secondord((/u_vel(a,b,c,1),u_vel(a,b,c+1,1),u_vel(a,b,c+2,1),&
                     u_vel(a,b,c-1,1),u_vel(a,b,c-2,1)/),(/state(a,b,c,1),state(a,b,c,4),&
                     state(a,b,c,18),state(a,b,c,7),state(a,b,c,21)/),w_vel(a,b,c,1),upwind,downwind,central)
            dwdx = derivative_windage_secondord((/w_vel(a,b,c,1),w_vel(a+1,b,c,1),w_vel(a+2,b,c,1),&
                     w_vel(a-1,b,c,1),w_vel(a-2,b,c,1)/),(/state(a,b,c,1),state(a,b,c,2),&
                     state(a,b,c,16),state(a,b,c,5),state(a,b,c,19)/),u_vel(a,b,c,1),upwind,downwind,central)

            holder(a,b,c,g) = dudz-dwdx
          end do
        end do
      end do
    case(3) !yz shear
      v_vel => mfv_vel(lvl)%dataptr(mfi)
      w_vel => mfw_vel(lvl)%dataptr(mfi)
      state => mfstate(lvl)%dataptr(mfi)

      do c = tux_lo(3)-reach,tux_hi(3)+reach
        do b = tux_lo(2)-reach,tux_hi(2)+reach
          do a = tux_lo(1)-reach,tux_hi(1)+reach

            dvdz = derivative_windage_secondord((/v_vel(a,b,c,1),v_vel(a,b,c+1,1),v_vel(a,b,c+2,1),&
                     v_vel(a,b,c-1,1),v_vel(a,b,c-2,1)/),(/state(a,b,c,1),state(a,b,c,4),&
                     state(a,b,c,18),state(a,b,c,7),state(a,b,c,21)/),w_vel(a,b,c,1),upwind,downwind,central)
            dwdy = derivative_windage_secondord((/w_vel(a,b,c,1),w_vel(a,b+1,c,1),w_vel(a,b+2,c,1),&
                     w_vel(a,b-1,c,1),w_vel(a,b-2,c,1)/),(/state(a,b,c,1),state(a,b,c,3),&
                     state(a,b,c,17),state(a,b,c,6),state(a,b,c,20)/),v_vel(a,b,c,1),upwind,downwind,central)

            holder(a,b,c,g) = dvdz-dwdy
          end do
        end do
      end do
    case(4) !shear mag
      u_vel => mfu_vel(lvl)%dataptr(mfi)
      v_vel => mfv_vel(lvl)%dataptr(mfi)
      w_vel => mfw_vel(lvl)%dataptr(mfi)
      state => mfstate(lvl)%dataptr(mfi)

      do c = tux_lo(3)-reach,tux_hi(3)+reach
        do b = tux_lo(2)-reach,tux_hi(2)+reach
          do a = tux_lo(1)-reach,tux_hi(1)+reach

            dudy = derivative_windage_secondord((/u_vel(a,b,c,1),u_vel(a,b+1,c,1),u_vel(a,b+2,c,1),&
                     u_vel(a,b-1,c,1),u_vel(a,b-2,c,1)/),(/state(a,b,c,1),state(a,b,c,3),&
                     state(a,b,c,17),state(a,b,c,6),state(a,b,c,20)/),v_vel(a,b,c,1),upwind,downwind,central)
            dudz = derivative_windage_secondord((/u_vel(a,b,c,1),u_vel(a,b,c+1,1),u_vel(a,b,c+2,1),&
                     u_vel(a,b,c-1,1),u_vel(a,b,c-2,1)/),(/state(a,b,c,1),state(a,b,c,4),&
                     state(a,b,c,18),state(a,b,c,7),state(a,b,c,21)/),w_vel(a,b,c,1),upwind,downwind,central)
            dvdx = derivative_windage_secondord((/v_vel(a,b,c,1),v_vel(a+1,b,c,1),v_vel(a+2,b,c,1),&
                     v_vel(a-1,b,c,1),v_vel(a-2,b,c,1)/),(/state(a,b,c,1),state(a,b,c,2),&
                     state(a,b,c,16),state(a,b,c,5),state(a,b,c,19)/),u_vel(a,b,c,1),upwind,downwind,central)
            dvdz = derivative_windage_secondord((/v_vel(a,b,c,1),v_vel(a,b,c+1,1),v_vel(a,b,c+2,1),&
                     v_vel(a,b,c-1,1),v_vel(a,b,c-2,1)/),(/state(a,b,c,1),state(a,b,c,4),&
                     state(a,b,c,18),state(a,b,c,7),state(a,b,c,21)/),w_vel(a,b,c,1),upwind,downwind,central)
            dwdx = derivative_windage_secondord((/w_vel(a,b,c,1),w_vel(a+1,b,c,1),w_vel(a+2,b,c,1),&
                     w_vel(a-1,b,c,1),w_vel(a-2,b,c,1)/),(/state(a,b,c,1),state(a,b,c,2),&
                     state(a,b,c,16),state(a,b,c,5),state(a,b,c,19)/),u_vel(a,b,c,1),upwind,downwind,central)
            dwdy = derivative_windage_secondord((/w_vel(a,b,c,1),w_vel(a,b+1,c,1),w_vel(a,b+2,c,1),&
                     w_vel(a,b-1,c,1),w_vel(a,b-2,c,1)/),(/state(a,b,c,1),state(a,b,c,3),&
                     state(a,b,c,17),state(a,b,c,6),state(a,b,c,20)/),v_vel(a,b,c,1),upwind,downwind,central)

            holder(a,b,c,g) = sqrt((dudy-dvdx)**2 + (dudz-dwdx)**2 + (dvdz-dwdy)**2)
          end do
        end do
      end do

  end select

end if


!      END DO

END SUBROUTINE
!      IF (dimensions == 2) THEN
!        IF (choice == 1) THEN
!          IF (bounding_method == -1 .OR. bounding_method == -2) THEN
!            q = 0.0D0
!          ELSE
!            DO a = 1,link_number
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = w_vel(link(6,a))
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            in_values(1) = w_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            in_values(5) = w_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE
!            vort_dwdy = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!            in_values(2) = v_vel(link(7,a))
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            in_values(1) = v_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            in_values(5) = v_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE
!            vort_dvdz = 0.0D0
!          END IF
!
!          q(a) = vort_dwdy-vort_dvdz
!
!            END DO
!          END IF
!
!        ELSE IF (choice == 2) THEN
!          IF (bounding_method == -1 .OR. bounding_method == -3) THEN
!            q = 0.0D0
!          ELSE
!          DO a = 1,link_number
!
!          IF (deriv_form(3,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = u_vel(link(7,a))
!            in_values(4) = u_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(7,a))
!            in_values(1) = u_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(4,a))
!            in_values(5) = u_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE
!            vort_dudz = 0.0D0
!          END IF
!
!          IF (deriv_form(1,1,a)) THEN
!            in_values(2) = w_vel(link(5,a))
!            in_values(4) = w_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(5,a))
!            in_values(1) = w_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(2,a))
!            in_values(5) = w_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE
!            vort_dwdx = 0.0D0
!          END IF
!
!          q(a) = vort_dudz-vort_dwdx
!
!          END DO
!          END IF
!        ELSE IF (choice == 3) THEN
!          IF (bounding_method == -2 .OR. bounding_method == -3) THEN
!            q = 0.0D0
!          ELSE
!          DO a = 1,link_number
!          IF (deriv_form(1,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(1) = v_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            in_values(5) = v_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE
!            vort_dvdx = 0.0D0
!          END IF
!
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = u_vel(link(6,a))
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            in_values(1) = u_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            in_values(5) = u_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE
!            vort_dudy = 0.0D0
!          END IF
!
!          q(a) = vort_dvdx-vort_dudy
!
!          END DO
!
!          END IF
!        ELSE IF (choice == 4) THEN
!!
!! vort mag for x-y plane is just ABS(w-vort)
!!
!          IF (bounding_method == -1) THEN
!
!          DO a = 1,link_number
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = w_vel(link(6,a))
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            in_values(1) = w_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            in_values(5) = w_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE
!            vort_dwdy = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!            in_values(2) = v_vel(link(7,a))
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            in_values(1) = v_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            in_values(5) = v_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE
!            vort_dvdz = 0.0D0
!          END IF
!
!          q(a) = ABS(vort_dwdy-vort_dvdz)
!          END DO
!!
!! vort mag for x-z plane is just ABS(y-vort)
!!
!          ELSE IF (bounding_method == -2) THEN
!
!          DO a = 1,link_number
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = w_vel(link(6,a))
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            in_values(1) = w_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            in_values(5) = w_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE
!            vort_dwdy = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!            in_values(2) = v_vel(link(7,a))
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            in_values(1) = v_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            in_values(5) = v_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE
!            vort_dvdz = 0.0D0
!          END IF
!
!          q(a) = ABS(vort_dwdy-vort_dvdz)
!          END DO
!!
!! vort mag for y-z plane is just ABS(x-vort)
!!
!          ELSE IF (bounding_method == -3) THEN
!
!            DO a = 1,link_number
!
!          IF (deriv_form(1,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(1) = v_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            in_values(5) = v_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE
!            vort_dvdx = 0.0D0
!          END IF
!
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = u_vel(link(6,a))
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            in_values(1) = u_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            in_values(5) = u_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE
!            vort_dudy = 0.0D0
!          END IF
!
!          q(a) = ABS(vort_dvdx-vort_dudy)
!
!          END DO
!
!            END IF
!          END IF
!
!
!        ELSE IF (dimensions == 3) THEN
!!
!! x-vorticity
!! dw _ dv
!! dy   dz
!!
!        IF (choice == 1) THEN
!
!          DO a = 1,link_number
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = w_vel(link(6,a))
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            in_values(1) = w_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            in_values(5) = w_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE
!            vort_dwdy = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!            in_values(2) = v_vel(link(7,a))
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            in_values(1) = v_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            in_values(5) = v_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!
!          ELSE
!            vort_dvdz = 0.0D0
!          END IF
!
!          q(a) = vort_dwdy-vort_dvdz
!
!          END DO
!!
!! y-vorticity
!! du _ dw
!! dz   dx
!!
!        ELSE IF (choice == 2) THEN
!
!          DO a = 1,link_number
!
!          IF (deriv_form(3,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = u_vel(link(7,a))
!            in_values(4) = u_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(7,a))
!            in_values(1) = u_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(4,a))
!            in_values(5) = u_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE
!            vort_dudz = 0.0D0
!          END IF
!
!          IF (deriv_form(1,1,a)) THEN
!            in_values(2) = w_vel(link(5,a))
!            in_values(4) = w_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(5,a))
!            in_values(1) = w_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(2,a))
!            in_values(5) = w_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE
!            vort_dwdx = 0.0D0
!          END IF
!
!          q(a) = vort_dudz-vort_dwdx
!
!          END DO
!!
!! z-vorticity
!! dv _ du
!! dx   dy
!!
!        ELSE IF (choice == 3) THEN
!
!          DO a = 1,link_number
!
!          IF (deriv_form(1,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(1) = v_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            in_values(5) = v_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!
!          ELSE
!            vort_dvdx = 0.0D0
!          END IF
!
!
!          IF (deriv_form(2,1,a)) THEN
!            in_values(2) = u_vel(link(6,a))
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            in_values(1) = u_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            in_values(5) = u_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!
!          ELSE
!            vort_dudy = 0.0D0
!          END IF
!
!          q(a) = vort_dvdx-vort_dudy
!
!          END DO
!
!!
!!
!! vort_mag = SQRT(x-vort^2 + y-vort^2 + z-vort^2)
!!
!!
!        ELSE IF (choice == 4) THEN
!
!          DO a = 1,link_number
!
!          IF (deriv_form(1,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!            in_values(2) = w_vel(link(5,a))
!            in_values(4) = w_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            in_values(1) = v_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(5,a))
!            in_values(1) = w_vel(link(5,link(5,a)))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            in_values(5) = v_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(2,a))
!            in_values(5) = w_vel(link(2,link(2,a)))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(5,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE IF (deriv_form(1,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dvdx,in_values)
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(2,a))
!            CALL x_dir_derivative(a,vort_dwdx,in_values)
!
!          ELSE
!            vort_dvdx = 0.0D0
!            vort_dwdx = 0.0D0
!          END IF
!
!          IF (deriv_form(2,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = u_vel(link(6,a))
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!            in_values(2) = w_vel(link(6,a))
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,2,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            in_values(1) = u_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            in_values(1) = w_vel(link(6,link(6,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,3,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            in_values(5) = u_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            in_values(5) = w_vel(link(3,link(3,a)))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,4,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!            in_values(3) = w_vel(a)
!            in_values(2) = w_vel(link(6,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE IF (deriv_form(2,5,a)) THEN
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dudy,in_values)
!            in_values(3) = w_vel(a)
!            in_values(4) = w_vel(link(3,a))
!            CALL y_dir_derivative(a,vort_dwdy,in_values)
!
!          ELSE
!            vort_dudy = 0.0D0
!            vort_dwdy = 0.0D0
!          END IF
!
!          IF (deriv_form(3,1,a)) THEN
!!            in_values(3) = q(a)
!            in_values(2) = v_vel(link(7,a))
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!            in_values(2) = u_vel(link(7,a))
!            in_values(4) = u_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,2,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            in_values(1) = v_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(7,a))
!            in_values(1) = u_vel(link(7,link(7,a)))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,3,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            in_values(5) = v_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(4,a))
!            in_values(5) = u_vel(link(4,link(4,a)))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,4,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(2) = v_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!            in_values(3) = u_vel(a)
!            in_values(2) = u_vel(link(7,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE IF (deriv_form(3,5,a)) THEN
!            in_values(3) = v_vel(a)
!            in_values(4) = v_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dvdz,in_values)
!            in_values(3) = u_vel(a)
!            in_values(4) = u_vel(link(4,a))
!            CALL z_dir_derivative(a,vort_dudz,in_values)
!
!          ELSE
!            vort_dudz = 0.0D0
!            vort_dvdz = 0.0D0
!          END IF
!
!
!            q(a) = SQRT((vort_dwdy-vort_dvdz)**2.0D0 + (vort_dudz-vort_dwdx)**2.0D0 +&
!              (vort_dvdx - vort_dudy)**2.0D0)
!
!          END DO
!
!          END IF
!
!        END IF
