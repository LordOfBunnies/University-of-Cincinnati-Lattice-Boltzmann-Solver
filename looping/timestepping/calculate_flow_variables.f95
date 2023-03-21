SUBROUTINE calculate_flow_variables(level_step,fine,lvl)
!
! This calculates the flow variables based on values of fi and gi
!
! Called by: loop
! Calls: interlevel_interpolater
! External calls: amrex_mfiter_build,amrex_mfiter_destroy
!
USE precise
use amrex_amr_module
use amrex_base_module
use quick_calcs
USE mpi
USE startup
USE constants
use linkwise
USE freestream_values
USE amr_info_holder
use amr_processes, only: amr_max_lvl,self,commune
use grid_data, only: dir
use mpi
IMPLICIT NONE
REAL(KIND=dp) :: rho_sum,u_sum,v_sum,w_sum,g_sum
real(kind=dp) :: old_temp,old_u,old_v,old_w,old_rho
real(kind=dp) :: new_pxx,new_pyy,new_pzz,new_pxy,new_pxz,new_pyz
real(kind=dp) :: loc_fieq(dir),loc_fi_diff(dir)
real(kind=dp) :: dummy_lambda_z,dummy_w,dummy_pi_xz,dummy_pi_yz,&
  dummy_pi_zz,dummy_zeta_z,energy_sum,v_mag
real(kind=dp) :: tmp_chi,tmp_zeta_x,tmp_zeta_y,tmp_zeta_z,tmp_pi_xx,tmp_pi_yy,&
  tmp_pi_zz,tmp_pi_xy,tmp_pi_xz,tmp_pi_yz,tmp_lambda_x,tmp_lambda_y,tmp_lambda_z
INTEGER :: i,a,b,c,cur_lvl_boxes,ier,reach
integer,intent(in) :: lvl,fine
logical :: level_step(0:amr_max_lvl)

TYPE(amrex_mfiter) :: mfi
type(amrex_box) :: lox
!
!
!
!
!
if (self == 0) then
  write(*,*) 'calculating flow variables',lvl,level_step
end if
!
if (.not. shifted) then
if (lvl == 0) then
  reach = nghosts/2
else
  if (mod(timestep(lvl),2) == 0) then
    reach = nghosts/2
  else
    reach = 0
  end if
end if
else
if (lvl == 0) then
  reach = nghosts/2
else
  if (mod(timestep(lvl),2) == 0) then
    reach = nghosts/2
  else
    reach = 0
  end if
end if
end if
!
! Must step through all valid levels, this is found elsewhere
!
call mffi(lvl)% fill_boundary( amrex_geom(lvl),1,dir)
call mfgi(lvl)% fill_boundary( amrex_geom(lvl),1,dir)
!  write(*,*) 'Flow variable calculation on level ',lvl

  IF (level_step(lvl)) THEN
    cur_lvl_boxes = 0
!
! Build the memory piece to be able to do data pointers to the flow data
!
    CALL amrex_mfiter_build(mfi,mffi(lvl),tiling=.false.)
!
! Cycle through all the blocks on a given level
!
    do while (mfi%next())
      cur_lvl_boxes = cur_lvl_boxes+1
      lox = mfi%validbox()
!
!
!
      rho => mfrho(lvl)%dataptr(mfi)
      u_vel => mfu_vel(lvl)%dataptr(mfi)
      v_vel => mfv_vel(lvl)%dataptr(mfi)
      temp => mftemp(lvl)%dataptr(mfi)
      press => mfpress(lvl)%dataptr(mfi)
      state => mfstate(lvl)%dataptr(mfi)
!
!
!
      fi => mffi(lvl)%dataptr(mfi)
      gi => mfgi(lvl)%dataptr(mfi)
      alpha => mfalpha(lvl)%dataptr(mfi)
!
! LMs
!
      chi => mfchi(lvl)%dataptr(mfi)
      zeta_x => mfzeta_x(lvl)%dataptr(mfi)
      zeta_y => mfzeta_y(lvl)%dataptr(mfi)
      pi_xx => mfpi_xx(lvl)%dataptr(mfi)
      pi_yy => mfpi_yy(lvl)%dataptr(mfi)
      pi_xy => mfpi_xy(lvl)%dataptr(mfi)
      lambda_x => mflambda_x(lvl)%dataptr(mfi)
      lambda_y => mflambda_y(lvl)%dataptr(mfi)
!
! 3D values
!
      IF (dimensions == 3) THEN
        zeta_z => mfzeta_z(lvl)%dataptr(mfi)
        w_vel => mfw_vel(lvl)%dataptr(mfi)
        pi_xz => mfpi_xz(lvl)%dataptr(mfi)
        pi_zz => mfpi_zz(lvl)%dataptr(mfi)
        pi_yz => mfpi_yz(lvl)%dataptr(mfi)
        lambda_z => mflambda_z(lvl)%dataptr(mfi)

      else
        dummy_zeta_z = 0.0D0
        dummy_w = 0.0D0
        dummy_pi_xz = 0.0D0
        dummy_pi_yz = 0.0D0
        dummy_pi_zz = 0.0D0
        dummy_lambda_z = 0.0D0
      END IF
!
! Cycle through all nodes in the block
!
      do c = lox%lo(3)-reach,lox%hi(3)+reach
        do b = lox%lo(2)-reach,lox%hi(2)+reach
          do a = lox%lo(1)-reach,lox%hi(1)+reach
!
! Skip null nodes
!
            if (state(a,b,c,1) < 0) cycle
!
! Skip overlapped ghost nodes which get their values using fill_boundary
!
! May be able to skip 51k-100k because the information will come from the fine grid
!
            if (state(a,b,c,1) >= 1000 .and. state(a,b,c,1) < 51000) cycle
!
! conditions for 3-D terms
!
            if (dimensions == 3) then
              dummy_zeta_z = zeta_z(a,b,c,1)
              dummy_w = w_vel(a,b,c,1)
              dummy_pi_xz = pi_xz(a,b,c,1)
              dummy_pi_zz = pi_zz(a,b,c,1)
              dummy_pi_yz = pi_yz(a,b,c,1)
              dummy_lambda_z = lambda_z(a,b,c,1)
            end if
!
!            if (a == 257 .and. b == 1 .and. c == 0) then
!              write(*,*) 'pre-calcs flow vals',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!                w_vel(a,b,c,1),temp(a,b,c,1)
!              write(*,*) 'beep beep',state(a,b,c,1:dir),a,b,c
!              write(*,*) 'JINGOISM JINGOISM',fi(a,b,c,1:dir)
!              write(*,*) 'jumpin jehosephat',gi(a,b,c,1:dir)
!              write(*,*) 'pre-LMs',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),dummy_pi_zz,pi_xy(a,b,c,1),dummy_pi_xz,&
!              dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z
!            end if
!
!            if (a == 64 .and. b == 1 .and. c == 1) then
!              write(*,*) 'pre-calcs flow vals',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!                w_vel(a,b,c,1),temp(a,b,c,1)
!              write(*,*) 'beep beep',state(a,b,c,1:dir),a,b,c
!              write(*,*) 'JINGOISM JINGOISM',fi(a,b,c,1:dir)
!              write(*,*) 'jumpin jehosephat',gi(a,b,c,1:dir)
!              write(*,*) 'pre-LMs',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),dummy_pi_zz,pi_xy(a,b,c,1),dummy_pi_xz,&
!              dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z
!            end if
!
!            if (a == 479 .and. b == 239 .and. c == 72) then
!              write(*,*) 'state of the error',state(a,b,c,1:dir),lox%lo,lox%hi,a,b,c
!              write(*,*) 'flow variables',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!                dummy_w,temp(a,b,c,1)
!              write(*,*) 'lagrange mults',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
!                pi_xx(a,b,c,1),pi_yy(a,b,c,1),dummy_pi_zz,pi_xy(a,b,c,1),dummy_pi_xz,&
!                dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z
!            end if
!
! Initialize the sums to 0
!
            rho_sum = 0.0D0
            u_sum = 0.0D0
            v_sum = 0.0D0
            w_sum = 0.0D0
            energy_sum = 0.0D0
            g_sum = 0.0D0
!
! sum up fi and gi to find density, velocities, and temperature
!
            do i = 1,dir
              rho_sum = rho_sum + fi(a,b,c,i)
              u_sum = u_sum + rcx(i)*fi(a,b,c,i)
              v_sum = v_sum + rcy(i)*fi(a,b,c,i)
              energy_sum = energy_sum + cmag(i)*fi(a,b,c,i)+gi(a,b,c,i)
              g_sum = g_sum + gi(a,b,c,i)
            end do
!
            rho(a,b,c,1) = rho_sum
            u_vel(a,b,c,1) = u_sum/rho_sum
            v_vel(a,b,c,1) = v_sum/rho_sum
            if (dimensions ==3 ) then
              do i = 1,dir
                w_sum = w_sum + rcz(i)*fi(a,b,c,i)
              end do
              dummy_w = w_sum/rho_sum
            end if
!
            v_mag = u_vel(a,b,c,1)**2 + v_vel(a,b,c,1)**2 + dummy_w**2
            temp(a,b,c,1) = (energy_sum - rho(a,b,c,1)*v_mag)/(2*const_vol*rho(a,b,c,1))
            press(a,b,c,1) = rho(a,b,c,1)*(temp(a,b,c,1) + v_mag/2.0D0)

            new_pxx = rho(a,b,c,1)*(temp(a,b,c,1) + u_vel(a,b,c,1)**2)
            new_pyy = rho(a,b,c,1)*(temp(a,b,c,1) + v_vel(a,b,c,1)**2)
            new_pzz = rho(a,b,c,1)*(temp(a,b,c,1) + w_vel(a,b,c,1)**2)
            new_pxy = rho(a,b,c,1)*u_vel(a,b,c,1)*v_vel(a,b,c,1)
            new_pxz = rho(a,b,c,1)*u_vel(a,b,c,1)*dummy_w
            new_pyz = rho(a,b,c,1)*v_vel(a,b,c,1)*dummy_w

            !write(*,*) 'current location',a,b,c,state(a,b,c,1),lvl,self
            !write(*,*) rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
!            if (a == 70 .and. b == 32 .and. c == 32) then
!              write(*,*) 'sending values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!              dummy_w,temp(a,b,c,1),new_pxx,new_pyy,new_pzz,new_pxy,new_pxz,new_pyz
!            end if
            if (timestep(lvl) <= 1 .and. any(state(a,b,c,1:dir)<-1000) .and. shifted) then

              tmp_chi = su_chi
              tmp_zeta_x = su_zeta_x
              tmp_zeta_y = su_zeta_y
              tmp_zeta_z = su_zeta_z
              tmp_pi_xx = su_pixx
              tmp_pi_yy = su_piyy
              tmp_pi_zz = su_pizz
              tmp_pi_xy = su_pixy
              tmp_pi_xz = su_pixz
              tmp_pi_yz = su_piyz
              tmp_lambda_x = su_lambda_x
              tmp_lambda_y = su_lambda_y
              tmp_lambda_z = su_lambda_z

              CALL calculate_lagrangian_mults(rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
                dummy_w,temp(a,b,c,1),new_pxx,new_pyy,new_pzz,new_pxy,new_pxz,new_pyz,&
                tmp_chi,tmp_zeta_x,tmp_zeta_y,tmp_zeta_z,&
                tmp_pi_xx,tmp_pi_yy,tmp_pi_zz,tmp_pi_xy,tmp_pi_xz,&
                tmp_pi_yz,tmp_lambda_x,tmp_lambda_y,tmp_lambda_z,loc_fieq,a,b,c)

              chi(a,b,c,1) = tmp_chi
              zeta_x(a,b,c,1) = tmp_zeta_x
              zeta_y(a,b,c,1) = tmp_zeta_y
!              zeta_z(a,b,c,1) = tmp_zeta_z
              dummy_zeta_z = tmp_zeta_z
              pi_xx(a,b,c,1) = tmp_pi_xx
              pi_yy(a,b,c,1) = tmp_pi_yy
!              pi_zz(a,b,c,1) = tmp_pizz
              dummy_pi_zz = tmp_pi_zz
              pi_xy(a,b,c,1) = tmp_pi_xy
!              pi_xz(a,b,c,1) = tmp_pixz
!              pi_yz(a,b,c,1) = tmp_piyz
              dummy_pi_xz = tmp_pi_xz
              dummy_pi_yz = tmp_pi_yz
              lambda_x(a,b,c,1) = tmp_lambda_x
              lambda_y(a,b,c,1) = tmp_lambda_y
!              lambda_z(a,b,c,1) = tmp_lambda_z
              dummy_lambda_z = tmp_lambda_z
            else

              CALL calculate_lagrangian_mults(rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
                dummy_w,temp(a,b,c,1),new_pxx,new_pyy,new_pzz,new_pxy,new_pxz,new_pyz,&
                chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
                pi_xx(a,b,c,1),pi_yy(a,b,c,1),dummy_pi_zz,pi_xy(a,b,c,1),dummy_pi_xz,&
                dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z,loc_fieq,a,b,c)

            end if
!
!            if (any(isnan(loc_fieq))) then
!              write(*,*) 'local fieq',loc_fieq,a,b,c
!              write(*,*) 'local fis',fi(a,b,c,1:dir),a,b,c
!              write(*,*) 'local stuff',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
!            end if
!
            loc_fi_diff = loc_fieq(1:dir)-fi(a,b,c,1:dir)
!
!            if (a == 368 .and. b == 264 .and. c == 255 .and. lvl == 3) then
!              write(*,*) 'how the fluff',fi(a,b,c,1:dir),a,b,c,self
!              write(*,*) 'fieqs after',loc_fieq
!              write(*,*) 'flow variables afterward',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),&
!                dummy_w,temp(a,b,c,1)
!              write(*,*) 'lagrange mults afterward',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
!                pi_xx(a,b,c,1),pi_yy(a,b,c,1),dummy_pi_zz,pi_xy(a,b,c,1),dummy_pi_xz,&
!                dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z
!
!!              write(*,*) 'oh my fluffs',state(a,b,c,1:dir),lvl,a,b,c,self
!            end if

            if (turbulent) then
              CALL calculate_alpha(fi(a,b,c,1:dir),loc_fi_diff(1:dir),alpha(a,b,c,1),temp(a,b,c,1),&
                state(a,b,c,1),lvl,a,b,c)
            else
              alpha(a,b,c,1) = 2.0D0
            end if

!           if (a == 128 .and. b == 1 .and. c == 54) then
!              write(*,*) 'post calc vals',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!                temp(a,b,c,1)
!              write(*,*) 'post calc-LMs',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),dummy_pi_zz,pi_xy(a,b,c,1),dummy_pi_xz,&
!              dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z
!            end if

            IF (dimensions == 3) THEN
              w_vel(a,b,c,1) = dummy_w
              pi_zz(a,b,c,1) = dummy_pi_zz
              pi_xz(a,b,c,1) = dummy_pi_xz
              pi_yz(a,b,c,1) = dummy_pi_yz
              zeta_z(a,b,c,1) = dummy_zeta_z
              lambda_z(a,b,c,1) = dummy_lambda_z
            END IF

          END DO
        END DO
      END DO
    end do

!    n_boxes_lvl(lvl) = cur_lvl_boxes
!
! Destroy the mfiter so we can move on.
!
    call amrex_mfiter_destroy(mfi)
    !call amrex_mfiter_destroy(mfj)

!  ELSE
!    CYCLE
  END IF


END SUBROUTINE
!            if (a == 97 .and. b == 53 .and. c == 49) then
!              write(*,*) 'After calculating',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!                temp(a,b,c,1),alpha(a,b,c,1)
!            end if

!END DO
!call mpi_barrier(commune,ier)
!call fine_to_coarse(lvl,fine)

!            if (lvl ==3 .and. timestep(lvl) == 1) then
!              write(*,*) 'flow variables',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!                temp(a,b,c,1),fi(a,b,c,1:dir)
!            end if
            !write(*,*) 'flow vars',rho_sum,u_sum,v_sum,w_sum,energy_sum,v_mag,g_sum,a,b,c,lvl,self
            !do j = 1,dir
            !  v_mag = SQRT(u_sum**2.0 + v_sum**2.0 + w_sum**2.0)
            !end do

!            if (rho(a,b,c,1) < 0.6D0 .or. rho(a,b,c,1) > 1.5D0) then
!              write(*,*) 'large changes in flow variables',rho(a,b,c,1),u_vel(a,b,c,1),&
!                v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),state(a,b,c,1),a,b,c,lvl,timestep(lvl)
!              write(*,*) 'state values',state(a,b,c,1:dir),a,b,c,lvl
!              write(*,*) 'fi values',fi(a,b,c,1:dir),a,b,c,lvl
!              !write(*,*) 'fout values',fout(a,b,c,1:dir)
!              write(*,*) 'gi values',gi(a,b,c,1:dir),a,b,c,lvl
!              !write(*,*) 'gout values',gout(a,b,c,1:dir)
!              write(*,*) 'LMs',chi(a,b,c,1),zeta_x(a,b,c,1),zeta_y(a,b,c,1),dummy_zeta_z,&
!                pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_xy(a,b,c,1),dummy_pi_zz,dummy_pi_xz,&
!                dummy_pi_yz,lambda_x(a,b,c,1),lambda_y(a,b,c,1),dummy_lambda_z
!            end if

!
!
!
!            if (abs(v_vel(a,b,c,1)) >= v_inf .or. abs(dummy_w) >= v_inf) then
!              write(*,*) 'weird velocities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),dummy_w,&
!                temp(a,b,c,1),a,b,c,lvl
!              !write(*,*) 'location',a,b,c,i,lvl
!              write(*,*) 'states',state(a,b,c,1:dir),a,b,c,lvl
!              write(*,*) 'fi',fi(a,b,c,1:dir),a,b,c,lvl
!              write(*,*) 'gi',gi(a,b,c,1:dir),a,b,c,lvl
!
!            end if
!
!
!
!              write(*,*) 'local numbers ',fi(a,b,c,i),gi(a,b,c,i),a,b,c,i

              !c_mag2 = REAL(cx(direct)**2)+REAL(cy(direct)**2)+REAL(cz(direct)**2)
              !temp_sum = temp_sum + gi(a,b,c,i)

!              if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) <0.0D0) then
!                write(*,*) 'negative fi',fi(a,b,c,i),state(a,b,c,1:dir),lox%lo,lox%hi,a,b,c,&
!                  i,lvl,self,timestep(lvl)
!                write(*,*) 'negative gi',gi(a,b,c,i),state(a,b,c,1:dir),lox%lo,lox%hi,a,b,c,&
!                  i,lvl,self,timestep(lvl)
!              end if

!            if (abs(old_temp - temp(a,b,c,1)) > 0.05) then
!            if (temp(a,b,c,1) < 0.5D0 .or. temp(a,b,c,1) > 1.0D0) then
!              write(*,*) 'large temperature change',old_temp,temp(a,b,c,1),old_rho,rho(a,b,c,1),&
!                old_u,u_vel(a,b,c,1),old_v,v_vel(a,b,c,1),old_w,w_vel(a,b,c,1)
!              write(*,*) 'why is it dropping',old_temp,temp(a,b,c,1),g_sum,lvl,self
!              write(*,*) 'local gis',gi(a,b,c,1:dir)
!              write(*,*) 'local states',state(a,b,c,1:dir)
!              write(*,*) 'nearby temps',temp(a,b,c,1),temp(a+1,b,c,1),temp(a-1,b,c,1),temp(a,b+1,c,1),&
!                temp(a,b-1,c,1),temp(a,b,c+1,1),temp(a,b,c-1,1)
!            end if

            !temp(a,b,c,1) = (energy_sum-rho_sum*v_mag)/(2*rho_sum*const_vol)
!            if (a == 289 .and. b == 127 .and. c == 122 .and. lvl == 3) then
!              old_temp = temp(a,b,c,1)
!              temp(a,b,c,1) = g_sum/(2*const_vol - dimensions)
!              write(*,*) 'why is it dropping',old_temp,temp(a,b,c,1),g_sum,state(a,b,c,1),self
!              write(*,*) 'local gis',gi(a,b,c,1:dir)
!            end if
!            old_temp = temp(a,b,c,1)


!            if (a == 1 .and. b == 16 .and. c == 31) then
!              write(*,*) 'Before calculating',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!                temp(a,b,c,1)
!            end if
