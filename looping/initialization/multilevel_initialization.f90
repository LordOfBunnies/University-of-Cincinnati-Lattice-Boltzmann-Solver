SUBROUTINE multilevel_initialization
!
! Start with a basic multigrid solution around the area of interest at the level
! specified in the input file by the user.
!
! Called by: init_fi
! Calls:
!
use precise
USE amrex_base_module
use amrex_amr_module
use mpi
use amr_info_holder
use amr_processes, only: amr_max_lvl,self,commune
use constants
use linkwise
use freestream_values
use grid_data
use mpi
IMPLICIT NONE
integer :: lvl,i,a,b,c
integer :: fine,ier
!real(kind=dp) :: rest_chi,rest_zeta_x,rest_zeta_y,rest_zeta_z,rest_pixx,rest_pixy,&
!  rest_piyy,rest_pizz,rest_pixz,rest_piyz,rest_lambda_x,rest_lambda_y,rest_lambda_z
real(kind=dp) :: rest_rho
real(kind=dp) :: rest_fieq(dir),rest_gieq(dir)

type(amrex_mfiter) :: mfi
type(amrex_box) :: nox

fine = amrex_get_finest_level()

    write(*,*) 'initialization values ',self,rho_inf,fs_u,fs_v,fs_w,temperature,fs_alpha,self
!    write(*,*) 'initial LMs ',self,fs_chi,fs_zeta_x,fs_zeta_y,fs_zeta_z,fs_pixx,&
!      fs_piyy,fs_pizz,fs_pixy,fs_pixz,fs_piyz,fs_lambda_x,fs_lambda_y,fs_lambda_z
!    write(*,*) 'initial fi ',fs_fieq
!    write(*,*) 'initial gi ',fs_gieq

rest_press = stag_density*stag_temperature

!initial LMs
if (.not. shifted) then
  rest_chi=-2.4D0
  rest_zeta_x=0.0D0
  rest_zeta_y=0.0D0
  rest_zeta_z=0.0D0
  rest_pixx=-0.40D0
  rest_piyy=-0.40D0
  rest_pizz=-0.40D0
  rest_pixy=0.0D0
  rest_pixz=0.0D0
  rest_piyz=0.0D0
  rest_lambda_x=0.0D0
  rest_lambda_y=0.0D0
  rest_lambda_z=0.0D0
else
    if (primary_flow_dir == 1) then
    rest_zeta_x = -mach/4.0D0
    rest_pixx = -0.2D0*1.5

    rest_zeta_y = 0.0D0
    rest_zeta_z = 0.0D0
    rest_pizz = -0.2D0
    rest_piyy = -0.2D0
    else if (primary_flow_dir == 2) then
    rest_zeta_y = -mach/4.0D0
    rest_piyy = -0.2D0*1.5

    rest_zeta_x = 0.0D0
    rest_zeta_z = 0.0D0
    rest_pixx = -0.2D0
    rest_pizz = -0.2D0
    else if (primary_flow_dir == 3) then
    rest_zeta_z = -mach/4.0D0
    rest_pizz = -0.2D0*1.5

    rest_zeta_y = 0.0D0
    rest_zeta_x = 0.0D0
    rest_pixx = -0.2D0
    rest_piyy = -0.2D0
    else if (primary_flow_dir == 4) then
    rest_zeta_x = mach/4.0D0
    rest_pixx = -0.2D0*1.5

    rest_zeta_y = 0.0D0
    rest_zeta_z = 0.0D0
    rest_pizz = -0.2D0
    rest_piyy = -0.2D0
    else if (primary_flow_dir == 5) then
    rest_zeta_y = mach/4.0D0
    rest_piyy = -0.2D0*1.5

    rest_zeta_x = 0.0D0
    rest_zeta_z = 0.0D0
    rest_pixx = -0.2D0
    rest_pizz = -0.2D0
    else if (primary_flow_dir == 6) then
    rest_zeta_z = mach/4.0D0
    rest_pizz = -0.2D0*1.5

    rest_zeta_y = 0.0D0
    rest_zeta_x = 0.0D0
    rest_pixx = -0.2D0
    rest_piyy = -0.2D0
    end if

  rest_chi=-2.4D0
  !rest_pixx=-0.280D0
  !rest_piyy=-0.280D0
  !rest_pizz=-0.280D0
  rest_pixy=0.0D0
  rest_pixz=0.0D0
  rest_piyz=0.0D0
  rest_lambda_x = -0.005D0
  rest_lambda_y = 0.0D0
  rest_lambda_z = 0.0D0
end if

call calculate_lagrangian_mults(stag_density,0.0D0,0.0D0,0.0D0,stag_temperature,rest_press,&
  rest_press,rest_press,0.0D0,0.0D0,0.0D0,&
  rest_chi,rest_zeta_x,rest_zeta_y,rest_zeta_z,rest_pixx,&
  rest_piyy,rest_pizz,rest_pixy,rest_pixz,rest_piyz,rest_lambda_x,rest_lambda_y,rest_lambda_z,&
  rest_fieq,1,1,1)

rest_gieq = rest_fieq*(2*const_vol - dimensions)*temperature

write(*,*) 'rest fieq',rest_fieq
write(*,*) 'rest gieq',rest_gieq
write(*,*) 'rest values',stag_density,stag_temperature,rest_press,rest_chi,rest_zeta_x,rest_zeta_y,rest_zeta_z,rest_pixx,&
  rest_piyy,rest_pizz,rest_pixy,rest_pixz,rest_piyz,rest_lambda_x,rest_lambda_y,rest_lambda_z

do lvl = 0,fine

  call mpi_barrier(commune,ier)

  call amrex_mfiter_build(mfi,mffi(lvl),tiling=.false.)

  do while (mfi%next())
    nox=mfi%validbox()

    write(*,*) 'noxious',nox%lo,nox%hi,lvl,self

    state => mfstate(lvl)%dataptr(mfi)
    fi => mffi(lvl)%dataptr(mfi)
    gi => mfgi(lvl)%dataptr(mfi)
    fout => mffout(lvl)%dataptr(mfi)
    gout => mfgout(lvl)%dataptr(mfi)
    alpha => mfalpha(lvl)%dataptr(mfi)

    rho => mfrho(lvl)%dataptr(mfi)
    u_vel => mfu_vel(lvl)%dataptr(mfi)
    v_vel => mfv_vel(lvl)%dataptr(mfi)
    temp => mftemp(lvl)%dataptr(mfi)
    press => mfpress(lvl)%dataptr(mfi)
!
!
!
    chi => mfchi(lvl)%dataptr(mfi)
    zeta_x => mfzeta_x(lvl)%dataptr(mfi)
    zeta_y => mfzeta_y(lvl)%dataptr(mfi)
    pi_xx => mfpi_xx(lvl)%dataptr(mfi)
    pi_yy => mfpi_yy(lvl)%dataptr(mfi)
    pi_xy => mfpi_xy(lvl)%dataptr(mfi)
    lambda_x => mflambda_x(lvl)%dataptr(mfi)
    lambda_y => mflambda_y(lvl)%dataptr(mfi)
!    write(*,*) 'faddle'
    IF (dimensions == 3) THEN
!      write(*,*) 'blarg?'
      zeta_z => mfzeta_z(lvl)%dataptr(mfi)
!              dummy_zeta = zeta_z(a,b,c,1)
      w_vel => mfw_vel(lvl)%dataptr(mfi)
!              dummy_w = w_vel(a,b,c,1)
!
      pi_xz => mfpi_xz(lvl)%dataptr(mfi)
      pi_zz => mfpi_zz(lvl)%dataptr(mfi)
      pi_yz => mfpi_yz(lvl)%dataptr(mfi)
      lambda_z => mflambda_z(lvl)%dataptr(mfi)
!      write(*,*) 'twaddle'

    end if

!    write(*,*) 'boxes and loxes and noxes!!! ',self,lvl,nox%lo,nox%hi
!    write(*,*) 'pre-assignement values for rho and u',rho(1,1,1,1),u_vel(1,1,1,1)

    do c = nox%lo(3)-nghosts_mid,nox%hi(3)+nghosts_mid
      do b = nox%lo(2)-nghosts_mid,nox%hi(2)+nghosts_mid
        do a = nox%lo(1)-nghosts_mid,nox%hi(1)+nghosts_mid

!            if (a == 113 .and. b == 49 .and. c == 49 .and. lvl == 2) then
!              write(*,*) 'node status, at init',state(a,b,c,1),a,b,c,lvl
!            end if

          if (state(a,b,c,1) < -1000) then
            rho(a,b,c,1) = rho_inf
            u_vel(a,b,c,1) = 0.0D0
            v_vel(a,b,c,1) = 0.0D0
            temp(a,b,c,1) = temperature
            press(a,b,c,1) = pressure

            alpha(a,b,c,1) = fs_alpha

            chi(a,b,c,1) = rest_chi
            zeta_x(a,b,c,1) = rest_zeta_x
            zeta_y(a,b,c,1) = rest_zeta_y
            pi_xx(a,b,c,1) = rest_pixx
            pi_yy(a,b,c,1) = rest_piyy
            pi_xy(a,b,c,1) = rest_pixy
            lambda_x(a,b,c,1) = rest_lambda_x
            lambda_y(a,b,c,1) = rest_lambda_y
!
!
!
            if (dimensions == 3) then
              w_vel(a,b,c,1) = 0.0D0

              zeta_z(a,b,c,1) = rest_zeta_z
              pi_xz(a,b,c,1) = rest_pixz
              pi_zz(a,b,c,1) = rest_pizz
              pi_yz(a,b,c,1) = rest_piyz
              lambda_z(a,b,c,1) = rest_lambda_z
            end if
            !write(*,*) 'local rho value ',self,rho(1,1,1,1)
!            if (a == 51 .and. b == 37 .and. c == 2 .and. lvl == 1) then
!              write(*,*) 'ARGH!!! in init',state(a,b,c,1),state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!                    state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),i,nox%lo,nox%hi,self
!            end if
            !write(*,*) 'local u value ',self,u_vel(1,1,1,1)

            do i = 1,dir
              fi(a,b,c,i) = rest_fieq(i)
              gi(a,b,c,i) = rest_gieq(i)
              fout(a,b,c,i) = rest_fieq(i)
              gout(a,b,c,i) = rest_gieq(i)
            end do

          else

            rho(a,b,c,1) = rho_inf
            u_vel(a,b,c,1) = fs_u
            v_vel(a,b,c,1) = fs_v
            temp(a,b,c,1) = temperature
            press(a,b,c,1) = pressure

            alpha(a,b,c,1) = fs_alpha

            chi(a,b,c,1) = fs_chi
            zeta_x(a,b,c,1) = fs_zeta_x
            zeta_y(a,b,c,1) = fs_zeta_y
            pi_xx(a,b,c,1) = fs_pixx
            pi_yy(a,b,c,1) = fs_piyy
            pi_xy(a,b,c,1) = fs_pixy
            lambda_x(a,b,c,1) = fs_lambda_x
            lambda_y(a,b,c,1) = fs_lambda_y
!
!
!
            if (dimensions == 3) then
              w_vel(a,b,c,1) = fs_w

              zeta_z(a,b,c,1) = fs_zeta_z
              pi_xz(a,b,c,1) = fs_pixz
              pi_zz(a,b,c,1) = fs_pizz
              pi_yz(a,b,c,1) = fs_piyz
              lambda_z(a,b,c,1) = fs_lambda_z
            end if
            !write(*,*) 'local rho value ',self,rho(1,1,1,1)
!            if (a == 51 .and. b == 37 .and. c == 2 .and. lvl == 1) then
!              write(*,*) 'ARGH!!! in init',state(a,b,c,1),state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!                    state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),i,nox%lo,nox%hi,self
!            end if
            !write(*,*) 'local u value ',self,u_vel(1,1,1,1)

            do i = 1,dir
              fi(a,b,c,i) = fs_fieq(i)
              gi(a,b,c,i) = fs_gieq(i)
              fout(a,b,c,i) = fs_fieq(i)
              gout(a,b,c,i) = fs_gieq(i)
            end do

          end if

        end do
      end do
    end do

!    write(*,*) 'local rho value ',self,lvl,rho(nox%lo(1),nox%lo(2),nox%lo(3),1),rho_inf
!    write(*,*) 'local u value ',self,lvl,u_vel(nox%lo(1),nox%lo(2),nox%lo(3),1),fs_u
!    write(*,*) 'local fi values ',self,lvl,fi(nox%lo(1),nox%lo(2),nox%lo(3),1),&
!      fi(nox%lo(1),nox%lo(2),nox%lo(3),dir),&
!      fs_fieq(1),fs_fieq(dir)
  end do

  call amrex_mfiter_destroy(mfi)
end do


END SUBROUTINE


!mfrho = rho_inf
!IF (primary_flow_dir == 1) THEN
!  mfu_vel = v_inf
!  mfv_vel = 0.0D0
!  mfw_vel = 0.0D0
!ELSE IF (primary_flow_dir == 2) THEN
!  mfu_vel = 0.0D0
!  mfv_vel = v_inf
!  mfw_vel = 0.0D0
!ELSE IF (primary_flow_dir == 3) THEN
!  mfu_vel = 0.0D0
!  mfv_vel = 0.0D0
!  mfw_vel = v_inf
!ELSE IF (primary_flow_dir == 4) THEN
!  mfu_vel = -v_inf
!  mfv_vel = 0.0D0
!  mfw_vel = 0.0D0
!ELSE IF (primary_flow_dir == 5) THEN
!  mfu_vel = 0.0D0
!  mfv_vel = -v_inf
!  mfw_vel = 0.0D0
!ELSE IF (primary_flow_dir == 6) THEN
!  mfu_vel = 0.0D0
!  mfv_vel = 0.0D0
!  mfw_vel = -v_inf
!END IF
!mftemp = temperature
!mfpress = pressure
!
!mfalpha = fs_alpha
!mffi = fs_fieq
!mffout = fs_fieq
!mfgi = fs_gieq
!mfgout = fs_gieq
!
!mfchi = fs_chi
!mfzeta_x = fs_zeta_x
!mfzeta_y = fs_zeta_y
!mfpi_xx = fs_pixx
!mfpi_xy = fs_pixy
!mfpi_yy = fs_piyy
!mflambda_x = fs_lambda_x
!mflambda_y = fs_lambda_y
!
!
!if (dimensions == 3) then
!  mfzeta_z = fs_zeta_z
!  mfpi_xz = fs_pixz
!  mfpi_zz = fs_pizz
!  mfpi_yz = fs_piyz
!  mflambda_z = fs_lambda_z
!end if

!do lvl = 0,fine
!
!
!
!end do
