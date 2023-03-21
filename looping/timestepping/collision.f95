SUBROUTINE collision(level_step,fine,lvl)
!
! This is the collisions step of LBM. This will definitely take the most time
! out of all the steps as the Lagrangian Multipliers must be calculated here and
! fi's for every direction must be found.
!
! Called by: loop
! Calls: calculate_lagrangian_mults,collide,find_weights
! External calls: amrex_mfiter_build,amrex_mfiter_destroy
!
use amrex_base_module
use mpi
use amrex_amr_module
use precise
USE grid_data
use quick_calcs
use amr_info_holder
use amr_processes
USE linkwise
use constants
IMPLICIT NONE
INTEGER :: fine,i,lvl
REAL(kind=dp) :: dummy_zeta_z,dummy_w,dummy_pi_xz,&
  dummy_pi_yz,dummy_pi_zz,dummy_lambda_z,vel_mag_2
REAL(kind=dp) :: lead_gieq_const,loc_mu,loc_kappa
REAL(kind=dp) :: beta_1,beta_2,second_coeff,loc_coords(dimensions)
REAL(kind=dp),allocatable :: loc_fieq(:),loc_gieq(:),loc_fi_diff(:),q_bar_dir(:),&
  contr_heat_flux(:),loc_gi_diff(:),old_fout(:)
real(kind=dp) :: rho_out,u_vel_out,v_vel_out,w_vel_out,temp_out,u_fieq,v_fieq,rho_fieq,&
  w_fieq,temp_fieq,pxx_eq,pyy_eq,pzz_eq,pxy_eq,pxz_eq,pyz_eq
real(kind=dp) :: loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz
LOGICAL :: level_step(0:amr_max_lvl)
INTEGER :: a,b,c,reach,null_counter

type(amrex_mfiter) :: mfi
type(amrex_box) :: box
!
! This is the workhorse subroutine where huge amounts of processing get done.
!
! As tempting as it is, the streaming and collision steps cannot be combined. This is
! because of the MPI method used by Amrex.  As the information must be shared between
! boxes, it has to be kept in its original box.  Then a fill_boundary call in streaming,
! fills the ghost nodes in so streaming can move the PDFs as needed.
!
!
if (self == 0) then
  write(*,*) 'starting collision step',lvl,self
end if
!
!
!
!
if (allocated(loc_fieq)) deallocate(loc_fieq)
allocate(loc_fieq(dir))
!
if (allocated(q_bar_dir)) deallocate(q_bar_dir)
allocate(q_bar_dir(dir))
!
if (allocated(loc_gieq)) deallocate(loc_gieq)
allocate(loc_gieq(dir))
!
if (allocated(contr_heat_flux)) deallocate(contr_heat_flux)
allocate(contr_heat_flux(dir))
!
if (allocated(loc_fi_diff)) deallocate(loc_fi_diff)
allocate(loc_fi_diff(dir))
!
if (allocated(loc_gi_diff)) deallocate(loc_gi_diff)
allocate(loc_gi_diff(dir))
!
allocate(old_fout(dir))

lead_gieq_const = 2.0D0*const_vol - dimensions

!DO lvl = fine,0,-1
!  IF (level_step(lvl)) THEN

CALL amrex_mfiter_build(mfi,mffi(lvl),tiling=.false.)
!call amrex_mfiter_build(mfj,mfrho(lvl),tiling=.true.)
if (.not. shifted) then
if (lvl == 0) then
  reach = nghosts/2
else
  if (mod(timestep(lvl),2) == 0) then
    reach = nghosts/2
  else
    reach = nghosts
  end if
end if
else
if (lvl == 0) then
  reach = nghosts/2
else
  if (mod(timestep(lvl),2) == 0) then
    reach = nghosts/2
  else
    reach = nghosts
  end if
end if
end if
!
! Go through all the boxes on a given level
!
do while(mfi%next())
!
  rho => mfrho(lvl)%dataptr(mfi)
  u_vel => mfu_vel(lvl)%dataptr(mfi)
  v_vel => mfv_vel(lvl)%dataptr(mfi)
  temp => mftemp(lvl)%dataptr(mfi)
  IF (dimensions == 2) THEN
    dummy_w = 0.0D0
  END IF
!
! Point at the LM data
!

  chi => mfchi(lvl)%dataptr(mfi)
  zeta_x => mfzeta_x(lvl)%dataptr(mfi)
  zeta_y => mfzeta_y(lvl)%dataptr(mfi)
  pi_xx => mfpi_xx(lvl)%dataptr(mfi)
  pi_yy => mfpi_yy(lvl)%dataptr(mfi)
  pi_xy => mfpi_xy(lvl)%dataptr(mfi)
  lambda_x => mflambda_x(lvl)%dataptr(mfi)
  lambda_y => mflambda_y(lvl)%dataptr(mfi)
  IF (dimensions == 3) THEN
    zeta_z => mfzeta_z(lvl)%dataptr(mfi)
!              dummy_zeta = zeta_z(a,b,c,1)
    w_vel => mfw_vel(lvl)%dataptr(mfi)
!              dummy_w = w_vel(a,b,c,1)
!
    pi_xz => mfpi_xz(lvl)%dataptr(mfi)
    pi_zz => mfpi_zz(lvl)%dataptr(mfi)
    pi_yz => mfpi_yz(lvl)%dataptr(mfi)
    lambda_z => mflambda_z(lvl)%dataptr(mfi)
  END IF
!
! Assign the big arrays to their pointers
!
  fi => mffi(lvl)%dataptr(mfi)
  fout => mffout(lvl)%dataptr(mfi)
  gi => mfgi(lvl)%dataptr(mfi)
  gout => mfgout(lvl)%dataptr(mfi)
  alpha => mfalpha(lvl)%dataptr(mfi)
  state => mfstate(lvl)%dataptr(mfi)
! The box allows us to get the ranges on the arrays
  box = mfi%validbox()
!
! Get the ranges on the box
!
!
! Sweep through those and do the collisions
!
  DO c = box%lo(3)-reach,box%hi(3)+reach
    DO b = box%lo(2)-reach,box%hi(2)+reach
      DO a = box%lo(1)-reach,box%hi(1)+reach

        if (state(a,b,c,1) <0) cycle
        if (state(a,b,c,1) >= 1000 .and. state(a,b,c,1) < 51000) cycle
        if (state(a,b,c,1) >= 200 .and. state(a,b,c,1) < 299) cycle
!
! Either set up the necessary terms for 3D or set them to 0 for 2D. This is passing
! to subroutines remains the same and some calculations don't have to change.
!
        if (dimensions == 3) then
          dummy_w = w_vel(a,b,c,1)

          dummy_zeta_z = zeta_z(a,b,c,1)
          dummy_lambda_z = lambda_z(a,b,c,1)
          dummy_pi_zz = pi_zz(a,b,c,1)
          dummy_pi_xz = pi_xz(a,b,c,1)
          dummy_pi_yz = pi_yz(a,b,c,1)
        else
          dummy_zeta_z = 0.0D0
          dummy_w = 0.0D0
          dummy_pi_xz = 0.0D0
          dummy_pi_yz = 0.0D0
          dummy_pi_zz = 0.0D0
          dummy_lambda_z = 0.0D0
        end if
        old_fout(1:dir) = fout(a,b,c,1:dir)
        vel_mag_2 = u_vel(a,b,c,1)**2 + v_vel(a,b,c,1)**2 + &
          dummy_w**2
!
!
!
       u_fieq = 0.0D0
       v_fieq = 0.0D0
       w_fieq = 0.0D0
       rho_fieq = 0.0D0
       temp_fieq = 0.0D0
       pxx_eq = 0.0D0
       pyy_eq = 0.0D0
       pzz_eq = 0.0D0
       pxy_eq = 0.0D0
       pxz_eq = 0.0D0
       pyz_eq = 0.0D0
       loc_pxx = 0.0D0
       loc_pyy = 0.0D0
       loc_pzz = 0.0D0
       loc_pxy = 0.0D0
       loc_pxz = 0.0D0
       loc_pyz = 0.0D0
!
! Create local arrays for fi, fieq,gi, and gieq as well as find K to find little_q
!
        if (.not. shifted) then
          DO i = 1,dir
!
! find your local fieq and gieq
!
            select case(dir)
              case (19)
                loc_fieq(i) = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                    zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                    pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                    pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
              case (39)
                loc_fieq(i) = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                    zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                    pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                    pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
            end select

            loc_fi_diff(i) = loc_fieq(i)-fi(a,b,c,i)
            loc_gieq(i) = lead_gieq_const*temp(a,b,c,1)*loc_fieq(i)
            loc_gi_diff(i) = loc_gieq(i)-gi(a,b,c,i)

          END DO
        else
          DO i = 1,dir
!
! find your local fieq and gieq
!
            select case(dir)
              case (19)
                loc_fieq(i) = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                    zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                    pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                    pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
              case (39)
                loc_fieq(i) = fieq_comp_39_shifted(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                    zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                    pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                    pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
            end select

            loc_fi_diff(i) = loc_fieq(i)-fi(a,b,c,i)
            loc_gieq(i) = lead_gieq_const*temp(a,b,c,1)*loc_fieq(i)
            loc_gi_diff(i) = loc_gieq(i)-gi(a,b,c,i)
          END DO
        end if

!
! Find the heat flux and first order contracted heat flux
!
        CALL find_q_bar(loc_fi_diff,loc_fieq,q_bar_dir,&
          u_vel(a,b,c,1),v_vel(a,b,c,1),dummy_w,temp(a,b,c,1),lvl)
        CALL find_little_q_bar(loc_gi_diff,u_vel(a,b,c,1),v_vel(a,b,c,1),&
          dummy_w,temp(a,b,c,1),contr_heat_flux,lvl)
!        q_bar_dir = 0.0D0
!        contr_heat_flux = 0.0D0
!
! Get local values for thermal conductivity (kappa) and viscosity (mu)
! Then use those to find beta 1 and beta 2
!
        CALL quick_sutherlands(loc_mu, temp(a,b,c,1), .false.)
        CALL quick_therm_cond(loc_kappa, rho(a,b,c,1), temp(a,b,c,1) )
!
! Calculate what is essentially a sub-grid scale model for the flow.
!
!        beta_1 = 1.0D0/((2.0D0*loc_mu)/(rho(a,b,c,1)*temp(a,b,c,1))+1.0D0)
!        beta_2 = 2.0D0/((2.0D0*loc_kappa)/(rho(a,b,c,1)*const_press*temp(a,b,c,1))+1.0D0)
!
!        if (timestep(fine) < 10) then
!          beta_1 = 0.9D0 + 0.01D0*timestep(fine)
!        end if
        beta_1 = beta_1_calculator(loc_mu,rho(a,b,c,1),temp(a,b,c,1),lvl,timestep(amr_max_lvl))
        beta_2 = beta_2_calculator(loc_kappa,rho(a,b,c,1),temp(a,b,c,1),lvl,timestep(amr_max_lvl))
!        beta_1 = 0.51D0
!        beta_2 = 0.51D0
!        second_coeff = 2.0D0*(beta_1-beta_2)
!
! Solve the collisions for fout
!

        DO i = 1,dir
!
!          fout(a,b,c,i) = fi(a,b,c,i) + alpha(a,b,c,1)*beta_1*loc_fi_diff(i) - second_coeff*q_bar_dir(i)
          fout(a,b,c,i) = fi(a,b,c,i) + alpha(a,b,c,1)*beta_1*(loc_fi_diff(i)-q_bar_dir(i)) + beta_2*q_bar_dir(i)
!
!
!          gout(a,b,c,i) = gi(a,b,c,i) + alpha(a,b,c,1)*beta_1*loc_gi_diff(i) - second_coeff*contr_heat_flux(i)
          gout(a,b,c,i) = gi(a,b,c,i) + alpha(a,b,c,1)*beta_1*(loc_gi_diff(i) - contr_heat_flux(i)) + beta_2*contr_heat_flux(i)
!
!          if (self == 0) then
!          if (fout(a,b,c,i) < 0.0D0 .or. gout(a,b,c,i) < 0.0D0 .or. fout(a,b,c,i) >1.0D0 .or. &
!              gout(a,b,c,i) > 1.0D0 .or. fi(a,b,c,i) <0.0D0 .or. gi(a,b,c,i) < 0.0D0 .or. &
!              fi(a,b,c,i) > 1.0D0 .or. gi(a,b,c,i) > 1.0D0 &
!              .or. isnan(fi(a,b,c,i)) .or. isnan(gi(a,b,c,i)) .or. isnan(fout(a,b,c,i)) &
!              .or. isnan(gout(a,b,c,i))) then
!
!            !loc_coords = get_real_coords(amrex_geom(lvl),(/a,b,c/))
!
!           !if (state(a,b,c,1) < 1000) then
!            !if (self == 0) then
!            write(*,*) 'collisional values for negative fi or gi',fi(a,b,c,i),fout(a,b,c,i),loc_fieq(i),loc_fi_diff(i),&
!              gi(a,b,c,i),gout(a,b,c,i),loc_gieq(i),loc_gi_diff(i),&
!              i,lvl,timestep(lvl),time(lvl),box%lo,box%hi,a,b,c,self
!            write(*,*) 'node information, negative fi or gi',state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i)),&
!              state(a,b,c,1:dir),a,b,c,lvl
!            write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
!              alpha(a,b,c,1),a,b,c
!           write(*,*) 'more components',alpha(a,b,c,1),beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const,a,b,c
!!                write(*,*) 'local fi diffs',loc_fi_diff,a,b,c
!!                write(*,*) 'local gi diffs',loc_gi_diff,a,b,c
!            write(*,*) 'local fis',fi(a,b,c,1:dir),a,b,c
!            write(*,*) 'fi diffs',loc_fi_diff(1:dir),a,b,c
!            write(*,*) 'local fieqs',loc_fieq,a,b,c
!            write(*,*) 'fouts',fout(a,b,c,1:dir),a,b,c
!            write(*,*) 'local gis',gi(a,b,c,1:dir),a,b,c
!            write(*,*) 'local gieqs',loc_gieq,a,b,c
!            write(*,*) 'local gi diffs',loc_gi_diff(1:dir),a,b,c
!            write(*,*) 'gouts',gout(a,b,c,1:dir),a,b,c
!            write(*,*) 'Q bar',q_bar_dir,sum(q_bar_dir(1:dir)),a,b,c
!            write(*,*) 'Contracted heat flux',contr_heat_flux,sum(contr_heat_flux(1:dir)),a,b,c
!            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!                      zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                      pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!                      pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),a,b,c
!            end if
!          end if
!          if (isnan(fi(a,b,c,i)) .or. isnan(gi(a,b,c,i)) .or. isnan(fout(a,b,c,i)) &
!              .or. isnan(gout(a,b,c,i))) then
!
!            !loc_coords = get_real_coords(amrex_geom(lvl),(/a,b,c/))
!
!           !if (state(a,b,c,1) < 1000) then
!            !if (self == 0) then
!            write(*,*) 'collisional values for negative fi or gi',fi(a,b,c,i),fout(a,b,c,i),loc_fieq(i),loc_fi_diff(i),&
!              gi(a,b,c,i),gout(a,b,c,i),loc_gieq(i),loc_gi_diff(i),&
!              i,lvl,timestep(lvl),time(lvl),box%lo,box%hi,a,b,c,self
!            write(*,*) 'node information, negative fi or gi',state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i)),&
!              state(a,b,c,1:dir),a,b,c,lvl
!            write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
!              alpha(a,b,c,1),a,b,c
!           write(*,*) 'more components',alpha(a,b,c,1),beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const,a,b,c
!!                write(*,*) 'local fi diffs',loc_fi_diff,a,b,c
!!                write(*,*) 'local gi diffs',loc_gi_diff,a,b,c
!            write(*,*) 'local fis',fi(a,b,c,1:dir),a,b,c
!            write(*,*) 'fi diffs',loc_fi_diff(1:dir),a,b,c
!            write(*,*) 'local fieqs',loc_fieq,a,b,c
!            write(*,*) 'fouts',fout(a,b,c,1:dir),a,b,c
!            write(*,*) 'local gis',gi(a,b,c,1:dir),a,b,c
!            write(*,*) 'local gieqs',loc_gieq,a,b,c
!            write(*,*) 'local gi diffs',loc_gi_diff(1:dir),a,b,c
!            write(*,*) 'gouts',gout(a,b,c,1:dir),a,b,c
!            write(*,*) 'Q bar',q_bar_dir,sum(q_bar_dir(1:dir)),a,b,c
!            write(*,*) 'Contracted heat flux',contr_heat_flux,sum(contr_heat_flux(1:dir)),a,b,c
!            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!                      zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                      pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!                      pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),a,b,c
!            end if
!          end if

!          end if
          if (a == 1 .and. b == 32 .and. c == 32 .and. i == 39) then
          !if (state(a,b,c,1) < 51000 .and. state(a,b,c,1) >1000 .and. lvl == 3) then
            write(*,*) 'node oddities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
              temp(a,b,c,1),alpha(a,b,c,1),a,b,c,timestep(lvl),self
            !write(*,*) 'out velocities and stuff',rho_out,u_vel_out,v_vel_out,w_vel_out,temp_out
!            write(*,*) 'pressure terms',loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,pxx_eq,pyy_eq,pzz_eq,&
!              pxy_eq,pxz_eq,pyz_eq
            write(*,*) 'more components',alpha(a,b,c,1),beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const,a,b,c
            write(*,*) 'f-terms',fout(a,b,c,i),fi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_fi_diff(i),q_bar_dir(i),beta_2
            write(*,*) 'g-terms',gout(a,b,c,i),gi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_gi_diff(i),contr_heat_flux(i),beta_2
            write(*,*) 'oddball states',state(a,b,c,1:dir),box%lo,box%hi,a,b,c
            write(*,*) 'oddball fis',fi(a,b,c,1:dir)
            write(*,*) 'oddball fouts',fout(a,b,c,1:dir)
            write(*,*) 'fi diffs',loc_fi_diff
            write(*,*) 'fieqs',loc_fieq
            write(*,*) 'nearby states',state(a+1,b,c,1),state(a,b+1,c,1),state(a,b,c+1,1), state(a-1,b,c,1),&
              state(a,b-1,c,1),state(a,b,c-1,1)
            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
            write(*,*) 'oddball gis',gi(a,b,c,1:dir)
            write(*,*) 'oddball gouts',gout(a,b,c,1:dir)
            write(*,*) 'gi diffs',loc_gi_diff
            write(*,*) 'gieqs',loc_gieq
            write(*,*) 'Q bars',q_bar_dir,sum(q_bar_dir),elbm_epsilon(lvl)
            write(*,*) 'contr heat flux',contr_heat_flux,sum(contr_heat_flux),elbm_epsilon(lvl)
          end if

!          if (a == 293 .and. b == 105 .and. c == 127 .and. i == 39) then
!            write(*,*) 'node oddities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!              temp(a,b,c,1),alpha(a,b,c,1),a,b,c,self
!            !write(*,*) 'out velocities and stuff',rho_out,u_vel_out,v_vel_out,w_vel_out,temp_out
!!            write(*,*) 'pressure terms',loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,pxx_eq,pyy_eq,pzz_eq,&
!!              pxy_eq,pxz_eq,pyz_eq
!            write(*,*) 'more components',alpha(a,b,c,1),beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const,a,b,c
!            write(*,*) 'f-terms',fout(a,b,c,i),fi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_fi_diff(i),q_bar_dir(i),beta_2
!            write(*,*) 'g-terms',gout(a,b,c,i),gi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_gi_diff(i),contr_heat_flux(i),beta_2
!            write(*,*) 'oddball states',state(a,b,c,1:dir),box%lo,box%hi,a,b,c
!            write(*,*) 'oddball fis',fi(a,b,c,1:dir)
!            write(*,*) 'oddball fouts',fout(a,b,c,1:dir)
!            write(*,*) 'fi diffs',loc_fi_diff
!            write(*,*) 'fieqs',loc_fieq
!            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!            write(*,*) 'oddball gis',gi(a,b,c,1:dir)
!            write(*,*) 'oddball gouts',gout(a,b,c,1:dir)
!            write(*,*) 'gi diffs',loc_gi_diff
!            write(*,*) 'gieqs',loc_gieq
!            write(*,*) 'Q bars',q_bar_dir,sum(q_bar_dir),elbm_epsilon(lvl)
!            write(*,*) 'contr heat flux',contr_heat_flux,sum(contr_heat_flux),elbm_epsilon(lvl)
!          end if

!          if (a == 296 .and. b == 105 .and. c == 3 .and. i == 39) then
!            write(*,*) 'node oddities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!              temp(a,b,c,1),alpha(a,b,c,1),a,b,c,self
!            !write(*,*) 'out velocities and stuff',rho_out,u_vel_out,v_vel_out,w_vel_out,temp_out
!!            write(*,*) 'pressure terms',loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,pxx_eq,pyy_eq,pzz_eq,&
!!              pxy_eq,pxz_eq,pyz_eq
!            write(*,*) 'more components',alpha(a,b,c,1),beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const,a,b,c
!            write(*,*) 'f-terms',fout(a,b,c,i),fi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_fi_diff(i),q_bar_dir(i),beta_2
!            write(*,*) 'g-terms',gout(a,b,c,i),gi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_gi_diff(i),contr_heat_flux(i),beta_2
!            write(*,*) 'oddball states',state(a,b,c,1:dir),box%lo,box%hi,a,b,c
!            write(*,*) 'oddball fis',fi(a,b,c,1:dir)
!            write(*,*) 'oddball fouts',fout(a,b,c,1:dir)
!            write(*,*) 'fi diffs',loc_fi_diff
!            write(*,*) 'fieqs',loc_fieq
!            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!            write(*,*) 'oddball gis',gi(a,b,c,1:dir)
!            write(*,*) 'oddball gouts',gout(a,b,c,1:dir)
!            write(*,*) 'gi diffs',loc_gi_diff
!            write(*,*) 'gieqs',loc_gieq
!            write(*,*) 'Q bars',q_bar_dir,sum(q_bar_dir),elbm_epsilon(lvl)
!            write(*,*) 'contr heat flux',contr_heat_flux,sum(contr_heat_flux),elbm_epsilon(lvl)
!          end if

!          if (a == 345 .and. b == 249 .and. c == 3 .and. i == 39) then
!            write(*,*) 'node oddities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!              temp(a,b,c,1),alpha(a,b,c,1),a,b,c,self
!            !write(*,*) 'out velocities and stuff',rho_out,u_vel_out,v_vel_out,w_vel_out,temp_out
!!            write(*,*) 'pressure terms',loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,pxx_eq,pyy_eq,pzz_eq,&
!!              pxy_eq,pxz_eq,pyz_eq
!            write(*,*) 'more components',alpha(a,b,c,1),beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const,a,b,c
!            write(*,*) 'f-terms',fout(a,b,c,i),fi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_fi_diff(i),q_bar_dir(i),beta_2
!            write(*,*) 'g-terms',gout(a,b,c,i),gi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_gi_diff(i),contr_heat_flux(i),beta_2
!            write(*,*) 'oddball states',state(a,b,c,1:dir),box%lo,box%hi,a,b,c
!            write(*,*) 'oddball fis',fi(a,b,c,1:dir)
!            write(*,*) 'oddball fouts',fout(a,b,c,1:dir)
!            write(*,*) 'fi diffs',loc_fi_diff
!            write(*,*) 'fieqs',loc_fieq
!            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!            write(*,*) 'oddball gis',gi(a,b,c,1:dir)
!            write(*,*) 'oddball gouts',gout(a,b,c,1:dir)
!            write(*,*) 'gi diffs',loc_gi_diff
!            write(*,*) 'gieqs',loc_gieq
!            write(*,*) 'Q bars',q_bar_dir,sum(q_bar_dir),elbm_epsilon(lvl)
!            write(*,*) 'contr heat flux',contr_heat_flux,sum(contr_heat_flux),elbm_epsilon(lvl)
!          end if

!          if (a == 114 .and. b == 64 .and. c == 64 .and. i == 39) then
!            write(*,*) 'node oddities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!              temp(a,b,c,1),alpha(a,b,c,1),a,b,c,self
!            !write(*,*) 'out velocities and stuff',rho_out,u_vel_out,v_vel_out,w_vel_out,temp_out
!!            write(*,*) 'pressure terms',loc_pxx,loc_pyy,loc_pzz,loc_pxy,loc_pxz,loc_pyz,pxx_eq,pyy_eq,pzz_eq,&
!!              pxy_eq,pxz_eq,pyz_eq
!            write(*,*) 'f-terms',fout(a,b,c,i),fi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_fi_diff(i),q_bar_dir(i),beta_2
!            write(*,*) 'g-terms',gout(a,b,c,i),gi(a,b,c,i),alpha(a,b,c,1),beta_1,loc_gi_diff(i),contr_heat_flux(i),beta_2
!            write(*,*) 'oddball states',state(a,b,c,1:dir),box%lo,box%hi,a,b,c
!            write(*,*) 'oddball fis',fi(a,b,c,1:dir)
!            write(*,*) 'oddball fouts',fout(a,b,c,1:dir)
!            write(*,*) 'fi diffs',loc_fi_diff
!            write(*,*) 'fieqs',loc_fieq
!            write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!            write(*,*) 'oddball gis',gi(a,b,c,1:dir)
!            write(*,*) 'oddball gouts',gout(a,b,c,1:dir)
!            write(*,*) 'gi diffs',loc_gi_diff
!            write(*,*) 'gieqs',loc_gieq
!            write(*,*) 'Q bars',q_bar_dir,sum(q_bar_dir),elbm_epsilon(lvl)
!            write(*,*) 'contr heat flux',contr_heat_flux,sum(contr_heat_flux),elbm_epsilon(lvl)
!          end if

        END DO

!        if (a == 108 .and. b == 72 .and. c == 38 .and. lvl == 1) then
!          write(*,*) 'state police!',state(a,b,c,1:dir),amrex_geom(lvl)%dx
!        end if
!
!
!        if (state(a,b,c,1) >= 200 .and. state(a,b,c,1) <= 299) then
!          fi(a,b,c,1:dir) = fout(a,b,c,1:dir)
!          gi(a,b,c,1:dir) = gout(a,b,c,1:dir)
!        end if

      END DO
    END DO
  END DO

end do
!
! Destroy the mfiter and move on to streaming
!
call amrex_mfiter_destroy(mfi)
!    call mffout(lvl)%fill_boundary(amrex_geom(lvl),1,dir,.true.)

    !call amrex_mfiter_destroy(mfj)
!  ELSE
!    CYCLE
!  END IF
!END DO

DEALLOCATE(loc_fieq)
DEALLOCATE(loc_gieq)
DEALLOCATE(loc_fi_diff)
deallocate(loc_gi_diff)
deallocate(q_bar_dir)
deallocate(contr_heat_flux)
END SUBROUTINE


!
!              if (lvl ==3 ) then
!                write(*,*) 'fi information',loc_fieq(i),loc_fi_diff(i),fi(a,b,c,i)
!              end if
!                    loc_gi_diff(i) = loc_gieq(i)-loc_gi(i)
!            P_xx = rho(a,b,c,1)*temp(a,b,c,1) + rho(a,b,c,1)*u_vel(a,b,c,1)**2
!            P_yy = rho(a,b,c,1)*temp(a,b,c,1) + rho(a,b,c,1)*v_vel(a,b,c,1)**2
!            P_xy = rho(a,b,c,1)*v_vel(a,b,c,1)*u_vel(a,b,c,1)
!            if (dimensions == 3) then
!              P_zz = rho(a,b,c,1)*temp(a,b,c,1) + rho(a,b,c,1)*w_vel(a,b,c,1)**2
!              P_xz = rho(a,b,c,1)*w_vel(a,b,c,1)*u_vel(a,b,c,1)
!              P_yz = rho(a,b,c,1)*v_vel(a,b,c,1)*w_vel(a,b,c,1)
!            end if
!
! Find all 13 Lagrange multipliers
!

!                  little_q_bar(1) = (2.0D0*rho(a,b,c,1)*const_vol*temp(a,b,c,1) + &
!                    2.0D0*rho(a,b,c,1)*temp(a,b,c,1) + 2.0D0*rho(a,b,c,1)*vel_mag_2)*u_vel(a,b,c,1)
!                  little_q_bar(2) = (2.0D0*rho(a,b,c,1)*const_vol*temp(a,b,c,1) + &
!                    2.0D0*rho(a,b,c,1)*temp(a,b,c,1) + 2.0D0*rho(a,b,c,1)*vel_mag_2)*v_vel(a,b,c,1)
!                  IF (dimensions == 3) THEN
!                    little_q_bar(3) = (2.0D0*rho(a,b,c,1)*const_vol*temp(a,b,c,1) + &
!                      2.0D0*rho(a,b,c,1)*temp(a,b,c,1) + 2.0D0*rho(a,b,c,1)*vel_mag_2)*w_vel(a,b,c,1)
!                  END IF

!              if (state(a,b,c,i) >=0) then
!                fout(a+cx(i),b+cy(i),c+cz(i),i) = fi(a,b,c,i)+beta_1*(loc_fi_diff(i)-q_bar_dir(i)) - &
!                  beta_2*q_bar_dir(i)
!                gout(a+cx(i),b+cy(i),c+cz(i),i) = gi(a,b,c,i)+beta_1*(loc_gieq(i)-gi(a,b,c,i)+contr_heat_flux(i)) - &
!                  beta_2*contr_heat_flux(i)
!              else if (state(a,b,c,i) <= -100 .and. state(a,b,c,i) >= -199) then
!                call inlet_bc(fout(a,b,c,i),i,state(a,b,c,i))
!              else if (state(a,b,c,i) <= -200 .and. state(a,b,c,i) >= -299) then
!                call outlet_bc(fout(a,b,c,i),i,state(a,b,c,i))
!              else if (state(a,b,c,i) < 0 .and. state(a,b,c,i) > -100) then
!                call bc_switchboard(fout(a,b,c,i),i)


!              end if
!                    fistar = loc_fieq(i)+q
               !CALL collide(fi(a,b,c,i),fout(a,b,c,i),gi(a,b,c,i),gout(a,b,c,i),&
               !  loc_fi_diff(i),&
               !  rho(a,b,c,i),u_vel(a,b,c,i),v_vel(a,b,c,i),dummy_w_vel,temp(a,b,c,i),&
               !  i,vel_mag,alpha(a,b,c,1),beta_1,beta_2)

!              if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0 .or. fi(a,b,c,i) >1.0D0 .or. &
!                  gi(a,b,c,i) > 1.0D0) then
!                write(*,*) 'collisional values',fi(a,b,c,i),fout(a,b,c,i),loc_fieq(i),&
!                  gi(a,b,c,i),gout(a,b,c,i),loc_gieq(i),i,lvl,timestep(lvl)
!                write(*,*) 'node information',state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i))
!                write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
!                write(*,*) 'Q bar',q_bar_dir
!                write(*,*) 'Contracted heat flux',contr_heat_flux
!                write(*,*) 'more components',beta_1,beta_2,loc_mu,loc_kappa,lead_gieq_const
!                write(*,*) 'collision LMs',chi(a,b,c,1),zeta_x(a,b,c,1),&
!                          zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                          pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!                          pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!              end if

!          if (a == 289 .and. b == 132 .and. c == 125 .and. lvl == 3 .and. i == 39) then
!            write(*,*) 'node oddities',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!              temp(a,b,c,1),a,b,c,self
!            write(*,*) 'oddball states',state(a,b,c,1:dir)
!            write(*,*) 'oddball fis',fi(a,b,c,1:dir)
!            write(*,*) 'oddball fouts',fout(a,b,c,1:dir)
!            write(*,*) 'fieqs',loc_fieq
!            write(*,*) 'oddball gis',gi(a,b,c,1:dir)
!            write(*,*) 'oddball gouts',gout(a,b,c,1:dir)
!            write(*,*) 'gieqs',loc_gieq
!            write(*,*) 'Q bars',q_bar_dir
!            write(*,*) 'contr heat flux',contr_heat_flux
!          end if
!          if (a == 277 .and. b == 132 .and. c == 125 .and. lvl == 3 .and. i == 39) then
!            write(*,*) 'node oddities -3',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
!              temp(a,b,c,1),a,b,c,self
!            write(*,*) 'oddball states -3',state(a,b,c,1:dir)
!            write(*,*) 'oddball fis -3',fi(a,b,c,1:dir)
!            write(*,*) 'oddball fouts -3',fout(a,b,c,1:dir)
!            write(*,*) 'fieqs -3',loc_fieq
!            write(*,*) 'oddball gis -3',gi(a,b,c,1:dir)
!            write(*,*) 'oddball gouts -3',gout(a,b,c,1:dir)
!            write(*,*) 'gieqs -3',loc_gieq
!            write(*,*) 'Q bars -3',q_bar_dir
!            write(*,*) 'contr heat flux -3',contr_heat_flux
!          end if
!
! Save the new zeta_z and w_vel to the arrays before it's destroyed
!
!            IF (dimensions == 3) THEN
!              zeta_z(a,b,c,1) = dummy_zeta
!              w_vel(a,b,c,1) = dummy_w
!              pi_zz(a,b,c,1) = dummy_pi_xx
!              pi_xz(a,b,c,1) = dummy_pi_xz
!              pi_yz(a,b,c,1) = dummy_pi_yz
!              zeta_z(a,b,c,1) = dummy_zeta_z
!              lambda_z(a,b,c,1) = dummy_lambda_z
!            END IF

!                    loc_fieq(i) = rho(a,b,c,i)*EXP(chi(a,b,c,1) + zeta_x(a,b,c,1)*rcx(i) +&
!                      zeta_y(a,b,c,1)*rcy(i) + dummy_zeta_z*rcz(i) + pi_xx(a,b,c,1)*rcx(i)**2 +&
!                      pi_yy(a,b,c,1)*rcy(i)**2 + pi_xy(a,b,c,1)*rcx(i)*rcy(i) + dummy_pi_zz*rcz(i)**2 +&
!                      dummy_pi_xz*rcx(i)*rcy(i) + dummy_pi_yz*rcy(i)*rcz(i) + &
!                      lambda_x(a,b,c,1)*cmag(i)*rcx(i) + lambda_y(a,b,c,1)*cmag(i)*rcy(i) + &
!                      dummy_lambda_z*cmag(i)*rcz(i))


!        u_fieq = u_fieq/rho_fieq
!        v_fieq = v_fieq/rho_fieq
!        w_fieq = w_fieq/rho_fieq
!        temp_fieq = (temp_fieq - rho_fieq*(u_fieq**2+v_fieq**2+w_fieq**2))/(2*const_vol*rho_fieq)
!        rho_out = 0.0D0
!        u_vel_out= 0.0D0
!        v_vel_out = 0.0D0
!        w_vel_out = 0.0D0
!        temp_out = 0.0D0
!
!
!
!          if (fi(a,b,c,i) < 0.0D0 .or. fi(a,b,c,i) > 1.0D0 .or. loc_fieq(i) < 0.0D0 .or. &
!            loc_fieq(i) > 1.0D0 .or. fout(a,b,c,i) < 0.0D0 .or. fout(a,b,c,i) > 1.0D0) then
!            write(*,*) 'bleh, really fi? ',fi(a,b,c,i),loc_fieq(i),fout(a,b,c,i),&
!              fout(a,b,c,opp(i)),state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i)),&
!              i,lvl,timestep(lvl)
!          end if
!!
!!
!!
!          if (gi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) > 1.0D0 .or. loc_gieq(i) < 0.0D0 .or. &
!            loc_gieq(i) > 1.0D0 .or. gout(a,b,c,i) < 0.0D0 .or. gout(a,b,c,i) > 1.0D0) then
!            write(*,*) 'bleh, really gi? ',gi(a,b,c,i),loc_gieq(i),gout(a,b,c,i),&
!              gout(a,b,c,opp(i)),state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i)),&
!              i,lvl,timestep(lvl)
!          end if
!
!          u_fieq = u_fieq + rcx(i)*loc_fieq(i)
!          v_fieq = v_fieq + rcy(i)*loc_fieq(i)
!          w_fieq = w_fieq + rcz(i)*loc_fieq(i)
!          rho_fieq = rho_fieq + loc_fieq(i)
!          temp_fieq = temp_fieq + cmag(i)*loc_fieq(i)+loc_gieq(i)
!          pxx_eq = pxx_eq + rcx(i)**2*loc_fieq(i)
!          pyy_eq = pyy_eq +rcy(i)**2*loc_fieq(i)
!          pzz_eq = pzz_eq +rcz(i)**2*loc_fieq(i)
!          pxy_eq = pxy_eq +rcx(i)*rcy(i)*loc_fieq(i)
!          pxz_eq = pxz_eq +rcx(i)*rcz(i)*loc_fieq(i)
!          pyz_eq = pyz_eq +rcy(i)*rcz(i)*loc_fieq(i)
!
!          loc_pxx = loc_pxx + rcx(i)**2*fi(a,b,c,i)
!          loc_pyy = loc_pyy + rcy(i)**2*fi(a,b,c,i)
!          loc_pzz = loc_pzz + rcz(i)**2*fi(a,b,c,i)
!          loc_pxy = loc_pxy + rcx(i)*rcy(i)*fi(a,b,c,i)
!          loc_pxz = loc_pxz + rcx(i)*rcz(i)*fi(a,b,c,i)
!          loc_pyz = loc_pyz + rcy(i)*rcz(i)*fi(a,b,c,i)


!          if (abs(old_fout(i)-fout(a,b,c,i)) > 0.05D0) then
!            write(*,*) 'things are going wrong',old_fout(i),fout(a,b,c,i),i,fi(a,b,c,i),rho(a,b,c,1),u_vel(a,b,c,1),&
!              v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),loc_fieq(i),loc_fi_diff(i),state(a,b,c,1),state(a,b,c,i)
!
!          end if
!          rho_out = rho_out + fout(a,b,c,i)
!          u_vel_out = u_vel_out + rcx(i)*fout(a,b,c,i)
!          v_vel_out = v_vel_out + rcy(i)*fout(a,b,c,i)
!          w_vel_out = w_vel_out + rcz(i)*fout(a,b,c,i)
!          temp_out = temp_out + gout(a,b,c,i)
