SUBROUTINE streaming(level_step,fine,lvl)
  !
  ! This is the streaming step where the fouts move in their direction
  ! and become the fis for the next timestep
  !
  ! Called by: loop
  ! Calls:
  ! External calls:
  !
  use precise
  use constants
  use amrex_base_module
  use amrex_amr_module
  use amr_info_holder
  use amr_processes
  use quick_calcs
  use bc_info
  !use mpi
  use freestream_values
  use linkwise
  use grid_data, only: dir
  IMPLICIT NONE
  integer :: fine,max_level,coarsest_lvl_run,lvl
  integer :: a,b,c,i,j,reach,k,ier,loc_cx,loc_cy,loc_cz,state_pass(7)
  logical :: level_step(:)
  real(kind=dp) :: loc_fieq,loc_gieq,loc_fieq_opp,loc_gieq_opp,k_gieq,pass_fieq(dir)
  real(kind=dp),allocatable :: u_pass(:),v_pass(:),w_pass(:),temp_pass(:)
  real(kind=dp),allocatable :: u_local(:),v_local(:),w_local(:),temp_local(:)
  !real(kind=dp),contiguous,pointer :: fi_temp(:,:,:,:)
  real(kind=dp),contiguous,pointer :: u_tgt_temp(:,:,:,:),v_tgt_temp(:,:,:,:),&
    w_tgt_temp(:,:,:,:),temp_tgt_temp(:,:,:,:),rho_tgt_temp(:,:,:,:)
  real(kind=dp),contiguous,pointer :: chi_temp(:,:,:,:),zeta_x_temp(:,:,:,:),zeta_y_temp(:,:,:,:),&
    zeta_z_temp(:,:,:,:),pi_xx_temp(:,:,:,:),pi_yy_temp(:,:,:,:),pi_zz_temp(:,:,:,:),pi_xy_temp(:,:,:,:),&
    pi_xz_temp(:,:,:,:),pi_yz_temp(:,:,:,:),lambda_x_temp(:,:,:,:),lambda_y_temp(:,:,:,:),&
    lambda_z_temp(:,:,:,:)
  real(kind=dp) :: pxx_temp,pyy_temp,pzz_temp,pxy_temp,pxz_temp,pyz_temp
  real(kind=dp) :: dummy_fout,dummy_gout
  integer,contiguous,pointer :: node_pt(:,:,:,:)!,state_temp(:,:,:,:)


  type(amrex_mfiter) :: mfi
  type(amrex_box) :: wax
  type(amrex_multifab) :: mfu_tgt_temp,mfv_tgt_temp,mfw_tgt_temp,mftemp_tgt_temp,mfrho_tgt_temp
  type(amrex_multifab) :: mf_chi_temp,mf_zeta_x_temp,mf_zeta_y_temp,mf_zeta_z_temp,mf_pi_xx_temp,&
    mf_pi_yy_temp,mf_pi_zz_temp,mf_pi_xy_temp,mf_pi_xz_temp,mf_pi_yz_temp,mf_lambda_x_temp,&
    mf_lambda_y_temp,mf_lambda_z_temp
  type(amrex_imultifab) :: imf_node_done
  !
  !
  !
  !
  !
  k_gieq = 2.0D0*const_vol-dimensions
  !
  !
  !
  if (allocated(u_pass)) deallocate(u_pass)
  if (allocated(v_pass)) deallocate(v_pass)
  if (allocated(w_pass)) deallocate(w_pass)
  if (allocated(temp_pass)) deallocate(temp_pass)

  allocate(u_pass(dir))
  allocate(v_pass(dir))
  allocate(w_pass(dir))
  allocate(temp_pass(dir))
  !
  ! targets
  ! 1 = u velocity
  ! 2 = v velocity
  ! 3 = w_velocity
  ! 4 = temperature
  !
  ! In case of 2-D, 3 = temperature
  !
  !if (dimensions == 3) then
  !  num_targets = 4
  !  allocate(mf_targets(num_targets))
  !else
  !  num_targets = 3
  !  allocate(mf_targets(num_targets))
  !end if

  !do lvl = fine,coarsest_lvl_run,-1
  if (self == 0) then
    write(*,*) 'beginning streaming on level ',lvl,time!,elbm_epsilon
  end if
  !
  ! Build the multifabs so things can get saved.
  !
  call amrex_multifab_build(mfrho_tgt_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfu_tgt_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfv_tgt_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mfw_tgt_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mftemp_tgt_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)

  call amrex_multifab_build(mf_chi_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)

  call amrex_multifab_build(mf_zeta_x_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_zeta_y_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_zeta_z_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)

  call amrex_multifab_build(mf_pi_xx_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_pi_yy_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_pi_zz_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_pi_xy_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_pi_xz_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_pi_yz_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)

  call amrex_multifab_build(mf_lambda_x_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_lambda_y_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  call amrex_multifab_build(mf_lambda_z_temp,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  !
  ! save whether a node has been done before and therefore the target values exists
  !
  call amrex_imultifab_build(imf_node_done,mffi(lvl)%ba,mffi(lvl)%dm,1,nghosts_mid,node_based_3d)
  !
  !
  !
  !  write(*,*) 'wiggle',self,lvl
  call mffout(lvl)% fill_boundary( amrex_geom(lvl),1,dir)
  !  write(*,*) 'waggle',self,lvl
  call mfgout(lvl)% fill_boundary(amrex_geom(lvl),1,dir)
  !  write(*,*) 'woggle',self,lvl
  !
  call mfrho(lvl)% fill_boundary(amrex_geom(lvl))
  call mfu_vel(lvl)% fill_boundary(amrex_geom(lvl))
  call mfv_vel(lvl)% fill_boundary(amrex_geom(lvl))
  call mftemp(lvl)% fill_boundary(amrex_geom(lvl))
  !
  call mfchi(lvl)% fill_boundary(amrex_geom(lvl))
  call mfzeta_x(lvl)% fill_boundary(amrex_geom(lvl))
  call mfzeta_y(lvl)% fill_boundary(amrex_geom(lvl))
  call mfpi_xx(lvl)% fill_boundary(amrex_geom(lvl))
  call mfpi_xy(lvl)% fill_boundary(amrex_geom(lvl))
  call mfpi_yy(lvl)% fill_boundary(amrex_geom(lvl))
  call mflambda_x(lvl)% fill_boundary(amrex_geom(lvl))
  call mflambda_y(lvl)% fill_boundary(amrex_geom(lvl))

  if (dimensions == 3) then
    call mfw_vel(lvl)% fill_boundary(amrex_geom(lvl))
    call mfzeta_z(lvl)% fill_boundary(amrex_geom(lvl))
    call mfpi_xz(lvl)% fill_boundary(amrex_geom(lvl))
    call mfpi_yz(lvl)% fill_boundary(amrex_geom(lvl))
    call mfpi_zz(lvl)% fill_boundary(amrex_geom(lvl))
    call mflambda_z(lvl)% fill_boundary(amrex_geom(lvl))

    if (allocated(u_local)) deallocate(u_local)
    if (allocated(v_local)) deallocate(v_local)
    if (allocated(w_local)) deallocate(w_local)
    if (allocated(temp_local)) deallocate(temp_local)

    allocate(u_local(7))
    allocate(v_local(7))
    allocate(w_local(7))
    allocate(temp_local(7))

  end if

  !
  ! Build the iterator so things can point
  !
  call amrex_mfiter_build(mfi,mffout(lvl),tiling=.false.)
  if (.not. shifted) then
    if (lvl == 0) then
      reach = nghosts/2
    !  else if (mod(timestep(lvl),2) == 1) then
    !    reach = nghosts/2
    else
      if (mod(timestep(lvl),2) == 0) then
        reach = nghosts
      else
        reach = nghosts/2
      end if
    end if
  else
    if (lvl == 0) then
      reach = nghosts/2
    !  else if (mod(timestep(lvl),2) == 1) then
    !    reach = nghosts/2
    else
      if (mod(timestep(lvl),2) == 0) then
        reach = nghosts
      else
        reach = nghosts/2
      end if
    end if
  end if
  !
  !  do i = 1,dir
  !    call mffout(lvl)%fill_boundary(amrex_geom(lvl),.true.)
  !  end do
  !
  do while(mfi%next())
    wax=mfi%validbox()
    !
    ! Target multifabs
    !
    !    write(*,*) 'location of multifabs',loc(mfu_tgt_temp),loc(mfv_tgt_temp),loc(imf_node_done)
    state => mfstate(lvl)%dataptr(mfi)
    !    state => mfstate(lvl)%dataptr(mfi)

    fi => mffi(lvl)%dataptr(mfi)
    fout => mffout(lvl)%dataptr(mfi)
    gi => mfgi(lvl)%dataptr(mfi)
    gout => mfgout(lvl)%dataptr(mfi)
    !
    !
    !
    rho => mfrho(lvl)%dataptr(mfi)
    u_vel => mfu_vel(lvl)%dataptr(mfi)
    v_vel => mfv_vel(lvl)%dataptr(mfi)
    temp => mftemp(lvl)%dataptr(mfi)
    alpha => mfalpha(lvl)%dataptr(mfi)
    !
    ! interpolative targets
    !
    rho_tgt_temp => mfrho_tgt_temp%dataptr(mfi)
    u_tgt_temp => mfu_tgt_temp%dataptr(mfi)
    v_tgt_temp => mfv_tgt_temp%dataptr(mfi)
    w_tgt_temp => mfw_tgt_temp%dataptr(mfi)
    temp_tgt_temp => mftemp_tgt_temp%dataptr(mfi)
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
    IF (dimensions == 3) THEN
      zeta_z => mfzeta_z(lvl)%dataptr(mfi)
      w_vel => mfw_vel(lvl)%dataptr(mfi)
      pi_xz => mfpi_xz(lvl)%dataptr(mfi)
      pi_zz => mfpi_zz(lvl)%dataptr(mfi)
      pi_yz => mfpi_yz(lvl)%dataptr(mfi)
      lambda_z => mflambda_z(lvl)%dataptr(mfi)
    END IF

    chi_temp => mf_chi_temp%dataptr(mfi)
    zeta_x_temp => mf_zeta_x_temp%dataptr(mfi)
    zeta_y_temp => mf_zeta_y_temp%dataptr(mfi)
    pi_xx_temp => mf_pi_xx_temp%dataptr(mfi)
    pi_yy_temp => mf_pi_yy_temp%dataptr(mfi)
    pi_xy_temp => mf_pi_xy_temp%dataptr(mfi)
    lambda_x_temp => mf_lambda_x_temp%dataptr(mfi)
    lambda_y_temp => mf_lambda_y_temp%dataptr(mfi)
    IF (dimensions == 3) THEN
      zeta_z_temp => mf_zeta_z_temp%dataptr(mfi)
      pi_xz_temp => mf_pi_xz_temp%dataptr(mfi)
      pi_zz_temp => mf_pi_zz_temp%dataptr(mfi)
      pi_yz_temp => mf_pi_yz_temp%dataptr(mfi)
      lambda_z_temp => mf_lambda_z_temp%dataptr(mfi)
    END IF

    ! Since I can't do a logical mf, this will have to do
    node_pt => imf_node_done%dataptr(mfi)

    node_pt = 0
    !
    ! Loop through every node and direction Streaming fi values after  0.10675898826350119        8.1199848080719297E-002   7.5472603156772552E-002   7.1668243153084440E-002   6.2170421641303292E-002   7.3268019629481429E-002   7.1668294228699539E-002   3.5441048288652778E-002   2.7491027990083756E-002   2.9666618238734462E-002   3.5050647688395019E-002   3.5441048288652778E-002   2.7491027990083756E-002   2.9666744847538281E-002   3.5051319948574793E-002   2.5008288600439885E-002   2.0096507867446163E-002   2.1678354059375664E-002   1.7217542334600794E-002   2.1748385100336495E-002   2.1678648834210366E-002   4.1424512614056325E-003   4.5172756304959465E-003   3.4426717920899553E-003   2.8438668527215297E-003   4.2415902925789662E-003   3.8213442677570769E-003   4.5175267597818173E-003   3.8411752171805228E-003   3.8411522603449388E-003   3.4426717920899553E-003   4.2491662319996994E-003   4.2488426053475100E-003   2.8732290182369130E-003   2.2148247961806494E-003   2.9518222640950906E-003   2.3833651353159153E-003   2.6694021731949800E-003   3.1549611430556107E-003

    !
    do i = 1,dir
      do c = wax%lo(3)-reach,wax%hi(3)+reach
        do b = wax%lo(2)-reach,wax%hi(2)+reach
          do a = wax%lo(1)-reach,wax%hi(1)+reach
            !            write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
            !            if (state(a,b,c,1) >= 1000 .and. state(a,b,c,1) <51000 .and. lvl == 3) then
            !              write(*,*) 'Streaming fi values before',fi(a,b,c,1:dir)
            !              write(*,*) 'Streaming gi values before',gi(a,b,c,1:dir)
            !              write(*,*) 'values at streaming time',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),&
            !                temp(a,b,c,1),alpha(a,b,c,1)
            !              write(*,*) 'Streaming state before',state(a,b,c,1:dir),self
            !            end if
            !if (state(a,b,c,1) < 0 .and. state(a,b,c,1) >= -1000) cycle

            !select case (state(a,b,c,1))
              !case(0:100000000)
            if (state(a,b,c,1) >= 0) then
!              if (isnan(fi(a,b,c,i)) .or. isnan(fout(a,b,c,i)) .or. isnan(gi(a,b,c,i)) &
!                .or. isnan(gout(a,b,c,i)) .or. fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) <0.0D0 .or. &
!                fout(a,b,c,i) <0.0D0 .or. gout(a,b,c,i) < 0.0D0) then
!                write(*,*) 'shit crap fuck fuck fuck',a,b,c,i,state(a,b,c,1),state(a,b,c,i),&
!                  wax%lo,wax%hi,lvl,self
!                write(*,*) 'tigtig',state(a,b,c,1:dir)
!                write(*,*) 'lunlun',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
!
!                write(*,*) 'bunbun',fi(a,b,c,1:dir)
!                write(*,*) 'funfun',gi(a,b,c,1:dir)
!                write(*,*) 'munmun',fout(a,b,c,1:dir)
!                write(*,*) 'sunsun',gout(a,b,c,1:dir)
!
!                write(*,*) 'jubjub',chi(a,b,c,1),zeta_x(a,b,c,1),&
!                  zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                  pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!                  pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!              end if



!                            if (a == 268 .and. b == 22 .and. c == 176 .and. i == 39) then
!                              write(*,*) 'streaming states',state(a,b,c,1:dir),a,b,c,i,lvl
!                              write(*,*) 'streaming fis',fi(a,b,c,1:dir),a,b,c
!                              write(*,*) 'streaming gis',gi(a,b,c,1:dir),a,b,c
!                            end if
!              !!!
!                            if (a == 267 .and. b == 21 .and. c == 176 .and. i == 39) then
!                              write(*,*) 'streaming states',state(a,b,c,1:dir),a,b,c,&
!                                wax%lo,wax%hi,i,lvl
!                              write(*,*) 'streaming fis',fi(a,b,c,1:dir),a,b,c
!                              write(*,*) 'streaming gis',gi(a,b,c,1:dir),a,b,c
!                            end if
              !!
!              if (a == 139 .and. b == 127 .and. c == 120 .and. i == 39) then
!                write(*,*) 'streaming states',state(a,b,c,1:dir),state(a-1,b,c,1),&
!                  state(a-2,b,c,1),a,b,c,i
!                write(*,*) 'streaming fis',fi(a,b,c,1:dir)
!                write(*,*) 'streaming gis',gi(a,b,c,1:dir)
!                write(*,*) 'streaming fouts',fout(a,b,c,1:dir)
!                write(*,*) 'streaming gouts',gout(a,b,c,1:dir)
!              end if
!              if (a == 374 .and. b == 256 .and. c == 256 .and. i == 39) then
!                write(*,*) 'streaming states',state(a,b,c,1:dir),a,b,c,i
!                write(*,*) 'streaming fis',fi(a,b,c,1:dir)
!                write(*,*) 'streaming gis',gi(a,b,c,1:dir)
!                write(*,*) 'streaming fouts',fout(a,b,c,1:dir)
!                write(*,*) 'streaming gouts',gout(a,b,c,1:dir)
!
!              end if
!              if (a == 377 .and. b == 256 .and. c == 256 .and. i == 39) then
!                write(*,*) 'streaming states',state(a,b,c,1:dir),a,b,c,i
!                write(*,*) 'streaming fis',fi(a,b,c,1:dir)
!                write(*,*) 'streaming gis',gi(a,b,c,1:dir)
!                write(*,*) 'streaming fouts',fout(a,b,c,1:dir)
!                write(*,*) 'streaming gouts',gout(a,b,c,1:dir)
!              end if

              if (i == 1 .and. shifted .and. state(a,b,c,1) >= 0) then
                if (state(a+cx(i),b+cy(i),c+cz(i),1) >= 0) then
                  go to 1000
                else if (state(a+cx(i),b+cy(i),c+cz(i),1) == -1001 .or. &
                  state(a+cx(i),b+cy(i),c+cz(i),1) == -666) then
                  cycle
                !                else if
                !
                !                else if
                !
                !                else if
                !
                !                else if
                !
                !                else if

                else
                  cycle
                end if

              end if
              !
              ! To ensure pure outflow is handled correctly.
              !
              select case (state(a,b,c,i))
                case (0:1100000000)
1000              fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                  gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                !
                case(-12:-11)
                  !                call bc_switchboard(fout(a,b,c,i),fout(a,b,c,opp(i)),gout(a,b,c,i),&
                  !                  gout(a,b,c,opp(i)),i,state(a,b,c,i),j,state(a,b,c,1:dir),loc_cx,loc_cy,loc_cz)
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if

                  call inviscid_edge(fout(a,b,c,i),fout(a,b,c,opp(i)),gout(a,b,c,i),&
                    gout(a,b,c,opp(i)),i,state(a,b,c,i),state(a,b,c,1:dir),j,loc_cx,loc_cy,loc_cz)
                  !if (state(a+loc_cx,b+loc_cy,c+loc_cz,j) >= 0) then !XXXX figure out weird corner cases later
                  fi(a+loc_cx,b+loc_cy,c+loc_cz,j) = fout(a,b,c,i)
                  gi(a+loc_cx,b+loc_cy,c+loc_cz,j) = gout(a,b,c,i)
                  !else
!                  write(*,*) 'switchboard operator speaking!',state(a,b,c,1),state(a,b,c,i),&
!                    state(a,b,c,erp(i)),a,b,c,i

                  !end if

                case(-10:-1)
!                  write(*,*) 'switchboard operator speaking!',state(a,b,c,1),state(a,b,c,i),&
!                    state(a,b,c,erp(i)),a,b,c,i
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if
                  call bc_switchboard(fout(a,b,c,i),fout(a,b,c,opp(i)),gout(a,b,c,i),&
                    gout(a,b,c,opp(i)),i,state(a,b,c,i))
                  fi(a,b,c,opp(i)) = fout(a,b,c,i)
                  gi(a,b,c,opp(i)) = gout(a,b,c,i)
                  if (fi(a,b,c,opp(i)) < 0.0D0) then
                    write(*,*) 'negative fi in baseline bcs'
                  end if
                case(-199:-100)
!                  write(*,*) 'inlets',state(a,b,c,1),state(a,b,c,i),&
!                    state(a,b,c,erp(i)),a,b,c,i
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if
                  call inlet_bc(fout(a,b,c,i),gout(a,b,c,i),i,state(a,b,c,i))
                  fi(a,b,c,opp(i)) = fout(a,b,c,i)
                  gi(a,b,c,opp(i)) = gout(a,b,c,i)
                  if (fi(a,b,c,opp(i)) < 0.0D0) then
                    write(*,*) 'negative fi in inlets'
                  end if
                case(-299:-200)
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if
                  call outlet_bc(fout(a,b,c,i),fout(a,b,c,opp(i)),gout(a,b,c,i),gout(a,b,c,opp(i)),&
                    i,state(a,b,c,i),state(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),1),&
                    fi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),opp(i)),&
                    gi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),opp(i)),a,b,c)
                  fi(a,b,c,opp(i)) = fout(a,b,c,i)
                  gi(a,b,c,opp(i)) = gout(a,b,c,i)
                  if (fi(a,b,c,opp(i)) < 0.0D0 .or. gi(a,b,c,opp(i)) < 0.0D0) then
                    write(*,*) 'negative fi in outlet bc'
                    write(*,*) 'states for negative fi',state(a,b,c,1:dir),a,b,c,i
                    write(*,*) 'fis',fi(a,b,c,1:dir)
                    write(*,*) 'gis',gi(a,b,c,1:dir)
                  end if
                case(-399:-300)
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if
                  call inlet_bc(dummy_fout,dummy_gout,i,state(a,b,c,i))
                  fi(a,b,c,i) = fout(a,b,c,i)
                  gi(a,b,c,i) = gout(a,b,c,i)
                !                if (fi(a,b,c,opp(i)) < 0.0D0) then
                !                  write(*,*) 'negative fi in outlet bc'
                !                end if
                case(-499:-400)
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if
                  call outlet_bc(dummy_fout,fout(a,b,c,opp(i)),dummy_gout,gout(a,b,c,opp(i)),&
                    i,state(a,b,c,i),state(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),1),&
                    fi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),opp(i)),&
                    gi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),opp(i)),a,b,c)
                  fi(a,b,c,i) = fout(a,b,c,i)
                  gi(a,b,c,i) = gout(a,b,c,i)
                !                if (fi(a,b,c,opp(i)) < 0.0D0) then
                !                  write(*,*) 'negative fi in outlet bc'
                !                end if

                case(-666)
                  !                if (state(a,b,c,1) >= 200 .and. state(a,b,c,1) <= 299) then
                  !                  !write(*,*) 'roughest puff'
                  !                  fi(a,b,c,i) = fout(a,b,c,i)
                  !                  gi(a,b,c,i) = gout(a,b,c,i)
                  !                else
                  cycle
                  !end if
                case(-1000)
                  if (shifted) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)
                  end if
                  call freestream_edge_bc(fout(a,b,c,i),gout(a,b,c,i),i,state(a,b,c,i))
                  fi(a,b,c,opp(i)) = fout(a,b,c,i)
                  gi(a,b,c,opp(i)) = gout(a,b,c,i)
                  if (fi(a,b,c,opp(i)) < 0.0000001D0) then
                    write(*,*) 'negative fi in freestream bc'
                  end if
                case(-1100000000:-1001)

                  if (state(a,b,c,1) >= 1000 .and. state(a,b,c,1) < 51000) cycle

                  if (shifted .and. state(a+cx(i),b+cy(i),c+cz(i),1) >=0) then
                    fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                    gi(a+cx(i),b+cy(i),c+cz(i),i) = gout(a,b,c,i)

                  end if

                  if (node_pt(a,b,c,1) == 1) then !if the node already has target information
                    !
                    ! Find local fieq
                    !
                    if (any(state(a,b,c,1:15) < -1000)) then !If adjacent to a wall
                      if (.not. shifted) then
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),&
                              pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
                              lambda_z_temp(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),&
                              pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
                              lambda_z_temp(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                          case (39)
                            loc_fieq = fieq_comp_39(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq
                            loc_fieq_opp = fieq_comp_39(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                        end select
                      else
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),&
                              pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
                              lambda_z_temp(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),&
                              pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
                              lambda_z_temp(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                          case (39)
                            loc_fieq = fieq_comp_39_shifted(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq
                            loc_fieq_opp = fieq_comp_39_shifted(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                        end select
                      end if
                      !
                      ! Call the ibb subroutine.  It requires lots of input.
                      !
                      ! fout becomes the missing population, and is associated with loc_fieq
                      ! fi is a known population but is modified in startup_transient_damping, it's associated with loc_fieq_opp
                      !
                      ! These are passed to the IBB BC because during startup, the values can shift by a factor of 20 or more.
                      !   They are damped in order to prevent negative fi and gi values during collision.  The opposite direction
                      !   is modified to prevent divergence and non-physical behavior as startup is inherently non-physical itself.
                      !
                      do j = 1,7
                        if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= 0) then
                          u_local(j) = u_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          v_local(j) = v_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          w_local(j) = w_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          temp_local(j) = temp(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -11) then
                          u_local(j) = u_vel(a,b,c,1)
                          v_local(j) = v_vel(a,b,c,1)
                          w_local(j) = w_vel(a,b,c,1)
                          temp_local(j) = temp(a,b,c,1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -1000) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) <= -100 .and. &
                          state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= -499) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else
                          u_local(j) = 0.0D0
                          v_local(j) = 0.0D0
                          w_local(j) = 0.0D0
                          temp_local(j) = temperature
                        end if
                      end do
                      !                    if (any(u_local < very_small_number)) then
                      !
                      !                      write(*,*) 'bad u_local',state(a,b,c,1),state(a+1,b,c,1),state(a,b+1,c,1),state(a,b,c+1,1),&
                      !                        state(a-1,b,c,1),state(a,b-1,c,1),state(a,b,c-1,1),a,b,c
                      !
                      !                    end if
                      !                  if (any(state(a,b,c,1:15) < -1000)) then
                      !                    write(*,*) 'moving target thingy',rho_tgt_temp(a,b,c,1),u_tgt_temp(a,b,c,1),&
                      !                      v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1)

                      if (.not. shifted) then
                        call interpolative_bounceback_bc(fout(a,b,c,i),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& !
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_tgt_temp(a,b,c,1),v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1),&
                          opp(i),lvl,state(a,b,c,1:dir), rho_tgt_temp(a,b,c,1),&
                          u_local,v_local,w_local,temp_local,fine,a,b,c)
                      else

                        state_pass = (/state(a,b,c,1),state(a+1,b,c,1),state(a,b+1,c,1),&
                          state(a,b,c+1,1),state(a-1,b,c,1),state(a,b-1,c,1),state(a,b,c-1,1)/)

                        call shifted_ibb(fout(a,b,c,opp(i)),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& !
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_tgt_temp(a,b,c,1),v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1),&
                          i,lvl,state(a,b,c,1:dir), rho_tgt_temp(a,b,c,1),&
                          u_local,v_local,w_local,temp_local,state_pass,fine,a,b,c) !XXXX finish

                      end if
                    !
                    !
                    !
                    else !If not adjacent
                      if (.not. shifted) then
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                          case (39)
                            loc_fieq = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                        end select
                      else
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                          case (39)
                            loc_fieq = fieq_comp_39_shifted(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_39_shifted(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                        end select
                      end if
                      !                    u_local = (/u_vel(a,b,c,1),u_vel(a+1,b,c,1),u_vel(a,b+1,c,1),u_vel(a,b,c+1,1),u_vel(a-1,b,c,1),&
                      !                       u_vel(a,b-1,c,1),u_vel(a,b,c-1,1)/)
                      !                    v_local = (/v_vel(a,b,c,1),v_vel(a+1,b,c,1),v_vel(a,b+1,c,1),v_vel(a,b,c+1,1),v_vel(a-1,b,c,1),v_vel(a,b-1,c,1),&
                      !                      v_vel(a,b,c-1,1)/)
                      !                    w_local = (/w_vel(a,b,c,1),w_vel(a+1,b,c,1),w_vel(a,b+1,c,1),w_vel(a,b,c+1,1),w_vel(a-1,b,c,1),&
                      !                      w_vel(a,b-1,c,1),w_vel(a,b,c-1,1)/)
                      !                    temp_local = (/temp(a,b,c,1),temp(a+1,b,c,1),temp(a,b+1,c,1),temp(a,b,c+1,1),&
                      !                      temp(a-1,b,c,1),temp(a,b-1,c,1),temp(a,b,c-1,1)/)
                      !
                      !
                      !
                      do j = 1,7
                        if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= 0) then
                          u_local(j) = u_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          v_local(j) = v_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          w_local(j) = w_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          temp_local(j) = temp(a+cx(j),b+cy_27(j),c+cz_27(j),1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -11) then
                          u_local(j) = u_vel(a,b,c,1)
                          v_local(j) = v_vel(a,b,c,1)
                          w_local(j) = w_vel(a,b,c,1)
                          temp_local(j) = temp(a,b,c,1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -1000) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) <= -100 .and. &
                          state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= -499) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else
                          u_local(j) = 0.0D0
                          v_local(j) = 0.0D0
                          w_local(j) = 0.0D0
                          temp_local(j) = temperature
                        end if
                      end do


                      if (.not. shifted) then
                        call interpolative_bounceback_bc(fout(a,b,c,i),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& !
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
                          opp(i),lvl,state(a,b,c,1:dir), rho(a,b,c,1),u_local,v_local,w_local,temp_local,fine,a,b,c)

                      else

                        state_pass = (/state(a,b,c,1),state(a+1,b,c,1),state(a,b+1,c,1),&
                          state(a,b,c+1,1),state(a-1,b,c,1),state(a,b-1,c,1),state(a,b,c-1,1)/)

                        call shifted_ibb(fout(a,b,c,opp(i)),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& !
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
                          i,lvl,state(a,b,c,1:dir), rho(a,b,c,1),&
                          u_local,v_local,w_local,temp_local,state_pass,fine,a,b,c)
                      end if
                    end if
                  !
                  !                  fi(a,b,c,opp(i)) = fout(a,b,c,i)
                  !                  gi(a,b,c,opp(i)) = gout(a,b,c,i)


                    !write(*,*) 'IBB outputs ',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),i
                  !
                  ! Since the target values have not been found yet, find them, do the
                  ! bounceback, and make the target flag true
                  !
                  else !if the node does not already have target information
                    if (any(state(a,b,c,1:15) < -1000)) then !If adjacent to a wall, find target values

                      do j=1,dir
                        u_pass(j) = u_vel(a+cx(j),b+cy(j),c+cz(j),1)
                        v_pass(j) = v_vel(a+cx(j),b+cy(j),c+cz(j),1)
                        w_pass(j) = w_vel(a+cx(j),b+cy(j),c+cz(j),1)
                        temp_pass(j) = temp(a+cx(j),b+cy(j),c+cz(j),1)

                        !write(*,*) 'passing information',u_pass(j),v_pass(j),temp_pass(j)
                      end do
                      !
                      ! Find the velocity targets and temperature target the first time.
                      !
                      !                  write(*,*) 'temp mf locs',loc(u_tgt_temp(a,b,c,1)),loc(v_tgt_temp(a,b,c,1)),loc(node_pt(a,b,c,1))
                      !                  write(*,*) 'multifabs from modules',loc(chi(a,b,c,1)),loc(rho(a,b,c,1)),loc(fi(a,b,c,1))
                      call target_finder(state(a,b,c,1:dir),u_pass,v_pass,&
                        w_pass,temp_pass,fi(a,b,c,1:dir),rho_tgt_temp(a,b,c,1),&
                        u_tgt_temp(a,b,c,1),v_tgt_temp(a,b,c,1),&
                        w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1),a,b,c)

                      chi_temp(a,b,c,1) = chi(a,b,c,1)
                      zeta_x_temp(a,b,c,1) = zeta_x(a,b,c,1)
                      zeta_y_temp(a,b,c,1) = zeta_y(a,b,c,1)
                      zeta_z_temp(a,b,c,1) = zeta_z(a,b,c,1)
                      pi_xx_temp(a,b,c,1) = pi_xx(a,b,c,1)
                      pi_yy_temp(a,b,c,1) = pi_yy(a,b,c,1)
                      pi_zz_temp(a,b,c,1) = pi_zz(a,b,c,1)
                      pi_xy_temp(a,b,c,1) = pi_xy(a,b,c,1)
                      pi_xz_temp(a,b,c,1) = pi_xz(a,b,c,1)
                      pi_yz_temp(a,b,c,1) = pi_yz(a,b,c,1)
                      lambda_x_temp(a,b,c,1) = lambda_x(a,b,c,1)
                      lambda_y_temp(a,b,c,1) = lambda_y(a,b,c,1)
                      lambda_z_temp(a,b,c,1) = lambda_z(a,b,c,1)
                      !
                      !
                      !
                      pxx_temp = rho_tgt_temp(a,b,c,1)*temp_tgt_temp(a,b,c,1)+rho_tgt_temp(a,b,c,1)*u_tgt_temp(a,b,c,1)**2
                      pyy_temp = rho_tgt_temp(a,b,c,1)*temp_tgt_temp(a,b,c,1)+rho_tgt_temp(a,b,c,1)*v_tgt_temp(a,b,c,1)**2
                      pzz_temp = rho_tgt_temp(a,b,c,1)*temp_tgt_temp(a,b,c,1)+rho_tgt_temp(a,b,c,1)*w_tgt_temp(a,b,c,1)**2
                      pxy_temp = rho_tgt_temp(a,b,c,1)*u_tgt_temp(a,b,c,1)*v_tgt_temp(a,b,c,1)
                      pxz_temp = rho_tgt_temp(a,b,c,1)*u_tgt_temp(a,b,c,1)*w_tgt_temp(a,b,c,1)
                      pyz_temp = rho_tgt_temp(a,b,c,1)*v_tgt_temp(a,b,c,1)*w_tgt_temp(a,b,c,1)
                      !
                      !
                      !
                      if (timestep(lvl) <= 4 .and. shifted) then

                        chi_temp(a,b,c,1) = su_chi
                        zeta_x_temp(a,b,c,1) = su_zeta_x
                        zeta_y_temp(a,b,c,1) = su_zeta_y
                        zeta_z_temp(a,b,c,1) = su_zeta_z
                        pi_xx_temp(a,b,c,1) = su_pixx
                        pi_yy_temp(a,b,c,1) = su_piyy
                        pi_zz_temp(a,b,c,1) = su_pizz
                        pi_xy_temp(a,b,c,1) = su_pixy
                        pi_xz_temp(a,b,c,1) = su_pixz
                        pi_yz_temp(a,b,c,1) = su_piyz
                        lambda_x_temp(a,b,c,1) = su_lambda_x
                        lambda_y_temp(a,b,c,1) = su_lambda_y
                        lambda_z_temp(a,b,c,1) = su_lambda_z

                        call calculate_lagrangian_mults(rho_tgt_temp(a,b,c,1),u_tgt_temp(a,b,c,1),&
                          v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),&
                          temp_tgt_temp(a,b,c,1),pxx_temp,pyy_temp,pzz_temp,pxy_temp,pxz_temp,pyz_temp,&
                          chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                          pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),pi_xy_temp(a,b,c,1),&
                          pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
                          lambda_z_temp(a,b,c,1),pass_fieq,a,b,c)



                      else
                        call calculate_lagrangian_mults(rho_tgt_temp(a,b,c,1),u_tgt_temp(a,b,c,1),&
                          v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),&
                          temp_tgt_temp(a,b,c,1),pxx_temp,pyy_temp,pzz_temp,pxy_temp,pxz_temp,pyz_temp,&
                          chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                          pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),pi_xy_temp(a,b,c,1),&
                          pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
                          lambda_z_temp(a,b,c,1),pass_fieq,a,b,c)
                      end if
                      !                  loc_fieq = pass_fieq(i)
                      !                  loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq !XXXX change back to loc_fieq
                      !
                      !select case (dir)
                      if (.not. shifted) then
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq
                            loc_fieq_opp = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
                              pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp

                          case (39)
                            loc_fieq = fieq_comp_39(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_39(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),i)
                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                            !write(*,*) 'local equilibrium values',loc_fieq,loc_gieq
                        end select
                      else
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq
                            loc_fieq_opp = fieq_comp_19(rho_tgt_temp(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
                              pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp

                          case (39)
                            loc_fieq = fieq_comp_39_shifted(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_39_shifted(rho_tgt_temp(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
                              zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
                              pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
                              pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),&
                              lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),i)
                            loc_gieq_opp = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq_opp
                            !write(*,*) 'local equilibrium values',loc_fieq,loc_gieq
                        end select
                      end if

                      !                  if (loc_fieq < 0.0D0 .or. loc_fieq > 1.0D0 .or. loc_gieq < 0.0D0 .or. &
                      !                      loc_gieq > 1.0D0) then
                      !                    write(*,*) 'Bad local equilibrium values',loc_fieq,loc_gieq,&
                      !                        rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                      !                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                      !                        pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
                      !                        pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz_temp(a,b,c,1),&
                      !                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
                      !                  end if
                      do j = 1,7
                        if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= 0) then
                          u_local(j) = u_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          v_local(j) = v_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          w_local(j) = w_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          temp_local(j) = temp(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -11) then
                          u_local(j) = u_vel(a,b,c,1)
                          v_local(j) = v_vel(a,b,c,1)
                          w_local(j) = w_vel(a,b,c,1)
                          temp_local(j) = temp(a,b,c,1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -1000) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) <= -100 .and. &
                          state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= -499) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else
                          u_local(j) = 0.0D0
                          v_local(j) = 0.0D0
                          w_local(j) = 0.0D0
                          temp_local(j) = temperature
                        end if
                      end do
                      !                    u_local = (/u_vel(a,b,c,1),u_vel(a+1,b,c,1),u_vel(a,b+1,c,1),u_vel(a,b,c+1,1),u_vel(a-1,b,c,1),&
                      !                       u_vel(a,b-1,c,1),u_vel(a,b,c-1,1)/)
                      !                    v_local = (/v_vel(a,b,c,1),v_vel(a+1,b,c,1),v_vel(a,b+1,c,1),v_vel(a,b,c+1,1),v_vel(a-1,b,c,1),v_vel(a,b-1,c,1),&
                      !                      v_vel(a,b,c-1,1)/)
                      !                    w_local = (/w_vel(a,b,c,1),w_vel(a+1,b,c,1),w_vel(a,b+1,c,1),w_vel(a,b,c+1,1),w_vel(a-1,b,c,1),&
                      !                      w_vel(a,b-1,c,1),w_vel(a,b,c-1,1)/)
                      !                    temp_local = (/temp(a,b,c,1),temp(a+1,b,c,1),temp(a,b+1,c,1),temp(a,b,c+1,1),&
                      !                      temp(a-1,b,c,1),temp(a,b-1,c,1),temp(a,b,c-1,1)/)
                      !
                      if (.not. shifted) then
                        call interpolative_bounceback_bc(fout(a,b,c,i),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& ! XXXX change back to loc_fieq
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_tgt_temp(a,b,c,1),v_tgt_temp(a,b,c,1),&
                          w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1),opp(i),lvl,state(a,b,c,1:dir),&
                          rho_tgt_temp(a,b,c,1),u_local,v_local,w_local,temp_local,fine,a,b,c)
                      else
                        state_pass = (/state(a,b,c,1),state(a+1,b,c,1),state(a,b+1,c,1),&
                          state(a,b,c+1,1),state(a-1,b,c,1),state(a,b-1,c,1),state(a,b,c-1,1)/)

                        call shifted_ibb(fout(a,b,c,opp(i)),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& !
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_tgt_temp(a,b,c,1),v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1),&
                          i,lvl,state(a,b,c,1:dir), rho_tgt_temp(a,b,c,1),&
                          u_local,v_local,w_local,temp_local,state_pass,fine,a,b,c)
                      end if
                    !
                    !
                    !
                    else !If not adjacent to a wall, skip the target finding process
                      if (.not. shifted) then
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                          case (39)
                            loc_fieq = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                        end select
                      else
                        select case (dir)
                          case (19)
                            loc_fieq = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))

                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                          case (39)
                            loc_fieq = fieq_comp_39_shifted(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),opp(i))
                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq

                            loc_fieq_opp = fieq_comp_39_shifted(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
                              zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
                              pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
                              pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
                            loc_gieq_opp = k_gieq*temp(a,b,c,1)*loc_fieq_opp
                        end select
                      end if
                      !
                      !
                      !
                      do j = 1,7
                        if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= 0) then
                          u_local(j) = u_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          v_local(j) = v_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          w_local(j) = w_vel(a+cx_27(j),b+cy_27(j),c+cz_27(j),1)
                          temp_local(j) = temp(a+cx(j),b+cy(j),c+cz(j),1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -11) then
                          u_local(j) = u_vel(a,b,c,1)
                          v_local(j) = v_vel(a,b,c,1)
                          w_local(j) = w_vel(a,b,c,1)
                          temp_local(j) = temp(a,b,c,1)
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) == -1000) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else if (state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) <= -100 .and. &
                          state(a+cx_27(j),b+cy_27(j),c+cz_27(j),1) >= -499) then
                          u_local(j) = fs_u
                          v_local(j) = fs_v
                          w_local(j) = fs_w
                          temp_local(j) = temperature
                        else
                          u_local(j) = 0.0D0
                          v_local(j) = 0.0D0
                          w_local(j) = 0.0D0
                          temp_local(j) = temperature
                        end if
                      end do
                      !                    u_local = (/u_vel(a,b,c,1),u_vel(a+1,b,c,1),u_vel(a,b+1,c,1),u_vel(a,b,c+1,1),u_vel(a-1,b,c,1),&
                      !                       u_vel(a,b-1,c,1),u_vel(a,b,c-1,1)/)
                      !                    v_local = (/v_vel(a,b,c,1),v_vel(a+1,b,c,1),v_vel(a,b+1,c,1),v_vel(a,b,c+1,1),v_vel(a-1,b,c,1),v_vel(a,b-1,c,1),&
                      !                      v_vel(a,b,c-1,1)/)
                      !                    w_local = (/w_vel(a,b,c,1),w_vel(a+1,b,c,1),w_vel(a,b+1,c,1),w_vel(a,b,c+1,1),w_vel(a-1,b,c,1),&
                      !                      w_vel(a,b-1,c,1),w_vel(a,b,c-1,1)/)
                      !                    temp_local = (/temp(a,b,c,1),temp(a+1,b,c,1),temp(a,b+1,c,1),temp(a,b,c+1,1),&
                      !                      temp(a-1,b,c,1),temp(a,b-1,c,1),temp(a,b,c-1,1)/)
                      !
                      !
                      !
                      if (.not. shifted) then
                        call interpolative_bounceback_bc(fout(a,b,c,i),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& ! XXXX change back to loc_fieq
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_vel(a,b,c,1),v_vel(a,b,c,1),&
                          w_vel(a,b,c,1),temp(a,b,c,1),opp(i),lvl,state(a,b,c,1:dir),&
                          rho(a,b,c,1),u_local,v_local,w_local,temp_local,fine,a,b,c)
                      else
                        state_pass = (/state(a,b,c,1),state(a+1,b,c,1),state(a,b+1,c,1),&
                          state(a,b,c+1,1),state(a-1,b,c,1),state(a,b-1,c,1),state(a,b,c-1,1)/)

                        call shifted_ibb(fout(a,b,c,opp(i)),fi(a,b,c,i),fi(a,b,c,opp(i)),loc_fieq,loc_fieq_opp,& !
                          fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),gout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),&
                          gout(a,b,c,i),gi(a,b,c,i),gi(a,b,c,opp(i)),loc_gieq,loc_gieq_opp,&
                          u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
                          i,lvl,state(a,b,c,1:dir), rho(a,b,c,1),&
                          u_local,v_local,w_local,temp_local,state_pass,fine,a,b,c)
                      end if
                    end if

!                  if (abs(fi(a,b,c,opp(i)) - fout(a,b,c,i)) > 0.01D0 .or. &
!                      abs(gi(a,b,c,opp(i))-gout(a,b,c,i)) > 0.06D0) then
!                    write(*,*) 'what the fluff',fi(a,b,c,opp(i)),fi(a,b,c,i),fout(a,b,c,i),loc_fieq,gi(a,b,c,opp(i)),&
!                      gi(a,b,c,i),gout(a,b,c,i),loc_gieq,u_tgt_temp(a,b,c,1),i
!                  end if
!                  if (isnan(fout(a,b,c,i)) .or. isnan(fout(a,b,c,i))) then
!                    write(*,*) 'its a nan nan nan world',fout(a,b,c,i),gout(a,b,c,i),loc_fieq,loc_gieq,u_tgt_temp(a,b,c,1),&
!                      v_tgt_temp(a,b,c,1),w_tgt_temp(a,b,c,1),temp_tgt_temp(a,b,c,1)
!                    write(*,*) 'whats up with LMs',chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),zeta_y_temp(a,b,c,1),&
!                    zeta_z_temp(a,b,c,1),&
!                    pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),pi_xy_temp(a,b,c,1),&
!                    pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),&
!                    lambda_z_temp(a,b,c,1)
!                  end if

                    node_pt(a,b,c,1) = 1

                  end if !targer conditional
                  fi(a,b,c,opp(i)) = fout(a,b,c,i)
                  gi(a,b,c,opp(i)) = gout(a,b,c,i)
              end select

            else
              cycle
            end if !if the node is a fluid node

          if (fi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) < 0.0D0) then
            write(*,*) 'stream information, negative fi or gi',state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i)),&
              state(a,b,c,1:dir),a,b,c,i,lvl
            write(*,*) 'stream values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
              alpha(a,b,c,1),a,b,c
            write(*,*) 'stream fis',fi(a,b,c,1:dir),a,b,c
            write(*,*) 'stream fouts',fout(a,b,c,1:dir),a,b,c
            write(*,*) 'stream gis',gi(a,b,c,1:dir),a,b,c
            write(*,*) 'stream gouts',gout(a,b,c,1:dir),a,b,c

          end if

          !            if (a == 67 .and. b == 37 .and. c == 32 .and. i == 39) then
          !              write(*,*) 'Streaming fi values after',fi(a,b,c,1:dir)
          !              write(*,*) 'Streaming gi values after',gi(a,b,c,1:dir)
          !              write(*,*) 'Streaming states',state(a,b,c,1:dir),self
          !            end if

          end do !x loop
        end do !y loop
      end do !z loop
    end do !direction loop
    !write(*,*) 'Next box please!',self,lvl
  end do



  call amrex_mfiter_destroy(mfi)

  call amrex_multifab_destroy(mfrho_tgt_temp)
  call amrex_multifab_destroy(mfu_tgt_temp)
  call amrex_multifab_destroy(mfv_tgt_temp)
  call amrex_multifab_destroy(mfw_tgt_temp)
  call amrex_multifab_destroy(mftemp_tgt_temp)

  call amrex_multifab_destroy(mf_chi_temp)
  call amrex_multifab_destroy(mf_zeta_x_temp)
  call amrex_multifab_destroy(mf_zeta_y_temp)
  call amrex_multifab_destroy(mf_zeta_z_temp)

  call amrex_multifab_destroy(mf_pi_xx_temp)
  call amrex_multifab_destroy(mf_pi_yy_temp)
  call amrex_multifab_destroy(mf_pi_zz_temp)
  call amrex_multifab_destroy(mf_pi_xy_temp)
  call amrex_multifab_destroy(mf_pi_xz_temp)
  call amrex_multifab_destroy(mf_pi_yz_temp)

  call amrex_multifab_destroy(mf_lambda_x_temp)
  call amrex_multifab_destroy(mf_lambda_y_temp)
  call amrex_multifab_destroy(mf_lambda_z_temp)

  call amrex_imultifab_destroy(imf_node_done)

  call mpi_barrier(commune,ier)

  !  write(*,*) 'Done streaming on level',lvl
  !end do

  deallocate(u_pass)
  deallocate(v_pass)
  deallocate(w_pass)
  deallocate(temp_pass)

!write(*,*) 'Done streaming'

END SUBROUTINE
!            if (fi(a,b,c,opp(i)) < 0.0D0 .and. i>1 .and. state(a,b,c,1) >=0) then
!              write(*,*) 'fi opp after streaming',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),fout(a,b,c,i),gout(a,b,c,i),&
!                state(a,b,c,1),state(a,b,c,opp(i)),&
!                rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl,timestep(lvl)
!              write(*,*) 'relevant values',state(a,b,c,1),state(a,b,c,i),state(a-cx(i),b-cy(i),c-cz(i),1),&
!                state(a-cx(i),b-cy(i),c-cz(i),i)
!            end if
!!            if (fi(a,b,c,i) < 0.0D0 .and. state(a,b,c,1) >=0 )then
!!              write(*,*) 'fi after streaming',fi(a,b,c,i),gi(a,b,c,i),state(a,b,c,1),state(a,b,c,i),&
!!                rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl,timestep(lvl)
!!              write(*,*) 'locations and bounds',a,b,c,wax%lo,wax%hi
!!              write(*,*) 'weird fis',fi(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i),fout(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i)
!!              write(*,*) 'relevant values',state(a,b,c,1),state(a,b,c,i),state(a,b,c,opp(i)),&
!!                state(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),1),&
!!                state(a+cx(opp(i)),b+cy(opp(i)),c+cz(opp(i)),i)
!!            end if
!            if (fi(a,b,c,i) < 0.0D0 .or. fi(a,b,c,i) > 1.0D0 .or. fout(a,b,c,1) < 0.0D0 .or. &
!                 fout(a,b,c,1) > 1.0D0) then
!              write(*,*) 'fi bad after streaming',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),fout(a,b,c,i),gout(a,b,c,i),&
!                state(a,b,c,1),state(a,b,c,opp(i)),&
!                rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl,timestep(lvl)
!              write(*,*) 'relevant values',state(a,b,c,1),state(a,b,c,i),state(a-cx(i),b-cy(i),c-cz(i),1),&
!                state(a-cx(i),b-cy(i),c-cz(i),i)
!            end if
!            if (gi(a,b,c,i) < 0.0D0 .or. gi(a,b,c,i) > 1.0D0 .or. gout(a,b,c,1) < 0.0D0 .or. &
!                 gout(a,b,c,1) > 1.0D0) then
!              write(*,*) 'gi bad after streaming',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),fout(a,b,c,i),gout(a,b,c,i),&
!                state(a,b,c,1),state(a,b,c,opp(i)),&
!                rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl,timestep(lvl)
!              write(*,*) 'relevant values',state(a,b,c,1),state(a,b,c,i),state(a-cx(i),b-cy(i),c-cz(i),1),&
!                state(a-cx(i),b-cy(i),c-cz(i),i)
!
!            end if

!                    case (39)
!                      loc_fieq = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                        pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
!                        pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz_temp(a,b,c,1),&
!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)

!            if (state(a,b,c,1) < -300) cycle
!            if (state(a+cx(i),b+cy(i),c+cz(i),1) >=0) then
!              fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
!            else if (state(a+cx(i),b+cy(i),c+cz(i),1) <=-1001) then
!              call interpolative_bounceback_bc()
!            else if (state(a+cx(i),b+cy(i),c+cz(i),1) ==-1000) then
!              call freestream_bc()
!            else if (state(a+cx(i),b+cy(i),c+cz(i),1) <=-200) then
!              call
!            else if (state(a+cx(i),b+cy(i),c+cz(i),1) <=-100) then
!              call
!            else if (state(a+cx(i),b+cy(i),c+cz(i),1) ) then
!              call
!            end if
!
! Find the loc fi_eq value for the given direction.
!
!                      select case (dir)
!                        case (19)
!                          loc_fieq = fieq_comp_19(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!                            zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                            pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),pi_xy(a,b,c,1),pi_xz(a,b,c,1),&
!                            pi_yz(a,b,c,1),lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
!
!                            loc_gieq = k_gieq*temp(a,b,c,1)*loc_fieq
!                        case (39)
!                          loc_fieq = fieq_comp_39(rho(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
!                            zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
!                            pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),pi_xy_temp(a,b,c,1),&
!                            pi_xz_temp(a,b,c,1),pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),&
!                            lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),i)
!
!                            loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq
!
!                            !write(*,*) 'local equilibrium values',loc_fieq,loc_gieq
!                      end select
!
! Call the ibb subroutine.  It requires lots of input
!
!                  if (loc_fieq > 0.08D0 .and. i > 7) then
!                    write(*,*) 'fieq calculation values',loc_fieq,&
!                        rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
!                        chi(a,b,c,1),zeta_x(a,b,c,1),&
!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                        pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
!                        pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz_temp(a,b,c,1),&
!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!                    !write(*,*)
!                  end if
!                  if (loc_fieq > 0.08D0 .and. i > 7) then
!                    write(*,*) 'fieq calculation values',loc_fieq,&
!                        rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),&
!                        chi(a,b,c,1),zeta_x(a,b,c,1),&
!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                        pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
!                        pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz_temp(a,b,c,1),&
!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1)
!                    !write(*,*)
!                  end if
!
!            if (a == 51 .and. b == 37 .and. c == 2 .and. lvl == 1) then
!              write(*,*) 'ARGH!!! in streaming',state(a,b,c,1),state(a,b,c,2),state(a,b,c,3),state(a,b,c,4),&
!                    state(a,b,c,5),state(a,b,c,6),state(a,b,c,7),i,wax%lo,wax%hi,self
!            end if


!                  if (fi(a,b,c,opp(i)) < 0.0D0) then
!                    write(*,*) 'fi check after ibb',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),state(a,b,c,opp(i)),&
!                      rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl,timestep(lvl)
!                  end if
!                  if (gi(a,b,c,opp(i)) < 0.0D0) then
!                    write(*,*) 'gi check after ibb',fi(a,b,c,opp(i)),gi(a,b,c,opp(i)),state(a,b,c,opp(i)),&
!                      rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1),a,b,c,i,lvl
!                  end if
                  !write(*,*) 'IBB outputs',fi(a,b,c,opp(i)),gi(a,b,c,opp(i))


!                  if (abs(fi(a,b,c,opp(i)) - fout(a,b,c,i)) > 0.03D0 .or. &
!                      abs(gi(a,b,c,opp(i))-gout(a,b,c,i)) > 0.06D0) then
!                    write(*,*) 'what the fluff',fi(a,b,c,opp(i)),fi(a,b,c,i),fout(a,b,c,i),gi(a,b,c,opp(i)),&
!                      gi(a,b,c,i),gout(a,b,c,i),i
!                  end if


!              if (fout(a,b,c,i) < 0.0D0 .or. fi(a,b,c,i) <0.0D0 .or. fout(a,b,c,i) > 1.0D0 .or. &
!                  fi(a,b,c,i) > 1.0D0 .or. temp(a,b,c,1) > 1.0D0) then
!                write(*,*) 'fi check on streaming',fi(a,b,c,i),fout(a,b,c,i),gi(a,b,c,i),gout(a,b,c,i),&
!                  state(a,b,c,1),state(a,b,c,i),a,b,c,i,lvl
!                write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
!                write(*,*) 'node check for bad stuff',state(a,b,c,1:dir),a,b,c,i,lvl,timestep(lvl),wax%lo,wax%hi
!              end if
!              if (gout(a,b,c,i) < 0.00 .or. gi(a,b,c,i) <0.0D0 .or. gout(a,b,c,i) > 1.0D0 .or. &
!                  gi(a,b,c,i) > 1.0D0  .or. temp(a,b,c,1) > 1.0D0) then
!                write(*,*) 'gi check on streaming',fout(a,b,c,i),gout(a,b,c,i),&
!                  state(a,b,c,1),state(a,b,c,i),a,b,c,i,lvl,wax%lo,wax%hi
!                write(*,*) 'local values',rho(a,b,c,1),u_vel(a,b,c,1),v_vel(a,b,c,1),w_vel(a,b,c,1),temp(a,b,c,1)
!                write(*,*) 'node check for bad stuff',state(a,b,c,1:dir),a,b,c,i,lvl,timestep(lvl)
!              end if




                !else
                  !fi(a+cx(i),b+cy(i),c+cz(i),i) = fout(a,b,c,i)
                !end if
!              case(-99:-1)
!                call bc_switchboard(fout(a,b,c,i),fout(a,b,c,opp(i)),gout(a,b,c,i),&
!                  gout(a,b,c,opp(i)),i,state(a,b,c,i))
!                fi(a,b,c,opp(i)) = fout(a,b,c,i)
!                gi(a,b,c,opp(i)) = gout(a,b,c,i)
!              case(-199:-100)
!                call inlet_bc(fout(a,b,c,i),gout(a,b,c,i),i,state(a,b,c,i))
!                fi(a,b,c,opp(i)) = fout(a,b,c,i)
!                gi(a,b,c,opp(i)) = gout(a,b,c,i)
!              case(-299:-200)
!                call outlet_bc(fout(a,b,c,i),gout(a,b,c,i),i,state(a,b,c,i))
!               fi(a,b,c,opp(i)) = fout(a,b,c,i)
!                gi(a,b,c,opp(i)) = gout(a,b,c,i)
!              case(-666)
!                cycle
!              case(-1000)
!                call freestream_bc(fout(a,b,c,i),gout(a,b,c,i),i,state(a,b,c,i))
!                fi(a,b,c,opp(i)) = fout(a,b,c,i)
!                gi(a,b,c,opp(i)) = gout(a,b,c,i)
!              case(-1000000000:-1001)
!                !cycle
!              end select

!
!              if (state(a,b,c,1) >= 200 .and. state(a,b,c,1) <= 299 .and. &
!                  state(a+cx(i),b+cy(i),c+cz(i),i) >= 0 ) then
!                fi(a,b,c,1:dir) = fout(a,b,c,1:dir)
!                gi(a,b,c,1:dir) = gout(a,b,c,1:dir)
!              end if
!                      loc_fieq = fieq_comp_19(rho(a,b,c,1),chi_temp(a,b,c,1),zeta_x_temp(a,b,c,1),&
!                        zeta_y_temp(a,b,c,1),zeta_z_temp(a,b,c,1),&
!                        pi_xx_temp(a,b,c,1),pi_yy_temp(a,b,c,1),pi_zz_temp(a,b,c,1),&
!                        pi_xy_temp(a,b,c,1),pi_xz_temp(a,b,c,1),&
!                        pi_yz_temp(a,b,c,1),lambda_x_temp(a,b,c,1),lambda_y_temp(a,b,c,1),lambda_z_temp(a,b,c,1),i)
!
!                      loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq
!                      loc_fieq = fieq_comp_39(rho(a,b,c,1),chi(a,b,c,1),zeta_x(a,b,c,1),&
!                        zeta_y(a,b,c,1),zeta_z(a,b,c,1),&
!                        pi_xx(a,b,c,1),pi_yy(a,b,c,1),pi_zz(a,b,c,1),&
!                        pi_xy(a,b,c,1),pi_xz(a,b,c,1),pi_yz(a,b,c,1),&
!                        lambda_x(a,b,c,1),lambda_y(a,b,c,1),lambda_z(a,b,c,1),i)
!
!                      loc_gieq = k_gieq*temp_tgt_temp(a,b,c,1)*loc_fieq









!                    u_local = (/u_vel(a,b,c,1),u_vel(a+1,b,c,1),u_vel(a,b+1,c,1),u_vel(a,b,c+1,1),u_vel(a-1,b,c,1),&
!                       u_vel(a,b-1,c,1),u_vel(a,b,c-1,1)/)
!                    v_local = (/v_vel(a,b,c,1),v_vel(a+1,b,c,1),v_vel(a,b+1,c,1),v_vel(a,b,c+1,1),v_vel(a-1,b,c,1),v_vel(a,b-1,c,1),&
!                      v_vel(a,b,c-1,1)/)
!                    w_local = (/w_vel(a,b,c,1),w_vel(a+1,b,c,1),w_vel(a,b+1,c,1),w_vel(a,b,c+1,1),w_vel(a-1,b,c,1),&
!                      w_vel(a,b-1,c,1),w_vel(a,b,c-1,1)/)
!                    temp_local = (/temp(a,b,c,1),temp(a+1,b,c,1),temp(a,b+1,c,1),temp(a,b,c+1,1),&
!                      temp(a-1,b,c,1),temp(a,b-1,c,1),temp(a,b,c-1,1)/)
