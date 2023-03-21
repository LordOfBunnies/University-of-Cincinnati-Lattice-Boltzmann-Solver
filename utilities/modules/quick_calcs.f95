module quick_calcs
!
! Putting these in an external procedure for error check and optional arguments
!   to be able to work.
!
!
!
!
!
use precise
use constants
implicit none
save
private
public :: quick_sutherlands,quick_therm_cond



contains
!
!
!
!
SUBROUTINE quick_sutherlands(loc_mu,loc_temp,inter_lvl,dim_mu,dim_temp)
!
! A quick Sutherland's Law calculation to find the viscosity at the local node
!
!
! Called by: collision
! Calls:
!
IMPLICIT NONE
REAL(kind=dp),INTENT(IN) :: loc_temp
REAL(kind=dp),INTENT(OUT) :: loc_mu
real(kind=dp),optional :: dim_mu,dim_temp
REAL(kind = dp),parameter :: const_s = 110.15,temp_ref = 273.15,visc_ref = 1.716D-5
REAL(kind=dp) :: unit_mu,unit_temp
logical,optional :: inter_lvl



unit_temp = loc_temp*t_ref*gama

unit_mu = visc_ref*((unit_temp/temp_ref)**1.5)*((temp_ref+&
  const_s)/(unit_temp+const_s))

if (present(inter_lvl)) then
  if (inter_lvl) then
    loc_mu = (unit_mu*2)/mu_ref
  else
    loc_mu = unit_mu/mu_ref
  end if
else
  loc_mu = unit_mu/mu_ref
end if

if (present(dim_mu)) then
  dim_mu = unit_mu
end if
if (present(dim_temp)) then
  dim_temp = unit_temp
end if

END SUBROUTINE
!
!
!
!
!
SUBROUTINE quick_therm_cond(loc_kappa,loc_rho,loc_temp,startup)
!
! This calculates the heat conductivity for a broad range of conditions
!
!
! Called by: collision,basic_fluids_calc
! Calls: error_out
!
IMPLICIT NONE
REAL(kind=dp),INTENT(IN) :: loc_rho,loc_temp
REAL(kind=dp),INTENT(OUT) :: loc_kappa
REAL(kind=dp) :: unit_temp,unit_kappa,unit_rho,unit_p
REAL(kind=dp) :: ak,bk,ck,dk,ek,kek,kappa_bottom
logical,optional :: startup
!INTEGER ::

!mfp =
!c_bar = 3.0D0*kb*temper/molec_mass_air

if (present(startup)) then
  kappa_bottom = 1.0D0
else
  kappa_bottom = kappa_ref
end if

unit_temp = loc_temp*t_ref*gama
unit_rho = loc_rho*rho_ref
unit_p = unit_temp*unit_rho*gas_const_R
!write(*,*) 'temps',unit_temp,loc_temp,t_ref,gama

if (unit_temp <= 500.0D0) then
  loc_kappa = kappa_convert*(5.9776D-6 * unit_temp**1.5 / (unit_temp +194.4))/kappa_bottom
else if (unit_temp > 500.0D0 .AND. unit_temp <= 2250.0D0) THEN
  unit_rho = loc_rho*rho_ref

  unit_p = unit_rho*gas_const_R*unit_temp/101300.0D0
!
! Set up the coefficients for the most likely band
!
  if (unit_p < 10.0D0 .AND. unit_p > 0.1D0) THEN
    ak = 0.334316D0
    bk = 3.28202D0
    ck = 11.9939D0
    dk = 20.0944D0
    ek = 4.62882D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else if (unit_p > 10.0D0 .AND. unit_p < 100.0D0) then
    ak = 0.413573D0
    bk = 3.83393D0
    ck = 13.1885D0
    dk = 20.7305D0
    ek = 4.27728D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else if (unit_p > 100.0D0 .AND. unit_p < 1000.0D0) then
    ak = 0.208749D0
    bk = 1.92122D0
    ck = 6.58813D0
    dk = 10.7630D0
    ek = -1.27699D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else if (unit_p > 0.10D0 .AND. unit_p < 0.010D0) then
    ak = 1.05928D0
    bk = 10.0924D0
    ck = 35.6709D0
    dk = 56.1818D0
    ek = 24.9670D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else
    WRITE(*,*) 'Parameters are way out of whack, something went wrong'
    write(11,*) 'Parameters are way out of whack, something went wrong'
    call error_out
  end if
!
! Find the local value of kappa, the thermal conductivity
!

else if (unit_temp > 2250.0D0 .AND. unit_temp <= 4250.0D0) THEN

  if (unit_p < 10.0D0 .AND. unit_p > 0.1D0) THEN
    ak = 10.9992D0
    bk = 387.106D0
    ck = 387.282D0
    dk = 5.48304D0
    ek = -1.20106D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else if (unit_p > 10.0D0 .AND. unit_p < 100.0D0) then
    ak = 82.1184D0
    bk = 308.927D0
    ck = 423.174D0
    dk = 250.668D0
    ek = 47.5889D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else if (unit_p > 100.0D0 .AND. unit_p < 1000.0D0) then
    ak = 38.8677D0
    bk = 123.284D0
    ck = 144.224D0
    dk = 72.8083D0
    ek = 0.684807D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else if (unit_p > 0.10D0 .AND. unit_p < 0.010D0) then
    ak = 101.351D0
    bk = 490.653D0
    ck = 868.620D0
    dk = 666.792D0
    ek = 180.596D0

    kek = LOG(unit_temp/10000.0D0)
    unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
    loc_kappa = unit_kappa/kappa_bottom
  else
    WRITE(*,*) 'Parameters are way out of whack, something went wrong'
    write(11,*) 'Parameters are way out of whack, something went wrong'
    call error_out
  end if
!  kek = LOG(unit_temp/10000.0D0)
!  unit_kappa = EXP(ak*kek**4 + bk*kek**3 + ck*kek**2 + dk*kek +ek)
!  loc_kappa = unit_kappa/kappa_ref
else
  WRITE(*,*) 'Temperature diverging, exiting',loc_temp,unit_temp,loc_kappa,loc_rho
  CALL error_out
end if

!write(*,*) 'boing boing'

END SUBROUTINE







end module
