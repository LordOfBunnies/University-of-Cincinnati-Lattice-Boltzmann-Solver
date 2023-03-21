SUBROUTINE calculate_alpha(loc_fi,loc_fi_diff,loc_alpha,loc_temp,&
  loc_state,lvl,a,b,c)
!
! Find the sub-grid-scale modeling for the ELBM
!
! This subroutine does a Newton Raphson iteration based on maximizing entropy.
!
! The basic equation is:
! H(fi) = sum( fi * ln(fi))
!
! H(fi) = H(f_alpha)  where  f_alpha = fi + alpha*(fieq-fi)
! H(fi) - H(f_alpha) = 0
!
! Take the derivatives with respective to alpha
!
! loc_fi_diff = fieq - fi
! delta = fi/(fi - fieq) = fi/( -loc_fi_diff)
!
!
!
! Called by: collision
! Calls:
!
USE precise
use timing
USE constants
use grid_data
use linkwise
use quick_calcs
use amr_processes, only: self,amr_max_lvl
use amr_info_holder, only: dt,timestep!,elbm_epsilon
IMPLICIT NONE
real(kind=dp),intent(in) :: loc_fi(a:a,b:b,c:c,1:dir),loc_fi_diff(dir),loc_temp
real(kind=dp) :: loc_alpha,new_alpha,loc_delta,alpha_abso_min,alpha_start,alpha_sum_min
real(kind=dp) :: chi_sum_2,chi_sum_3,chi_sum_1
real(kind=dp) :: first_term_1,first_term_2,first_term_3,first_term_4,second_term_2
real(kind=dp) :: second_term_1,second_term_coeff,second_term_2_1,second_term_2_2,second_term_2_3
real(kind=dp) :: third_term,first_total,second_total,third_total
real(kind=dp) :: chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
real(kind=dp) :: beta_mfp,alpha_lo,alpha_hi,mfp,mft,dim_temp,dim_mu,loc_mu
real(kind=dp) :: derp_fi,derp_fi_diff
!      REAL(kind=dp),ALLOCATABLE :: fi_diff(:),wwfi(:),wwfieq(:),fieq(:)
integer :: i!,iter_counter
integer,intent(in) :: a,b,c,lvl,loc_state
logical :: nonzero_denom(dir)
!
!
!
!
! Define delta so alpha_max and alpha_min can be found
!
!alpha_lo = -1.0D9
!alpha_hi = 1.0D9
alpha_lo = 1.0D0
alpha_hi = 2.125D0
alpha_start = 1.0D-6

!if (timestep(amr_max_lvl) < startup_timesteps) then
!  alpha_abso_min = (alpha_lo - alpha_start)/real(startup_timesteps,kind=dp)*timestep(amr_max_lvl)
!else
!  alpha_abso_min = alpha_lo
!end if

!if (any(isnan(loc_fi)) .or. any(isnan(loc_fi_diff))) then
!  write(*,*) 'PANIC PANIC PANIC!',a,b,c
!end if
!write(*,*) 'incoming fi values',loc_fi
!write(*,*) 'incoming fi diff vlaues',loc_fi_diff

!do i = 1,dir
!  if (abs(loc_fi_diff(i)) <= very_small_number) then
!    nonzero_denom(i) = .false.
!  else
!    nonzero_denom(i) = .true.
!    loc_delta = -loc_fi(a,b,c,i)/loc_fi_diff(i)
!    if (loc_delta < 0) then
!      if (loc_delta > alpha_lo ) then
!        alpha_lo = loc_delta
!      end if
!    else
!      if (loc_delta < alpha_hi ) then
!        alpha_hi = loc_delta
!      end if
!    end if
!  end if
!end do

!if (alpha_lo < 0.10D0) alpha_lo = 0.10D0
!if (alpha_lo < 1.0D0) alpha_lo = 1.0D0
!if (alpha_hi > 2.50D0) alpha_hi = 2.50D0
!if (alpha_hi < 2.0D0) then
!
!!  if (alpha_hi < alpha_abso_min) then
!!    loc_alpha = alpha_abso_min
!!  else
!    loc_alpha = alpha_hi
!!  end if
!!  if (loc_alpha < 0.00001D0) then
!!  write(*,*) 'bad alpha hi',alpha_hi,alpha_lo
!!  write(*,*) 'fis for bad alpha hi',loc_fi(a,b,c,1:dir)
!!  write(*,*) 'fi diffs for bad alpha hi',loc_fi_diff
!!  end if
!  return
!end if
!if (alpha_hi < 1.10D0) alpha_hi = 1.1D0
!if (alpha_lo > alpha_hi) alpha_lo = alpha_hi - 0.000001D0
!if (alpha_lo < 0.50D0) alpha_lo = 0.50D0
!if (alpha_hi > 2.1D0) alpha_hi = 2.1D0
!
!
!
!
!
!upper_alpha = 2.5D0
!lower_alpha = 0.1D0


call quick_sutherlands(loc_mu,loc_temp,.false.,dim_mu,dim_temp)

mfp = dim_mu/p_ref*sqrt(pi*k_boltz*dim_temp/(2*molec_weight))
mft = mfp/sqrt(3*k_boltz*dim_temp/molec_weight)

beta_mfp = dt(lvl)/(2*mft + dt(lvl))
!
first_term_1 = 0.0D0
first_term_2 = 0.0D0
first_term_3 = 0.0D0
first_term_4 = 0.0D0
!
second_term_1 = 0.0D0
second_term_coeff = 0.0D0
second_term_2 = 0.0D0
second_term_2_1 = 0.0D0
second_term_2_2 = 0.0D0
second_term_2_3 = 0.0D0

third_term = 0.0D0

first_total = 0.0D0
second_total = 0.0D0
third_total = 0.0D0

alpha_sum_min = 0.0D0

chi_sum_1 = 0.0D0
chi_sum_2 = 0.0D0
chi_sum_3 = 0.0D0

do i = 1,dir

    if (loc_fi(a,b,c,i) < 0.0D0) then
      derp_fi = abs(loc_fi(a,b,c,i))
    else
      derp_fi = loc_fi(a,b,c,i)
    end if

!    if (derp_fi < 1.0D-6) then
!      derp_fi = 1.0D-6
!    end if

    chi_sum_1 = chi_sum_1 + derp_fi*(2.0D0*(loc_fi_diff(i)/derp_fi)**2/(2+loc_fi_diff(i)/derp_fi))


      chi_sum_2 = chi_sum_2 + derp_fi*(loc_fi_diff(i)/derp_fi)**2

    if (loc_fi_diff(i)/derp_fi <= 0.0D0 .and. loc_fi_diff(i)/derp_fi >= -1.0D0) then
    chi_sum_3 = chi_sum_3 + derp_fi*(loc_fi_diff(i)/derp_fi)**3
    end if
end do

if (abs(chi_sum_3) < plus_epsilon) then
  loc_alpha = 2.0D0
  return
  !chi_sum_3 = minus_epsilon
end if

alpha_sum_min = (chi_sum_2 - sqrt(chi_sum_2**2 - 8.0D0*chi_sum_3*chi_sum_1))/(2.0D0*chi_sum_3)

!write(*,*) 'stuffs', alpha_sum_min,chi_sum_1,chi_sum_2,chi_sum_3

if (alpha_sum_min > alpha_hi) then
  alpha_sum_min = alpha_hi - 1.0D-9
end if

if (timestep(amr_max_lvl) < startup_timesteps) then
  alpha_sum_min = (alpha_sum_min - alpha_start)/real(startup_timesteps,kind=dp)*timestep(amr_max_lvl)
end if
!
!
!
!if (timestep(amr_max_lvl) < startup_timesteps) then
  do i = 1,dir

    if (loc_fi(a,b,c,i) < 0.0D0) then
      derp_fi = abs(loc_fi(a,b,c,i))
    else
      derp_fi = loc_fi(a,b,c,i)
    end if

!    if (derp_fi < 1.0D-6) then
!      derp_fi = 1.0D-6
!    end if

    chi_1 = loc_fi_diff(i)/derp_fi
    chi_2 = chi_1**2
    chi_3 = chi_1**3

    chi_f_2 = chi_2*derp_fi
    chi_f_3 = chi_3*derp_fi
    chi_f_4 = chi_2**2*derp_fi
    chi_f_5 = chi_2*chi_3*derp_fi
    chi_f_6 = chi_3**2*derp_fi

    if (chi_1 <= 0.0D0 .and. chi_1 >= -1.0D0) then
      first_term_1 = first_term_1 + chi_f_3
      first_term_2 = first_term_2 + chi_f_4
      first_term_3 = first_term_3 + chi_f_5
      first_term_4 = first_term_4 + chi_f_6
    end if

    second_term_1 = second_term_1 + chi_f_2

!    second_term_2 = second_term_2 + chi_f_3*(2.0D0/(4 + alpha_lo*chi_1) + &
!      1.0D0/(4 + 2*alpha_lo*chi_1) + 2.0D0/(4 + 3*alpha_lo*chi_1))

    if (chi_1 >= 0.0D0 ) then
      second_term_2 = second_term_2 + chi_f_3*(2.0D0/(4 + alpha_sum_min*chi_1) + &
        1.0D0/(4 + 2*alpha_sum_min*chi_1) + 2.0D0/(4 + 3*alpha_sum_min*chi_1))
    end if
    third_total = third_total - (60*chi_f_2 + 60*chi_f_3 + 11*chi_f_4)/(60 + 90*chi_1 + &
      36*chi_2 + 3*chi_3)

  end do

!else
!do i = 1,dir
!
!  chi_1 = loc_fi_diff(i)/loc_fi(a,b,c,i)
!  chi_2 = chi_1**2
!  chi_3 = chi_1**3
!
!  chi_f_2 = chi_2*loc_fi(a,b,c,i)
!  chi_f_3 = chi_3*loc_fi(a,b,c,i)
!  chi_f_4 = chi_2**2*loc_fi(a,b,c,i)
!  chi_f_5 = chi_2*chi_3*loc_fi(a,b,c,i)
!  chi_f_6 = chi_3**2*loc_fi(a,b,c,i)
!
!  first_term_1 = first_term_1 + chi_f_3
!  first_term_2 = first_term_2 + chi_f_4
!  first_term_3 = first_term_3 + chi_f_5
!  first_term_4 = first_term_4 + chi_f_6
!
!  second_term_1 = second_term_1 + chi_f_2
!
!!  second_term_2 = second_term_2 + chi_f_3*(2.0D0/(4 + alpha_lo*chi_1) + &
!!    1.0D0/(4 + 2*alpha_lo*chi_1) + 2.0D0/(4 + 3*alpha_lo*chi_1))
!  second_term_2 = second_term_2 + chi_f_3*(2.0D0/(4 + alpha_abso_min*chi_1) + &
!    1.0D0/(4 + 2*alpha_abso_min*chi_1) + 2.0D0/(4 + 3*alpha_abso_min*chi_1))
!
!  third_total = third_total - (60*chi_f_2 + 60*chi_f_3 + 11*chi_f_4)/(60 + 90*chi_1 + &
!    36*chi_2 + 3*chi_3)
!
!!  write(*,*) 'Quadratic stuff, chis',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!!  write(*,*) 'Quadratic constants',beta_mfp,alpha_hi,alpha_lo,mfp,mft,dim_mu,dim_temp,molec_weight,&
!!    k_boltz,p_ref,t_ref,loc_mu,loc_temp,self
!end do
!
!end if
!
first_total = -beta_mfp**2 * (first_term_1/6 - alpha_hi*beta_mfp*first_term_2/12 + &
  alpha_hi**2*beta_mfp**2*first_term_3/20 - alpha_hi**3*beta_mfp**3*first_term_4/5 )

!second_total = second_term_1/2 - 2*alpha_lo*beta_mfp**2*second_term_2/15
second_total = second_term_1/2 - 2*alpha_sum_min*beta_mfp**2*second_term_2/15
!
!
if (abs(first_total) < very_small_number) then
  loc_alpha = 2.0D0
else
  new_alpha = (-second_total + sqrt(second_total**2 - 4*first_total*third_total))/(2*first_total)
  if (new_alpha < alpha_sum_min) then
!    loc_alpha = alpha_lo
    loc_alpha = alpha_sum_min
  else if (new_alpha > alpha_hi) then
    loc_alpha = alpha_hi
  else
    loc_alpha = new_alpha
  end if
end if
!
!
!
if (isnan(new_alpha)) then
  write(*,*) 'Quadratic stuff, first term',first_term_1,first_term_2,first_term_3,first_term_4
  write(*,*) 'Quadratic stuff, second term',second_term_1,second_term_2
  write(*,*) 'Quadratic stuff, third term',third_total
  write(*,*) 'Alpha constants',beta_mfp,alpha_hi,alpha_lo,alpha_abso_min,alpha_sum_min,&
    mfp,mft,dim_mu,dim_temp,molec_weight,&
    k_boltz,p_ref,t_ref,loc_mu,loc_temp,loc_state,a,b,c,lvl,self
  write(*,*) 'Quadratic equation issue, shutting down',first_total,second_total,third_total
  write(*,*) 'fis for bad alpha',loc_fi
  write(*,*) 'fi_diffs for bad alpha',loc_fi_diff
  write(*,*) 'da fluff!',a,b,c,self,lvl
  call error_out
end if

!if (a == 345 .and. b == 249 .and. c == 3) then
!  write(*,*) 'alpha information middle,center',first_total,second_total,third_total,loc_alpha,new_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'alpha subcomponents middle,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information middle,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, middle,center',loc_fi
!  write(*,*) 'fi diffs, middle,center',loc_fi_diff
!end if

!if (a == 347 .and. b == 251 .and. c == 3) then
!  write(*,*) 'alpha information middle,center',first_total,second_total,third_total,loc_alpha,new_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'alpha subcomponents middle,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information middle,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, middle,center',loc_fi
!  write(*,*) 'fi diffs, middle,center',loc_fi_diff
!end if

!if (a == 349 .and. b == 253 .and. c == 3) then
!  write(*,*) 'alpha information middle,center',first_total,second_total,third_total,loc_alpha,new_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'alpha subcomponents middle,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information middle,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, middle,center',loc_fi
!  write(*,*) 'fi diffs, middle,center',loc_fi_diff
!end if

END SUBROUTINE
!real(kind=dp) :: !max_diff,f_alpha,f_alpha_prev,mod_alpha,alpha_prev
!REAL(kind=dp) :: hf,hf_prime,hf_double,old_alpha,hf_prev
!real(kind=dp) :: alpha_min,alpha_max,loc_delta
!real(kind=dp) :: a1,a2,a3,a4,a_rat,base_a

!delta = -loc_fi/loc_fi_diff
!diff_max = minloc(delta,MASK = delta < 0.)
!diff_min = maxloc(delta,MASK = delta > 0.)
!IF (diff_min(1) <= 2.) THEN
!  loc_alpha = diff_min(1)
!  write(*,*) 'local alpha',loc_alpha
!  RETURN
!END IF
!
! Pre-define a maximum difference so the loop is entered properly
!

!write(*,*) 'Quadratic equation issue, shutting down',first_total,second_total,third_total,new_alpha

!if (a == 277 .and. b == 124 .and. c == 132) then
!  write(*,*) 'alpha information middle,center',first_total,second_total,third_total,loc_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'alpha subcomponents middle,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information middle,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, middle,center',loc_fi
!  write(*,*) 'fi diffs, middle,center',loc_fi_diff
!end if

!if (loc_alpha < 0.1D0) then
!  write(*,*) 'alpha weirdness',first_total,second_total,third_total,loc_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'weird alpha fis',loc_fi
!  write(*,*) 'weird alpha fi diffs',loc_fi_diff
!
!end if



!write(*,*)  'alpha information front,center',first_total,second_total,third_total,loc_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'alpha subcomponents front,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!  write(*,*) 'fis, front,center',loc_fi
!  write(*,*) 'fi diffs, front,center',loc_fi_diff

!if (a == 58 .and. b == 67 .and. c == 5) then
!  write(*,*) 'alpha information front,center',first_total,second_total,third_total,loc_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,self,a,b,c
!  write(*,*) 'alpha subcomponents front,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information front,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, front,center',loc_fi
!  write(*,*) 'fi diffs, front,center',loc_fi_diff
!end if



!if (a == 72 .and. b == 32 .and. c == 32) then
!  write(*,*) 'alpha information near,center',first_total,second_total,third_total,loc_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,a,b,c
!  write(*,*) 'alpha subcomponents near,center',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information near,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, near,center',loc_fi
!  write(*,*) 'fi diffs, near,center',loc_fi_diff
!end if
!if (a == 73 .and. b == 32 .and. c == 32) then
!  write(*,*) 'alpha information back',first_total,second_total,third_total,loc_alpha,&
!    alpha_lo,alpha_hi,beta_mfp,mfp,mft,a,b,c
!  write(*,*) 'alpha subcomponents back',first_term_1,first_term_2,first_term_3,first_term_4,second_term_1,&
!    second_term_2
!!  write(*,*) 'chi information near,center',chi_1,chi_2,chi_3,chi_f_2,chi_f_3,chi_f_4,chi_f_5,chi_f_6
!  write(*,*) 'fis, back',loc_fi
!  write(*,*) 'fi diffs, back',loc_fi_diff
!end if
!
!
! Set initial conditions
!
!alpha_prev = loc_alpha + 0.01D0
!!
!!
!!
!! Define delta so alpha_max and alpha_min can be found
!!
!alpha_min = -1.0D9
!alpha_max = 1.0D9
!!write(*,*) 'incoming fi values',loc_fi
!!write(*,*) 'incoming fi diff vlaues',loc_fi_diff
!
!do i = 1,dir
!  if (loc_fi_diff(i) <= very_small_number .and. loc_fi_diff(i) >= -very_small_number) then
!    nonzero_denom(i) = .false.
!  else
!    nonzero_denom(i) = .true.
!    loc_delta = -loc_fi(a,b,c,i)/loc_fi_diff(i)
!    if (loc_delta < 0) then
!      if (loc_delta > alpha_min ) then
!        alpha_min = loc_delta
!      end if
!    else
!      if (loc_delta < alpha_max ) then
!        alpha_max = loc_delta
!      end if
!    end if
!  end if
!end do
!
!!write(*,*) nonzero_denom
!
!if (alpha_max < 2) then
!  loc_alpha = alpha_max
!!  write(*,*) 'alpha maxed out',loc_alpha
!!  write(*,*) 'local fi',loc_fi
!!  write(*,*) 'local fi diff',loc_fi_diff
!  return
!end if
!
!if (any(nonzero_denom)) then
!!  write(*,*) 'scoot'
!!
!  max_diff = 1.0D0
!  iter_counter = 0
!!
!! Newton-Raphson iteration for the subgrid scale model factor
!!
!  do while (max_diff > 1.0D-7)
!    hf_prime = 0.0D0
!    hf = 0.0D0
!    hf_double = 0.0D0
!    hf_prev = 0.0D0
!    old_alpha = loc_alpha
!!
!!    write(*,*) 'snoot',max_diff,loc_alpha,iter_counter,self
!!    write(*,*) 'local fi',loc_fi,iter_counter,self
!!    write(*,*) 'local fi diff',loc_fi_diff,iter_counter,self
!
!  iter_counter = iter_counter+1
!    DO i = 1,dir
!!
!! Calculate things that are required every time so as not to repeat the same calcs
!!
!      f_alpha = loc_fi(a,b,c,i)+loc_alpha*loc_fi_diff(i)
!      f_alpha_prev = loc_fi(a,b,c,i)+alpha_prev*loc_fi_diff(i)
!!
!! Calculate hf and hf_alpha
!!
!!      hf = hf + loc_fi(a,b,c,i)*LOG(loc_fi(a,b,c,i)) - f_alpha*LOG(f_alpha)
!!      hf_prev = hf_prev + loc_fi(a,b,c,i)*LOG(loc_fi(a,b,c,i)) - f_alpha_prev*LOG(f_alpha_prev)
!      hf = hf + loc_fi(a,b,c,i)*LOG(loc_fi(a,b,c,i)/elbm_epsilon(lvl)) - f_alpha*LOG(f_alpha/elbm_epsilon(lvl))
!!
!!      hf_prime = hf_prime - loc_fi_diff(i)*(1.0D0+LOG(f_alpha))
!      hf_prime = hf_prime - loc_fi_diff(i)*(1.0D0 + LOG(f_alpha/elbm_epsilon(lvl)))
!
!!      hf_double = hf_double - loc_fi_diff(i)**2/f_alpha
!      write(*,*) 'internal values',f_alpha,hf,hf_prime,loc_fi(a,b,c,i),loc_fi_diff(i),loc_alpha,iter_counter,a,b,c,self
!!      if (self == 0) then
!!      write(*,*) 'alpha calculations',hf,hf_prev,loc_fi(a,b,c,i),loc_fi_diff(i),f_alpha,f_alpha_prev,&
!!        LOG(loc_fi(a,b,c,i)),LOG(f_alpha),log(f_alpha_prev),loc_alpha,alpha_prev,iter_counter,i,self,a,b,c
!!      end if
!      write(*,*) 'log values and limits',loc_fi(a,b,c,i)*LOG(loc_fi(a,b,c,i)),f_alpha*LOG(f_alpha),&
!        alpha_max,alpha_min,iter_counter,a,b,c,self
!      if (f_alpha < 0.0D0 ) then
!        write(*,*) 'shit shit shit shit'
!      end if
!
!    END DO
!!
!! Do the Newton-Raphson iteration
!!
!!    if (hf_prime <= very_small_number .and. hf_prime >= -very_small_number) then
!!      loc_alpha = 2.0D0
!!      write(*,*) 'alpha, fi = fieq',loc_alpha,hf,hf_prime,hf_double,alpha_max,alpha_min,max_diff,iter_counter,self
!!!      write(*,*) 'local fi',loc_fi
!!!      write(*,*) 'local fi diff',loc_fi_diff
!!      return
!!    end if
!!
!! Standard Newton-Raphson iteration
!!
!    new_alpha = loc_alpha - hf/hf_prime
!!    mod_alpha = loc_alpha - (hf*hf_prime)/(hf_prime**2 - hf*hf_double)
!!    new_alpha = loc_alpha - hf*(loc_alpha-alpha_prev)/(hf-hf_prev)
!    write(*,*) 'alpha junk',loc_alpha,new_alpha,old_alpha,hf,hf_prime,hf_double,a,b,c
!!    write(*,*) 'alpha junk',loc_alpha,new_alpha,old_alpha,hf,hf_prev,a,b,c
!!    alpha_prev = loc_alpha
!!    loc_alpha = new_alpha
!!
!!
!!
!    if (new_alpha > alpha_max) then
!      loc_alpha = alpha_max
!!      write(*,*) 'alpha above maximum alpha',loc_alpha,new_alpha,alpha_max,iter_counter,self
!      return
!    else if (new_alpha < alpha_min) then
!      loc_alpha = alpha_min
!!      write(*,*) 'alpha below minimum alpha',loc_alpha,new_alpha,alpha_min,iter_counter,self
!      return
!    end if
!
!
!!    max_diff = ABS(alpha_prev - loc_alpha)
!    max_diff = ABS(new_alpha - loc_alpha)
!    loc_alpha = new_alpha
!    write(*,*) 'normal alpha calc',loc_alpha,new_alpha,old_alpha,hf,hf_prime,hf_double,&
!      alpha_max,alpha_min,max_diff,iter_counter,self
!!    write(*,*) 'normal alpha calc',loc_alpha,new_alpha,old_alpha,hf,hf_prev,&
!!      alpha_max,alpha_min,max_diff,iter_counter,self
!  end do
!!
!! if fi=fieq
!!
!else
!!  write(*,*) 'blarg!!!'
!  loc_alpha = 2.0D0
!!  write(*,*) 'fi = fieq',loc_alpha
!  return
!end if

!a1 = 0.0D0
!a2 = 0.0D0
!a3 = 0.0D0
!a4 = 0.0D0
!do i = 1,dir
!  a_rat = loc_fi_diff(i)/loc_fi(a,b,c,i)
!  base_a = loc_fi_diff(i)**2/loc_fi(a,b,c,i)
!  a1 = a1 + base_a
!  a2 = a2 + base_a*a_rat
!  a3 = a3 + base_a*a_rat**2
!  a4 = a4 + base_a*a_rat**3
!end do
!
!a1 = a1/2.0D0
!a2 = -a2/6.0D0
!a3 = a3/12.0D0
!a4 = -a4/20.0D0
!
!loc_alpha = 2.0D0 - 4.0D0/a1*(a2 + 4.0D0*a2**2/a1 - 2.0D0*a3 + 20.0D0*a2/a1*(a3-a2**2/a1) - 4.0D0*a4)

!write(*,*) 'New alpha value',loc_alpha
