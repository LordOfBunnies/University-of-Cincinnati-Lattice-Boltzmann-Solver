subroutine target_finder(loc_state,u_pass,v_pass,w_pass,temp_pass,&
  loc_fi,rho_tgt,u_tgt,v_tgt,w_tgt,temp_tgt,a,b,c)
!
! Finds target values for a given node and stores them in a temporary multifab
!
!
! Called by: interpolative_bounceback
! Calls:
!
use precise
use constants
use grid_data
use linkwise
use freestream_values
implicit none
real(kind=dp),intent(out) :: rho_tgt,u_tgt,v_tgt,w_tgt,temp_tgt
!real(kind=dp) :: squirrel_1,squirrel_2,squirrel_3,squirrel_4
real(kind=dp) :: u_pass(dir),v_pass(dir),w_pass(dir),temp_pass(dir),&
  latt_dist,loc_fi(a:a,b:b,c:c,1:dir)
integer :: i,solid_count
integer,intent(in) :: loc_state(a:a,b:b,c:c,1:dir),a,b,c

solid_count = 0
rho_tgt = 0.0D0
u_tgt = 0.0D0
v_tgt = 0.0D0
w_tgt = 0.0D0
temp_tgt = 0.0D0

!write(*,*) 'Passed info, u ',u_pass,a,b,c
!write(*,*) 'Passed info, v ',v_pass,a,b,c
!write(*,*) 'Passed info, w ',w_pass,a,b,c
!write(*,*) 'Passed info, temp ',temp_pass,a,b,c
!write(*,*) 'Passed states', loc_state(a,b,c,1:dir),a,b,c

if (.not. shifted) then
  do i = 1,dir

    if (loc_state(a,b,c,i) <= -1001) then
      latt_dist = ABS(real(loc_state(a,b,c,i)+1001.0D0,kind=dp)/1.0D9)
!    write(*,*) 'state info ',loc_state(a,b,c,i),i,a,b,c
!    write(*,*) 'local lattice distance',latt_dist,solid_count,a,b,c,i

      solid_count = solid_count+1
      if (loc_state(a,b,c,opp(i)) >= 0) then
        u_tgt = u_tgt + (latt_dist * u_pass(opp(i))) / (latt_dist+1.0D0)

        v_tgt = v_tgt + (latt_dist * v_pass(opp(i))) / (latt_dist+1.0D0)

        w_tgt = w_tgt + (latt_dist * w_pass(opp(i))) / (latt_dist+1.0D0)

!          temp_tgt = temp_tgt + (latt_dist*temp_pass(opp(i)) + temperature) &
!            / (latt_dist+1.0D0)
        temp_tgt = temp_tgt+(latt_dist*temp_pass(opp(i)) + temp_pass(opp(i)))/ (latt_dist+1.0D0)

      else
        u_tgt = u_tgt + (latt_dist * u_pass(1)) / (latt_dist+1.0D0)

        v_tgt = v_tgt + (latt_dist * v_pass(1)) / (latt_dist+1.0D0)

        w_tgt = w_tgt + (latt_dist * w_pass(1)) / (latt_dist+1.0D0)

        temp_tgt = temp_tgt + (latt_dist*temp_pass(1) + temperature) / (latt_dist+1.0D0)
      end if

!    write(*,*) 'in loop targets',u_tgt,v_tgt,w_tgt,temp_tgt
!    write(*,*) 'tgt opps ',u_pass(opp(i)),v_pass(opp(i)),w_pass(opp(i)),temp_pass(opp(i)),&
!      loc_state(a,b,c,i),loc_state(a,b,c,opp(i))
!    if (isnan(temp_tgt)) call error_out()
    end if

    if (loc_state(a,b,c,opp(i)) < -1000 ) then
      rho_tgt = rho_tgt + loc_fi(a,b,c,opp(i))
    else
      rho_tgt = rho_tgt + loc_fi(a,b,c,i)
    end if

  end do

  u_tgt = u_tgt/solid_count
  v_tgt = v_tgt/solid_count
  w_tgt = w_tgt/solid_count
  temp_tgt = temp_tgt/solid_count

!write(*,*) 'local target values',u_tgt,v_tgt,w_tgt,temp_tgt,rho_tgt,a,b,c!,loc_state
!
! Shifted lattice
!
else
  do i = 1,dir

    if (loc_state(a,b,c,i) <= -1001) then
      latt_dist = ABS(real(loc_state(a,b,c,i)+1001.0D0,kind=dp)/1.0D9)
!    write(*,*) 'state info ',loc_state(a,b,c,i),i,a,b,c
!    write(*,*) 'local lattice distance',latt_dist,solid_count,a,b,c,i

      solid_count = solid_count+1
      if (loc_state(a,b,c,i) <= -1001 .and. loc_state(a,b,c,erp(i)) == -666) then
        u_tgt = u_tgt + (latt_dist * u_pass(i)) / (latt_dist+1.0D0)

        v_tgt = v_tgt + (latt_dist * v_pass(i)) / (latt_dist+1.0D0)

        w_tgt = w_tgt + (latt_dist * w_pass(i)) / (latt_dist+1.0D0)

        temp_tgt = temp_tgt + (latt_dist*temp_pass(i) + temp_pass(i))/(latt_dist+1.0D0)

      else
        u_tgt = u_tgt + (latt_dist * u_pass(5)) / (latt_dist+1.0D0)

        v_tgt = v_tgt + (latt_dist * v_pass(5)) / (latt_dist+1.0D0)

        w_tgt = w_tgt + (latt_dist * w_pass(5)) / (latt_dist+1.0D0)

        temp_tgt = temp_tgt + (latt_dist*temp_pass(5) + temp_pass(5)) / (latt_dist+1.0D0)
      end if

!    write(*,*) 'in loop targets',u_tgt,v_tgt,w_tgt,temp_tgt,latt_dist,i,loc_state(a,b,c,i),a,b,c
!    write(*,*) 'tgt opps ',u_pass(opp(i)),v_pass(opp(i)),w_pass(opp(i)),temp_pass(opp(i)),&
!      loc_state(a,b,c,i),loc_state(a,b,c,opp(i))
!    if (isnan(temp_tgt)) call error_out()
    end if

    if (loc_state(a,b,c,opp(i)) < -1000 ) then
      rho_tgt = rho_tgt + loc_fi(a,b,c,erp(i))
    else
      rho_tgt = rho_tgt + loc_fi(a,b,c,i)
    end if

  end do

  u_tgt = u_tgt/solid_count
  v_tgt = v_tgt/solid_count
  w_tgt = w_tgt/solid_count
  temp_tgt = temp_tgt/solid_count
end if

if (isnan(temp_tgt)) then
  write(*,*) 'target finder state values ',loc_state,a,b,c
  write(*,*) 'target values',u_tgt,v_tgt,w_tgt,temp_tgt,solid_count

  call error_out()
end if



end subroutine

!    write(*,*) u_tgt,latt_dist,u_pass(opp(i))
!    squirrel_1 = latt_dist*v_pass(opp(i))
!    squirrel_2 = latt_dist+1
!    squirrel_3 = squirrel_1/squirrel_2
!    squirrel_4 = v_tgt+squirrel_3
!    write(*,*) 'Squirrel 1 ',squirrel_1,' Squirrel 2',squirrel_2,' Squirrel 3',squirrel_3
!    write(*,*) 'Squirrel target',squirrel_4
!    write(*,*) loc(u_tgt),loc(v_tgt),loc(w_tgt),loc(temp_tgt)

