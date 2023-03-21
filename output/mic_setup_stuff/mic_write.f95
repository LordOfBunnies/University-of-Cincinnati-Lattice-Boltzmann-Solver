SUBROUTINE mic_write(a)
!
! This does microphone writes, it is given the timestep and current mic
!
! Called by: lb_output
! Calls: error_out
!
USE output_data
USE linkwise
USE precise
USE constants
USE nml_output
USE freestream_values
use mpi
use amrex_amr_module
use amr_info_holder
use amr_processes, only: commune,magic_box
IMPLICIT NONE
INTEGER :: i,a,unit_file,ier,b,output_counter,fine
integer :: mic_cur_lvl
REAL(KIND=dp) :: loc_rho,loc_u,loc_v,loc_w,loc_press,loc_temp
REAL(KIND=dp),ALLOCATABLE :: data_output(:)
CHARACTER(len=64) :: mic_file,file_path,junk
logical :: inside
!      CHARACTER(len=8) :: p_num,temp_num,rho_num,u_vel_num,v_vel_num,&
!        w_vel_num
type(amrex_mfiter) :: mfi
type(amrex_box) :: fox

fine = amrex_get_finest_level()

if (fine >= mic_lvl(a)) then
    mic_cur_lvl = mic_lvl(a)
  else
    mic_cur_lvl = amrex_get_finest_level()
end if
inside = .false.
junk = '(I5.5)'
unit_file = 29

call amrex_mfiter_build(mfi,mfrho(mic_cur_lvl),tiling=.false.)
do while (mfi%next())
  fox = mfi%validbox()
  inside = magic_box(fox%lo,fox%hi,mic_nodes(1,1,a))
  if (inside) then

!      mic_directory = 'microphones/'
 !MOD(a,100)+100
!      mic_file = 'Microphone_'//CHAR(a)//'.txt'
    WRITE(mic_file,junk) a
    file_path = 'microphone'//TRIM(mic_file)//'.txt'

    OPEN(UNIT = unit_file, FILE=file_path, FORM = 'FORMATTED',&
      ACTION = 'WRITE',POSITION='APPEND',IOSTAT=ier)
    IF (ier /=0) THEN
      WRITE(*,*) 'Something went wrong opening microphone file ',&
        mic_file, ' for microphone ', a
      CALL error_out
    END IF
!
! Initialize the character output variables.
!
    output_counter = 0

!
! Averaging for local pressure
!
    IF (mic_type(1,a)) THEN
      loc_press = 0.0D0
      press => mfpress(mic_cur_lvl)%dataptr(mfi)
      output_counter = output_counter +1
      DO b = 1,8
        loc_press = loc_press + mic_weights(b,a)*press(mic_nodes(1,b,a),&
          mic_nodes(2,b,a),mic_nodes(3,b,a),1)
!        IF (mic_links(b,a) /= 0) THEN
!          loc_press = loc_press + mic_weights(b,a)*rho(mic_links(b,a))&
!            *gas_const_R*temperature
!        END IF
      END DO

!        p_num = CHAR(loc_press)
    END IF
!
! Averaging for local temperatures
!
!      output_counter = 0
    IF (mic_type(2,a)) THEN
      loc_temp = 0.0D0
      temp => mftemp(mic_cur_lvl)%dataptr(mfi)
      output_counter = output_counter +1
      DO b = 1,8
        loc_temp = loc_temp + mic_weights(b,a)*temp(mic_nodes(1,b,a),&
          mic_nodes(2,b,a),mic_nodes(3,b,a),1)
!        IF (mic_links(b,a) /= 0) THEN
!          loc_temp = loc_temp + mic_weights(b,a)*temperature
!        END IF
      END DO
    END IF
!        temp_num = CHAR(loc_temp)
!
! Averaging for local density
!
    IF (mic_type(3,a)) THEN
      loc_rho = 0.0D0
      rho => mfrho(mic_cur_lvl)%dataptr(mfi)
      output_counter = output_counter +1
      DO b = 1,8
        loc_rho = loc_rho + mic_weights(b,a)*rho(mic_nodes(1,b,a),&
          mic_nodes(2,b,a),mic_nodes(3,b,a),1)
!        IF (mic_links(b,a) /= 0) THEN
!          loc_rho = loc_rho + mic_weights(b,a)*rho(mic_links(b,a))
!        END IF
      END DO

    !        rho_num = CHAR(loc_rho)
    END IF
    !
    ! Averaging for local u-velocity
    !
    IF (mic_type(4,a)) THEN
      loc_u = 0.0D0
      u_vel => mfu_vel(mic_cur_lvl)%dataptr(mfi)
      output_counter = output_counter +1
      DO b = 1,8
        loc_u = loc_u + mic_weights(b,a)*u_vel(mic_nodes(1,b,a),&
          mic_nodes(2,b,a),mic_nodes(3,b,a),1)
!        IF (mic_links(b,a) /= 0) THEN
!          loc_u = loc_u + mic_weights(b,a)*u_vel(mic_links(b,a))
!        END IF
      END DO

    !        u_vel_num = CHAR(loc_u)
    END IF
    !
    ! Averaging for local v-velocity
    !
    IF (mic_type(5,a)) THEN
      loc_v = 0.0D0
      v_vel => mfv_vel(mic_cur_lvl)%dataptr(mfi)
      output_counter = output_counter +1
      DO b = 1,8
        loc_v = loc_v + mic_weights(b,a)*v_vel(mic_nodes(1,b,a),&
          mic_nodes(2,b,a),mic_nodes(3,b,a),1)
!        IF (mic_links(b,a) /= 0) THEN
!          loc_v = loc_v + mic_weights(b,a)*v_vel(mic_links(b,a))
!        END IF
      END DO
    !        v_vel_num = CHAR(loc_v)
    END IF
    !
    ! Averaging for local w-velocity
    !
    IF (mic_type(6,a)) THEN
      loc_w = 0.0D0
      w_vel => mfw_vel(mic_cur_lvl)%dataptr(mfi)
      output_counter = output_counter +1
      DO b = 1,8
        loc_w = loc_w + mic_weights(b,a)*w_vel(mic_nodes(1,b,a),&
          mic_nodes(2,b,a),mic_nodes(3,b,a),1)
!        IF (mic_links(b,a) /= 0) THEN
!          loc_w = loc_w + mic_weights(b,a)*w_vel(mic_links(b,a))

!        END IF
      END DO
    !        w_vel_num = CHAR(loc_w)
    END IF
    if (allocated(data_output)) deallocate(data_output)
    ALLOCATE(data_output(output_counter))

    output_counter = 0
    DO i = 1,6
      IF (mic_type(i,a)) THEN
        output_counter = output_counter + 1
        SELECT CASE(i)
        CASE(1)
          data_output(output_counter) = loc_press
        CASE(2)
          data_output(output_counter) = loc_temp
        CASE(3)
          data_output(output_counter) = loc_rho
        CASE(4)
          data_output(output_counter) = loc_u
        CASE(5)
          data_output(output_counter) = loc_v
        CASE(6)
          data_output(output_counter) = loc_w
        END SELECT


      END IF
    END DO

    !      WRITE(unit_file) p_num,' ',temp_num,' ',rho_num,' ',&
    !        u_vel_num,' ', v_vel_num,' ',w_vel_num
    IF (output_counter == 1) THEN
      WRITE(unit_file,100) data_output(1)
    ELSE IF (output_counter == 2) THEN
      WRITE(unit_file,101) data_output(1),data_output(2)
    ELSE IF (output_counter == 3) THEN
      WRITE(unit_file,102) data_output(1),data_output(2),data_output(3)
    ELSE IF (output_counter == 4) THEN
      WRITE(unit_file,103) data_output(1),data_output(2),data_output(3),&
        data_output(4)
    ELSE IF (output_counter == 5) THEN
      WRITE(unit_file,104) data_output(1),data_output(2),data_output(3),&
        data_output(4),data_output(5)
    ELSE IF (output_counter == 6) THEN
      WRITE(unit_file,105) data_output(1),data_output(2),data_output(3),&
        data_output(4),data_output(5),data_output(6)

    END IF


    CLOSE(unit_file)
  else
    CYCLE
  end if
end do
call amrex_mfiter_destroy(mfi)
call mpi_barrier(commune,ier)


 100  FORMAT (F12.5)
 101  FORMAT (F12.5,' ',F12.5)
 102  FORMAT (F12.5,' ',F12.5,' ',F12.5)
 103  FORMAT (F12.5,' ',F12.5,' ',F12.5,' ',F12.5)
 104  FORMAT (F12.5,' ',F12.5,' ',F12.5,' ',F12.5,' ',F12.5)
 105  FORMAT (F12.5,' ',F12.5,' ',F12.5,' ',F12.5,' ',F12.5,' ',F12.5)
 106  FORMAT (3F9.7,I8,2F9.7,I8)

END SUBROUTINE
