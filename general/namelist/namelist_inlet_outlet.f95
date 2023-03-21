SUBROUTINE namelist_inlet_outlet
!
! Subroutine loads the inlet and outlet information
!
! Called by: namelist_global
! Calls: error_out
!
USE nml_inlet_outlet
USE precise
IMPLICIT NONE
INTEGER :: inlet_flow_dir,outlet_flow_dir
INTEGER :: istat,i,in_type,out_type
REAL(KIND=dp) :: inlet_val_1,inlet_val_2,inlet_val_3,inlet_val_4,&
  inlet_val_5,outlet_val_1,outlet_val_2,outlet_val_3,outlet_val_4,&
  outlet_val_5
CHARACTER(len=80) :: inlet_name,outlet_name
!
!Declare the NAMELISTs for all the inlets and outlets
!
NAMELIST /INLET_OUTLET/ num_inlets,num_outlets

NAMELIST /INLET_DATA/ inlet_val_1,inlet_val_2,inlet_val_3,inlet_val_4,&
  inlet_val_5,inlet_name,in_type,inlet_flow_dir

NAMELIST /OUTLET_DATA/ outlet_val_1,outlet_val_2,outlet_val_3,outlet_val_4,&
  outlet_val_5,outlet_name,out_type,outlet_flow_dir
! Initialize the first round of namelist reads
num_inlets = 1
num_outlets = 1
! Read the numbers of inlets and outlets
READ(2,NML=INLET_OUTLET,IOSTAT=istat)
IF (istat /= 0) THEN
  WRITE(*,*) 'Error reading inlet/outlet data, exiting.'
  WRITE(11,*) 'Error reading inlet/outlet data, exiting.'
  WRITE(*,*) istat
  CALL error_out
END IF
!
! Allocate all the arrays stored in the module
!
ALLOCATE(inlet_type(2,num_inlets))
ALLOCATE(outlet_type(2,num_outlets))
ALLOCATE(in_name(num_inlets))
ALLOCATE(out_name(num_outlets))
ALLOCATE(inlet_flow_val(5,num_inlets))
ALLOCATE(outlet_flow_val(5,num_outlets))

inlet_side = .FALSE.
outlet_side = .FALSE.
!
! Make sure there's at least one inlet and go from there
! Loops through the inlet data reading as many as there are.
!
! Inlet/outlet types and value meanings
! type 1 = momentum bc, val_1 = rho, val_2 = velocity, temp = freestream
! type 2 = pressure bc, val_1 = pressure, val_2 = temp, rho = calculated
! type 3 = freestream bc, values mean nothing and freestream velocity and rho used
! type 4 = temperature bc, val_1 = velocity, val_2 = temperature, rho = freestream
! type 5 = periodic (flow from outlet goes to the paired outlet on the other side)
!                      (outlet only)
! type 6 = pure outflow (outlet only), nodes on that side are designated as fluid (assuming they are)
!            and allowed to propagate back into the fluid.  The nodes calculate the new fout, and that
!            becomes the new fi WITHOUT MOVING for these nodes.
!
! type 8 = freestream edge bc, essentially a far field condition flowing in the same
!            direction as the primary flow direction
!
! type 10 = viscous wall
! type 11 = inviscid wall
!
!
!
! Inlet flow direction (same for outlets)
! 1 = +x-direction
! 2 = +y-direction
! 3 = +z-direction
! 4 = -x-direction
! 5 = -y-direction
! 6 = -z-direction
!
IF (num_inlets < 1) THEN
  WRITE(11,*) 'Too few inlets, you need at least 1 inlet, exiting.'
  CALL error_out
ELSE
  DO i =1,num_inlets
    inlet_val_1 = 0.0D0
    inlet_val_2 = 0.0D0
    inlet_val_3 = 0.0D0
    inlet_val_4 = 0.0D0
    inlet_val_5 = 0.0D0
    inlet_name = 'Basic Inlet'
    in_type = 1
    inlet_flow_dir = 1

    READ(2,NML=INLET_DATA,IOSTAT=istat)
    IF (istat /= 0) THEN
      WRITE(*,*) 'Error reading inlet data at inlet ',i,&
        ', exiting.'
      WRITE(11,*) 'Error reading inlet data at inlet ',i,&
        ', exiting.'
      WRITE(*,*) istat
      CALL error_out
    ELSE
      inlet_flow_val(1,i) = inlet_val_1
      inlet_flow_val(2,i) = inlet_val_2
      inlet_flow_val(3,i) = inlet_val_3
      inlet_flow_val(4,i) = inlet_val_4
      inlet_flow_val(5,i) = inlet_val_5
      in_name(i) = inlet_name
      inlet_type(1,i) = in_type
      inlet_type(2,i) = inlet_flow_dir
      IF (inlet_type(2,i) < 0 .OR. inlet_type(2,i) > 6) THEN
        WRITE(*,*) 'Inlet flow direction out of range, must be between 1 and 6.'
        CALL error_out
      END IF


  inlet_side(inlet_flow_dir) = .TRUE.
!
! This will do the output to the output file, essentially reiterating
! everything the user puts in for the inlets.
!
      WRITE(11,101) i,in_name(i)

      IF (inlet_type(1,i) == 1) THEN
        WRITE(11,104) i
      ELSE IF (inlet_type(1,i) == 2) THEN
        WRITE(11,105) i
      ELSE IF (inlet_type(1,i) == 3) THEN
        WRITE(11,106) i
      ELSE IF (inlet_type(1,i) == 4) THEN
        WRITE(11,107) i
      ELSE IF (inlet_type(1,i) == 5) THEN
        WRITE(11,109) i
        WRITE(11,*) 'Only outlets can be periodic boundary conditions.'
        WRITE(11,*) 'Change this to an outlet and run program again',&
          ', exiting'
        CALL error_out
      else if (inlet_type(1,i) == 7) then
        write(11,122) i


      ELSE IF (inlet_type(1,i) == 10) THEN
        WRITE(11,110) i
      ELSE IF (inlet_type(1,i) == 11) THEN
        WRITE(11,109) i

      END IF

      WRITE(11,102) i,inlet_flow_val(1,i)
      WRITE(11,103) i,inlet_flow_val(2,i)
      WRITE(11,103) i,inlet_flow_val(3,i)
      WRITE(11,103) i,inlet_flow_val(4,i)
      WRITE(11,103) i,inlet_flow_val(5,i)

      IF (inlet_type(2,i) == 1) THEN
        WRITE(11,*) 'Inlet is in the positive x-direction.'
      ELSE IF (inlet_type(2,i) == 2) THEN
        WRITE(11,*) 'Inlet is in the positive y-direction.'
      ELSE IF (inlet_type(2,i) == 3) THEN
        WRITE(11,*) 'Inlet is in the positive z-direction.'
      ELSE IF (inlet_type(2,i) == 4) THEN
        WRITE(11,*) 'Inlet is in the negative x-direction.'
      ELSE IF (inlet_type(2,i) == 5) THEN
        WRITE(11,*) 'Inlet is in the negative y-direction.'
      ELSE IF (inlet_type(2,i) == 6) THEN
        WRITE(11,*) 'Inlet is in the negative z-direction.'
      END IF
    END IF

    WRITE(11,*) ''

  END DO
END IF
 101  FORMAT ('Inlet number ',I4,' is called ',A32)
 102  FORMAT ('Inlet number ',I4,' has inlet value 1 = ',&
              F8.4)
 103  FORMAT ('Inlet number ',I4,' has inlet value 2 = ',&
              F8.4)
 104  FORMAT ('Inlet number ',I4,' is a momentum',&
                ' boundary condition.')
 105  FORMAT ('Inlet number ',I4,' is a pressure',&
                ' boundary condition.')
 106  FORMAT ('Inlet number ',I4,' is a freestream',&
                ' boundary condition.')
 107  FORMAT ('Inlet number ',I4,' is a temperature and',&
                ' velocity boundary condition.')
 108  FORMAT ('Inlet number ',I4,' is a periodic',&
                ' boundary condition.')
 109  FORMAT ('Inlet number ',I4,' is an invisid wall',&
                ' boundary condition.')
 110  FORMAT ('Inlet number ',I4,' is a viscous wall',&
                ' boundary condition.')
 122  format ('Inlet number ',I4,' directional flow boundary condition.')
! 101  FORMAT ()
!
! Make sure there's at least one outlet and go on
! Loops through the outlet data reading as many as there are.
!
      IF (num_outlets < 1) THEN
        WRITE(11,*) 'Too few outlets, you need at least 1 outlet, exiting.'
        CALL error_out
      ELSE
        DO i =1,num_outlets
          outlet_val_1 = 0.0D0
          outlet_val_2 = 0.0D0
          outlet_val_3 = 0.0D0
          outlet_val_4 = 0.0D0
          outlet_val_5 = 0.0D0
          outlet_name = 'Basic Outlet'
          out_type = 1
          outlet_flow_dir = 4
 
          READ(2,NML=OUTLET_DATA,IOSTAT=istat)
          IF (istat /= 0) THEN
            WRITE(*,*) 'Error reading outlet data at outlet ',i,', exiting.'
            WRITE(11,*) 'Error reading outlet data at outlet ',i,', exiting.'
            WRITE(*,*) istat
            CALL error_out
          ELSE
            outlet_flow_val(1,i) = outlet_val_1
            outlet_flow_val(2,i) = outlet_val_2
            outlet_flow_val(3,i) = outlet_val_3
            outlet_flow_val(4,i) = outlet_val_4
            outlet_flow_val(5,i) = outlet_val_5
            out_name(i) = outlet_name
            outlet_type(1,i) = out_type
            outlet_type(2,i) = outlet_flow_dir

            IF (outlet_type(2,i) < 0 .OR. outlet_type(2,i) > 6) THEN
              WRITE(*,*) 'Outlet flow direction out of range, must be between 1 and 6.'
              CALL error_out
            END IF


            outlet_side(outlet_flow_dir) = .TRUE.
!
! This will do the output to the output file, essentially reiterating
! everything the user puts in for the outlets.
!
            WRITE(11,201) i,out_name(i)

            IF (outlet_type(1,i) == 1) THEN
              WRITE(11,204) i
            ELSE IF (outlet_type(1,i) == 2) THEN
              WRITE(11,205) i
            ELSE IF (outlet_type(1,i) == 3) THEN
              WRITE(11,206) i
            ELSE IF (outlet_type(1,i) == 4) THEN
              WRITE(11,207) i
            ELSE IF (outlet_type(1,i) == 5) THEN
              WRITE(11,208) i
            else if (outlet_type(1,i) == 6) then
              write(11,222) i
            END IF

            IF (outlet_type(2,i) == 1) THEN
              WRITE(11,*) 'Outlet is in the positive x-direction.'
            ELSE IF (outlet_type(2,i) == 2) THEN
              WRITE(11,*) 'Outlet is in the positive y-direction.'
            ELSE IF (outlet_type(2,i) == 3) THEN
              WRITE(11,*) 'Outlet is in the positive z-direction.'
            ELSE IF (outlet_type(2,i) == 4) THEN
              WRITE(11,*) 'Outlet is in the negative x-direction.'
            ELSE IF (outlet_type(2,i) == 5) THEN
              WRITE(11,*) 'Outlet is in the negative y-direction.'
            ELSE IF (outlet_type(2,i) == 6) THEN
              WRITE(11,*) 'Outlet is in the negative z-direction.'
            END IF

          END IF

          WRITE(11,202) i,outlet_flow_val(1,i)
          WRITE(11,203) i,outlet_flow_val(2,i)
          WRITE(11,203) i,outlet_flow_val(3,i)
          WRITE(11,203) i,outlet_flow_val(4,i)
          WRITE(11,203) i,outlet_flow_val(5,i)

!          WRITE(*,202) i,outlet_flow_val(1,i)
!          WRITE(*,203) i,outlet_flow_val(2,i)
!          WRITE(*,203) i,outlet_flow_val(3,i)
!          WRITE(*,203) i,outlet_flow_val(4,i)
!          WRITE(*,203) i,outlet_flow_val(5,i)

          WRITE(11,*) ''

        END DO
      END IF  

      DO i = 1,6
        IF (inlet_side(i) .AND. outlet_side(i)) THEN
          WRITE(11,*) 'Inlet and outlet are on the same side of the',&
            'bounding box. Exiting gracefully.'
          CALL error_out
        END IF
      END DO
      WRITE(11,*) ''
 201  FORMAT ('Outlet number ',I4,' is called ',A32)
 202  FORMAT ('Outlet number ',I4,' has outlet value 1 = ',&
              F8.4)
 203  FORMAT ('Outlet number ',I4,' has outlet value 2 = ',&
              F8.4)
 204  FORMAT ('Outlet number ',I4,' is a momentum',&
                ' boundary condition.')
 205  FORMAT ('Outlet number ',I4,' is a pressure',&
                ' boundary condition.')
 206  FORMAT ('Outlet number ',I4,' is a freestream',&
                ' boundary condition.')
 207  FORMAT ('Outlet number ',I4,' is a temperature and',&
                ' velocity boundary condition.')
 208  FORMAT ('Outlet number ',I4,' is a periodic',&
                ' boundary condition.')
! 209  FORMAT ('Outlet number ',I4,' is an invisid wall',&
!                ' boundary condition.')
! 210  FORMAT ('Outlet number ',I4,' is a viscous wall',&
!                ' boundary condition.')
 222  format ('Outlet number ',I4,' is a pure outlet boundary condition.')

      END SUBROUTINE
