      SUBROUTINE ggm_output(a,t)
!
! Output ggm to a CGNS file that was prepared and started before.
!
! Calls: cg_open_f,cg_sol_write_f,cg_field_write_f,cg_error_exit_f
!          cg_close_f
! Called by: ggm_calculator
!
      USE precise
      USE constants
      USE cgns
      USE ggm_stuff
      USE freestream_values
      IMPLICIT NONE
      INTEGER :: a,t
      INTEGER :: index_base,index_zone,index_file,&
        ier,index_flow,index_field
!      REAL(KIND=dp),ALL ::
      CHARACTER(len=32) :: filename,char_save_number,format_junk,&
        solut_name,format_junk_ts,solut_save_number

        index_base = 1
!        index_file = 1
        index_zone = 1
        index_flow = 1
        index_field = 1


      format_junk = '(I3.3)'
      format_junk_ts = '(I9.9)'
      WRITE(solut_save_number,format_junk_ts) t

      WRITE(char_save_number,format_junk) a

!        WRITE(*,*) SIZE(x_vertices),SIZE(y_vertices),SIZE(z_vertices)

      filename = 'ggm_output_'//TRIM(char_save_number)//'.cgns'
      solut_name = 'Timestep_'//TRIM(solut_save_number)
!      basename = 'outputbase_'//TRIM(char_save_number)
!      zonename = 'outputzone_'//TRIM(char_save_number)

      CALL cg_open_f(filename,CG_MODE_MODIFY,index_file,ier)
      IF (ier /= CG_OK) CALL cg_error_exit_f
!
! Write the solution for this timestep
!
!        index_file = 1
!        WRITE(*,*) ''
!        WRITE(*,*) index_file,index_base,index_zone,solut_name,&
!          index_flow,ier
      CALL cg_sol_write_f(index_file,index_base,index_zone,solut_name,&
        CellCenter,index_flow,ier)
      IF (ier /= CG_OK) CALL cg_error_exit_f
!        index_file = 1
!
! Write density
!
!      IF (save_data(1,a)) THEN
       CALL cg_field_write_f(index_file,index_base,index_zone,index_flow,&
         RealDouble,'GGM',ggm,index_field,ier)
!        WRITE(11,*) rho_local
       IF (ier /= CG_OK) CALL cg_error_exit_f
!      END IF

      CALL cg_close_f(index_file,ier)
      IF (ier /= CG_OK) CALL cg_error_exit_f

      END SUBROUTINE
