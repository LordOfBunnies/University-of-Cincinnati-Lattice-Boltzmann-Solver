      SUBROUTINE gause_elim_and_newton_raph(jacobian,em,deter,&
        cramer_1,cramer_2,cramer_3,cramer_4,cramer_5)
!
!
!
!
!
!
!
      USE precise
      USE constants
      
      IMPLICIT NONE
      REAL(kind=dp) :: jacobian(em,em)
      REAL(kind=dp) :: deter,cramer_1,cramer_2,cramer_3,cramer_4,cramer_5 
      INTEGER :: em,i,j,k

      DO i = 0,em






        SELECT CASE (i)
          CASE (0)

      END DO


      END SUBROUTINE
