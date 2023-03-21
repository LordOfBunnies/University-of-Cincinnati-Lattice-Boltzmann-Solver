SUBROUTINE calculate_q_criterion
!
! Calculates the q-criterion, a vortex identification factor
!
! Calls: calculate_vort,calculate_shear_stress
! Called by: populate_q
!
!      USE ggm_stuff
!      USE precise
!      USE constants
!      USE derivative_form
!      USE linkwise
IMPLICIT NONE
INTEGER :: a,choice

!
! Q-criteria is a method of determining vortex locations
! In essence, it is the strain rate tensor + vorticity tensor
! Q_ij = S_ij+Omega_ij
! S_ij = .5*(du_i/dx_j+du_j/dx_i)
! Omega_ij = .5*(du_i/dx_j-du_j/dx_i)
!
! Essentially, it ensure that the flow is rotating and not simply
! straining the local space.
!
!      ALLOCATE(q2(link_number))
!
!      choice = 4
!      CALL calculate_vort(choice)
!
!      q2 = q
!
!      choice = 4
!      CALL calculate_shear_stress(choice)
!
!      DO a = 1, link_number
!
!
!        q(a) = q(a)-q2(a)
!
!
!      END DO

END SUBROUTINE
