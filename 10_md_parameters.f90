!-----------------------------------------------------
PROGRAM main
!-----------------------------------------------------
  USE ions_base, ONLY: amass, nsp
  USE dynamics_module, ONLY: dt
  USE constants, ONLY: amu_ry
  IMPLICIT NONE
  INTEGER :: isp

  CALL prepare_all()

  DO isp = 1,nsp
    WRITE(*,*) 'amass = ', amass(isp)*amu_ry*2.d0
  ENDDO
  WRITE(*,*) 'dt = ', dt

END PROGRAM