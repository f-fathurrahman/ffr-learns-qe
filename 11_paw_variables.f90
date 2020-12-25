PROGRAM main
  USE kinds, ONLY: DP
  IMPLICIT NONE 

  CALL prepare_all()
  CALL test_paw_variables()

END PROGRAM 

SUBROUTINE test_paw_variables()
  USE paw_variables, ONLY: total_core_energy, okpaw
  USE uspp, ONLY: okvan

  WRITE(*,*) 'okpaw = ', okpaw
  WRITE(*,*) 'okvan = ', okvan
  WRITE(*,*) 'total_core_energy = ', total_core_energy

END SUBROUTINE